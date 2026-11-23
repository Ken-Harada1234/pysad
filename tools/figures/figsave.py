# -*- coding: utf-8 -*-
"""図の保存 (折れ線系は SVG、そうでなければ PNG)。図のスクリプト (tools/figures と example の figures_make_*.py) が使う。

保存の方針 (`save_figure`):
- 図を SVG (文字はパスにして、フォントのない環境でも崩れない) と PNG の両方に描いて大きさを比べる。
  SVG が PNG 以下、または 150 KB 以下なら SVG だけを置く (同名の古い PNG は消す)。そうでなければ PNG だけを置く。
- 点が非常に多い図 (散布図・画像) は `vector=False` を渡す。SVG を試さず、PNG だけを置く。
- PNG は、画素を変えずに小さくする (zlib 最大、フィルターは画像ごとに最良のものを選ぶ)。
- 環境変数 `PYSAD_FIG_SVG_ONLY=1` のときは、PNG にする図を書かない (SVG にできるかの確認だけをしたいとき)。
"""
import io
import os
import pathlib
import struct
import zlib

import matplotlib

matplotlib.rcParams["svg.fonttype"] = "path"     # 文字をパスにする (フォントがなくても同じ見た目)
matplotlib.rcParams["svg.hashsalt"] = "pysad"    # ID を固定して、同じ図なら同じ SVG にする
SVG_LIMIT = 150 * 1024


def _chunk(tag, data):
    c = struct.pack(">I", len(data)) + tag + data
    return c + struct.pack(">I", zlib.crc32(tag + data) & 0xFFFFFFFF)


def optimize_png(data):
    """PNG のバイト列を、画素が同じまま小さく作り直す。小さくならなければ元のまま返す。
    RGBA で不透明なら RGB に、色が 256 以下なら 8 ビットパレットにする (どちらも画素は同じ)。"""
    import numpy as np
    from PIL import Image
    im = Image.open(io.BytesIO(data))
    im.load()
    if im.mode not in ("RGB", "RGBA", "L"):
        return data
    a = np.array(im)
    variants = []
    if im.mode == "RGBA" and (a[..., 3] == 255).all():
        a = a[..., :3]
    if a.ndim == 3 and a.shape[2] == 3:
        flat = a.reshape(-1, 3)
        key = (flat[:, 0].astype(np.int64) << 16) | (flat[:, 1].astype(np.int64) << 8) | flat[:, 2]
        uk = np.unique(key)
        if len(uk) <= 256:
            pal = np.stack([(uk >> 16) & 255, (uk >> 8) & 255, uk & 255], 1).astype(np.uint8)
            variants.append((np.searchsorted(uk, key).reshape(a.shape[:2]).astype(np.uint8), 3, pal))
    variants.append((a, {1: 0, 3: 2, 4: 6}[1 if a.ndim == 2 else a.shape[2]], None))
    best = data
    for arr, ctype, pal in variants:
        bpp = 1 if arr.ndim == 2 else arr.shape[2]
        h = arr.shape[0]
        rb = arr.reshape(h, -1)
        for f in (0, 1, 2):
            if f == 0:
                filt = rb
            elif f == 1:
                filt = rb.copy()
                filt[:, bpp:] = (rb[:, bpp:].astype(np.int16) - rb[:, :-bpp]).astype(np.uint8)
            else:
                filt = rb.copy()
                filt[1:] = (rb[1:].astype(np.int16) - rb[:-1]).astype(np.uint8)
            raw = b"".join(bytes([f]) + filt[y].tobytes() for y in range(h))
            z = zlib.compress(raw, 9)
            out = b"\x89PNG\r\n\x1a\n" + _chunk(b"IHDR", struct.pack(">IIBBBBB", arr.shape[1], h, 8, ctype, 0, 0, 0))
            if pal is not None:
                out += _chunk(b"PLTE", pal.tobytes())
            out += _chunk(b"IDAT", z) + _chunk(b"IEND", b"")
            if len(out) < len(best):
                best = out
    return best


def shrink_svg(svg):
    """SVG を小さくする: 数値を小数 2 桁に丸める (1/100 ポイントの単位で、見た目は変わらない)。作成情報・コメント・字下げを除く。"""
    import re
    t = svg.decode("utf-8")
    t = re.sub(r"(?<![\w.])(-?\d+\.\d{3,})", lambda m: ("%.2f" % float(m.group(1))).rstrip("0").rstrip("."), t)
    t = re.sub(r"<metadata>.*?</metadata>\s*", "", t, flags=re.S)       # 作成情報
    t = re.sub(r"<!--.*?-->\s*", "", t, flags=re.S)                      # コメント
    t = re.sub(r"\n[ \t]+", "\n", t)                                     # 字下げ
    t = re.sub(r"[ \t]+\n", "\n", t)
    return t.encode("utf-8")


def save_figure(fig, path, dpi=115, vector=True, png_max_kb=None, **kw):
    """図を保存する。path は .png で与える。SVG にしたときは .svg を書き、書いたパスを返す。
    png_max_kb を与えると、PNG がそれを超える間 dpi を 10 ずつ下げる (60 まで)。"""
    path = pathlib.Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    svg_path = path.with_suffix(".svg")
    kw.setdefault("bbox_inches", "tight")
    while True:
        b = io.BytesIO()
        fig.savefig(b, format="png", dpi=dpi, **kw)
        png = optimize_png(b.getvalue())
        if png_max_kb is None or len(png) <= png_max_kb * 1024 or dpi <= 60:
            break
        dpi -= 10
    svg = None
    if vector:
        b = io.BytesIO()
        fig.savefig(b, format="svg", metadata={"Date": None}, **kw)
        svg = shrink_svg(b.getvalue())
    import matplotlib.pyplot as plt
    plt.close(fig)
    if svg is not None and (len(svg) <= len(png) or len(svg) <= SVG_LIMIT):
        svg_path.write_bytes(svg)
        if path.exists():
            path.unlink()
        print(svg_path, "SVG", len(svg) // 1024, "KB (PNG なら", len(png) // 1024, "KB)")
        return svg_path
    if os.environ.get("PYSAD_FIG_SVG_ONLY") == "1":
        print(path, "PNG のまま (SVG", "-" if svg is None else len(svg) // 1024, "KB / PNG", len(png) // 1024, "KB) 書かない")
        return None
    path.write_bytes(png)
    if svg_path.exists():
        svg_path.unlink()
    print(path, "PNG", len(png) // 1024, "KB", "(SVG", "-" if svg is None else "%d KB" % (len(svg) // 1024), ")")
    return path
