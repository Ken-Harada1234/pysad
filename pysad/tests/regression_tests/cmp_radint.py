# -*- coding: utf-8 -*-
"""cmp_radint.py SAD.result PYSAD.out [tol] — RADINT の出力ブロック (見出し行 '****' から 'LINEAR CHROMATICITY GY' まで) を
ブロックごと・行ごとに突き合わせる。数値は相対 tol (既定 1e-6。SAD の表示が 7 桁) で比べ、両方とも |x|<1e-10 (エミッタンス行は 1e-14) なら雑音として一致扱い。
行の文字部分 (数を除いた形) は完全一致を要求する。不一致の行だけ表示し、最後に 一致/不一致 の数を出す。"""
import re
import sys

NUM = re.compile(r'[-+]?(?:\d+\.\d*|\.\d+|\d+)(?:[eEdD][-+]?\d+)?')


def blocks(path):
    return blocks_text(open(path, encoding='utf-8', errors='replace').read())


def blocks_text(text):
    lines = text.split('\n')
    out, cur = [], None
    for l in lines:
        if l.startswith('*' * 20):
            if cur is None:
                cur = [l]
                continue
        if cur is not None:
            cur.append(l)
            if 'LINEAR CHROMATICITY GY' in l:
                out.append(cur)
                cur = None
    return out


def fnum(s):
    return float(s.replace('D', 'e').replace('d', 'e'))


def main(a, b, tol=1e-6, quiet=False):
    return compare(blocks(a), blocks(b), tol, quiet)


def compare(A, B, tol=1e-6, quiet=False):
    """A, B: blocks() の結果 (SAD, pysad)。戻り (一致行数, 不一致行数, 最大相対差)。tol はブロックごとのリストでもよい"""
    ok = bad = 0
    worst = 0.0
    if len(A) != len(B):
        if not quiet:
            print(f'ブロック数が違う: SAD {len(A)} / pysad {len(B)}')
        bad += abs(len(A) - len(B))
    for k, (ba, bb) in enumerate(zip(A, B), 1):
        tolk = tol[k - 1] if isinstance(tol, (list, tuple)) else tol
        ba = [l for l in ba if l.strip()]
        bb = [l for l in bb if l.strip()]
        if len(ba) != len(bb):
            bad += 1
            if not quiet:
                print(f'[{k}] 行数が違う: {len(ba)} / {len(bb)}')
        for la, lb in zip(ba, bb):
            na, nb = NUM.findall(la), NUM.findall(lb)
            ta, tb = ' '.join(NUM.sub('#', la).split()), ' '.join(NUM.sub('#', lb).split())
            good = (ta == tb and len(na) == len(nb))
            if good:
                for x, y in zip(na, nb):
                    x, y = fnum(x), fnum(y)
                    # 0 同然の量 (非結合の V 面など。SAD は 1e-35〜1e-11 の丸め雑音) は雑音として一致扱い。エミッタンスは桁が小さいので除く
                    thr = 1e-14 if 'EMITTANCE' in la else 1e-10
                    if abs(x) < thr and abs(y) < thr:
                        continue
                    r = abs(x - y) / max(abs(x), 1e-300)
                    worst = max(worst, r)
                    if r > tolk:
                        good = False
            if good:
                ok += 1
            else:
                bad += 1
                if not quiet:
                    print(f'[{k}] SAD  : {la}\n[{k}] pysad: {lb}')
    if not quiet:
        print(f'一致 {ok} 行 / 不一致 {bad} 行 (最大相対差 {worst:.2e}, 許容 {tol})')
    return ok, bad, worst


if __name__ == '__main__':
    sys.exit(1 if main(sys.argv[1], sys.argv[2], float(sys.argv[3]) if len(sys.argv) > 3 else 1e-6)[1] else 0)
