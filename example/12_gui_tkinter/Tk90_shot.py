# -*- coding: utf-8 -*-
"""Tk 9.0 (Python 3.14) の GUI の見た目を撮る。PIL なしで動く (PrintWindow + 自前の PNG 書き出し)。

使い方 (リポジトリ直下から):
  <py314 の venv>/Scripts/python.exe example/12_gui_tkinter/Tk90_shot.py example/12_gui_tkinter/fn_Canvas.sad out.png
  py -3.12 example/12_gui_tkinter/Tk90_shot.py ...          (Python 3.12 = Tk 8.6 の窓と比べる)
.sad の最後の Exit; を外して pysad で走らせ、tkinter のルート窓を PrintWindow で取る。標準出力に Python・Tk の版と窓の大きさ、窓の中の部品の数を出す。
"""
import ctypes
import os
import struct
import sys
import time
import zlib
from ctypes import wintypes

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
# 環境変数 PYSAD_SAD_PACKAGES (K64 版 SAD の Packages のフォルダ) を指定して実行する
sys.path.insert(0, ROOT)


def write_png(path, w, h, rgb):
    raw = b''.join(b'\x00' + rgb[y * w * 3:(y + 1) * w * 3] for y in range(h))
    def chunk(t, d):
        c = struct.pack('>I', len(d)) + t + d
        return c + struct.pack('>I', zlib.crc32(t + d) & 0xffffffff)
    with open(path, 'wb') as f:
        f.write(b'\x89PNG\r\n\x1a\n' + chunk(b'IHDR', struct.pack('>IIBBBBB', w, h, 8, 2, 0, 0, 0)) + chunk(b'IDAT', zlib.compress(raw, 9)) + chunk(b'IEND', b''))


def printwindow(root):
    u, g = ctypes.windll.user32, ctypes.windll.gdi32
    hwnd = u.GetParent(root.winfo_id()) or root.winfo_id()
    r = wintypes.RECT()
    u.GetWindowRect(hwnd, ctypes.byref(r))
    w, h = r.right - r.left, r.bottom - r.top
    hdc = u.GetDC(0)
    mdc = g.CreateCompatibleDC(hdc)
    bmp = g.CreateCompatibleBitmap(hdc, w, h)
    g.SelectObject(mdc, bmp)
    ok = u.PrintWindow(hwnd, mdc, 2)

    class BIH(ctypes.Structure):
        _fields_ = [('biSize', wintypes.DWORD), ('biWidth', ctypes.c_long), ('biHeight', ctypes.c_long), ('biPlanes', wintypes.WORD), ('biBitCount', wintypes.WORD),
                    ('biCompression', wintypes.DWORD), ('biSizeImage', wintypes.DWORD), ('x', ctypes.c_long), ('y', ctypes.c_long), ('c', wintypes.DWORD), ('i', wintypes.DWORD)]
    bi = BIH(ctypes.sizeof(BIH), w, -h, 1, 32, 0, 0, 0, 0, 0, 0)
    buf = ctypes.create_string_buffer(w * h * 4)
    g.GetDIBits(mdc, bmp, 0, h, buf, ctypes.byref(bi), 0)
    b = bytes(buf)
    rgb = bytearray(w * h * 3)
    rgb[0::3] = b[2::4]
    rgb[1::3] = b[1::4]
    rgb[2::3] = b[0::4]
    return ok, w, h, bytes(rgb)


def count_widgets(w):
    n = 1
    for c in w.winfo_children():
        n += count_widgets(c)
    return n


def main(sad, png):
    import tkinter
    print('python', sys.version.split()[0], 'Tk', tkinter.TkVersion)
    src = open(sad, encoding='utf-8').read()
    i = src.rstrip().rfind('Exit;')
    src = src[:i]
    from pysad.ffscmd import run_sad_string
    run_sad_string(src)
    root = tkinter._default_root
    root.update()
    root.attributes('-topmost', True)
    root.lift()
    root.update()
    time.sleep(1.0)
    root.update()
    ok, w, h, rgb = printwindow(root)
    write_png(png, w, h, rgb)
    nonwhite = sum(1 for k in range(0, len(rgb), 3) if rgb[k] != 255 or rgb[k + 1] != 255 or rgb[k + 2] != 255)
    print('window', w, 'x', h, 'PrintWindow ok =', ok, 'widgets', count_widgets(root), 'non-white pixels', nonwhite, '->', png)


if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2])
