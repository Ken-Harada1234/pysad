# -*- coding: utf-8 -*-
"""12_gui_tkinter の例: pysad が開いた窓を Windows の画面から切り取って PNG にする。

使い方: py -3.12 example/tools/gui_shot.py example/12_gui_tkinter/fn_X.sad example/12_gui_tkinter/fn_X_window.png
.sad の最後の Exit; を外した台本を pysad で走らせ、tkinter のルート窓の位置と大きさを取り、PIL の ImageGrab で画面から切り取る。
(pysad では Window[] が別の窓を作らず、部品はルート窓 `.` の下に並ぶ。docs/22 の既知の差。)
gs の窓 (WSLg の X11) は撮らない (撮る道具が WSL に無い)。
"""
import os, sys, time
ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
# 環境変数 PYSAD_SAD_PACKAGES (K64 版 SAD の Packages のフォルダ) を指定して実行する
sys.path.insert(0, ROOT)


def shot(sad, png, wait=1.0):
    src = open(sad, encoding='utf-8').read()
    i = src.rstrip().rfind('Exit;')
    src = src[:i]
    from pysad.ffscmd import run_sad_string
    run_sad_string(src)
    import tkinter
    root = tkinter._default_root
    root.update()
    root.attributes('-topmost', True)
    root.lift(); root.update(); time.sleep(wait); root.update()
    x, y, w, h = root.winfo_rootx(), root.winfo_rooty(), root.winfo_width(), root.winfo_height()
    img = None
    try:
        from PIL import ImageGrab
        img = ImageGrab.grab(bbox=(x - 2, y - 2, x + w + 2, y + h + 2))
    except OSError:
        img = _printwindow(root)            # 画面の切り取りが使えない環境 (デスクトップ無し) では PrintWindow
    img.save(png, optimize=True)
    print('saved', png, img.size, (x, y, w, h))


def _printwindow(root):
    import ctypes
    from ctypes import wintypes
    from PIL import Image
    u, g = ctypes.windll.user32, ctypes.windll.gdi32
    hwnd = u.GetParent(root.winfo_id()) or root.winfo_id()
    r = wintypes.RECT(); u.GetWindowRect(hwnd, ctypes.byref(r))
    w, h = r.right - r.left, r.bottom - r.top
    hdc = u.GetDC(0); mdc = g.CreateCompatibleDC(hdc); bmp = g.CreateCompatibleBitmap(hdc, w, h); g.SelectObject(mdc, bmp)
    ok = u.PrintWindow(hwnd, mdc, 2)
    class BIH(ctypes.Structure):
        _fields_ = [('biSize', wintypes.DWORD), ('biWidth', ctypes.c_long), ('biHeight', ctypes.c_long), ('biPlanes', wintypes.WORD), ('biBitCount', wintypes.WORD),
                    ('biCompression', wintypes.DWORD), ('biSizeImage', wintypes.DWORD), ('x', ctypes.c_long), ('y', ctypes.c_long), ('c', wintypes.DWORD), ('i', wintypes.DWORD)]
    bi = BIH(ctypes.sizeof(BIH), w, -h, 1, 32, 0, 0, 0, 0, 0, 0)
    buf = ctypes.create_string_buffer(w * h * 4)
    g.GetDIBits(mdc, bmp, 0, h, buf, ctypes.byref(bi), 0)
    img = Image.frombuffer('RGBA', (w, h), buf, 'raw', 'BGRA', 0, 1).convert('RGB')
    print('PrintWindow ok=%s size=%s' % (ok, img.size))
    return img


if __name__ == '__main__':
    shot(sys.argv[1], sys.argv[2])
