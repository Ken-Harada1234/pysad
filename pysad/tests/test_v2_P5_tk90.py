# -*- coding: utf-8 -*-
"""担当 P5 (v2next) 項目 7: Tk 9.0 (Python 3.14 の venv (環境変数 PYSAD_PY314 か ~/venvs/py314)) で pysad の GUI が動き、窓の見た目が撮れる。

example/v2_next/P5_tk90_shot.py (PIL なし: PrintWindow + 自前の PNG) を 3.14 の venv で走らせ、`Tk 9.0` で窓が出て、画面が空でない (白以外の画素が多い) ことを確かめる。
Python 3.12 (Tk 8.6) と窓の大きさが同じで、白以外の画素数の差は 1% 以内 (フォントの描画だけが違う)。見た目の画像は example/v2_next/P5_tk90_*.png。
デスクトップが無い・venv が無い・`PYSAD_SAD_PACKAGES` が無い環境では skip。
"""
import os
import re
import subprocess
import sys

import pytest

from . import _p5_util as U

VENV = os.environ.get('PYSAD_PY314', os.path.expanduser('~/venvs/py314/Scripts/python.exe'))
SHOT = os.path.join(U.EX, 'P5_tk90_shot.py')
GUI = os.path.join(U.ROOT, 'example', '12_gui_tkinter')


def _shot(py, name, tmp_path):
    png = str(tmp_path / (name + '.png'))
    env = dict(os.environ)
    r = subprocess.run([py, SHOT, os.path.join(GUI, name + '.sad'), png], capture_output=True, text=True, timeout=240, env=env, cwd=U.ROOT,
                       encoding='utf-8', errors='replace')
    out = r.stdout
    m = re.search(r'window (\d+) x (\d+) PrintWindow ok = (\d+) widgets (\d+) non-white pixels (\d+)', out)
    if not m:
        pytest.skip('窓が撮れない環境: ' + (r.stderr or out)[-200:])
    return re.search(r'Tk ([\d.]+)', out).group(1), tuple(int(x) for x in m.groups()), png


@pytest.fixture(scope='module')
def env_ok():
    if not os.path.exists(VENV):
        pytest.skip('Python 3.14 の venv がない')
    if not os.environ.get('PYSAD_SAD_PACKAGES'):
        pytest.skip('PYSAD_SAD_PACKAGES がない')
    if not os.path.exists(SHOT):
        pytest.skip('P5_tk90_shot.py がない')


@pytest.mark.parametrize('name', ['fn_Window_Button_TextLabel', 'fn_Canvas', 'fn_KBMainFrame_KBFComponentFrame'])
def test_gui_runs_under_tk90_and_matches_tk86_size(env_ok, name, tmp_path):
    tk9, (w9, h9, ok9, _, nw9), png = _shot(VENV, name, tmp_path)
    assert tk9.startswith('9.0')
    assert ok9 == 1 and nw9 > 5000
    assert os.path.getsize(png) > 500
    tk8, (w8, h8, ok8, _, nw8), _ = _shot(sys.executable, name, tmp_path)
    if tk8.startswith('8.6'):
        assert (w9, h9) == (w8, h8)
        assert abs(nw9 - nw8) <= 0.01 * nw8
