# -*- coding: utf-8 -*-
"""PF の履歴ラティスの解析の共通部 (pysad を直接呼ぶ)。run(lat, mode) -> (ffs, 画面表示)"""
import os, sys, io, re
HERE = os.path.dirname(os.path.abspath(__file__)); EX = os.path.dirname(os.path.dirname(HERE)); ROOT = os.path.dirname(EX)
sys.path.insert(0, ROOT)
LATS = [('pf_1981_460nmrad', 1981, 460), ('pf_1986_128nmrad', 1986, 128), ('pf_1997_036nmrad', 1997, 36), ('pf_2005_036nmrad', 2005, 36)]
def run(lat, mode='k64', body='USE RING;\ncell calc emit;\n'):
    os.environ['PYSAD_COMPAT'] = mode
    from pysad import compat
    compat.set_mode(mode)
    from pysad import ffscmd
    src = open(os.path.join(EX, 'pf', 'lattice', lat + '.sad'), encoding='utf-8').read()
    buf = io.StringIO()
    ffs, I = ffscmd.run_sad_string(src + '\nON LOG;\nON RAD FLUC COD RFSW RADCOD;\nFFS;\n$FORM = "15.10"; PageWidth = 1999;\n' + body + 'Exit;\n', out=buf)
    return ffs, buf.getvalue()
def emit_table(txt):
    """emit の表の画面表示から値を取り出す"""
    g = lambda pat: float(re.search(pat, txt).group(1))
    d = {}
    d['U0'] = g(r'Energy loss per turn U0 =\s*([-\d.E+]+) MV')
    d['Vc'] = g(r'Effective voltage\s+Vc =\s*([-\d.E+]+) MV')
    d['alpha'] = g(r'Momentum compact\. alpha =\s*([-\d.E+]+)')
    d['h'] = g(r'Effective harmonic #\s+h =\s*([-\d.E+]+)')
    d['f0'] = g(r'f0 =\s*([-\d.E+]+) Hz')
    d['fs'] = g(r'Synchrotron frequency\s+=\s*([-\d.E+]+) Hz')
    d['bucket'] = g(r'Bucket height\s+dV/P0 =\s*([-\d.E+]+)')
    m = re.search(r'Damping time \(sec\):\s*X :\s*([-\d.E+]+)\s+Y :\s*([-\d.E+]+)\s+Z :\s*([-\d.E+]+)', txt)
    d['taux'], d['tauy'], d['tauz'] = map(float, m.groups())
    m = re.search(r'Damping partition number:\s*X :\s*([-\d.E+]+)\s+Y :\s*([-\d.E+]+)\s+Z :\s*([-\d.E+]+)', txt)
    d['Jx'], d['Jy'], d['Jz'] = map(float, m.groups())
    d['ex'] = g(r'Emittance X\s+=\s*([-\d.E+]+) m')
    d['sige'] = g(r'Energy spread\s+=\s*([-\d.E+]+)')
    d['sigz'] = g(r'Bunch Length\s+=\s*([-\d.E+]+) mm')
    return d
