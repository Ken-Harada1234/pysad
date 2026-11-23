# POL + RAD + TRPT + 加速 (track6 の経路) のスピンの歳差が、断熱加速の解析式 Φ = aγ1θ1 + aγ2θ2 になる。
# 1 GeV、扇形の偏向磁石 (θ1 = 0.3、θ2 = 0.2) の間に加速する空洞 (VOLT = 0・0.2・0.5 GV、PHI = -PI/2)。縦向きのスピンの (sx, sz) = (sin Φ, cos Φ)。
# 使い方: py -3.12 example/16_spin_polarization/spin_accel.py   (v2 モードだけ。k64 は SADError)
import io, math, os, sys
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..'))
from pysad import ffscmd, emit6d, compat
compat.set_mode('v2')
SRC = '''MOMENTUM = 1 GEV;
DRIFT D = (L=0.5);
BEND B1 = (L=1.0 ANGLE=0.3 E1=0 E2=0) B2 = (L=1.0 ANGLE=0.2 E1=0 E2=0);
CAVI RF = (L=0 VOLT=%g FREQ=5e8 PHI=-PI/2);
MARK ST = (BX=5 BY=5);
LINE BL = (ST D B1 D RF D B2 D);
FFS;
USE BL;
TRPT; INS; RFSW; RAD; NOFLUC; CALC;
Exit;
'''
a = 0.00115965218
for volt in (0.0, 2e8, 5e8):
    ffs2 = ffscmd.run_sad_string(SRC % volt, out=io.StringIO())[0]
    ffs2.flags['RAD'] = True
    ffs2.flags['FLUC'] = False
    ffs2.flags['POL'] = True
    n = len(ffs2.els) + 1
    cols = [[0.0]] * 6 + [[math.atan2(1.0, 0.0)], [0.0], [1.0]]
    r = ffs2.fn_track([[1.0, cols], float(n)])[1]
    ang, sy = r[6][0], r[7][0]
    sx, sz = math.sqrt(1 - sy * sy) * math.cos(ang), math.sqrt(1 - sy * sy) * math.sin(ang)
    ph = emit6d._Phys(ffs2.momentum(), 10.0)
    h1 = ph.h0
    h2 = h1 + volt / ph.amass
    phi = a * (h1 * 0.3 + h2 * 0.2)
    print(f'VOLT = {volt:.1e} V: γ1 = {h1:.2f}、γ2 = {h2:.2f}、Φ = {phi:.6f} rad: 追跡 (sx, sz) = ({sx:+.9f}, {sz:+.9f})、式 (±sinΦ, cosΦ) = ({math.sin(phi):.9f}, {math.cos(phi):.9f})')
