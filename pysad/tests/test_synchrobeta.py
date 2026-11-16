"""SynchroBetaEmittance (temits.f)。参照は SAD K64 (gs): 12 セルのリング + CAVI (VOLT=2e6 HARM=64)、CELL; RAD; CALC のあと。
pysad は運動量のずれごとの閉軌道・1 周写像と、(アクション, 方位角モード) の連立方程式の反復を直訳した。gs と相対 1e-5 以内で一致する (小さな差は閉軌道の収束の履歴など)。
SynchroBetaEmittance は Packages/emit.n (PYSAD_SAD_PACKAGES が必要)。"""
import io
import os
import re

import pytest

PK = os.environ.get('PYSAD_SAD_PACKAGES')
needs_pk = pytest.mark.skipif(not PK or not os.path.isdir(PK), reason='環境変数 PYSAD_SAD_PACKAGES (SAD の Packages) が必要')

SRC = r"""MOMENTUM = 3GEV;
MARK ST=();
DRIFT D1 = (L=1) D2 = (L=0.5);
QUAD QF = (L=0.3 K1=0.7) QD = (L=0.3 K1=-0.6);
BEND B = (L=1.2 ANGLE=0.2617993877991494 E1=0.1308996938995747 E2=0.1308996938995747);
CAVI RFC = (VOLT=2e6 HARM=64);
LINE CELL1 = (QF D2 B D2 QD D1);
LINE RING = (ST 12*CELL1 RFC);
FFS;
PageWidth=1999;
USE RING; CELL; RFSW; RAD; CALC;
"""


def _run(body):
    from pysad.ffscmd import run_sad_string
    buf = io.StringIO()
    run_sad_string(SRC + body + 'Exit;\n', out=buf)
    return buf.getvalue()


def _nums(line):
    return [float(x) for x in re.findall(r'[-+]?\d*\.?\d+(?:e[-+]?\d+)?', line)]


@needs_pk
def test_synchrobeta_matches_gs():
    out = _run('r=SynchroBetaEmittance[{0.005,0.02,0.005},AzimuthalModes->5];\nPrint["R ",r];\n')
    line = [l for l in out.splitlines() if l.startswith('R ')][0]
    got = _nums(line[2:])
    assert len(got) == 24
    # gs: {nus, emitx, emity, emitxp, emityp, conv} × 4 (emity は 1e-18 の雑音)
    want = [(.005, 3.1670642157266124e-6, 3.1666436597976484e-6, 6.842685666702896e-10),
            (.01, 3.1697582462604716e-6, 3.1694024756716114e-6, 5.434330033433681e-10),
            (.015, 3.174348521962957e-6, 3.174101140498108e-6, 4.221629816394806e-10),
            (.02, 3.179470474788337e-6, 3.179375065752695e-6, 3.967387730426545e-10)]
    for k, (nus, ex, exp_, conv) in enumerate(want):
        row = got[6 * k:6 * k + 6]
        assert abs(row[0] - nus) < 1e-12
        assert abs(row[1] - ex) < 3e-5 * ex
        assert abs(row[3] - exp_) < 3e-5 * exp_
        assert abs(row[5] - conv) < 3e-3 * conv


@needs_pk
def test_synchrobeta_default_azimuthal_modes():
    out = _run('r=SynchroBetaEmittance[0.01];\nPrint["R ",r];\n')
    row = _nums([l for l in out.splitlines() if l.startswith('R ')][0][2:])
    assert abs(row[1] - 2.663127335876036e-6) < 3e-5 * 2.66e-6
    assert abs(row[3] - 2.6626116449334464e-6) < 3e-5 * 2.66e-6
