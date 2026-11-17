"""
ffscmd.py — SAD の FFS コマンド層 (tffsa.f の直訳的な再現) と FFS 系組み込み関数

SAD の FFS では各文の先頭語が
  tflag (RAD/NORAD/RFSW/COD/...) → tgetfv (FIT 条件: NX 7.5 2 など) → terror (DELX/DELY/DK/DTHETA/...)
  → 要素名・パターン (tfgetv: "QF K1 0.5", "SDA MAX 0") → コマンド (USE/CELL/INS/CALC/GO/FIT/FREE/FIX/FITP/
  EMIT/DISP/TYPE/SAVE/RESET/COUPLE/MAXITERATION/STOP/...)
の順に解釈される。ここではそれを `FFS.execute(statement)` で行う。
SADScript からは Twiss[], LINE[], FFS[], Element[], Emittance[], TrackParticles[], SeedRandom[] などを
`FFS.install(interp)` で sadscript.Interp に登録する。

光学計算は pysad.lattice (CELL: twiss_table, INS: ins_twiss)、マッチングは pysad.match6 (tffsmatch の直訳)、
エミッタンスは pysad.emit6d、追跡は pysad.track6 を使う。
"""
from __future__ import annotations

import fnmatch
import functools
import math
import re
import sys
from typing import Any, Dict, List, Optional

import numpy as np

from . import lattice as _lat
from .lattice import TwissRow
from .match6 import Matcher, _FUNS, _sad_get, _sad_set, _PER_LENGTH, PI2
from .sadrandom import SadRandom
from .sadscript import _is_pattern, Interp, Sym, Rule, Null, Expr, SADError, S, fmt_real, is_num, _Lit, truthy
from . import compat as _compat



def _compat_v2():
    from . import compat
    return compat.v2()

class _DefaultMomentum(float):
    """MOMENTUM が未設定で、K64 の既定 (30 GeV) のままであることの印 (v2 の警告用。B67)"""

# SAD の要素型コード (idtype)
# Twiss["*"] の並び (tffs.f tffitcode: mfitax=1 .. mfitzpy=ntwissfun=28)
_TWISS_ALL = ['AX', 'BX', 'NX', 'AY', 'BY', 'NY', 'EX', 'EPX', 'EY', 'EPY', 'R1', 'R2', 'R3', 'R4', 'DETR',
              'DX', 'DPX', 'DY', 'DPY', 'DZ', 'DDP', 'AZ', 'BZ', 'NZ', 'ZX', 'ZPX', 'ZY', 'ZPY']
# 要素種別番号 (tfstk.f: icSOL=20, icMULT=22)
_TYPECODE = {'Drift': 1, 'Bend': 2, 'Quadrupole': 4, 'Sextupole': 6, 'Octupole': 8, 'Decapole': 10,
             'Dodecapole': 12, 'Solenoid': 20, 'Multipole': 22, 'Cavity': 31, 'TCavity': 32, 'Monitor': 42, 'Marker': 41}
_SADTYPE2CLS = {'BEND': 'Bend', 'QUAD': 'Quadrupole', 'SEXT': 'Sextupole', 'OCT': 'Octupole', 'MULT': 'Multipole',
                'CAVI': 'Cavity', 'MONI': 'Monitor', 'MARK': 'Marker', 'DRIFT': 'Drift', 'SOL': 'Solenoid'}


def _cls(el) -> str:
    """要素の SAD 上の型 (parser がドリフトに変換した ANGLE=0 BEND などは元の型)"""
    t = getattr(el, '_sad_type', None)
    return _SADTYPE2CLS.get(t, type(el).__name__) if t else type(el).__name__
_TYPENAME = {'Drift': 'DRIFT', 'Bend': 'BEND', 'Quadrupole': 'QUAD', 'Sextupole': 'SEXT', 'Octupole': 'OCT',
             'Decapole': 'DECA', 'Dodecapole': 'DODECA', 'Multipole': 'MULT', 'Solenoid': 'SOL', 'Cavity': 'CAVI',
             'TCavity': 'TCAVI', 'Monitor': 'MONI', 'Marker': 'MARK'}
# 要素型ごとのキーワード (SAD の順: 2 番目が既定の変数)
_KEYS = {'Drift': ['L'], 'Bend': ['L', 'ANGLE', 'K1', 'E1', 'E2', 'ROTATE', 'DROTATE', 'K0', 'DX', 'DY', 'F1', 'FRINGE', 'FB1', 'FB2'],
         'Quadrupole': ['L', 'K1', 'DX', 'DY', 'ROTATE', 'F1', 'F2', 'FRINGE', 'K0', 'DISFRIN'],
         'Sextupole': ['L', 'K2', 'DX', 'DY', 'ROTATE', 'K0'], 'Octupole': ['L', 'K3', 'DX', 'DY', 'ROTATE'],
         'Decapole': ['L', 'K4', 'DX', 'DY', 'ROTATE'], 'Dodecapole': ['L', 'K5', 'DX', 'DY', 'ROTATE'],
         'Multipole': ['L', 'K1', 'K0', 'K2', 'K3', 'DX', 'DY', 'ROTATE', 'ANGLE'], 'Solenoid': ['L', 'BZ', 'DX', 'DY', 'BOUND'],
         'Cavity': ['L', 'VOLT', 'FREQ', 'HARM', 'PHI', 'DX', 'DY'], 'Monitor': ['L', 'DX', 'DY'],
         'Marker': ['L', 'AX', 'BX', 'NX', 'AY', 'BY', 'NY', 'EX', 'EPX', 'EY', 'EPY', 'DX', 'DPX', 'DY', 'DPY', 'DZ', 'DDP',
                    'EMITX', 'EMITY', 'EMITZ', 'DP', 'SIGZ', 'SIGE', 'OFFSET']}
_ATTR = {'ROTATE': 'ROTATION', 'DROTATE': 'DROTATE'}       # SAD キーワード → pysad 属性 (BEND の ROTATE だけは Bend.ROTATE: _attr_of)


class _AttrMap:
    """_amap(el).get(key, key) の代わり。BEND の ROT (ROTATE) は設計ロールの属性 Bend.ROTATE に入れる (他の型は ROTATION)。
    DROTATE は DROT の属性 (BEND/MULT)。"""
    def __init__(self, el):
        self.bend = type(el).__name__ == 'Bend'

    def get(self, key, default=None):
        if key == 'ROTATE' and self.bend:
            return 'ROTATE'
        return _ATTR.get(key, default)

    def __contains__(self, key):
        return key in _ATTR


def _amap(el):
    return _AttrMap(el)
# FFS のフラグ (tffs.f fname/sino)。値は initb1.f の既定
# RFSW の既定は True (track.f:146-148: FFS 進入時に trpt=.false., rfsw=.true. に強制)
_FLAGS_DEFAULT = {'TRPT': False, 'CELL': False, 'RFSW': True, 'RAD': False, 'FLUC': True, 'RADCOD': False, 'RADTAPER': False,
                  'CALC6D': False, 'KEEPEXP': True, 'CALEXP': True, 'EMIOUT': False, 'CODPLOT': False, 'COD': True, 'INTRA': False,
                  'POL': False, 'RADPOL': True, 'LWAKE': False, 'TWAKE': False, 'WSPAC': False, 'SPAC': False, 'SELFCOD': False,
                  'CONV': False, 'STABLE': True, 'GEOCAL': True, 'RADLIGHT': False, 'PHOTONS': False, 'LOSSMAP': False, 'SORG': False,
                  'INTRES': True, 'HALFRES': True, 'SUMRES': True, 'DIFFRES': False, 'FFSPRMPT': False, 'CONVCASE': True,
                  'PRSVCASE': False, 'SUS': True, 'K64': True, 'GAUSS': False, 'FIXSEED': True, 'BIPOL': True, 'PSPAC': False,
                  'ORBITCAL': True, 'CALOPT': True, 'DAPERT': False, 'IDEAL': False, 'FOURIER': False, 'TRACKSIZ': False,
                  'SIMULATE': True, 'ABSW': True, 'JITTER': True, 'TRGAUSS': True, 'BARYCOD': True, 'WAKEOPT': False,
                  'LOG': False, 'ECHO': False, 'DEBUG': False, 'CTIME': False, 'SMEAR': True}
# sino (反対語): word → (flag, 値)。NO+sino は逆
_SINO = {'RING': 'TRPT', 'INS': 'CELL', 'DAMPONLY': 'FLUC', 'CALC4D': 'CALC6D', 'UNSTABLE': 'STABLE', 'GEOFIX': 'GEOCAL',
         'LEGACY': 'K64', 'UNIFORM': 'GAUSS', 'MOVESEED': 'FIXSEED', 'UNIPOL': 'BIPOL', 'ORBONLY': 'CALOPT', 'REAL': 'IDEAL',
         'OPERATE': 'SIMULATE', 'RELW': 'ABSW', 'QUIET': 'JITTER', 'TRUNI': 'TRGAUSS'}
_FLAG_ALT = {k: (v, False) for k, v in _SINO.items()}
_FLAG_ALT.update({'NO' + k: (v, True) for k, v in _SINO.items()})
_ELEM_KEYWORDS = {'L', 'K0', 'K1', 'K2', 'K3', 'K4', 'K5', 'ANGLE', 'E1', 'E2', 'ROTATE', 'DROTATE', 'DX', 'DY', 'DZ', 'VOLT', 'FREQ',
                  'HARM', 'PHI', 'BZ', 'F1', 'F2', 'FRINGE', 'DISFRIN', 'SK1', 'SK2', 'AX', 'BX', 'AY', 'BY', 'EX', 'EPX', 'EY', 'EPY'}
_COMMANDS = {'USE', 'VISIT', 'BYE', 'STOP', 'EXIT', 'QUIT', 'END', 'ABORT', 'SAVE', 'RESET', 'FREE', 'FIX', 'FIT', 'FITP', 'FITPOINTS',
             'CAL', 'CALC', 'CALCULATE', 'GO', 'TYPE', 'DISP', 'DISPLAY', 'EMIT', 'MAXI', 'MAXITERATION', 'COUP', 'COUPLE',
             'DECOUP', 'DECOUPLE', 'INDEP', 'INDEPENDENT', 'DUMP', 'VARY', 'VARS', 'VARIABLES', 'REJ', 'REJECT', 'ORG', 'ORIGIN',
             'STAT', 'STATUS', 'SHOW', 'DRAW', 'GEO', 'MCOD', 'TCOD', 'MON', 'REC', 'RECOVER', 'INC', 'INCLUDE', 'REPEAT',
             'UNTIL', 'IF', 'ELSE', 'ENDIF', 'CHRO', 'CHROMATICITY', 'RADINT', 'BEAM', 'BEAMSIZE', 'EXPAND', 'SCALE',
             'ALIGN', 'ATT', 'ATTRIBUTE', 'GCUT', 'DP', 'DP0', 'MATRIX', 'MAT', 'RENUM', 'RENUMB', 'RENUMBE', 'RENUMBER', 'REV', 'REVERSE', 'ZAP', 'DIMAD',
             'SYNCHROB', 'INTRA', 'CODPLOT', 'NOCODPLOT', 'ON', 'OFF', 'REF', 'REFE', 'REFER', 'REFERE', 'REFEREN', 'REFERENC', 'REFERENCE', 'PRINT', 'MEASURE', 'TRC', 'TRCOD', 'DELC', 'DELCO', 'DELCOR',
             'ATTR', 'VAR', 'EXEC', 'EXECUTE', 'PRI', 'SPLIT', 'OUTPUT', 'OUT', 'OUTP', 'APPEND', 'APP', 'CLOSE', 'CLO', 'TERMINATE', 'TERM',
             'INPUT', 'IN', 'READ', 'SUSPEND', 'SUSP', 'RESUME', 'RES'}
_ERRCMDS = {'DELX', 'DELY', 'DTHETA', 'DL', 'DK', 'DDK', 'DBEAM', 'DUMP'}


class _ZeroRow:
    """REFERENCE 未設定のときの参照光学 (twiss(:,-1,:) は 0)"""
    def __getattr__(self, name):
        return 0.0


def _shortest(v: float) -> str:
    """SAD の TYPE の値の表現: 最短の往復可能表現で、先頭の 0 は書かない (.413, -.6422919019387, 1.0604857646104)"""
    if v == int(v) and abs(v) < 1e15:
        return str(int(v))
    r = repr(float(v))
    if 'e' in r:
        return r
    if r.startswith('0.'):
        return r[1:]
    if r.startswith('-0.'):
        return '-' + r[2:]
    return r


class _NullOut:
    def write(self, s):
        pass

    def flush(self):
        pass


class _RecOut:
    """FFS[] の出力を転送しながら記録する (エラー行の検出用)"""
    def __init__(self, target):
        self.target = target
        self.text = ''

    def write(self, s):
        self.text += s
        self.target.write(s)

    def flush(self):
        if hasattr(self.target, 'flush'):
            self.target.flush()


def _env_get(env, name):
    try:
        return env.get(name)
    except Exception:
        return None


def _fr(v, form):
    """EMIT 報告の autofg(v, 'w.d'): fmt_real と同じだが NaN は幅 w に右詰め (SAD の報告は桁をそろえる)"""
    v = float(v)
    return 'NaN'.rjust(int(form.split('.')[0])) if v != v else fmt_real(v, form)


def _ff(v, w, d):
    """Fortran の fw.d (NaN は 'NaN')"""
    return 'NaN'.rjust(w) if v != v else f'{v:{w}.{d}f}'


def _unit_ext_twiss():
    """ri = 単位行列の tfetwiss (28 個)。閉軌道が見つからなかった EMIT の報告用"""
    from .calc6d import tfetwiss
    return np.asarray(tfetwiss(np.eye(6), np.zeros(6), True), dtype=float)


def _f107(v):
    """Fortran の f10.7 (NaN は 'NaN'。gfortran と同じ)"""
    return _ff(v, 10, 7)


def _line9013(vals, labs=('X :', 'Y :', 'Z :')):
    """temitprint 9013: 10x, 3(a, 1p g14.6, 3x) の 'X :' 'Y :' 'Z :' 行 (行末の 3x は出力されない)"""
    from .sadscript import fmt_g

    def g(v):
        v = float(v)
        return 'NaN'.rjust(14) if v != v else fmt_g(v, 14, 6)
    return '          ' + '   '.join(f'{lab}{g(v)}' for lab, v in zip(labs, vals)) + '\n'


def _gwd(v, w, d):
    """Fortran (gfortran) の 1p gw.d: 0.1<=|x|<10^d は小数 d-k 桁 (k = 整数部の桁数) の F + 4 空白、0 は '0.' + d-1 桁の 0 + 4 空白、それ以外は E 形式"""
    v = float(v)
    if v != v:
        return 'NaN'.rjust(w)
    a = abs(v)
    if a == 0.0:
        s = '0.' + '0' * (d - 1)
    elif 0.1 <= a < 10.0**d:
        k = len(str(int(a))) if a >= 1.0 else 0
        s = f'{v:.{d - k}f}'
    else:
        return f'{v:.{d}E}'.rjust(w)                  # 1p の E は小数 d 桁 (gs: 4.22829E-04)
    return (s + '    ').rjust(w)


def _g157(v):
    """Fortran の 1p g15.7 (gfortran): 0.1<=|x|<1e7 は小数 7-k 桁 (k = 整数部の桁数) の F + 4 空白、それ以外は 1.2345678E-02"""
    v = float(v)
    if v != v:
        return 'NaN'.rjust(15)
    a = abs(v)
    if a == 0.0:
        s = '0.000000'
    elif 0.1 <= a < 1e7:
        k = len(str(int(a))) if a >= 1.0 else 0
        s = f'{v:.{7 - k}f}'
    else:
        return f'{v:.7E}'.rjust(15)
    return (s + '    ').rjust(15)


def _abbrev(word: str, full: str) -> bool:
    """SAD の abbrev: 'CAL_CULATE' は CAL, CALC, ... CALCULATE を許す"""
    head, _, tail = full.partition('_')
    if not word.startswith(head):
        return False
    rest = word[len(head):]
    return (head + tail).startswith(word) and len(rest) <= len(tail)


def _squeeze_space(text: str, n: int):
    """tfdisp.f tfsqueezespace: 直前も空白の空白を先頭から n 個まで除く。戻り (除いて末尾を空白で埋めた文字列, 除いた個数)"""
    if n <= 0 or not text:
        return text, 0
    out = [text[0]]
    nsq = 0
    for i in range(1, len(text)):
        if text[i - 1] == ' ' and text[i] == ' ' and nsq < n:
            nsq += 1
            continue
        out.append(text[i])
    return ''.join(out) + ' ' * nsq, nsq


@functools.lru_cache(maxsize=4096)
def _wild_rx(p: str):
    """_wild_match の正規表現 (p は大文字化済み。コンパイル結果をキャッシュ)"""
    rx = ''
    i = 0
    while i < len(p):
        c = p[i]
        if c == '*':
            rx += '.*'
        elif c == '%':
            rx += '.'
        elif c == '{':
            j = p.find('}', i)
            if j < 0:
                rx += re.escape(c)
            else:
                body = p[i + 1:j]
                if body.startswith('^'):
                    rx += '[^' + re.escape(body[1:]) + ']'
                else:
                    rx += '[' + re.escape(body) + ']'
                i = j
        elif c == '|':
            rx += '|'
        else:
            rx += re.escape(c)
        i += 1
    return re.compile(rx)


def _wild_match(name: str, pat: str) -> bool:
    """SAD の tmatch (src/tmatch.f): * (任意)、% (1 文字)、{abc} (文字集合)、{^abc} (否定)、| (または)。大文字で比較。
    ? は SAD では普通の文字 (ワイルドカードではない: Element["NAME","Q?"] は一致しない — probe13 で確認)"""
    return _wild_rx(pat.upper()).fullmatch(name.upper()) is not None


def _mark_offset_value(el) -> float:
    """MARK の OFFSET キーワード (属性か _raw)"""
    v = getattr(el, 'OFFSET', None)
    if v is None:
        v = (getattr(el, '_raw', None) or {}).get('OFFSET', 0.0)
    try:
        return float(v or 0.0)
    except (TypeError, ValueError):
        return 0.0


class _RadLight:
    """RADLIGHT (tlinit.f tlinit / tlstore / tlresult): 粒子の軌道 {t, x, y, z} (GEO 座標) を要素ごとに記録する。
    k64: 記録するのは区間の始めと、L≠0 の BEND 以外 (型 < CAVI) の要素の後だけ (tturn.f:555。BEND と曲げのある MULT の中は記録しない: K64 のまま)。
         横方向の軸は 0 (tlinit が geo0 の軸 1・2 を書かない) で、BEND で向きも変わらない (docs/21 B65)。
    v2 (B65 の修正): 横方向の軸を幾何から取り、BEND では tlstore の dphi≠0 の式 (円弧) で位置と向きを進める。BEND の中は
         1/γ の角度あたり RADLIGHTDIV 点 (既定 16) に刻んで記録する (放射場の時間構造は 1/γ の角度で決まる)。
    geo = ls の位置の幾何 (3×4: 列 0..2 が x, y, z の軸、列 3 が位置)"""

    def __init__(self, geo, ph, npart, v2=False, div=16.0):
        g = np.array(geo, dtype=float)
        self.geo0 = g.copy()                 # tlinit: geo0(:,3)=geo(:,3)、位置 = geo(:,4)。軸 1・2 は書かない (common block の初期値 0 のまま: K64 の癖。
        self.v2 = v2
        if not v2:
            self.geo0[:, :2] = 0.0           # gs で確かめた: 横方向の座標は軌道に入らず、進行方向 (geo の z 軸) に沿った点だけが記録される)
        self.t0 = 0.0
        self.v0 = ph.p0 / ph.h0 * 299792458.0
        self.gamma = ph.h0
        self.div = div
        self.npart = npart
        self.pts = [[[], [], [], []] for _ in range(npart)]
        self._bprog = {}

    def store(self, tr, dl, dphi=0.0, rho0=0.0, theta=0.0, S=None, n=None, kid=None):
        """tlstore (keep=True)。dl は要素の長さ (開始点は 0)。dphi≠0 (v2 だけ) は偏向角 dphi・曲率半径 rho0・回転 theta の円弧"""
        g0 = self.geo0
        v0 = self.v0
        if dphi == 0.0:
            gpos = g0[:, 3] + dl * g0[:, 2]
            gt = self.t0 + dl / v0
            geoi = g0[:, :3].copy()
        else:
            c, s_ = math.cos(theta), math.sin(theta)
            ex = c * g0[:, 0] - s_ * g0[:, 1]            # 偏向面の軸 (ROTATE だけ回した x)
            ey = s_ * g0[:, 0] + c * g0[:, 1]
            ez = g0[:, 2]
            sp0, cp0 = math.sin(dphi), math.cos(dphi)
            r1 = rho0 * sp0
            r2 = rho0 * sp0**2 / (1.0 + cp0) if cp0 >= 0.0 else rho0 * (1.0 - cp0)
            gpos = g0[:, 3] + (r1 * ez - r2 * ex)
            ex1 = cp0 * ex + sp0 * ez
            ez1 = -sp0 * ex + cp0 * ez
            geoi = np.empty((3, 3))
            geoi[:, 0] = c * ex1 + s_ * ey               # 回転を戻した出口の軸 (粒子の x, y はこの軸の座標)
            geoi[:, 1] = -s_ * ex1 + c * ey
            geoi[:, 2] = ez1
            gt = self.t0 + abs(rho0 * dphi) / v0
        if S is None:
            S, n, kid = tr.S, (tr.n if tr.S is not None else 0), tr.kid
        dvfs = tr.ph.dvfs
        for s in range(n):
            i = int(kid[s])
            x, y, z, dv = S[s, 0], S[s, 2], S[s, 4], S[s, 6]
            self.pts[i][0].append(gt - z / (v0 * (1.0 - dv + dvfs)))
            self.pts[i][1].append(float(gpos[0] + x * geoi[0, 0] + y * geoi[0, 1]))
            self.pts[i][2].append(float(gpos[1] + x * geoi[1, 0] + y * geoi[1, 1]))
            self.pts[i][3].append(float(gpos[2] + x * geoi[2, 0] + y * geoi[2, 1]))
        self.geo0[:, :3] = geoi
        self.geo0[:, 3] = gpos
        self.t0 = gt

    def bend(self, tr, el, l):
        """v2: BEND (ANGLE≠0、L≠0) を n 個に刻み、刻みの終わりごとに記録する。粒子の本当の追跡は要素全体を 1 回 (tr.segment) で、
        中の点は状態の写しを刻んだ BEND (端の角は最初と最後だけ) で進めて作る (放射の揺らぎなし)。
        刻みと全体の一致は ずれ・F1・FRINGE がなければ 1e-11 以下、あれば 1e-6 m 程度 (中の点だけ)"""
        import copy as _copy
        from .track6n import Prog, turn
        ang = float(el.ANGLE)
        n = max(1, int(math.ceil(abs(ang) * self.gamma * self.div)))
        rho0 = float(el.L) / ang
        theta = float(getattr(el, 'ROTATE', 0.0) or 0.0)
        if n > 1 and tr.n > 0:
            key = (id(el), n)
            progs = self._bprog.get(key)
            if progs is None:
                progs = []
                for k in (0, 1):                 # 最初の刻み (入口の端)、中の刻み
                    e = _copy.copy(el)
                    e.L = el.L / n
                    e.ANGLE = ang / n
                    e.E1 = el.E1 if k == 0 else 0.0
                    e.E2 = 0.0
                    if hasattr(e, 'K0'):
                        e.K0 = float(el.K0 or 0.0) / n
                    progs.append(Prog([e], tr.ph))
                self._bprog[key] = progs
            S = tr.S.copy()
            n0 = tr.n
            kid = tr.kid.copy()
            rs = tr.rs.copy()
            with np.errstate(all='ignore'):
                for k in range(n - 1):
                    turn(progs[0 if k == 0 else 1], S, n0, tr.ph, tr.trf0, tr.rfsw, tr.rad, False, rs, trpt=tr.trpt)
                    self.store(tr, 0.0, ang / n, rho0, theta, S=S, n=n0, kid=kid)
        tr.segment(l, l + 1)
        self.store(tr, 0.0, ang / n, rho0, theta)

    def result(self):
        return [[[float(v) for v in row] for row in pt] for pt in self.pts]


class FFS:
    def __init__(self, script_env, interp: Interp, out=None):
        self.env = script_env           # parser.SADEnvironment (要素定義・ビームライン)
        self.I = interp
        self.out = out or sys.stdout
        self.flags = dict(_FLAGS_DEFAULT)
        # MAIN 部の `ON/OFF RFSW` は MAIN の整数フラグ $RFSW$ (track.f:80 isynch) で、MAIN の TRACK だけが使う。
        # FFS に入るとき rfsw=.true. に固定される (track.f:149) ので FFS の RFSW には届かない (<work-item>/mflag*.sad)
        self.main_flags = {'RFSW': False}
        self.lat = None
        self._spread_pending = set()     # LINE[key,"NAME.n"]= で設定した (ファミリー, key): 次の CALC で代表の値に揃える
        self.name = ''
        self.els: List = []
        self.nlat = 1
        self._last_rows = None
        self.rows: Optional[List[TwissRow]] = None
        self._beam_cache = None            # LINE["SIG.."] 用の平衡ビーム行列 (遅延計算)
        self.geo = None                    # geo.GeoResult (tfgeo の結果)
        self.geo_dirty = True              # CALC / EMIT / USE / ORG で再計算が必要 (GEOCAL のとき)
        self.org_k, self.org_g, self.org_chi = 1, [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]   # ORG コマンドの状態 (chi は rad)
        self._geo_warned = set()
        self._visit_stack: List[str] = []   # VISIT した LINE の履歴 (BYE で戻る)
        self._visit_syms: List[dict] = []   # VISIT で退避した FitFunction・FitValue・ElementValues・OpticsProlog・OpticsEpilog の値 (BYE で戻す)
        self._newline_count = 0            # USE BeamLine[...] が作る LINE の通し番号 (L000001$ ...)
        self.touschek = None               # INTRA ON の emit の結果 (touschek.touschek_table の戻り)
        self.calc6d_result = None          # 直近の CALC が CALC6D だったときの calc6d.Calc6DResult
        self.matcher: Optional[Matcher] = None
        self.couples: Dict[str, tuple] = {}
        self.emit: Optional[dict] = None
        self.track_results: List[Optional[dict]] = []    # MAIN の TRACK ごとの結果 (maintrack.run_track の戻り。検証用)
        self.trf0 = 0.0
        from . import track6 as _t6
        self.rng = _t6.RNG               # 追跡 (track6/track6n/dapert) と同じ乱数 (SAD は 1 つの乱数列)
        self.rng.restore(17)
        interp.rng = self.rng
        self.stopped = False
        self._stop_count = 0
        self.last_match_print = ''
        g = interp.globals
        g.setdefault('DP', 0.01)
        g.setdefault('DP0', 0.0)
        g.setdefault('NPARA', 1.0)
        g.setdefault('DAPWIDTH', 51.0)
        g.setdefault('TURNS', -1.0)          # gs: FFS に入った直後の TURNS は -1
        g.setdefault('NP', 1000.0)
        g.setdefault('GCUT', 1e300)
        g.setdefault('FSHIFT', 0.0)
        g.setdefault('MINCOUP', 1e-2)
        g.setdefault('XANGLE', 0.0)
        g.setdefault('LOSSAMPL', 1.0)
        for k, v in (('OMEGA0', 0.0), ('DTSYNCH', 0.0), ('PHICAV', 0.0), ('EFFVCRATIO', 1.0), ('EFFVC', 1.0), ('EFFRFFREQ', 0.0),
                     ('EMITXE', 0.0), ('EMITYE', 0.0), ('EMITZE', 0.0), ('EMITXC', 0.0), ('EMITYC', 0.0), ('EMITZC', 0.0),
                     ('PBUNCH', 1e10), ('NBUNCH', 1.0)):
            g.setdefault(k, v)                         # initb1.f の既定値
        g.setdefault('LOSSDZ', 100.0)
        g.setdefault('CONVERGENCE', 1e-9)
        g.setdefault('MAXITERATION', 40.0)
        g.setdefault('ExponentOfResidual', 2.0)
        g.setdefault('OffMomentumWeight', 1.0)
        g.setdefault('MatchingResidual', 0.0)
        g.setdefault('NetResidual', 0.0)
        g.setdefault('UnstableLevel', 0.0)
        for k in ('EMITX', 'EMITY', 'EMITZ', 'SIGE', 'SIGZ', 'XIX', 'XIY', 'DPM'):
            g.setdefault(k, 0.0)
        g.setdefault('CHARGE', 1.0)
        g.setdefault('MASS', 0.51099895e6)          # SAD K64 の elmass (tfstk.f, CODATA 2018)
        if 'MOMENTUM' not in g:
            m = _env_get(script_env, 'MOMENTUM')
            g['MOMENTUM'] = float(m) if isinstance(m, (int, float)) else 30e9      # K64 initb1.f: RsetGL('MOMENTUM',30.d9)
            if not isinstance(m, (int, float)) and _compat.v2():
                g['MOMENTUM'] = _DefaultMomentum(30e9)      # v2 (B67): 既定のままかを見分ける印 (値は同じ 30e9。代入されると普通の実数になる)
        self._env_mom_seen = _env_get(script_env, 'MOMENTUM')
        self.install(interp)

    # ------------------------------------------------------------------ 組み込み関数の登録
    def install(self, I: Interp):
        I.ffs_handler = self.execute
        I.ffs_words = self.is_word
        B = I.builtins
        B['Twiss'] = self.fn_twiss
        B['LINE'] = self.fn_line
        B['Master'] = self.fn_master

        B['RenumberElement'] = self.fn_renumber_element
        B['ResponseMatrix'] = self.fn_response_matrix
        from . import corr as _corr
        _corr.install(I)
        from . import pkgload as _pkgload
        _pkgload.register(I)
        from .optics_fn import OpticsFunctions
        self._optfn = OpticsFunctions(self)
        self._optfn.install(B)                 # TransferMatrix, BeamMatrix, OrbitGeo, GeoBase, SymplecticJ, SurvivedParticles, AccelerateParticles
        B['FFS$SHOW'] = lambda a: self.ffs_show_list()
        B['FFS$FlagQ'] = lambda a: self.flag_query(a[0] if a else '')          # ?FLAG (SADScript の構文、sadscript.parse_prefix)
        def _rdims(a, name, draw):
            # Random[n]・Random[n, m, …]: n×m×… の配列 (小数は切り捨て、0 は {}、負は wrongnum)。リストの引数は wrongtype (gs)
            if not all(is_num(x) for x in a):
                raise SADError('General::wrongtype: Argument must be Real numbers in %s[%s]' % (name, ','.join(self.I.to_str(x, None, True) for x in a)))
            dims = [int(x) for x in a]
            if any(d < 0 for d in dims):
                raise SADError('General::wrongnum:  Number must be positive or 0 in %s[%s]' % (name, ','.join(self.I.to_str(x) for x in a)))
            tot = 1
            for d in dims:
                tot *= d
            flat = draw(tot) if tot > 0 else []
            flat = flat if isinstance(flat, list) else [flat]

            def build(k, off):
                if k == len(dims):
                    return flat[off]
                stride = 1
                for d in dims[k + 1:]:
                    stride *= d
                return [build(k + 1, off + t * stride) for t in range(dims[k])]
            return build(0, 0) if tot > 0 or len(dims) == 1 else build(0, 0) if dims[0] > 0 else []

        def _parabolarandom(a):
            if not a:
                return self.rng.parabola(1)
            return _rdims(a, 'ParabolaRandom', lambda n: self.rng.parabola(n))
        B['ParabolaRandom'] = _parabolarandom
        B['Element'] = self.fn_element
        B['SynchroBetaEmittance1'] = self.fn_synchrobeta
        B['FFS'] = self.fn_ffs
        B['Emittance'] = self.fn_emittance
        B['TrackParticles'] = self.fn_track
        B['SeedRandom'] = self.fn_seedrandom
        def _gaussrandom(a):
            if not a:
                return self.rng.gauss(1, self.gcut())
            return _rdims(a, 'GaussRandom', lambda n: self.rng.gauss(n, self.gcut()))
        B['GaussRandom'] = _gaussrandom

        def _random(a):
            if not a:
                return self.rng.uniform0()
            return _rdims(a, 'Random', lambda n: [self.rng.uniform0() for _ in range(n)])
        B['Random'] = _random
        B['BeamLineName'] = lambda a: self.name
        B['CalculateOptics'] = self.fn_calculate_optics
        B['DynamicApertureSurvey'] = self.fn_dasurvey
        B['TouschekLifetime'] = self.fn_touschek
        B['FFSCalc'] = lambda a: (self.calc(), Null)[1]
        I.set_hooks = getattr(I, 'set_hooks', {})
        I.set_call_hooks = {'LINE': self.set_line, 'Element': self.set_element}
        # Set[LINE[...], v] を扱えるよう Set を包む
        orig_set = B['Set']

        def _set(args):
            lhs = args[0]
            if isinstance(lhs, Expr) and isinstance(lhs.head, Sym) and lhs.head.name in I.set_call_hooks:
                v = I.eval(args[1])
                r = I.set_call_hooks[lhs.head.name](I.eval_args(lhs.args), v)
                return r if (lhs.head.name == 'Element' and r is not None) else v
            return orig_set(args)
        B['Set'] = _set
        # Element/SetElement/FLAG/ExtractBeamLine/BeamLine/PrintBeamLine/WriteBeamLine/GetMAIN (sadelem.py)
        from . import sadelem
        sadelem.install(self)

    # ------------------------------------------------------------------ ユーティリティ
    def write(self, s: str):
        self.out.write(s)
        self.out.flush()

    def gcut(self):
        g = self.I.get('GCUT', 1e300)
        return float(g) if is_num(g) and g > 0 else 1e300

    def dp(self):
        return float(self.I.get('DP', 0.01))

    def dp0(self):
        return float(self.I.get('DP0', 0.0))

    def momentum(self, warn: bool = True):
        """設計運動量 [eV/c]。v2 (B67): MAIN でも FFS でも MOMENTUM が設定されないまま、運動量を使う計算
        (EMIT・Emittance・追跡・RADINT・CALC6D・DA ほか) に入ったら、1 度だけ警告する (K64 は黙って 30 GeV)"""
        m = self.I.get('MOMENTUM', 30e9)
        if warn and type(m) is _DefaultMomentum and not getattr(self, '_mom_warned', False) and _compat.v2():
            self._mom_warned = True
            self.write(' ???-FFS-Warning-MOMENTUM is not set: using the default 30 GeV (K64 initb1.f). '
                       'Set it as "MOMENTUM = 3 GEV;" in MAIN or "MOMENTUM = 3e9;" in FFS. ' + chr(10))
        return float(m)

    def is_word(self, word: str) -> bool:
        w = word.upper()
        if w in _COMMANDS or w in _ERRCMDS or w in self.flags or (w.startswith('NO') and w[2:] in self.flags) or w in _FLAG_ALT:
            return True
        if _compat.v2() and (_abbrev(w, 'MEA_SURE') or _abbrev(w, 'TRC_OD') or _abbrev(w, 'LOGM_EASURE')
                             or _abbrev(w, 'ORBITJ_ITTER') or _abbrev(w, 'PUTM_EASURE') or _abbrev(w, 'GET_MEASURE')):
            return True                                 # v2: MEASURE / TRCOD / LOGMEASURE / ORBITJITTER / PUTMEASURE / GETMEASURE の省略形 (measure.py)
        if w == 'WAKE' and _compat.v2():
            return True                                 # v2: 旧式の WAKE コマンド (_wake_cmd_v2)
        if w in _FUNS or (w[:-1] in _FUNS and w[-1:] in ('M', 'I')) or w.endswith('SCALE') and w[:-5] in _FUNS:
            return True
        if any(c in w for c in '*%{|?'):
            return bool(self.els)
        return self._family(w) is not None or self._position(w) is not None

    def _nidx(self):
        """要素名の索引 (USE ごとに作り直す。要素列の同一性と長さで無効化)。性能: positions/_family/_position/elname を O(1) に"""
        key = (id(self.els), len(self.els))
        c = getattr(self, '_nidx_cache', None)
        if c is not None and c['key'] == key:
            return c
        names = [el.name.upper() for el in self.els]
        fam: Dict[str, List[int]] = {}
        ords = []
        for i, n in enumerate(names):
            lst = fam.setdefault(n, [])
            lst.append(i + 1)
            ords.append(len(lst))
        c = {'key': key, 'names': names, 'fam': fam, 'ords': ords, 'pat': {}}
        # tfinimult: BEND・QUAD・薄肉 (SEXT…) の INDEX 属性が 2 以上なら、同名の要素の番号はその値から始まる (QF の INDEX=5 → QF.5 QF.6 …)
        renum = {}
        for n, lst in fam.items():
            base = self._index_base(self.els[lst[0] - 1]) if len(lst) > 1 else 1
            if base > 1:
                for j, pos in enumerate(lst):
                    ords[pos - 1] = base + j
                renum[n] = {base + j: pos for j, pos in enumerate(lst)}
        if renum:
            c['renum'] = renum
        self._nidx_cache = c
        return c

    def _index_base(self, el) -> int:
        """要素の INDEX 属性 (K64 の kwINDX: BEND・QUAD・SEXT・OCT・DECA・DODECA) による番号の初期値。無ければ 1 (tfinit.f tfinimult: max(1, int(INDEX)))"""
        if _TYPENAME.get(_cls(el), '') not in ('BEND', 'QUAD', 'SEXT', 'OCT', 'DECA', 'DODECA'):
            return 1
        raw = getattr(el, '_raw', None) or {}
        v = raw.get('INDEX')
        if v is None:
            v = (self.env._elements.get(el.name, {}) or {}).get('INDEX', 0.0)
        try:
            return max(1, int(float(v or 0.0)))
        except (TypeError, ValueError):
            return 1

    def _family(self, name):
        lst = self._nidx()['fam'].get(name.upper())
        return lst[0] - 1 if lst else None

    def _position(self, name):
        """NAME.n → 1 始まりの位置。$$$ → nlat"""
        n = name.upper()
        if n in ('$$$', '***'):
            return self.nlat
        if n == '^^^':
            return 1
        fam, _, num = n.partition('.')
        if not num.isdigit():
            return None
        k = int(num)
        ix = self._nidx()
        rn = ix.get('renum', {}).get(fam)
        if rn is not None:
            return rn.get(k)                     # RenumberElement で付け替えた番号
        lst = ix['fam'].get(fam)
        if lst and 1 <= k <= len(lst):
            return lst[k - 1]
        return None

    def fn_renumber_element(self, args):
        """RenumberElement[name, {n1, n2, ..}] (tfrenumber_element.f): 名前が name の要素の連番 (NAME.n の n) を n1 n2 .. に付け替える。
        足りない分は最大の番号から続ける。USE のやり直しで元に戻る"""
        if len(args) != 2:
            raise SADError('General::narg: Number of arguments is expected to be 2 in RenumberElement')
        name, nums = args
        if not isinstance(name, str) or any(c in name for c in '*%{<|.'):
            raise SADError('General::wrongtype: Argument must be Name of component family for #1')
        if not isinstance(nums, list):
            raise SADError('General::wrongtype: Argument must be List of positive sequence numbers for #2')
        ix = self._nidx()
        lst = ix['fam'].get(name.upper())
        if not lst:
            raise SADError('General::wrongtype: Argument must be Name of component family for #1')
        if len(lst) < 2 or not nums:
            return Null
        fwd = [int(x) for x in nums[:min(len(nums), len(lst))]]
        last = 0
        for j, v in enumerate(fwd):
            if v < 1 or v in fwd[:j] and v <= last:
                raise SADError('General::wrongtype: Argument must be List of positive sequence numbers for #2')
            last = max(last, v)
        newmap = {}
        for j, pos in enumerate(lst):
            if j < len(fwd):
                num = fwd[j]
            else:
                last += 1
                num = last
            ix['ords'][pos - 1] = num
            newmap[num] = pos
        ix.setdefault('renum', {})[name.upper()] = newmap
        ix['pat'].clear()
        return Null

    def _renumber_cmd(self, words, i):
        """RENUMBER comp (tfinit.f tffsrenumber/tfinimult): comp を先頭に、周回順に同名の要素の番号 (NAME.n の n) を 1 から付け直す。
        comp が無ければ ' ???-FFS-Error-Missing origin component for RENUM_BER '。番号が 1 台だけのものは 0 (LINE["MULT"]) で名前に番号が付かない。
        要素の INDEX 属性 (kwINDX) が 2 以上なら、番号はその値から始まる (_index_base)"""
        name = words[i + 1] if i + 1 < len(words) else None
        k0 = None
        if name is not None and not any(c in name for c in '*%{|?'):
            pos = self.positions(name)
            k0 = pos[0] if pos else None
        if k0 is None:
            self.write(' ???-FFS-Error-Missing origin component for RENUM_BER ' + chr(10))
            return 1
        ix = self._nidx()
        n = self.nlat - 1
        if k0 > n:
            k0 = 1
        count: Dict[str, int] = {}
        renum: Dict[str, Dict[int, int]] = {}
        for t in range(n):
            pos = (k0 - 1 + t) % n + 1
            nm = ix['names'][pos - 1]
            if nm not in count:                  # INDEX 属性による番号の初期値 (tfinimult。gs: INDEX=5 の QF は RENUMBER 後も QF.5 から)
                count[nm] = (self._index_base(self.els[pos - 1]) - 1) if len(ix['fam'].get(nm, ())) > 1 else 0
            count[nm] = c = count[nm] + 1
            ix['ords'][pos - 1] = c
            renum.setdefault(nm, {})[c] = pos
        ix['renum'] = renum
        ix['pat'].clear()
        return 2

    def positions(self, pat: str) -> List[int]:
        """パターン (または名前) に合う位置 (1 始まり) のリスト。tflinenamestk 相当"""
        p = pat.upper()
        if p in ('$$$', '***'):
            return [self.nlat]
        ix = self._nidx()
        if p.endswith('.*') and not any(c in p[:-2] for c in '*%{|'):
            return list(ix['fam'].get(p[:-2], ()))
        if any(c in p for c in '*%{|?'):
            hit = ix['pat'].get((p, self.nlat))
            if hit is None:
                rx = _wild_rx(p)
                names, ords = ix['names'], ix['ords']
                hit = []
                cache = {}
                for i, nm in enumerate(names):
                    m = cache.get(nm)
                    if m is None:
                        m = cache[nm] = rx.fullmatch(nm) is not None
                    if m or rx.fullmatch(f'{nm}.{ords[i]}') is not None:
                        hit.append(i + 1)
                if rx.fullmatch('$$$') is not None:
                    hit.append(self.nlat)
                ix['pat'][(p, self.nlat)] = hit
            return list(hit)
        k = self._position(p)
        if k is not None:
            return [k]
        i = self._family(p)
        return [i + 1] if i is not None else []

    def _ordinal(self, i):
        return self._nidx()['ords'][i]

    def families(self, pat: str) -> List[str]:
        """パターンに合うファミリー名 (tfgetlineps mode 1)"""
        seen = []
        for k in self.positions(pat):
            if k <= len(self.els):
                nm = self.els[k - 1].name
                if nm not in seen:
                    seen.append(nm)
        return seen

    def instances(self, fam):
        lst = self._nidx()['fam'].get(fam.upper(), ())
        els = self.els
        return [els[i - 1] for i in lst if els[i - 1].name == fam]

    def elname(self, k):
        """SAD の elname: 同名が 1 台だけなら NAME、複数なら NAME.n"""
        if k >= self.nlat:
            return '$$$'
        nm = self.els[k - 1].name
        if len(self.instances(nm)) == 1:
            return nm
        return f'{nm}.{self._ordinal(k - 1)}'

    # ------------------------------------------------------------------ 要素値 (SAD 単位)
    def get_value(self, el, key: str):
        key = key.upper()
        tn = type(el).__name__
        if getattr(el, '_sad_type', None) and not hasattr(el, _amap(el).get(key, key)):
            return float(self.env._elements.get(el.name, {}).get(key, 0.0) or 0.0)
        if tn == 'Marker':
            if key == 'DDP' and self.els and el is self.els[0]:
                return self.dp0()           # SAD track.f: DP0 は LINE["DDP",1] の別名
            raw = self.env._elements.get(el.name, {})
            v = raw.get(key)
            if v is None:
                v = getattr(el, _amap(el).get(key, key), 0.0)
            return float(v or 0.0)
        if key in ('K0', 'K1', 'K2', 'K3', 'K4', 'K5', 'ANGLE', 'L', 'E1', 'E2', 'VOLT', 'FREQ', 'HARM', 'PHI', 'BZ', 'F1', 'F2', 'FRINGE', 'DISFRIN', 'DX', 'DY', 'ROTATE', 'DROTATE', 'K1S', 'SK1'):
            if key in ('E1', 'E2') and tn == 'Bend':
                return self._e_ratio(el, key)           # SAD の E1/E2 は ANGLE に対する比 (tffs.f:1261)
            if key == 'HARM':
                return float(getattr(el, '_harm', 0) or 0)
            v = _sad_get(el, key if key not in _ATTR else _amap(el).get(key, key)) if key in ('K1', 'K2', 'K3', 'K4', 'K5') else float(getattr(el, _amap(el).get(key, key), 0.0) or 0.0)
            v = v / getattr(el, '_errk', 1.0) if (key in ('K1', 'K2', 'K3', 'K4', 'K5', 'ANGLE', 'VOLT') and key == self.default_key(el)) else v
            raw = getattr(el, '_raw', None)
            if raw is not None and key in raw and is_num(raw[key]) and abs(v - raw[key]) <= 1e-13 * abs(raw[key]):
                return float(raw[key])              # 読み込んだ値のまま (L で割って掛け戻す丸めを消す)
            return v
        if key == 'INDEX' and not hasattr(el, 'INDEX'):
            return float((getattr(el, '_raw', None) or {}).get('INDEX', 0.0) or 0.0)     # 番号の初期値 (kwINDX)。gs: LINE["INDEX","QF"] は定義の値
        return float(getattr(el, _amap(el).get(key, key), 0.0) or 0.0)

    @staticmethod
    def _e_attr(el, key):
        """定義の E1/E2 に対応する pysad の属性名 (逆向き通過では入れ替わっている)"""
        if getattr(el, '_rev', False):
            return 'E2' if key == 'E1' else 'E1'
        return key

    def _e_ratio(self, el, key):
        """BEND の E1/E2 (SAD の値 = ANGLE に対する比)。pysad の属性は rad (= 比·ANGLE + AE)。
        ANGLE=0 のときは属性に比が残らないので _raw の値 (set_value が更新する) を返す。"""
        raw = getattr(el, '_raw', None) or {}
        ang = float(el.ANGLE)
        if ang == 0.0:
            return float(raw.get(key, 0.0) or 0.0)
        ae = float(getattr(el, 'AE' + key[1], 0.0) or 0.0)
        return (float(getattr(el, self._e_attr(el, key))) - ae) / ang

    @staticmethod
    def _e_as_fraction(el):
        raw = getattr(el, '_e_fraction', None)
        return bool(raw)

    def set_value(self, el, key: str, v: float):
        key = key.upper()
        tn = type(el).__name__
        if tn == 'Marker':
            raw = self.env._elements.setdefault(el.name, {'TYPE': 'MARK'})
            raw[key] = v
            if key == 'DDP' and self.els and el is self.els[0]:
                self.I.globals['DP0'] = float(v)    # DP0 ↔ LINE["DDP",1] (track.f の別名定義)
            if hasattr(el, _amap(el).get(key, key)):
                setattr(el, _amap(el).get(key, key), v)
            return
        if key in ('K1', 'K2', 'K3', 'K4', 'K5'):
            _sad_set(el, key, v * (getattr(el, '_errk', 1.0) if key == self.default_key(el) else 1.0))
        elif key == 'ANGLE' and tn == 'Bend':
            # E1/E2 は ANGLE に対する比なので、ANGLE を変えても比を保つ (SAD の psi1 = E1·ANGLE + AE1)
            r = {k: self._e_ratio(el, k) for k in ('E1', 'E2')}
            raw = getattr(el, '_raw', None)
            if raw is not None:
                for k in r:
                    raw[k] = r[k]
            setattr(el, key, v * getattr(el, '_errk', 1.0))
            for k in r:
                ae = float(getattr(el, 'AE' + k[1], 0.0) or 0.0)
                setattr(el, self._e_attr(el, k), r[k] * el.ANGLE + ae)
        elif key in ('E1', 'E2') and tn == 'Bend':
            raw = getattr(el, '_raw', None)
            if raw is not None:
                raw[key] = v
            ae = float(getattr(el, 'AE' + key[1], 0.0) or 0.0)
            setattr(el, self._e_attr(el, key), v * el.ANGLE + ae)
        elif key in ('ANGLE', 'VOLT'):
            setattr(el, key, v * getattr(el, '_errk', 1.0))
        elif key == 'L' and any(t == tn for t, _ in _PER_LENGTH):
            # SAD の K1.. は積分強度: L を変えても保たれる (Element["L","QF"]=.5 で LINE["K1","QF"] は変わらない)
            keep = {k: _sad_get(el, k) for t, k in _PER_LENGTH if t == tn}
            setattr(el, 'L', v)
            for k, val in keep.items():
                _sad_set(el, k, val)
        elif key == 'HARM':
            el._harm = int(v)
        else:
            setattr(el, _amap(el).get(key, key), v)
        # 幾何 (geo.py の _val は el._raw を先に見る): L/ROTATE/ANGLE/DX/DY/DZ を FFS (fit の変数を含む) で変えたら _raw も更新
        # (CALC のたびに tfgeo が計算し直す。11/05: L を足した — fit の free D L で幾何 GX.. の微分が 0 になっていた)
        _r = getattr(el, '_raw', None)
        if isinstance(_r, dict) and key in ('L', 'ROTATE', 'ANGLE', 'DX', 'DY', 'DZ') and key in _r:
            _r[key] = v * (getattr(el, '_errk', 1.0) if key == 'ANGLE' else 1.0)
        self.invalidate()

    def default_key(self, el) -> Optional[str]:
        vk = getattr(el, '_vary_key', None)
        if vk:
            return vk                            # VARY で変えた既定のキー
        tn = type(el).__name__
        if tn == 'Drift':
            return 'L'
        keys = _KEYS.get(tn)
        if tn == 'Marker' or not keys or len(keys) < 2:
            return None
        return keys[1]

    def set_family(self, fam: str, key: str, v: float):
        for el in self.instances(fam):
            self.set_value(el, key, v)
        for slave, (master, co) in self.couples.items():
            if master == fam:
                for el in self.instances(slave):
                    self.set_value(el, key, v * co)

    def invalidate(self):
        if self.rows is not None:
            self._last_rows = self.rows         # 光学の計算 (qcod) が失敗したとき SAD は直前の twiss を残す
        self.rows = None
        self.geo_dirty = True
        self._beam_cache = None
        if self.lat is not None:
            self.lat._twiss_cache = None

    # ------------------------------------------------------------------ 光学計算
    def _check_sol_optics(self, calc_only=False, emit=False):
        """SOL 領域の Twiss (CALC/INS) は solmap (追跡と同じ写像の数値微分)、Emittance[]・EMIT は solemit (放射は bz 付きの tradke)。
        CALC6D・RADCOD つきの CALC は未対応なので拒否する。SOL を無視して黙って違う値を出さない"""
        if self.els is not None:
            from . import und as _und
            _und.annotate(self.els, self.momentum(warn=False))       # v2: 磁場のある UND の写像 (k64 では何もしない)
            from . import phsrot as _phsrot
            _phsrot.annotate(self.els, self.momentum(warn=False))    # v2: PHSROT の写像
        if emit or calc_only:
            # (v2next B1: CALC6D・RADCOD つきの CALC も solmap で行う。calc6d._turn が _solmap の要素を solemit で通す)
            from . import solmap as _sm
            if self.els is not None and _sm.annotate(self.els, self.momentum(warn=False)):
                return
        for el in self.els:
            if type(el).__name__ == 'Solenoid' and (float(getattr(el, 'BZ', 0.0) or 0.0) != 0.0
                                                    or float(getattr(el, 'BOUND', 0.0) or 0.0) != 0.0):
                raise SADError('SOL (BZ≠0 か BOUND≠0) を含む CALC6D・RADCOD つきの CALC は未対応 (追跡 TrackParticles、CALC/INS の Twiss、Emittance[]・EMIT は対応)')

    def _apply_element_values(self):
        """ElementValues = {"K1"["QD"] -> expr, …} (Couple.n): CALC の前に要素値を設定する (K64 の tfinitvar → SetCoupledElements)。
        Packages (Couple.n) が無いときは何もできないので、ElementValues があれば SADError にする"""
        ev = self.I.globals.get('ElementValues')
        if not isinstance(ev, list) or not ev or getattr(self, '_in_ev', False):
            return
        from . import pkgload
        if pkgload.packages_dir() is None:
            raise SADError('ElementValues: SAD の Packages (Couple.n) が必要 (PYSAD_SAD_PACKAGES)')
        self._in_ev = True
        try:
            self.I.eval(Sym('SetCoupledElements'))
        finally:
            self._in_ev = False

    def calc(self, quiet: bool = False, ini: bool = False):
        """CELL/INS の光学計算。self.rows[l-1] = 要素 l の入口 (SAD の twiss(l))。
        CELL は不安定 (|cos μ|>1 等) でも SAD の qcell1 と同じく Twiss を作って続行する (判定は self.lat._qstat)。quiet=True は qcod の報告を出さない"""
        self._check_sol_optics(calc_only=True)
        if self.lat is None:
            raise SADError('No beam line (USE first)')
        self._apply_element_values()
        if self.flags.get('CALEXP', True):
            self.spread_pending()               # CALEXP: tffsadjustvar (NOCALEXP のときは成分ごとの値をそのまま使う)
        delta = self.dp0()
        self.calc6d_result = None
        self.geo_dirty = True              # tffscalc: tfgeo(.not. fitflg) — CALC のたびに幾何を軸から計算し直す (GEOCAL)
        if (not self.flags['CELL']) and self.flags['CALC6D']:
            # qcell (qcell.f:62) は calc6d なら CELL/INS に関わらず qcell6d に入る。INS (cell=.false.) の枝は閉軌道も tecalc も使わず、
            # 先頭 MARK の 28 成分 (AX..ZPY) と軌道 (DX..DDP) を初期条件に 6D で伝搬する (calc6d.calc6d_ins)。MARK には書き戻さない
            from .calc6d import calc6d_ins
            try:
                res = calc6d_ins(self.lat, self.momentum(), self._mark_tw28(delta), rfsw=self.flags['RFSW'],
                                 radcod=self.flags['RADCOD'], trf0=self.trf0, rad=self.flags['RAD'], fshift=self._fshift(), **self._mass_charge())
            except (ValueError, RuntimeError, FloatingPointError, ZeroDivisionError, np.linalg.LinAlgError) as e:
                self.write(' *****qcod---> Overflow & closed orbit not found            1' + chr(10))
                if self.rows is None:
                    self.rows = self._last_rows        # SAD は直前の twiss をそのまま残す
                if self.rows is None:
                    self.rows = _lat.nan_rows(self.lat)   # SAD は壊れた値で続行する (pysad は NaN)
                return self.rows
            self.rows = res.rows
            self._z_stale = self._zcols(res.rows[0])
            self.calc6d_result = res
            return self.rows
        sc_ = self._sc_opts() if _compat.v2() else None
        sc4d = bool(sc_ is not None and self.flags['CELL'] and not self.flags['CALC6D'] and self.flags['RFSW'])
        if self.flags['CELL'] and (self.flags['CALC6D'] or sc4d):
            # v2next B1: WSPAC (空間電荷) の CELL CALC は 6D の経路 (temit の反復でビーム行列を自己無撞着にし、空間電荷のレンズを含む T から
            # Twiss を作る)。4 次元の CALC (CALC6D なし) も RFSW ON のときはこの経路 (縦方向の列は 4D の既定値)。K64 (gs) は WSPAC ありの
            # CALC で入口以外の Twiss が NaN になる (k64 モードは今までどおり空間電荷を入れない)
            # qcell6d (qcell.f:404): CELL + CALC6D は temit の経路 (閉軌道 + 一周写像 + tecalc の正規形) で拡張 Twiss 28 成分を作る。
            # DP0≠0 のオフモメンタムも 6D の閉軌道に任せる (SAD も同じ)
            from .calc6d import calc6d
            prev = self.rows[0] if self.rows else None        # 不安定時の初期値 (SAD の直前の twiss(1))
            try:
                # v2 (C2): RADTAPER を CALC6D にも入れる (Emittance[] と同じ閉軌道・同じ係数)。K64 の qcell6d は rt=radtaper を渡すが、
                # gs は 'qcell6d---> Unstable optics' を出して入口の軌道がテーパーなしの値のまま (動作しない)。k64 モードは従来どおりテーパーなし
                res = calc6d(self._radtaper_lattice(), self.momentum(), rfsw=self.flags['RFSW'], radcod=self.flags['RADCOD'],
                             rad=self.flags['RAD'], twiss0=prev, out=self, fshift=self._fshift(), sc=sc_, **self._mass_charge())
            except (ValueError, RuntimeError, FloatingPointError, ZeroDivisionError, np.linalg.LinAlgError) as e:
                self.write(' *****qcod---> Overflow & closed orbit not found            1' + chr(10))
                if self.rows is None:
                    self.rows = self._last_rows        # SAD は直前の twiss をそのまま残す
                if self.rows is None:
                    self.rows = _lat.nan_rows(self.lat)   # SAD は壊れた値で続行する (pysad は NaN)
                return self.rows
            self.rows = res.rows
            if sc4d:
                for r_ in res.rows:            # 4 次元の CALC: 縦方向の列は既定値 (v2 B11)
                    for k_, v_ in zip(self._ZCOLS, (0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0)):
                        setattr(r_, k_, v_)
                self._z_stale = None
            else:
                self._z_stale = self._zcols(res.rows[0])
            self.calc6d_result = res
            self.trf0 = res.trf0
            if self.flags['CELL']:
                self._twmov_to_mark(res.rows)
                # track.f: DP0 は LINE["DDP",1] の別名。CELL CALC6D (RADCOD ON なら放射で動く閉軌道の δ) の後は DP0 が先頭の DDP になり、
                # 次の 4 次元の CALC はそのエネルギーで光学を計算する (gs で確かめた: RADCOD の CALC6D のあとの CALC のチューンが動く)
                d0_ = float(res.rows[0].delta_cod)
                if math.isfinite(d0_):
                    self.I.globals['DP0'] = d0_
                    self._dp0_from_c6d = d0_          # (RADTAPER の閉軌道の δ には使わない: _radtaper_lattice)
            return self.rows
        # 4 次元の CALC (qtwiss) では空洞はドリフト (RFSW に依らない: SAD で確認)
        _lat.RFSW = False
        if ini and not self.flags['CELL']:
            io = self.initial_orbits()             # CALC/GO の後の Twiss は InitialOrbits の中心 (idp=0) の軌道 (wake があればそのキック込み)
            if io is not None:
                self.rows = self.ini_optics(io, delta)[0]
                return self.rows
        try:
            if self.flags['CELL']:
                # tffscalc.f:93 twmov(..,.true.): 閉軌道探索 (qcod) の初期値は先頭 MARK の DX DPX DY DPY (直前の CELL の周期解か、
                # スクリプトが LINE["DX",..]= で書いた値)。収束判定が緩いので、初期値によって Twiss が 1e-7 程度変わる (SAD も同じ)
                m0 = self.els[0]
                cod_start = None
                if type(m0).__name__ == 'Marker':
                    raw0 = self.env._elements.get(m0.name, {})
                    cod_start = np.array([float(raw0.get(k, 0.0) or 0.0) for k in ('DX', 'DPX', 'DY', 'DPY')] + [0.0, delta])
                self.lat._mark_init = self._mark_init_raw()
                rows = list(self.lat.twiss_table(delta=delta, cod=cod_start,
                                                 sad_stab={'intres': bool(self.flags['INTRES']), 'halfres': bool(self.flags['HALFRES'])}))
                q = self.lat._qstat
                if q is None or (q['stabx'] and q['staby'] and not q['over']):
                    self._twmov_to_mark(rows)          # tffscalc: 全 idp が安定のときだけ MARK に書き戻す
            else:
                rows = self._ins_rows(delta)
        except (ValueError, RuntimeError, FloatingPointError, ZeroDivisionError) as e:
            if not quiet:
                self.write(' *****qcod---> Overflow & closed orbit not found            1' + chr(10))
            if self.rows is None:
                self.rows = self._last_rows            # SAD は直前の twiss をそのまま残す
            if self.rows is None:
                self.rows = _lat.nan_rows(self.lat)       # SAD は壊れた値で続行する (pysad は NaN)
            return self.rows
        # 4 次元の CALC は先頭行の縦方向の列 (AZ BZ NZ ZX ZPX ZY ZPY = 22〜28 列目) を作り直さない (tffscalc.f:90 は ntfun=15 列だけ入れ直す)。
        # 各行は qtwiss1 (qtwiss.f:536) で直前の行から写されるので、直前の 6D の CALC の先頭行の値が全行に残る (bench_level29c)
        zs = getattr(self, '_z_stale', None)
        if _compat.v2():
            zs = None              # v2 (B11): 4 次元の CALC は縦方向の列を既定値 (AZ=0, BZ=1, 他 0) にする (古い 6D の値を残さない)
        if zs is not None:
            for r in rows:
                for k, v in zip(self._ZCOLS, zs):
                    setattr(r, k, v)
        rows = self._rows_pos(self._mark_offsets(rows))
        self.rows = rows
        return self.rows

    _ZCOLS = ('alfz', 'betz', 'muz', 'zx', 'zpx', 'zy', 'zpy')
    # 先頭 MARK の 28 成分のキーワード (twmov の ntwissfun 全部。tffitcode の順)
    _MARK28 = ('AX', 'BX', 'PSIX', 'AY', 'BY', 'PSIY', 'EX', 'EPX', 'EY', 'EPY', 'R1', 'R2', 'R3', 'R4', 'DETR',
               'DX', 'DPX', 'DY', 'DPY', 'DZ', 'DDP', 'AZ', 'BZ', 'PSIZ', 'ZX', 'ZPX', 'ZY', 'ZPY')

    def _mark_tw28(self, delta):
        """INS + CALC6D の初期条件: 先頭 MARK の 28 成分 (tffscalc.f:93 twmov(1,twiss,..,.true.) = cmp%value(1:ntwissfun))。
        未設定のキーは 0、BX/BY/BZ は 1。DDP は DP0 (LINE["DDP",1] の別名)"""
        m0 = self.els[0]
        raw = self.env._elements.get(m0.name, {}) if type(m0).__name__ == 'Marker' else {}
        tw = []
        for k in self._MARK28:
            d = 1.0 if k in ('BX', 'BY', 'BZ') else 0.0
            tw.append(float(raw.get(k, d) or 0.0) if raw.get(k, None) is not None else d)
        tw[20] = float(delta)
        return tw

    def _zcols(self, row):
        return tuple(getattr(row, k) for k in self._ZCOLS)

    # SAD の twmov(1,twiss,nlat,ndim,.false.) が先頭 MARK に書く量: 周期解の入口 Twiss (TwissRow の属性) と閉軌道
    _TWMOV_ROW = (('AX', 'alfx'), ('BX', 'betx'), ('AY', 'alfy'), ('BY', 'bety'), ('EX', 'etx'), ('EPX', 'etpx'),
                  ('EY', 'ety'), ('EPY', 'etpy'), ('R1', 'r1'), ('R2', 'r2'), ('R3', 'r3'), ('R4', 'r4'), ('DETR', 'detr'),
                  ('DX', 'x_cod'), ('DPX', 'px_cod'), ('DY', 'y_cod'), ('DPY', 'py_cod'), ('DZ', 'z_cod'), ('DDP', 'delta_cod'))

    def _twmov_to_mark(self, rows):
        """SAD tffscalc.f:119 の直訳: CELL CALC が安定 (ここに来た = 閉軌道と周期解が求まった) なら、周期解の入口 Twiss と閉軌道を
        先頭要素 (MARK) の値に書き戻す (twmov(..,.false.))。次の INS CALC は tffscalc.f:93 の twmov(..,.true.) でこの値を入口条件に
        読む (CELL→INS 切替の直前の周期解が初期値になる)。PSIX=PSIY=0 (twiss(1,0,mfitnx)=0)。先頭が MARK でなければ何もしない"""
        m0 = self.els[0] if self.els else None
        if m0 is None or type(m0).__name__ != 'Marker' or not rows:
            return
        r0 = rows[0]
        vals = {k: float(getattr(r0, a, 0.0) or 0.0) for k, a in self._TWMOV_ROW}
        if not all(math.isfinite(v) for v in vals.values()):
            return
        vals['PSIX'] = 0.0
        vals['PSIY'] = 0.0
        # twmov(..,.false.) は 28 成分全部を書く。縦方向の列 (AZ BZ PSIZ ZX ZPX ZY ZPY) は 6D の CALC のときはその値、
        # 4D の CALC のときは直前の 6D の先頭行の値 (作り直されない: _z_stale、まだなければ AZ=0 BZ=1 他 0)
        zs = getattr(self, '_z_stale', None)
        if _compat.v2() and not self.flags.get('CALC6D'):
            zs = None              # v2 (B11): 4 次元の CALC のあとは MARK にも既定値を書く
        if zs is None:
            zs = (0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0)
        for k, v in zip(('AZ', 'BZ', 'PSIZ', 'ZX', 'ZPX', 'ZY', 'ZPY'), zs):
            vals[k] = float(v)
        vals['PSIZ'] = 0.0
        raw = self.env._elements.setdefault(m0.name, {'TYPE': 'MARK'})
        raw.update(vals)

    def _ins_start(self):
        """INS の入口条件: 先頭 MARK の Twiss (AX..DETR) と軌道 (DX DPX DY DPY DZ)。DDP は delta = DP0 と同じ"""
        m0 = self.els[0]
        raw = self.env._elements.get(m0.name, {}) if type(m0).__name__ == 'Marker' else {}
        g = lambda k, d=0.0: float(raw.get(k, d) or 0.0)
        tw0 = dict(ax=g('AX'), bx=g('BX', 1.0), ay=g('AY'), by=g('BY', 1.0), ex=g('EX'), epx=g('EPX'), ey=g('EY'), epy=g('EPY'),
                   r1=g('R1'), r2=g('R2'), r3=g('R3'), r4=g('R4'), detr=g('DETR'))
        return tw0, [g('DX'), g('DPX'), g('DY'), g('DPY'), g('DZ')]

    def _ins_rows(self, delta, cod0=None):
        tw0, cod_m = self._ins_start()
        # 初期軌道 = MARK の DX DPX DY DPY DZ (DDP は delta = DP0 と同じ: DP0 は LINE["DDP",1] の別名)
        cod0 = cod_m if cod0 is None else [float(v) for v in cod0[:5]]
        tws = _lat.ins_twiss(self.els, tw0, cod0, delta, acc=self._trpt_acc())
        return self._rows_from_tws(tws)

    def _rows_from_tws(self, tws):
        rows = []
        s = 0.0
        for i, tw in enumerate(tws):
            c = tw['cod']
            rows.append(TwissRow(s=s, name=('*' if i == 0 else self.els[i - 1].name),
                                 alfx=tw['ax'], betx=tw['bx'], mux=tw['nx'], alfy=tw['ay'], bety=tw['by'], muy=tw['ny'],
                                 etx=tw['ex'], etpx=tw['epx'], ety=tw.get('ey', 0.0), etpy=tw.get('epy', 0.0),
                                 r1=tw.get('r1', 0.0), r2=tw.get('r2', 0.0), r3=tw.get('r3', 0.0), r4=tw.get('r4', 0.0),
                                 detr=tw.get('detr', 0.0),
                                 x_cod=float(c[0]), px_cod=float(c[1]), y_cod=float(c[2]), py_cod=float(c[3]),
                                 z_cod=float(c[4]), delta_cod=float(c[5])))
            if i < len(self.els):
                s += self.els[i].L
        return rows

    # ------------------------------------------------------------------ InitialOrbits / WAKEOPT (tffsinitialcond, tffswake)
    def initial_orbits(self):
        """tffsinitialcond: シンボル InitialOrbits = {{x, px, y, py, z, dp}, ...} (6 個) か {{28 個の Twiss}, ...} のリスト。
        未定義 (シンボルのまま、Null) なら None。形が違えば SADError ('Error in InitialOrbits')。
        軌道は idp = nfr1..nfr (0 を除く) の順に入る (nfr=(n+1)/2、nfr1=−(n/2))。0 は中心の軌道 (ずれなし)。
        戻り: dict(nfr, nfam1, nfam, uini={idp: 28 個の配列}, has_tw={idp: 28 個形式で BX>0}, iuid={idp: 1 始まりの番号})"""
        v = self.I.get('InitialOrbits')
        if isinstance(v, Sym) or v is Null or v is None:
            return None
        if not isinstance(v, list) or not v:
            raise SADError('Error in InitialOrbits')
        n = len(v)
        nfr = (n + 1) // 2
        nfr1 = -(n // 2)
        uini, has_tw, iuid = {}, {}, {}
        j = 0
        for i in range(nfr1, nfr + 1):
            if i == 0:
                uini[0] = np.zeros(len(_TWISS_ALL))
                has_tw[0] = False
                iuid[0] = 0
                continue
            e = v[j]
            j += 1
            if not (isinstance(e, list) and all(is_num(x) for x in e)):
                raise SADError('Error in InitialOrbits')
            u = np.zeros(len(_TWISS_ALL))
            if len(e) == 6:
                u[15:21] = [float(x) for x in e]
                has_tw[i] = False
            elif len(e) == len(_TWISS_ALL):
                u[:] = [float(x) for x in e]
                has_tw[i] = u[1] > 0.0
            else:
                raise SADError('Error in InitialOrbits')
            uini[i] = u
            iuid[i] = j
        nfam = nfr
        nfam1 = -nfam if n % 2 == 0 else 1 - nfam
        return dict(nfr=nfr, nfam1=nfam1, nfam=nfam, uini=uini, has_tw=has_tw, iuid=iuid)

    def wake_optics_on(self):
        """tffsa.f:1186: wake = (TWAKE | LWAKE) & TRPT & WAKEOPT"""
        return bool((self.flags['TWAKE'] or self.flags['LWAKE']) and self.flags['TRPT'] and self.flags['WAKEOPT'])

    def _wake_opt_tables(self):
        """tffssetupwake: 各要素 i (1..nlat−1) の WakeFunction[Longitudinal/Transverse, 要素名]。{i: (縦, 横, DX, DY)}。
        LWAKE/TWAKE の ON/OFF は見ない (キックの側で見る)。どちらの表も無い要素は含めない"""
        out = {}
        for i in range(1, self.nlat):
            nm = self.elname(i)
            tabs = []
            for tag in ('Longitudinal', 'Transverse'):
                try:
                    r = self.I.eval(Expr(S('WakeFunction'), [S(tag), _Lit(nm)]))
                except SADError:
                    r = None
                a = None
                if isinstance(r, list) and r and all(isinstance(q, list) and len(q) == 2 and all(is_num(x) for x in q) for q in r):
                    a = np.array(r, dtype=float)
                tabs.append(a)
            if tabs[0] is None and tabs[1] is None:
                continue
            el = self.els[i - 1]
            out[i] = (tabs[0], tabs[1], float(getattr(el, 'DX', 0.0) or 0.0), float(getattr(el, 'DY', 0.0) or 0.0))
        return out

    def ini_optics(self, ini, dp0):
        """tffscalc の InitialOrbits の部分 (INS): 中心 (idp=0) と各軌道の Twiss と軌道を計算する。
        WAKEOPT (wake_optics_on) で WakeFunction を持つ要素があれば、要素の前後で全軌道が互いにキックを与え合う (wakeopt.wake_orbits)。
        そうでなければ各軌道は独立 (qcell1(fam=.true.) = 中心の入口の Twiss (28 個形式なら自分の Twiss) と、中心の入口の軌道 + ずれ)。
        戻り: {idp: TwissRow のリスト (位置 1..nlat)}"""
        if self.flags['CELL']:
            raise SADError('InitialOrbits: CELL は未対応 (INS のみ)')
        if self._trpt_acc() is not None:
            raise SADError('InitialOrbits: TRPT で加速要素があるときは未対応')
        tw0, cod_m = self._ins_start()
        tw_keys = ('ax', 'bx', 'nx', 'ay', 'by', 'ny', 'ex', 'epx', 'ey', 'epy', 'r1', 'r2', 'r3', 'r4', 'detr')
        slots = {}
        for idp in range(ini['nfam1'], ini['nfam'] + 1):
            u = ini['uini'][idp]
            if idp != 0 and ini['has_tw'][idp]:
                tw = {k: float(u[a]) for a, k in enumerate(tw_keys)}
            else:
                tw = dict(tw0)
                tw['nx'] = tw['ny'] = 0.0
            cod =np.array(cod_m + [dp0], dtype=float) + (u[15:21] if idp != 0 else 0.0)
            slots[idp] = (tw, cod)
        wake = None
        if self.wake_optics_on():
            tables = self._wake_opt_tables()
            if tables:
                wake = tables
        out = {}
        if wake is None:
            for idp, (tw, cod) in slots.items():
                tw = dict(tw)
                if idp != 0:
                    tw['nx'] = tw['ny'] = 0.0             # tffscalc.f:188
                out[idp] = self._rows_from_tws(_lat.ins_twiss(self.els, tw, cod[:5], float(cod[5])))
            return out
        from .wakeopt import wake_orbits
        g = self.I.globals
        amass = float(g.get('MASS', 0.51099895e6) or 0.51099895e6)
        res = wake_orbits(self.els, slots, wake, sigz=float(g.get('SIGZ', 0.0) or 0.0),
                          nbunch=float(g.get('NBUNCH', 1.0) or 1.0), charge=float(g.get('CHARGE', 1.0) or 1.0),
                          pbunch=float(g.get('PBUNCH', 0.0) or 0.0), lwake=bool(self.flags['LWAKE']),
                          twake=bool(self.flags['TWAKE']), gammab=self._gammab_at, amass=amass)
        for idp, d in res.items():
            out[idp] = self._rows_from_tws([d[k] for k in range(1, self.nlat + 1)])
        return out

    def _gammab_at(self, k, fr=0.0):
        """gammab(k) = 基準運動量 p/mc (単位 MASS)。TRPT かつ RFSW で加速要素があれば tfgeo1 の分布 (trptacc)、なければ一定。
        tftwiss.f: rlist(ifgamm+ia-1)*(1-v) + rlist(ifgamm+min(nlat-1,ia))*v"""
        amass = float(self.I.get('MASS', 0.51099895e6) or 0.51099895e6)
        p0 = self.momentum() / amass
        acc = self._trpt_acc()
        if acc is None:
            return p0
        gb = acc['gb']
        k = max(1, min(self.nlat, int(k)))
        a = float(gb[k - 1])
        b = float(gb[min(self.nlat - 1, k)])
        return a * (1.0 - fr) + b * fr

    def _trpt_acc(self):
        """TRPT (輸送線) で RFSW ON のときの加速 (trptacc.make_acc)。qcav は rfsw = rfsw .and. trpt なので RING では None。
        基準運動量 p0 = MOMENTUM/MASS、trf0 = 0 (qtwiss.f:87: trpt なら trf0=0)、omega0 は全長から (tphyzp)"""
        if not (self.flags['TRPT'] and self.flags['RFSW']):
            return None
        from . import trptacc
        g = self.I.globals
        amass = float(g.get('MASS', 0.51099895e6) or 0.51099895e6)
        charge = float(g.get('CHARGE', 1.0) or 1.0)
        p0 = self.momentum() / amass
        circ = sum(el.L for el in self.els)
        h0 = math.sqrt(1.0 + p0**2)
        omega0 = 2.0 * math.pi * 299792458.0 * p0 / h0 / circ if circ != 0.0 else 0.0
        return trptacc.make_acc(self.els, p0, charge=charge, amass=amass, omega0=omega0, trf0=0.0,
                                radcod=bool(self.flags['RADCOD']), rfsw=True)

    def _mark_offsets(self, rows):
        """qcell.f:318-345: OFFSET≠0 の MARK は位置 xp = tffselmoffset(l) (前後の要素の中の分数位置、MARK を指すなら辿る) の
        Twiss を自分の行にする (twiss(l,:) = twiss(lx,:) か qtwissfrac1(lx, fr))。行の s はその位置"""
        nlat = self.nlat
        for l in range(1, nlat):
            el = self.els[l - 1]
            if type(el).__name__ != 'Marker':
                continue
            xp = self._elm_offset(l)
            if xp == float(l) or not (1.0 <= xp < float(nlat)):
                continue
            lx = int(xp)
            fr = xp - lx
            if fr == 0.0:
                r = _lat._clone_rows([rows[lx - 1]])[0]
                r.name = rows[l - 1].name
            else:
                r = self.row_frac(lx, fr, rows=rows)
                r.name = rows[l - 1].name
            rows[l - 1] = r
        return rows

    def _mark_offset(self, l):
        """tffs.f tffsmarkoffset: MARK の OFFSET (逆向きなら (OFFSET−.5)·(−1)+.5)。nlat は 0"""
        if l >= self.nlat:
            return 0.0
        el = self.els[l - 1]
        off = _mark_offset_value(el)
        if off == 0.0:
            return 0.0
        if getattr(el, '_rev', False):
            off = (off - 0.5) * -1.0 + 0.5
        return off

    def _elm_offset(self, l):
        """tffs.f tffselmoffset: 要素 l の表示・fit の位置 xp (OFFSET 付き MARK は offset+l、それが MARK なら辿る、最大 256 回)"""
        xp = float(l)
        if l == self.nlat:
            return xp
        lm = l
        nm = 0
        while True:
            off = self._mark_offset(lm)
            if off == 0.0:
                break
            xp = off + lm
            if 1.0 <= xp <= float(self.nlat):
                lx = int(xp)
                if type(self.els[lx - 1]).__name__ == 'Marker' if lx < self.nlat else False:
                    nm += 1
                    if nm < 256:
                        lm = lx
                        continue
                    self.write(f'?Recursive OFFSET in {self.elname(l)}' + chr(10))
            break
        return xp

    def _rows_pos(self, rows):
        """SAD の twiss の s (pos) は tfgeo の pos: BEND/MULT のフリンジの長さ補正と OFFSET MARK の位置を含む。
        それらが無ければ要素長の累積と同じなので触らない"""
        need = any(_mark_offset_value(el) != 0.0 for el in self.els if type(el).__name__ == 'Marker')
        if not need:
            for el in self.els:
                if type(el).__name__ in ('Bend', 'Multipole') and float(getattr(el, 'FRINGE', 0.0) or 0.0) != 0.0 \
                        and any(float(getattr(el, k, 0.0) or 0.0) for k in ('F1', 'FB1', 'FB2')):
                    need = True
                    break
        if not need:
            return rows
        pos = self.ensure_geo().pos
        for i, r in enumerate(rows):
            r.s = float(pos[i])
        return rows

    def ensure_geo(self):
        """tfgeo: 要素の幾何 (geo.compute_geo)。GEOCAL (既定) なら CALC/USE/ORG などの後に計算し直し、
        GEOFIX なら最初の USE 直後の結果を保持する (tffsa.f:114-118: USE 直後は GEOFIX でも一度だけ計算)"""
        if self.lat is None:
            raise SADError('No beam line (USE first)')
        if self.geo is None or (self.geo_dirty and self.flags['GEOCAL']):
            from . import geo as _geo

            def warn(m):
                if m not in self._geo_warned:
                    self._geo_warned.add(m)
                    self.write(f'[pysad] 幾何: {m}' + chr(10))
            gx, gy, gz = self.org_g
            self.geo = _geo.compute_geo(self.els, gx0=gx, gy0=gy, gz0=gz, chi0=tuple(self.org_chi), iorg=self.org_k, warn=warn, momentum_eV=self.momentum(warn=False))
            self.geo_dirty = False
        return self.geo

    def _org(self, words, i):
        """ORG [要素 dgx dgy dgz dchi1 dchi2 dchi3] (tffsa.f ORI_GIN, chi は度)。要素名なしは初期設定に戻す
        (次の語はコマンドとして処理される)。値は getvad と同じく「+/- で始まる語は前の式の続き」と読む
        (`ORG B3 1.5 -2.5 0.5 ..` は gx=1.5-2.5=-1.0)。chi の欠落は ?ORG-Missing parameters を出すが、それまでの値で設定する"""
        nxt = words[i + 1] if i + 1 < len(words) else ''
        pos = []
        if nxt and not any(c in nxt for c in '*%{|?') and not re.match(r'^[-+.0-9]', nxt):
            pos = self.positions(nxt)
        if not pos:
            self.org_k, self.org_g, self.org_chi = 1, [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]
            self.geo_dirty = True
            self.write(' ???-FFS-Info-ORG-Reset to the initial setting. ' + chr(10))
            return 1
        # 式のグループ化 (各グループは語のリスト)。+/- などで始まる語、または演算子で終わる語の次の語は前の式の続き
        groups = []
        for w in words[i + 2:]:
            # v2 (B26): 空白で区切った +/- で始まる語は別の値 (`ORG B3 1.5 -2.5 ..` は gx=1.5, gy=-2.5)。K64 は前の式の続き (1.5-2.5)
            if groups and (w[0] in ('*/^' if _compat.v2() else '+-*/^') or groups[-1][-1][-1] in '+-*/^('):
                groups[-1].append(w)
            else:
                groups.append([w])
        vals = []
        consumed = 2
        for g in groups:
            if len(vals) >= 6:
                break
            v, _ = self._getva([''.join(g)], 0)
            if v is None or isinstance(v, list):
                break
            vals.append(float(v))
            consumed += len(g)
        if len(vals) < 3:
            self.write(' ???-FFS-Error-?ORG-Missing parameters: ORG dgx dgy dgz dchi1 dchi2 dch3. ' + chr(10))
            return consumed
        g0 = vals[:3]
        chi = (vals[3:6] + [0.0, 0.0, 0.0])[:3]
        self.org_k = pos[0]
        self.org_g = g0
        self.org_chi = [c * math.pi / 180.0 for c in chi]
        self.geo_dirty = True
        if len(vals) < 6:
            self.write(' ???-FFS-Error-?ORG-Missing parameters: ORG dgx dgy dgz dchi1 dchi2 dch3. ' + chr(10))
        return consumed

    def ensure_rows(self):
        if self.rows is None:
            self.calc()
        return self.rows

    def row(self, k):
        """位置 k (1 始まり) の入口の Twiss"""
        rows = self.ensure_rows()
        return rows[k - 1]

    @staticmethod
    def frac_element(el, fr):
        """SAD qfraccomp (qtwiss.f) の直訳: 要素 el の入口から割合 fr (rx1=0, rx2=fr) までの部分要素を作る。
        f1=1 (入口側あり), f2=0 (出口側なし)。pysad の K1/K2/K3 は長さあたりなので L を縮めれば積分値は r 倍になる。
        BEND: ANGLE·r, K0·r, E1 はそのまま (E1·f1/r × ANGLE·r = E1·ANGLE), E2=0,
              FRMD=0 なら F1=0、FRMD≠0 なら FRMD=-f1-2f2=-1。
        共通: fr0=FRMD (0 なら 3) について fr0∈{0,2} なら FB1 (QUAD では F1) を 0、fr0∈{0,1} なら FB2 (F2) を 0、
              FRMD = 1 (fr0∈{3,1} のとき入口フリンジ) に置き換える。
        CAVI: VOLT·r。MULT: K·r。逆向き要素 (ori=.false.) の E1/E2 の入れ替えは未対応。"""
        import copy
        part = copy.copy(el)
        r = fr
        part.L = el.L * r
        tn = type(el).__name__
        g = lambda k: float(getattr(el, k, 0.0) or 0.0)
        if tn == 'Drift' or tn == 'Marker':
            return part
        if tn == 'Bend':
            if g('FRINGE') == 0.0:
                part.F1 = 0.0
            else:
                part.FRINGE = -1.0
            part.E2 = 0.0
            part.ANGLE = el.ANGLE * r
            if hasattr(el, 'K0'):
                part.K0 = g('K0') * r
            return part
        if tn == 'Multipole':
            knl = getattr(el, 'KnL', None)
            if isinstance(knl, dict):
                part.KnL = {k: v * r for k, v in knl.items()}
            ksl = getattr(el, 'KsL', None)
            if isinstance(ksl, dict):
                part.KsL = {k: v * r for k, v in ksl.items()}
            for k in ('VOLT', 'DVOLT', 'W1'):
                if hasattr(el, k):
                    setattr(part, k, g(k) * r)
            if g('ANGLE') != 0.0:
                part.ANGLE = g('ANGLE') * r
                part.E2 = 0.0
                return part
        elif tn == 'Cavity':
            for k in ('VOLT', 'DVOLT', 'RANV', 'V1'):
                if hasattr(el, k):
                    setattr(part, k, g(k) * r)
        elif tn == 'Quadrupole':
            pass
        else:
            return part
        # 共通部 (QUAD/CAVI/MULT(ANGLE=0))。QUAD の FB1/FB2 は F1/F2
        fr0 = g('FRINGE')
        fb1, fb2 = ('F1', 'F2') if tn == 'Quadrupole' else ('FB1', 'FB2')
        if fr0 == 0.0 or fr0 == 2.0:
            if hasattr(part, fb1):
                setattr(part, fb1, 0.0)
        if fr0 == 0.0 or fr0 == 1.0:
            if hasattr(part, fb2):
                setattr(part, fb2, 0.0)
        if fr0 == 0.0:
            fr0 = 3.0
        part.FRINGE = 1.0 if fr0 in (3.0, 1.0) else 0.0
        return part

    def row_frac(self, k, fr, rows=None):
        """要素 k の中 (入口から fr の割合) の Twiss。SAD qtwissfrac の直訳: qfraccomp で部分要素を作り、
        入口の閉軌道から qtwiss (閉軌道まわりの転送行列) で伝搬する。rows: 使う Twiss の行 (既定は ensure_rows)"""
        if fr == 0.0 or k >= self.nlat:
            return self.row(k) if rows is None else rows[k - 1]
        res6 = getattr(self, 'calc6d_result', None)
        if rows is None and res6 is not None and self.rows is res6.rows:
            from .calc6d import twiss_frac            # CALC6D の分数位置 (qtwissfrac の calc6d 分岐)
            return twiss_frac(res6, k, fr)
        el = self.els[k - 1]
        r0 = self.row(k) if rows is None else rows[k - 1]
        part = self.frac_element(el, fr)
        delta = self.dp0()
        cod_in = np.array([r0.x_cod, r0.px_cod, r0.y_cod, r0.py_cod, r0.z_cod, r0.delta_cod], dtype=float)
        _lat.RFSW = False
        _, cod_out, el_mats = _lat._one_turn_cod_matrix([part], delta, cod_in)
        M = el_mats[0][0]
        c = cod_out
        tw0 = dict(ax=r0.alfx, bx=r0.betx, nx=r0.mux, ay=r0.alfy, by=r0.bety, ny=r0.muy, ex=r0.etx, epx=r0.etpx, ey=r0.ety, epy=r0.etpy,
                   r1=r0.r1, r2=r0.r2, r3=r0.r3, r4=r0.r4, detr=r0.detr)
        tw = _lat._transport_twiss(tw0, M)
        return TwissRow(s=r0.s + part.L, name=el.name, alfx=tw['ax'], betx=tw['bx'], mux=tw['nx'], alfy=tw['ay'], bety=tw['by'],
                        muy=tw['ny'], etx=tw['ex'], etpx=tw['epx'], ety=tw.get('ey', 0.0), etpy=tw.get('epy', 0.0),
                        r1=tw.get('r1', 0.0), r2=tw.get('r2', 0.0), r3=tw.get('r3', 0.0), r4=tw.get('r4', 0.0), detr=tw.get('detr', 0.0),
                        x_cod=float(c[0]), px_cod=float(c[1]), y_cod=float(c[2]), py_cod=float(c[3]), z_cod=float(c[4]), delta_cod=float(c[5]))

    _TWKEY = {'AX': 'alfx', 'BX': 'betx', 'NX': 'mux', 'AY': 'alfy', 'BY': 'bety', 'NY': 'muy', 'EX': 'etx', 'EPX': 'etpx',
              'EY': 'ety', 'EPY': 'etpy', 'R1': 'r1', 'R2': 'r2', 'R3': 'r3', 'R4': 'r4', 'DETR': 'detr',
              'DX': 'x_cod', 'DPX': 'px_cod', 'DY': 'y_cod', 'DPY': 'py_cod', 'DZ': 'z_cod', 'DDP': 'delta_cod',
              'AZ': 'alfz', 'BZ': 'betz', 'NZ': 'muz', 'ZX': 'zx', 'ZPX': 'zpx', 'ZY': 'zy', 'ZPY': 'zpy', 'S': 's', 'LENG': 's'}

    _TWKEY_ALL = frozenset(list(_TWKEY) + ['GMX', 'GMY', 'GMZ'])

    def twiss_value(self, key: str, k: int, fr: float = 0.0):
        key = key.upper()
        r = self.row_frac(k, fr) if fr else self.row(k)
        if key in self._TWKEY:
            return float(getattr(r, self._TWKEY[key], 0.0))
        if key == 'GMX':
            return (1 + r.alfx**2) / r.betx            # γ 関数 (SAD の GMX/GMY)
        if key == 'GMY':
            return (1 + r.alfy**2) / r.bety
        if key in ('GX', 'GY', 'GZ', 'CHI1', 'CHI2', 'CHI3'):
            # K64 では Twiss["GX"] は光学関数名ではない (幾何は LINE["GX"..]) — tftwiss.f
            raise SADError(f'General::wrongval:  #1 ({key}) is must be to be name of optical function in Twiss["{key}"]')
        if key.startswith('SIG') and key[3:4] != 'E' and not key.startswith('SIGMA'):
            return self.beam_size_value(key[3:], k, fr)          # tftwiss.f:270: SIG.. は tfline (ビームサイズ) へ
        if len(key) > 1 and key[0] in 'RD' and key[1:] in self._TWKEY_ALL:
            # tftwiss.f: R… は参照光学 (twiss(:,-1,:))、D… は現在 - 参照 (BX BY BZ は比)。小数の位置は不可
            if fr:
                raise SADError('General::wrongval:  #2 is must be to be integer position for reference optics')
            ref = getattr(self, '_ref_rows', None)
            r0 = ref[k - 1] if ref is not None else _ZeroRow()
            rv = self._disp_val(r0, key[1:]) / (PI2 if False else 1.0)
            if key[0] == 'R':
                return rv
            cv = self._disp_val(r, key[1:])
            if key[1:] in ('BX', 'BY', 'BZ'):
                return (cv - rv) / rv if rv else float('nan')
            return cv - rv
        if key in ('GAMMA', 'GAMMABETA', 'S') or key.startswith('SIZE'):
            return self.line_value(key, k, fr)                # tftwiss.f:271: Twiss[…] はこれらを tfline (LINE と同じ) へ回す
        if key == 'PEX':
            return r.etx
        if key == 'PEPX':
            return r.etpx
        if key == 'PEY':
            return r.ety
        if key == 'PEPY':
            return r.etpy
        raise SADError(f'General::wrongval:  #1 ({key}) is must be to be name of optical function in Twiss["{key}"]')

    # ------------------------------------------------------------------ 組み込み: Twiss[], LINE[]
    def _resolve_loc(self, loc):
        """位置引数 → [(k, fr)]。裸のシンボル (`LINE["K1", Q1A]` の Q1A: 未定義のシンボル) は SAD の tfgetstrns と同じく名前として扱う"""
        if isinstance(loc, Sym):
            loc = loc.name
        if is_num(loc):
            v = float(loc)
            k = int(math.floor(v))
            fr = v - k
            if k < 0:
                k = self.nlat + k + 1
            while k > self.nlat:             # tflinestk: nlat を超える番号は周回する (k − nlat)
                k -= self.nlat
            return [(max(1, min(self.nlat, k)), fr)]
        if isinstance(loc, str):
            s = loc.strip()
            m = re.match(r'^(.*?)([+-][0-9.]+)$', s)
            off, fr = 0, 0.0
            if m and m.group(1) and not any(c in m.group(1) for c in '*%{|'):
                s = m.group(1)
                v = float(m.group(2))
                off = int(math.floor(v))
                fr = v - off
            return [(k + off, fr) for k in self.positions(s)]
        if isinstance(loc, list):
            out = []
            for x in loc:
                out.extend(self._resolve_loc(x))
            return out
        raise SADError(f'bad location {loc!r}')

    def fn_twiss(self, args):
        if not args:
            raise SADError('Twiss: keyword required')
        key = args[0]
        if isinstance(key, list):
            return [self.fn_twiss([k] + list(args[1:])) for k in key]
        key = key.upper()
        if key == 'LENGTH':
            if len(args) > 1:
                raise SADError(f'General::narg:  Number of arguments is expected to be 1 in Twiss["LENGTH",{self.I.to_str(args[1], None, True)}]')
            return float(self.nlat)
        if key in ('*', 'ALL', 'RALL', 'DALL'):
            # tftwiss.f: Twiss["*"] / Twiss["*", loc] は ntwissfun=28 個の光学関数 (tffs.f mfitax..mfitzpy の順) のリスト。RALL は参照光学、DALL は現在 - 参照
            pre = key[0] if key in ('RALL', 'DALL') else ''
            nm_ = lambda k2: (pre + k2) if pre else k2
            if len(args) == 1:
                return [[self.twiss_value(nm_(k2), k) for k2 in _TWISS_ALL] for k in range(1, self.nlat + 1)]
            locs = self._resolve_loc(args[1])
            vals = [[self.twiss_value(nm_(k2), k, fr) for k2 in _TWISS_ALL] for k, fr in locs]
            # tftwiss.f tflinestk: 一致が 1 個 (isp == isp0+1) ならパターンの形に関わらずスカラー側 (リストでない)
            return vals[0] if len(vals) == 1 and not isinstance(args[1], list) else vals
        if len(args) == 1:
            return [self.twiss_value(key, k) for k in range(1, self.nlat + 1)]
        locs = self._resolve_loc(args[1])
        vals = [self.twiss_value(key, k, fr) for k, fr in locs]
        # tftwiss.f: tflinestk の一致数が 1 ならスカラー (`***`、`RFC*` など 1 個だけ一致するパターンも)。0 個は {}
        if len(vals) == 1 and not isinstance(args[1], list):
            return vals[0]
        return vals

    def line_value(self, key: str, k: int, fr: float = 0.0):
        key = key.upper()
        if len(key) > 4 and key.endswith('$SUM'):
            # tffs.f tfkeyv: KEY$SUM は KEY の値。値が実数のリスト (SPLIT した要素の分割ごとの値) ならその和 (pysad は SPLIT を持たないので値そのもの)
            v = self.line_value(key[:-4], k, fr)
            return float(sum(v)) if isinstance(v, list) and all(isinstance(t, (int, float)) for t in v) else v
        if key == 'NAME':
            return self.elname(k)
        if key == 'ELEMENT':
            return '$$$' if k >= self.nlat else self.els[k - 1].name
        if key == 'POSITION':
            return float(k)
        if key.startswith('SIG') and key[3:4] != 'E' and not key.startswith('SIGMA'):
            # tftwiss.f:982: LINE["SIG"//key, l] = ビーム行列 (key 空 → 6x6、X → σx、XY → Σxy 等)。tfsize: CALC6D なら temit の
            # 平衡ビーム行列、そうでなければ tfbeam (EMIT の EMITX/EMITY/SIGE と Twiss から)
            return self.beam_size_value(key[3:], k, fr)
        if key in ('GX', 'GY', 'GZ', 'GCHI1', 'GCHI2', 'GCHI3', 'GEO'):
            return self.ensure_geo().line_value(key, k, fr)
        if key in ('OGX', 'OGY', 'OGZ', 'OCHI1', 'OCHI2', 'OCHI3', 'OGEO'):
            # tforbitgeo: 閉軌道を乗せた幾何 (閉軌道は CALC の twiss の DX/DPX/DY/DPY)
            rw = self.row_frac(k, fr) if fr else self.row(k)
            return self.ensure_geo().orbit_value(key, k, (rw.x_cod, rw.px_cod, rw.y_cod, rw.py_cod), fr)
        if key in ('S', 'LENG'):
            try:
                return float(self.ensure_geo().pos_frac(k, fr))       # tfgeo の pos (フリンジによる長さ補正込み)
            except Exception:
                r = self.ensure_rows()
                s0 = r[k - 1].s
                s1 = r[min(k, self.nlat - 1)].s if k < self.nlat else s0
                return s0 * (1 - fr) + s1 * fr
        if key == 'TYPE':
            return 0.0 if k >= self.nlat else float(_TYPECODE.get(_cls(self.els[k - 1]), 0))
        if key == 'TYPENAME':
            return '' if k >= self.nlat else _TYPENAME.get(_cls(self.els[k - 1]), '')
        if key == 'DK':
            return 0.0 if k >= self.nlat else getattr(self.els[k - 1], '_errk', 1.0)
        if key == 'MULT':
            # tfinimult: 同名が複数ある要素は周回順の番号 (NAME.n の n)、1 台だけ (と $$$) は 0
            if k >= self.nlat or len(self.instances(self.els[k - 1].name)) < 2:
                return 0.0
            return float(self._ordinal(k - 1))
        if key == 'DIR':
            # tftwiss.f:942 direlc: 反転配置 (LINE の `-`、cmp%ori=.false.) は -1
            return -1.0 if (0 < k < self.nlat and getattr(self.els[k - 1], '_rev', False)) else 1.0
        if key in ('GAMMA', 'GAMMABETA'):
            # tftwiss.f:837-842: gammab(k) (基準運動量 p/mc、TRPT の加速があればその分布) を fr で補間。GAMMA は sqrt(1+gammab²)
            gb = self._gammab_at(k, fr)
            return gb if key == 'GAMMABETA' else math.sqrt(1.0 + gb * gb)
        if k >= self.nlat:
            return 0.0
        el = self.els[k - 1]
        if key in ('EMITX', 'EMITY', 'EMITZ', 'SIGMAZ', 'SIGZ', 'SIGE', 'DP') and type(el).__name__ != 'Marker':
            return 0.0
        if key == 'SIGMAZ':
            key = 'SIGZ'
        v = self.get_value(el, key)
        ek = getattr(el, '_errk', 1.0)
        if ek != 1.0 and key in ('K1', 'K2', 'K3', 'K4', 'K5', 'ANGLE', 'VOLT') and key == self.default_key(el):
            v = v * ek
        return v

    def fn_synchrobeta(self, args):
        """SynchroBetaEmittance1[nus0, nus1, dnus, mphi] (tfemit.f tfemits → temits.f)。Packages/emit.n の SynchroBetaEmittance が呼ぶ。
        戻り: {{nus, emitx, emity, emitxp, emityp, conv}, …} (pysad/synchrobeta.py)"""
        if len(args) != 4 or not all(is_num(a) for a in args):
            raise SADError('General::wrongtype: "[{nusstart,nusstop,nusstep},mphi]"')
        if self.lat is None:
            raise SADError('No beam line (USE first)')
        from .synchrobeta import synchro_beta_emittance
        res = synchro_beta_emittance(self.lat, self.momentum(), args[0] * 2.0 * math.pi, args[1] * 2.0 * math.pi,
                                     args[2] * 2.0 * math.pi, int(args[3]))
        if _compat.v2() and int(args[3]) >= 2 and int(args[3]) < 32:
            # D1: 方位角モード数 (AzimuthalModes) について収束していない値を黙って返さない。大きなリング A (3 GeV、周長 603 m、νs = 0.004) で
            # 既定の 9 や 3〜8 の値 (6.18e-9) が 12 で 3.2e-9、20 で 1.5e-9、45 で 7.9e-10 と動き続け、Emittance[] (1.55e-10)・追跡と合わなかった。
            # 1.5 倍のモード数でもう一度解き、emitx が 20% 以上違えば警告する (値はそのまま返す)
            try:
                m0 = int(args[3])
                m1 = min(32, max(m0 + 3, int(math.ceil(1.5 * m0))))
                res1 = synchro_beta_emittance(self.lat, self.momentum(), args[0] * 2.0 * math.pi, args[0] * 2.0 * math.pi, 0.0, m1)
                a, b2 = abs(res[0][1]), abs(res1[0][1])
                if b2 > 0.0 and abs(a - b2) > 0.2 * max(a, b2):
                    self.write(f' (pysad) SynchroBetaEmittance: AzimuthalModes->{m0} の emitx {res[0][1]:.4g} は AzimuthalModes->{m1} の値 {res1[0][1]:.4g} と'
                               f' {abs(a - b2) / max(a, b2) * 100:.0f}% 違う (モード数について収束していない。Emittance[] と追跡で確かめること)' + chr(10))
            except Exception:
                pass
        return [[float(v) for v in row] for row in res]

    def fn_line(self, args):
        key = args[0]
        if isinstance(key, list):
            return [self.fn_line([k] + list(args[1:])) for k in key]
        key = key.upper()
        if key == 'LENGTH':
            if len(args) > 1:
                if _compat.v2():
                    # v2 (B24): K64 のソース (tftwiss.f) も LENGTH は引数 1 個だけ (要素数を返す)。エラーのままにして、位置の取り方を知らせる
                    raise SADError(f'General::narg:  Number of arguments is expected to be 1 in LINE["LENGTH","{args[1]}"]'
                                   f' (LINE["LENGTH"] is the number of components; use LINE["POSITION",loc] for the index, LINE["LENG",loc] or LINE["S",loc] for the orbit length)')
                raise SADError(f'General::narg:  Number of arguments is expected to be 1 in LINE["LENGTH","{args[1]}"]')
            return float(self.nlat)
        if key == 'EXPAND':
            return Null
        if len(args) == 1:
            return [self.line_value(key, k) for k in range(1, self.nlat + 1)]
        locs = self._resolve_loc(args[1])
        vals = [self.line_value(key, k, fr) for k, fr in locs]
        # tftwiss.f tfline: m == 1 ならスカラー、それ以外 (0 個を含む) はリスト
        if len(vals) == 1 and not isinstance(args[1], list):
            return vals[0]
        return vals


    def fn_master(self, args):
        """Master[] (pgmast.f): 各要素の「マスター」。分割された要素 (同じ icomp が続くもの) の先頭は最後の分割の位置、残りは -先頭。
        BEND/QUAD/SEXT/OCT/MULT(22)/31/SOL 以外は 0。pysad は SPLIT を持たないので、対象の要素は自分自身 (長さ 0 の要素が続いても同じ)"""
        if len(args) > 1:
            raise SADError('General::narg: Number of arguments is expected to be 1 in Master')
        n = self.nlat
        out = [0.0] * n
        for i in range(1, n):
            id_ = int(self.line_value('TYPE', i))
            if (id_ > 8 and id_ != 31 and id_ != 22 and id_ != 20) or id_ == 1:
                continue
            out[i - 1] = float(i)
        return out

    def fn_response_matrix(self, args):
        """ResponseMatrix[mon, st, stt, opt, periodic, kind (, tune)] (pgrmat.f)。mon がリストなら行ごと"""
        from . import corr as _corr
        if len(args) < 6 or len(args) > 7:
            raise SADError('General::narg: Number of arguments is expected to be 6 or 7 in ResponseMatrix')
        mon = args[0]
        if isinstance(mon, list):
            return [self.fn_response_matrix([m] + list(args[1:])) for m in mon]
        st, stt, opt, periodic, kind = args[1], args[2], args[3], args[4], args[5]
        if not isinstance(opt, list) or not opt or not isinstance(opt[0], list):
            raise SADError('General::wrongtype: Argument must be List for #4')
        if not is_num(periodic):
            raise SADError('General::wrongtype: Argument must be Real for #5')
        if not isinstance(kind, str):
            raise SADError('General::wrongtype: Argument must be String for #6')
        tune = [0.0, 0.0]
        if len(args) == 7:
            t = args[6]
            if not isinstance(t, list) or len(t) != 2:
                raise SADError('General::wrongleng: #7 must be a list of length 2')
            tune = [float(t[0]), float(t[1])]
        locs = self._resolve_loc(mon)
        if not locs:
            raise SADError('General::wrongval: ResponseMatrix: invalid element for #1')
        im = int(locs[0][0])
        rows = opt[0]

        def rot(js):
            if js < 1 or js > self.nlat:
                return None
            if int(self.line_value('TYPE', js)) == 2:
                return -float(self.line_value('ROTATE', js))
            return None
        scalar = not isinstance(st, list)
        stl = [st] if scalar else st
        sttl = [stt] if not isinstance(stt, list) else stt
        if len(sttl) < len(stl):
            raise SADError('General::equalleng: #2 and #3')
        ar = _corr.response_row(im, stl, sttl, rows, rot, bool(periodic), kind, tune[0], tune[1])
        if self.flags['TRPT']:
            # pgrmat.f:354・498: TRPT (加速のある輸送路) は sqrt(gammab(補正磁石)/gammab(モニタ)) を掛ける (断熱減衰:
            # 規格化した座標の応答に、基準運動量の比の平方根が掛かる)。加速がなければ 1
            from .trptacc import gammab_profile
            from .emit6d import ELMASS as _EM
            gb = gammab_profile(self.els, self.momentum() / _EM)
            for i_ in range(len(stl)):
                js_ = int(sttl[i_]) + 1               # pgrmat1: ループの後の js は stt(i)+1 (補正磁石の終端の次の位置)
                if 1 <= js_ <= len(gb) and 1 <= im <= len(gb):
                    ar[i_] = ar[i_] * math.sqrt(gb[js_ - 1] / gb[im - 1])
        return float(ar[0]) if scalar else [float(v) for v in ar]

    def set_line(self, args, v):
        key = args[0].upper()
        locs = self._resolve_loc(args[1])
        if not locs:
            # SAD (tftwiss.f): 位置が解決できない LINE[...]= は General::invset で文を打ち切る (<work-perf>/err_test.sad)
            raise SADError('General::invset:  Invalid set to LINE in LINE["%s","%s"]' % (args[0], args[1]))
        vals = v if isinstance(v, list) else [v] * len(locs)
        for (k, fr), val in zip(locs, vals):
            if k < self.nlat:
                el = self.els[k - 1]
                if isinstance(args[1], str) and '.' not in args[1] and not any(c in args[1] for c in '*%{|'):
                    self.set_family(el.name, key, float(val))      # LINE["K1","QF"]=v はファミリー全体
                else:
                    self.set_value(el, key, float(val))
                    # LINE["K0","BPI6.3"]=v は読み出しではその 1 台だけ変わるが、次の CALC でファミリー全体が
                    # 代表 (.1) の値に揃えられる (SAD tfsetv/tffsadjust。bench_level23: .3 だけ設定すると CALC 後に 0 に戻る)
                    self._spread_pending.add((el.name, key))
        self.invalidate()

    def spread_pending(self, fams=None):
        """CALC の前に、LINE[key,"NAME.n"]= で設定されたファミリーの値を代表 (最初の 1 台) の値で全台に揃える。
        SAD の tftouch/tffsadjust と同じく、触ったキーは USE まで残り、CALC のたびに揃える
        (後から DELX/DELY 等で非代表に入れた誤差や tfadjst の共有値は、次の CALC で代表の値に戻る: bench_level31a の R11)"""
        pend = list(self._spread_pending)
        for fam, key in pend:
            if fams is not None and fam not in fams:
                continue
            inst = self.instances(fam)
            if len(inst) < 2:
                continue
            v = self.get_value(inst[0], key)
            for el in inst[1:]:
                self.set_value(el, key, v)

    def fn_element(self, args):
        if not args:
            return Null
        key = args[0]
        if isinstance(key, Sym):
            key = key.name
        if isinstance(key, str) and key.upper() in ('EXPAND', 'UPDATE'):
            return Null
        if len(args) >= 2 and isinstance(args[1], (str, Sym)):
            fam = args[1].name if isinstance(args[1], Sym) else args[1]       # 裸のシンボルも名前 (tfgetstrns)
            i = self._family(fam)
            if i is None:
                raw = self.env._elements.get(fam.upper(), {})
                return float(raw.get(key.upper(), 0.0) or 0.0)
            return self.get_value(self.els[i], key)
        return Null

    def set_element(self, args, v):
        fam = args[1]
        self.set_family(fam.upper(), args[0].upper(), float(v))

    def fn_ffs(self, args):
        cmd = args[0]
        if not isinstance(cmd, str):
            return Null
        if cmd.strip().upper() == 'SHOW':
            return self.ffs_show_list()                 # FFS["SHOW"] = FFS$SHOW[] (tfshow.f ret=.true.)
        if cmd.strip().upper() in ('VAR', 'VARS', 'VARIABLES'):
            return self.ffs_var_list()                  # FFS["VAR"] = tfvars (ret=.true.)
        unit = args[1] if len(args) > 1 else 0.0        # SAD: 出力ユニット省略時は 0 (表示しない)
        saved = self.out
        target = saved
        if is_num(unit) and int(unit) == 0:
            target = _NullOut()
        elif is_num(unit) and int(unit) != 6 and int(unit) in self.I.files:
            target = self.I.files[int(unit)]
        rec = _RecOut(target)
        self.out = rec
        self._measure_ret = None                        # v2 の MEASURE が入れる戻り値 (trackb.f ret: {np, 重心, 分散行列})
        try:
            self.run_text(cmd)
        finally:
            self.out = saved
        # tffsa.f (FFS[] の経路): コマンドがエラー (?...、???-FFS-Error-) を出したら、出力ユニットに依らずエラー行を表示して
        # FFS::error を投げ、呼び出し元の文 (Do など) を打ち切る (<work-perf>/err_test.sad)
        # 打ち切るのは SAD が ' ???-FFS-Error-' で出すもの (?No variable. / ?Missing value / ?Missing number / ?Missing amount of error)。
        # DTHETA/disp/fit の未定義パターンや ???-Emittance[] の表示は打ち切らない (<work-perf>/err_test2.sad)
        errs = [ln for ln in rec.text.split(chr(10)) if '???-FFS-Error' in ln
                or ln.lstrip().startswith(('?No variable', '?Missing value', '?Missing number', '?Missing amount'))]
        if errs:
            if isinstance(target, _NullOut):
                for ln in errs:
                    saved.write(ln + chr(10))
            cmds = [c for c in cmd.split(';') if c.strip()]
            raise SADError('FFS::error:  Error in FFS command %d in FFS["%s"]' % (max(1, len(cmds)), cmd))
        if getattr(self, '_measure_ret', None) is not None:
            ret_, self._measure_ret = self._measure_ret, None
            return ret_                                 # FFS["MEASURE …"] (v2 だけ。k64 では _measure_ret は常に None)
        return Null


    def ffs_var_list(self):
        """FFS["VAR"] (tfvars.f): 変数ごとに {名前, キー, Now, Previous, Saved(要素の現在値), Minimum, Maximum, 結合先の要素名, 係数}"""
        m = self.matcher_()
        rows = []
        for v in m.variables:
            el0 = v['els'][0]
            val = float(self.get_value(el0, v['key']))
            lo, hi = v['vlim'][0], v['vlim'][1]
            rows.append([v.get('label', v['family']), v['key'], float(v.get('val', val)), float(v.get('val', val)), val, float(lo), float(hi),
                         self.elname_of(el0), 1.0])
        return rows

    def elname_of(self, el):
        """要素オブジェクトの '.n' 付き名前 (最初の 1 台なら NAME.1)"""
        for k in range(1, self.nlat):
            if self.els[k - 1] is el:
                return self.elname(k)
        return el.name

    def ffs_show_list(self):
        """FFS$SHOW[] / FFS["SHOW"] (tfinitn.f 1025): 現在の fit 条件のリスト。各要素は SAD の出力どおり
        {comp1, comp2 ("" なら 1 点), fun (maxfit なら名前+'M'), goal, 運動量点の数 |mfitp|, scale (NX/NY/NZ は 2π、CHI は π/180、他 1)}"""
        m = self.matcher_()
        conds = []
        for c in m.conditions:
            if c['mfitp'] == 0:
                continue
            kp, kp1, kf = c['kp'], c['kp1'], c['kf']
            n1 = self.elname(kp)
            n2 = self.elname(kp1) if kp != kp1 else ''
            fun = kf + ('M' if c['mfitp'] < 0 else '')
            sc = PI2 if kf in ('NX', 'NY', 'NZ') else (math.pi / 180.0 if kf.startswith('CHI') else 1.0)
            goal = float(c['fmin']) / sc if c['fmin'] == c['fmax'] else [float(c['fmin']) / sc, float(c['fmax']) / sc]
            conds.append([n1, n2, fun, goal, float(abs(c['mfitp']) - 1), sc])       # pysad の mfitp は SAD の npoints+1
        return conds

    def _definition_of(self, name):
        """?名前 = Definition[名前]: Hold[{名前:=名前, 定義 …}]。K64 は定義を内部のハッシュ表の順に返す (再現できない) ので、pysad は
        書いた順 (記号の値、具体的な添字の定義、パターンの定義の順) で返す。Plots.n の Delete[?contl, {{1,1},{1,-1}}] は先頭の `名前:=名前` と
        末尾のパターンの定義 (contl[_][_]) を除く"""
        from .sadscript import _Lit
        I = self.I
        dv = [d for d in I.downvalues.get(name, []) if d[0] is not None]
        sv = list(I.subvalues.get(name, []))
        if not dv and not sv and name not in I.globals:
            return None

        def un(p):
            return p.v if isinstance(p, _Lit) else p

        lit, pat = [], []
        for pats, body, delayed in dv:
            lhs = Expr(S(name), [un(p) for p in pats])
            (pat if any(_is_pattern(un(p)) for p in pats) else lit).append(Expr(S('SetDelayed'), [lhs, body]))
        for hp, ap, body, delayed in sv:
            lhs = Expr(Expr(S(name), [un(p) for p in hp]), [un(p) for p in ap])
            (pat if any(_is_pattern(un(p)) for p in list(hp) + list(ap)) else lit).append(Expr(S('SetDelayed'), [lhs, body]))
        first = [Expr(S('SetDelayed'), [S(name), S(name)])]
        return Expr(S('Hold'), [first + lit[::-1] + pat[::-1]])

    def flag_query(self, name) -> float:
        """?FLAG: フラグの値 (True → 1.0、False → 0.0)。NOxxx と sino (RING/INS/UNIFORM …) も受ける"""
        w = str(name).upper()
        neg = False
        if w in _FLAG_ALT:
            f, val = _FLAG_ALT[w]
            return 1.0 if bool(self.flags.get(f, False)) == val else 0.0
        if w.startswith('NO') and w[2:] in self.flags:
            w, neg = w[2:], True
        if w not in self.flags:
            # フラグでなければ `?名前` は定義の問い合わせ (Information)。定義 1 つのリストを返す (KBFrame.n の ?ReturnToSAD)
            d = self._definition_of(str(name))
            return d if d is not None else [Expr(S('Hold'), [Null])]
        v = bool(self.flags[w])
        return 1.0 if (v != neg) else 0.0

    def _mass_charge(self):
        """tfsetparam.f: charge = CHARGE、amass = MASS (eV)。既定 (電子) のときは何も渡さない (従来と同じ)"""
        g = self.I.globals
        out = {}
        try:
            am = float(g.get('MASS', 0.51099895e6))
            ch = float(g.get('CHARGE', 1.0))
        except (TypeError, ValueError):
            return out
        if am != 0.51099895e6 and am > 0.0:
            out['amass'] = am
        if ch != 1.0:
            out['charge'] = ch
        return out

    def _set_omega0(self):
        """tphyzp: OMEGA0 = 2π c p0/h0 / C (C = 全長)"""
        if self.lat is None:
            return
        circ = sum(el.L for el in self.els)
        if circ <= 0.0:
            return
        amass = float(self.I.get('MASS', 0.51099895e6) or 0.51099895e6)
        p0 = self.momentum(warn=False) / amass
        self.I.globals['OMEGA0'] = 2.0 * math.pi * 299792458.0 * p0 / math.sqrt(1.0 + p0 * p0) / circ

    def run_text(self, text: str):
        """FFS["..."] / Get の中の FFS コマンド列"""
        return self.I.run(text, 'ffs')

    def fn_seedrandom(self, args):
        """SeedRandom[n] は種 (n を uint32 に切り捨て)、SeedRandom["SAD"] は乱数の種類の選択 (状態は変えない。"SAD" 以外は System::error)、
        SeedRandom[] は状態 {"SAD", seed, filled, next の下位語, 上位語} を返す。tfRandom_.c の癖 (gs で確認): リストは先頭が文字列なら種類の選択だけ
        (状態は変わらず、続く数は読まない)、数のリストは要素の値を未初期化のまま詰める (gs では 0 になり種 1)"""
        import struct
        if len(args) > 1:
            raise SADError('General::narg:  Number of arguments is expected to be 0 or 1 in SeedRandom[%s]' % ','.join(self.I.to_str(x, None, True) for x in args))
        if args:
            x = args[0]
            if is_num(x):
                self.rng.restore(int(x) & 0xFFFFFFFF if int(x) >= 0 else int(x))
            elif isinstance(x, str) or (isinstance(x, list) and x and isinstance(x[0], str)):
                name = x if isinstance(x, str) else x[0]
                if name != 'SAD':
                    raise SADError('System::error:  No such random plugin in SeedRandom[%s]' % self.I.to_str(x, None, True))
            elif isinstance(x, list):
                nums = [v for v in x if is_num(v)]
                if any(not is_num(v) and not isinstance(v, str) for v in x):
                    raise SADError('General::wrongtype:  Argument must be Real in SeedRandom[%s]' % self.I.to_str(x, None, True))
                if nums:
                    self.rng.restore(0)
            elif x is not Null:
                raise SADError('General::wrongtype:  Argument must be String or Real or List or Null in SeedRandom[%s]' % self.I.to_str(x, None, True))
        lo, hi = struct.unpack('<II', struct.pack('<d', float(self.rng.next)))
        return ['SAD', float(self.rng.seed % (1 << 32)), float(1 if self.rng.filled else 0), float(lo), float(hi)]

    # ------------------------------------------------------------------ Emittance[] / TrackParticles[]
    def fn_emittance(self, args):
        """Emittance[opt] (Packages/emit.n): 規則のリストを返す。mode は Matrix→3, Orbit→2, OneTurnInformation→1, 既定 0
        (emit_extra.emittance_rules)。SaveEMIT->False (既定) のとき EMITX/EMITY/EMITZ/SIGE/SIGZ は Emittance[] の前の値に戻す
        (emit.n: `If[~(SaveEMIT/.op),EMITX=ex0;...]`)。EnergyLossU0 は eV (SAD で確認)。INTRA ON なら TouschekTable も作る"""
        from .emit_extra import emittance_rules, rules_to_sad_list
        if self.lat is None:
            raise SADError('No beam line')
        self.spread_pending()
        opts = {}
        for a in args:
            for r in (a if isinstance(a, list) else [a]):
                if isinstance(r, Rule) and isinstance(r.lhs, Sym):
                    opts[r.lhs.name] = r.rhs

        def flag(name):
            v = opts.get(name, 0.0)
            return bool(is_num(v) and float(v) != 0.0)
        matrix, orbit, onet = flag('Matrix'), flag('Orbit'), flag('OneTurnInformation')
        mode = 3 if matrix else (2 if orbit else (1 if onet else 0))
        g = self.I.globals
        keep = {k: g.get(k) for k in ('EMITX', 'EMITY', 'EMITZ', 'SIGE', 'SIGZ')}
        self._check_sol_optics(emit=True)
        e = x = None
        try:
            if self.flags['INTRA']:
                e = self.run_emit(quiet=True)          # IBS の反復と TouschekTable (グローバル)
            _lat.RFSW = self.flags['RFSW']
            self._sync_emitdiv()
            x = emittance_rules(self._radtaper_lattice(), self.momentum(), rfsw=self.flags['RFSW'], radcod=self.flags['RADCOD'],
                                one_turn=onet, orbit=orbit, matrix=matrix, beam_matrices=False, fshift=self._fshift(),
                                pol=bool(self.flags['POL']), sc=self._sc_opts(), **self._mass_charge())
        except Exception:
            x = None
        if x is None:
            # 閉軌道が見つからない等: EMIT と同じ NaN の報告 (run_emit が ???-Emittance[]-closed orbit not found. を出す)
            e = self.run_emit(quiet=True)
            res = self._emit_rules_nan(e)
        else:
            em = x['_emit']
            self.trf0 = float(em['trf0'])
            sc_x = x.get('_sc')
            if sc_x is not None:
                g['EMITXC'], g['EMITYC'], g['EMITZC'] = float(sc_x['EMITXC']), float(sc_x['EMITYC']), float(sc_x['EMITZC'])     # tintraconv: rsetgl1('EMITXC', eemx)
            if self.flags['INTRA'] and e is not None:
                # 最終パスの tecalc の値 (IBS 込み) を Emittances などに反映
                x['Emittances'] = np.array([e['emittance_x'], e['emittance_y'], e['emittance_z']])
                x['MomentumSpread'] = e['momentum_spread']
                x['BunchLength'] = e['bunch_length']
            res = [self._rule_value(k, v) for k, v in rules_to_sad_list(x, mode)]
            if 'Emittance' in opts and not flag('Emittance') and mode == 0:
                # emit.n: mode = If[Emittance, 0, -1, 0]。Emittance->False は mode -1 で、DampingRate 以降 (放射の量) の規則を返さない
                cut = {'DampingRate', 'Emittances', 'MomentumSpread', 'BunchLength', 'TuneShiftByRadiation', 'NominalSpinTune', 'PolarizationTime',
                       'SpinTune', 'Polarization', 'Polarizations', 'PolarizationVector'}
                res = [r for r in res if not (isinstance(r, Rule) and isinstance(r.lhs, Sym) and r.lhs.name in cut)]
            if flag('SaveEMIT') and not self.flags['INTRA']:
                em3 = x['Emittances']
                g['EMITX'], g['EMITY'], g['EMITZ'] = float(em3[0]), float(em3[1]), float(em3[2])
                g['SIGE'], g['SIGZ'] = float(x['MomentumSpread']), float(x['BunchLength'])
            self._beam_cache = None
        if not flag('SaveEMIT'):
            for k, v in keep.items():
                if v is not None:
                    g[k] = v
        return res

    @staticmethod
    def _rule_value(k, v):
        """emit_extra の値 → SADScript の値 (numpy → リスト、bool → 1/0、ExtendedTwissParameters → 規則のリスト)"""
        if k == 'ExtendedTwissParameters':
            return Rule(S(k), [Rule(S(n), float(x)) for n, x in v])
        if isinstance(v, (bool, np.bool_)):
            v = 1.0 if v else 0.0
        elif isinstance(v, np.ndarray):
            v = v.tolist()
        elif isinstance(v, (tuple, list)):
            v = [x.tolist() if isinstance(x, np.ndarray) else x for x in v]
            v = [float(x) if isinstance(x, (int, float, np.floating)) else x for x in v]
        elif isinstance(v, (int, float, np.floating)):
            v = float(v)
        return Rule(S(k), v)

    def _emit_rules_nan(self, e):
        """閉軌道が見つからなかったときの規則 (NaN)。EMIT の失敗と同じ値で、既定 mode の規則だけ"""
        nan = float('nan')
        nlat = self.nlat
        return [Rule(S('Stable'), 0.0), Rule(S('Region'), [1.0, float(nlat)]), Rule(S('Tunes'), list(e['tunes'])),
                Rule(S('EnergyLossU0'), e['energy_loss_eV']), Rule(S('RfVoltageVc'), e['rf_voltage_eff']),
                Rule(S('EquilibriumPosition'), e['trf0']), Rule(S('MomentumCompaction'), e['momentum_compaction']),
                Rule(S('OrbitDilation'), 0.0), Rule(S('BucketHeight'), e['bucket_height']), Rule(S('HarmonicNumber'), nan),
                Rule(S('OrbitAtExit'), list(map(float, e['closed_orbit']))),
                Rule(S('DampingRate'), list(e['damping_rate'])),
                Rule(S('Emittances'), [e['emittance_x'], e['emittance_y'], e['emittance_z']]),
                Rule(S('MomentumSpread'), e['momentum_spread']), Rule(S('BunchLength'), e['bunch_length'])]

    def beam_size_value(self, suffix: str, k: int, fr: float = 0.0):
        """Twiss["SIG"//suffix, l] / LINE["SIG"//suffix, l] の値 (tftwiss.f:982-1005 + tfbeamkey): suffix 空 → 6×6 のリスト、
        'X' 'PX' 'Y' 'PY' 'Z' 'DP' → sqrt(Σ_ii)、'XY' 等 → Σ_ij。ビーム行列は beam_matrix_twiss (CALC6D なら beam_matrix_at)"""
        from .emit_extra import beam_size_at
        kk = min(k, self.nlat)
        if self.flags.get('CALC6D'):
            if fr:
                raise SADError('LINE["SIG..."] at a fractional position is not supported with CALC6D')
            B = self.beam_matrix_at(kk)
        else:
            B = self.beam_matrix_twiss(kk, fr)
        try:
            v = beam_size_at({'BeamMatrices': {kk - 1: B}}, kk, suffix)
        except KeyError:
            raise SADError('General::wrongval:  "coordinate" is must be to be "X, PX, Y, PY, Z, DP"')
        return v.tolist() if hasattr(v, 'tolist') else float(v)

    def beam_matrix_twiss(self, k: int, fr: float = 0.0) -> np.ndarray:
        """tfbeam (tfbeam.f、非 CALC6D) の直訳: 位置 k (+fr) の物理座標のビーム行列 Σ (6×6)。
        pysad の仕様 (11/05 ユーザー指示): SAD はグローバル EMITX/EMITY/SIGE をそのまま使うが、ここでは **EMIT を先に実行していること** を
        条件にし、EMIT が決めたエミッタンス・エネルギー広がりと現在の Twiss からビームサイズを作る。
        基準座標 (ノーマルモード) のビーム行列と物理座標への射影 (qtent の R1-R4 の変換、MINCOUP の下限、gauss の pfi=SIGE²) は SAD のまま。
        加速 (gammab の比 r) は 1 とする。"""
        if self.emit is None:
            raise SADError('ビームサイズ (SIG..) には EMIT の事前実行が必要です (pysad の仕様: EMITX/EMITY/SIGE は EMIT の結果を使う)')
        r = self.row_frac(k, fr) if fr else self.row(k)
        g = self.I.globals
        emx = float(g.get('EMITX', 0.0) or 0.0)
        emy = float(g.get('EMITY', 0.0) or 0.0)
        sige = float(g.get('SIGE', 0.0) or 0.0)
        mincoup = float(g.get('MINCOUP', 1e-2) or 0.0)
        pfi = sige**2                                   # gauss=.true. (EMIT 後)
        emmin = (emx + emy) * mincoup
        emx = max(emmin, emx)
        emy = max(emmin, emy)
        ax, bx, ay, by = r.alfx, r.betx, r.alfy, r.bety
        ex, epx, ey, epy = r.etx, r.etpx, r.ety, r.etpy
        B = np.zeros((6, 6))
        B[0, 0] = bx * emx + ex**2 * pfi
        B[0, 1] = B[1, 0] = -ax * emx + ex * epx * pfi
        B[0, 2] = B[2, 0] = ex * ey * pfi
        B[0, 3] = B[3, 0] = ex * epy * pfi
        B[1, 1] = (1.0 + ax**2) / bx * emx + epx**2 * pfi
        B[1, 2] = B[2, 1] = epx * ey * pfi
        B[1, 3] = B[3, 1] = epx * epy * pfi
        B[2, 2] = by * emy + ey**2 * pfi
        B[2, 3] = B[3, 2] = -ay * emy + ey * epy * pfi
        B[3, 3] = (1.0 + ay**2) / by * emy + epy**2 * pfi
        B[0, 5] = B[5, 0] = ex * pfi
        B[1, 5] = B[5, 1] = epx * pfi
        B[2, 5] = B[5, 2] = ey * pfi
        B[3, 5] = B[5, 3] = epy * pfi
        B[5, 5] = pfi
        # qtent (disp=.false.): 基準座標 → 物理座標の 4×4 (R1-R4)。theta=0 なので qchg は恒等
        r1, r2, r3, r4 = r.r1, r.r2, r.r3, r.r4
        detr = r1 * r4 - r2 * r3
        cc = math.sqrt(max(0.0, 1.0 - detr))
        T = np.eye(6)
        if r.detr < 1.0:
            T[0:4, 0:4] = [[cc, 0.0, r4, -r2], [0.0, cc, -r3, r1], [-r1, -r2, cc, 0.0], [-r3, -r4, 0.0, cc]]
        else:
            T[0:4, 0:4] = [[r4, -r2, cc, 0.0], [-r3, r1, 0.0, cc], [cc, 0.0, -r1, -r2], [0.0, cc, -r3, -r4]]
        return T @ B @ T.T                              # tmulbs: Σ' = T Σ Tᵀ

    def _sc_opts(self):
        """WSPAC ON のとき Emittance の反復 (temit.f の calint1 = wspac or intra) に渡す dict。OFF なら None"""
        if not self.flags['WSPAC']:
            return None
        g = self.I.globals
        return dict(pbunch=float(g.get('PBUNCH', 0.0) or 0.0), mincoup=float(g.get('MINCOUP', 1e-2) or 0.0), intra=bool(self.flags['INTRA']))

    def beam_matrix_at(self, k):
        """LINE["SIG..", k] 用: 要素 k の入口の平衡ビーム行列 Σ (6x6)。emit_extra の BeamMatrices (temit.f:214-222) を遅延計算してキャッシュ"""
        key = (self.flags['RFSW'], self.flags['RADCOD'], id(self.lat))
        cache = getattr(self, '_beam_cache', None)
        if cache is None or cache[0] != key:
            from .emit_extra import emittance_rules
            _lat.RFSW = self.flags['RFSW']
            self._sync_emitdiv()
            x = emittance_rules(self._radtaper_lattice(), self.momentum(), rfsw=self.flags['RFSW'], radcod=self.flags['RADCOD'],
                                one_turn=False, orbit=False, matrix=False, beam_matrices=True, fshift=self._fshift(), **self._mass_charge())
            cache = (key, x['BeamMatrices'])
            self._beam_cache = cache
        return cache[1][k - 1]

    def _sync_emitdiv(self):
        """大域変数 EMITDIV・EMITDIVB・EMITDIVQ (放射の分割数の倍率。既定 1) を emit6d に渡す (K64 の tfsetparam)"""
        from . import emit6d as _e6
        for i, k in enumerate(('EMITDIV', 'EMITDIVB', 'EMITDIVQ')):
            try:
                v = float(self.I.get(k, 1.0))          # MAIN の `EMITDIV = 20;` も FFS の `EMITDIV=20;` も同じ変数
            except (TypeError, ValueError):
                v = 1.0
            _e6.EMITDIV[i] = v if (v == v and v > 0.0) else 1.0

    def _fshift(self):
        """グローバル FSHIFT (Orbit dilation)。emit の閉軌道・報告で RFSW ON のとき dv に足される"""
        return float(self.I.globals.get('FSHIFT', 0.0) or 0.0)

    def _radtaper_on(self):
        f = self.flags
        return bool(f['RAD'] and f['RADCOD'] and f['RADTAPER'] and f['RFSW'])

    def _radtaper_lattice(self):
        """RADTAPER (radtaper.py): RAD・RADCOD・RADTAPER・RFSW が ON なら、閉軌道の局所の運動量に磁石の強さを合わせた要素の列の
        Lattice を返す (Emittance[]・EMIT 用)。Twiss の DDP (閉軌道の δ) と、追跡用の状態 self._taper を作る。そうでなければ self.lat"""
        self._taper = None
        if not self._radtaper_on():
            return self.lat
        import copy as _copy
        from . import radtaper as _rt
        d0 = self.dp0()
        if _compat.v2() and d0 == getattr(self, '_dp0_from_c6d', None):
            d0 = 0.0          # v2: CALC6D が DP0 (= LINE["DDP",1] の別名) に書いた閉軌道の δ は、テーパーの閉軌道の探索の δ にしない
                              # (k64 は gs と同じに DP0 を使う: RADCOD の CALC6D のあとの RADTAPER の CALC6D のチューンが gs と 1e-8)
        r = _rt.solve(self.els, self.momentum(), dp0=d0, fshift=self._fshift())
        if not r['found']:
            self.write(' ???-RADTAPER-closed orbit with taper not found; taper is not applied.' + chr(10))
            return self.lat
        self._taper = dict(ddp=r['ddp'], shape=r['shape'], u=r['u'], dptaper=r['dptaper'], key=(id(self.els), len(self.els)))
        if self.rows is not None and len(self.rows) == len(r['ddp']):
            for row, d in zip(self.rows, r['ddp']):
                row.delta_cod = float(d)              # tturne.f: twiss(l, mfitddp) = cod(6)
        lat = _copy.copy(self.lat)
        lat.elements = r['els']
        lat._twiss_cache = None
        lat._z_into_trf0 = True                   # 閉軌道の z を trf0 に移す (emit6d._z_into_trf0)
        return lat

    def run_emit(self, quiet=False):
        self._check_sol_optics(emit=True)
        from .emit6d import emittance6d
        if self.lat is None:
            raise SADError('No beam line')
        _lat.RFSW = self.flags['RFSW']
        self._sync_emitdiv()
        try:
            gg = self.I.globals
            floor = tuple(float(gg.get(k, 0.0) or 0.0) for k in ('EMITXE', 'EMITYE', 'EMITZE'))
            lat_ = self._radtaper_lattice()
            e = emittance6d(lat_, self.momentum(), rfsw=self.flags['RFSW'], radcod=self.flags['RADCOD'],
                            fshift=self._fshift(), emit_floor=floor, pol=bool(self.flags['POL']), **self._mass_charge())
        except Exception:
            # SAD: "???-Emittance[]-closed orbit not found." と表示して NaN の報告を出し、続行する
            self.write(' ???-Emittance[]-closed orbit not found.' + chr(10))
            nan = float('nan')
            e = dict(emittance_x=nan, emittance_y=nan, emittance_z=nan, momentum_spread=nan, bunch_length=nan,
                     momentum_compaction=nan, rf_voltage_eff=0.0, trf0=nan, bucket_height=0.0,
                     closed_orbit=np.full(6, nan), energy_loss_eV=nan, tunes=[nan, nan, nan],
                     damping_rate=[nan, nan, nan], cd=np.full(6, complex(nan, nan)), omegaz=nan, heff=0.0,
                     # tecalc まで来ないので Extended Twiss は初期値の ri (単位行列) から (SAD の報告: BX=BY=BZ=1 で他は 0)
                     ext_twiss=_unit_ext_twiss())
            self.emit = e
            # temit.f: 閉軌道が見つからなくても EMITX.. に NaN が入る (以後の重みは emminv 由来になる)
            g = self.I.globals
            for k in ('EMITX', 'EMITY', 'EMITZ', 'SIGE', 'SIGZ'):
                g[k] = nan
            if not quiet:
                self.print_emit(e)
            return e
        self.emit = e
        self.trf0 = e['trf0']
        g = self.I.globals
        # tturne.f:60-64 (update=.true.): DTSYNCH=trf0, PHICAV=trf0·wrfeff, EFFVCRATIO=vcalpha (RING では 1), EFFVC=vceff, EFFRFFREQ=wrfeff·c/2π
        wrf = float(e.get('wrfeff', 0.0) or 0.0)
        g['DTSYNCH'] = float(e['trf0']) if e['trf0'] == e['trf0'] else 0.0
        g['PHICAV'] = g['DTSYNCH'] * wrf
        g['EFFVCRATIO'] = 1.0
        g['EFFVC'] = float(e.get('rf_voltage_eff', 0.0) or 0.0)
        g['EFFRFFREQ'] = wrf * 299792458.0 / (2.0 * math.pi)
        self._set_omega0()
        g['EMITX'] = e['emittance_x']
        g['EMITY'] = e['emittance_y']
        g['EMITZ'] = e['emittance_z']
        g['SIGE'] = e['momentum_spread']
        g['SIGZ'] = e['bunch_length']
        if self.flags['INTRA'] or self.flags['WSPAC']:
            # temit.f: INTRA ON なら一周ごとに tintrb で IBS の拡散を足し、収束後に Touschek の損失率表を作る。
            # WSPAC ON なら tintrb の twspace のレンズ (空間電荷) を一周の写像に入れ、ビーム行列を反復して自己無撞着にする (v2next B1)。
            # 値は最終パスの tecalc (EMITX 等) と収束値 (EMITXC 等)。TouschekTable はグローバル
            from .touschek import touschek_table
            try:
                tt = touschek_table(self.lat, self.momentum(), e, float(g.get('PBUNCH', 0.0) or 0.0),
                                    nbunch=float(g.get('NBUNCH', 1.0) or 1.0), mincoup=float(g.get('MINCOUP', 1e-2)),
                                    rfsw=self.flags['RFSW'], radcod=self.flags['RADCOD'], verbose=False,
                                    intra=bool(self.flags['INTRA']), wspac=bool(self.flags['WSPAC']), **self._mass_charge())
            except Exception as ex:
                self.write(f' ???-Emittance[]-{"INTRA" if self.flags["INTRA"] else "WSPAC"} failed: {ex}' + chr(10))
                tt = None
            if tt is not None:
                self.touschek = tt
                g['EMITX'] = tt['emittance_x']
                g['EMITY'] = tt['emittance_y']
                g['EMITZ'] = tt['emittance_z']
                g['SIGE'] = tt['momentum_spread']
                g['SIGZ'] = tt['bunch_length']
                g['EMITXC'] = tt['EMITXC']
                g['EMITYC'] = tt['EMITYC']
                g['EMITZC'] = tt['EMITZC']
                g['TouschekTable'] = tt['table']
                e0 = e
                e = dict(e)
                if self.flags['WSPAC']:
                    e.update({k: v for k, v in tt['tecalc'].items() if k not in ('orbit_entrance', 'orbit_exit', 'closed_orbit', 'trf0')})
                e.update(emittance_x=tt['emittance_x'], emittance_y=tt['emittance_y'], emittance_z=tt['emittance_z'],
                         momentum_spread=tt['momentum_spread'], bunch_length=tt['bunch_length'])
                self.emit = e
                if not quiet:
                    # temit.f: 1 周目 (零電流) の報告 → Touschek の表 (7301 の分岐) → 収束後の temitprint (intpri=F) の IBS ブロック
                    self.print_emit(e0)
                    if tt.get('table'):
                        self._print_touschek(tt['table'])
                    if self.flags['WSPAC']:
                        self._print_emit_wspac(e, tt)       # 空間電荷つきの全報告 (temitprint: intpri .or. wspac)
                    else:
                        self._print_emit_intra(tt)
                return e
        if not quiet:
            self.print_emit(e)
        return e

    def _print_emit_wspac(self, e, tt):
        """temitprint (intpri=F, wspac): 'Parameters with space charge:' (INTRA も ON なら 'intrabeam scattering + space charge:') の後に、
        空間電荷つきの一周の写像の報告 (Extended Twiss から Damping partition number まで) と、9103 (PBUNCH/1e10, MINCOUP %)、9102"""
        w = self.write
        g = self.I.globals
        w(('    Parameters with intrabeam scattering + space charge:' if self.flags['INTRA'] else '    Parameters with space charge:') + '\n')
        e2 = dict(e)
        e2.update(emittance_x=tt['EMITXC'], emittance_y=tt['EMITYC'], emittance_z=tt['EMITZC'])
        self.print_emit(e2, before9102=lambda: w('Particles/bunch/1e10   =' + _fr(float(tt.get('pbunch', g.get('PBUNCH', 0.0)) or 0.0) / 1e10, '11.8')
                                                 + '    ' + ' Minimum coupling       =' + _fr(float(g.get('MINCOUP', 1e-2) or 0.0) * 100.0, '11.8') + ' %  \n'),
                        no_orbit=True)

    def _print_emit_intra(self, tt):
        """temitprint (intpri=F, intra): '   Parameters with intrabeam scattering:' + 9103 (PBUNCH/1e10, MINCOUP %) + 9102
        (eemx/eemy/eemz = EMITXC/EMITYC/EMITZC、sige/sigz と Beam tilt/size は収束した最終パスの tecalc から)"""
        w = self.write
        w('    Parameters with intrabeam scattering:\n\n')
        g = self.I.globals
        w('Particles/bunch/1e10   =' + _fr(float(tt.get('pbunch', g.get('PBUNCH', 0.0)) or 0.0) / 1e10, '11.8') + '    '
          + ' Minimum coupling       =' + _fr(float(g.get('MINCOUP', 1e-2) or 0.0) * 100.0, '11.8') + ' %  \n')
        tec = tt.get('tecalc') or {}
        e2 = dict(tec)
        e2.update(emittance_x=tt['EMITXC'], emittance_y=tt['EMITYC'], emittance_z=tt['EMITZC'],
                  momentum_spread=tt['momentum_spread'], bunch_length=tt['bunch_length'])
        self._print_emit_9102(e2)

    def _print_touschek(self, table):
        """temit.f 9104/9105/9131 (pri かつ caltouck): INTRA ON の EMIT 報告の末尾。運動量アクセプタンスごとの Touschek 寿命
        (1/(touckl·tf) を autofg '9.6' の先頭 9 文字で) と、振幅 (Nx, Nz) / (Ny, Nz) ごとの寿命/100s の表 (nint(.01/touckm) を 0..999 に抑える)"""
        w = self.write
        dp, rate = table[1]
        n = len(rate)
        w('\n')
        for iii in range((n - 1) // 5 + 1):
            w(' Momentum acceptance:  ' + ''.join(f'{(5 * iii + i + 1) * 0.2:8.1f}  ' for i in range(1, 6)) + '  %\n')
            cells = []
            for i in range(1, 6):
                tt = float(rate[5 * iii + i - 1]) if 5 * iii + i <= n else 0.0
                cells.append(_fr(1.0 / tt, '9.6')[:9] if tt != 0.0 else '   ---')
            w(' Touschek lifetime:    ' + ''.join(c + ' ' for c in cells) + ' sec\n')
        w('\n')
        tampl, mx, my = table[2]
        head = ('  Nz: 4.6 5.3   6 6.8 7.7 8.8  10  11  13  15  17  19  21  24  28  32  36  41  46  53'
                '  60  68  77  88 100\n')
        bar = '!---1---1---1---1---1---1---1---1---1---1---1---1---1---1---1---1---1---1---1---1---1---1---1---1---1\n'

        def cell(v):
            if v == 0.0 or v != v:
                return 999
            return max(0, min(999, int(math.floor(0.01 / v + 0.5))))

        for label, col, m in (('Nx', 0, mx), ('Ny', 1, my)):
            w(" Touschek lifetime/100s for aperture 2J%s/(N%s**2 emitx') + 2Jz/(Nz**2 emitz) < 1:\n" % (label[1], label[1]))
            w(head)
            w('  ' + label + bar)
            for k in range(len(m)):
                w(f'{int(tampl[col][k]):4d}!' + ''.join(f'{cell(float(v)):4d}' for v in m[k]) + '\n')

    def print_emit(self, e, before9102=None, no_orbit=False):
        p0 = self.momentum()
        circ = self.lat.circumference
        # temitprint: f0 = omega0/2π = c β / C (β = p/h)。以前は β を落としていて SAD と 4e-8 ずれ、桁も丸めていた
        pn = p0 / 0.51099895e6
        f0 = 2.99792458e8 * (pn / math.sqrt(1.0 + pn * pn)) / circ
        u0 = e['energy_loss_eV']
        tun = e['tunes']
        dmp = e['damping_rate']
        w = self.write
        c_in = e.get('orbit_entrance')
        c_out = e.get('orbit_exit', e['closed_orbit'])
        if not no_orbit:
            w('\n    Closed orbit:\n')
            w('                  x     px/p0         y     py/p0         z     dp/p0 \n')
            # tput.f: 各値は ' '//autofg(v,'9.6') (10 桁)。Entrance は tcod の cod、Exit は最後の一周の出口 (temit.f:960-962)
            if c_in is not None:
                w('  Entrance :' + ''.join(' ' + _fr(float(v), '9.6') for v in c_in) + '\n')
            w('      Exit :' + ''.join(' ' + _fr(float(v), '9.6') for v in c_out) + '\n\n')
        else:
            w('\n')
        # temitprint: rxrxi = rx·tinv6(rx) - 1 の絶対値の和が somin=1e-8 を超えたら偏差を出す (不安定な光学など T が悪条件のとき)。
        # list-directed の real*8 (17 桁・指数 3 桁) で出る
        rxm = e.get('rx')
        if rxm is not None:
            from .emit6d import _tinv6
            with np.errstate(all='ignore'):
                so = float(np.sum(np.abs(rxm @ _tinv6(rxm) - np.eye(6))))
            if so > 1e-8:
                w(' *** Deviation from symplectic matrix = ' + '   ' + f'{so:.16E}'.replace('E-', 'E-0').replace('E+', 'E+0') + chr(10))
        self._print_emit_exttwiss(e)
        w(f'Design momentum      P0 ={_fr(p0 / 1e9, "10.7")} GeV Revolution freq.     f0 ={_fr(f0, "10.7")} Hz \n')
        w(f'Energy loss per turn U0 ={_fr(u0 / 1e6, "10.7")} MV  Effective voltage    Vc ={_fr(e["rf_voltage_eff"] / 1e6, "10.7")} MV \n')
        w(f'Equilibrium position dz ={_fr(e["trf0"] * 1e3, "10.7")} mm  Momentum compact. alpha ={_fr(e["momentum_compaction"], "10.5")}\n')
        # temitprint 9101: Orbit dilation dl = -dleng (dleng = FSHIFT·circ, tfgeo.f:134)、h = wrfeff c/ω0、fs = ωz/2π
        nan = float('nan')
        dleng = float(self.I.globals.get('FSHIFT', 0.0) or 0.0) * circ
        w(f'Orbit dilation       dl ={_fr(0.0 - dleng * 1e3, "10.7")} mm  Effective harmonic #  h ={_fr(e.get("heff", nan), "10.7")}\n')
        w(f'Bucket height     dV/P0 ={_fr(e["bucket_height"], "10.7")}     Synchrotron frequency   ={_fr(e.get("omegaz", nan) / (2.0 * math.pi), "10.7")} Hz \n\n')
        # 9012: 3(f10.7,10x)。行末の 10x は出力されない。Imag.tune = Re log λ / 2π (安定なら 0 の雑音)
        cd = e.get('cd')
        imt = [float(cd[3 + k].real) / (2.0 * math.pi) for k in range(3)] if cd is not None else [nan] * 3
        w('  Imag.tune:' + (' ' * 10).join(_f107(v) for v in imt) + '\n')
        w('  Real tune:' + (' ' * 10).join(_f107(float(v)) for v in tun) + '\n\n')
        w('    Damping per one revolution:\n')
        w(_line9013([float(v) for v in dmp]))
        # 9013: 減衰時間 = -2π/ω0/Re cd、放射によるチューンシフト = Im cd/2π、分配数 J = Re cd · 4/ΣRe cd (tecalc)
        om0 = e.get('omega0') or 2.0 * math.pi * f0
        if cd is not None and om0:
            rd3 = [float(cd[k].real) for k in range(3)]
            w('    Damping time (sec):\n')
            w(_line9013([(-2.0 * math.pi / om0 / v) if v != 0.0 else nan for v in rd3]))
            w('    Tune shift due to radiation:\n')
            w(_line9013([float(cd[k].imag) / (2.0 * math.pi) for k in range(3)]))
            sr = sum(rd3)
            jj = [v * 4.0 / sr for v in rd3] if sr != 0.0 else [0.0, 0.0, 0.0]
            w('    Damping partition number:\n')
            w('          X :' + _ff(jj[0], 10, 4) + ' ' * 7 + 'Y :' + _ff(jj[1], 10, 4) + ' ' * 7 + 'Z :' + _ff(jj[2], 10, 4) + '\n')
        w('\n')
        if before9102 is not None:
            before9102()
        self._print_emit_9102(e)

    def _print_emit_9102(self, e):
        """temitprint 9102: Emittance X/Y/Z, Energy spread, Bunch Length, Beam tilt, Beam size xi/eta, Nominal spin tune, Polarization time"""
        w = self.write
        pn = self.momentum() / 0.51099895e6
        w(f'Emittance X            ={_fr(e["emittance_x"], "11.8")} m   Emittance Y            ={_fr(e["emittance_y"], "11.8")} m\n')
        w(f'Emittance Z            ={_fr(e["emittance_z"], "11.8")} m   Energy spread          ={_fr(e["momentum_spread"], "11.8")}\n')
        # temitprint (temit.f:1126-1141, 9102): Beam tilt と主軸方向の Beam size は入口の平衡ビーム行列から。
        # 閉軌道が見つからなかった (beam_matrix なし) ときは NaN
        if e.get('beam_matrix') is not None:
            from .emit6d import beam_tilt_sizes
            btilt, sigx, sigy = beam_tilt_sizes(np.asarray(e['beam_matrix']))
        else:
            btilt = 0.0                     # SAD は NaN でも Beam tilt を 0 と出す (btilt=0 の初期値のまま)
            sigx = sigy = float('nan')
        w(f'Bunch Length           ={_fr(e["bunch_length"] * 1e3, "11.8")} mm  Beam tilt              ={_fr(btilt, "11.8")} rad\n')
        w(f'Beam size xi           ={_fr(sigx * 1e3, "11.8")} mm  Beam size eta          ={_fr(sigy * 1e3, "11.8")} mm\n')
        # 9102 の最後の行: Nominal spin tune = h0·GSPIN (params(ipnnup))、Polarization time (POL OFF は 0: calpol でなければ params=0)
        from .emit_extra import GSPIN
        h0 = e.get('h0') or math.sqrt(1.0 + pn * pn)
        nnu = h0 * GSPIN if h0 is not None else float('nan')
        sp = e.get('spin')
        taup = sp['taup'] / 60.0 if sp is not None else 0.0
        w(f'Nominal spin tune      ={_fr(nnu, "11.7")}     Polarization time      ={_fr(taup, "11.7")} min\n\n')
        if sp is not None:
            self._print_emit_spin(e, sp)

    def _print_emit_spin(self, e, sp):
        """temitprint の calpol 部 (temit.f:1188-1225)。1p3g15.7 / 1p6g15.7 は Fortran の G 編集 (_g157)"""
        w = self.write
        w('   Polarization vector at entrance:\n')
        w(_line9013(sp['sps'][:, 0], ('x :', 'y :', 'z :')))
        if self.I.globals.get('EMIOUT') or self.flags.get('EMIOUT'):
            w(' \n   Spin precession matrix:\n')
            w('                 sx             sy             sz\n')
            for lab, row in zip(('sx', 'sy', 'sz'), sp['srot1']):
                w('        ' + lab + ''.join(_g157(v) for v in row) + '\n')
            w(' \n   One-turn depolarization vectors:\n')
            w('                 x              px             y              py             z              pz\n')
            srot = sp['srot']
            for lab, row in zip(('gx', 'gy', 'gz'), srot[:, 3:9]):
                w('        ' + lab + ''.join(_g157(v) for v in row) + '\n')
            srotn = srot[:, 3:9] @ e['normal_r'] if e.get('normal_r') is not None else srot[:, 3:9]
            w(' \n                 X              Px             Y              Py             Z              Pz\n')
            for lab, row in zip(('gx', 'gy', 'gz'), srotn):
                w('        ' + lab + ''.join(_g157(v) for v in row) + '\n')
        w('\nSpin tune              =' + _fr(sp['nup'], '11.8') + '     Equil. polarization    ='
          + _fr(sp['equpol'][0] * 100.0, '11.8') + ' %\n\n')

    def _print_emit_exttwiss(self, e):
        """temitprint 9001: Extended Twiss Parameters (params(iptws0+i) = tfetwiss(ri, cod) を autofg '9.6' で)。
        SAD の癖: EPX の欄に mfitey (EY) を渡している (引数の並びのミス)ので EPX は EY と同じ値が出る"""
        tw = e.get('ext_twiss')
        if tw is None:
            return
        from .emit_extra import TWISS_FUNS
        v = {k: _fr(float(tw[i]), '9.6') for i, k in enumerate(TWISS_FUNS)}
        sp = ' '
        w = self.write
        w('    Extended Twiss Parameters\n')
        w('AX:' + v['AX'] + ' BX:' + v['BX'] + sp * 26 + '  ZX:' + v['ZX'] + '  EX:' + v['EX'] + '\n')
        w(sp * 11 + 'PSIX:' + v['PSIX'] + sp * 26 + ' ZPX:' + v['ZPX'] + ' EPX:' + v['EY'] + '\n')
        w('R1:' + v['R1'] + ' R2:' + v['R2'] + ' AY:' + v['AY'] + ' BY:' + v['BY'] + '  ZY:' + v['ZY'] + '  EY:' + v['EY'] + '\n')
        w('R3:' + v['R3'] + ' R4:' + v['R4'] + sp * 12 + 'PSIY:' + v['PSIY'] + ' ZPY:' + v['ZPY'] + ' EPY:' + v['EPY'] + '\n')
        w(sp * 51 + '  AZ:' + v['AZ'] + '  BZ:' + v['BZ'] + '\n')
        w(sp * 65 + 'PSIZ:' + v['PSIZ'] + '\n')
        w('    Units: B(X,Y,Z), E(X,Y), R2: m | PSI(X,Y,Z): radian | ZP(X,Y), R3: 1/m\n\n')

    def fn_track(self, args):
        """TrackParticles[{ls, beam}, ld, nt, nend] (tftrack.f:40-115, 233-300)。
        ls: 始点、ld: 終点 (既定 nlat)、nt..nend: 周の番号 (mt = nend−nt+1 周、既定 1)。
        ls, ld は nlat−1 を法にして正規化 (ls∈[1,nlat]、ld∈[2,nlat])。ld<=ls なら mt+1。
        mt>1 なら ls→nlat、(mt−2) 周、1→ld。mt==1 なら ls→ld。戻りの先頭は正規化した ld (merge(ld,ls,normal))。
        RFSW OFF では tturn ごと (区間ごと) に z=0 (tturn.f:142-145)。空洞の周波数は全周長 (omega0) で決まるので区間でも全周長を渡す"""
        from .track6n import track_particles_fast as track_particles
        if not args or len(args) > 4:
            raise SADError('General::narg: Number of arguments is expected to be "1, 2, 3, or 4" in TrackParticles')
        beam = args[0]
        if not (isinstance(beam, list) and len(beam) >= 2):
            raise SADError('TrackParticles: {n, {x,px,y,py,z,dp,flag}} expected')
        if self.lat is None:
            raise SADError('No beam line (USE first)')
        nlat = self.nlat

        def loc_of(a):
            r = self._resolve_loc(a)
            if len(r) != 1:
                raise SADError('General::wrongtype: Argument must be one location in TrackParticles')
            return int(r[0][0])

        ls = loc_of(beam[0])
        ld = loc_of(args[1]) if len(args) >= 2 and args[1] is not None else nlat
        nt = nend = 1
        if len(args) >= 3 and is_num(args[2]):
            nt = nend = int(args[2])
        if len(args) >= 4 and is_num(args[3]):
            nend = int(args[3])
            if nend < nt:
                raise SADError('General::wrongval: "#4 >= #3"')
        mt = nend - nt + 1
        ls = (ls - 1) % (nlat - 1) + 1 if ls > 0 else nlat - (1 - ls) % (nlat - 1)
        ld = (ld - 2) % (nlat - 1) + 2 if ld > 1 else nlat + 1 - (2 - ld) % (nlat - 1)
        if ld <= ls:
            mt += 1
        cols = beam[1]
        if not cols or not isinstance(cols[0], list):
            raise SADError('TrackParticles: beam columns expected')
        npart = len(cols[0])
        _split = self._npara_chunks(npart, ls, ld, nlat)
        if _split is not None:
            return self._fn_track_npara(args, beam, _split)
        x = np.array([[cols[j][i] for j in range(6)] for i in range(npart)], dtype=float)
        pol = bool(self.flags['POL'])
        spin_in = None
        if pol:
            # tftrack.f: POL ON のときは 9〜11 列 ({x,px,y,py,z,dp, 角度, sy/|s|, flag}…)。スピン (sx, sy, sz) = (st·cos, sy, st·sin)、st = sqrt(1−sy²)
            if not (9 <= len(cols) <= 11):
                raise SADError('General::wrongtype: Beam must be a list of 9 to 11 lists with POL in TrackParticles')
            spin_in = []
            for i in range(npart):
                sy_ = float(cols[7][i])
                st_ = 1.0 + (-sy_ * sy_) / (1.0 + math.sqrt(1.0 - sy_ * sy_))
                ang = float(cols[6][i])
                spin_in.append([st_ * math.cos(ang), sy_, st_ * math.sin(ang)])
            flags = list(cols[8])
        else:
            flags = list(cols[6]) if len(cols) > 6 else [1.0] * npart
        els = self.els
        # v2 (C2) AUTOPHI: 空洞の位相の基準にする設計軌道の到着時刻 (Twiss の DZ) を追跡に渡す (tturn.f: gettwiss(mfitdz, l))
        from . import track6n as _n6
        _n6.AUTOPHI_DZ = None
        if _compat.v2() and any(float(getattr(e, 'AUTOPHI', 0.0) or 0.0) != 0.0 for e in els):
            rows_ = self.rows if self.rows is not None else self.calc(quiet=True)
            _n6.AUTOPHI_DZ = [float(r.z_cod) for r in rows_]
        tp_ = getattr(self, '_taper', None)
        if tp_ is not None and self._radtaper_on() and tp_['key'] == (id(self.els), len(self.els)):
            # RADTAPER (tturn.f): 直前の Emittance[]・EMIT が作った閉軌道の δ (Twiss の DDP) と dptaper で磁石の強さを変える
            from . import radtaper as _rt
            els = _rt.tapered(self.els, tp_['shape'], tp_['u'], tracking=True)
        from .emit6d import _emit_circ
        circ = _emit_circ(els)                  # K64 は SOL (BOUND≠0) の長さを周長に入れない (omega0・空洞の位相)
        radlight = bool(self.flags['RADLIGHT'])
        kw = dict(trf0=self.trf0, rfsw=self.flags['RFSW'], rad=self.flags['RAD'] or radlight, fluc=self.flags['FLUC'],
                  fshift=self._fshift(), circ=circ)
        links6 = _n6.bend_links(els) if any(_n6._bend_type(e) for e in els) else None      # tfgeo の iprev・inext (F1 の分割)
        kw.update(self._mass_charge())          # tfsetparam: charge = CHARGE, amass = MASS (電子以外も追える)
        g = self.I.globals
        from . import beambeam as _bbm
        _bbm.PHOTON_SINK = None
        photons_on = bool(self.flags['PHOTONS'])
        photons_py = False
        if photons_on:
            if not (self.flags['RAD'] and self.flags['FLUC']):
                raise SADError('Track::photons:  PHTONS requires RAD & FULC in TrackParticles')
            # v2next B3: BEND の光子 (tbend・tbendi の tphrec) は track6 (Python の追跡) で記録する。
            # D1: QUAD・MULT・SEXT・OCT…・磁場つきドリフトの光子も track6 で記録する (BEAMBEAM と同じラティスでも、BEAMBEAM は track6 の中で通す)。
            # BEAMBEAM だけが光子を出すラティスは今までどおり高速追跡 (numba)。SOL 領域内の光子は未対応 (記録しない)
            _ph_types = ('Quadrupole', 'Multipole', 'Sextupole', 'Octupole', 'Decapole', 'Dodecapole')
            photons_py = any((type(e).__name__ == 'Bend' and float(getattr(e, 'ANGLE', 0.0) or 0.0) != 0.0 and e.L != 0.0)
                             or (type(e).__name__ in _ph_types and e.L != 0.0 and not any(type(x).__name__ in ('Coord', 'TCavity', 'Map') for x in els))
                             for e in els)
            gg = self.ensure_geo().geo
            for i, e in enumerate(els):
                if type(e).__name__ == 'BeamBeam':
                    e._photon = (i + 1, np.array(gg[i], dtype=float))      # bbinit の tsetphotongeo (l_track = この要素)
            _bbm.PHOTON_SINK = []                # tphotoninit
            _bbm.PHOTON_OUT = self.write
        try:
            if photons_py:
                raise NotImplementedError('PHOTONS (BEND): 光子は Python の追跡 (track6) で記録する')
            if pol and self.flags['RAD']:
                raise NotImplementedError('POL: スピンは Python の追跡 (track6) で追う')
            from .track6n import TrackRun
            tr = TrackRun(els, x, flags, self.momentum(), alost=float(g.get('LOSSAMPL', 1.0) or 0.0),
                          zlost=float(g.get('LOSSDZ', 100.0) or 0.0), kturn0=nt,
                          trpt=bool(self.flags['TRPT'] or radlight), **kw)
            if _compat.k64() and not getattr(self, '_tracked_once', False):   # v2 (B33): 最初の tapert でも落とす
                tr.l_track0 = 0                    # セッションの最初の TrackParticles: l_track=0 (gs で確かめた: 最初から LOSSAMPL の外の粒子は最初の tapert で落ちない)
            self._tracked_once = True
            if True:
                # LOSSMAP: l_track (最後に処理した要素の番号) はセッションを通じて残る (K64 の tmacro の変数。最初は 0)。
                # v2 は最初の値を「始点の 1 つ前の要素」(周回で数える) にする
                lt_ = getattr(self, '_l_track', None)
                if lt_ is None:
                    lt_ = 0 if not _compat_v2() else (ls - 1 if ls > 1 else nlat - 1)
                tr.lmap = [lt_]
        except NotImplementedError:
            tr = None                            # numba の表に無い要素 (SOL 等): 損失なしの track6 で追う (以前どおり)
        from . import track6n as _t6n
        if tr is not None and any(type(e).__name__ == 'Map' for e in els):
            _t6n.MAP_HOOK = lambda l, nt_, S_, n_, kid_, lost_: self._map_track(tr, npart, l, nt_, S_, n_, kid_, lost_)
        else:
            _t6n.MAP_HOOK = None
        wake_on = (self.flags['LWAKE'] or self.flags['TWAKE']) and self.flags['TRPT']      # tftrack.f:120
        if wake_on:
            table = self._wake_table()
            if table:
                if tr is None:
                    raise SADError('WakeFunction: この要素の並びは高速追跡 (TrackRun) で扱えないため未対応')
                from .wake import WakeKick, kick_factor
                fw = kick_factor(float(g.get('CHARGE', 1.0) or 1.0), float(g.get('PBUNCH', 0.0) or 0.0),
                                 float(g.get('NBUNCH', 1.0) or 1.0), npart, tr.ph.amass)
                tr.set_wake(WakeKick(tr.ph.p0, tr.ph.h0, fw, int(float(g.get('NBUNCH', 1.0) or 1.0))), table)
        if self.flags.get('SPAC') or self.flags.get('PSPAC') or self.flags.get('WSPAC'):
            # v2: 空間電荷 (spch.py: SPAC + DRIFT・MULT の RADIUS、PSPAC + SPCH 要素、WSPAC)。k64 モードでは今までどおり入れない
            from . import compat as _cmp_spch, spch as _spch     # (この関数の _compat は上の try の中の局所変数)
            _sc = _spch.make(self, tr, npart) if _cmp_spch.v2() else None
            if _sc is not None:
                tr.spch = _sc

        from . import beambeam as _bbm
        if any(type(e).__name__ == 'BeamBeam' for e in els):
            _bbm.LUM_CB = self._luminosity_hook   # bbstrhl.f: symbol LuminosityMonitor に値、LuminosityMonitor[lum] を呼ぶ
        rl = None
        if radlight:
            if tr is None:
                raise SADError('RADLIGHT: この要素の並びは高速追跡 (TrackRun) で扱えない')
            rl = _RadLight(self.ensure_geo().geo[ls - 1], tr.ph, npart, v2=_compat_v2(),
                           div=float(g.get('RADLIGHTDIV', 16.0) or 16.0))

        # track6 (Python) の追跡: APERT (tapert) で落ちた粒子は alive6 が False になる (LOSSAMPL/LOSSDZ の判定は numba 版だけ。zlost は APERT の |z| 判定)
        from .track6 import _is_apert as _is_ap
        alive6 = np.array([float(f) != 0.0 for f in flags])
        zlost6 = float(g.get('LOSSDZ', 100.0) or 0.0)
        # スピン追跡 (POL + RAD) の損失: tturn の tapert (LOSSAMPL・LOSSDZ、最初の要素・15 要素ごと・区間の終わり) を track6 の経路でも行う。
        # 落ちた要素と周は lost6 に入る (LOSSMAP)
        alost6 = float(g.get('LOSSAMPL', 1.0) or 0.0)
        lost6 = {}
        kt6 = [nt]
        lt6 = [getattr(self, '_l_track', None)]
        first6 = [False]
        if lt6[0] is None:
            # セッションの最初の追跡: K64 は l_track=0 で、最初の tapert は落とさない (B33)。v2 は始点の 1 つ前の要素として落とす
            lt6[0] = (ls - 1 if ls > 1 else nlat - 1) if _compat_v2() else 0
            first6[0] = not _compat_v2()

        def seg(x, a, b):                        # tturn(ls=a, ld=b): 要素 a..b−1 (els[a−1:b−1])
            if b <= a:
                return x
            if rl is not None:
                rl.store(tr, 0.0)                  # tturn.f:61 tlstore(…, dl=0): 区間の始めの点
                for l in range(a, b):
                    e = els[l - 1]
                    tn = type(e).__name__
                    if rl.v2:
                        # v2 (B65): BEND は円弧で進めて中も記録、曲げのある MULT は出口だけ、ほかは長さのある要素すべて (CAVI も)
                        ang_ = float(getattr(e, 'ANGLE', 0.0) or 0.0)
                        if tn == 'Bend' and ang_ != 0.0 and e.L != 0.0:
                            rl.bend(tr, e, l)
                            continue
                        tr.segment(l, l + 1)
                        if e.L != 0.0:
                            if tn == 'Multipole' and ang_ != 0.0:
                                rl.store(tr, 0.0, ang_, float(e.L) / ang_, float(getattr(e, 'ROTATION', 0.0) or 0.0))
                            else:
                                rl.store(tr, float(e.L))
                        continue
                    tr.segment(l, l + 1)
                    if tn == 'Bend' or (tn == 'Multipole' and float(getattr(e, 'ANGLE', 0.0) or 0.0) != 0.0):
                        continue                   # tturn.f:555: BEND と曲げのある MULT は記録しない (K64 のまま)
                    if _TYPECODE.get(tn, 99) < 31 and e.L != 0.0:
                        rl.store(tr, float(e.L))
                return x
            if tr is not None:
                tr.segment(a, b)
                return x
            if spin_in is not None or photons_py:
                if spin_in is not None and self.flags['TRPT'] and self.flags['RFSW'] and (any(
                        float(getattr(e_, 'VOLT', 0.0) or 0.0) != 0.0 for e_ in els[a - 1:b - 1] if type(e_).__name__ in ('Multipole', 'TCavity'))
                        or (any(float(getattr(e_, 'VOLT', 0.0) or 0.0) != 0.0 for e_ in els[a - 1:b - 1] if type(e_).__name__ == 'Cavity')
                            and (not _compat.v2() or any(type(e_).__name__ == 'Solenoid' for e_ in els[a - 1:b - 1])))):
                    # v2next B1: スピン追跡 (POL、track6 の経路) は TRPT の加速で変わる基準運動量を扱わない (RF の蹴りは入るが運動量の規格化が
                    # 入口のまま)。黙って δ の違う値を返さず、加速のある輸送路は POL なしの高速追跡 (SOL 領域内も対応) を使わせる
                    raise SADError('POL (スピン追跡) + TRPT + RFSW の加速 (VOLT≠0 の CAVI・MULT・TCAVI) は未対応 '
                                   '(基準運動量が変わる。POL なしの TrackParticles は対応)')
                from .track6 import track_particles as _tp6
                _res6 = _tp6(els[a - 1:b - 1], x, 1, self.momentum(), spins=spin_in, radpol=bool(self.flags['RADPOL']),
                             alive=alive6, zlost=zlost6, alost=alost6, l0=a - 1, kturn0=kt6[0], lt_prev=lt6[0],
                             lost_info=lost6, nodrop_first=first6[0], links=links6, links_off=a - 1,
                             photons=(_bbm.PHOTON_SINK if photons_py else None),
                             photon_geo=(self.ensure_geo() if photons_py else None), trpt=bool(self.flags['TRPT']), **kw)
                lt6[0] = b - 1
                first6[0] = False
                self._l_track = lt6[0]
                self._tracked_once = True
                if spin_in is None:
                    return _res6
                xo, sp_ = _res6
                spin_in[:] = [list(v) for v in sp_]
                return xo
            if any(_is_ap(e) for e in els[a - 1:b - 1]):
                from .track6 import track_particles as _tp6
                return _tp6(els[a - 1:b - 1], x, 1, self.momentum(), alive=alive6, zlost=zlost6, links=links6, links_off=a - 1, **kw)
            return track_particles(els[a - 1:b - 1], x, 1, self.momentum(), links=links6, links_off=a - 1, **kw)

        if mt > 1:
            x = seg(x, ls, nlat)
            if tr is not None:
                tr.next_turn()
            kt6[0] += 1
            ls = 1
            mt -= 1
        for _ in range(mt - 1):
            x = seg(x, 1, nlat)
            if tr is not None:
                tr.next_turn()
            kt6[0] += 1
        normal = True
        if ld <= ls:
            normal = False
        else:
            x = seg(x, ls, ld)
        if photons_on:
            self.I.globals['PhotonList'] = list(_bbm.PHOTON_SINK or [])         # tphotonlist
            _bbm.PHOTON_SINK = None
        alive = [fl if a else 0.0 for fl, a in zip(flags, alive6)]
        if tr is not None:
            x, alive = tr.result()               # tftrack.f ktfresetparticles: 番号順、alive (落ちた粒子は 0、座標は落ちた時点)
        res = [[float(x[i][j]) for i in range(npart)] for j in range(6)] + [list(alive)]
        lm_rows = None
        if tr is not None and getattr(tr, 'lmap', None) is not None:
            self._l_track = tr.lmap[0]
        if self.flags['LOSSMAP']:
            # tftrack.f ktfresetparticles: flag の後に 2 行 (落ちた要素の番号 l_track、落ちた周。落ちていなければ 0 と nend+1) を足す
            if tr is not None and getattr(tr, 'lmap', None) is not None:
                lm_rows = [list(r_) for r_ in tr.lossmap(nend, nlat, v2=_compat_v2())]
            elif spin_in is not None:
                # スピン追跡 (track6 の経路): lost6 = {粒子: (要素, 周)}
                lm_rows = [[float(lost6[i][0]) if i in lost6 else 0.0 for i in range(npart)],
                           [float(lost6[i][1]) if i in lost6 else float(nend + 1) for i in range(npart)]]
                self._l_track = lt6[0]
            else:
                raise SADError('LOSSMAP: この要素の並び・フラグは高速追跡 (TrackRun) で扱えず、損失の位置を記録できない')
            res = res + lm_rows
        if rl is not None:
            return [[float(ld if normal else ls), res], rl.result()]
        if spin_in is not None:
            # ktfresetparticles: 7 列目 = atan(sz, sx)、8 列目 = sy/|s|、9 列目 = flag
            ang = [math.atan2(sv[2], sv[0]) for sv in spin_in]
            syn = [sv[1] / math.sqrt(sv[0]**2 + sv[1]**2 + sv[2]**2) for sv in spin_in]
            res = res[:6] + [ang, syn, list(alive)] + (lm_rows or [])
        return [float(ld if normal else ls), res]

    def _npara_chunks(self, npart, ls, ld, nlat):
        """tftrack.f:124-195 の NPARA (fork による並列) の粒子の分け方。RAD と FLUC の追跡だけで結果 (乱数列) が変わる。
        npara = min(NPARA, 256, npart/9 + 1, ne·npart/3000 + 1) (ne = ld − ls、≤ 0 なら + nlat)。npara > 1 なら
        np1 = npart/npara ずつ、子プロセス k (k = 1..) が [ipn, ipn+np1) を乱数の種 + 2k で追い、親が残りを元の種で追う。
        戻り: [(始点, 終点, kseed), ...] (最後が親)。分割しない (npara = 1 か対象外のフラグ) なら None。
        pysad の既定は NPARA = 1 (gs の既定は計算機の設定で 4 などになり、gs の既定と比べるときは NPARA=1; か NPARA=4; を揃える)"""
        if getattr(self, '_in_npara', False):
            return None
        try:
            npara = int(float(self.I.globals.get('NPARA', 1.0) or 1.0))
        except (TypeError, ValueError):
            return None
        fl = self.flags
        if npara <= 1 or not (fl['RAD'] and fl['FLUC']) or fl['RADLIGHT'] or fl['POL'] or fl['PHOTONS'] or fl['LWAKE'] or fl['TWAKE']:
            return None
        ne = ld - ls
        if ne <= 0:
            ne += nlat
        npara = min(npara, 256, npart // 9 + 1, ne * npart // 3000 + 1)
        if npara <= 1:
            return None
        np1 = npart // npara
        ipn = 0
        kseed = 0
        chunks = []
        while ipn + np1 < npart and len(chunks) < npara:
            kseed += 2
            chunks.append((ipn, ipn + np1, kseed))
            ipn += np1
        chunks.append((ipn, npart, 0))
        return chunks

    def _fn_track_npara(self, args, beam, chunks):
        """NPARA > 1 の TrackParticles: 子 (種 + kseed) を順に追い、最後に親 (元の種) を追う。親の乱数の状態だけが残る (fork した子の消費は戻らない)"""
        from . import track6 as _t6
        rng = _t6.RNG
        s0 = rng.seed
        cols = beam[1]
        outs = []
        self._in_npara = True
        try:
            for (a, b, ks) in chunks:
                rng.restore(s0 + ks)                # SeedRandom[vn + kseed] (tfaddseed)。親 (kseed = 0) は元の種
                sub = [[c[i] for i in range(a, b)] for c in cols]
                outs.append(self.fn_track([[beam[0], sub]] + list(args[1:])))
        finally:
            self._in_npara = False
        nrow = len(outs[0][1])
        res = [[v for o in outs for v in o[1][j]] for j in range(nrow)]
        return [outs[0][0], res]

    def _wake_table(self):
        """tffssetupwake: 各要素 i について WakeFunction[Longitudinal, 要素名] と WakeFunction[Transverse, 要素名] を評価し、
        {i: (縦の表, 横の表, dx, dy)} を返す (表は (行, 2) の配列か None。数のリスト {{z, w}, ...} でなければ None)。
        LWAKE が OFF なら縦は使わず、TWAKE が OFF なら横は使わない (txwake の lwl, lwt)。
        CAVI / MULT(VOLT≠0) の要素に表があるときは SADError (tcav.f / tmulti.f の中での呼び出しは未対応)"""
        out = {}
        in_sol = False                  # tsol の領域 (BOUND≠0 の SOL から次の BOUND≠0 の SOL まで)
        for i in range(1, self.nlat):
            el0 = self.els[i - 1]
            is_sol = type(el0).__name__ == 'Solenoid' and float(getattr(el0, 'BOUND', 0.0) or 0.0) != 0.0
            region_inner = in_sol and not is_sol          # tsol の中の要素 (回転 ROTATE を正しく使う: tsol.f:628-631)
            if is_sol:
                if in_sol:
                    in_sol = False
                else:
                    in_sol = True
            nm = self.elname(i)
            tabs = []
            for tag, on in (('Longitudinal', self.flags['LWAKE']), ('Transverse', self.flags['TWAKE'])):
                try:
                    r = self.I.eval(Expr(S('WakeFunction'), [S(tag), _Lit(nm)]))
                except SADError:
                    r = None
                a = None
                if isinstance(r, list) and r and all(isinstance(q, list) and len(q) == 2 and all(is_num(v) for v in q) for q in r):
                    a = np.array(r, dtype=float)
                tabs.append(a if on else None)
            if tabs[0] is None and tabs[1] is None:
                continue
            el = self.els[i - 1]
            tn = type(el).__name__
            if tn == 'Cavity':
                out[i] = (tabs[0], tabs[1], 0.0, 0.0, 0.0)      # tcav.f は txwake に dx=dy=theta=0 を渡す
                continue
            if tn == 'Multipole' and (float(getattr(el, 'VOLT', 0.0) or 0.0) + float(getattr(el, 'DVOLT', 0.0) or 0.0)) != 0.0:
                out[i] = (tabs[0], tabs[1], 0.0, 0.0, 0.0)      # tmulti.f tmultiacc は txwake に dx=dy=theta=0 を渡す
                continue
            rot_ = float(getattr(el, 'ROTATE' if tn == 'Bend' else 'ROTATION', 0.0) or 0.0)      # SAD の ROTATE (Bend だけ Bend.ROTATE、他は ROTATION)
            if region_inner or (is_sol and not in_sol):
                # tsol.f: dx=DX、dy=DY、theta=ROTATE
                out[i] = (tabs[0], tabs[1], float(getattr(el, 'DX', 0.0) or 0.0), float(getattr(el, 'DY', 0.0) or 0.0), rot_)
            else:
                # tturn.f:263-264 の癖: dy=ROTATE の値、theta は未初期化 (0 とする)
                out[i] = (tabs[0], tabs[1], float(getattr(el, 'DX', 0.0) or 0.0), rot_, 0.0)
        return out

    def fn_calculate_optics(self, args):
        """CalculateOptics[lb, le, inicond, cell, idim (, geocal)] (tfoptics.f → qcell1)。
        区間 lb..le の Twiss を inicond (28 個: Twiss["*",loc] の並び) から計算し、
        {各位置 lb..le の 28 個の Twiss のリスト, {stabx, staby, stabz}, {tracex, tracey, tracez}} を返す。
        cell≠0: 区間を周期系として周期解 (閉軌道 + Twiss) を求める (inicond は閉軌道探索の初期値と位相 PSIX/PSIY のみ使う)。
        cell=0: inicond を位置 lb の入口条件として伝搬する。idim=2 は 4 次元 (qtwiss)、3 は 6 次元 (qcell6d)。
        計算した Twiss は大域の Twiss 配列に書き込まれる (呼出し後の Twiss[] に反映される)。
        tracex/y: 区間の転送行列の対角成分の和 (位置 lb と le の Twiss から x11+x22)。tracez は 0。
        stabz は SAD では未初期化 (呼出しごとに 0/1 が不定)。ここでは 1 とする。"""
        if len(args) not in (5, 6):
            raise SADError('General::narg: Number of arguments is expected to be "5 or 6" in CalculateOptics')
        a_lb, a_le, ini, a_cell, a_idim = args[:5]
        if not is_num(a_idim):
            raise SADError('General::wrongtype: Argument must be Real(IDIM) for #5 in CalculateOptics')
        idim = int(a_idim)
        if idim not in (2, 3):
            raise SADError('General::wrongnum: IDIM must be 2 or 3 in CalculateOptics')
        if not is_num(a_cell):
            raise SADError('General::wrongtype: Argument must be Real(CELL) for #4 in CalculateOptics')
        cell = a_cell != 0
        if not isinstance(ini, list):
            raise SADError('General::wrongtype: Argument must be List(inicond) for #3 in CalculateOptics')
        if len(ini) != len(_TWISS_ALL):
            raise SADError('General::wrongleng: Length of #3 must be 28 in CalculateOptics')
        if not all(is_num(x) for x in ini):
            raise SADError('General::wrongtype: Argument must be List of Reals for #3 in CalculateOptics')
        if self.lat is None:
            raise SADError('No beam line (USE first)')
        locs = []
        for a in (a_lb, a_le):
            r = self._resolve_loc(a)
            if len(r) != 1:
                raise SADError('General::wrongtype: Argument must be one location in CalculateOptics')
            locs.append(r[0][0])
        lb, le = locs
        if not 1 <= lb <= le <= self.nlat:
            raise SADError('General::wrongval: bad region in CalculateOptics')
        try:
            # 直前の値の変更で Twiss 配列が無効なら再計算するが、SAD の CalculateOptics は Twiss 配列を再計算しない (報告も出さない) ので黙って行う
            rows = list(self.rows if self.rows is not None else self.calc(quiet=True))
        except SADError:
            # 直前の値の変更で Twiss 配列が無効 (SAD は古い値を保持) かつ現在の光学が計算できない: s だけ持つ行を作る
            rows, s_ = [], 0.0
            for i in range(self.nlat):
                rows.append(TwissRow(s=s_, name='*' if i == 0 else self.els[i - 1].name))
                if i < len(self.els):
                    s_ += self.els[i].L
        ini = [float(x) for x in ini]
        delta = ini[20]
        els = self.els[lb - 1:le - 1]
        sub = _lat.Lattice(els, name='CALCOPT') if els else None
        stab = [True, True]
        qq = None
        new = None
        if idim == 3:
            new, stab = self._calcopt_6d(lb, le, ini, cell, rows)
        elif cell:
            if sub is None:
                raise SADError('General::wrongval: empty region for cell=1 in CalculateOptics')
            cod_start = np.array(ini[15:20] + [delta], dtype=float)
            try:
                # 不安定な周期系でも SAD の qcell1 と同じく続行する (|cos μ|>1 は 1/cos μ に置き換え、INTRES/HALFRES で stabx/staby を判定)
                sub_rows = sub.twiss_table(delta=delta, cod=cod_start,
                                           sad_stab={'intres': bool(self.flags['INTRES']), 'halfres': bool(self.flags['HALFRES'])})
                qq = sub._qstat
            except (ValueError, RuntimeError, FloatingPointError, ZeroDivisionError, np.linalg.LinAlgError, SADError):
                sub_rows = None
            if sub_rows is None or not all(math.isfinite(getattr(r, k)) for r in sub_rows for k in ('betx', 'bety', 'mux', 'muy')):
                self.write(' *****qcod---> Overflow & closed orbit not found            1' + chr(10))
                return [[[float(getattr(rows[k - 1], self._TWKEY[k2], 0.0)) for k2 in _TWISS_ALL] for k in range(lb, le + 1)],
                        [0.0, 0.0, 1.0], [0.0, 0.0, 0.0]]
            new = [TwissRow(**{**sub_rows[i].__dict__, 's': rows[lb - 1 + i].s, 'name': rows[lb - 1 + i].name,
                               'mux': sub_rows[i].mux + ini[2], 'muy': sub_rows[i].muy + ini[5]})
                   for i in range(len(sub_rows))]
        else:
            tw0 = dict(ax=ini[0], bx=ini[1], nx=ini[2], ay=ini[3], by=ini[4], ny=ini[5], ex=ini[6], epx=ini[7], ey=ini[8],
                       epy=ini[9], r1=ini[10], r2=ini[11], r3=ini[12], r4=ini[13], detr=ini[14])
            # N3 (pysad の不具合の修正、モードによらない): TRPT かつ RFSW の加速 (CALC の _ins_rows と同じ trptacc) を区間に入れる
            acc = self._trpt_acc()
            if acc is not None:
                o = lb - 1
                sidx = {i - o for i in acc['idx'] if o <= i < le - 1}
                acc = ({'gb': acc['gb'][o:], 'idx': sidx, 'step': (lambda el, i, cod, _s=acc['step'], _o=o: _s(el, i + _o, cod))}
                       if sidx else None)
            tws = _lat.ins_twiss(els, tw0, ini[15:20], delta, acc=acc)
            new = []
            for i, tw in enumerate(tws):
                c = tw['cod']
                new.append(TwissRow(s=rows[lb - 1 + i].s, name=rows[lb - 1 + i].name,
                                    alfx=tw['ax'], betx=tw['bx'], mux=tw['nx'], alfy=tw['ay'], bety=tw['by'], muy=tw['ny'],
                                    etx=tw['ex'], etpx=tw['epx'], ety=tw.get('ey', 0.0), etpy=tw.get('epy', 0.0),
                                    r1=tw.get('r1', 0.0), r2=tw.get('r2', 0.0), r3=tw.get('r3', 0.0), r4=tw.get('r4', 0.0),
                                    detr=tw.get('detr', 0.0),
                                    x_cod=float(c[0]), px_cod=float(c[1]), y_cod=float(c[2]), py_cod=float(c[3]),
                                    z_cod=float(c[4]), delta_cod=float(c[5])))
        for i, r in enumerate(new):
            rows[lb - 1 + i] = r
        self.rows = rows
        # qcell1 の末尾: 位置 lb と le の Twiss から区間の転送行列の対角和 (x11+x22)
        b, e = rows[lb - 1], rows[le - 1]
        trace = []
        for (ab, bb, pb), (ae, be_, pe) in (((b.alfx, b.betx, b.mux), (e.alfx, e.betx, e.mux)),
                                           ((b.alfy, b.bety, b.muy), (e.alfy, e.bety, e.muy))):
            d = pe - pb
            c_, s_ = math.cos(d), math.sin(d)
            trace.append(math.sqrt(be_ / bb) * (c_ + ab * s_) + math.sqrt(bb / be_) * (c_ - ae * s_))
        st = [bool(stab[0]), bool(stab[1])]
        if cell:
            # qcell1: |trace| > 2 のとき INTRES/HALFRES が真なら不安定、偽なら 4/trace (fam=False)
            for j in (0, 1):
                if trace[j] > 2.0:
                    if self.flags['INTRES']:
                        st[j] = False
                    else:
                        trace[j] = 4.0 / trace[j]
                if trace[j] < -2.0:
                    if self.flags['HALFRES']:
                        st[j] = False
                    else:
                        trace[j] = 4.0 / trace[j]
        if cell and qq is not None and idim != 3:
            st = [bool(qq['stabx']), bool(qq['staby'])]       # qcell1 の判定そのもの (qmdiag の stab、codfnd、over も含む)
            trace = [qq['tracex'], qq['tracey']]
        out = [[self.twiss_value(k2, k) for k2 in _TWISS_ALL] for k in range(lb, le + 1)]
        stz, trz = 1.0, 0.0
        if _compat.v2() and idim == 3:
            # v2 (B16): idim=3 は縦方向の対角和と安定性も返す (K64 は未初期化)。idim=2 は縦を計算しないので stabz=1、tracez=0 と決める
            try:
                d = e.muz - b.muz
                trz = (math.sqrt(e.betz / b.betz) * (math.cos(d) + b.alfz * math.sin(d))
                       + math.sqrt(b.betz / e.betz) * (math.cos(d) - e.alfz * math.sin(d)))
                stz = 1.0 if (math.isfinite(trz) and abs(trz) <= 2.0) else 0.0
            except (ValueError, ZeroDivisionError, AttributeError):
                stz, trz = 0.0, 0.0
        return [out, [1.0 if st[0] else 0.0, 1.0 if st[1] else 0.0, stz], [trace[0], trace[1], trz]]

    def _calcopt_6d(self, lb, le, ini, cell, rows):
        """CalculateOptics の idim=3 (qcell6d): 区間全体 (lb=1, le=nlat) の周期 6 次元解のみ"""
        if not (cell and lb == 1 and le == self.nlat):
            raise SADError('CalculateOptics: idim=3 は cell≠0 で全周 (1, LINE["LENGTH"]) のみ実装')
        save = (self.flags['CALC6D'], self.flags['CELL'])
        self.flags['CALC6D'] = True
        self.flags['CELL'] = True
        try:
            new = self.calc()
        finally:
            self.flags['CALC6D'], self.flags['CELL'] = save
        return list(new), [True, True]

    def fn_dasurvey(self, args):
        """DynamicApertureSurvey[range, nt, opt...] (Packages/dapert.n): 既定 Output->0, ExpandElementValues->True,
        ReferenceOrbit->Default, PhaseX->0, PhaseY->0, PhaseZ->-Pi/2, Axes->"ZX", Damping->False, DampEnough->0"""
        from .dapert import dynamic_aperture
        self.spread_pending()
        if len(args) < 2:
            raise SADError('DynamicApertureSurvey[range, nturn, options]')
        rng = args[0]
        nt = int(args[1])
        opts = {}
        for a in args[2:]:
            for r in (a if isinstance(a, list) else [a]):
                if isinstance(r, Rule) and isinstance(r.lhs, Sym):
                    opts[r.lhs.name] = r.rhs
        lfno = opts.get('Output', 0.0)
        lfno = int(lfno) if is_num(lfno) else 0
        phi = [float(opts.get('PhaseX', 0.0)), float(opts.get('PhaseY', 0.0)), float(opts.get('PhaseZ', -math.pi / 2))]
        axes = opts.get('Axes', 'ZX')
        axes = axes.upper() if isinstance(axes, str) else 'ZX'
        modes = {'ZX': (3, 1), 'ZY': (3, 2), 'XY': (1, 2), 'XZ': (1, 3), 'YX': (2, 1), 'YZ': (2, 3)}
        if axes not in modes:
            raise SADError(f'Wrong Axes {axes}')
        iv1, iv2 = modes[axes]
        iv3 = 6 - iv1 - iv2
        damping = opts.get('Damping', 0.0)
        if truthy(damping):
            raise SADError('DynamicApertureSurvey: Damping->True は未対応')
        if not (isinstance(rng, list) and len(rng) == 3):
            raise SADError('DynamicApertureSurvey: range must be a list of 3')
        r1 = [float(x) for x in rng[iv1 - 1]]
        n1 = min(201, len(r1))
        g = self.I.globals
        dapw = float(g.get('DAPWIDTH', 51.0))
        turns0, np0 = g.get('TURNS', 1.0), g.get('NP', 1000.0)
        g['TURNS'] = float(nt)
        g['NP'] = float(min(4000, n1 * min(dapw, 51) * 2))
        out = self.write if lfno == 6 else (lambda s: None)
        if lfno == 6:
            self.write(' Turns =' + self.I.to_str(float(nt)) + ' Maximum number of particles =' + self.I.to_str(g['NP']) + chr(10))
        if self.emit is None:
            self.run_emit(quiet=True)
        emx = math.sqrt(abs(float(g.get('EMITX', 0.0))) + abs(float(g.get('EMITY', 0.0))))
        emz = math.sqrt(abs(float(g.get('EMITZ', 0.0)))) if self.flags['RFSW'] else abs(float(g.get('SIGE', 0.0)))
        range3 = [[float(rng[iv2 - 1][0]), float(rng[iv2 - 1][1])], [float(rng[iv3 - 1][0]), float(rng[iv3 - 1][1])]]
        score, result = dynamic_aperture(self.els, range3, r1[:n1], nt, phi, iv1, iv2, self.momentum(), self.emit,
                                         emx, emz, self.trf0, self.flags['RFSW'], self.flags['RAD'], fluc=self.flags['FLUC'],
                                         alost=float(g.get('LOSSAMPL', 1.0)), zlost=float(g.get('LOSSDZ', 100.0)),
                                         out=out, progress=int(float(g.get('PYSAD$PROGRESS', 0.0))), fshift=self._fshift())
        g['ResultOfDAPERT'] = result
        g['NP'] = np0
        g['TURNS'] = turns0
        return [score, [rng, result[1]]]

    def fn_touschek(self, args):
        """TouschekLifetime[nx, ny, nz] (Packages/emit.n): TouschekTable (INTRA ON の emit/Emittance[] が作る) をスプライン補間して寿命 [s]。
        nx, ny は INF (= sadscript.INF) か実数、sige は SIGE グローバル。"""
        from .touschek import touschek_lifetime
        tt = getattr(self, 'touschek', None)
        if tt is None or tt.get('table') is None:
            raise SADError('Touschek::undef:  TouschekTable Undefined; Do Emittance[] with INTRA;')
        if len(args) != 3:
            raise SADError('TouschekLifetime[nx, ny, nz]')
        nx, ny, nz = [(float(a) if is_num(a) else math.inf) for a in args]
        sige = float(self.I.get('SIGE', 0.0))
        return float(touschek_lifetime(tt['table'], nx, ny, nz, sige))

    # ------------------------------------------------------------------ MAIN レベルの命令 (FFS; の前)
    def main_exec(self, c: dict):
        """parser.parse_main_commands の 1 命令を実行する (ON/OFF, EMIT, SET, TRACK)。
        MAIN の EMIT は USE=名前 でラインを選び EMIT を実行 (報告も出す)。TRACK は maintrack.run_track。"""
        from .maintrack import run_track
        k = c['cmd']
        if k in ('ON', 'OFF'):
            for x in c['flags']:                 # MAIN 部の `ON RAD FLUC;` / `OFF RAD;` は有効 (FFS 内の ON/OFF は無効: _exec_word)
                X = x.upper()
                if X in _FLAG_ALT:
                    name, val = _FLAG_ALT[X]
                    val = val if k == 'ON' else (not val)
                elif X in self.flags:
                    name, val = X, (k == 'ON')
                else:
                    break
                if name == 'RFSW':               # $RFSW$ は MAIN の TRACK 用。FFS の RFSW は常に ON で始まる (track.f:149)
                    self.main_flags['RFSW'] = val
                    if _compat.v2():
                        self.flags['RFSW'] = val     # v2 (B30): MAIN の ON/OFF RFSW を、ほかのフラグと同じく FFS に届ける
                else:
                    self.flags[name] = val
        elif k == 'EMIT':
            if c.get('use') and c['use'].upper() != self.name.upper():
                self.use(c['use'])
            elif self.lat is None:
                self.write(' ???-EMIT-Missing beamline in USE.' + chr(10))
                return
            self.run_emit()
        elif k == 'SET':
            self.I.globals[c['name']] = float(c['value'])      # PHSPLOTS は描画の出力番号 (値だけ保持して使わない)
        elif k == 'TRACK':
            res = run_track(self, c, self.write)
            self.track_results.append(res)

    # ------------------------------------------------------------------ コマンドの実行
    def execute(self, stmt: str):
        if self.stopped:
            return
        words = self._split_words(stmt)
        i = 0
        while i < len(words):
            w = words[i]
            W = w.upper()
            n = self._exec_word(words, i)
            if n is None:
                break
            i += max(1, n)

    @staticmethod
    def _split_words(stmt: str) -> List[str]:
        """空白区切り。ただし {..} [..] "..." の中は分けない"""
        out = []
        buf = ''
        depth = 0
        q = False
        for c in stmt:
            if q:
                buf += c
                if c == '"':
                    q = False
                continue
            if c == '"':
                q = True
                buf += c
            elif c in '{[(':
                depth += 1
                buf += c
            elif c in '}])':
                depth -= 1
                buf += c
            elif c.isspace() and depth <= 0:
                if buf:
                    out.append(buf)
                    buf = ''
            else:
                buf += c
        if buf:
            out.append(buf)
        return out

    def _eval_word(self, w: str):
        try:
            return self.I.run(w, 'script')
        except Exception:
            raise SADError(f'?Cannot evaluate {w}')

    def _getva(self, words, i):
        """words[i] を数値として評価。(値, 消費数) — 評価できなければ (None, 0)"""
        if i >= len(words):
            return None, 0
        try:
            v = self.I.run(words[i], 'script')
        except Exception:
            return None, 0
        if is_num(v):
            return float(v), 1
        if isinstance(v, list) and all(is_num(x) for x in v):
            return v, 1
        return None, 0

    def _exec_word(self, words, i) -> Optional[int]:
        """words[i] から始まるコマンドを実行し、消費した語数を返す。None なら文の終わり"""
        w = words[i]
        W = w.upper()
        rest = words[i + 1:]
        # ---- ON/OFF: SAD K64 の FFS では無効 (bench_level29d)。語そのものは何もせず `Out[1]:= USE` を出すだけで、後ろの語は
        # 普通に次の文として読まれる。つまり `ON INTRA;` は `INTRA;` (フラグを立てる) と同じで、`OFF INTRA;` も INTRA を **立てる**
        # (切るには `NOINTRA;`)。`ON NORAD NOFLUC;` は RAD と FLUC を切る。MAIN 部の ON/OFF は有効 (main_exec)
        if W in ('ON', 'OFF'):
            self.write('Out[1]:= USE' + chr(10))
            return 1
        # ---- フラグ
        if W in _FLAG_ALT:
            name, val = _FLAG_ALT[W]
            self.flags[name] = val
            return 1
        if W in self.flags:
            self.flags[W] = True
            return 1
        if W.startswith('NO') and W[2:] in self.flags:
            self.flags[W[2:]] = False
            return 1
        # ---- FIT 条件 (tgetfv)
        fun = W
        maxfit = False
        if W not in _FUNS and W[:-1] in _FUNS and W.endswith('M'):
            fun, maxfit = W[:-1], True
        if fun in _FUNS:
            return self._fit_condition(fun, maxfit, words, i)
        if W.endswith('I') and W[:-1] in _FUNS:
            v, n = self._getva(words, i + 1)
            if v is None:
                self.write(f'?Missing number for {w}\n')
                return 1
            self._set_initial(W[:-1], float(v))
            return 1 + n
        if W.endswith('SCALE') and W[:-5] in _FUNS:
            v, n = self._getva(words, i + 1)
            return 1 + n
        # ---- 誤差
        if W in _ERRCMDS:
            return self._error_cmd(W, words, i)
        # ---- コマンド
        if W in ('STOP', 'EXIT', 'END', 'QUIT', 'ABORT'):
            self.stopped = True
            self.I.stopped = True
            return None
        if W == 'USE' or W == 'VISIT':
            nm = None
            expnd, nskip = True, 0
            if rest and _abbrev(rest[0].upper(), 'NOEXP_AND'):       # tffsa.f: USE [NOEXPAND|EXPAND] line
                expnd, nskip = False, 1
            elif rest and _abbrev(rest[0].upper(), 'EXP_AND'):
                nskip = 1
            if nskip:
                rest = rest[1:]
            plain = bool(rest) and re.fullmatch(r'[A-Za-z0-9_$.]+', rest[0]) is not None
            if rest and rest[0] in self.env._beamlines:
                nm = rest[0]                         # gs: LINE 名は書かれたとおり (USE ring と USE RING は別。ring が無ければ Missing beamline)
            elif rest and not plain:
                nm = self._use_beamline_expr(' '.join(rest))
            if nm is None:
                # tffsa.f: ビームラインが未定義なら 'Missing beamline in USE.' (前の USE はそのまま)
                if rest and not plain:
                    self.write(' ???-FFS-Error-Missing beamline in USE. ' + chr(10))
                else:
                    self.write('???General::wrongtype:  Argument must be BeamLine[ ... ]:\n ???-FFS-Error-Missing beamline in USE. \n')
                return 1 + nskip if not rest else len(words) - i
            if W == 'VISIT':
                self._visit_stack.append(self.name)
                # tffsa.f: VISIT は tfblocksym で 5 つの記号の定義を隠す (BYE で戻る)。gs で確かめた: 訪問先では未定義
                self._visit_syms.append(self._block_syms())
            else:
                # USE は tfunblocksym で 5 つの記号を消す (gs で確かめた)。VISIT の中の USE は、その訪問先の分だけ
                self._block_syms()
            self.use(nm, expand=expnd)
            self.flags['CALEXP'] = expnd            # tffsa.f: calexp=expnd
            return len(words) - i
        if W == 'BYE':
            if self._visit_stack:
                def _restore(saved):
                    self._block_syms()
                    self.I.globals.update(saved.get('g', {}))
                    self.I.downvalues.update(saved.get('d', {}))
                if i + 1 < len(words) and words[i + 1].upper() == 'ALL':
                    prev = self._visit_stack[0]
                    saved = self._visit_syms[0] if self._visit_syms else {}
                    self._visit_stack = []
                    self._visit_syms = []
                    self.use(prev)
                    _restore(saved)
                    return 2
                self.use(self._visit_stack.pop())
                if self._visit_syms:
                    _restore(self._visit_syms.pop())
            else:
                self.write(' ???-FFS-Error-BYE without VISIT. ' + chr(10))
            return 1
        if _abbrev(W, 'CAL_CULATE'):
            n = 1
            while i + n < len(words) and (words[i + n].upper().rstrip('-') in _FUNS or _abbrev(words[i + n].upper(), 'NOEXP_AND')
                                          or _abbrev(words[i + n].upper(), 'EXP_AND')):
                n += 1
            # tffsa.f: CALC [funs] の直後の語が NOEXPAND なら calexp=.false.、EXPAND なら .true. (KEEPEXP なので以後の CALC/GO に残る)
            if n > 1 and _abbrev(words[i + n - 1].upper(), 'NOEXP_AND'):
                self.flags['CALEXP'] = False
            elif n > 1 and _abbrev(words[i + n - 1].upper(), 'EXP_AND') and words[i + n - 1].upper() not in _FUNS:
                self.flags['CALEXP'] = True
            elif not self.flags.get('KEEPEXP', True):
                self.flags['CALEXP'] = not self._tfvcomp()      # NOKEEPEXP (tffsa.f:553): calexp=.not. tfvcomp()
            self.go(fit=False)
            return n
        if W == 'GO':
            n = 1
            if i + n < len(words) and _abbrev(words[i + n].upper(), 'NOEXP_AND'):
                self.flags['CALEXP'] = False
                n += 1
            elif i + n < len(words) and _abbrev(words[i + n].upper(), 'EXP_AND'):
                self.flags['CALEXP'] = True
                n += 1
            elif not self.flags.get('KEEPEXP', True):
                self.flags['CALEXP'] = not self._tfvcomp()      # NOKEEPEXP (tffsa.f:569)
            self.go(fit=True)
            return n
        if W == 'EMIT':
            self.flags['GAUSS'] = True            # SAD: EMIT コマンドは gauss=.true. にする (Emittance[] はしない。11/06 に SAD で確認)
            self.run_emit()
            return 1
        if W in ('FIT',):
            return self._fit_location(words, i)
        if _abbrev(W, 'FITP_OINTS'):
            v, n = self._getva(words, i + 1)
            if v is not None:
                self.matcher_().fitp(int(v))
            return 1 + n
        if W in ('FREE', 'FIX'):
            return self._freefix(W == 'FREE', words, i)
        if _abbrev(W, 'MAXI_TERATION'):
            v, n = self._getva(words, i + 1)
            if v is not None:
                self.I.globals['MAXITERATION'] = float(v)
            return 1 + n
        if _abbrev(W, 'T_YPE'):
            self.type_cmd(rest)
            return len(words) - i
        if _abbrev(W, 'DISP_LAY'):
            return self.disp_cmd(words, i)
        if W == 'SAVE' or _abbrev(W, 'EXPAND'):
            n = 1
            pats = []
            while i + n < len(words) and (words[i + n].upper() == 'ALL' or self.families(words[i + n])):
                pats.append(words[i + n])
                n += 1
            self.save_cmd(pats)
            return n
        if W == 'RESET':
            self.reset_cmd()
            return 1
        if _abbrev(W, 'REJ_ECT'):
            m = self.matcher_()
            n = 1
            while i + n < len(words):
                w2 = words[i + n].upper()
                if w2 in ('TOTAL', 'TOTALFIT'):
                    m.conditions = []
                elif w2 in _FUNS or w2[:-1] in _FUNS:
                    kf = w2 if w2 in _FUNS else w2[:-1]
                    kp, kp1 = min(m._mfpnt, m._mfpnt1), max(m._mfpnt, m._mfpnt1)
                    m.conditions = [c for c in m.conditions if not (c['kp'] == kp and c['kp1'] == kp1 and c['kf'] == kf)]
                else:
                    break
                n += 1
            return n
        if _abbrev(W, 'DELC_OR'):
            return self._delcor(words, i)
        if W == 'WAKE':
            return self._wake_cmd(words, i)
        if _abbrev(W, 'COUP_LE'):
            return self._couple(words, i)
        if _abbrev(W, 'DECOUP_LE') or _abbrev(W, 'INDEP_ENDENT'):
            self.couples = {}
            return 1
        if W in ('GCUT',):
            v, n = self._getva(words, i + 1)
            if v is not None:
                self.I.globals['GCUT'] = float(v)
            return 1 + n
        if W in ('DP', 'DP0'):
            v, n = self._getva(words, i + 1)
            if v is not None:
                self.I.globals[W] = float(v)
            return 1 + n
        if W in ('CODPLOT', 'NOCODPLOT'):
            self.flags['CODPLOT'] = W == 'CODPLOT'
            return 1
        if W in ('MCOD', 'TCOD', 'MON'):
            # この K64 ビルドの SAD では MON/MCOD/TCOD は FFS コマンドではなく式として評価されて何もしない (SAD で確認)。
            # 後ろの語は普通に次の文として読まれるので、`MON 0 0 0 0 PM*` は `PM*` の値が無く ?Missing value for PM* になる
            # (<work-userscripts>/esrf correct_cod, print_emittance: SAD は FFS::error で Do を抜ける)
            if _compat.v2():
                # v2 (B22): K64 にこのコマンドは無い。文の残りを式として読んで別のエラー (Missing value for PM*) で止まるのを避け、
                # 無いことを知らせて文全体を読み飛ばす (軌道補正は CorrectOrbit パッケージ: docs/34)
                self.write(f' ???-FFS-Warning-{W} is not a command of SAD K64 (ignored): {" ".join(words[i:])} ' + chr(10))
                return len(words) - i
            return 1
        if W == 'ORG' or _abbrev(W, 'ORI_GIN'):
            return self._org(words, i)
        if _abbrev(W, 'CHRO_MATICITY'):
            self.chro_cmd()
            return 1
        if W == 'RADINT':
            self.radint_cmd()
            return 1
        if _abbrev(W, 'BEAM_SIZE'):
            # tffsa.f BEAM_SIZE: tfsetparam + tfsize(.true.)。pysad の仕様: ビームサイズは EMIT の結果から作るので、EMIT 済みかだけ確かめる
            if self.emit is None:
                self.write('?BEAMSIZE: EMIT を先に実行してください (pysad: ビームサイズは EMIT の EMITX/EMITY/SIGE と Twiss から)' + chr(10))
            return 1
        if _abbrev(W, 'RENUM_BER'):
            return self._renumber_cmd(words, i)
        if _abbrev(W, 'REF_ERENCE'):
            # tfsetref: twiss(:,-1,:)=twiss(:,0,:) (直前の CALC の光学を参照光学にする)
            import copy as _copy
            self._ref_rows = [_copy.copy(r) for r in self.ensure_rows()]
            return 1
        if (_abbrev(W, 'MEA_SURE') or _abbrev(W, 'TRC_OD') or _abbrev(W, 'LOGM_EASURE')
                or (_compat.v2() and (_abbrev(W, 'ORBITJ_ITTER') or _abbrev(W, 'PUTM_EASURE') or _abbrev(W, 'GET_MEASURE')))):
            # v2: MEASURE / TRCOD / LOGMEASURE と ORBITJITTER / PUTMEASURE / GETMEASURE (measure.py)。k64 は従来どおり (gs は MEASURE が無出力、TRCOD は SIGSEGV: B42)
            # (_compat はモジュールの先頭で import 済み。この関数の中で import し直すと、前の方の _compat の参照が UnboundLocalError になる)
            if _compat.v2():
                from . import measure as _ms
                r_ = _ms.ffs_command(self, W, words, i)
                if r_ is not None:
                    return r_
        from . import ffsextra as _fx
        if W in _fx.COMMANDS:
            r_ = _fx.run(self, W, words, i)
            if r_ is not None:
                return r_
        if W in ('DRAW', 'GEO', 'SHOW', 'STAT', 'STATUS', 'VARS', 'VARIABLES', 'VARY',
                 'BEAM', 'MATRIX', 'MAT', 'ALIGN', 'ATT', 'ATTRIBUTE',
                 'INC', 'INCLUDE', 'REC', 'RECOVER', 'REJ', 'REJECT', 'SCALE', 'RENUM', 'REV', 'REVERSE', 'ZAP', 'DIMAD',
                 'SYNCHROB', 'REF', 'REFERENCE', 'MEASURE', 'TRC', 'TRCOD', 'PRINT'):
            self.write(f'[pysad] FFS command {W} は未対応 (無視): {" ".join(words[i:])}\n')
            return len(words) - i
        # ---- 要素値 (tfgetv): "QF K1 0.5", "S* 0", "SDA MAX 0"
        m_ = re.fullmatch(r'(.+)\.(\d+)', w)
        if m_ and not any(c in w for c in '*%{|') and self.families(m_.group(1)):
            # `QF.3 K1 0.9;`: 要素名の定義済みシンボル + 実数なので、gs (tfprint) は文全体を式として評価するだけで何も設定しない
            # (成分ごとの設定は LINE["K1","QF.3"]=0.9)
            return len(words) - i
        if self.families(w) or self._position(w) is not None:
            return self._element_value(words, i)
        if re.fullmatch(r'[-+]?(?:\d+\.?\d*|\.\d+)(?:[EeDd][-+]?\d+)?', w):
            return 1                                    # tfprint: 数値の語は式として評価されるだけ (SAD は Out[n]:= 値 を出す)
        self.write(f'?Undefined command or element: {w}\n')
        return len(words) - i

    # ------------------------------------------------------------------ USE
    def _use_beamline_expr(self, text: str):
        """USE / VISIT の引数が名前でなければ SADScript の式 (BeamLine[...] や Join[a,-a]) として評価し、
        BeamLine なら新しい LINE (名前 L%06d$) を MAIN に作って名前を返す。先頭が MARK でなければ FFS::firstmark エラー。"""
        try:
            v = self.I.run(text, 'script')
        except (SADError, ValueError, TypeError, KeyError):
            return None
        if not (isinstance(v, Expr) and v.head == S('BeamLine')) or not v.args:
            return None
        refs = []
        for x in v.args:
            neg = isinstance(x, Expr) and x.head == S('Minus') and len(x.args) == 1
            y = x.args[0] if neg else x
            n = y.name if isinstance(y, Sym) else y
            if not isinstance(n, str):
                return None
            n = n.upper()
            if n not in self.env._elements and n not in self.env._beamlines:
                self.write(f'???MAIN::wrongtype:  Wrong type of symbol {n} in MAIN:' + chr(10))
                return None
            refs.append(('REVERSE', n) if neg else n)
        first = refs[0]
        first_name = first[1] if isinstance(first, tuple) else first
        first_el = self.env._elements.get(first_name)
        if first_el is None or first_el.get('TYPE', '').upper() != 'MARK':
            self.write('???FFS::firstmark: The 1st component must be MARK:' + chr(10) + first_name[:1] + chr(10) + '^' + chr(10))
            return None
        self._newline_count += 1
        nm = f'L{self._newline_count:06d}$'
        self.env.define_beamline(nm, refs)
        return nm

    def _block_syms(self):
        """tfblocksym / tfunblocksym: FitFunction・FitValue・ElementValues・OpticsProlog・OpticsEpilog の定義 (値と := の規則) を取り出して消す。取り出したものを返す"""
        I = self.I
        return {'g': {k: I.globals.pop(k) for k in _BLOCKED_SYMS if k in I.globals},
                'd': {k: I.downvalues.pop(k) for k in _BLOCKED_SYMS if k in I.downvalues}}

    def use(self, name: str, expand: bool = True):
        name = self.env.resolve_beamline(name) or name.upper()      # LINE 名は書かれたとおり (LinacPFBT)
        # tffsa.f: ilist(2,idval(iuse)) <= 0 .or. expnd のときだけ expnln (展開済みの LINE を NOEXPAND で USE すると成分ごとの値が残る)
        cache = self.__dict__.setdefault('_built_lines', {})
        lat = cache.get(name) if not expand else None
        if lat is None:
            lat = self.env.build_lattice(name)
        cache[name] = lat
        self.lat = lat
        self.name = name
        self.els = lat.elements
        self.nlat = len(self.els) + 1
        self._install_map_hooks()
        for el in self.els:
            if not hasattr(el, '_errk'):
                el._errk = 1.0
        self.rows = None
        self.matcher = None
        self._spread_pending = set()
        self._disp_id1, self._disp_id2 = 1, self.nlat      # tffsa.f:108 id1=1, id2=nlat (USE のたびに)
        self._z_stale = None                              # 縦方向の列 (AZ..ZPY) の直前の 6D の値 (USE で初期化)
        self.couples = {}
        self.emit = None
        self.touschek = None
        self.geo, self.geo_dirty = None, True
        self.org_k, self.org_g, self.org_chi = 1, [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]
        # trf0 (同期位置の時間) は K64 では tmacro の大域変数で、temit (EMIT・Emittance[]) と CALC6D だけが更新する。USE では 0 に戻らない
        # (gs で確かめた: EMIT のあと USE しても、追跡の RF の位相は EMIT の trf0 のまま。次の EMIT まで残る)。最初の値は 0
        self.trf0 = getattr(self, 'trf0', 0.0)
        # 先頭 MARK の EMITX/EMITY/EMITZ/SIGZ/SIGE/DDP → グローバル (tffsa.f:135-139 の
        # {EMITX,EMITY,EMITZ,SIGZ,SIGE}=LINE[{..},1]; DP0=LINE["DDP",1])。これは FFS に入った最初の 1 回だけで、
        # 以後の USE (chguse=.true. で 101 へ戻る経路) では実行されない: 直前の emit の結果 (失敗時は NaN も) と DP0 が残る
        g = self.I.globals
        if not getattr(self, '_globals_loaded', False):
            self._globals_loaded = True
            m0 = self.els[0]
            raw = self.env._elements.get(m0.name, {}) if type(m0).__name__ == 'Marker' else {}
            for k, gk in (('EMITX', 'EMITX'), ('EMITY', 'EMITY'), ('EMITZ', 'EMITZ'), ('SIGZ', 'SIGZ'), ('SIGE', 'SIGE')):
                # MARK に書かれていない (0 の) ものは、MAIN で代入した値を残す (gs で確かめた: MAIN の `EMITY=2E-11;` は MARK に EMITY が
                # 無ければ FFS でも 2E-11、MARK の EMITX=5E-9 は MAIN の EMITX=1E-9 より優先)。MAIN で代入していなければ従来どおり 0
                g[gk] = float(raw.get(k, 0.0) or g.get(gk, 0.0) or 0.0)
            g['DP0'] = float(raw.get('DDP', 0.0) or 0.0)
            if raw.get('DP'):
                g['DP'] = float(raw['DP'])
        # MAIN 側の MOMENTUM は、値が変わったときだけ FFS のグローバルへ写す (gs: FFS 内の MOMENTUM=… は USE で戻らない)
        m = _env_get(self.env, 'MOMENTUM')
        if isinstance(m, (int, float)) and m > 0 and m != getattr(self, '_env_mom_seen', None):
            self._env_mom_seen = m
            g['MOMENTUM'] = float(m)
        self._set_omega0()

    # ------------------------------------------------------------------ MAP (temap.f、ExternalMap)
    def _install_map_hooks(self):
        """MAP 要素があるときだけ ExternalMap["OPTICS"/"EMIT"/"TRACK"] へのコールバックを入れる (なければ外す)"""
        from . import emit6d as _e6, track6n as _t6n
        has = False
        for i, el in enumerate(self.els):
            if type(el).__name__ == 'Map':
                el._pos = i + 1
                has = True
        if has:
            _lat.MAP_HOOKS['OPTICS'] = self._map_optics
            _e6.MAP_EMIT_HOOK = self._map_emit
        else:
            _lat.MAP_HOOKS.pop('OPTICS', None)
            _e6.MAP_EMIT_HOOK = None
        _t6n.MAP_HOOK = None

    def _extmap(self, *args):
        """ExternalMap[args] を評価。定義がなければ None (K64: 評価されずに返る → 何もしない)"""
        if 'ExternalMap' not in self.I.downvalues:
            return None
        r = self.I.apply_value(S('ExternalMap'), list(args))
        if isinstance(r, Expr):
            return None
        return r

    @staticmethod
    def _mat66(m):
        a = np.array(m, dtype=float)
        return a if a.shape == (6, 6) else None

    def _map_optics(self, l, cod):
        r = self._extmap('OPTICS', float(l), [float(v) for v in cod])
        if not (isinstance(r, list) and len(r) == 2):
            return None
        c1 = np.array(r[0], dtype=float)
        M = self._mat66(r[1])
        if c1.shape != (6,) or M is None:
            self.write("ExternalMap(OPTICS) of %d element did not return {cod(6), trans(6,6)}.\n" % l)
            return None
        return c1, M

    def _map_emit(self, l, cod):
        r = self._extmap('EMIT', float(l), [float(v) for v in cod])
        if not isinstance(r, list) or len(r) not in (2, 4):
            if r is not None:
                self.write("ExternalMap(EMIT) of %d element did not return {cod(6), trans(6,6)} or {cod(6), trans(6,6), dtrans(6,6), dbeam(6,6)}.\n" % l)
            return None
        c1 = np.array(r[0], dtype=float)
        M = self._mat66(r[1])
        if c1.shape != (6,) or M is None:
            self.write("ExternalMap(EMIT) of %d element did not return {cod(6), trans(6,6)} or {cod(6), trans(6,6), dtrans(6,6), dbeam(6,6)}.\n" % l)
            return None
        if len(r) == 2:
            return c1, M, None, None
        dM, dB = self._mat66(r[2]), self._mat66(r[3])
        if dM is None or dB is None:
            return None
        return c1, M, dM, dB

    def _map_track(self, tr, npart, l, nt, S_, n, kid, lost_turn):
        """temap: ExternalMap["TRACK", l, nt, x] (x は (7, 粒子数) の実数のリスト: x px*(1+δ) y py*(1+δ) z δ 生存フラグ、粒子の番号順)。
        同じ形のリストが返れば座標を置き換え、フラグ 0 の粒子は落とす (落ちた粒子の復活は未対応)。返らなければ何もしない"""
        if 'ExternalMap' not in self.I.downvalues:
            return n
        X = [[0.0] * npart for _ in range(7)]
        for s in range(n):
            i = int(kid[s])
            pr = 1.0 + S_[s, 5]
            X[0][i], X[1][i], X[2][i], X[3][i] = float(S_[s, 0]), float(S_[s, 1] * pr), float(S_[s, 2]), float(S_[s, 3] * pr)
            X[4][i], X[5][i], X[6][i] = float(S_[s, 4]), float(S_[s, 5]), 1.0
        r = self._extmap('TRACK', float(l), float(nt), X)
        if not (isinstance(r, list) and len(r) == 7 and all(isinstance(v, list) and len(v) == npart for v in r)):
            return n
        from .track6n import drop_particles, MAXTURN_LOST
        drop = np.zeros(n, dtype=bool)
        for s in range(n):
            i = int(kid[s])
            if float(r[6][i]) == 0.0:
                drop[s] = True
                lost_turn[i] = nt
                continue
            self._map_set(tr, S_, s, r, i)
        n = drop_particles(S_, n, kid, drop)
        # 落ちている粒子 (x(j) が生き返る): K64 の doinject。返された座標で生き返らせる (flag≠0 なら)
        for i in range(npart):
            if float(r[6][i]) == 0.0:
                continue
            if i in set(int(kid[q]) for q in range(n)):
                continue
            s = next((q for q in range(n, S_.shape[0]) if int(kid[q]) == i), None)
            if s is None:
                raise SADError('ExternalMap["TRACK"]: 最初から追跡していない粒子 (flag=0) の復活は未対応')
            self._map_set(tr, S_, s, r, i)
            lost_turn[i] = MAXTURN_LOST
            for k in range(S_.shape[1]):
                S_[s, k], S_[n, k] = S_[n, k], S_[s, k]
            kid[s], kid[n] = kid[n], kid[s]
            n += 1
        return n

    @staticmethod
    def _map_set(tr, S_, s, r, i):
        g = float(r[5][i])
        pr = 1.0 + g
        S_[s, 0], S_[s, 2], S_[s, 4], S_[s, 5] = float(r[0][i]), float(r[2][i]), float(r[4][i]), g
        S_[s, 1], S_[s, 3] = float(r[1][i]) / pr, float(r[3][i]) / pr
        S_[s, 6] = tr.ph.dv(g)


    def matcher_(self) -> Matcher:
        if self.lat is None:
            raise SADError('No beam line (USE first)')
        if self.matcher is not None:
            self.matcher.bipol = bool(self.flags['BIPOL'])        # tffsmatch.f:669: UNIPOL (NOBIPOL) は変数の符号を変えない (gs で確かめた)
            self.matcher.absweit = bool(self.flags['ABSW'])       # RELW (NOABSW): tweigh の重み (gs で確かめた)
        if self.matcher is None:
            self.matcher = Matcher(self.lat, cell=self.flags['CELL'], out=lambda s: self.write(s + '\n'))
            self.matcher.bipol = bool(self.flags['BIPOL'])
            self.matcher.absweit = bool(self.flags['ABSW'])
            self.matcher.ffs = self
            self.matcher._calc_optics = self._matcher_optics
            self.matcher.geo_fn = self._matcher_geo
            self.matcher.fam_setup = self._fam_setup
            self.matcher.calc_optics_fixed = self._matcher_optics_fixed
            self.matcher.ini_fn = self.initial_orbits
            self.matcher.ini_optics = self.ini_optics
            # tffscalc.f:64/398: シンボル OpticsProlog/OpticsEpilog を評価する (:= の定義があればそれが走る。未定義なら何もしない)
            self.matcher.optics_prolog = lambda: self.I.run('OpticsProlog', 'script')
            self.matcher.optics_epilog = lambda: self.I.run('OpticsEpilog', 'script')
        return self.matcher

    def _fam_setup(self):
        """tffamsetup の入力: MatchingAmplitude ({{dp, xamp, yamp}, …})、NFAMP、em = |EMITX|+|EMITY|、先頭 MARK の Twiss (twmov(1,·))"""
        amp = self.I.get('MatchingAmplitude')
        if not isinstance(amp, list) or not amp:
            return None
        nfp = self.I.get('NFAMP')
        nfp = float(nfp) if is_num(nfp) else 0.0
        g = self.I.globals
        em = abs(float(g.get('EMITX', 0.0) or 0.0)) + abs(float(g.get('EMITY', 0.0) or 0.0))
        m0 = self.els[0]
        raw = self.env._elements.get(m0.name, {}) if type(m0).__name__ == 'Marker' else {}
        gv = lambda k, d=0.0: float(raw.get(k, d) or 0.0)
        tw = dict(ax=gv('AX'), bx=gv('BX', 1.0), ay=gv('AY'), by=gv('BY', 1.0),
                  r1=gv('R1'), r2=gv('R2'), r3=gv('R3'), r4=gv('R4'), detr=gv('DETR'))
        if tw['bx'] <= 0.0:
            tw['bx'] = 1.0
        if tw['by'] <= 0.0:
            tw['by'] = 1.0
        return dict(amp=amp, nfp=nfp, em=em, tw=tw)

    def _matcher_optics_fixed(self, delta, cod0):
        """fam 点の光学 (閉軌道を探さない): CELL は Matcher の既定 (与えた軌道のまわりの周期解)、INS はその軌道からの伝搬"""
        if self.lat is None:
            return None
        _lat.RFSW = False
        m = self.matcher
        if self.flags['CELL']:
            return Matcher._calc_optics_fixed_default(m, delta, cod0)
        m._last_qstat = None
        try:
            return self._ins_rows(delta, cod0=cod0)
        except Exception:
            return None

    def _matcher_geo(self):
        """Matcher の幾何 (fit 条件 GX/GY/GZ/CHI1-3 用): 変数を反映した現在の要素から ORG 込みで作り直す (tffscalc の tfgeo)"""
        from . import geo as _geo
        gx, gy, gz = self.org_g
        self.geo = _geo.compute_geo(self.els, gx0=gx, gy0=gy, gz0=gz, chi0=tuple(self.org_chi), iorg=self.org_k,
                                    warn=lambda m: None, momentum_eV=self.momentum(warn=False))
        self.geo_dirty = False
        return self.geo

    def fit_fun_value(self, kf: str, kp: int) -> float:
        """条件表示用の関数値: 幾何 (GX..CHI3) は geo から、それ以外は twiss_value"""
        if kf in ('GX', 'GY', 'GZ', 'CHI1', 'CHI2', 'CHI3'):
            g = self.ensure_geo()
            if kf[0] == 'G':
                return float(g.gxyz(kp)['XYZ'.index(kf[1])])
            return float(g.chi(kp)[int(kf[3]) - 1])
        return self.twiss_value(kf, kp)

    def _matcher_optics(self, delta=None, cod0=None):
        """Matcher._calc_optics の差し替え: CELL/INS と δ (オフモメンタム点) に対応"""
        if self.lat is None:
            return None
        self._check_sol_optics(calc_only=True)       # SOL 領域の要素に写像を付ける (付けないと _pack_el が拒否して qcod の誤報が出る)
        d = self.dp0() if delta is None else delta
        _lat.RFSW = False
        m = self.matcher
        if self.flags['CELL']:
            # CELL: Matcher._calc_optics (不安定でも SAD の qcell1 と同じく続行し、stabx/staby/over/codfnd を返す)
            rows = Matcher._calc_optics(m, d, cod0)
        else:
            m._last_qstat = None
            try:
                rows = self._ins_rows(d)
            except Exception:
                return None
            if any(not np.isfinite(r.betx) or not np.isfinite(r.bety) for r in (rows[0], rows[-1])):
                return None
        if rows is None:
            return None
        if delta is None or delta == self.dp0():
            self.rows = rows
        return rows

    def _mark_init_raw(self):
        """SAD の先頭 MARK に残っている Twiss (twmov の「以前の値」)。不安定な CELL CALC で更新されない値の元になる"""
        m0 = self.els[0] if self.els else None
        if m0 is None or type(m0).__name__ != 'Marker':
            return None
        raw = self.env._elements.get(m0.name, {})
        g = lambda k, d=0.0: float(raw.get(k, d) or 0.0)
        return dict(ax=g('AX'), bx=g('BX', 1.0), ay=g('AY'), by=g('BY', 1.0), ex=g('EX'), epx=g('EPX'), ey=g('EY'),
                    epy=g('EPY'), r1=g('R1'), r2=g('R2'), r3=g('R3'), r4=g('R4'), detr=g('DETR'))

    # ------------------------------------------------------------------ FIT / FREE / FIX / GO
    def _fit_location(self, words, i):
        """FIT [loc [loc1]]: 位置を設定 (tffsa)。位置がなければ (FIT; や FIT NX ...) mfpnt1=mfpnt で、
        ielmf が見つからない語に返す nlat ($$$) になる。条件は消さない (消すのは REJECT TOTAL)"""
        m = self.matcher_()
        if i + 1 >= len(words):
            m._mfpnt = m._mfpnt1 = self.nlat
            return 1
        a = words[i + 1]
        if a.upper() == 'ALL':
            m._mfpnt, m._mfpnt1 = 1, self.nlat
            return 2
        au = a.upper()
        if au in _FUNS or au[:-1] in _FUNS:
            # FIT DX 0 : 位置なし → ielmf が nlat を返し $$$ 扱い。語は次のコマンド (fit 関数) として残す
            m._mfpnt = m._mfpnt1 = self.nlat
            return 1
        pa = self.positions(a)
        if not pa:
            # tffsa.f: ielme が見つからない語 (現在のラインに無い要素名など) は mfpnt=mfpnt1=nlat ($$$) にして
            # go to 12 → tfprint で式として評価 (記号なら何も起きない) され、語は消費される。
            # 続く "DX 0" は $$$ の条件になる (pf2 bump3_r の MKSTRSS、calc_bump_r の B10)
            m._mfpnt = m._mfpnt1 = self.nlat
            return 2
        m._mfpnt = pa[0]
        m._mfpnt1 = pa[0]
        n = 2
        if i + 2 < len(words) and words[i + 2].upper() not in _FUNS and words[i + 2].upper()[:-1] not in _FUNS:
            pb = self.positions(words[i + 2])
            if pb:
                m._mfpnt1 = pb[0]
                n = 3
        return n

    def _fit_condition(self, fun, maxfit, words, i):
        m = self.matcher_()
        if i + 1 >= len(words):
            self.write(f'?Missing value for {fun}\n')
            return 1
        vw = words[i + 1]
        n = 2
        if vw == '*':
            value = None
        elif vw.startswith('@'):
            raise SADError('FIT @ (MARK value) is not supported')
        else:
            val = self._eval_word(vw)
            if isinstance(val, list) and len(val) == 2 and all(is_num(x) for x in val):
                value = (float(val[0]), float(val[1]))
            elif is_num(val):
                value = float(val)
            elif isinstance(val, list) and len(val) == 1:
                value = float(val[0])
            else:
                raise SADError(f'?Bad value for {fun}: {vw}')
        npts = None
        v2, k = self._getva(words, i + 2)
        if v2 is not None and is_num(v2):
            npts = int(v2)
            n += k
        if value is None:
            kp, kp1 = min(m._mfpnt, m._mfpnt1), max(m._mfpnt, m._mfpnt1)
            for c in m.conditions:
                if c['kp'] == kp and c['kp1'] == kp1 and c['kf'] == fun:
                    value = (c['fmin'] / (PI2 if fun in ('NX', 'NY', 'NZ') else 1.0), c['fmax'] / (PI2 if fun in ('NX', 'NY', 'NZ') else 1.0))
                    break
            else:
                self.write(f'No default value of {fun}\n')
                return n
        m.fit(fun + ('M' if maxfit else ''), value, npts)
        return n

    def _set_initial(self, fun, v):
        """BXI 1.0 など: 先頭要素 (MARK) の初期 Twiss"""
        if self.els:
            self.set_value(self.els[0], fun, v)

    def _freefix(self, free: bool, words, i):
        m = self.matcher_()
        n = 1
        items = []
        while i + n < len(words):
            w = words[i + n]
            W = w.upper()
            if W in _ELEM_KEYWORDS and items and items[-1][1] is None and self._kw_of_pattern(items[-1][0], W):
                # tffsfreefix: 直前のパターンに合った要素の型のキーワード (tfkwrd) に一致する語だけがキーワード。
                # 合わなければ (FREE K1 K0 K2 K0 の K2 など) 次のパターン
                items[-1] = (items[-1][0], W)
                n += 1
                continue
            if W in _COMMANDS or W in _FUNS or W in _ERRCMDS:
                break
            for p in w.split(','):
                if p.strip():
                    items.append((p.strip(), None))
            n += 1
        for pat, key in items:
            if free:
                if '.' in pat:
                    self._free_components(pat, key)          # tffsfreefix: 語に . があれば成分 (QF.3 など) の変数
                    continue
                for fam in self.families(pat):
                    if fam in self.couples:
                        continue
                    k_ = (key or self.default_key(self.instances(fam)[0]) or '').upper()
                    comp_ = [v for v in m.variables if v.get('comp') and v['family'] == fam and v['key'] == k_]
                    if comp_:
                        self.write(' ???-FFS-Error-A component has been already used as variable: ' + comp_[0]['label'] + ' ' + chr(10))
                        return n
                    self._free_family(fam, key)
            else:
                self._fix_vars(pat, key)
        return n

    def _sort_vars(self, m):
        """変数の並び (tffsfreefix: 要素 (ivarele) の順。同じ要素の変数は、後から FREE したものが先 = LOOP_I_2 が既存の変数の前に挿入する)"""
        pos = {}
        for i, e in enumerate(self.els):
            pos.setdefault(e.name, i)
        for v in m.variables:
            if 'seq' not in v:
                m._vseq = getattr(m, '_vseq', 0) + 1
                v['seq'] = m._vseq
        m.variables.sort(key=lambda v: (pos.get(v['family'], 0), -v['seq']))

    def _free_components(self, pat: str, key):
        """成分ごとの FREE (tffsfreefix の comp 分岐。変数は ivcomp = 成分の位置)。K64 の規則:
          - 既定のキー (ival: QUAD の K1、BEND の ANGLE、DRIFT の L など) の変数は、同じ名前の要素すべて (成分は先頭の成分に結合: icomp) の値。
            先頭の成分 (QF.1) を指定すると要素全体の変数と同じ (結合している成分を Info で知らせる)。先頭でない成分 (QF.3) の変数は、
            結合しているので値を持たず、フィットに効かない (gs で確かめた)。pysad はこの変数を作らず、Info だけ出す
          - 既定でないキー (BEND の K1 など) の変数は、その成分だけの値 (tfsetv の ievar2: ivcomp == その成分)"""
        m = self.matcher_()
        poss = self.positions(pat)
        found = False
        for k in poss:
            if k >= self.nlat:
                continue
            el = self.els[k - 1]
            fam = el.name
            tn = type(el).__name__
            if tn in ('Marker', 'Monitor'):
                continue
            dk = (self.default_key(el) or '').upper()
            kk = (key or dk).upper()
            if not kk:
                continue
            found = True
            insts = self.instances(fam)
            if kk == dk:
                if el is insts[0]:
                    if any(v['family'] == fam and v['key'] == kk and not v.get('comp') for v in m.variables):
                        continue
                    self._free_family(fam, kk)
                    for e2 in insts[1:]:
                        self.write(' ???-FFS-Info-Component %s is coupled to %s .\n' % (self.elname_of(e2), self.elname_of(el)))
                else:
                    self.write(' ???-FFS-Info-Component %s is coupled to %s .\n' % (self.elname_of(el), self.elname_of(insts[0])))
                continue
            if any(v['family'] == fam and v['key'] == kk and not v.get('comp') for v in m.variables):
                self.write(' ???-FFS-Error-Element already used as variable: ' + fam + ' ' + chr(10))
                return
            if any(v.get('comp') == k and v['key'] == kk for v in m.variables):
                continue
            v0 = self.get_value(el, kk)
            m.variables.append(dict(family=fam, key=kk, els=[el], type=tn, val=v0, val2=v0, comp=k, label=self.elname(k),
                                    vlim=list(getattr(el, '_vlim', {}).get(kk, [0.1, 1e10] if tn == 'Drift' and kk == 'L' else ([-1e10, 1e10] if kk == dk else [-1e99, 1e99])))))
        if not found:
            self.write(' ???-FFS-Error-Undefined element ' + pat + ' ' + chr(10))
        self._sort_vars(m)

    def _fix_vars(self, pat: str, key):
        """FIX: 成分の変数は、語に . があるか (QF.3)、ワイルドカードで成分名が合うときだけ消える (tffsfreefix の else 分岐: 要素全体の語は成分の変数に触れない)"""
        m = self.matcher_()
        comps = [v for v in m.variables if v.get('comp')]
        if not comps:
            m.fix(pat, key)
            return
        sel = set()
        for p in [x.strip() for x in pat.split(',') if x.strip()]:
            if '.' in p or any(c in p for c in '*%{|'):
                sel.update(self.positions(p))
        keep = [v for v in comps if not (v['comp'] in sel and (key is None or v['key'] == key.upper()))]
        m.variables = [v for v in m.variables if not v.get('comp')]
        pats = ','.join(x.strip() for x in pat.split(',') if x.strip() and '.' not in x)
        if pats:
            m.fix(pats, key)
        m.variables = m.variables + keep
        self._sort_vars(m)

    def set_var(self, v, val):
        """変数 v の値を要素に入れる。成分の変数はその成分だけ、ほかは同じ名前の要素すべて (と COUPLE の相手)"""
        if v.get('comp'):
            self.set_value(v['els'][0], v['key'], val)
        else:
            self.set_family(v['family'], v['key'], val)

    def _tfvcomp(self) -> bool:
        """tffsa.f tfvcomp: 成分ごとのフィット変数 (ivcomp≠0) で、既定以外のキーワードのものがあれば真。
        成分ごとの変数 (`comp` つき。既定のキー以外の成分の変数だけが作られる) があれば真"""
        return any(v.get('comp') for v in (self.matcher.variables if self.matcher else ()))

    def _kw_of_pattern(self, pat: str, kw: str) -> bool:
        """pat に合うファミリーのどれかの型 (SAD の型) のキーワード一覧に kw があるか (tffsfreefix の tfkwrd 照合)"""
        from .sadkw import KEYWORDS
        for fam in self.families(pat):
            el = self.instances(fam)[0]
            if kw in KEYWORDS.get(_TYPENAME.get(_cls(el), ''), ()):
                return True
        return False

    def _free_family(self, fam, key):
        m = self.matcher_()
        els = self.instances(fam)
        el0 = els[0]
        tn = type(el0).__name__
        k = (key or self.default_key(el0) or '').upper()
        if not k or tn in ('Marker', 'Monitor'):
            return
        if any(v['family'] == fam and v['key'] == k for v in m.variables):
            return
        m.variables.append(dict(family=fam, key=k, els=els, type=tn, val=self.get_value(el0, k), val2=self.get_value(el0, k),
                                vlim=list(getattr(el0, '_vlim', {}).get(k, [0.1, 1e10] if tn == 'Drift' and k == 'L' else
                                                ([-1e10, 1e10] if k == (self.default_key(el0) or '').upper() else [-1e99, 1e99])))))      # gs: 既定のキーは ±1e10、ほかは ±1e99
        self._sort_vars(m)

    def go(self, fit: bool):
        m = self.matcher_()
        m.cell = self.flags['CELL']
        m.dp0 = self.dp0()
        m.dpmax = self.dp()
        m.itmax = int(float(self.I.get('MAXITERATION', 40)))
        m.rsconv = float(self.I.get('CONVERGENCE', 1e-9))
        m.wexponent = max(0.1, min(4.0, float(self.I.get('ExponentOfResidual', 2.0))))
        m.offmw = float(self.I.get('OffMomentumWeight', 1.0))
        m.nlat = self.nlat
        m.emx = float(self.I.get('EMITX', 0.0))
        m.emy = float(self.I.get('EMITY', 0.0))
        m.coumin = float(self.I.get('MINCOUP', 1e-2))
        # FFS$NumericalDerivative (SAD のグローバル): True なら数値微分、既定は解析微分 (match6d)
        nd = self.I.get('FFS$NumericalDerivative', 0.0)
        m.numerical_derivative = bool(is_num(nd) and float(nd) != 0.0)
        m.couples_active = bool(self.couples)
        m.xix = float(self.I.get('XIX', 0.0)) * PI2
        m.xiy = float(self.I.get('XIY', 0.0)) * PI2
        m.dpm2 = float(self.I.get('DPM', 0.0))
        m.resflags = (self.flags['INTRES'], self.flags['HALFRES'], self.flags['SUMRES'], self.flags['DIFFRES'])
        m.fit_value = self._fit_value_hook if 'FitValue' in self.I.downvalues else None
        m.fit_weight = self._fit_weight_hook if 'FitWeight' in self.I.downvalues else None
        m.var_range = self._var_range_hook if 'VariableRange' in self.I.downvalues else None
        m.var_weight = self._var_weight_hook if 'VariableWeight' in self.I.downvalues else None
        m.dlim = False
        m.calexp = bool(self.flags.get('CALEXP', True))
        if not m.calexp:
            self.write(' ???-FFS-Info-Element values are not expanded. ' + chr(10))
        self.lat._mark_init = self._mark_init_raw()
        m.fit_function = self._fit_function_hook if ('FitFunction' in self.I.downvalues or self.I.has('FitFunction')) else None
        # 変数の set/get を FFS の errk/couple 付きに
        for v in m.variables:
            v['val'] = self.get_value(v['els'][0], v['key'])
            if fit:
                v['val2'] = v['val']             # VARIABLES の Previous / RECOVER の交換相手: GO の前の値
        self._patch_matcher_setters(m)
        try:
            r = m.go(calc_only=not fit)
        except RuntimeError as e:
            # 光学が計算できない (不安定・オーバーフロー): SAD は報告を出して文を終えるので、スクリプトを止めずにエラーの文にする
            self.write(' *****qcod---> Overflow & closed orbit not found            1' + chr(10))
            raise SADError(str(e))
        g = self.I.globals
        g['MatchingResidual'] = m.matching_residual
        g['NetResidual'] = float(m.residual.r)           # tffscalc: NetResidual=r%r (nstab を含まない)、MatchingResidual=res2r(r)
        # SAD K64 の tffscalc は rlist(isl) (isl=ktfsymbolz('UnstableLevel',14)-4: 名前は 13 文字なのに 14 を渡している) に nstab を書くので、
        # UnstableLevel 自体は更新されず常に 0 のまま (不安定なまま終わっても MatchingResidual は 10^n·nstab+r で、UnstableLevel は 0。bench_level28)
        g['UnstableLevel'] = 0.0
        if _compat.v2():
            g['UnstableLevel'] = float(m.residual.nstab)     # v2 (B15): 不安定な点の数 (MatchingResidual = 10^n·nstab + r の nstab)
        self.rows = None
        self.calc(quiet=True, ini=True)
        if not getattr(m, 'no_variable', False):
            self.print_fit_table()               # tffscalc の '?No variable.' (error=.true.) のときは Matched. の行だけ (条件表なし)
        self.print_unstable()
        return r

    def print_unstable(self):
        """tffsa.f: CELL で、最後の tffscalc の optstat のどれかが不安定なら 'Horizontal' / 'Vertical' / 'Horizontal/Vertical' + ' Unstable.'"""
        m = self.matcher
        if m is None or not self.flags['CELL']:
            return
        st = m.last_stab
        s = ''
        if any(not sx for sx, _ in st.values()):
            s = 'Horizontal'
        if any(not sy for _, sy in st.values()):
            s = 'Horizontal/Vertical' if s else 'Vertical'
        if s:
            self.write(' ' + s + ' Unstable.' + chr(10))

    def _patch_matcher_setters(self, m):
        ffs = self

        def apply_vars():
            for v in m.variables:
                ffs.set_var(v, v['val'])
            ffs.rows = None
        m._apply_vars = apply_vars

    def _fit_value_hook(self, name, kf, dp, goal, vnow, *extra):
        fam = name.split('.')[0]
        gval = goal if not isinstance(goal, tuple) else list(goal)
        # SAD の呼び出し形: FitValue[comp, fn, {id, dp}, goal, now] / 2 点条件 FitValue[comp1, comp2, fn, {id, dp}, goal, now1, now2]
        # (id は MatchingAmplitude/InitialOrbits の軌道番号、通常の条件では 0)
        fid = float(getattr(self.matcher, '_cur_fam_id', 0) or 0)        # kfam (fam 点の番号、通常の点は 0)
        if extra:
            r = self.I.eval(Expr(S('FitValue'), [_Lit(fam), _Lit(extra[0].split('.')[0]), _Lit(kf), _Lit([fid, dp]),
                                                 _Lit(gval), _Lit(extra[1]), _Lit(vnow)]))
        else:
            r = self.I.eval(Expr(S('FitValue'), [_Lit(fam), _Lit(kf), _Lit([fid, dp]), _Lit(gval), _Lit(vnow)]))
        if r is Null:
            return None
        if is_num(r):
            return float(r)
        if isinstance(r, list) and len(r) == 2:
            return (float(r[0]), float(r[1]))
        if isinstance(r, Expr):
            return goal
        return None

    def _var_range_hook(self, fam, key, x):
        """VariableRange[element, keyword, v_] := 式 (tffsvarfun id=1)。実数 → そのまま (0 は範囲外)、{min,max} → タプル、他は None"""
        try:
            r = self.I.eval(Expr(S('VariableRange'), [_Lit(fam), _Lit(key), _Lit(x)]))
        except SADError:
            self.write(f'Error in VariableRange {fam} {key}' + chr(10))
            return None
        if is_num(r) and not isinstance(r, complex):
            return float(r)
        if isinstance(r, list) and len(r) == 2 and all(is_num(t) for t in r):
            return (float(r[0]), float(r[1]))
        return None

    def _var_weight_hook(self, fam, key, gw):
        """VariableWeight[element, keyword, v_] := 式: 実数が返ればそれを重みに使う"""
        try:
            r = self.I.eval(Expr(S('VariableWeight'), [_Lit(fam), _Lit(key), _Lit(gw)]))
        except SADError:
            self.write(f'Error in VariableWeight {fam} {key}' + chr(10))
            return None
        return float(r) if is_num(r) and not isinstance(r, complex) else None

    def _luminosity_hook(self, lum):
        self.I.globals['LuminosityMonitor'] = float(lum)
        if 'LuminosityMonitor' in self.I.downvalues:
            try:
                self.I.eval(Expr(S('LuminosityMonitor'), [_Lit(float(lum))]))
            except SADError:
                pass

    def _fit_function_hook(self):
        """FitFunction (:= で遅延定義されることが多い) を評価して実数/実数リストを返す。他は無視"""
        try:
            r = self.I.eval(S('FitFunction'))
        except Exception:
            self.write('Error in FitFunction ' + chr(10))
            return None
        if is_num(r):
            return float(r)
        if isinstance(r, list):
            flat = [x for x in r if is_num(x)]
            return flat if len(flat) == len(r) else None
        return None

    def _fit_weight_hook(self, name, kf, dp, w):
        fam = name.split('.')[0]
        r = self.I.eval(Expr(S('FitWeight'), [_Lit(fam), _Lit(kf), _Lit([0.0, dp]), _Lit(w)]))
        return float(r) if is_num(r) else w

    def _show_orbit_header(self, m, fitted):
        """tfshow.f:70-170 (InitialOrbits): 軌道ごとの列の幅 lf、表示する軌道 jshow (残差の大きいもの nn 本、dp の昇順)、
        DP / Orbit ID / Res. の 3 行を書く。戻り: dict(form, jshow, lf, lfs)"""
        ini = m.ini
        nfam, nfam1 = ini['nfam'], ini['nfam1']
        lw = max(79, min(256, int(self.I.get('PageWidth', 131.0))))
        lcf = 20
        lf = max(6, min(10, (lw - lcf) // (2 * nfam + 2)))
        fmm = any(c['fmin'] != c['fmax'] for c in fitted.values())
        if lf > 6:
            lfs = max(6, min(8, lf, lw - lf * (2 * nfam + 1) - lcf - 2))
            if fmm:
                lfs += 4
            nn = (lw - lcf - 2 - lfs) // lf
        else:
            lfs = max(6, min(8, lf, lw - lf * (2 * nfam + 1) - lcf))
            if fmm:
                lfs += 4
            nn = (lw - lcf - lfs) // lf
        form = f'{lf}.{min(6, lf - 3)}'
        nn = min(nn, nfam - nfam1 + 1)
        res = getattr(m, '_res1', {})
        dps = m._dps_show
        jshow = []
        for i in range(nfam1, nfam + 1):
            for k, jk in enumerate(jshow):
                if res.get(jk, 0.0) <= res.get(i, 0.0):
                    jshow.insert(k, i)
                    break
            else:
                jshow.append(i)
            jshow = jshow[:nn]
        if 0 not in jshow:
            jshow[nn - 1] = 0
        for i in range(1, nn):
            ji = jshow[i]
            for j in range(i):
                if dps[jshow[j]] > dps[ji]:
                    del jshow[i]
                    jshow.insert(j, ji)
                    break
        pre = ' ' * 11 if lfs > 6 else ' ' * 10
        mid = ' ' * 7 if lfs > 6 else ' ' * 6

        def line(label, cells):
            vout = label.ljust(lfs + 3)
            self.write((pre + vout + mid + ''.join(cells)).rstrip() + chr(10))
        line(' DP ', [fmt_real(dps[j] + m.dp0, form) for j in jshow])
        line(' Orbit ID ', [f'{ini["iuid"][j]:6d}'.ljust(lf) for j in jshow])
        line(' Res.', [fmt_real(res.get(j, 0.0), form) for j in jshow])
        return dict(form=form, jshow=jshow, lf=lf, lfs=lfs)

    def print_fit_table(self):
        """CALC/GO 後の条件表示 (SAD の twelm 風): 位置ごとに、$$$ は AX BX NX AY BY NY LENG を全部、
        他は fit した関数だけ。書式: 名前 a10, 'f '/'  ', 関数 a5, 目標 a10 (S 書式), 点数 1 桁 or '#', 1x, 現在値 (10.6)"""
        m = self.matcher
        if m is None or not self.rows:
            return
        fitted = {}
        for c in m.conditions:
            if c['mfitp'] != 0:
                fitted[(c['kp'], c['kf'])] = c
        locs = sorted({kp for kp, _ in fitted} | {self.nlat})
        multi = getattr(m, 'ini', None) is not None
        cols = None
        if multi:
            cols = self._show_orbit_header(m, fitted)
        cells = []
        for kp in locs:
            funs = ['AX', 'BX', 'NX', 'AY', 'BY', 'NY', 'LENG'] if kp == self.nlat else []
            funs += [kf for (p, kf) in fitted if p == kp and kf not in funs]
            funs.sort(key=lambda f: _FUNS.index(f))
            for kf in funs:
                v = self.fit_fun_value(kf, kp)
                sc = PI2 if kf in ('NX', 'NY', 'NZ') else (math.pi / 180.0 if kf.startswith('CHI') else 1.0)
                c = fitted.get((kp, kf))
                vs = fmt_real(v / sc, "10.6")
                if multi and kf in _FUNS and _FUNS.index(kf) <= _FUNS.index('TRY') and kf not in ('LENG',):
                    # tfshow: k <= mfittry (軌道ごとの値) は表示する軌道の数だけ並べる (幅 lf)
                    vs = ''.join(fmt_real(m._tgfun(m._rows_idp[j], kf, kp) / sc, cols['form']) for j in cols['jshow'])
                # 名前の 11 桁目の 'f' は tfshow の kpa==mfpnta .and. kpb==mfpnta1 (現在の FIT 位置の行)
                fm = 'f' if (kp == m._mfpnt and kp == m._mfpnt1) else ' '
                if c is not None:
                    goal = fmt_real(c['fmin'] / sc, 'S10.6')
                    cells.append(f'{self.elname(kp):10s}{fm} {kf:5s}{goal:10s}{abs(c["mfitp"]) - 1:1d} {vs}')
                else:
                    cells.append(f'{self.elname(kp):10s}{fm} {kf:5s} #######  # {vs}')
        if m.cell and not m.cellstab:
            # tfshow: CELL が不安定 (cellstab=偽) なら 1/$$$ の BMAGX, BMAGY (tgetbmagu2: 先頭と終点の Twiss の不一致) を追加で出す
            rw = self.rows
            nm = f'{self.elname(1)}/{self.elname(self.nlat)}'
            for kf, ka, kb in (('BMAGX', 'alfx', 'betx'), ('BMAGY', 'alfy', 'bety')):
                if any(c['kp'] == 1 and c['kp1'] == self.nlat and c['kf'] == kf and c['mfitp'] != 0 for c in m.conditions):
                    continue
                a0, b0 = getattr(rw[0], ka), getattr(rw[0], kb)
                a1, b1 = getattr(rw[self.nlat - 1], ka), getattr(rw[self.nlat - 1], kb)
                with np.errstate(all='ignore'):
                    v = float(0.5 * (np.float64(b0) / b1 + np.float64(b1) / b0 + (a0 * b1 - a1 * b0)**2 / b0 / b1))
                cells.append(f'{nm:11s} {kf:5s} #######  # {fmt_real(v, "10.6")}')
        # tfshow.f:70,418-433: 行の長さ lw = max(79, min(256, 出力の行長 (PageWidth)))。セルを ' ' でつないで、次を足すと lw を超えるなら改行する
        lw = max(79, min(256, int(self.I.get('PageWidth', 131.0))))
        line = ''
        for c in cells:
            c = c.rstrip()
            if line and len(line) + len(c) + 1 > lw:
                self.write(line + chr(10))
                line = ''
            line = c if not line else line + ' ' + c
        if line:
            self.write(line + chr(10))

    # ------------------------------------------------------------------ 要素値コマンド (tfgetv)
    def _element_value(self, words, i):
        pat = words[i]
        fams = self.families(pat)
        if not fams:
            self.write(f'?Undefined element {pat}\n')
            return 1
        n = 1
        key = None
        mode = None
        if i + n < len(words):
            W = words[i + n].upper()
            if _abbrev(W, 'R_ELATIVE'):
                mode = 'rel'
                n += 1
            elif W == 'MAX':
                mode = 'max'
                n += 1
            elif W == 'MIN':
                mode = 'min'
                n += 1
            elif W in ('MAXABS', 'ABSMAX', 'MAXMIN', 'MINMAX'):
                mode = 'maxmin'
                n += 1
            elif W in _ELEM_KEYWORDS or W in _KEYS.get(type(self.instances(fams[0])[0]).__name__, []):
                key = W
                n += 1
        v, k = self._getva(words, i + n)
        if v is None:
            self.write(f'?Missing value for {pat}\n')
            return n
        n += k
        for fam in fams:
            el0 = self.instances(fam)[0]
            kk = key or self.default_key(el0)
            if kk is None:
                self.write(f'?No default keyword for {fam}\n')
                continue
            if mode in ('max', 'min', 'maxmin'):
                vl = getattr(el0, '_vlim', None)
                if vl is None:
                    vl = {}
                    for el in self.instances(fam):
                        el._vlim = vl
                lim = vl.setdefault(kk, [0.1 if type(el0).__name__ == 'Drift' and kk == 'L' else -1e10, 1e10])
                cur = self.get_value(el0, kk)
                if mode == 'max':
                    lim[1] = v
                    cur = min(v, cur)
                elif mode == 'min':
                    lim[0] = v
                    cur = max(v, cur)
                else:
                    lim[0], lim[1] = -abs(v), abs(v)
                    cur = min(abs(v), max(-abs(v), cur))
                self.set_family(fam, kk, cur)
                if self.matcher:
                    for var in self.matcher.variables:
                        if var['family'] == fam and var['key'] == kk:
                            var['vlim'] = list(lim)
            elif mode == 'rel':
                self.set_family(fam, kk, self.get_value(el0, kk) * (1.0 + v))
            else:
                self.set_family(fam, kk, v)
        self.invalidate()
        return n

    # ------------------------------------------------------------------ 誤差 (terror)
    def _error_cmd(self, W, words, i):
        if W == 'DUMP':
            self._dump(words[i + 1] if i + 1 < len(words) else '*')
            return len(words) - i
        if W == 'DBEAM':
            return len(words) - i
        n = 1
        uni = cohe = add = False
        err = None
        while i + n < len(words):
            w = words[i + n].upper()
            if _abbrev(w, 'U_NIFORM'):
                uni = True
            elif _abbrev(w, 'R_ANDOM'):
                uni = False
            elif _abbrev(w, 'C_OHERENT'):
                cohe = True
            elif _abbrev(w, 'INC_OHERENT'):
                cohe = False
            elif _abbrev(w, 'A_DD'):
                add = True
            elif _abbrev(w, 'P_UT'):
                add = False
            else:
                v, k = self._getva(words, i + n)
                if v is None:
                    self.write(f'?Missing amount of error {W}\n')
                    return n
                err = float(v)
                n += k
                break
            n += 1
        if err is None:
            self.write(f'?Missing amount of error {W}\n')
            return n
        errg = 0.0
        while i + n < len(words):
            w = words[i + n]
            Wu = w.upper()
            if _abbrev(Wu, 'U_NIFORM'):
                uni = True
            elif _abbrev(Wu, 'R_ANDOM'):
                uni = False
            elif _abbrev(Wu, 'C_OHERENT'):
                cohe = True
                errg = err * self.rng.gauss(1, self.gcut())
            elif _abbrev(Wu, 'INC_OHERENT'):
                cohe = False
            elif _abbrev(Wu, 'A_DD'):
                add = True
            elif _abbrev(Wu, 'P_UT'):
                add = False
            else:
                if Wu in _COMMANDS or Wu in _ERRCMDS or Wu in _FUNS:
                    break
                if cohe and errg == 0.0:
                    errg = err * self.rng.gauss(1, self.gcut())
                hit = False
                matched = set(k for k in self.positions(w) if k < self.nlat)
                abb = any(c in w for c in '*%{}')
                if matched:
                    master = self._master_array()
                    for i1 in range(1, self.nlat):
                        if i1 not in matched:
                            continue
                        el_i = self.els[i1 - 1]
                        id_i = _TYPECODE.get(_cls(el_i), 0)
                        # terror.f: ドリフトは自分だけ、それ以外は master(i) = 同じ要素 (間にゼロ長要素だけ) の最後。
                        # master(i) <= 0 (グループの 2 番目以降 = スレーブ、マーカー等) は飛ばす
                        if id_i == 1:
                            ie = i1
                        else:
                            ie = master[i1]
                            if ie <= 0:
                                continue
                        # 乱数はグループ (先頭 i) ごとに 1 回だけ引く (スレーブには同じ誤差が入る)
                        if uni:
                            eg = err
                        elif not cohe:
                            eg = err * self.rng.gauss(1, self.gcut())
                        else:
                            eg = errg
                        hit = True
                        for j1 in range(i1, ie + 1):
                            if (id_i != 20 and master[j1] != 0) or id_i == 1:
                                self._apply_error(W, self.els[j1 - 1], eg, add)
                        if not abb:
                            break
                    if not uni:
                        self._tfadjst()
                if not hit:
                    self.write(f'?Undefined element {w}\n')
            n += 1
        self.invalidate()
        return n

    # ------------------------------------------------------------------ WAKE (tffsa.f:972)
    def _wake_cmd(self, words, i):
        """WAKE: K64 の tffsa.f は tfwake.f (Fortran の旧 WAKE 構文) を呼ばず、SADScript の WAKECOMMAND[] を評価する (init.n で
        Packages/WAKECOMMAND.n を AutoLoad。この配布には無い)。gs で確かめた: ファイルが無いと fileopen のエラーを出し、WAKECOMMAND[] は True にならず
        その行の残りを捨てる。ここも同じ (WAKECOMMAND.n があればそれが評価される)"""
        if _compat_v2() and i + 1 < len(words) and (_abbrev(words[i + 1].upper(), 'L_ONGITUDINAL')
                                                   or _abbrev(words[i + 1].upper(), 'T_RANSVERSE')):
            return self._wake_cmd_v2(words, i)
        try:
            ok = self.I.run('WAKECOMMAND[]', 'script')
        except SADError:
            ok = None
        if not (isinstance(ok, Sym) and ok.name == 'True'):
            return max(1, len(words) - i)
        return 1

    def _wake_cmd_v2(self, words, i):
        """v2 (B66): 旧式の WAKE コマンド (K64 の tfwake.f の構文。K64 では呼ばれない死んだコード)。
          WAKE L_ONGITUDINAL|T_RANSVERSE 要素 s1 w1 s2 w2 ... ENDWAKE
        名前 (ワイルドカード可) に合う CAVI に、表 {{s, w}, ...} を WakeFunction[Longitudinal|Transverse, 要素名] として登録する
        (追跡は WakeFunction と同じ: TRPT と LWAKE / TWAKE が ON のとき TrackParticles で効く。s: m、w: V/C (縦)、V/C/m (横))。
        表は行をまたいで書ける (tfwake.f と同じ: 1 行に何組でも、注釈 `!` 可。Interp._wake_table_open が ENDWAKE まで 1 文にする)。
        表の途中に `;` は書けない (tfwake.f も `;` で ?Missing ENDWAKE.)。ENDWAKE の後の `;` はあってもなくてもよい"""
        tag = 'Longitudinal' if _abbrev(words[i + 1].upper(), 'L_ONGITUDINAL') else 'Transverse'
        if i + 2 >= len(words) or words[i + 2].upper() == 'ENDWAKE':
            self.write('?Missing element for WAKE. (書き方: WAKE L|T 要素 s1 w1 s2 w2 ... ENDWAKE;)' + chr(10))
            return len(words) - i
        pat = words[i + 2]
        n = 3
        vals = []
        while True:
            if i + n >= len(words):
                self.write('?Missing ENDWAKE. (WAKE %s: 表は ENDWAKE で閉じる。表の途中に ; は書けない。%d 個の値まで読んだ)'
                           % (pat, len(vals)) + chr(10))
                return len(words) - i
            if words[i + n].upper() == 'ENDWAKE':
                n += 1
                break
            v, k = self._getva(words, i + n)
            if v is None or isinstance(v, list) or isinstance(v, complex):
                self.write(('?Missing value of wake.' if len(vals) % 2 else '?Missing ENDWAKE.')
                           + ' (WAKE %s: %d 個の値の後の "%s" は数として読めない。ENDWAKE を忘れていないか)'
                           % (pat, len(vals), words[i + n][:30]) + chr(10))
                return len(words) - i
            vals.append(float(v))
            n += k
        if len(vals) % 2:
            self.write('?Missing value of wake. (WAKE %s: 値が %d 個 (奇数)。s と w の組で書く)' % (pat, len(vals)) + chr(10))
            return n
        if not vals:
            self.write('?No wake data. (WAKE %s: WAKE と ENDWAKE の間に s w の組が無い)' % pat + chr(10))
            return n
        allp = [k for k in sorted(self.positions(pat)) if k < self.nlat]
        hit = [k for k in allp if type(self.els[k - 1]).__name__ == 'Cavity']
        if not hit:
            if allp:
                self.write('?WAKE: %s は CAVI ではない (%s)。旧式の WAKE は CAVI だけ'
                           % (pat, _cls(self.els[allp[0] - 1])) + chr(10))
            else:
                self.write('?Undefined element ' + pat + ' (WAKE は CAVI だけ)' + chr(10))
            return n
        tab = '{' + ','.join('{%r,%r}' % (vals[2 * j], vals[2 * j + 1]) for j in range(len(vals) // 2)) + '}'
        for k in hit:
            self.I.run('WakeFunction[%s,"%s"]=%s' % (tag, self.elname(k), tab), 'script')
        return n

    # ------------------------------------------------------------------ DELCOR (tcorr.f, tvcorr.f)
    def _delcor(self, words, i):
        """DELCOR [ADD] corlength delx dely element [element1 ...]: リングに沿った相関のあるランダムなずれ (tcorr.f)。
        corlength=0 は ランダムウォーク (両端を引き戻す)、≠0 は相関長 corlength のガウス過程 (tvcorr.f: 周波数空間で exp(-(k·2π cl/2L)²/4)
        の重みを付けた複素ガウス乱数を逆 FFT)。要素の中心の位置 s での値 (周期的に線形補間) × (delx, dely) を DX, DY に入れる (ADD は足す)"""
        n = 1
        add = False
        if i + 1 < len(words) and _abbrev(words[i + 1].upper(), 'A_DD'):
            add = True
            n += 1
        vals = []
        for _ in range(3):
            v, k = self._getva(words, i + n)
            if v is None:
                self.write('Syntax: DELCOR corlength delx dely element [element1...]\n')
                return max(n, len(words) - i)
            vals.append(float(v))
            n += k
        cl, dx, dy = vals
        els = self.els
        nlat = self.nlat
        pos = [0.0] * (nlat + 1)
        for k in range(1, nlat):
            pos[k + 1] = pos[k] + float(getattr(els[k - 1], 'L', 0.0) or 0.0)
        al = pos[nlat]
        if not self.flags['CELL']:
            al *= 2.0
        nc = 128
        while nc <= 32768 and nc * 0.5 < al:
            nc *= 2
        cv = self._tvcorr(cl, al, nc)
        master = self._master_array()
        design = lambda el, key: float(self.env._elements.get(el.name, {}).get(key, 0.0) or 0.0)
        exist1 = False
        while i + n < len(words):
            w = words[i + n]
            if w.upper() in _COMMANDS or w.upper() in _ERRCMDS:
                break
            abb = any(c in w for c in '*%{}')
            exist = False
            for i1 in sorted(self.positions(w)):
                if i1 >= nlat:
                    continue
                el_i = els[i1 - 1]
                id_ = _TYPECODE.get(_cls(el_i), 0)
                if id_ == 1 or (id_ > 8 and id_ not in (20, 31, 41)):
                    continue
                exist = True
                ie = i1 if id_ == 41 else master[i1]
                if ie <= 0:
                    continue
                s = 0.5 * (pos[i1] + pos[ie + 1])
                ns = (int(nc * s / al) + nc // 2) % nc + 1
                ns1 = ns % nc + 1
                wt = (nc * s / al + nc // 2 - (ns - 1)) % 1.0
                delx = (cv[ns - 1].real * (1.0 - wt) + cv[ns1 - 1].real * wt) * dx
                dely = (cv[ns - 1].imag * (1.0 - wt) + cv[ns1 - 1].imag * wt) * dy
                for j in range(i1, ie + 1):
                    if master[j] != 0 or id_ == 41:
                        ej = els[j - 1]
                        if id_ == 41:
                            # K64 のまま: MARK には latt+15, +17 に書く。値の添字では ky_DETR_Mark=15、ky_DPX_Mark=17 (DX=16、DY=18 ではない)。
                            # INS の始点が読む raw の DETR (R の行列式) と DPX (初期の傾き) に入る
                            # v2 (B40): 意図どおり MARK の DX・DY (ky_DX_Mark=16、ky_DY_Mark=18) に書く
                            raw = self.env._elements.setdefault(ej.name, {})
                            kx_, ky_ = ('DX', 'DY') if _compat_v2() else ('DETR', 'DPX')
                            ox = float(raw.get(kx_, 0.0) or 0.0)
                            oy = float(raw.get(ky_, 0.0) or 0.0)
                            raw[kx_] = ox + delx if add else delx
                            raw[ky_] = oy + dely if add else dely
                            continue
                        if id_ == 31 and _compat_v2():
                            # v2 (B40・B62): CAVI も他の磁石と同じく DX・DY に、同じ符号のずれ (+delx, +dely) を入れる。
                            # K64 は HARM・PHI に (設計の DX・DY) - (delx, dely) を書く (添字のずれ)
                            ej.DX = (float(getattr(ej, 'DX', 0.0) or 0.0) + delx) if add else delx
                            ej.DY = (float(getattr(ej, 'DY', 0.0) or 0.0) + dely) if add else dely
                            continue
                        if id_ == 31:
                            # K64 (tcorr.f): CAVI は latt+3, +4 = HARM と PHI に (DX, DY ではなく) 書く。add でなければ設計の DX・DY (idvalc+13, +14) から
                            # delx・dely を引いた値、add なら -delx・-dely を足す (K64 の添字のずれ。gs で確かめた: docs/21 B62)
                            dxo = float(getattr(ej, 'DX', 0.0) or 0.0)
                            dyo = float(getattr(ej, 'DY', 0.0) or 0.0)
                            if add:
                                ej._harm = float(getattr(ej, '_harm', 0) or 0) - delx          # HARM は (整数でなく) 実数のまま持つ
                                ej.PHI = float(getattr(ej, 'PHI', 0.0) or 0.0) - dely
                            else:
                                ej._harm = dxo - delx
                                ej.PHI = dyo - dely
                            continue
                        if not (hasattr(ej, 'DX') or id_ == 41):
                            continue
                        ej.DX = (float(getattr(ej, 'DX', 0.0) or 0.0) + delx) if add else delx
                        ej.DY = (float(getattr(ej, 'DY', 0.0) or 0.0) + dely) if add else dely
                if not abb:
                    break
            self._tfadjst()
            exist1 = exist1 or exist
            n += 1
            if not exist and not abb:
                if not exist1:
                    self.write('?Undefined element ' + w + '\n')
                break
        self.invalidate()
        return n

    def _tvcorr(self, x0, al, n):
        """tvcorr.f: 相関のある複素ガウス過程 cv (n 点)"""
        g = lambda: float(self.rng.gauss(1, self.gcut()))
        cv = np.zeros(n, dtype=complex)
        if x0 == 0.0:
            w = math.sqrt(al / n)
            for k in range(1, n):
                re_ = w * g()
                im_ = w * g()
                cv[k] = cv[k - 1] + complex(re_, im_)
            cd = cv[n - 1] / n
            for k in range(1, n):
                cv[k] = cv[k] - cd * k
            return cv
        x = x0 / al
        s = 0.0
        x1 = 0.25 * (x * 2.0 * math.pi) ** 2
        for k in range(2, n // 2 + 1):             # Fortran の添字 i = k (1 始まり)
            w = math.exp(-((k - 1) ** 2) * x1)
            re_ = w * g(); im_ = w * g()
            cv[k - 1] = complex(re_, im_)
            re_ = w * g(); im_ = w * g()
            if _compat_v2():
                # v2 (B62): 波数 +(k-1) の相手は -(k-1) (Fortran の添字 n-k+2)。K64 は cv(n-k+1) = 波数 -k に置き、
                # cv(n) (波数 -1) が未設定 (ktaloc の古い内容) で、2 回目以降の DELCOR が再現できない
                cv[n - k + 1] = complex(re_, im_)
            else:
                cv[n - k] = complex(re_, im_)
            s += w * w
        cv[0] = complex(g(), g())
        if _compat_v2():
            # v2 (C2): 実部 (DX に使う)・虚部 (DY に使う) それぞれの分散を 1 にする。Var(Re) = Σ_k w_k² = 2Σw² + 1。
            # K64 は複素数の絶対値の 2 乗平均を 1 にする (4Σw² + 2) ので、DX の rms が delx/√2 になる
            s = 1.0 / math.sqrt(2.0 * s + 1.0)
        else:
            s = 1.0 / math.sqrt(4.0 * s + 2.0)
        cv = np.fft.ifft(cv) * n * s           # tcftr (conj=.false.) = exp(+2πi jk/n) の和
        return cv

    def _master_array(self):
        """tmast.f: master[i] (1 始まり、nlat-1 まで)。磁石 (BEND/QUAD/SEXT/OCT/DECA/DODECA/MULT/CAVI/TCAV) の先頭 i は
        同じ要素 (名前が同じ) の最後の位置 ie、同じ要素の 2 番目以降 j は -i、それ以外 (ドリフト・ソレノイド・マーカー・
        別の要素) は 0。同じ要素の間には長さ 0 の別要素 (マーカー・L=0 の BEND 等) だけを許す"""
        n = self.nlat
        master = [0] * (n + 1)
        els = self.els
        ie = 0
        for i in range(1, n):
            if i <= ie:
                continue
            id_ = _TYPECODE.get(_cls(els[i - 1]), 0)
            if id_ == 20:
                continue
            if id_ == 1 or (id_ > 12 and id_ not in (22, 31, 32)):
                continue
            ie = i
            name = els[i - 1].name
            for j in range(i + 1, n):
                if els[j - 1].name == name:
                    ie = j
                    master[j] = -i
                elif float(getattr(els[j - 1], 'L', 0.0) or 0.0) != 0.0:
                    break
            master[i] = ie
        return master

    def _tfadjst(self):
        """tfadjst.f: 隙間なく並んだ磁石 (BEND/QUAD/SEXT/OCT/MULT/CAVI。間にゼロ長要素だけ) の鎖は、鎖の先頭の磁石の DX/DY の
        ずれ (設計値からの差) を鎖の全磁石が共有する"""
        n = self.nlat
        els = self.els
        pos = [0.0] * (n + 1)
        for k in range(1, n):
            pos[k + 1] = pos[k] + float(getattr(els[k - 1], 'L', 0.0) or 0.0)
        magnet = (2, 4, 6, 8, 22, 31)

        def mag(k):
            return _TYPECODE.get(_cls(els[k - 1]), 0) in magnet and hasattr(els[k - 1], 'DX')

        def design(el, key):
            return float(self.env._elements.get(el.name, {}).get(key, 0.0) or 0.0)

        i0 = 0
        for i in range(1, n):
            if i > i0 and mag(i):
                el = els[i - 1]
                s0 = pos[i]
                sx = float(el.DX or 0.0) - design(el, 'DX')
                sy = float(el.DY or 0.0) - design(el, 'DY')
                for j in range(i, n):
                    if pos[j] != s0:
                        break
                    if mag(j):
                        i0 = j
                        s0 = pos[j + 1]
                        ej = els[j - 1]
                        ej.DX = sx + design(ej, 'DX')
                        ej.DY = sy + design(ej, 'DY')

    def _apply_error(self, W, el, eg, add):
        tn = type(el).__name__
        raw = self.env._elements.get(el.name, {})
        if W == 'DK':
            key = self.default_key(el)
            if key is None or tn in ('Marker', 'Solenoid', 'Drift'):
                return
            design = self.get_value(el, key)
            el._errk = (el._errk if add else 1.0) + eg
            self.set_value(el, key, design)
        elif W == 'DL':
            el.L = max(0.0, (el.L if add else float(raw.get('L', el.L) or 0.0)) + eg)
        elif W == 'DDK':
            # terror.f: K0 のキーワードがある要素 (BEND・MULT・TCAVI) だけ。QUAD などの K0 属性 (pysad の内部) には入れない
            from .sadkw import KEYWORDS
            if hasattr(el, 'K0') and 'K0' in KEYWORDS.get(_TYPENAME.get(tn, ''), ()):
                el.K0 = (el.K0 if add else float(raw.get('K0', 0.0) or 0.0)) + eg
        elif W == 'DELX':
            if hasattr(el, 'DX'):
                el.DX = (el.DX if add else float(raw.get('DX', 0.0) or 0.0)) + eg
        elif W == 'DELY':
            if hasattr(el, 'DY'):
                el.DY = (el.DY if add else float(raw.get('DY', 0.0) or 0.0)) + eg
        elif W == 'DTHETA':
            if tn in ('Marker', 'Monitor', 'Drift'):
                return
            # terror.f: DTHETA は DROT のある要素 (BEND, MULT) では DROT、なければ ROT (QUAD, SEXT...) に入る
            if tn in ('Bend', 'Multipole'):
                attr, kw = 'DROTATE', 'DROTATE'
            else:
                attr, kw = 'ROTATION', 'ROTATE'
            base = float(getattr(el, attr, 0.0) or 0.0) if add else float(raw.get(kw, 0.0) or 0.0)
            setattr(el, attr, base + eg)

    def _dump(self, pat):
        self.write('  element     DL(mm)    DK          DDK      DTHETA(mrad) DELX(mm) DELY(mm)\n')
        pat = pat[:20]                  # terror.f: word1 は character*20 (長い `A|B|C…` は 20 文字で切れる)
        if not self._dump_pass(pat):
            self._dump_pass('*')        # どの要素にも合わなければ全要素 (all=.true.)

    def _dump_pass(self, pat):
        exist = False
        for k in range(1, self.nlat):
            el = self.els[k - 1]
            if not _wild_match(el.name, pat) and not _wild_match(self.elname(k), pat):
                continue
            exist = True
            raw = self.env._elements.get(el.name, {})
            dl = el.L - float(raw.get('L', el.L) or 0.0)
            dk = getattr(el, '_errk', 1.0) - 1.0
            ddk = float(getattr(el, 'K0', 0.0) or 0.0) - float(raw.get('K0', 0.0) or 0.0)
            if type(el).__name__ in ('Bend', 'Multipole'):
                dth = float(getattr(el, 'DROTATE', 0.0) or 0.0) - float(raw.get('DROTATE', 0.0) or 0.0)
            else:
                dth = float(getattr(el, 'ROTATION', 0.0) or 0.0) - float(raw.get('ROTATE', 0.0) or 0.0)
            dx = float(getattr(el, 'DX', 0.0) or 0.0) - float(raw.get('DX', 0.0) or 0.0)
            dy = float(getattr(el, 'DY', 0.0) or 0.0) - float(raw.get('DY', 0.0) or 0.0)
            if any((dl, dk, ddk, dth, dx, dy)):
                self.write(f' {self.elname(k):10s}{dl * 1e3:10.5f}{_gwd(dk, 12, 5)}{_gwd(ddk, 12, 5)}{dth * 1e3:10.5f}{dx * 1e3:10.5f}{dy * 1e3:10.5f}\n')
        return exist

    # ------------------------------------------------------------------ COUPLE
    def _couple(self, words, i):
        if i + 3 > len(words) - 1 + 1 and len(words) - i < 4:
            if len(words) - i == 2 and words[i + 1] == '*':
                return 2
            self.write('Missing coefficient for COUP_LE.\n')
            return len(words) - i
        slave, master = words[i + 1], words[i + 2]
        v, k = self._getva(words, i + 3)
        if v is None:
            self.write('Missing coefficient for COUP_LE.\n')
            return 3
        sf = self.families(slave)
        mf = self.families(master)
        if not sf or not mf:
            self.write(f'Undefined element for COUP_LE: {slave} {master}\n')
            return 4
        self.couples[sf[0]] = (mf[0], float(v))
        if self.matcher:
            self.matcher.variables = [x for x in self.matcher.variables if x['family'] != sf[0]]
        key = self.default_key(self.instances(mf[0])[0])
        if key:
            self.set_family(mf[0], key, self.get_value(self.instances(mf[0])[0], key))
        return 4

    # ------------------------------------------------------------------ TYPE / DISP / SAVE / RESET
    # tftype.f ltyp の順 (型のグループの出力順)
    _TYPE_ORDER = ['DRIFT', 'BEND', 'QUAD', 'SEXT', 'OCT', 'DECA', 'DODECA', 'MULT', 'SOL', 'CAVI', 'TCAVI',
                   'MAP', 'BEAM', 'SPCH', 'INS', 'COORD', 'APERT', 'MONI', 'MARK']

    def _type_value(self, el, tname: str, key: str, l: int):
        """tftype.f tftyp1 の値: MARK の光学関数は先頭要素 (l=1) なら要素自身の値、それ以外は直前の CALC の twiss。
        EMITX/EMITY/EMITZ は emx/emy/emz、SIGMAZ は sigzs、SIGE は sizedp、DP は dpmax。他は要素の値 (誤差なしの設計値)"""
        try:
            if tname == 'MARK':
                key = {'PSIX': 'NX', 'PSIY': 'NY', 'PSIZ': 'NZ'}.get(key, key)   # MARK のキーワード名 → 光学関数名
                if key in _TWISS_ALL:
                    if l == 1:
                        return self.get_value(el, key)
                    if self.rows is None:
                        return 0.0                       # CALC 前は twiss 配列が 0 (SAD)
                    try:
                        return self.twiss_value(key, l)
                    except Exception:
                        return 0.0
                if key in ('EMITX', 'EMITY', 'EMITZ'):
                    return float(self.I.get(key, 0.0) or 0.0)
                if key == 'SIGMAZ':
                    return float(self.I.get('SIGZ', 0.0) or 0.0)
                if key == 'SIGE':
                    return float(self.I.get('SIGE', 0.0) or 0.0)
                if key == 'DP':
                    return self.dp()
            return self.get_value(el, key)
        except Exception:
            return None

    def type_cmd(self, pats):
        """TYPE [pattern ...]: SAD の tftype の書式。パターンごとに、型 (tftype.f ltyp 順) のグループの前に ' ;' を出し、
        値が 0 でないキーワード (と既定の変数 ival: BEND の ANGLE、QUAD の K1、DRIFT の L など。これは 0 でも出す) を
        SAD のキーワード順 (sadkw.KEYWORDS) で並べる。ROTATE は度で ' DEG' を付ける。値は最短表現、値の後は 5 の倍数の幅まで空白
        (最後は ' )')、全部の後に ' ;'"""
        from .sadkw import KEYWORDS, DEFAULT_KEY
        if not pats:
            pats = ['*']
        for pat in pats:
            groups = {}
            for fam in self.families(pat):
                el = self.instances(fam)[0]
                groups.setdefault(_TYPENAME.get(_cls(el), 'DRIFT'), []).append(fam)
            lpw = 131 - 2 + 2                    # itfgetrecl の幅 (twbuf の右端は lpw、値は lpw-2)
            lpw = 131
            st = {'buff': '', 'ip': 0}

            def twbuf(word, unit, ilmgn, irmgn, itab, icmd):
                """twbuf.f の直訳 (行の組み立て。単語はタブ位置 itab の倍数に置く。幅を超えたら改行)"""
                if icmd == 0:
                    st['ip'] = 0
                    return
                if icmd == -1:
                    if st['ip'] > ilmgn - 1:
                        self.write(st['buff'][:st['ip']] + chr(10))
                    st['ip'] = 0
                    return
                while True:
                    ip = st['ip']
                    ip1 = (max(ip, ilmgn - 1) + itab - 1) // itab * itab
                    if ip1 >= irmgn:
                        self.write(st['buff'][:ip] + chr(10))
                        st['ip'] = 0
                        continue
                    lw = len(word.rstrip())
                    ln = lw + len(unit) + 1
                    if ln > 1:
                        if ip1 + ln > irmgn:
                            self.write(st['buff'][:ip] + chr(10))
                            st['ip'] = 0
                            continue
                        buf = st['buff'][:ip].ljust(ip1 + 1) + word[:lw] + unit
                        st['buff'] = buf
                        st['ip'] = ip1 + ln
                    return

            for tname in self._TYPE_ORDER:
                fams = groups.get(tname)
                if not fams:
                    continue
                mulc = tname != 'DRIFT'              # kytbl(kwMAX, lt) > ky_MAX_DRFT: DRIFT 以外は 1 要素ごとに改行、DRIFT は詰めて並べる
                self.write(' ;' + chr(10))
                twbuf(' ', '', 1, lpw, 0, 0)
                twbuf(tname, '', 1, lpw, 7, 1)
                for fam in fams:
                    el = self.instances(fam)[0]
                    l = self.els.index(el) + 1
                    ival = DEFAULT_KEY.get(tname, '')
                    pairs = []
                    for k in KEYWORDS.get(tname, ['L']):
                        v = self._type_value(el, tname, k, l)
                        if v is None or not is_num(v):
                            continue
                        v = float(v)
                        if v != 0.0 or k == ival:
                            if k == 'ROTATE':
                                deg = float(f'{v * 180.0 / math.pi:.14g}')
                                pairs.append((k, deg, ' DEG'))
                            else:
                                pairs.append((k, v, ''))
                    nm = f'{fam:8s}'
                    if not pairs:
                        twbuf(nm, '=()', 7, lpw, 7, 1)
                    else:
                        for n, (k, v, unit) in enumerate(pairs):
                            vout = f'{k} ={_shortest(v)}'
                            if n == 0:
                                twbuf(nm + '=(' + vout, unit.rstrip(), 7, lpw - 2, 7, 1)
                            else:
                                twbuf(vout, unit.rstrip(), 10, lpw - 2, 5, 1)
                        twbuf(')', '', 10, lpw, 1, 1)
                    if mulc:
                        twbuf(' ', '', 10, lpw, 0, -1)
                if not mulc:
                    twbuf(' ', '', 10, lpw, 0, -1)
        self.write(' ;' + chr(10))

    # ------------------------------------------------------------------ DISP (tfdisp.f)
    def _eval_real(self, text: str):
        """式を評価して実数ならその値、そうでなければ None (tfevalb / tfeval の ktfrealq)"""
        try:
            v = self.I.run(text, 'script')
        except Exception:
            return None
        return float(v) if is_num(v) else None

    def _ielma(self, i: int) -> int:
        """tffs.f ielma (trpt=.false.): 1..nlat に収める (0 以下は nlat+i+1)"""
        n = self.nlat
        return min(n, i) if i > 0 else max(1, n + i + 1)

    def _ielmh(self, name: str, iord: int):
        """ielm.f ielmh: iord=0 は 1 台だけの要素かファミリーの代表 (最初)、iord=n>0 は同名が複数あるときの n 番目 (1 台だけの要素の NAME.1 は見つからない)"""
        pos = [k + 1 for k, el in enumerate(self.els) if el.name.upper() == name.upper()]
        if not pos:
            return None
        if iord == 0:
            return pos[0]
        if len(pos) > 1 and 1 <= iord <= len(pos):
            return pos[iord - 1]
        return None

    def _ielmf(self, word: str):
        """ielm.f ielmf: 位置の記法 `NAME`, `NAME.n`, `NAME.n+2`, `$$$`, `***`, `^^^`。戻り (位置 or None, 未定義のときのメッセージ or None)"""
        lw = len(word)
        d = word.find('.')
        ip, im = word.find('+'), word.find('-')
        cands = [x for x in (ip, im) if x >= 0]
        ipm = min(cands) if cands else -1
        if d >= 0:
            if 0 <= ipm < d:
                ln, iord = ipm, 0
            else:
                v = self._eval_real(word[d + 1:ipm if ipm >= 0 else lw])
                if v is None:
                    return None, None
                ln, iord = d, int(v)
        else:
            iord = 0
            ln = ipm if ipm >= 0 else lw
        name = word[:ln]
        ioff = 0
        if ipm >= 0:
            v = self._eval_real(word[ipm:lw])
            if v is None:
                return None, None
            ioff = int(v)
        if name in ('$$$', '***'):
            k = self.nlat
        elif name == '^^^':
            k = 1
        else:
            k = self._ielmh(name, iord)
            if k is None:
                return None, ('?Undefined location ' + word if d >= 0 else None)
        return self._ielma(k + ioff), None

    def _ielmex(self, word: str):
        """ielm.f ielmex: 数式ならその値 (v>=0: int(v+0.499)、負: nlat+1+v+0.5) を 1..nlat に収めた位置、そうでなければ名前。戻り (位置 or None, メッセージ or None)"""
        if not word.strip():
            return None, None
        if word[0] != '^':
            v = self._eval_real(word)
            if v is not None:
                iv = int(v + 0.499) if v >= 0.0 else int(self.nlat + 1 + v + 0.5)
                return max(1, min(self.nlat, iv)), None
        return self._ielmf(word)

    def _temat(self, k: int, pat: str) -> bool:
        """temat.f (パターンのとき): 標準名 (NAME / NAME.n) に合う、または 1 台だけの要素なら NAME.1 に合う、複数ある要素の代表 (最初) ならファミリー名に合う"""
        if _wild_match(self.elname(k), pat):
            return True
        if k == self.nlat:
            return False
        el = self.els[k - 1]
        inst = self.instances(el.name)
        if len(inst) == 1:
            return _wild_match(el.name + '.1', pat)
        if self.els.index(inst[0]) == k - 1:
            return _wild_match(el.name, pat)
        return False

    _DISP_MODES = (('A_CCELERATION', 'A'), ('G_EOMETRY', 'G'), ('GA_MMA', 'GA'), ('OG_EOMETRY', 'OG'), ('O_RBIT', 'O'),
                   ('R_MATRIX', 'R'), ('B_EAM', 'B'), ('P_HYSICAL', 'P'), ('D_UMPOPTICS', 'D'))

    def disp_cmd(self, words, i):
        """DISP [パターン] [モード] [loc1 [loc2]]: tfdisp.f の直訳 (表示範囲)。
        語を 1 つずつ読み、パターン (| < { * % を含む) は表示する要素名の絞り込み (temat)、ALL は範囲 1..nlat、モード語 (R_MATRIX O_RBIT E_XTREMUM ...)、
        それ以外は位置 (数式・NAME・NAME.n・NAME.n+2・$$$・^^^) として 1 つ目を loc1、続く 1 語が位置なら loc2 (そうでなければ loc2=loc1)。
        範囲 (id1, id2) は FFS に保持され (USE で 1..nlat に戻る)、位置を省くと前回の範囲が使われる。
        loc1==loc2 はその要素から $$$ まで (刻み nlat-loc2)、loc2<loc1 は nlat を超えて周回する。
        戻り値: 消費した語数 (DISP 自身を含む)。位置でない語に出会ったらそこで表示せずに止め、その語は次のコマンドとして読まれる。"""
        nlat = self.nlat
        id1, id2 = getattr(self, '_disp_id1', 1), getattr(self, '_disp_id2', nlat)
        seldis = False
        word1 = None
        mdisp = None
        dpeak = False
        icolm, dref = 0, False
        n = 1
        while True:
            w = words[i + n] if i + n < len(words) else ' '
            W = w.upper()
            n += 1
            if w != ' ' and W[:3] != '***' and any(c in W for c in '|<{*%'):
                seldis, word1 = True, W
                continue
            if W == 'ALL':
                id1, id2 = 1, nlat
                continue
            hit = False
            for full, code in self._DISP_MODES:
                if _abbrev(W, full):
                    mdisp = code
                    hit = True
                    break
            if hit:
                continue
            if W == 'Z':
                mdisp = 'Z'
                continue
            if _abbrev(W, 'E_XTREMUM'):
                dpeak = True
                continue
            if _abbrev(W, 'RE_FERENCE'):
                icolm = -1
                continue
            if _abbrev(W, 'DRE_FERENCE'):
                dref = True
                continue
            break
        k1, msg = self._ielmex(W)
        if k1 is not None:
            id1 = k1
            w2 = words[i + n] if i + n < len(words) else ' '
            k2, msg2 = self._ielmex(w2.upper())
            if k2 is not None:
                id2 = k2
                n += 1
            else:
                id2 = id1
                if msg2:
                    self.write(' ???-FFS-Error-' + msg2 + chr(10) + 'Out[1]:= USE' + chr(10))
                    return len(words) - i
        elif w.strip():
            if msg:
                self.write(' ???-FFS-Error-' + msg + chr(10) + 'Out[1]:= USE' + chr(10))
                return len(words) - i
            return n - 1                  # 位置でない語: 表示せずに戻り、その語は次のコマンド
        else:
            n -= 1                         # 文の終わり (位置なし): 前回の範囲で表示
        self._disp_id1, self._disp_id2 = id1, id2
        if mdisp in ('G', 'OG'):
            self._disp_geo(id1, id2, mdisp, seldis, word1)
            return n
        self._disp_show(id1, id2, mdisp, seldis, word1, dpeak, icolm, dref)
        return n

    def _disp_geo(self, id1, id2, mdisp, seldis, word1):
        """DISP G_EOMETRY / OG_EOMETRY (tfdisp.f:322-370): Element, Gx Gy Gz (12.6), s (12.6), Length (10.6), Value (10.7),
        Chi1-3 (度, 12.6), #。OG は閉軌道 (DX DPX DY DPY) を乗せた幾何 (tforbitgeo)。ヘッダは 66 行ごと"""
        from .geo import tforbitgeo, tfchi
        nlat = self.nlat
        rows = self.ensure_rows()
        g = self.ensure_geo()
        if id1 == id2:
            id3, idstep = nlat, max(nlat - id2, 1)
        elif id2 < id1:
            id3, idstep = id2 + nlat, 1
        else:
            id3, idstep = id2, 1
        hc = ' ' if mdisp == 'G' else 'O'
        header = (' Element             ' + hc + 'Gx        ' + hc + 'Gy          ' + hc + 'Gz    '
                  + '        s        Length     Value' + '       ' + hc + 'Chi1       ' + hc + 'Chi2       ' + hc + 'Chi3      #' + chr(10))
        f = fmt_real

        def tz(v, form):
            s_ = f(v, form)
            if v == 0.0:
                d = s_.find('.')
                if d != -1 and d < len(s_) - 2:
                    s_ = s_[:d + 2] + ' ' * (len(s_) - d - 2)
            return s_

        deg = math.pi / 180.0
        lines = 0
        for l1 in range(id1, id3 + 1, idstep):
            k = (l1 - 1) % nlat + 1
            if seldis and not self._temat(k, word1):
                continue
            if lines % 66 == 0:
                self.write(header)
            r = rows[k - 1]
            if k < nlat:
                el = self.els[k - 1]
                key = self.default_key(el)
                tn = type(el).__name__
                if key and tn not in ('Marker', 'Monitor'):
                    val = self.get_value(el, key)
                    vout = ' ' + tz(val, '10.7')
                else:
                    vout = ' 0'
                L = el.L
                name = self.elname(k)
                rev = bool(getattr(el, '_rev', False))
            else:
                vout, L, name, rev = ' 0', 0.0, '$$$', False
            gk = g.geo[k - 1]
            if mdisp == 'OG':
                gk = tforbitgeo(gk, (r.x_cod, r.px_cod, r.y_cod, r.py_cod))
            buff = ''.join(tz(gk[i, 3], '12.6') for i in range(3))
            buff += tz(float(g.pos[k - 1]), '12.6') + tz(L, '10.6') + vout.ljust(11)
            buff += ''.join(tz(tfchi(gk, i) / deg, '12.6') for i in (1, 2, 3))
            self.write(('-' if rev else ' ') + name.ljust(18) + buff[:105] + f' {(k - 1) % 10000 + 1:5d}' + chr(10))
            lines += 1

    def _disp_val(self, row, key):
        """DISP の 1 行の値 (twiss(l,icol,mf))。GMX/GMY/GMZ は (1+α²)/β (tgetgm)"""
        if key == 'GMX':
            return (1.0 + row.alfx ** 2) / row.betx
        if key == 'GMY':
            return (1.0 + row.alfy ** 2) / row.bety
        if key == 'GMZ':
            return (1.0 + row.alfz ** 2) / row.betz if row.betz else 0.0
        return float(getattr(row, self._TWKEY[key], 0.0))

    def _disp_show(self, id1, id2, mdisp, seldis, word1, dpeak, icolm=0, dref=False):
        """tfdisp.f の表。mdisp: None (標準)・R・O・P・GA・Z・A・D・B。icolm=-1 は REFERENCE (参照光学)、dref は DREFERENCE (現在 - 参照、
        BX BY BZ GMX GMY GMZ は比)"""
        nlat = self.nlat
        rows = self.ensure_rows()
        ref = getattr(self, '_ref_rows', None)
        if (icolm == -1 or dref) and ref is None:
            ref = [_ZeroRow() for _ in rows]
        if mdisp == 'B' and (dref or icolm != 0):
            self.write(' ???-FFS-Info-REF and DREF not implemented for DISP B. ' + chr(10))
            return
        if id1 == id2:
            id3, idstep = nlat, max(nlat - id2, 1)
        elif id2 < id1:
            id3, idstep = id2 + nlat, 1
        else:
            id3, idstep = id2, 1
        # tfdisp.f:398 `if(mod(lines,nlc) == 0) write header` (nlc=66): ヘッダは 66 行ごとに繰り返す
        left, mid, right = ('   AX      BX      NX      EX      EPX  ', ' Element   Length   Value      s(m)    ',
                            '   AY      BY      NY      EY      EPY    DetR     #')
        if mdisp == 'R':
            mid = ' Element    R1     R2     R3     R4    '
        elif mdisp == 'P':
            mid = ' Element    PEX    PEPX   PEY    PEPY  '
        elif mdisp == 'GA':
            left = '   AX      BX      NX      GMX     GMX*L'
            right = '   AY      BY      NY      GMY     GMY*L  DetR     #'
        elif mdisp == 'O':
            mid = ' Element    DX     DPX    DY     DPY   '
        elif mdisp == 'B':
            left = '   AXp     BXp    EMITXp   EXp     EPXp '
            mid = ' Element   Sigx(mm)  Sigy(mm)  Rot(deg)'
            right = '   AYp     BYp    EMITYp   EYp     EPYp   sigp     #'
        elif mdisp == 'A':
            mid = ' Element    p(GeV)emitx(m)emity(m) DDP '
            right = '   AY      BY      NY      EY      EPY     DZ      #'
        elif mdisp == 'Z':
            left = '   AZ      BZ      NZ      DZ     DDP   '
            right = '   ZX      ZPX     ZY      ZPY     GMZ    DetR     #'
        header = left + mid + right
        if (dref or icolm == -1) and mdisp != 'B':
            h = list(header.ljust(131))

            def put(pos, txt):
                h[pos - 1:pos - 1 + len(txt)] = list(txt)
            if dref:
                for pos in (3, 11, 19, 27, 35, 82, 90, 98, 106, 114, 121):
                    put(pos, 'd')
                put(14, '/BX')
                put(93, '/BY')
                if mdisp == 'Z':
                    put(14, '/BZ')
                    put(93, 'X  ')
                    put(114, ' ')
                elif mdisp in ('O', 'R', 'P'):
                    for pos in (52, 59, 66, 73):
                        put(pos, 'd')
                elif mdisp == 'A':
                    put(121, ' d')
            else:
                for pos in (6, 14, 22, 30, 39, 85, 93, 101, 106, 118, 126):
                    put(pos, 'R')
                if mdisp == 'Z':
                    put(118, ' ')
                elif mdisp in ('O', 'R', 'P'):
                    for pos in (55, 63, 69, 77):
                        put(pos, 'R')
                elif mdisp == 'A':
                    put(125, 'R ')
            header = ''.join(h).rstrip()
        f = fmt_real

        def tz(v, form):
            """tdtrimz: 値が 0 なら '.' の 1 桁後ろから空白 (`.00000` → `.0`)"""
            s = f(v, form)
            if v == 0.0:
                d = s.find('.')
                if d != -1 and d < len(s) - 2:
                    s = s[:d + 2] + ' ' * (len(s) - d - 2)
            return s

        cur = ref if icolm == -1 else rows

        def tdv(k, key, form, a=None):
            """tdispv"""
            sc = PI2 if key in ('NX', 'NY', 'NZ') else 1.0
            if dref:
                c, r0 = self._disp_val(rows[k - 1], key), self._disp_val(ref[k - 1], key)
                v = ((c - r0) / r0 if r0 else float('nan')) if key in ('BX', 'BY', 'BZ', 'GMX', 'GMY', 'GMZ') else c - r0
            else:
                v = self._disp_val(cur[k - 1], key) / sc
                if a is not None:
                    v = v * a
            return tz(v, form)

        def bxby(k):
            r = cur[min(k, nlat) - 1]
            return r.betx, r.bety

        if dpeak:
            bx0, by0 = bxby(max(1, id1 - 1))
            bx1, by1 = bxby(id1)
        lines = 0
        for l1 in range(id1, id3 + 1, idstep):
            k = (l1 - 1) % nlat + 1
            if seldis and not self._temat(k, word1):
                continue
            if dpeak:
                bx2, by2 = bxby(min(k + 1, nlat))
                peak = ((bx1 >= bx0 and bx1 >= bx2) or (by1 >= by0 and by1 >= by2)
                        or (bx1 <= bx0 and bx1 <= bx2) or (by1 <= by0 and by1 <= by2))
                bx0, by0, bx1, by1 = bx1, by1, bx2, by2
                if not peak:
                    continue
            r = rows[k - 1]
            key = None
            if k < nlat:
                el = self.els[k - 1]
                key = self.default_key(el)
                val = self.get_value(el, key) if key else 0.0
                L = el.L
                name = self.elname(k)
                rev = bool(getattr(el, '_rev', False))
            else:
                val, L, name = 0.0, 0.0, '$$$'
                rev = False
            if mdisp == 'D':
                # moded: 最初の行の前に scale の行、各行は 名前 (12 桁)・型番号 (I3)・28 個の光学関数 (12.8、scale で割る)
                if lines == 0:
                    self.write((' ' * 15 + ''.join(tz(PI2 if kk in ('NX', 'NY', 'NZ') else 1.0, '12.8').ljust(12) for kk in _TWISS_ALL)).rstrip() + chr(10))
                idt = 0 if k >= nlat else int(_TYPECODE.get(_cls(self.els[k - 1]), 0))
                vals = ''.join(tz(self._disp_val(cur[k - 1], kk) / (PI2 if kk in ('NX', 'NY', 'NZ') else 1.0), '12.8').ljust(12)
                               for kk in _TWISS_ALL)
                self.write((name.ljust(12) + f'{idt:3d}' + vals).rstrip() + chr(10))
                lines += 1
                continue
            if lines % 66 == 0:
                self.write(header + chr(10))
            # tfdisp.f:438-448: 反転配置 (cmp%ori=.false.、LINE の `-`) の要素は ' -'//name (行頭側の空白 2 つのうち 2 つ目が `-` になる)。
            # LINE["NAME"] / Twiss["NAME"] は elname で名前だけを返す (`-` は付かない)
            #   forward: lname<=8 → '  '//name を 10 桁、それより長い名前は ' '//name。reversed: ' -'//name を max(10, lname+2) 桁
            if rev:
                bname = (' -' + name).ljust(max(10, len(name) + 2))
            elif len(name) <= 8:
                bname = ('  ' + name).ljust(10)
            else:
                bname = ' ' + name
            if mdisp == 'Z':
                head = ''.join(tdv(k, kk, '8.5') for kk in ('AZ', 'BZ', 'NZ', 'DZ', 'DDP'))
            elif mdisp == 'GA':
                head = ''.join(tdv(k, kk, '8.5') for kk in ('AX', 'BX', 'NX', 'GMX')) + tdv(k, 'GMX', '8.5', L)
            elif mdisp == 'B':
                B = self.beam_matrix_twiss(k)
                sigpp = B[5, 5]
                etaxp, etapxp = B[0, 5] / sigpp, B[1, 5] / sigpp
                sxx, sxpx, spxpx = B[0, 0] - etaxp ** 2 * sigpp, B[0, 1] - etaxp * etapxp * sigpp, B[1, 1] - etapxp ** 2 * sigpp
                emixp = math.sqrt(sxx * spxpx - sxpx ** 2)
                head = ''.join((tz(-sxpx / emixp, '8.5'), tz(sxx / emixp, '8.5'), tz(emixp, '8.5'), tz(etaxp, '8.5'), tz(etapxp, '8.5')))
            else:
                head = ''.join(tdv(k, kk, '8.5') for kk in ('AX', 'BX', 'NX', 'EX', 'EPX'))
            if mdisp == 'R':
                mid_s = ''.join(tdv(k, kk, '7.4') for kk in ('R1', 'R2', 'R3', 'R4'))
            elif mdisp == 'O':
                mid_s = ''.join(tdv(k, kk, '7.4') for kk in ('DX', 'DPX', 'DY', 'DPY'))
            elif mdisp == 'P':
                if dref:
                    pe = [self._disp_val(rows[k - 1], kk) - self._disp_val(ref[k - 1], kk) for kk in ('EX', 'EPX', 'EY', 'EPY')]
                else:
                    pe = [self._disp_val(cur[k - 1], kk) for kk in ('EX', 'EPX', 'EY', 'EPY')]
                mid_s = ''.join(f(v, '7.4') for v in pe)
            elif mdisp == 'B':
                mid_s = (tz(math.sqrt(B[0, 0]) * 1e3, '10.8') + tz(math.sqrt(B[2, 2]) * 1e3, '10.8')
                         + tz(math.atan2(-2.0 * B[0, 2], B[0, 0] - B[2, 2]) * 90.0 / math.pi, '9.4'))
            elif mdisp == 'A':
                amass = float(self.I.get('MASS', 0.51099895e6) or 0.51099895e6)
                rr = self._gammab_at(1) / self._gammab_at(k)
                emx = float(self.I.get('EMITX', 0.0) or 0.0)
                emy = float(self.I.get('EMITY', 0.0) or 0.0)
                mid_s = (tz(self._gammab_at(k) * amass / 1e9, '8.6') + tz(emx * rr, '7.5') + tz(emy * rr, '7.5')
                         + tz(cur[k - 1].delta_cod, '7.4'))
            else:
                # Length は 8.5 (8 桁)、Value は vout=10.7 (10 桁。負の値も 7 桁出る)、s(m) は 11.6 (tfdisp.f の buff(51:58), (59:68), (69:79))
                ls = '  .0    ' if L == 0.0 else f(L, '8.5')
                # tfdisp.f:293-302: 既定キー (nelvx%ival) のある要素は tdtrimz(値,'10.7') (0 でも ' .0')、MARK/キー無しは ' 0'
                vs = ' 0        ' if (k == nlat or not key) else tz(val, '10.7').ljust(10)
                mid_s = f'{ls}{vs}{tz(r.s, "11.6")}'
            mid_s = mid_s.ljust(29) if mdisp in ('R', 'O', 'P') else mid_s
            if mdisp == 'Z':
                tail_cols = ''.join(tdv(k, kk, '8.5') for kk in ('ZX', 'ZPX', 'ZY', 'ZPY', 'GMZ')) + tdv(k, 'DETR', '7.4')
            elif mdisp == 'GA':
                tail_cols = ''.join(tdv(k, kk, '8.5') for kk in ('AY', 'BY', 'NY', 'GMY')) + tdv(k, 'GMY', '8.5', L) + tdv(k, 'DETR', '7.4')
            elif mdisp == 'B':
                etayp, etapyp = B[2, 5] / sigpp, B[3, 5] / sigpp
                syy, sypy, spypy = B[2, 2] - etayp ** 2 * sigpp, B[2, 3] - etayp * etapyp * sigpp, B[3, 3] - etapyp ** 2 * sigpp
                emiyp = math.sqrt(syy * spypy - sypy ** 2)
                tail_cols = (tz(-sypy / emiyp, '8.5') + tz(syy / emiyp, '8.5') + tz(emiyp, '8.5') + tz(etayp, '8.5') + tz(etapyp, '8.5')
                             + tz(math.sqrt(sigpp), '7.4'))
            else:
                tail_cols = (''.join(tdv(k, kk, '8.5') for kk in ('AY', 'BY', 'NY', 'EY', 'EPY'))
                             + (tdv(k, 'DZ', '7.4') if mdisp == 'A' else tdv(k, 'DETR', '7.5')))
            tail = mid_s + tail_cols + f'{k:5d}'
            # tfdisp.f:572-575: 名前が 10 桁を超えるときは、後ろの欄 (buff(51:131)) の連続空白を lname-10 個、足りなければ前の欄 (buff(1:40)) からも詰めて桁を保つ
            tail_sq, nsq = _squeeze_space(tail, len(bname) - 10)
            lb = len(tail) - nsq
            head_sq, nsq2 = _squeeze_space(head, lb + len(bname) - 91)
            self.write(head_sq[:len(head) - nsq2] + bname + tail_sq[:lb] + chr(10))
            lines += 1

    def radint_cmd(self):
        """RADINT: INTGRL.f の直訳 (pysad/radint.py)。直前の CALC の Twiss (DP0 の行) を使う。"""
        from . import radint as _radint
        rows = self.ensure_rows()
        chrom = None
        if _compat.v2():
            # v2 (R2): LINEAR CHROMATICITY は運動量を ±1e-4 ずらした光学のチューンの差分 dν/dδ (六極・偏向磁石の寄与を含む)。
            # K64 の式 -I9/4π は磁石の中だけの ∮γ ds で、色収差にならない
            chrom = self._chrom_fd()
            if chrom is None:
                self.write(' RADINT: off-momentum optics not found; GX/GY below are the K64 integral (not a chromaticity)' + chr(10))
        for line in _radint.report(self, rows, chrom=chrom):
            self.write(line + chr(10))

    def _chrom_fd(self, h=1.0e-4):
        """一次の色収差 (dνx/dδ, dνy/dδ): DP0±h の 4 次元の光学 (CELL は周期解、そうでなければ先頭の条件からの伝搬) の位相の差分。
        求まらなければ None。状態 (rows、安定判定) は変えない"""
        d0 = self.dp0()
        q0 = getattr(self.lat, '_qstat', None)
        rf0 = _lat.RFSW
        _lat.RFSW = False
        try:
            nu = []
            for sgn in (1.0, -1.0):
                d = d0 + sgn * h
                if self.flags['CELL']:
                    rows = list(self.lat.twiss_table(delta=d, sad_stab={'intres': bool(self.flags['INTRES']),
                                                                        'halfres': bool(self.flags['HALFRES'])}))
                    q = self.lat._qstat
                    if q is not None and not (q['stabx'] and q['staby'] and not q['over']):
                        return None
                else:
                    rows = self._ins_rows(d)
                last = rows[len(self.els)] if len(rows) > len(self.els) else rows[-1]
                nu.append((last.mux / (2.0 * math.pi), last.muy / (2.0 * math.pi)))
            gx = (nu[0][0] - nu[1][0]) / (2.0 * h)
            gy = (nu[0][1] - nu[1][1]) / (2.0 * h)
            if not (math.isfinite(gx) and math.isfinite(gy)):
                return None
            return gx, gy
        except Exception:
            return None
        finally:
            _lat.RFSW = rf0
            self.lat._qstat = q0

    def chro_cmd(self):
        """CHROMATICITY: tfchro.f の直訳。四極・六極・多極 (MULT) の寄与を 1 要素ずつ表示し、最後に合計。
        現在の Twiss (直前の CALC、dp0 の行) の要素入口/出口の α β と分散 EX から閉形式で求める。BEND は寄与に数えない (SAD も同じ)。
        四極 / MULT の K1:  -(v β + (1+α²)/β L + α(i+1) - α(i)) / 2   (水平)、  -(-v β + (1+α²)/β L + Δα) / 2 (垂直)
        六極 / MULT の K2:   ±(β(i) η(i) + β(i+1) η(i+1)) v / 2        (水平 +、垂直 -)"""
        rows = self.ensure_rows()
        self.write(' Element           X Chromaticity   Y Chromaticity' + chr(10))
        gx = gy = 0.0
        for i in range(1, self.nlat):
            el = self.els[i - 1]
            cls = _cls(el)
            if cls not in ('Quadrupole', 'Sextupole', 'Multipole'):
                continue
            r0, r1 = rows[i - 1], rows[i]
            ali = float(el.L)
            xix = xiy = 0.0
            if cls in ('Quadrupole', 'Multipole'):
                v = self.get_value(el, 'K1')
                xix = -(v * r0.betx + (1.0 + r0.alfx ** 2) / r0.betx * ali + r1.alfx - r0.alfx) * 0.5
                xiy = -(-v * r0.bety + (1.0 + r0.alfy ** 2) / r0.bety * ali + r1.alfy - r0.alfy) * 0.5
            if cls in ('Sextupole', 'Multipole'):
                v = self.get_value(el, 'K2')
                if cls == 'Sextupole':
                    # Fortran は代入 (xiy=-(..)*v*.5) なので v=0 のとき -0.000 と表示される。足し込む MULT の場合は -0.0 にならない
                    xix = (r0.betx * r0.etx + r1.betx * r1.etx) * v * 0.5
                    xiy = -(r0.bety * r0.etx + r1.bety * r1.etx) * v * 0.5
                else:
                    xix = xix + (r0.betx * r0.etx + r1.betx * r1.etx) * v * 0.5
                    xiy = xiy - (r0.bety * r0.etx + r1.bety * r1.etx) * v * 0.5
            self.write(f' {self.elname(i):12s} {xix:17.3f}{xiy:17.3f}' + chr(10))
            gx += xix
            gy += xiy
        self.write(f' {"TOTAL":12s} {gx:17.3f}{gy:17.3f}' + chr(10))

    def save_cmd(self, pats):
        """SAVE [ALL] [pattern..]: 現在値を設計値に (誤差の基準、RESET の戻り先)"""
        fams = []
        for p in pats or ['*']:
            if p.upper() == 'ALL':
                continue
            fams.extend(self.families(p))
        # tfsave: LINE[key,"NAME.n"]= で触れた成分の値は、代表 (最初の 1 台) の値で全台に揃えられる (SAVE も EXPAND も同じ)
        self.spread_pending(set(fams or self.families('*')))
        rows = self.ensure_rows() if self.rows is not None else None
        for fam in fams or self.families('*'):
            el0 = self.instances(fam)[0]
            tn = type(el0).__name__
            raw = self.env._elements.setdefault(fam, {'TYPE': _TYPENAME.get(tn, 'DRIFT')})
            for k in _KEYS.get(tn, ['L']):
                if tn == 'Marker' and k in self._TWKEY and rows is not None:
                    pos = self.els.index(el0) + 1
                    if pos != 1:
                        raw[k] = self.twiss_value(k, pos)
                    continue
                v = self.get_value(el0, k)
                if v != 0.0 or k in raw:
                    raw[k] = v
            for el in self.instances(fam):
                el._saved = {k: self.get_value(el, k) for k in _KEYS.get(tn, ['L'])}
        if self.matcher:
            for v in self.matcher.variables:
                v['val'] = self.get_value(v['els'][0], v['key'])

    def reset_cmd(self):
        for el in self.els:
            sv = getattr(el, '_saved', None)
            if sv:
                for k, v in sv.items():
                    try:
                        self.set_value(el, k, v)
                    except Exception:
                        pass
            el._errk = 1.0
        self.invalidate()


# ---------------------------------------------------------------------------- 実行
_BLOCKED_SYMS = ('FitFunction', 'FitValue', 'ElementValues', 'OpticsProlog', 'OpticsEpilog')     # tffsa.f の tfblocksym / tfunblocksym の対象


def _locate_ffs(src: str):
    """MAIN 部の `FFS` 文を探す (gs で確かめた読み方)。無ければ None。
      - `FFS` は行の先頭か、同じ行の前の文の `;` の直後に書ける (`LINE RING = (...); FFS;`)。`FFS` と `ffs` だけ
        (gs と同じ。`Ffs`・`FFs` は gs の MAIN の構文エラーで、FFS 部に入れない)。
      - 文は次の `;` まで。`FFS USE RING;` / `FFS USE=RING;` は FFS 部の最初に USE RING をしたのと同じ。
      - FFS 部は、その `;` のある行の次の行から。`;` の後ろの同じ行の残り (`FFS;NPARA=1;` の NPARA=1;) は FFS 部では実行されない
        (gs は FFS を抜けた後に MAIN の文として読む)。
      - `;` の無い `FFS` だけの行は従来どおり、次の行から FFS 部 (gs は次の `;` までを FFS の引数として読む)。
    戻り値: {'main_end' (MAIN 部の終わり = FFS の語の位置), 'ffs_start' (FFS 部の始まり), 'use' (名前 | None), 'bad' (読めない引数), 'dropped' (捨てた残り)}"""
    from . import parser as _parser
    b = _parser._blank_comments(src)
    n = len(b)
    depth, q = 0, False
    i = 0
    at_start = True                       # 文の始まり (行の先頭の空白の後、または深さ 0 の ; の後)
    line_start = True                     # 行の先頭 (空白だけが前にある)
    while i < n:
        c = b[i]
        if q:
            if c == '"':
                q = False
            i += 1
            continue
        if c == '\n':
            line_start = True
            i += 1
            continue
        if c in ' \t\r':
            i += 1
            continue
        if ((depth == 0 and at_start) or line_start) and b[i:i + 3] in ('FFS', 'ffs') and (i + 3 >= n or b[i + 3] in ' \t\r\n;'):
            eol = b.find('\n', i)
            eol = n if eol < 0 else eol
            own = b[i + 3:eol]
            if (depth == 0 and at_start) or re.fullmatch(r'\s*;?\s*', own):      # 文の途中の行頭の FFS は、行が FFS だけのとき (従来の規則)
                if ';' not in own:
                    if own.strip():
                        k = b.find(';', eol)
                    else:
                        k = -1
                    if k < 0:
                        args, tail_end, dropped = own, eol, ''
                    else:
                        tail_end = b.find('\n', k)
                        tail_end = n if tail_end < 0 else tail_end
                        args, dropped = b[i + 3:k], b[k + 1:tail_end]
                else:
                    k = i + 3 + own.index(';')
                    args, tail_end, dropped = b[i + 3:k], eol, b[k + 1:eol]
                args = ' '.join(args.split())
                m = re.fullmatch(r'USE\s*=?\s*([A-Za-z_$][\w$.]*)', args, re.I)
                return {'main_end': i, 'ffs_start': min(n, tail_end + 1), 'use': m.group(1) if m else None,
                        'bad': '' if (m or not args) else args, 'dropped': dropped.strip()}
        at_start = False
        line_start = False
        if c == '"':
            q = True
        elif c in '([{':
            depth += 1
        elif c in ')]}':
            depth = max(0, depth - 1)
        elif c == ';' and depth == 0:
            at_start = True
        i += 1
    return None


def split_main_ffs(src: str):
    """FFS の文で MAIN 部と FFS 部に分ける (_locate_ffs)。MAIN 部の ON/OFF 行は [(flag, bool)] で返す"""
    loc = _locate_ffs(src)
    if loc is None:
        main, rest = src.split('\n'), []
    else:
        main = src[:loc['main_end']].split('\n')
        rest = src[loc['ffs_start']:].split('\n')
        if loc['use']:
            rest.insert(0, 'USE ' + loc['use'] + ';')
    flags = []
    main2 = []
    for ln in main:
        s = ln.split('!')[0].strip().rstrip(';').strip()
        m = re.match(r'^(ON|OFF)\s+(.+)$', s, re.I)
        if m:
            for f in m.group(2).split():
                flags.append((f.upper(), m.group(1).upper() == 'ON'))
            continue
        main2.append(ln)
    return '\n'.join(main2), '\n'.join(rest), flags


def run_sad(paths: List[str], out=None, cwd: Optional[str] = None):
    """複数ファイルを連結して (SAD の `gs lattice script` の慣習) 実行する"""
    from . import parser as _parser
    src = ''
    for p in paths:
        src += open(p, encoding='utf-8', errors='replace').read() + '\n'
    return run_sad_string(src, out=out, cwd=cwd or os.path.dirname(os.path.abspath(paths[-1])))


def main_part(src: str) -> str:
    """FFS の文より前 (MAIN 部) の原文 (ON/OFF も残す)"""
    loc = _locate_ffs(src)
    return src if loc is None else src[:loc['main_end']]


def run_sad_string(src: str, out=None, cwd: Optional[str] = None):
    from . import parser as _parser
    out = out or sys.stdout
    main, ffs_src, flags = split_main_ffs(src)
    if _locate_ffs(src) is None and re.search(r'(?m)(?:^|;)[ \t]*(?!FFS\b|ffs\b)[Ff][Ff][Ss]\b[ \t]*(?:;|$)', _parser._blank_comments(src)):
        # gs: `Ffs;`・`FFs;` (FFS と ffs 以外の綴り) は MAIN の構文エラーで、FFS 部には入れない (何も実行されない)
        out.write('LALR(1) parser: syntax error\n')
        return None, Interp(out=out)
    # MAIN 部の一般の代入 `NAME = 式;` (gs と同じく FFS に入ってからも見える。関数形・PI の値は要素定義でも使えるように数値に直す)
    assigns = _parser.parse_main_assignments(main)
    main = _parser.rewrite_main_assignments(main, assigns)
    sc = _parser.load_sad_string(main)
    I = Interp(out=out)
    if cwd:
        I.cwd = cwd
    consts = set(I.globals)           # E・Pi などの定数 (MAIN の `E = …` は gs でも FFS の E を変えない)
    ffs = FFS(sc.env, I, out=out)
    for a in assigns:
        if a['kind'] != 'VAR' or a['name'] in consts:
            continue
        if a['value'] is not None:    # MAIN の式として読めない右辺 (gs でも構文エラー) は代入しない
            I.assign(a['name'], float(a['value']))
    loc = _locate_ffs(src)
    if loc is not None and loc['bad']:
        out.write('  *** $FFS *** syntax error:undefined variable $%s (FFS の後に書けるのは ; か USE 名前; だけ)\n' % loc['bad'].split()[0])
    if loc is not None and loc['dropped']:
        out.write(' (pysad) FFS の文と同じ行の残り "%s" は実行しない (gs と同じ。FFS 部は次の行から)\n' % loc['dropped'][:60])
    from . import emit6d as _e6
    _e6.EMITDIV[:] = [1.0, 1.0, 1.0]  # 前の台本の EMITDIV を持ち越さない (Emittance[]・EMIT のたびに _sync_emitdiv が入れ直す)
    I.continue_on_error = True        # SAD と同じく、エラーの文は表示して次の文へ進む
    I.echo_out = True                 # SAD と同じく、; のない最上位の式の値を Out[n]:= で表示する
    # MAIN に TRACK 文があるときは、MAIN の命令 (ON/OFF, EMIT, `PHSPLOTS=..;`, TRACK, STOP) を書かれた順に実行する
    # (ON/OFF が EMIT/TRACK の前後で効き方が変わるため。TRACK の無いスクリプトは従来どおり ON/OFF を先に反映する)
    cmds = _parser.parse_main_commands(main_part(src))
    if any(c['cmd'] == 'TRACK' for c in cmds):
        try:
            for c in cmds:
                if c['cmd'] == 'STOP':
                    break
                ffs.main_exec(c)
            if not ffs.stopped:
                I.run(ffs_src, 'ffs')
        except SystemExit:
            pass
        return ffs, I
    for a in assigns:
        # MAIN の `EMITDIV = 20;`、`PBUNCH = 2E10;`、`EMITX = 1E-9;` など (_MAIN_SET_NAMES) を FFS に渡す (gs で確かめた: FFS に入っても値が残る)。
        # NPARA・DP・DP0・CHARGE は gs でも FFS の入口で既定値に戻るので渡さない
        U = a['name'].upper()
        if a['kind'] == 'SET' and a['value'] is not None and U not in ('NPARA', 'DP', 'DP0', 'CHARGE'):
            I.globals[U] = float(a['value'])
    for f, v in flags:
        if f in _FLAG_ALT:
            f, v = _FLAG_ALT[f][0], (_FLAG_ALT[f][1] if v else not _FLAG_ALT[f][1])
        if f == 'RFSW':                   # MAIN の $RFSW$ は FFS に届かない (track.f:149 rfsw=.true.; main_exec と同じ)
            ffs.main_flags['RFSW'] = v
            if _compat.v2():
                ffs.flags['RFSW'] = v         # v2 (B30): FFS に届ける (K64 は FFS 入口で rfsw=.true. に固定)
        elif f in ffs.flags:
            ffs.flags[f] = v
    try:
        I.run(ffs_src, 'ffs')
    except SystemExit:
        pass
    return ffs, I


import os  # noqa: E402  (run_sad で使用)
