"""
geo.py — SAD の要素幾何 (src/tfgeo.f tfgeo/tfgeo1, src/tgrot.f module geolib) の直訳

SAD は各要素の「入口」について
    geo(3,4,nlat):  geo(:,1)=x 軸, geo(:,2)=y 軸, geo(:,3)=z 軸 (進行方向) の方向余弦 (実験室系の 3 成分)
                    geo(:,4)=入口の位置 (gx, gy, gz)
    pos(nlat):      入口の経路長 s
を持ち、FFS では LINE["GX"/"GY"/"GZ"/"GCHI1"/"GCHI2"/"GCHI3"/"GEO", loc] で参照する
(tftwiss.f tfline の 'GX'..'GCHI3' 分岐; Twiss["GX"] は K64 では
 "is must be to be name of optical function" になり使えない。CHI1..CHI3 は tfchi() で geo から求める)。

初期値 (tgrot.f geoini):  x 軸 = (0,-1,0), y 軸 = (0,0,-1), z 軸 = (1,0,0), 位置 = (0,0,0)。
すなわち実験室系の X が進行方向、ビーム座標の x は -Y 方向、y は -Z 方向 (右手系)。

角度の定義 (tgrot.f tfchi):
    CHI1 = atan2(z_y, z_x)        ... z 軸の水平方位角 (z_y = 0 なら 0)
    CHI2 = asinz(z_z)             ... z 軸の仰角
    CHI3 = atan2(x_z, -y_z)       ... ロール角 (x_z = 0 なら 0)
    ここで z_x = geo(1,3), z_y = geo(2,3), z_z = geo(3,3), x_z = geo(3,1), y_z = geo(3,2)。

要素ごとの進め方 (tfgeo1):
    ANGLE を持つ要素 (BEND, MULT) で ANGLE≠0:
        ROTATE=θ で (x,y) を回した x' = cosθ x − sinθ y, y' = sinθ x + cosθ y を使い
        ρ = L/φ, r1 = ρ sinφ, r2 = ρ (1−cosφ)   (xsincos で精度を確保)
        位置'  = 位置 + r1 z − r2 x'
        x''    =  cosφ x' + sinφ z,   z' = −sinφ x' + cosφ z
        θ≠0 なら x_new = cosθ x'' + sinθ y',  y_new = −sinθ x'' + cosθ y',  θ=0 なら y_new = y
    ANGLE=0 またはその他の要素: 位置' = 位置 + z·L、軸はそのまま
    COORD: CHI1..3, DX..DZ, DIR による回転・並進 (下記 _coord)
    MARK で OFFSET≠0: 後処理で位置 i+OFFSET の幾何をコピー (tfgeo)
    SOL (tsgeo.f: ソレノイド内の軌道追跡で幾何を決める) と MAP (tgmap) は未移植 → 直進扱い (注意を出す)

pos(i+1) = pos(i) + L ただし BEND/MULT で FRINGE≠0 かつ F1/FB1/FB2≠0 のときは
    L -= ((φ fb1)^2 + (φ fb2)^2)/L/48 · sin(½(φ(1−E1−E2) − AE1 − AE2)) / sin(½φ)
    (fb1 = F1+FB1, fb2 = F1+FB2 (BEND); MULT は FB1, FB2) の補正を入れる。

最後に ORG (tffsa.f の ORI_GIN/ORG: iorgr, geo0=geoini を tfrotgeo で回したもの) の変換:
    dgo = geo(:,4,iorg) + geo0_x·z_iorg − geo0_y·x_iorg − geo0_z·y_iorg  (geo0(:,4) = (gx0,gy0,gz0))
    geo(:,4,:) -= dgo ;  geo(:,1:4,:) = (geo0(:,1:3) · geo(:,1:3,iorg)^T) · geo(:,1:4,:)
既定 (iorg=1, geo0=geoini) では恒等変換。

要素の値の取り方: SAD のキーワード名 (L, ANGLE, ROTATE, FRINGE, F1, FB1, FB2, E1, E2, AE1, AE2, OFFSET, ...) で
el._raw (parser が残した元の値) → el の属性 の順に探す (_val)。E1/E2 は _raw にあれば SAD と同じ ANGLE 比、
属性しかなければ pysad の rad 値 (E1·ANGLE+AE1) なので AE1/AE2 を重ねて数えない。

接続 (ffscmd.py には手を入れていない): ffscmd.FFS.line_value の 'GX','GY','GZ','GCHI1','GCHI2','GCHI3','GEO' を
    g = geo.compute_geo(self.els);  g.line_value(key, k, fr)
に振ればよい (docs/history/work_geo.md 参照)。
"""
from __future__ import annotations

import math
import sys
from dataclasses import dataclass
from typing import List, Optional, Sequence, Tuple

import numpy as np

# ---------------------------------------------------------------- 定数 (tgrot.f geoini, tfloor.f mathfun)
# geoini: Fortran の reshape (列優先) なので列 = 軸ベクトル
GEOINI = np.array([[0.0, 0.0, 1.0, 0.0],
                   [-1.0, 0.0, 0.0, 0.0],
                   [0.0, -1.0, 0.0, 0.0]])          # geo[:, j]: j=0 x 軸, 1 y 軸, 2 z 軸, 3 位置

_EPSIEEE = 5.0e-17
_XTHS1 = (_EPSIEEE * 362880.0 / 6.0) ** (1.0 / 6.0)
_XTHS2 = (_EPSIEEE * 6227020800.0 / 6.0e9) ** (1.0 / 10.0)
_XTHC1 = (_EPSIEEE * 40320.0 / 2.0) ** (1.0 / 6.0)
_XTHC2 = (_EPSIEEE * 479001600.0 / 2.0) ** (1.0 / 10.0)
_XTHCS1 = min(_XTHS1, _XTHC1)
_XTHCS2 = min(_XTHS2, _XTHC2)

# SAD の要素型コード (MACCODE) のうち幾何で意味を持つもの
_IC_BEND, _IC_SOL, _IC_MULT, _IC_MAP, _IC_COORD, _IC_BEAM, _IC_MARK = 2, 20, 22, 33, 35, 36, 41
_SADTYPE_CODE = {'DRIFT': 1, 'BEND': 2, 'QUAD': 4, 'SEXT': 6, 'OCT': 8, 'DECA': 10, 'DODECA': 12,
                 'UND': 18, 'WIG': 19, 'SOL': 20, 'ST': 21, 'MULT': 22, 'TEST': 30, 'CAVI': 31, 'TCAV': 32,
                 'MAP': 33, 'INS': 34, 'COORD': 35, 'BEAM': 36, 'PROT': 37, 'MARK': 41, 'MONI': 42, 'APRT': 43}
_CLS2SADTYPE = {'Drift': 'DRIFT', 'Bend': 'BEND', 'Quadrupole': 'QUAD', 'Sextupole': 'SEXT', 'Octupole': 'OCT',
                'Decapole': 'DECA', 'Dodecapole': 'DODECA', 'Multipole': 'MULT', 'Solenoid': 'SOL',
                'Cavity': 'CAVI', 'TCavity': 'TCAV', 'Monitor': 'MONI', 'Marker': 'MARK'}


# ---------------------------------------------------------------- 数学関数 (tfloor.f mathfun)
def xsincos(x: float) -> Tuple[float, float, float, float]:
    """tfloor.f xsincos の直訳: (sin x, x−sin x, cos x, cos x − 1) を小さい x でも桁落ちなく返す"""
    ax = abs(x)
    if ax < _XTHCS1:
        x2 = x * x
        xs = x * x2 * (1.0 / 6.0 - x2 * (1.0 / 120.0 - x2 / 5040.0))
        s = x - xs
        dc = -x2 * (0.5 - x2 * (1.0 / 24.0 - x2 / 720.0))
        c = 1.0 + dc
    elif ax < _XTHCS2:
        x2 = x * x
        xs = x * x2 * (1.0 / 6.0 - x2 * (1.0 / 120.0 - x2 * (
            1.0 / 5040.0 - x2 * (1.0 / 362880.0 - x2 / 39916800.0))))
        s = x - xs
        dc = -x2 * (0.5 - x2 * (1.0 / 24.0 - x2 * (1.0 / 720.0 - x2 * (1.0 / 40320.0 - x2 / 3628800.0))))
        c = 1.0 + dc
    else:
        s = math.sin(x)
        xs = x - s
        c = math.cos(x)
        dc = -s * s / (1.0 + c) if c > 0.0 else c - 1.0
    return s, xs, c, dc


def asinz(x: float) -> float:
    """tfloor.f asinz の直訳: |x|<0.0032 は級数、それ以外は [-1,1] に切って asin"""
    if abs(x) < 0.0032:
        x2 = x * x
        return x + x * x2 * (1.0 / 6.0 + x2 * 3.0 / 40.0)
    return math.asin(min(1.0, max(-1.0, x)))


def _atan(y: float, x: float) -> float:
    """Fortran の atan(y,x) (= atan2)"""
    return math.atan2(y, x)


# ---------------------------------------------------------------- geolib (tgrot.f) の直訳
def tfchi(geo: np.ndarray, i: int) -> float:
    """tgrot.f tfchi: geo (3×3 以上) から CHI_i (i=1,2,3) を求める"""
    if i == 2:
        return asinz(float(geo[2, 2]))
    if i == 1:
        if geo[1, 2] == 0.0:
            return 0.0
        return _atan(float(geo[1, 2]), float(geo[0, 2]))
    if geo[2, 0] == 0.0:
        return 0.0
    return _atan(float(geo[2, 0]), -float(geo[2, 1]))


def tfchitogeo(chi: Sequence[float]) -> np.ndarray:
    """tgrot.f tfchitogeo: CHI1..3 → 3×3 の軸ベクトル (geoini を回したもの)"""
    c1, c2, c3 = math.cos(chi[0]), math.cos(chi[1]), math.cos(chi[2])
    s1, s2, s3 = math.sin(chi[0]), math.sin(chi[1]), math.sin(chi[2])
    g = np.empty((3, 3))
    g[0, 0] = s1 * c3 - c1 * s2 * s3
    g[1, 0] = -c1 * c3 - s1 * s2 * s3
    g[2, 0] = c2 * s3
    g[0, 1] = s1 * s3 + c1 * s2 * c3
    g[1, 1] = -c1 * s3 + s1 * s2 * c3
    g[2, 1] = -c2 * c3
    g[0, 2] = c1 * c2
    g[1, 2] = s1 * c2
    g[2, 2] = s2
    return g


def tfchitorot(chi: Sequence[float]) -> np.ndarray:
    """tgrot.f tfchitorot: CHI1..3 → 回転行列"""
    c1, c2, c3 = math.cos(chi[0]), math.cos(chi[1]), math.cos(chi[2])
    s1, s2, s3 = math.sin(chi[0]), math.sin(chi[1]), math.sin(chi[2])
    r = np.empty((3, 3))
    r[0, 0] = c1 * c2
    r[1, 0] = c2 * s1
    r[2, 0] = s2
    r[0, 1] = -c3 * s1 + c1 * s2 * s3
    r[1, 1] = c1 * c3 + s1 * s2 * s3
    r[2, 1] = -c2 * s3
    r[0, 2] = -c1 * c3 * s2 - s1 * s3
    r[1, 2] = -c3 * s1 * s2
    r[2, 2] = c2 * c3
    return r


def tfrotgeo(geo: np.ndarray, chi: Sequence[float]) -> np.ndarray:
    """tgrot.f tfrotgeo: 軸ベクトル geo(3×3) を chi1 (y まわり) → chi2 → chi3 の順に回す"""
    c1, s1 = math.cos(chi[0]), math.sin(chi[0])
    c2, s2 = math.cos(chi[1]), math.sin(chi[1])
    c3, s3 = math.cos(chi[2]), math.sin(chi[2])
    g1 = geo[:, 0] * c1 + geo[:, 2] * s1
    out = np.empty((3, 3))
    out[:, 2] = -geo[:, 0] * s1 + geo[:, 2] * c1
    out[:, 1] = geo[:, 1] * c2 + out[:, 2] * s2
    out[:, 2] = -geo[:, 1] * s2 + out[:, 2] * c2
    out[:, 0] = g1 * c3 - out[:, 1] * s3
    out[:, 1] = g1 * s3 + out[:, 1] * c3
    return out


def tfderotgeo(geo: np.ndarray, chi: Sequence[float]) -> np.ndarray:
    """tgrot.f tfderotgeo: tfrotgeo の逆順 (chi3 → chi2 → chi1 を戻す)"""
    c1, s1 = math.cos(chi[0]), math.sin(chi[0])
    c2, s2 = math.cos(chi[1]), math.sin(chi[1])
    c3, s3 = math.cos(chi[2]), math.sin(chi[2])
    g1 = geo[:, 0] * c3 + geo[:, 1] * s3
    out = np.empty((3, 3))
    out[:, 1] = -geo[:, 0] * s3 + geo[:, 1] * c3
    g3 = out[:, 1] * s2 + geo[:, 2] * c2
    out[:, 1] = out[:, 1] * c2 - geo[:, 2] * s2
    out[:, 2] = g3 * c1 + g1 * s1
    out[:, 0] = -g3 * s1 + g1 * c1
    return out


def tgrot(geo0: np.ndarray, geo1: np.ndarray) -> np.ndarray:
    """tgrot.f tgrot: geo0 の軸系から見た geo1 の向きを (chi1, chi2, chi3) で返す"""
    a = np.zeros(3)
    a[1] = asinz(-float(geo0[0, 1] * geo1[0, 2] + geo0[1, 1] * geo1[1, 2] + geo0[2, 1] * geo1[2, 2]))
    a13 = float(geo0[0, 0] * geo1[0, 2] + geo0[1, 0] * geo1[1, 2] + geo0[2, 0] * geo1[2, 2])
    if a13 != 0.0:
        a[0] = -_atan(a13, float(geo0[0, 2] * geo1[0, 2] + geo0[1, 2] * geo1[1, 2] + geo0[2, 2] * geo1[2, 2]))
    a21 = float(geo0[0, 1] * geo1[0, 0] + geo0[1, 1] * geo1[1, 0] + geo0[2, 1] * geo1[2, 0])
    if a21 != 0.0:
        a[2] = _atan(-a21, float(geo0[0, 1] * geo1[0, 1] + geo0[1, 1] * geo1[1, 1] + geo0[2, 1] * geo1[2, 1]))
    return a


def tforbitgeo(geo: np.ndarray, cod: Sequence[float]) -> np.ndarray:
    """tgrot.f tforbitgeo: 設計幾何 geo(3×4) と閉軌道 (x,px,y,py) から軌道の幾何 (OGX.. / OCHI..) を作る"""
    chi2i = -asinz(float(cod[3]))
    c2, s2 = math.cos(chi2i), math.sin(chi2i)
    chi1i = -asinz(float(cod[1]) / c2)
    c1, s1 = math.cos(chi1i), math.sin(chi1i)
    og = np.empty((3, 4))
    og[:, 3] = geo[:, 3] + geo[:, 0] * cod[0] + geo[:, 1] * cod[2]
    og[:, 0] = c1 * geo[:, 0] + s1 * geo[:, 2]
    g1 = -s1 * geo[:, 0] + c1 * geo[:, 2]
    og[:, 2] = c2 * g1 - s2 * geo[:, 1]
    og[:, 1] = s2 * g1 + c2 * geo[:, 1]
    return og


# ---------------------------------------------------------------- 要素の値 (SAD キーワード名で取り出す)
def _sad_type(el) -> str:
    t = getattr(el, '_sad_type', None)
    if t:
        return str(t).upper()
    return _CLS2SADTYPE.get(type(el).__name__, 'DRIFT')


def _val(el, key: str, default: float = 0.0) -> float:
    """SAD のキーワード名で要素の値を取る: _raw (parser が残した元の値) → 属性 の順"""
    raw = getattr(el, '_raw', None)
    if isinstance(raw, dict) and key in raw:
        v = raw[key]
        if isinstance(v, (int, float)):
            return float(v)
    v = getattr(el, key, None)
    if isinstance(v, (int, float)):
        return float(v)
    return default


def _angle_of(el, sad_type: str) -> float:
    """kytbl(kwANGL,k)≠0 の型 (BEND, MULT) の ANGLE。それ以外は None 扱い (0)"""
    if sad_type in ('BEND', 'MULT'):
        return _val(el, 'ANGLE')
    return 0.0


def _has_angle_key(sad_type: str) -> bool:
    return sad_type in ('BEND', 'MULT')


def _edge_sum(el) -> Tuple[float, float, float, float]:
    """(E1, E2, AE1, AE2) を SAD の意味 (E1/E2 は ANGLE 比) で返す。
    _raw に E1/E2 がなく属性しかないときは pysad の rad 値 (= E1·ANGLE + AE1) なので ANGLE 比に戻し AE は 0 にする。"""
    raw = getattr(el, '_raw', None) or {}
    if 'E1' in raw or 'E2' in raw or 'AE1' in raw or 'AE2' in raw or not hasattr(el, 'E1'):
        return (_val(el, 'E1'), _val(el, 'E2'), _val(el, 'AE1'), _val(el, 'AE2'))
    phi = _val(el, 'ANGLE')
    if phi == 0.0:
        return (0.0, 0.0, float(getattr(el, 'E1', 0.0) or 0.0), float(getattr(el, 'E2', 0.0) or 0.0))
    return (float(getattr(el, 'E1', 0.0) or 0.0) / phi, float(getattr(el, 'E2', 0.0) or 0.0) / phi, 0.0, 0.0)


def _geo_length(el, sad_type: str) -> float:
    """tfgeo1 の ali: L に BEND/MULT のフリンジ補正 (FRINGE≠0 かつ F1/FB1/FB2≠0) を入れたもの"""
    ali = _val(el, 'L')
    if not _has_angle_key(sad_type):
        return ali
    phi = _val(el, 'ANGLE')
    if _val(el, 'FRINGE') != 0.0 and phi * ali != 0.0:
        if sad_type == 'BEND':
            fb1 = _val(el, 'F1') + _val(el, 'FB1')
            fb2 = _val(el, 'F1') + _val(el, 'FB2')
        else:
            fb1 = _val(el, 'FB1')
            fb2 = _val(el, 'FB2')
        if fb1 != 0.0 or fb2 != 0.0:
            e1, e2, ae1, ae2 = _edge_sum(el)
            ali = ali - ((phi * fb1) ** 2 + (phi * fb2) ** 2) / ali / 48.0 \
                * math.sin(0.5 * (phi * (1.0 - e1 - e2) - ae1 - ae2)) / math.sin(0.5 * phi)
    return ali


# ---------------------------------------------------------------- 1 要素の幾何の進め方 (tfgeo1 の本体)
def _advance(geo_i: np.ndarray, el, sad_type: str, ali: float) -> np.ndarray:
    """要素 el の入口の幾何 geo_i (3×4) から出口の幾何 (3×4) を作る (tfgeo1 の calgeo 部分)"""
    g1 = np.array(geo_i, dtype=float)        # k >= icBEAM 等: まずコピー (既定の分岐でも同じ)
    if _has_angle_key(sad_type):
        v = _val(el, 'ANGLE')
        ald = ali
        if v == 0.0:
            g1[:, 3] = geo_i[:, 2] * ald + geo_i[:, 3]
            g1[:, 0:3] = geo_i[:, 0:3]
            return g1
        theta = _val(el, 'ROTATE')
        if theta != 0.0:
            cost = math.cos(theta)
            sint = math.sin(theta)
            x = cost * geo_i[:, 0] - sint * geo_i[:, 1]
            y = sint * geo_i[:, 0] + cost * geo_i[:, 1]
        else:
            x = geo_i[:, 0].copy()
            y = np.zeros(3)
        rho0 = ald / v
        sp0, xsp0, cp0, dcp0 = xsincos(v)
        r1 = rho0 * sp0
        r2 = -rho0 * dcp0
        g1[:, 3] = geo_i[:, 3] + (r1 * geo_i[:, 2] - r2 * x)
        g1[:, 0] = cp0 * x + sp0 * geo_i[:, 2]
        g1[:, 2] = -sp0 * x + cp0 * geo_i[:, 2]
        if theta != 0.0:
            g1[:, 1] = -sint * g1[:, 0] + cost * y
            g1[:, 0] = cost * g1[:, 0] + sint * y
        else:
            g1[:, 1] = geo_i[:, 1]
        return g1
    if sad_type == 'COORD':
        return _coord(geo_i, el)
    if sad_type == 'MAP':
        # tgmap (MAP 要素) は未移植: 直進扱い
        pass
    # 既定: 直進
    g1[:, 3] = geo_i[:, 2] * ali + geo_i[:, 3]
    g1[:, 0:3] = geo_i[:, 0:3]
    return g1


def _coord(geo_i: np.ndarray, el) -> np.ndarray:
    """tfgeo1 の icCOORD 分岐: CHI1..3 の回転と DX/DY/DZ の並進 (DIR=0 なら入口系で並進してから回転、
    DIR≠0 なら回転してから出口系で並進 (DX は符号反転))"""
    c1, s1 = math.cos(_val(el, 'CHI1')), math.sin(_val(el, 'CHI1'))
    c2, s2 = math.cos(_val(el, 'CHI2')), math.sin(_val(el, 'CHI2'))
    c3, s3 = math.cos(_val(el, 'CHI3')), math.sin(_val(el, 'CHI3'))
    direct = _val(el, 'DIR') == 0.0
    g1 = np.empty((3, 4))
    if direct:
        dx, dy, dz = _val(el, 'DX'), _val(el, 'DY'), _val(el, 'DZ')
        r11 = c1 * c3 + s1 * s2 * s3
        r12 = -c2 * s3
        r13 = s1 * c3 - c1 * s2 * s3
        r21 = -s1 * s2 * c3 + c1 * s3
        r22 = c2 * c3
        r23 = c1 * s2 * c3 + s1 * s3
        r31 = -s1 * c2
        r32 = -s2
        r33 = c1 * c2
        g1[:, 3] = geo_i[:, 3] + dx * geo_i[:, 0] + dy * geo_i[:, 1] + dz * geo_i[:, 2]
    else:
        dx, dy, dz = -_val(el, 'DX'), _val(el, 'DY'), _val(el, 'DZ')
        r11 = c1 * c3 + s1 * s2 * s3
        r21 = c2 * s3
        r31 = -s1 * c3 + c1 * s2 * s3
        r12 = s1 * s2 * c3 - c1 * s3
        r22 = c2 * c3
        r32 = c1 * s2 * c3 + s1 * s3
        r13 = s1 * c2
        r23 = -s2
        r33 = c1 * c2
    g1[:, 0] = r11 * geo_i[:, 0] + r12 * geo_i[:, 1] + r13 * geo_i[:, 2]
    g1[:, 1] = r21 * geo_i[:, 0] + r22 * geo_i[:, 1] + r23 * geo_i[:, 2]
    g1[:, 2] = r31 * geo_i[:, 0] + r32 * geo_i[:, 1] + r33 * geo_i[:, 2]
    if not direct:
        g1[:, 3] = geo_i[:, 3] + dx * g1[:, 0] + dy * g1[:, 1] + dz * g1[:, 2]
    return g1


# ---------------------------------------------------------------- 部分要素 (tfgeofrac / qfracseg)
def _frac_values(el, sad_type: str, fr: float):
    """qfracseg (qtwiss.f) に合わせた部分要素の値: L·fr, ANGLE·fr, E2=0, フリンジキーの置き換え。
    幾何に効くのは L, ANGLE, ROTATE, FRINGE, F1/FB1/FB2, E1/E2/AE1/AE2 だけなので、それらを持つ軽い代理オブジェクトを返す"""

    class _Part:
        pass

    p = _Part()
    raw = {}
    raw['L'] = _val(el, 'L') * fr
    if _has_angle_key(sad_type):
        # qfraccomp (rx1=0, rx2=fr → f1=1, f2=0, r=fr): ANGLE·r, E1 → E1·f1/r = E1/r (ANGLE 比なので絶対角は不変),
        # E2 → E2·f2/r = 0, AE1/AE2 はそのまま。
        # BEND: FRMD=0 なら F1=0、FRMD≠0 なら FRMD=-f1-2f2=-1 (入口のみ)。FB1/FB2 はそのまま。
        # MULT (ANGLE≠0): FRMD=0 なら FB1=FB2=0、FRMD≠0 なら FRMD=-1。
        # 逆向き要素 (cmp%ori=.false.) では E1·f2/r=0, E2·f1/r=E2/r になるが pysad の要素は向きを持たないので未対応。
        raw['ANGLE'] = _val(el, 'ANGLE') * fr
        raw['ROTATE'] = _val(el, 'ROTATE')
        e1, e2, ae1, ae2 = _edge_sum(el)
        raw['E1'], raw['E2'], raw['AE1'], raw['AE2'] = e1 / fr, 0.0, ae1, ae2
        frmd = _val(el, 'FRINGE')
        if sad_type == 'BEND':
            raw['F1'] = 0.0 if frmd == 0.0 else _val(el, 'F1')
            raw['FRINGE'] = 0.0 if frmd == 0.0 else -1.0
            raw['FB1'] = _val(el, 'FB1')
            raw['FB2'] = _val(el, 'FB2')
        else:
            raw['FB1'] = 0.0 if frmd == 0.0 else _val(el, 'FB1')
            raw['FB2'] = 0.0 if frmd == 0.0 else _val(el, 'FB2')
            raw['FRINGE'] = 0.0 if frmd == 0.0 else -1.0
    p._raw = raw
    p._sad_type = sad_type
    return p


# ---------------------------------------------------------------- 結果
@dataclass
class GeoResult:
    """compute_geo の結果。geo[k-1] が SAD の geo(:,:,k) (k=1..nlat; nlat は要素数+1 で末尾 $$$)"""
    geo: np.ndarray          # (nlat, 3, 4)
    pos: np.ndarray          # (nlat,)  入口の s
    elements: Sequence       # 要素列 (部分要素の計算に使う)
    warnings: List[str]

    @property
    def nlat(self) -> int:
        return int(self.geo.shape[0])

    def gxyz(self, k: int) -> np.ndarray:
        """位置 k (1 始まり) の入口位置 (gx, gy, gz)"""
        return self.geo[k - 1, :, 3]

    def axes(self, k: int) -> np.ndarray:
        """位置 k の方向余弦 3×3 (列 0: x 軸, 1: y 軸, 2: z 軸)"""
        return self.geo[k - 1, :, 0:3]

    def chi(self, k: int) -> Tuple[float, float, float]:
        """位置 k の (CHI1, CHI2, CHI3) (tfchi)"""
        g = self.geo[k - 1]
        return tfchi(g, 1), tfchi(g, 2), tfchi(g, 3)

    def geo_frac(self, k: int, fr: float) -> np.ndarray:
        """tfgeofrac: 要素 k の入口から割合 fr の位置の幾何 (3×4)。fr=0 なら入口そのもの"""
        if fr == 0.0 or k >= self.nlat:
            return self.geo[min(k, self.nlat) - 1]
        el = self.elements[k - 1]
        st = _sad_type(el)
        if not (_val(el, 'L') != 0.0 or (_has_angle_key(st) and _val(el, 'ANGLE') != 0.0)):
            return self.geo[k - 1]              # qfracseg: 変化なし (chg=.false.)
        part = _frac_values(el, st, fr)
        return _advance(self.geo[k - 1], part, st, _geo_length(part, st))

    def pos_frac(self, k: int, fr: float) -> float:
        if fr == 0.0 or k >= self.nlat:
            return float(self.pos[min(k, self.nlat) - 1])
        el = self.elements[k - 1]
        st = _sad_type(el)
        part = _frac_values(el, st, fr)
        return float(self.pos[k - 1] + _geo_length(part, st))

    def line_value(self, key: str, k: int, fr: float = 0.0):
        """LINE["GX"/"GY"/"GZ"/"GCHI1"/"GCHI2"/"GCHI3"/"GEO", k+fr] (tftwiss.f tfline) と同じ値"""
        key = key.upper()
        gv = self.geo_frac(k, fr)
        if key == 'GX':
            return float(gv[0, 3])
        if key == 'GY':
            return float(gv[1, 3])
        if key == 'GZ':
            return float(gv[2, 3])
        if key == 'GCHI1':
            return tfchi(gv, 1)
        if key == 'GCHI2':
            return tfchi(gv, 2)
        if key == 'GCHI3':
            return tfchi(gv, 3)
        if key == 'GEO':
            return [[float(gv[0, 3]), float(gv[1, 3]), float(gv[2, 3])],
                    [tfchi(gv, 1), tfchi(gv, 2), tfchi(gv, 3)]]
        raise KeyError(key)

    def orbit_value(self, key: str, k: int, cod: Sequence[float], fr: float = 0.0):
        """LINE["OGX".."OCHI3"/"OGEO", k+fr]: 閉軌道 cod=(x,px,y,py) を乗せた幾何 (tforbitgeo)"""
        key = key.upper()
        og = tforbitgeo(self.geo_frac(k, fr), cod)
        m = {'OGX': ('G', 0), 'OGY': ('G', 1), 'OGZ': ('G', 2), 'OCHI1': ('C', 1), 'OCHI2': ('C', 2), 'OCHI3': ('C', 3)}
        if key == 'OGEO':
            return [[float(og[0, 3]), float(og[1, 3]), float(og[2, 3])], [tfchi(og, 1), tfchi(og, 2), tfchi(og, 3)]]
        kind, i = m[key]
        return float(og[i, 3]) if kind == 'G' else tfchi(og, i)


# ---------------------------------------------------------------- 本体 (tfgeo + tfgeo1)
def compute_geo(elements: Sequence, gx0: float = 0.0, gy0: float = 0.0, gz0: float = 0.0,
                chi0: Sequence[float] = (0.0, 0.0, 0.0), iorg: int = 1,
                geo_init: Optional[np.ndarray] = None, warn=None) -> GeoResult:
    """SAD tfgeo (GEOCAL 時, calgeo=.true.) の直訳。

    elements: pysad の要素列 (lat.elements / FFS.els)。要素数 n に対し nlat=n+1 で末尾が $$$。
    gx0, gy0, gz0, chi0, iorg: FFS の ORG コマンド (ORG 要素名 gx gy gz chi1 chi2 chi3; chi は度で与えるので
        ここではラジアンに直して渡す) に相当。既定 (iorg=1, 0, 0) では tfgeo の最終変換は恒等。
    geo_init: 先頭 (位置 1) の幾何 3×4。既定 geoini。
    warn: 警告の出力先 (callable)。None なら sys.stderr。
    戻り: GeoResult (geo[k-1] = SAD geo(:,:,k), pos[k-1] = pos(k))。"""
    warnings: List[str] = []
    n = len(elements)
    nlat = n + 1
    geo = np.zeros((nlat, 3, 4))
    pos = np.zeros(nlat)
    geo[0] = GEOINI if geo_init is None else np.asarray(geo_init, dtype=float)
    # ---- tfgeo1(1, nlat): 要素ごとに進める
    for i in range(n):
        el = elements[i]
        st = _sad_type(el)
        if st == 'SOL':
            # tsgeo (ソレノイド境界内の軌道で幾何を決める) は未移植。BOUND/GEO 付きの SOL は正しく出ない
            if _val(el, 'BOUND') != 0.0 or _val(el, 'GEO') != 0.0:
                warnings.append(f'SOL {getattr(el, "name", "?")} (位置 {i + 1}): tsgeo 未移植、直進扱い')
        ali = _geo_length(el, st)
        pos[i + 1] = pos[i] + ali
        geo[i + 1] = _advance(geo[i], el, st, ali)
    # ---- tfgeo: MARK の OFFSET (tffsmarkoffset)
    for i in range(n):
        el = elements[i]
        if _sad_type(el) != 'MARK':
            continue
        offset = _val(el, 'OFFSET')
        if offset == 0.0:
            continue
        xp = max(1.0, min(float(nlat), offset + float(i + 1)))
        lxp = int(xp)
        fr = xp - lxp
        if fr == 0.0:
            geo[i] = geo[lxp - 1]
            pos[i] = pos[lxp - 1]
        else:
            tmp = GeoResult(geo, pos, elements, warnings)
            gfrac = tmp.geo_frac(lxp, fr)
            pfrac = tmp.pos_frac(lxp, fr)
            geo[i] = gfrac
            pos[i] = pfrac
        if i == 0:
            pos = pos - pos[0]
        if i == n - 1:
            geo[nlat - 1] = geo[i]
            pos[nlat - 1] = pos[i]
    # ---- tfgeo: ORG による変換 (calgeo 時)
    geo0 = GEOINI.copy()
    geo0[:, 3] = (gx0, gy0, gz0)
    chi0 = np.asarray(chi0, dtype=float)
    if np.any(chi0 != 0.0):
        geo0[:, 0:3] = tfrotgeo(geo0[:, 0:3], -chi0)
    io = iorg - 1
    dgo = geo[io, :, 3] + geo0[0, 3] * geo[io, :, 2] - geo0[1, 3] * geo[io, :, 0] - geo0[2, 3] * geo[io, :, 1]
    geo[:, :, 3] -= dgo
    geo2 = geo0[:, 0:3] @ geo[io, :, 0:3].T
    geo = np.einsum('ab,kbc->kac', geo2, geo)
    if warn is None:
        warn = lambda s: print(s, file=sys.stderr)
    for w in warnings:
        warn(w)
    return GeoResult(geo, pos, elements, warnings)
