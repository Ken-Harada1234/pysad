"""
SADScript の光学・追跡の補助関数 (Packages/Optics.n, emit.n, dapert.n の Python 書き直し、SurvivedParticles は tftrack.f):
  TransferMatrix[from, to, opts]  Optics.n TransferMatrixPrimitive (4D: Twiss と R1-R4 から、CALC6D: NormalCoordinate6 = etwiss2ri)
  BeamMatrix[i]                   Optics.n: SymplecticInverse[NormalCoordinates[i]] . Diag(EMITX,EMITX,EMITY,EMITY,0,DP²) . ᵀ
  OrbitGeo[k]                     Optics.n: 設計幾何 + GeoBase(chi)·(軌道のずれ)
  GeoBase[{chi1,chi2,chi3}]       Optics.n: 回転角 → 基底ベクトルの行列
  SymplecticJ[n]                  emit.n
  SurvivedParticles[x]            tftrack.f tfsurvivedparticles: flag≠0 の粒子だけ
  AccelerateParticles[beam, mom, {n, nturn}, Synchronize->f]   dapert.n: 周ごとに MOMENTUM を変えながら TrackParticles
方針 (docs/15): .n のものは Python で書き直す (直訳ではない)。
"""
import math
from typing import List

import numpy as np

from .sadscript import SADError, Sym, Rule, Expr, Null, S, is_num

MFIDX = {'AX': 0, 'BX': 1, 'NX': 2, 'AY': 3, 'BY': 4, 'NY': 5, 'EX': 6, 'EPX': 7, 'EY': 8, 'EPY': 9,
         'R1': 10, 'R2': 11, 'R3': 12, 'R4': 13, 'DETR': 14, 'DX': 15, 'DPX': 16, 'DY': 17, 'DPY': 18, 'DZ': 19, 'DDP': 20,
         'AZ': 21, 'BZ': 22, 'NZ': 23, 'ZX': 24, 'ZPX': 25, 'ZY': 26, 'ZPY': 27}
ROWATTR = ('alfx', 'betx', 'mux', 'alfy', 'bety', 'muy', 'etx', 'etpx', 'ety', 'etpy', 'r1', 'r2', 'r3', 'r4', 'detr',
           'x_cod', 'px_cod', 'y_cod', 'py_cod', 'z_cod', 'delta_cod', 'alfz', 'betz', 'muz', 'zx', 'zpx', 'zy', 'zpy')


def row28(r) -> np.ndarray:
    """TwissRow → Twiss["*", k] の 28 成分"""
    return np.array([float(getattr(r, a, 0.0) or 0.0) for a in ROWATTR])


def symplectic_j(n: int) -> np.ndarray:
    """emit.n SymplecticJ[n]: {{0,1,..},{-1,0,..},{0,0,0,1,..},{0,0,-1,0,..},..}"""
    m = np.zeros((n, n))
    for i in range(0, (n // 2) * 2 - 1, 2):
        m[i, i + 1] = 1.0
        m[i + 1, i] = -1.0
    return m


def symplectic_inverse(m: np.ndarray) -> np.ndarray:
    """emit.n SymplecticInverse[m] = -J.mᵀ.J"""
    n = m.shape[0]
    J = symplectic_j(n)
    return -J @ m.T @ J


def rmatrix4(p: np.ndarray) -> np.ndarray:
    """TransferMatrixPrimitive の rmatrix (4×4、結合 R1-R4、DETR<1 で分岐)"""
    r1, r2, r3, r4 = p[10], p[11], p[12], p[13]
    a1 = math.sqrt(max(0.0, 1.0 - (r1 * r4 - r2 * r3)))
    if p[14] < 1.0:
        return np.array([[a1, 0.0, -r4, r2], [0.0, a1, r3, -r1], [r1, r2, a1, 0.0], [r3, r4, 0.0, a1]])
    return np.array([[-r4, r2, a1, 0.0], [r3, -r1, 0.0, a1], [a1, 0.0, r1, r2], [0.0, a1, r3, r4]])


def normal_coordinates4(p: np.ndarray) -> np.ndarray:
    """Optics.n NormalCoordinatesPrimitive (4D): bmatrix.hmatrix.rmatrix (6×6、物理 → 正規モード座標)"""
    r1, r2, r3, r4 = p[10], p[11], p[12], p[13]
    a1 = math.sqrt(max(0.0, 1.0 - (r1 * r4 - r2 * r3)))
    if p[14] < 1.0:
        rm = np.array([[a1, 0, -r4, r2, 0, 0], [0, a1, r3, -r1, 0, 0], [r1, r2, a1, 0, 0, 0], [r3, r4, 0, a1, 0, 0],
                       [0, 0, 0, 0, 1, 0], [0, 0, 0, 0, 0, 1]], dtype=float)
    else:
        rm = np.array([[-r4, r2, a1, 0, 0, 0], [r3, -r1, 0, a1, 0, 0], [a1, 0, r1, r2, 0, 0], [0, a1, r3, r4, 0, 0],
                       [0, 0, 0, 0, 1, 0], [0, 0, 0, 0, 0, 1]], dtype=float)
    sbx, sby = math.sqrt(p[1]), math.sqrt(p[4])
    bm = np.array([[1 / sbx, 0, 0, 0, 0, 0], [p[0] / sbx, sbx, 0, 0, 0, 0], [0, 0, 1 / sby, 0, 0, 0],
                   [0, 0, p[3] / sby, sby, 0, 0], [0, 0, 0, 0, 1, 0], [0, 0, 0, 0, 0, 1]], dtype=float)
    ex, epx, ey, epy = p[6], p[7], p[8], p[9]
    hm = np.array([[1, 0, 0, 0, 0, -ex], [0, 1, 0, 0, 0, -epx], [0, 0, 1, 0, 0, -ey], [0, 0, 0, 1, 0, -epy],
                   [epx, -ex, epy, -ey, 1, 0], [0, 0, 0, 0, 0, 1]], dtype=float)
    return bm @ hm @ rm


def normal_coordinates6(p: np.ndarray) -> np.ndarray:
    """NormalCoordinate6 (tfinitn.f 1023 = etwiss2ri): 拡張 Twiss 28 成分 → 物理 → 正規モード座標の 6×6"""
    from .calc6d import etwiss2ri
    ri = etwiss2ri(p)
    return np.asarray(ri[0] if isinstance(ri, tuple) else ri, dtype=float)


def transfer_matrix(p1: np.ndarray, p2: np.ndarray, g1: float, g2: float, nx_tot: float = 0.0, ny_tot: float = 0.0,
                    fold: bool = False, symp: bool = False, c6d: bool = False) -> np.ndarray:
    """Optics.n TransferMatrixPrimitive: 位置 1 と 2 の Twiss (28) と gammab から転送行列 (4×4 / 6×6)"""
    gr = 1.0 if symp else math.sqrt(g1 / g2)
    dpsix = p2[2] - p1[2]
    dpsiy = p2[5] - p1[5]
    if fold:
        if dpsix <= 0.0:
            dpsix += nx_tot
        if dpsiy <= 0.0:
            dpsiy += ny_tot
    cx, sx, cy, sy = math.cos(dpsix), math.sin(dpsix), math.cos(dpsiy), math.sin(dpsiy)
    if c6d:
        dpsiz = p2[23] - p1[23]
        cz, sz = math.cos(dpsiz), math.sin(dpsiz)
        rm1, rm2 = normal_coordinates6(p1), normal_coordinates6(p2)
        m6 = np.array([cx * rm1[0] + sx * rm1[1], -sx * rm1[0] + cx * rm1[1],
                       cy * rm1[2] + sy * rm1[3], -sy * rm1[2] + cy * rm1[3],
                       cz * rm1[4] + sz * rm1[5], -sz * rm1[4] + cz * rm1[5]])
        m = symplectic_inverse(rm2) @ m6
        if (1.0 - p1[14]) * (1.0 - p2[14]) <= 0.0:
            m = m[[2, 3, 0, 1, 4, 5]]
        return gr * m
    bx1, ax1, by1, ay1 = p1[1], p1[0], p1[4], p1[3]
    bx2, ax2, by2, ay2 = p2[1], p2[0], p2[4], p2[3]
    m4 = np.zeros((4, 4))
    m4[0, 0] = math.sqrt(bx2 / bx1) * (cx + ax1 * sx)
    m4[0, 1] = math.sqrt(bx2 * bx1) * sx
    m4[1, 0] = -((ax2 - ax1) * cx + (1 + ax2 * ax1) * sx) / math.sqrt(bx2 * bx1)
    m4[1, 1] = math.sqrt(bx1 / bx2) * (cx - ax2 * sx)
    m4[2, 2] = math.sqrt(by2 / by1) * (cy + ay1 * sy)
    m4[2, 3] = math.sqrt(by2 * by1) * sy
    m4[3, 2] = -((ay2 - ay1) * cy + (1 + ay2 * ay1) * sy) / math.sqrt(by2 * by1)
    m4[3, 3] = math.sqrt(by1 / by2) * (cy - ay2 * sy)
    return gr * (symplectic_inverse(rmatrix4(p2)) @ m4 @ rmatrix4(p1))


def beam_matrix(p: np.ndarray, emx: float, emy: float, dp: float, d: float, gauss: bool, c6d: bool) -> np.ndarray:
    """Optics.n BeamMatrix[i]: r = SymplecticInverse[NormalCoordinates[i]]、r.Diag({emx,emx,emy,emy,0,d·(dp か dp/√3)²}·d).rᵀ。
    d = gammab(1)/gammab(i)"""
    nc = normal_coordinates6(p) if c6d else normal_coordinates4(p)
    r = symplectic_inverse(nc)
    s = dp if gauss else dp / math.sqrt(3.0)
    diag = np.array([emx, emx, emy, emy, 0.0, d * s * s]) * d
    return r @ np.diag(diag) @ r.T


def geo_base(chi) -> np.ndarray:
    """Optics.n GeoBase[{chi1,chi2,chi3}]: r3.r2.r1.g0 の基底ベクトル (行 = x_g, y_g, z_g)"""
    c1, c2, c3 = (math.cos(float(v)) for v in chi)
    s1, s2, s3 = (math.sin(float(v)) for v in chi)
    return np.array([[c3 * s1 - c1 * s2 * s3, -c1 * c3 - s1 * s2 * s3, c2 * s3],
                     [c1 * c3 * s2 + s1 * s3, c3 * s1 * s2 - c1 * s3, -c2 * c3],
                     [c1 * c2, c2 * s1, s2]])


def orbit_geo(gxyz, chi, cod6) -> np.ndarray:
    """Optics.n OrbitGeo: g + Σ GeoBase(chi)[i]·v[i]、v = {DX,DY,0} + {DPX,DPY,pz}/(1+dp)·DZ"""
    dx, dpx, dy, dpy, dz, dp = (float(v) for v in cod6)
    pz = math.sqrt(max(0.0, (1.0 + dp)**2 - dpx**2 - dpy**2))
    v = np.array([dx + dpx / (1.0 + dp) * dz, dy + dpy / (1.0 + dp) * dz, pz / (1.0 + dp) * dz])
    M = geo_base(chi)
    return np.asarray(gxyz, dtype=float) + v[0] * M[0] + v[1] * M[1] + v[2] * M[2]


def survived_particles(x: list) -> list:
    """tftrack.f tfsurvivedparticles: {x, px, y, py, z, dp, flag[, lost..]} (7 か 9 列) から flag≠0 の粒子だけ。全部落ちていれば空のリスト ×列数"""
    if not isinstance(x, list) or len(x) not in (7, 9) or not all(isinstance(c, list) for c in x):
        raise SADError('General::wrongtype: "{loc,{{x1,..},{px1,..},{y1,..},{py1,..},{z1,..},{dp1,..},{flg1,..}}}"')
    fl = x[6]
    n = len(fl)
    if any(len(c) != n for c in x):
        raise SADError('General::wrongtype: "{loc,{{x1,..},{px1,..},{y1,..},{py1,..},{z1,..},{dp1,..},{flg1,..}}}"')
    keep = [i for i in range(n) if float(fl[i]) != 0.0]
    if len(keep) == n:
        return x
    out = [[c[i] for i in keep] for c in x]
    out[6] = [1.0] * len(keep)
    return out


def _opt(opts, name, default):
    for o in opts:
        if isinstance(o, Rule) and isinstance(o.lhs, Sym) and o.lhs.name == name:
            return o.rhs
        if isinstance(o, list):
            v = _opt(o, name, None)
            if v is not None:
                return v
    return default


def _truthy(v):
    return bool(v) and v is not Null and not (isinstance(v, Sym) and v.name == 'False')


class OpticsFunctions:
    """FFS に組み込む関数群 (ffscmd.install から登録)"""

    def __init__(self, ffs):
        self.ffs = ffs

    def install(self, B):
        B['TransferMatrix'] = self.fn_transfer_matrix
        B['BeamMatrix'] = self.fn_beam_matrix
        B['OrbitGeo'] = self.fn_orbit_geo
        B['GeoBase'] = self.fn_geo_base
        B['SymplecticJ'] = lambda a: symplectic_j(int(a[0])).tolist()
        B['SymplecticInverse'] = lambda a: symplectic_inverse(np.array(a[0], dtype=float)).tolist()
        B['SurvivedParticles'] = lambda a: survived_particles(a[0] if a else None)
        B['AccelerateParticles'] = self.fn_accelerate_particles
        B['NormalCoordinates'] = self.fn_normal_coordinates

    # ---- 位置
    def _one_loc(self, a, what):
        ffs = self.ffs
        if isinstance(a, Sym):
            a = a.name
        r = ffs._resolve_loc(a)
        if len(r) != 1:
            raise SADError(f'General::wrongtype: Argument must be one location in {what}')
        return r[0]

    def _p28(self, k, fr=0.0):
        ffs = self.ffs
        r = ffs.row_frac(k, fr) if fr else ffs.row(k)
        return row28(r)

    def fn_transfer_matrix(self, args):
        """TransferMatrix[from, to, Fold->, Symplectic->, Calc6D->]"""
        ffs = self.ffs
        if len(args) < 2:
            raise SADError('General::narg: Number of arguments is expected to be "2 or more" in TransferMatrix')
        opts = list(args[2:])
        param = None
        if opts and isinstance(opts[0], list) and opts[0] and isinstance(opts[0][0], list):
            param = opts.pop(0)                  # TransferMatrix[from, to, param_List, opts]: param[[{from,to}]]、g={1,1}
        fold = _truthy(_opt(opts, 'Fold', False))
        symp = _truthy(_opt(opts, 'Symplectic', False))
        c6d = _truthy(_opt(opts, 'Calc6D', False)) or bool(ffs.flags.get('CALC6D'))
        rows = ffs.ensure_rows()
        if param is not None:
            i1, i2 = int(args[0]), int(args[1])
            p1 = np.array(param[i1 - 1], dtype=float)
            p2 = np.array(param[i2 - 1], dtype=float)
            g1 = g2 = 1.0
        else:
            (k1, f1), (k2, f2) = self._one_loc(args[0], 'TransferMatrix'), self._one_loc(args[1], 'TransferMatrix')
            p1, p2 = self._p28(k1, f1), self._p28(k2, f2)
            g1, g2 = ffs._gammab_at(k1, f1), ffs._gammab_at(k2, f2)
        end = rows[ffs.nlat - 1]
        m = transfer_matrix(p1, p2, g1, g2, float(end.mux), float(end.muy), fold, symp, c6d)
        return m.tolist()

    def fn_normal_coordinates(self, args):
        """NormalCoordinates[i]: {物理 → 正規モード座標の 6×6, 軌道 (DX..DDP)}"""
        ffs = self.ffs
        k, fr = self._one_loc(args[0], 'NormalCoordinates')
        p = self._p28(k, fr)
        nc = normal_coordinates6(p) if ffs.flags.get('CALC6D') else normal_coordinates4(p)
        return [nc.tolist(), p[15:21].tolist()]

    def fn_beam_matrix(self, args):
        ffs = self.ffs
        if not args:
            raise SADError('General::narg: Number of arguments is expected to be "1" in BeamMatrix')
        k, fr = self._one_loc(args[0], 'BeamMatrix')
        p = self._p28(k, fr)
        g = ffs.I.globals
        emx = float(g.get('EMITX', 0.0) or 0.0)
        emy = float(g.get('EMITY', 0.0) or 0.0)
        dp = float(g.get('DP', 0.0) or 0.0)
        d = ffs._gammab_at(1) / ffs._gammab_at(k, fr)
        return beam_matrix(p, emx, emy, dp, d, bool(ffs.flags.get('GAUSS')), bool(ffs.flags.get('CALC6D'))).tolist()

    def fn_geo_base(self, args):
        a = args[0] if args else None
        if isinstance(a, list) and len(a) == 3 and all(is_num(v) for v in a):
            return geo_base(a).tolist()
        if isinstance(a, list):
            return [self.fn_geo_base([x]) for x in a]
        raise SADError('General::wrongtype: Argument must be {chi1, chi2, chi3} in GeoBase')

    def fn_orbit_geo(self, args):
        """OrbitGeo[location]: 設計幾何 LINE["GEO"] と軌道 (Twiss DX..DDP) から現在の軌道の幾何 {GX, GY, GZ}"""
        ffs = self.ffs
        locs = ffs._resolve_loc(args[0].name if isinstance(args[0], Sym) else args[0]) if args else [(k, 0.0) for k in range(1, ffs.nlat + 1)]
        geo = ffs.ensure_geo()
        out = []
        for k, fr in locs:
            gv = geo.line_value('GEO', k, fr)
            r = ffs.row_frac(k, fr) if fr else ffs.row(k)
            out.append(orbit_geo(gv[0], gv[1], (r.x_cod, r.px_cod, r.y_cod, r.py_cod, r.z_cod, r.delta_cod)).tolist())
        if len(out) == 1 and args and not isinstance(args[0], list):
            return out[0]
        return out

    def fn_accelerate_particles(self, args):
        """dapert.n AccelerateParticles[beam, mom, {n, nturn}, Synchronize->f] (HoldAll: mom は周ごとに n を入れて評価)。
        周ごとに p0=MOMENTUM, MOMENTUM=p1=mom, u=p0/p1, v=u·sqrt((p1²+m²)/(p0²+m²)) として
        b = TrackParticles[{1,{x, px·u, y, py·u, z/v, dp·u−(1−u), flag}}, 1][[2]]。戻りは {1, b}"""
        ffs = self.ffs
        I = ffs.I
        if len(args) < 3:
            raise SADError('General::narg: Number of arguments is expected to be "3 or more" in AccelerateParticles')
        beam = I.eval(args[0])
        spec = args[2]
        if not (isinstance(spec, Expr) and spec.head == S('List') and len(spec.args) == 2) and not (isinstance(spec, list) and len(spec) == 2):
            raise SADError('General::wrongtype: Argument must be {n, nturn} for #3 in AccelerateParticles')
        nsym, nturn_e = (spec.args if isinstance(spec, Expr) else spec)
        if not isinstance(nsym, Sym):
            raise SADError('General::wrongtype: #3 must be {n_Symbol, nturn} in AccelerateParticles')
        nturn = int(I.eval(nturn_e))
        opts = [I.eval(o) for o in args[3:]]
        synch = _opt(opts, 'Synchronize', None)
        ffs.execute('RING RFSW')
        b = beam[1]
        m2 = float(I.get('MASS') or 0.51099895e6)**2
        g = I.globals
        for n in range(1, nturn + 1):
            g[nsym.name] = float(n)
            p0 = float(g.get('MOMENTUM'))
            p1 = float(I.eval(args[1]))
            g['MOMENTUM'] = p1
            u = p0 / p1
            v = u * math.sqrt((p1 * p1 + m2) / (p0 * p0 + m2))
            if synch is not None and synch is not Null:
                I.eval(Expr(synch, [float(n), [1.0, b]]))
            bb = [list(b[0]), [float(x) * u for x in b[1]], list(b[2]), [float(x) * u for x in b[3]],
                  [float(x) / v for x in b[4]], [float(x) * u - (1.0 - u) for x in b[5]], list(b[6])]
            res = ffs.fn_track([[1.0, bb], 1.0])
            b = res[1]
        g.pop(nsym.name, None)
        return [1.0, b]
