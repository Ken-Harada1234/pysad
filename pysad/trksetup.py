# -*- coding: utf-8 -*-
"""Ring-tracking setup shared by the frequency map and the local momentum acceptance.
周回追跡の設定 (v2.3.0 Q1: 周波数マップと局所運動量アクセプタンスの共通部分)。

FFS の状態 (要素の並び、運動量、RF、フラグ、Emittance[] の結果) から、周回追跡に要るものを 1 つの束 (TrackSetup) にまとめる。
追跡は track6n.TrackRun (numba) を使う。これは SAD の TrackParticles と同じ損失判定 (tapert: 最初の要素・15 要素ごと・区間の終わりで
|x|,|y| < LOSSAMPL、|z| ≤ LOSSDZ、APERT 要素) をする。位置 loc (1 始まり、要素 loc の入口) から nturn 周する区間の分け方は、既定
では TrackParticles[{loc, beam}, loc, 1, nturn] の 1 回の呼び出しと同じ (loc→nlat、(nturn−1) 回の 1→nlat、1→loc)。周ごとの座標を記録するとき (周波数マップ) は
1 周ずつ loc から回す (1 周ずつ TrackParticles[{loc, beam}, loc] を呼ぶのと同じ)。

主な内容:
  TrackSetup.from_ffs(ffs, ...)     FFS の状態から作る (必要なら Emittance[] を走らせる)
  TrackSetup.closed_orbit(...)      位置 loc の閉軌道 (放射・RF を含む 6 次元。RFSW OFF なら δ を決めた 4 次元)
  TrackSetup.track(loc, coords, nturn, record=...)   位置 loc から nturn 周。周ごとの座標と損失の周を返す
  find_periodicity(els)             要素の並びの周期 (同じ セル が何回並ぶか)
  select_locations(...)             要素の境界の位置の選び方 (範囲・間隔)
  lower_priority()                  この処理の優先度を下げる (並列の子プロセス用)

単位: 座標は SAD と同じ (x, px, y, py, z, dp)。位置 loc の s は Twiss の行 (入口の s)。
"""
import math
import os

import numpy as np

from . import compat
from . import track6 as _t6
from .track6n import TrackRun, state_out, MAXTURN_LOST
from .track6n import make_parallel_run as _make_parallel_run, max_workers as _max_workers
from .emit6d import _emit_circ

_SKIP = ('Marker', 'Monitor', 'Cavity', 'Coord', 'Map')


def lower_priority():
    """Lower the priority of this process (BELOW_NORMAL on Windows, nice elsewhere).
    この処理の優先度を下げる (Windows: BELOW_NORMAL、それ以外: nice)。並列の子プロセスと、重い計算の本体で使う"""
    try:
        if os.name == 'nt':
            import ctypes
            k = ctypes.windll.kernel32
            k.SetPriorityClass(k.GetCurrentProcess(), 0x00004000)       # BELOW_NORMAL_PRIORITY_CLASS
        else:
            os.nice(10)
    except Exception:
        pass


class TrackSetup:
    """Ring-tracking setup holding els, nlat, p0, circ, trf0, rf/rad/fluc flags, loss limits, cod1, bucket and Twiss rows.
    周回追跡の設定。属性: els、nlat (= 要素の数 + 1)、p0 (運動量 eV)、circ、trf0、rfsw、rad、fluc、fshift、mass_charge、
    alost・zlost (LOSSAMPL・LOSSDZ)、cod1 (位置 1 の閉軌道)、bucket (RF バケットの情報)、rows (位置ごとの Twiss: s・βx・αx・βy・αy・ηx・η'x・μx・μy)"""

    def __init__(self, els, momentum_eV, trf0, rfsw, rad, fluc, fshift, mass_charge, alost, zlost, cod1, bucket, rows,
                 names, mode=None, bend_k1_edge=None, info=None, bend_k1_def=None):
        self.els = list(els)
        self.nlat = len(self.els) + 1
        self.p0 = float(momentum_eV)
        self.circ = _emit_circ(self.els)
        self.trf0 = float(trf0)
        self.rfsw = bool(rfsw)
        self.rad = bool(rad)
        self.fluc = bool(fluc)
        self.fshift = float(fshift)
        self.mass_charge = dict(mass_charge or {})
        self.alost = float(alost)
        self.zlost = float(zlost)
        self.cod1 = np.array(cod1, dtype=float)
        self.bucket = dict(bucket or {})
        self.rows = rows                  # [(s, betx, alfx, bety, alfy, etx, etpx), ...] 位置 1..nlat
        self.names = list(names)
        self.mode = mode if mode is not None else compat.mode()
        self.bend_k1_edge = bend_k1_edge if bend_k1_edge is not None else compat.bend_k1_edge()
        self.bend_k1_def = bend_k1_def if bend_k1_def is not None else compat.bend_k1_def()     # v2.3.0 U2 (B160 のオプション)
        self.info = dict(info or {})
        self.seed = None
        self.check_every = 0              # 0: SAD の tapert の判定のまま (最初の要素・15 要素ごと・区間の終わり)。k>0: k 要素ごとに区間を分けて、k 要素ごとに判定する
        self._orbits = {}

    def __setstate__(self, state):
        # 並列の worker (freqmap・lma。spawn で別プロセス) に K1 の定義 (compat.set_bend_k1_def) を引き継ぐ。unpickle した側の compat を設定の値にそろえる
        self.__dict__.update(state)
        d = state.get('bend_k1_def')
        if d is not None and compat.mode() == 'v2':
            try:
                compat.set_bend_k1_def(d)
            except Exception:
                pass

    # ------------------------------------------------------------------ 作る
    @classmethod
    def from_ffs(cls, ffs, rad=None, rfsw=None, fluc=False, alost=None, zlost=None, refit_orbit=True):
        """Build a `TrackSetup` from an FFS object of ffscmd (after `USE` and `CALC`).
        ffscmd の FFS オブジェクト (USE・CALC 済み) から作る。
        rad・rfsw: None なら FFS のフラグ (RAD・RFSW)。指定すると、そのフラグで Emittance[] (閉軌道・同期位相の決定) を走らせ直す
        (フラグはあとで元に戻す)。Emittance[] が未実行なら走らせる。alost・zlost: None なら LOSSAMPL・LOSSDZ (SAD の既定は 1 m と 100 m)"""
        if ffs.lat is None:
            raise ValueError('USE していない (ビームラインがない)')
        if ffs.rows is None:
            ffs.calc(quiet=True)
        fl = ffs.flags
        save = (fl['RAD'], fl['RFSW'], fl['RADCOD'])
        want_rad = bool(fl['RAD']) if rad is None else bool(rad)
        want_rf = bool(fl['RFSW']) if rfsw is None else bool(rfsw)
        try:
            fl['RAD'], fl['RFSW'] = want_rad, want_rf
            if want_rad:
                fl['RADCOD'] = True            # 放射つきの追跡の閉軌道は、放射を含む閉軌道 (RADCOD) から出発する
            need = (want_rf or want_rad) and (ffs.emit is None or save[0] != want_rad or save[1] != want_rf or (want_rad and not save[2]))
            if need:
                ffs.run_emit(quiet=True)
            em = ffs.emit if (want_rf or want_rad) else None
            trf0 = float(ffs.trf0) if em is not None and em.get('trf0') == em.get('trf0') else 0.0
            if not want_rad:
                trf0 = 0.0      # v2.3.0 S2: RAD OFF は放射損失がないので同期位相 0 (Emittance[] が出す trf0 = 損失 U0 を補う位相 φs を使うと、1 周ごとに正味の加速が入る。gs の CALC のあとの TrackParticles (trf0 = 0) と同じ)
        finally:
            fl['RAD'], fl['RFSW'], fl['RADCOD'] = save
        g = ffs.I.globals
        if alost is None:
            alost = float(g.get('LOSSAMPL', 1.0) or 0.0)
        if zlost is None:
            zlost = float(g.get('LOSSDZ', 100.0) or 0.0)
        cod1 = np.array(em['closed_orbit'], dtype=float) if (em is not None and want_rf) else np.zeros(6)
        if not np.all(np.isfinite(cod1)):
            cod1 = np.zeros(6)
        bucket = {}
        if em is not None and want_rf:
            wrf = float(em.get('wrfeff', 0.0) or 0.0)
            bucket = dict(height=float(em.get('bucket_height', 0.0) or 0.0), k_rf=wrf, trf0=trf0,
                          phis=trf0 * wrf, nus=float(em['tunes'][2]), U0=float(em.get('energy_loss_eV', 0.0)),
                          Vrf=float(em.get('rf_voltage_eff', 0.0) or 0.0),
                          alpha=float(em.get('momentum_compaction', 0.0) or 0.0))
        rows = [(float(r.s), float(r.betx), float(r.alfx), float(r.bety), float(r.alfy), float(r.etx), float(r.etpx), float(r.mux), float(r.muy))
                for r in ffs.rows]
        names = [e.name for e in ffs.els] + ['$$$']
        st = cls(ffs.els, ffs.momentum(), trf0, want_rf, want_rad, fluc, ffs._fshift(), ffs._mass_charge(), alost, zlost, cod1,
                 bucket, rows, names,
                 info=dict(rad=want_rad, rfsw=want_rf, fluc=bool(fluc), momentum_eV=float(ffs.momentum()), circ=None))
        st.info['circ'] = st.circ
        if em is not None:
            st.info['tunes_emit'] = [float(v) for v in em['tunes']]
        if refit_orbit and want_rf:
            st.cod1 = st.refine_orbit(st.cod1)
        return st

    @classmethod
    def from_source(cls, src, ring='RING', rad=True, rfsw=True, fluc=False, alost=None, zlost=None, extra=''):
        """Build a `TrackSetup` from SAD lattice source text (USE ring; CELL; CALC; Emittance); returns (TrackSetup, FFS object).
        SAD のラティス定義の文字列から作る (USE ring; CELL; CALC; Emittance)。extra: USE の前に実行する SAD の文 (例: 'MOMENTUM=3GEV;')。
        戻り (TrackSetup, FFS オブジェクト)"""
        import io
        from . import ffscmd
        flags = ('RFSW; ' if rfsw else 'NORFSW; ') + ('RAD; RADCOD; ' if rad else 'NORAD; ') + 'NOFLUC;'
        ffs, _ = ffscmd.run_sad_string(src + '\n' + extra + '\nFFS;\nUSE %s; CELL; CALC; %s\n' % (ring, flags), out=io.StringIO())
        return cls.from_ffs(ffs, rad=rad, rfsw=rfsw, fluc=fluc, alost=alost, zlost=zlost), ffs

    def copy(self, **kw):
        """Return a copy with some settings changed (rad, rfsw, fluc, alost, zlost, check_every), re-solving the closed orbit if needed.
        一部を変えた複製 (rad・rfsw・fluc・alost・zlost・check_every)。rad・rfsw を変えたときの閉軌道は、cod1 を出発点に求め直す"""
        import copy as _copy
        o = _copy.copy(self)
        o._orbits = {}
        o.info = dict(self.info)
        for k, v in kw.items():
            if k not in ('rad', 'rfsw', 'fluc', 'alost', 'zlost', 'check_every'):
                raise TypeError('copy: ' + k)
            setattr(o, k, bool(v) if k in ('rad', 'rfsw', 'fluc') else int(v) if k == 'check_every' else float(v))
        o.info.update({k: kw[k] for k in ('rad', 'rfsw', 'fluc') if k in kw})
        if 'rfsw' in kw and not o.rfsw:
            o.cod1 = np.zeros(6)
        return o

    # ------------------------------------------------------------------ 追跡
    def make_run(self, coords, alost=None, zlost=None):
        c = np.atleast_2d(np.asarray(coords, dtype=float))
        if self.seed is not None:
            _t6.RNG.restore(self.seed)
        return self._run_of(c, alost, zlost)

    def _run_of(self, c, alost=None, zlost=None):
        c = np.atleast_2d(np.asarray(c, dtype=float))
        return TrackRun(self.els, c, [1.0] * len(c), self.p0, alost=self.alost if alost is None else alost,
                        zlost=self.zlost if zlost is None else zlost, kturn0=1, trf0=self.trf0, rfsw=self.rfsw,
                        rad=self.rad, fluc=self.fluc, fshift=self.fshift, circ=self.circ, **self.mass_charge)

    def _seg(self, tr, a, b):
        """tturn(a, b)。check_every = k > 0 なら k 要素ずつの区間に分ける (区間ごとに tapert が最初の要素と終わりで判定する)"""
        k = self.check_every
        if k and b - a > k:
            for x in range(a, b, k):
                tr.segment(x, min(x + k, b))
        else:
            tr.segment(a, b)

    def one_turn(self, tr, loc):
        """Track one turn from location loc (entrance of element loc) as `tturn(loc, nlat)` then `tturn(1, loc)`, like `TrackParticles`.
        位置 loc (要素 loc の入口) から 1 周: tturn(loc, nlat)、tturn(1, loc) (TrackParticles と同じ)"""
        with np.errstate(all='ignore'):
            self._seg(tr, loc, self.nlat)
            if loc > 1:
                self._seg(tr, 1, loc)
        tr.next_turn()

    def run_call(self, tr, loc, nturn):
        """Track like one `TrackParticles[{loc, beam}, loc, 1, nturn]` call, with the same splitting into `tturn` segments.
        TrackParticles[{loc, beam}, loc, 1, nturn] (1 回の呼び出し) と同じ区間の分け方: tturn(loc, nlat)、(nturn−1) 回の tturn(1, nlat)、tturn(1, loc)。
        tapert (最初の要素・15 要素ごと・区間の終わり) の判定の位置は区間の始めから数えるので、1 周ずつ loc から回す one_turn とは、判定の位置が違う
        (損失ぎりぎりの粒子の判定が違いうる)。gs の 1 回の TrackParticles と同じ判定を得るにはこちら"""
        nlat = self.nlat
        with np.errstate(all='ignore'):
            self._seg(tr, loc, nlat)
            tr.next_turn()
            for _ in range(nturn - 1):
                if tr.n <= 0:
                    break
                self._seg(tr, 1, nlat)
                tr.next_turn()
            if loc > 1 and tr.n > 0:
                self._seg(tr, 1, loc)

    def track(self, loc, coords, nturn, record=False, alost=None, zlost=None, schedule=None, workers=1):
        """Track nturn turns from location loc; returns a dict with lost_turn, final coordinates and optionally per-turn records.
        位置 loc から nturn 周追跡する。戻り dict:
           lost_turn  (n,) 落ちた周 (1 始まり。0 は最後まで残った)
           final      (n, 6) 最後の座標 (落ちた粒子は落ちた時点の値)
           rec        (nturn, n, 6) 周ごとの座標 (record=True。落ちたあとは NaN。座標は位置 loc、1 周ごと)
        schedule: 'call' (既定、record=False のとき。TrackParticles[{loc, beam}, loc, 1, nturn] の 1 回の呼び出しと同じ区間の分け方)、
        'turn' (record=True のとき。1 周ずつ loc から回す = 1 周ずつ TrackParticles[{loc, beam}, loc] を呼ぶのと同じ)
        workers > 1: 粒子を塊に分けてスレッドで同時に追う (record=False のときだけ。CPU の論理コア数の 90% と 1 本あたり 16 粒子以上で抑える)。
        FLUC なしの結果は workers = 1 とビット単位で同じ。FLUC の乱数は塊ごとに種が違う"""
        c = np.atleast_2d(np.asarray(coords, dtype=float))
        n = len(c)
        sch = schedule or ('turn' if record else 'call')
        if record and sch != 'turn':
            raise ValueError('record=True は schedule=turn だけ')
        nw = 1 if record else _max_workers(workers, n)
        if nw > 1:
            if self.seed is not None:
                _t6.RNG.restore(self.seed)
            tr = _make_parallel_run(lambda c_, f_: self._run_of(c_, alost, zlost), c, [1.0] * n, nw)
        else:
            tr = self.make_run(c, alost, zlost)
        rec = np.full((nturn, n, 6), np.nan) if record else None
        if sch == 'call':
            self.run_call(tr, loc, nturn)
        else:
            for t in range(nturn):
                if tr.n <= 0:
                    break
                self.one_turn(tr, loc)
                if record and tr.n > 0:
                    rec[t, tr.kid[:tr.n]] = state_out(tr.S, tr.n)
        final, alive = tr.result()
        lt = np.array(tr.lost_turn[:n], dtype=np.int64)
        lost = np.where(lt == MAXTURN_LOST, 0, lt)
        return dict(lost_turn=lost, final=final, rec=rec, alive=np.array(alive, dtype=bool))

    # ------------------------------------------------------------------ 閉軌道
    def _map_batch(self, c, loc=1):
        tr = self.make_run(c, alost=1e30, zlost=1e30)
        self.one_turn(tr, loc)
        out, _a = tr.result()
        return out

    def refine_orbit(self, c0, dp=None, loc=1, tol=1e-12, itmax=8):
        """Compute the periodic solution (fixed point of the one-turn map) at loc by Newton's method (6D with RFSW ON, else 4D).
        位置 loc の周期解 (1 周の写像の不動点) を Newton 法で求める。RFSW ON (または放射つき) は 6 次元、RFSW OFF は δ = dp を決めた 4 次元。
        戻り: 閉軌道 (6)。残差 (1 周後との差) を self.info['orbit_resid'] に書く"""
        c = np.array(c0, dtype=float)
        six = self.rfsw
        if not six:
            c[4] = 0.0
            c[5] = 0.0 if dp is None else dp
        idx = list(range(6)) if six else [0, 1, 2, 3]
        h = 1e-8
        resid = None
        for _ in range(itmax):
            pts = [c.copy()]
            for k in idx:
                p = c.copy()
                p[k] += h
                pts.append(p)
            out = self._map_batch(np.array(pts), loc)
            f0 = out[0] - c
            resid = float(np.max(np.abs(f0[idx])))
            if resid < tol:
                break
            J = np.array([(out[1 + j] - out[0]) / h for j in range(len(idx))]).T[np.ix_(idx, idx)]
            dx = np.linalg.solve(J - np.eye(len(idx)), -f0[idx])
            for j, k in enumerate(idx):
                c[k] += dx[j]
        self.info['orbit_resid'] = resid
        return c

    def closed_orbit(self, locs):
        """Compute the closed orbit at the given locations by carrying one particle from cod1 at location 1 with radiation (no FLUC).
        位置 locs (1 始まりの整数のリスト) の閉軌道 (len(locs), 6)。位置 1 の cod1 から 1 粒子を放射 (FLUC なし) で運ぶ。
        RFSW OFF は δ=0 の 4 次元の閉軌道 (cod1 は δ=0)"""
        locs = [int(l) for l in locs]
        out = {}
        need = sorted(set(l for l in locs if l not in self._orbits))
        if need:
            save_fluc, self.fluc = self.fluc, False
            try:
                cur = 1
                c = self.cod1.copy()
                tr = self.make_run(c[None, :], alost=1e30, zlost=1e30)
                for l in need:
                    if l > cur:
                        with np.errstate(all='ignore'):
                            tr.segment(cur, l)
                        cur = l
                    self._orbits[l] = state_out(tr.S, tr.n)[0].copy() if l > 1 else c.copy()
            finally:
                self.fluc = save_fluc
        return np.array([self._orbits[l] for l in locs])

    def closed_orbit_dp(self, loc, dp):
        """Compute the 4D closed orbit at loc for delta = dp with RFSW OFF by Newton's method, solving at location 1 and carrying to loc.
        RFSW OFF: 位置 loc の δ = dp の閉軌道 (4 次元、Newton 法)。位置 1 で解いて loc へ運ぶ"""
        c1 = self.refine_orbit(np.zeros(6), dp=dp, loc=1)
        if loc == 1:
            return c1
        tr = self.make_run(c1[None, :], alost=1e30, zlost=1e30)
        with np.errstate(all='ignore'):
            tr.segment(1, loc)
        return state_out(tr.S, tr.n)[0]

    # ------------------------------------------------------------------ 位置
    def tunes_full(self):
        """Return the full tunes (nu_x, nu_y) including the integer part, from the Twiss phase advance over 2 pi.
        整数部を含む周回のチューン (νx, νy) (Twiss の位相の進み/2π)"""
        r = self.rows[-1]
        return r[7] / (2 * math.pi), r[8] / (2 * math.pi)

    def s_of(self, loc):
        return self.rows[loc - 1][0]

    def twiss_at(self, loc):
        r = self.rows[loc - 1]
        return dict(s=r[0], betx=r[1], alfx=r[2], bety=r[3], alfy=r[4], etx=r[5], etpx=r[6])

    # ------------------------------------------------------------------ バケット
    def bucket_zlimit(self, margin=1.0):
        """Return the |z| limit beyond which particles leave the RF bucket, (pi + 2|phi_s|)/k_rf; None with RFSW OFF or no RF.
        RF バケットから出た粒子を落とす |z| の限界: バケットの長さの半分 + 同期位相のずれ。(π + 2|φs|)/k_rf (バケットは安定固定点の周りに
        ψ = −π−2φs … π−2φs。φs = trf0·k_rf)。RFSW OFF か RF なしは None"""
        b = self.bucket
        if not self.rfsw or not b or b.get('k_rf', 0.0) <= 0.0:
            return None
        return margin * (math.pi + 2.0 * abs(b['phis'])) / b['k_rf']

    def nturn_for_synchrotron(self, nsynch=2.0, nmin=100):
        """Return the number of turns for nsynch synchrotron oscillation periods (period = 1/|nu_s|); None without RF.
        シンクロトロン振動 nsynch 周期分の周回数 (周期 = 1/|νs|)。RF なしは None"""
        b = self.bucket
        if not self.rfsw or not b:
            return None
        nus = abs(b.get('nus', 0.0))
        if nus <= 0.0:
            return None
        return max(int(nmin), int(math.ceil(nsynch / nus)))

    # ------------------------------------------------------------------ 条件の記録
    def conditions(self):
        b = self.bucket
        d = dict(compat_mode=self.mode, momentum_GeV=self.p0 * 1e-9, circumference_m=self.circ, nlat=self.nlat,
                 rfsw=self.rfsw, rad=self.rad, fluc=self.fluc, check_every=self.check_every, LOSSAMPL_m=self.alost, LOSSDZ_m=self.zlost,
                 orbit_residual=self.info.get('orbit_resid'))
        if b:
            d.update(bucket_height=b.get('height'), nu_s=b.get('nus'), phi_s_rad=b.get('phis'), U0_eV=b.get('U0'),
                     Vrf_eV=b.get('Vrf'), momentum_compaction=b.get('alpha'))
        return d


# ---------------------------------------------------------------------------------------------------------------------
# 周期性
# ---------------------------------------------------------------------------------------------------------------------

def _sig(el):
    items = []
    for k, v in el.__dict__.items():
        if k.startswith('_') or k == 'name':
            continue
        if isinstance(v, (int, float, np.floating)):
            items.append((k, round(float(v), 9)))
    return (type(el).__name__, tuple(sorted(items)))


def find_periodicity(els, rows=None):
    """Find the largest n for which the element sequence (ignoring zero-length markers, monitors, cavities) repeats n times.
    要素の並びの周期: マーカー・モニター・空洞 (L=0 の長さのない要素) を除いた並びが、同じ列の n 回の繰り返しになる最大の n を探す。
    戻り dict: nper (セルの数、1 なら周期的でない)、starts (各セルの最初の要素の位置 (1 始まり))、s_starts (その s)、cell_length (s の周期、
    rows (Twiss の行) を渡したとき)、ok (セルの長さが全部同じ)。同じ n でも、セルの中の要素の並び・磁場の値が全部同じことだけを調べる
    (誤差・補正は見ない。セルの鏡映の対称性は見ない)"""
    idx = [i for i, e in enumerate(els) if type(e).__name__ not in _SKIP]
    sigs = [_sig(els[i]) for i in idx]
    N = len(sigs)
    best = 1
    for n in range(N, 1, -1):
        if N % n:
            continue
        m = N // n
        if all(sigs[i] == sigs[i % m] for i in range(N)):
            best = n
            break
    m = N // best
    starts = [idx[k * m] + 1 for k in range(best)]
    out = dict(nper=best, starts=starts, nsig=N)
    if rows is not None:
        ss = [rows[l - 1][0] for l in starts]
        out['s_starts'] = ss
        circ = rows[-1][0]
        lens = [ss[k + 1] - ss[k] for k in range(best - 1)] + [circ - ss[-1] + ss[0]]
        out['cell_length'] = float(np.mean(lens))
        out['ok'] = bool(np.allclose(lens, lens[0], atol=1e-6))
    return out


def select_locations(setup, ls=1, le=None, step=1, ds=0.0, sort=True, unique_s=True):
    """Select starting locations for tracking (element entrances, 1-based) by range, step, minimum s spacing ds and unique s.
    追跡を始める位置 (要素の入口、1 始まり) の選び方。ls〜le (le は含む。既定は最後の要素)、step 要素おき、
    隣の選んだ位置と s が ds [m] 以上離れたものだけ (ds>0)、同じ s の位置 (長さ 0 の要素の列) は最初の 1 つだけ (unique_s)"""
    nel = setup.nlat - 1
    if le is None:
        le = nel
    le = min(int(le), nel)
    out = []
    last_s = None
    for l in range(max(1, int(ls)), le + 1, max(1, int(step))):
        s = setup.s_of(l)
        if unique_s and last_s is not None and abs(s - last_s) < 1e-9:
            continue
        if ds > 0.0 and last_s is not None and s - last_s < ds - 1e-12:
            continue
        out.append(l)
        last_s = s
    return out
