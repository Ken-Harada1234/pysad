# -*- coding: utf-8 -*-
"""SAD の台本から使う周波数マップ・局所運動量アクセプタンスの関数 (v2.3.0 Q1。K64 にない pysad の追加機能)。

  FrequencyMap[ls, nturn, X->{xmin, xmax, n}, Y->{...}, DP->{...}, opts]
      位置 ls (番号または名前) から (x, y) (または x–DP、y–DP) の格子を nturn 周追跡し、前半・後半のチューンから拡散 D を出す。
      X・Y・DP: 数 (固定)、{min, max, n} (n 点の等間隔)、または値のリスト (長さ ≠ 3 のとき)。単位は m・m・δ (閉軌道からのずれ)。
      オプション: Window->1 (NAFF の窓の次数)、Workers->1、Output->"file" (表)、Plot->"file.png" (図)、Order->4 (図の共鳴線の次数)。
      追跡の条件: v2.3.0 S2 から既定は RFSW OFF・RAD OFF・FLUC OFF (現在のフラグは使わない。使うには UseFlags->True。そのとき RF・放射・FLUC の警告が出る)。LOSSAMPL・LOSSDZ は現在の値。
      信頼できない・信頼度の低い条件 (周回数が少ない、失われた点が多い、作動点が共鳴に近い、など) の警告は標準出力に ` FrequencyMap: WARNING [コード] …` で出る。
      戻り: {X->…, Y->…, DP->…, NuX->…, NuY->…, NuX2->…, NuY2->…, Diffusion->…, LostTurn->…} (点ごとのリスト。NuX・NuY は前半、
      NuX2・NuY2 は後半の端数のチューン。失われた点は Diffusion = Indeterminate ではなく −99、LostTurn > 0)。

  MomentumAcceptance[locations, Turns->n, NSynch->2, ...]
      位置ごとの局所運動量アクセプタンス δ+(s)・δ−(s) (放射・RF つき)。locations は位置のリスト (番号・名前・パターン)。省略すると
      Range->{ls, le}、Step->k (k 要素おき)、Distance->ds (間隔 ds [m] 以上) で選ぶ要素の境界。Cell->True は、並びが n 回の繰り返しの
      リングの 1 周期 (最初のセル) の位置だけ。
      オプション: Turns (既定: NSynch 周期分、NSynch の既定は 2)、DPMax、Tolerance (1e-4)、Coarse (粗い走査の点数 12)、Sub (二分 1 回の点数 3。1 で二分法)、
      Bucket->True (RF バケットから出た粒子も失われたとする)、Workers (子プロセスの数、既定 2)、Aperture (|x|,|y| の限界。既定 LOSSAMPL)、CheckEvery (損失の判定の間隔、要素の数。**既定 1 = 要素ごと**。0 = SAD の tapert のまま (15 要素ごと、gs の TrackParticles と同じ判定。アパーチャが限界のとき δ を大きく見積もる))、
      Verify (限界の手前を n 点で再追跡して狭い損失の窓を探す。既定 0)、TurnCheck (限界の粒子を周回数のこの倍率で追跡し周回数依存を調べる。既定 2。0 で無効)、
      Output->"file" (表)、Plot->"file.png" (図)。
      戻り: {S->…, Loc->…, DPPlus->…, DPMinus->…} (δ− は負)。結果はセッションに残る (TouschekLifetimeLMA が使う)。

  TouschekLifetimeLMA[]    直前の MomentumAcceptance の結果を、INTRA ON の Emittance[] が作った TouschekTable に渡した Touschek 寿命 [s]
      (位置ごとの δ+(s)・δ−(s) から補間。Cell->True の結果は周期的に延ばす)。v2.3.0 T2: v2 の既定は Piwinski 型の式 (touschek.touschek_piwinski)、
      Formula->"K64" と k64 モードは K64 の表 (TouschekLifetime[INF, INF, nz] と同じ表)。
"""
import math
import os

import numpy as np

from .sadscript import Rule, Sym, S, SADError, is_num, truthy


def _opts(args):
    opts = {}
    rest = []
    for a in args:
        for r in (a if isinstance(a, list) else [a]):
            if isinstance(r, Rule) and isinstance(r.lhs, Sym):
                opts[r.lhs.name.upper()] = r.rhs
            else:
                rest.append(r)
    return opts, rest


def _axis(v, name):
    """X・Y・DP の指定 → 値のリスト"""
    if v is None:
        return [0.0]
    if is_num(v):
        return [float(v)]
    if isinstance(v, list):
        if len(v) == 3 and all(is_num(t) for t in v) and float(v[2]) >= 1 and float(v[2]) == int(v[2]) and float(v[1]) != float(v[0]) \
                and name != '_list':
            n = int(v[2])
            return [float(v[0])] if n == 1 else list(np.linspace(float(v[0]), float(v[1]), n))
        return [float(t) for t in v]
    raise SADError('FrequencyMap: %s must be a number, {min, max, n} or a list' % name)


def _num(opts, key, default):
    v = opts.get(key)
    return default if v is None or not is_num(v) else float(v)


def _str(opts, key):
    v = opts.get(key)
    return v if isinstance(v, str) else None


def install(ffs, B):
    state = {'lma': None, 'period': None}

    def _setup():
        from .trksetup import TrackSetup
        ffs.spread_pending() if hasattr(ffs, 'spread_pending') else None
        if ffs.lat is None:
            raise SADError('No beam line (USE first)')
        return TrackSetup.from_ffs(ffs)

    def _loc1(a):
        r = ffs._resolve_loc(a)
        if len(r) != 1:
            raise SADError('General::wrongtype: Argument must be one location in FrequencyMap')
        return int(r[0][0])

    def fn_frequency_map(args):
        from .freqmap import frequency_map
        opts, rest = _opts(args)
        if len(rest) < 2:
            raise SADError('FrequencyMap[location, nturn, X->{xmin, xmax, n}, Y->…, DP->…]')
        loc = _loc1(rest[0])
        nturn = int(rest[1])
        if truthy(opts.get('USEFLAGS')) is True:
            st = _setup()
        else:
            # v2.3.0 S2: 既定は RFSW OFF・RAD OFF・FLUC OFF (FFS に入ると RFSW は ON なので、台本に NORFSW を書き忘れると D が 0.03 → 0.57 の割合になる。S2 の記録)。現在のフラグを使うには UseFlags->True
            from .trksetup import TrackSetup
            ffs.spread_pending() if hasattr(ffs, 'spread_pending') else None
            if ffs.lat is None:
                raise SADError('No beam line (USE first)')
            st = TrackSetup.from_ffs(ffs, rad=False, rfsw=False)
            fl = ffs.flags
            if fl.get('RFSW') or fl.get('RAD'):      # FLUC は SAD の既定が ON (RAD OFF では効かない) なので、RFSW・RAD のどちらかが ON のときだけ知らせる
                ffs.write((' FrequencyMap: the current flags (RFSW=%s RAD=%s FLUC=%s) are not used: the map is made with RFSW OFF, RAD OFF, FLUC OFF. UseFlags->True uses them.\n') % (
                    'ON' if fl.get('RFSW') else 'OFF', 'ON' if fl.get('RAD') else 'OFF', 'ON' if fl.get('FLUC') else 'OFF'))
        res = frequency_map(st, loc, x=_axis(opts.get('X'), 'X'), y=_axis(opts.get('Y'), 'Y'), dp=_axis(opts.get('DP'), 'DP'), nturn=nturn,
                            window=int(_num(opts, 'WINDOW', 1)), workers=int(_num(opts, 'WORKERS', 1)))
        out = _str(opts, 'OUTPUT')
        if out:
            res.save(out)
        pl = _str(opts, 'PLOT')
        if pl:
            res.plot(pl, order=int(_num(opts, 'ORDER', 4)))
        ffs.write(' FrequencyMap: loc=%d (s=%.4f m) %d points, %d lost, nturn=%d, RFSW=%s RAD=%s\n' % (
            loc, st.s_of(loc), len(res), int((res.lost_turn > 0).sum()), res.cond['nturn'], 'ON' if st.rfsw else 'OFF', 'ON' if st.rad else 'OFF'))
        for (code, msg) in getattr(res, 'warnings', []):
            ffs.write(' FrequencyMap: WARNING [%s] %s\n' % (code, msg))
        D = np.where(np.isfinite(res.D), res.D, -99.0)
        return [Rule(S('X'), [float(v) for v in res.x]), Rule(S('Y'), [float(v) for v in res.y]), Rule(S('DP'), [float(v) for v in res.dp]),
                Rule(S('NuX'), [float(v) for v in res.nux1]), Rule(S('NuY'), [float(v) for v in res.nuy1]),
                Rule(S('NuX2'), [float(v) for v in res.nux2]), Rule(S('NuY2'), [float(v) for v in res.nuy2]),
                Rule(S('Diffusion'), [float(v) for v in D]), Rule(S('LostTurn'), [float(v) for v in res.lost_turn])]

    def fn_lma(args):
        from .lma import momentum_acceptance
        from . import trksetup as ts
        opts, rest = _opts(args)
        st = _setup()
        period = None
        if rest:
            locs = []
            for a in rest:
                locs.extend(int(k) for k, _fr in ffs._resolve_loc(a))
        else:
            rg = opts.get('RANGE')
            ls, le = 1, None
            if isinstance(rg, list) and len(rg) == 2:
                ls, le = _loc1(rg[0]), _loc1(rg[1])
            if truthy(opts.get('CELL')) is True:
                per = ts.find_periodicity(st.els, st.rows)
                from .lma import periodicity_report
                rep = periodicity_report(st, warn=False)
                if rep['periodic']:
                    ls, le = per['starts'][0], per['starts'][1] - 1
                    period = (per['cell_length'], per['s_starts'][0])
                else:
                    ffs.write(' MomentumAcceptance: the beam line is not periodic (%s); Cell->True ignored, the whole ring is used\n' % '; '.join(rep['messages']))
                    import warnings as _w
                    _w.warn('MomentumAcceptance: Cell->True ignored: ' + '; '.join(rep['messages']))
            locs = ts.select_locations(st, ls, le, int(_num(opts, 'STEP', 1)), _num(opts, 'DISTANCE', 0.0))
        nt = opts.get('TURNS')
        bk = opts.get('BUCKET')
        res = momentum_acceptance(st, locs, nturn=int(nt) if is_num(nt) else None, nsynch=_num(opts, 'NSYNCH', 2.0),
                                  dpmax=_num(opts, 'DPMAX', None) if is_num(opts.get('DPMAX')) else None,
                                  ncoarse=int(_num(opts, 'COARSE', 12)), nsub=int(_num(opts, 'SUB', 3)), tol=_num(opts, 'TOLERANCE', 1e-4),
                                  bucket=True if bk is None else truthy(bk) is True, workers=int(_num(opts, 'WORKERS', 2)),
                                  alost=_num(opts, 'APERTURE', None) if is_num(opts.get('APERTURE')) else None,
                                  check_every=int(_num(opts, 'CHECKEVERY', 1)) if is_num(opts.get('CHECKEVERY')) else None,
                                  verify=int(_num(opts, 'VERIFY', 0)), turncheck=_num(opts, 'TURNCHECK', 2.0))
        state['lma'] = res
        state['period'] = period
        out = _str(opts, 'OUTPUT')
        if out:
            res.save(out)
        pl = _str(opts, 'PLOT')
        if pl:
            res.plot(pl)
        ffs.write(' MomentumAcceptance: %d positions, nturn=%d, min dp+ = %.5f, min dp- = %.5f\n' % (
            len(res), res.cond['nturn'], res.dp_plus.min(), res.dp_minus.max()))
        return [Rule(S('S'), [float(v) for v in res.s]), Rule(S('Loc'), [float(v) for v in res.loc]),
                Rule(S('DPPlus'), [float(v) for v in res.dp_plus]), Rule(S('DPMinus'), [float(v) for v in res.dp_minus])]

    def fn_touschek_lma(args):
        from .lma import touschek_lifetime_local, touschek_lifetime_local_piwinski
        from . import compat as _compat
        opts, _rest = _opts(args)
        fm = opts.get('FORMULA')
        explicit = isinstance(fm, str)
        fm = fm.lower() if explicit else ('k64' if _compat.k64() else 'piwinski')      # v2.3.0 T2: v2 の既定は Piwinski 型の式 (docs/history/v23_T2.md)
        if fm not in ('k64', 'piwinski'):
            raise SADError('TouschekLifetimeLMA: Formula must be "K64" or "Piwinski"')
        tt = getattr(ffs, 'touschek', None)
        if tt is None or tt.get('table') is None:
            raise SADError('Touschek::undef:  TouschekTable Undefined; Do Emittance[] with INTRA;')
        if state['lma'] is None:
            raise SADError('TouschekLifetimeLMA: run MomentumAcceptance[] first')
        per = state['period']
        kw = dict(period=per[0] if per else None, s0=per[1] if per else 0.0)
        if fm == 'piwinski' and 'rest' in tt:
            r = touschek_lifetime_local_piwinski(tt, state['lma'], **kw)
        else:
            r = touschek_lifetime_local(tt['table'], state['lma'], **kw)
        return float(r['lifetime'])

    B['FrequencyMap'] = fn_frequency_map
    B['MomentumAcceptance'] = fn_lma
    B['TouschekLifetimeLMA'] = fn_touschek_lma
