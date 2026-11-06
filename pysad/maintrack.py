"""
maintrack.py — SAD の MAIN レベルの TRACK コマンド (toplvl → ActTra.f → track.f → tracka.f) の直訳

SAD の経路
    TRACK USE=RING NPART=21 TURNS=1000 NX=(..) NY=(..) NZ=(..);        (src/initbl.f が icACT として登録)
      ActTra.f : 引数 (USE, NPART, TURNS, NX, NY, NZ, 粒子種) を読み、粒子 j (0 始まり) の規格化振幅を
                 x = NX[j mod len(NX)], y = NY[j mod len(NY)], g(dp) = NZ[j mod len(NZ)]  (px = py = z = 0) とする。
                 リストが無ければ x = y = 1, g = 0。track() の後、NX=/NY=/NZ= のリストを autofg(S7.4) で 7 桁ずつ並べて表示。
      track.f  : フラグの一覧 (RFSW RADCOD RAD FLUC INTRA POL COD DAPER EMIOU K64 FOURI SMEAR) を表示し、tracka を呼ぶ。
                 追跡のあと `Tracking end:  CPU time = ...` を表示。TURNS=0 なら `Use EMIT_TANCE or Emittance[] within FFS.`
      tracka.f : emx = sqrt|EMITX+EMITY|、emz = sqrt|EMITZ| (RFSW) か |SIGE|。規格化座標 (x px y py z g) × (emx emx emx emx emz emz) に
                 emit の正規モード行列 r と閉軌道 codin を掛けて足して初期座標 (trackd.f tinip。TDXI/TDYI/TDZI を x,y,z に足す)。
                 周 n = 1..TURNS+1 について、周の頭で ri·(x − codin) (規格化座標) を zbuf(n, 3 面, 粒子) に (x+i px, y+i py, z+i g) と
                 して貯め、2J の 1 次・2 次和 (sjx, sjxjx …) を積算 (SMEAR ON のとき zbuf も貯める)。SMEAR ON のとき
                 n が nsmear の倍数か n=TURNS+1 なら tsmear (周ごとの 2J とそのゆらぎ) を表示。そのあと tturn で 1 周追跡
                 (n=TURNS+1 では追跡しない)。粒子が全部落ちたら終了。
                 追跡後、PHSPLOTS (TopDrawer 出力: 描画なので無視)、SMEAR ON なら wjfit (FOURIE OFF → wtunem: 面ごとに FFT の
                 ピークを 8 本、周波数の補間つき) を表示し、大域 `TunesByTracking` に結果を入れる。
      tapert   : 落ちた粒子は `P. id lost in n-th turn at l(name), amplitudes: x y z` を表示 (tapert_report_dropped)。
      tsmear.f / wjfit.f / wtunem.f / tcftr.f: smear_report / tune_table / wtunem 参照。

K64 では MAIN の EMIT と TRACK は動かない (EMIT は `Use EMIT_TANCE or Emittance[] within FFS.` と言うだけ、TRACK は emit の r が
未設定で SIGSEGV)。ここでは古典 SAD の意味 (EMIT が r, codin, EMITX.. を作り TRACK が使う) で実装し、追跡本体は
TrackParticles と同じ track6n.turn (SAD の tturn の直訳) を使う。FOURIE ON (wjfit の maxord>0) と PHSPLOTS の描画は未対応/無視。
"""
from __future__ import annotations

import math
import os
import time
from typing import List, Optional

import numpy as np

from . import track6 as _t6
from . import track6n as _tn
from .emit6d import _Phys

NFRE = 8                       # wjfit.f: nfre (面ごとに残すピークの数)
_MAXTURN = 2**29


# ---------------------------------------------------------------------------- Fortran の書式
def _fit(body: str, w: int) -> str:
    if len(body) > w:                       # gfortran: 幅が足りなければ先頭の 0 を落とす
        if body.startswith('0.'):
            body = body[1:]
        elif body.startswith('-0.'):
            body = '-' + body[2:]
    if len(body) > w:
        return '*' * w
    return body.rjust(w)


def fmt_g(x: float, w: int, d: int) -> str:
    """Fortran の 1P,Gw.d (gfortran と照合済み): 0.1 <= |x| < 10^d は F(w-4).(d-指数) + 4 空白、他は E 形式"""
    if x != x:
        return 'NaN'.rjust(w)
    if math.isinf(x):
        return ('Infinity' if x > 0 else '-Infinity').rjust(w)
    if x == 0.0:
        return _fit('0.' + '0' * (d - 1), w - 4) + ' ' * 4
    s = f'{abs(x):.{d - 1}e}'
    ex = int(s.split('e')[1])
    m = float(s)
    if 0.1 <= m < 10.0 ** d:
        nd = d - (ex + 1)
        body = f'{x:.{nd}f}' + ('.' if nd == 0 else '')
        return _fit(body, w - 4) + ' ' * 4
    mant, e = f'{x:.{d}E}'.split('E')
    en = int(e)
    es = f'E{en:+03d}' if abs(en) < 100 else f'{en:+04d}'
    return _fit(mant + es, w)


def fmt_f(x: float, w: int, d: int) -> str:
    if x != x:
        return 'NaN'.rjust(w)
    return _fit(f'{x:.{d}f}', w)


def ordinal(m: int) -> str:
    """tday.f の ord"""
    mm = m % 100
    if 3 < mm < 21:
        return 'th'
    mm = m % 10
    return {1: 'st', 2: 'nd', 3: 'rd'}.get(mm, 'th')


# ---------------------------------------------------------------------------- 解析 (tsmear / wtunem)
def smear_report(n: int, nsmear: int, sj, sjj, alive_ids, np0: int):
    """tsmear.f。sj/sjj は粒子 (ID 順) ごとの 2J の和と 2 乗和 (3 面)、alive_ids は生きている粒子の ID (0 始まり、昇順)。
    表示は ID 昇順 (kptbl(i,1) <= np の粒子)。戻り: (行のリスト, [(ID(1 始まり), [2Jx, 2Jy, 2Jz, σx, σy, σz])])"""
    out = []
    recs = []
    npa = len(alive_ids)
    if npa > 0:
        out.append(f' Smear at {n:6d}th turn:')
        out.append(' p. #   2Jx(m)      2Jy(m)      2Jz(m)    ' + '  sigma2Jx  sigma2Jy  sigma2Jz')
    aj = [0.0, 0.0, 0.0]
    for i in alive_ids:
        s = [sj[k][i] / nsmear for k in range(3)]
        sg = []
        for k in range(3):
            if s[k] > 0.0:
                sg.append(math.sqrt(abs(sjj[k][i] / nsmear - s[k] ** 2)) / s[k])
            else:
                sg.append(0.0)
        for k in range(3):
            aj[k] += s[k]
        out.append(f'{i + 1:6d}' + ''.join(fmt_g(v, 12, 4) for v in s) + ''.join(fmt_f(v, 10, 6) for v in sg))
        recs.append((int(i) + 1, s + sg))
    out.append(' Avrg ' + ''.join(fmt_g(v / npa if npa else float('nan'), 12, 4) for v in aj))
    return out, recs


def wtunem(z: np.ndarray, nfre: int = NFRE):
    """wtunem.f: z は (lsp, 3) の複素座標 (x+i px, y+i py, z+i g)。面ごとに FFT (tcftr: exp(+2πi nk/N)、正規化なし) のパワーを
    重み (.5, 1, .5) で平滑し、その極大の位置を放物線補間した周波数と振幅 (fa = r/lsp²) を大きい順に nfre 本。
    戻り f1[3][nfre], fa[3][nfre] (周波数は [-.5, .5] の tune 単位)"""
    lsp = z.shape[0]
    f1 = np.zeros((3, nfre))
    fa = np.zeros((3, nfre))
    for k in range(3):
        X = np.fft.ifft(z[:, k]) * lsp                  # tcftr(.., conj=.false.) は exp(+i...) の和
        P = X.real ** 2 + X.imag ** 2
        r = P + 0.5 * (np.roll(P, 1) + np.roll(P, -1))  # r(i) = P(i) + (P(i-1)+P(i+1))/2 (wtunem.f 66-76 行目の積算の整理)
        rm = np.roll(r, 1)
        rp = np.roll(r, -1)
        a = rm + rp - 2.0 * r
        cand = np.nonzero(a < 0.0)[0]
        fk = f1[k]
        ak = fa[k]
        for idx in cand:                                 # idx は 0 始まり (Fortran の i-1)
            aa = a[idx]
            xi = -(rp[idx] - rm[idx]) * 0.5 / aa
            if not (-0.5 <= xi < 0.5):
                continue
            ri = r[idx] + xi * (rp[idx] - rm[idx]) * 0.25
            fi = idx + xi
            merged = False
            for j in range(nfre):
                if abs(fk[j] - fi) < 1.0:
                    if ri > ak[j]:
                        ak[j] = ri
                        fk[j] = fi
                    merged = True
                    break
            if merged:
                continue
            for j in range(nfre):
                if ri > ak[j]:
                    for l in range(nfre - 1, j, -1):
                        ak[l] = ak[l - 1]
                        fk[l] = fk[l - 1]
                    ak[j] = ri
                    fk[j] = fi
                    break
        ak *= 1.0 / lsp ** 2
        fk *= 1.0 / lsp
        fk[fk > 0.5] -= 1.0
    return f1, fa


def tune_table(zbuf: np.ndarray, alive_mask: np.ndarray, lines: List[str]):
    """wjfit.f (maxord=0): 粒子ごとに wtunem の結果を表示し、TunesByTracking 用のリストを返す。
    zbuf は (np0, lsp, 3)、alive_mask は ID ごとの生存。"""
    res = []
    for i in range(zbuf.shape[0]):
        tag = ':     ' if alive_mask[i] else ',lost:'
        f1, fa = wtunem(zbuf[i])
        lines.append(f' P.{i + 1:3d}{tag}' + ' nux       ampl        nuy       ampl        nuz       ampl')
        for k in range(NFRE):
            lines.append(f'{k + 1:6d}  ' + ''.join(fmt_f(f1[j, k], 10, 6) + fmt_g(fa[j, k], 12, 4) for j in range(3)))
        res.append([1.0 if alive_mask[i] else 0.0] + [[[float(f1[j, k]), float(fa[j, k])] for k in range(NFRE)] for j in range(3)])
    return res


# ---------------------------------------------------------------------------- 本体
def initial_distribution(np0, nx, ny, nz, rmat, codin, emx, emz, tdx=0.0, tdy=0.0, tdz=0.0):
    """ActTra.f の規格化振幅 + tracka.f/tinip の変換: 物理座標 (np0, 6)。
    j 番目の粒子の振幅 = (x, px, y, py, z, g) = (NX[j mod n], 0, NY[j mod n], 0, 0, NZ[j mod n])。リストが無ければ x=y=1, g=0"""
    amp = np.zeros((np0, 6))
    for j in range(np0):
        amp[j, 0] = nx[j % len(nx)] if nx else 1.0
        amp[j, 2] = ny[j % len(ny)] if ny else 1.0
        amp[j, 5] = nz[j % len(nz)] if nz else 0.0
    scale = np.array([emx, emx, emx, emx, emz, emz])
    coords = (amp * scale) @ np.asarray(rmat, dtype=float).T + np.asarray(codin, dtype=float)
    coords[:, 0] += tdx
    coords[:, 2] += tdy
    coords[:, 4] += tdz
    return coords


def echo_lists(nx, ny, nz, fmt_real) -> List[str]:
    """ActTra.f の最後: NX= / NY= / NZ= のリストを autofg(x,'S7.4') を 7 桁ずつにして 70 桁を超えたら改行して表示"""
    lines = []
    for tag, lst in (('NX=', nx), ('NY=', ny), ('NZ=', nz)):
        mess = ''
        for v in (lst or []):
            mess += fmt_real(float(v), 'S7.4')[:7].ljust(7)
            if len(mess) >= 70:                        # Fortran: j = len+1 > 70
                lines.append(' ' + tag + mess)
                mess = ''
        if mess:
            lines.append(' ' + tag + mess)
    return lines


def run_track(ffs, spec: dict, write):
    """MAIN の TRACK を 1 回実行する。ffs: ffscmd.FFS (フラグ・emit・要素・大域変数を持つ)、write: 出力関数"""
    from .sadscript import fmt_real, SADError
    g = ffs.I.globals
    fl = ffs.flags
    name = spec.get('use')
    if name and name.upper() != ffs.name.upper():
        ffs.use(name)
    if ffs.lat is None:
        raise SADError('No beam line (TRACK USE=)')
    ffs.spread_pending()
    np0 = spec['npart'] if spec.get('npart') is not None else 0
    nturn = spec['turns'] if spec.get('turns') is not None else 0
    g['NP'] = float(np0)
    g['TURNS'] = float(nturn)
    ptype = spec.get('ptype')
    charge = -1.0 if ptype == 'ELEC' else 1.0
    amass = None
    if ptype == 'PROTON':
        amass = 938.272013e6
    # RFSW は MAIN の $RFSW$ (track.f:80, 107 isynch; FFS の RFSW とは別)。TURNS=-1 なら rfsw=.true. (track.f:160)
    rfsw = bool(ffs.main_flags.get('RFSW', False)) or nturn == -1
    rad, fluc, smearp = bool(fl['RAD']), bool(fl['FLUC']), bool(fl['SMEAR'])
    vals = [rfsw, fl['RADCOD'], fl['RAD'], fl['FLUC'], fl['INTRA'], fl['POL'], fl['COD'], fl['DAPERT'],
            fl['EMIOUT'], fl['K64'], fl.get('FOURIER', False), fl['SMEAR']]
    write(' RFSW RADCOD RAD  FLUC  INTRA POL   COD  DAPER EMIOU K64   FOURI SMEAR\n')
    write(' ' + ''.join(f'{"T" if v else "F":>3}   ' for v in vals) + '\n')
    if nturn == 0:
        write(' Use EMIT_TANCE or Emittance[] within FFS.\n')
        return None
    if np0 <= 0:
        raise SADError('TRACK: NPART must be positive')
    nturn = abs(nturn)
    t_start = time.time()
    # emit の結果 (正規モード行列 r と閉軌道 codin)。MAIN の EMIT が作る。まだなら黙って作る
    if ffs.emit is None:
        ffs.run_emit(quiet=True)
    emit = ffs.emit
    rmat = np.asarray(emit['normal_r'], dtype=float)
    ri = np.asarray(emit['normal_ri'], dtype=float)
    codin = np.asarray(emit['closed_orbit'], dtype=float)
    emx = math.sqrt(abs(float(g.get('EMITX', 0.0)) + float(g.get('EMITY', 0.0))))
    if rfsw:
        emz = math.sqrt(abs(float(g.get('EMITZ', 0.0))))
        if emz <= 0.0:
            emz = abs(float(g.get('SIGE', 0.0)))
    else:
        emz = abs(float(g.get('SIGE', 0.0)))
    coords = initial_distribution(np0, spec.get('nx'), spec.get('ny'), spec.get('nz'), rmat, codin, emx, emz,
                                  float(g.get('TDXI', 0.0)), float(g.get('TDYI', 0.0)), float(g.get('TDZI', 0.0)))
    els = ffs.els
    circ = sum(el.L for el in els)
    kw = {} if amass is None else {'amass': amass}
    ph = _Phys(ffs.momentum(), circ, charge=charge, **kw)
    ph.dvfs = ffs._fshift() if rfsw else 0.0                   # tsetdvfs
    try:
        prog = _tn.Prog(els, ph)
    except NotImplementedError as e:
        raise SADError(f'TRACK: unsupported element ({e})')
    S = _tn.make_state(coords, ph)
    n = S.shape[0]
    kid = np.arange(np0, dtype=np.int64)                 # スロット → 粒子 ID (落ちた粒子の入れ替えで変わる)
    lost_turn = np.full(np0, _MAXTURN, dtype=np.int64)
    log = np.zeros((np0, 6))
    nlog = np.zeros(1, dtype=np.int64)
    rs = np.array([_t6.RNG.seed], dtype=np.int64)
    flc = bool(rad and fluc)
    alost = float(g.get('LOSSAMPL', 1.0))
    zlost = float(g.get('LOSSDZ', 100.0))
    trf0 = ffs.trf0
    # zbuf の長さ nts (track.f): nturn+1 以上の最小の 2 のべき
    nts = 2 ** int(math.log(nturn + 1) / math.log(2.0))
    if nts < nturn + 1:
        nts *= 2
    if nts >= 2 * (nturn + 1):
        nts //= 2
    zbuf = np.zeros((np0, nts, 3), dtype=complex) if smearp else None
    nsmear = max(100, nturn // 50) // 100 * 100
    sj = [np.zeros(np0), np.zeros(np0), np.zeros(np0)]
    sjj = [np.zeros(np0), np.zeros(np0), np.zeros(np0)]
    names = [el.name for el in els]
    nlat = len(els) + 1
    smear_rec = []                                       # (n, ID, [2Jx 2Jy 2Jz σx σy σz]) 検証用の全桁データ
    loss_rec = []                                        # (ID, 周, l) 損失の記録
    progress = int(float(g.get('PYSAD$PROGRESS', 0.0))) or int(os.environ.get('PYSAD_PROGRESS', '0') or 0)
    with np.errstate(all='ignore'):
        for nn in range(1, nturn + 2):
            if n > 0:
                # 周の頭の規格化座標 ri·(x − codin) (tracka.f 7020): 生きている粒子だけ
                c = _tn.state_out(S, n) - codin
                nrm = c @ ri.T
                ids = kid[:n]
                z3 = np.stack([nrm[:, 0] + 1j * nrm[:, 1], nrm[:, 2] + 1j * nrm[:, 3], nrm[:, 4] + 1j * nrm[:, 5]], axis=1)
                if smearp:
                    zbuf[ids, nn - 1, :] = z3
                for k in range(3):
                    aj = nrm[:, 2 * k] ** 2 + nrm[:, 2 * k + 1] ** 2
                    sj[k][ids] += aj
                    sjj[k][ids] += aj * aj
                if smearp and (nn % nsmear == 0 or nn == nturn + 1):
                    nsm = (nn - 1) % nsmear + 1
                    if nsm > 1:
                        lns, recs = smear_report(nn, nsmear, sj, sjj, np.sort(ids), np0)
                        smear_rec.extend((nn, i1, v) for i1, v in recs)
                        for ln in lns:
                            write(ln + '\n')
                    for k in range(3):
                        sj[k][:] = 0.0
                        sjj[k][:] = 0.0
            if nn > nturn:
                break
            nlog[0] = 0
            n = _tn.turn(prog, S, n, ph, trf0, rfsw, rad, flc, rs, apert=(kid, lost_turn, nn, alost, zlost, log, nlog))
            for r in range(int(nlog[0])):
                pid, l, xv, yv, zv, t = log[r]
                l = int(l)
                t = int(t)
                loss_rec.append((int(pid) + 1, t, l))
                nm = names[l - 1] if 1 <= l < nlat else '$$$'
                write(f' P. {int(pid) + 1:5d} lost in {t:5d}{ordinal(t)} turn at {l:5d}({nm}), amplitudes:'
                      + ''.join(' ' + fmt_g(v, 13, 7) for v in (xv, yv, zv)) + '\n')
            if n <= 0:
                break
            if not rfsw:
                S[:n, 4] = 0.0
            if progress and (nn % progress == 0 or nn == nturn):
                print(f'    [track] turn {nn}/{nturn} alive {n}/{np0} ({time.time() - t_start:.1f} s)', flush=True)
    _t6.RNG.seed = int(rs[0])
    alive = np.zeros(np0, dtype=bool)
    alive[kid[:max(n, 0)]] = True
    tunes = None
    if smearp:
        lines: List[str] = []
        tunes = tune_table(zbuf, alive, lines)
        g['TunesByTracking'] = tunes
        for ln in lines:
            write(ln + '\n')
    write(' ' + 'Tracking'.ljust(20) + ' end:  CPU time =' + fmt_f(time.time() - t_start, 10, 3) + ' sec\n')
    for ln in echo_lists(spec.get('nx'), spec.get('ny'), spec.get('nz'), fmt_real):
        write(ln + '\n')
    final = np.full((np0, 6), np.nan)                    # 落ちた粒子は落ちた時点の座標 (limittrack で切った値) のまま末尾のスロットに残る
    final[kid] = _tn.state_out(S, np0)
    return {'alive': int(max(n, 0)), 'lost_turn': lost_turn, 'final': final, 'smear': smear_rec, 'tunes': tunes, 'losses': loss_rec,
            'emx': emx, 'emz': emz, 'codin': codin}
