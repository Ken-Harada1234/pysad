"""
emit_extra.py — SAD Emittance[] の戻り規則一式 (OneTurnInformation / Orbit / Matrix) の直訳
==========================================================================================

emit6d.py (temit/tecalc の直訳) が返すのは Emittances, MomentumSpread, BunchLength, Tunes,
DampingRate, EquilibriumBeamMatrix, NormalCoordinates, OrbitAtExit … の一部だけだった。
本モジュールは Packages/emit.n の Emittance[opt] が組み立てる規則のリストを **SAD の名前のまま**
dict で返す。既存ファイル (emit6d.py, calc6d.py, lattice.py) は変更せず、そこから private 関数を
借りて使う。

SAD 側の対応 (Packages/emit.n, tfemit.f, temit.f, tturne.f, tffs.f):
  mode = 3 (Matrix) / 2 (Orbit) / 1 (OneTurnInformation) / 0 (既定)
  e[[1]] = stab, e[[2]] = params(nparams) (setparams + tecalc), e[[3]] = setiamat の 7 要素
    (codin, trans(:,1:6), trans(:,7:12), ri, beamp(21), emitp(21), srot(3,3))
  e[[4]] = 各要素入口の cod (tturne1 の plot 分岐, nlat 個)
  e[[5]] = 各要素入口までの累積転送行列 trans(:,1:6) (同上)
  また CALC6D のとき LINE["SIG",l] / Twiss["SIGX",l] は tfsize → temit(plot, beamplt) が
  各要素入口に残す beamsize(:,l) (= 入口の平衡ビーム行列 Σ(s)) を返す。temit の最後の
  「plot .and. postcal」の一周 (beam(1:21)=emitp から tturne) がこれを作る (temit.f:214-222)。

規則の対応表 (キーは SAD の規則名そのまま):
  Stable                 |Re log λ|<1e-6 (3 モード) かつ閉軌道が見つかった
  Region                 {1, nlat} (Region 指定は未対応)
  Tunes                  Im(cd(4:6))/2π
  EnergyLossU0           u0·pgev [eV]
  RfVoltageVc            vceff [V]
  EquilibriumPosition    trf0 [m]
  MomentumCompaction     alphap
  OrbitDilation          dleng = FSHIFT·circ
  BucketHeight           bh
  HarmonicNumber         heff = wrfeff·c/ω0
  OrbitAtExit            一周後の cod
  DampingRate            Re cd(1:3)
  Emittances             (εx, εy, εz)
  MomentumSpread         σE,  BunchLength σz
  TuneShiftByRadiation   Im cd(1:3)/2π
  NominalSpinTune        h0·gspin (ipnnup)
  PolarizationTime, SpinTune, Polarization, Polarizations, PolarizationVector
                         偏極計算 (calpol) は未移植。SAD は params=0 で初期化するので 0 を返す
  OrbitAtEntrance        codin
  OneTurnTransferMatrix  T = trans(:,1:6) (放射なしのシンプレクティック部分。emit6d の one_turn_matrix は T+D)
  OneTurnDampingMatrix   D = trans(:,7:12)
  NormalCoordinates      ri (物理 → 正規モード)
  ExtendedTwissParameters  dict {AX, BX, PSIX, …, ZPY} = tfetwiss(ri, OrbitAtExit)
  OneTurnExcitation      beamp: beam=0 から一周したときの量子励起 (物理座標, 6×6 対称)
  EquilibriumBeamMatrix  emitp (6×6)
  SpinPrecessionMatrix   srot(:,1:3)。calpol 未移植なので単位行列 (srotinit の初期値)
  ClosedOrbit            [cod_1, …, cod_nlat] (各要素入口。nlat = len(elements)+1 は $$$)
  TransferMatrices       [T_1, …, T_nlat] (入口までの累積 T。T_1 = 単位行列)
  BeamMatrices           [Σ_1, …, Σ_nlat] (各要素入口の平衡ビーム行列 = SAD の LINE["SIG",l])
  IntrabeamExcitation    INTRA 未対応 (touschek.py 参照) → 入れない

未移植: SynchroBetaEmittance (temits.f)、Region (部分区間)、InitialOrbit/InitialBeamMatrix、TRPT 分岐。

接続方法 (既存ファイルは触らない):
  from pysad.emit_extra import emittance_rules
  e = emittance_rules(lat, 3.0e9, rfsw=True, radcod=False)   # Emittance[OneTurnInformation->True, Matrix->True] 相当
  e['Emittances'], e['ExtendedTwissParameters']['BX'], e['TransferMatrices'][k-1], e['BeamMatrices'][k-1]
  ffscmd.py の Emittance[] 実装 (emittance6d 呼び出し) をこちらに差し替えるには、
  emittance_rules(..., matrix=False) の戻り dict を SADScript の規則リストに変換すればよい
  (キー名が SAD と同じなので Rule[Symbol[key], value] を並べるだけ)。
"""

import math
from typing import List, Optional

import numpy as np

from . import lattice as _lat
from .elements import Element, Drift, Quadrupole, Bend
from .emit6d import (_Phys, _State, _tdrife0, _tbende, _tquade, _tcave, _tmulte, _tthine,
                     _rf_effective, _normal_form, _tsymp, _tinv6, tecalc, CVELOC, ELMASS)
from .calc6d import _closed_orbit, tfetwiss

# tfstk.f:269 gspin (電子の異常磁気能率 a)
GSPIN = 0.001159652181280002

# Packages/init.n FFS$TwissFuns の順 (= tffs.f の mfitax..mfitzpy の順)
TWISS_FUNS = ('AX', 'BX', 'PSIX', 'AY', 'BY', 'PSIY', 'EX', 'EPX', 'EY', 'EPY',
              'R1', 'R2', 'R3', 'R4', 'DETR', 'DX', 'DPX', 'DY', 'DPY', 'DZ', 'DDP',
              'AZ', 'BZ', 'PSIZ', 'ZX', 'ZPX', 'ZY', 'ZPY')

# emit.n の Emittance[] が返す規則の順 (mode=3 のとき。参考用)
RULE_ORDER = ('Stable', 'Region', 'Tunes', 'EnergyLossU0', 'RfVoltageVc', 'EquilibriumPosition',
              'MomentumCompaction', 'OrbitDilation', 'BucketHeight', 'HarmonicNumber', 'OrbitAtExit',
              'DampingRate', 'Emittances', 'MomentumSpread', 'BunchLength', 'TuneShiftByRadiation',
              'NominalSpinTune', 'PolarizationTime', 'SpinTune', 'Polarization', 'Polarizations',
              'PolarizationVector',
              'OrbitAtEntrance', 'OneTurnTransferMatrix', 'OneTurnDampingMatrix', 'NormalCoordinates',
              'ExtendedTwissParameters', 'OneTurnExcitation', 'EquilibriumBeamMatrix',
              'SpinPrecessionMatrix', 'ClosedOrbit', 'TransferMatrices')


# ---------------------------------------------------------------------------
# 一周追跡 (emit6d.tturne / calc6d._turn と同じ要素分岐 + 要素入口ごとに T, cod, beam を記録)
# ---------------------------------------------------------------------------

def _turn_record(elements: List[Element], cod0, ph: _Phys, trf0: float = 0.0,
                 rfsw: bool = False, radcod: bool = False,
                 beam0: Optional[np.ndarray] = None, record: bool = True) -> _State:
    """tturne1 (plot=.true.) 相当。
    beam0 : 入口のビーム行列 (6×6)。None なら 0 (temit の 1 周目 beam=beamin=0)。
            temit の最後の plot 一周は beam(1:21)=emitp から始める (temit.f:217)。
    record=True なら各要素入口 l=1..nlat で (T_{1→l}, cod_l, beam_l) を st.rec に積む
    (tturne1 の `if(plot)` の時点 = 要素 l を通す前。l=nlat は一周後 = tfsetplot(nlat))。"""
    st = _State(cod0, ph)
    st.trf0 = trf0
    st.rfsw = rfsw
    st.radcod = radcod
    if beam0 is not None:
        st.B = np.array(beam0, dtype=float).copy()
    st.setr0()
    rec = []
    for el in elements:
        if record:
            rec.append((st.T.copy(), st.cod.copy(), st.B.copy()))
        st.setr0()                       # tturne1: 要素ごとに tsetr0
        if isinstance(el, Drift):
            _tdrife0(st, el.L)
        elif isinstance(el, Bend):
            _tbende(st, el, enarad=True)
        elif isinstance(el, Quadrupole):
            _tquade(st, el, krad=True)
        elif type(el).__name__ in _lat._THIN_KORD:
            _tthine(st, el)          # tthine (放射なし)。emit6d.tturne と同じ
        elif type(el).__name__ == 'Cavity':
            _tcave(st, el)
        elif type(el).__name__ == 'Multipole':
            _tmulte(st, el, enarad=True)
        elif el.L == 0.0:
            continue          # MARK, MONI など
        else:
            raise NotImplementedError(f"emit_extra: 要素 {type(el).__name__} は未対応")
    if record:
        rec.append((st.T.copy(), st.cod.copy(), st.B.copy()))
    st.rec = rec
    return st


# ---------------------------------------------------------------------------
# Emittance[opt] の規則一式
# ---------------------------------------------------------------------------

def emittance_rules(lattice, momentum_eV: float, rfsw: bool = False, radcod: bool = False,
                    charge: float = 1.0, amass: float = ELMASS,
                    one_turn: bool = True, orbit: bool = True, matrix: bool = True,
                    beam_matrices: Optional[bool] = None, fshift: float = 0.0) -> dict:
    """
    SAD の Emittance[OneTurnInformation->one_turn, Orbit->orbit, Matrix->matrix] の戻り規則を
    dict (キー = SAD の規則名) で返す。emit.n と同じく Matrix ⊃ Orbit ⊃ OneTurnInformation。

    Parameters
    ----------
    lattice      : pysad.lattice.Lattice
    momentum_eV  : 運動量 [eV/c]
    rfsw, radcod : FFS フラグ (RAD/FLUC は Emittance[] に影響しない: temit は常に irad=12)
    one_turn / orbit / matrix : emit.n の mode (1 / 2 / 3)
    beam_matrices: True なら 'BeamMatrices' (各要素入口の Σ(s)) も計算する。None なら matrix と同じ

    Returns
    -------
    dict。'_emit' に emit6d.tecalc の dict、'_state' に一周後の _State も入れる (デバッグ用)。
    """
    els = list(lattice.elements)
    circ = sum(el.L for el in els)
    ph = _Phys(momentum_eV, circ, amass=amass, charge=charge)
    ph.dvfs = fshift if rfsw else 0.0          # tsetdvfs (FSHIFT = Orbit dilation)
    mode = 3 if matrix else (2 if orbit else (1 if one_turn else 0))
    if beam_matrices is None:
        beam_matrices = matrix

    # --- tcod 相当 (calc6d._closed_orbit = emit6d.emittance6d と同じ手順) ---
    rfsw0 = _lat.RFSW
    _lat.RFSW = rfsw
    try:
        codin, trf0, found = _closed_orbit(els, ph, rfsw, radcod, True)
        codin = np.array(codin, dtype=float)
        # --- temit の 1 周目: trans=1, beam=0, cod=codin → tecalc ---
        st = _turn_record(els, codin, ph, trf0=trf0, rfsw=rfsw, radcod=radcod,
                          beam0=None, record=(mode >= 2))
    finally:
        _lat.RFSW = rfsw0
    em = tecalc(st, circ, rfsw=rfsw)
    T = st.T.copy()
    D = st.D.copy()
    beamp = st.B.copy()                 # tecalc の beamp = 一周後の beam(1:21) (物理座標)
    cod = st.cod.copy()                 # 一周後の cod (OrbitAtExit)
    r = em['normal_r']
    ri = em['normal_ri']
    emitp = em['beam_matrix']

    # --- tecalc で params に入るが emit6d が返さないもの ---
    rx = T.copy()
    if not rfsw:
        rx[5, :5] = 0.0
        rx[5, 5] = 1.0
    _, ceig = _normal_form(rx)
    cc = np.zeros(6, dtype=complex)
    cc[0::2] = ceig[0::2]
    cc[1::2] = np.conj(cc[0::2])
    cd = np.zeros(6, dtype=complex)
    cd[3:6] = np.log(cc[0::2])
    stab = bool(np.all(np.abs(cd[3:6].real) < 1e-6)) and bool(found)
    rd = ri @ D @ r
    for k, i in enumerate((0, 2, 4)):
        cd[k] = complex((rd[i, i] + rd[i + 1, i + 1]) * 0.5,
                        (rd[i, i + 1] - rd[i + 1, i]) * 0.5) / cc[i]
    if rfsw:
        wrfeff, vceff = _rf_effective(st)
    else:
        wrfeff, vceff = 0.0, 0.0
    heff = wrfeff * CVELOC / ph.omega0 if vceff != 0.0 else 0.0

    nlat = len(els) + 1
    out = {
        'Stable': stab,
        'Region': (1, nlat),
        'Tunes': np.array(em['tunes'], dtype=float),
        'EnergyLossU0': em['energy_loss_eV'],
        'RfVoltageVc': em['rf_voltage_eff'],
        'EquilibriumPosition': em['trf0'],
        'MomentumCompaction': em['momentum_compaction'],
        'OrbitDilation': fshift * circ,             # dleng = FSHIFT·circ (tfgeo.f:134)
        'BucketHeight': em['bucket_height'],
        'HarmonicNumber': heff,
        'OrbitAtExit': cod,
        'DampingRate': np.array(em['damping_rate'], dtype=float),
        'Emittances': np.array([em['emittance_x'], em['emittance_y'], em['emittance_z']]),
        'MomentumSpread': em['momentum_spread'],
        'BunchLength': em['bunch_length'],
        'TuneShiftByRadiation': cd[0:3].imag / (2.0 * math.pi),
        'NominalSpinTune': ph.h0 * GSPIN,
        'PolarizationTime': 0.0,                    # calpol 未移植 (SAD は params=0 で初期化)
        'SpinTune': 0.0,
        'Polarization': 0.0,
        'Polarizations': np.zeros(3),
        'PolarizationVector': np.zeros(3),
        '_emit': em,
        '_state': st,
    }
    if mode >= 1:
        tw = tfetwiss(ri, cod, True)             # tecalc: params(iptwiss:) = tfetwiss(ri, cod, .true.)
        out.update({
            'OrbitAtEntrance': codin.copy(),
            'OneTurnTransferMatrix': T,
            'OneTurnDampingMatrix': D,
            'NormalCoordinates': ri.copy(),
            'ExtendedTwissParameters': {k: float(v) for k, v in zip(TWISS_FUNS, tw)},
            'OneTurnExcitation': beamp,
            'EquilibriumBeamMatrix': np.array(emitp, dtype=float),
            'SpinPrecessionMatrix': np.eye(3),       # srotinit の初期値 (calpol 未移植)
        })
    if mode >= 2:
        out['ClosedOrbit'] = [c for (_, c, _) in st.rec]
    if mode >= 3:
        out['TransferMatrices'] = [t for (t, _, _) in st.rec]
    if beam_matrices:
        # temit.f:214-222: plot .and. postcal → trans=1, cod=codin, beam=emitp でもう一周し、
        # tsetetwiss が各要素入口に beamsize(:,l) = beam を残す (LINE["SIG",l] の値)
        _lat.RFSW = rfsw
        try:
            st2 = _turn_record(els, codin, ph, trf0=trf0, rfsw=rfsw, radcod=radcod,
                               beam0=emitp, record=True)
        finally:
            _lat.RFSW = rfsw0
        out['BeamMatrices'] = [b for (_, _, b) in st2.rec]
        out['_state2'] = st2
    return out


def beam_size_at(e: dict, l: int, key: str = '') -> np.ndarray | float:
    """SAD の LINE["SIG"//key, l] 相当 (tftwiss.f:982-): key='' → 6×6、'X' 等 → sqrt(Σ_ii)、'XY' 等 → Σ_ij。
    l は 1 始まり (要素入口)。小数位置は未対応。"""
    B = e['BeamMatrices'][l - 1]
    idx = {'X': 0, 'PX': 1, 'Y': 2, 'PY': 3, 'Z': 4, 'DP': 5, 'P': 5}
    key = key.upper()
    if key == '':
        return B
    for k1 in sorted(idx, key=len, reverse=True):
        if key.startswith(k1):
            rest = key[len(k1):]
            i = idx[k1]
            if rest == '':
                return math.sqrt(B[i, i])
            j = idx[rest]
            return B[i, j]
    raise KeyError(key)


def rules_to_sad_list(e: dict, mode: int = 3) -> list:
    """SAD の Emittance[] と同じ順の [(名前, 値), …] に並べる (ffscmd/sadscript への接続用)。
    ExtendedTwissParameters は SAD と同じく {AX->v, …} の規則リスト (ここでは [(名前, 値)])。"""
    res = []
    for k in RULE_ORDER:
        if k not in e:
            continue
        if k in ('OrbitAtEntrance', 'OneTurnTransferMatrix') and mode < 1:
            break
        v = e[k]
        if k == 'ExtendedTwissParameters':
            v = [(kk, v[kk]) for kk in TWISS_FUNS]
        elif k == 'Region':
            v = list(v)
        res.append((k, v))
    return res
