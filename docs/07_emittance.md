# 07. Emittance[] / emit — SAD の temit と pysad の emit6d.py

> 最終更新: 2026-11-06 (現状のスナップショット)。行番号は 2026-10 時点。

照合先: `Packages/emit.n`, `src/tfemit.f`, `temit.f`, `tcod.f`, `tturne.f`, `tquade.f`, `tbende.f`, `tthine.f`, `tcave.f`, `tmulte.f`, `photons.f (tradke)`, `tffs.f (module temw)`, `tnorm.f`, `tdrife.f` と `pysad/emit6d.py`, `pysad/ffscmd.py (print_emit)`。

## 1. Emittance[opt] の定義 (Packages/emit.n)

`Override[opt, default]` で選択肢を決め `CalculateEmittance[mode, out, {iniorbit, inibeam}, reg]` (= tfemit → temit) を呼ぶ。既定: `OneTurnInformation->False, Orbit->False, Matrix->False, ExpandElementValues->True, Emittance->True, SaveEMIT->False, InitialOrbit->Null, InitialBeamMatrix->Null, Region->{1,-1}, Output->0`。mode は Matrix→3, Orbit→2, OneTurnInformation→1, Emittance→0, Emittance->False→-1 (固有値のみ)。

戻りは規則のリスト。`e[[2]]` = `params` (temw の ip* 番号)、`e[[3]]` = `setiamat` の 7 要素。

| 規則 | 由来 | 注意 |
|---|---|---|
| Stable | 3 モードの \|Re log λ\|<1e-6 かつ閉軌道あり | |
| Tunes | Im(cd(4:6))/2π | 小数部 (負もあり) |
| EnergyLossU0 / RfVoltageVc | u0·pgev / vceff | **eV** (報告の表示は MV) |
| EquilibriumPosition | trf0 [m] | 同期位置 dz |
| MomentumCompaction, OrbitDilation, BucketHeight, HarmonicNumber | alphap, dleng, bh, heff | |
| OrbitAtExit | cod (一周後) | |
| DampingRate / TuneShiftByRadiation | Re / Im (cd(1:3)) | 一周あたり |
| Emittances, MomentumSpread, BunchLength | eemx..z, sige, sigz | |
| OrbitAtEntrance, OneTurnTransferMatrix, OneTurnDampingMatrix, NormalCoordinates | codin, T=trans(:,1:6), D=trans(:,7:12), **ri** | mode≥1 |
| OneTurnExcitation / EquilibriumBeamMatrix | beamn (正規座標の励起) / emitp (物理座標) を 21 要素パックから 6×6 へ | mode≥1 |
| ExtendedTwissParameters | tfetwiss(ri, cod) | mode≥1 |
| ClosedOrbit | 各要素入口の cod | mode≥2 |
| TransferMatrices | 各要素入口までの**累積**行列 | mode≥3 |
| IntrabeamExcitation | | mode=3 かつ INTRA |

SaveEMIT->False (既定) では計算後に EMITX/EMITY/EMITZ/SIGE/SIGZ を元に戻す。

## 2. 計算の流れ (temit.f)

1. **tcod** (CELL かつ calcod): 6D 閉軌道。RFSW なら先に RFSW を切った 4D 解 (δ 固定) を求め、`cod(5)=asinz((u0·pgev−vcacc)/vceff)/wrfeff−trf0` を初期値に 6D Newton。各反復は `tturne(update=.true.)` で、末尾が `trf0 = asin(U0/Vc)/wrfeff` (trans(5,6)<0 なら π−asin) を決め直す。RADCOD では z のずれを trf0 へ移す (`trf0 += dz/v0`)。重み codw0=(1e-6,1e-7,1e-7,1e-8,1e-5,1e-6)、収束 1e-4。
   **RFSW OFF でも RADCOD ON なら放射込みの閉軌道**を使う (tcod の rfsw=F 枝、出口 δ=−U0/P0、z のずれは trf0 へ。docs/11 A32)。
2. `irad=12`、`tinitr12`、`tsetr0` の後 **tturne** で一周。`trans(6,12)` の 1:6 列 = シンプレクティック部 T、7:12 列 = 放射減衰 D (`tmuld6: D += tr1·(T+D)`)。`beam(1:21)` = 量子励起の 2 次モーメント。
3. **tecalc** で固有値分解と平衡解。4. lfno>0 なら `temitprint` (FFS の EMIT)。

### 2-1. tturne1 の要素ごとの扱い

| 要素 | 内容 |
|---|---|
| DRIFT | `tdrife0(dvon=.true.)`: 放射なし、`z −= L·dv(δ)` |
| BEND (tbende) | `ndiv = 1 (K1=0) / 1+int(\|phic/eps\|)`。enarad なら `nrad=1+int(\|L·crad·(h0 b)²/epsrad\|)`, `1+int(\|phib·h0·anrad\|/epsrad/1e6)` との max (×EMITDIV·EMITDIVB)。K1=0 は tbendebody0、K1≠0 は tbendecorr+tbendebody (tbendib の解析行列)。スライス間と出口エッジ後に tradke |
| QUAD (tquade) | tsolrote → tsetr0 → tqfrie (DISFRIN=0) → tqlfre (FRINGE=1,3; F1) → `tradke(f1in)` か tsetr0 → tsolque×ndiv (`itgetqraddiv`: cod 依存) → 出口側を逆順 |
| SEXT/OCT/… (tthine) | **放射なし** (`enarad=.false.`)。L/6 ドリフト → キック → 2L/3 → キック → L/6。ドリフトは `tdrife0(dvon=.true.)` なので **dv 項 `z −= L·dv(δ)` が付く** (落とすと ESRF (16b) で νs 2.8e-6, α_c 5.6e-6 の相対差、累積行列の (5,6) が六極 1 台ごとに 1.3e-9 ずれる) |
| MULT (tmulte) | キック分割数と itgetqraddiv の max (≤300)、スライスごとに tradke |
| CAVI (tcave) | `ndiv = 1+int(min(\|wL\|, sqrt((v(1/h0+1/h1))²/3/0.01)))`。RFSW: tcavfrie → [ドリフト → 薄いキック]×ndiv → ドリフト → tcavfrie。キック `dh = −va(sin φ + offset1)`, `φ = w t + PHI − trf0·w`, offset1 = sin(trf0 w) (RADCOD OFF) / 0 (ON)。RFSW OFF: 積算のみで長さ L のドリフト |

tturneg が `wrfeff = sqrt|ddvcacc/vcacc|`, `vceff = hypot(vcacc, dvcacc/wrfeff)` を作る。tturneg/tecalc は **RFSW に関係なく**空洞の積算量から Vc, wrfeff を作る (RFSW OFF でも Vc, bh, dz, εz, σz は 0 ではない)。

ROTATE と tchge: BEND は「tchge (入口) → 本体 → tchge (出口)」、薄い多極と L=0 の四極は tthine の先頭・末尾の tchge (pysad `_tbende` → `_tbende_core`、`_tthine`/`_tthine_gen`)。
th≠0 かつずれ≠0 の要素は SAD の tchge の癖 (入口の行列に回転が入らない、docs/11 A25) を `_tchge_jac` が再現するが、この要素を含むリングでは SAD の 6D 閉軌道 (tcod) の z が閉じない/収束が緩く、Emittance[] の δ・チューンは SAD と一致しない (CALC と追跡は一致)。bench_level25f、docs/history/work_rotate.md §3。

### 2-2. 放射 tradke (photons.f)

基準点 (codr0, transr) から現在点までの運動量ベクトルの角度 θ を外積 `xpa` から求め

```
anp = anrad·h1·θ,  uc = cuc·h1³/p0/al1·θ,  dg = −cave·anp·uc (cave=8/(15√3)),  u0 −= dg
dee = anp·uc²·cuu (cuu=11/27)      ← beam の (px,px),(py,py),(δ,·) に加算
```

`crad = r_e(c/amass)²/p0/1.5`, `cuc = 1.5 r_e/rcratio`, `anrad = 5/(2√3)·rcratio` (tphyzp)。`xpa=0` (軸上の四極) なら何も更新しない。RADCOD ON のときだけ cod を更新し `tesetdv` を取り直す。

### 2-3. tecalc

- RFSW OFF なら `rx(6,1:5)=0, rx(6,6)=1`。
- `teigenc` (= teigen.f: tbal の vx, thess の高速 Givens, tqr の二重シフト QR; pysad は `pysad/teigen.py` に直訳) → **tnorm**: 各モードのシンプレクティック積 s を 1 に正規化 (s<0 なら固有値を共役に)、z 成分最大のモードを 5,6 列、x 成分の大きい方を 1,2 列に並べ替え、`r(i,i+1)=0, r(i,i)>0` に回転。`r=tsymp(r)`, `ri=tinv6(r)`。
- `synchm = rfsw ∧ Im(log λz)≠0`。真: `alphap = −Im(cd6)|Im(cd6)|/(c·2π/ω0)/(dvcacc/pgev)`, `ωz = |Im cd6|ω0/2π`。偽: `alphap = −rx(5,6)/2π/c/p0·h0·ω0`, `ωz = sqrt|αp·2π·heff·Vc/pgev·cos φrf|·ω0/2π`。
- `φrf = asin(U0/Vc)`, `bh = sqrt|Vc/(π|αp| heff pgev)·(2cos φrf − (π−2φrf)U0/Vc)|`。
- `rd = ri·D·r`, 2×2 ブロックから `cd = ((rd_ii+rd_jj)/2 + i(rd_ij−rd_ji)/2)/λ`: Re = 減衰率、Im = チューンシフト、J = 4·Re cd/ΣRe cd。
- 平衡: beam を ri で正規座標へ (synchm でなければ z 行・列を 0)、`rdb` (rd + モード回転) から 21×21 の `btr` を組み、(Σ11±Σ22)/√2 の基底変換とスケーリングの後 `tsolva(btr, beam, emitn, 1e-8)` (SVD) で Σ = MΣMᵀ+B を解く。
- `εx = sign(sqrt|e11e22−e12²|, e22·charge)`、`emitp = r·emitn·rᵀ`、`σE = sqrt|emitp(6,6)|`。synchm なら `σz = sqrt|emitp(5,5)|`、εz も行列式。そうでなければ `σz = |αp| σE c p0/h0/ωz`, `εz = σz σE`。
- `rsetgl1` で EMITX/EMITY/EMITZ/SIGE/SIGZ に格納。
- **不安定な光学** (docs/history/work_emit_unstable.md A): 固有値が実数 (λ, 1/λ) でも teigen/tnorm はそのまま通り、`stab=(|Re cd(4:6)|<1e-6) ∧ fndcod` は戻り値 (Emittance[] の `Stable`) にしか効かない
  (INTRA ON のとき `Skip intrabeam because of unstable.`)。Emittance は不安定でも計算・格納され、temitprint は一周行列の偏差 `Σ|rx·tinv6(rx)−1|>1e-8` を ` *** Deviation from symplectic matrix = ` で出す。
  実固有値の対の並び (λ と 1/λ のどちらが先か) は teigen の返す順で、始点が鏡映対称 (a=d) のときは丸め雑音の符号で決まる (pysad `_order_real_pairs` は雑音のとき増大側が先)。

### 2-4. 報告 (temitprint / `FFS.print_emit`)

SAD の EMIT (temitprint) の行: Closed orbit (Entrance/Exit)、Extended Twiss、`P0, f0, U0 [MV], Vc [MV], dz [mm], alpha, dl, h, dV/P0, fs`、Imag./Real tune、Damping (一周・時間・チューンシフト・分配数)、`Emittance X/Y/Z, Energy spread, Bunch Length [mm], Beam tilt, Beam size`。
pysad `FFS.print_emit` はこの全行を同じ書式で出す (docs/history/work_emit_report.md):

- Closed orbit (Entrance と Exit。値は `' '//autofg(v,'9.6')`)、Extended Twiss (EPX の欄は SAD の癖で EY の値、docs/11 A30)。
- P0/f0 (f0 は `c·β/C`、β=p/h)、U0/Vc、dz/alpha、Orbit dilation `dl = -FSHIFT·周長` (RFSW ON のとき `dvfs=FSHIFT` が全要素の dv に入る: `_Phys.dvfs`)、h = `wrfeff c/ω0`、Bucket height、fs = `omegaz/2π`。
- Imag.tune = `Re log λ/2π`、Real tune、Damping (一周あたり Re cd、Damping time = `-2π/ω0/Re cd`、Tune shift = `Im cd/2π`、分配数 `J = Re cd · 4/Σ Re cd`)。
- Emittance X/Y/Z、Energy spread、Bunch Length、Beam tilt、Beam size xi/eta。Beam tilt = atan2(−2σxy, σxx−σyy)/2、xi/eta = 主軸方向の rms (temit.f:1126-1141、`emit6d.beam_tilt_sizes`)。
- Nominal spin tune = `h0·GSPIN`、Polarization time (POL OFF は 0)。spin/偏極 (calpol) の行は POL ON のときだけ出るので未実装 (POL OFF の形は SAD と同じ)。

書式: 実数は `fmt_real` (autofg 直訳)。`vout` は `autofg(x,'10.7')(1:10)` なので 7 桁目が**切り捨てに見える** (dz=20.61866962 mm → `20.618669`、f0=497069.3289 → `497069.32`) — `fmt_real(x,'10.7')`。
`Damping` 系は gfortran の `1p g14.6` (`sadscript.fmt_g(v,14,6)`)、tune は `3(f10.7,10x)` (行末の X 編集は出力されない)。NaN の報告 (閉軌道なし) は幅に右詰めの `NaN`。
ESRF の RFSW × RADCOD 4 通り (`bench_level26a`) の 148 行が SAD と一致 (0 同然の雑音の桁を除く)。

## 3. 内部の座標・変数

- 座標 (x, px/p0, y, py/p0, z, δ)。`tesetdv(δ)`: `dvemit = −(p1+p0)/h1/(p0h1+p1h0)·δ + dvfs`、`dvfs = FSHIFT` (RFSW のときのみ)。ドリフトは `tgetdv` の同じ式で `z −= L·dv`。
- `trf0` = 同期位置 [m]、空洞位相 `phis = trf0·w`。`PHICAV, DTSYNCH, EFFVC, EFFRFFREQ` に書き出される。
- temw の `r` (正規→物理), `ri` (物理→正規) を tecalc が更新。

## 4. pysad emit6d.py の構造

| 部品 | SAD 対応 | 要点 |
|---|---|---|
| `_Phys` | tphyzp | p0, h0, brhoz, crad, cuc, anrad, ω0, `dv(δ)` (dvfs=0) |
| `_State` | trans(6,12), beam, codr0/transr | `T, D, B(6×6), cod, u0, trf0, vc0/hvc0/vcacc/dvcacc/ddvcacc` |
| `_State.piece(f, jac=)` | tmultr5+tmulbs | COD 写像 f の Jacobian を T, D, B に掛ける。**解析ヤコビアン `jac(cod)`** (ドリフト・薄い多極・四極 (tqente/tqfrie/tsolque 本体)・BEND (tbfrie.f/tbendebody0 の trans + dv 項)) があればそれを、なければ**数値 Jacobian** (`lattice._num_jacobian`, 中心差分 eps=1e-7)。z 行 (5,1..4) は SAD と同じ**シンプレクティック条件** `t(5,j) = −t(1,j)t(2,6)+t(2,j)t(1,6)−t(3,j)t(4,6)+t(4,j)t(3,6)` で組み直す (空洞は `symp_z=False`)。`want_jac=False` は cod だけ追う一周、`want_db=False` は SAD の irad=6 (D・B なし) |
| `_jac_drift/_jac_thin/_jac_tqente/_jac_tqfrie/_jac_tsolque_body/_jac_quad_full/_jac_tqlfre`, `_dvz` | tdrife0, tthine, tqente, tqfrie, tsolque, tsolrote, tqlfre, `h0/h1emit**3` | 解析ヤコビアン。`z −= dv(δ)L` の δ 微分は `+h0/h1³·L` (`-dv'` と代数的に一致)。`_USE_ANALYTIC_JAC=False` で数値微分に戻る。ユニットテスト `tests/test_emit_jac.py` |
| `_State.tradke` | tradke | bzh=0 版の直訳 (radcod 分岐込み) |
| `_tdrife0, _tquade, _tbende, 薄い多極, _tcave, _tmulte, _tsteee` | 同名 | 分割数の式は SAD と同一。`_tbendecorr/_tbendebody` は解析行列ごと直訳 |
| `tturne(els, cod0, ph, trf0, rfsw, radcod)` | tturne1 | `_rf_effective` = tturneg、`_trf0_update` = tturne の update (mod は `math.fmod`) |
| `_normal_form` | teigenc+tnorm | numpy eig → 共役対 → z 成分最大を 3 番目、x 成分で 1,2 → tnorm と同じ正規化・回転。RF なし (λz=1) は `r(:,5)=e5, r(:,6)=(η,0,1)` を直接作る |
| `tecalc` | tecalc | btr の組み立て・`tsolva` (= `pysad/tsvd.py` の tsvd、特異値 1e-8) まで直訳 |
| `_tcod_radcod(els, ph, cod, rfsw)` | tcod.f (RADCOD ON) | **tcod の直訳**: RFSW ON は先に RFSW OFF で自分を呼び (`cod(5)=asin((U0 pgev−vcacc)/vceff)/wrfeff−trf0`)、本体は Newton (T−1 を `tsvd` 1e-8)、一周ごとに z の平均 dz を `codi(5)` から引いて `trf0 += dz/v0` (`sign(1,trf0)·mod(…)` で折り返し)、`r=Σ(dcod1/codw)²<1e-4` で収束、`fact` による減速 (fmin 1e-4)、`dpthre`。z と trf0 の最終値は反復の履歴で決まり、SAD と同じ分け方になる |
| `closed_orbit6d(els, ph, rfsw, radcod, close_z)` | tcod | RFSW+RADCOD、RFSW OFF+RADCOD ON (`_tcod_radcod(rfsw=False)`)、および **RFSW ON + RADCOD OFF (空洞あり) も `_tcod_radcod`** (後者は `radcod=False`: 各反復の tturne が trf0 を更新し、収束判定に (dtrf0/trw)² が入り、z を trf0 に移さない。FSHIFT≠0 でも SAD と ≤1.2e-11、bench_level30b)。初期値は CALC の 4D 閉軌道、z=0, dp=0。それ以外 (RFSW OFF + RADCOD OFF) は 4D 閉軌道 (δ 固定) で一周し、z の戻りが 0 になる δ を**割線法**で決める (tcod の直訳ではない。cod だけ追う一周 + `_qcod_accepts_start`)。RFSW: 1 周目 (trf0=0) で U0・Vc・T(5,6) の符号 → trf0。`calc6d._closed_orbit` (Emittance[]/CALC6D) も共用 |
| `emittance6d(lat, p_eV, rfsw, radcod, close_z)` | temit | `closed_orbit6d` → tturne (irad=12; tcod の最後の反復が完全な一周ならそれを使う) → `tecalc`。RFSW OFF でも tecalc は空洞の積算量から Vc・h・bh・ωz・σz・trf0 を出す (SAD と同じ) |

返り値 dict: `emittance_x/y/z` = Emittances、`momentum_spread`、`bunch_length`、`momentum_compaction`、`rf_voltage_eff` = RfVoltageVc、`trf0` = EquilibriumPosition、`bucket_height`、`closed_orbit` = OrbitAtExit (RFSW+RADCOD は tcod の入口の値)、`orbit_entrance`/`orbit_exit` (EMIT の Closed orbit の 2 行、emittance6d のみ)、`energy_loss_eV` = EnergyLossU0、`tunes`、`damping_rate`、`beam_matrix` = EquilibriumBeamMatrix、`one_turn_matrix` = T+D、`normal_r` = r、`normal_ri` = NormalCoordinates。

### 解析ヤコビアンの対象と理由

- 対象: ドリフト、薄い多極 (SEXT〜DODECA、ずれ・ロールは tchge の行列)、四極 (ずれ・ロール付き tsolrote、ハードエッジ tqfrie、線形フリンジ FRINGE=1,2,3 の tqlfre、K1<0 の −π/2 回転、tsolque 本体: `_jac_quad_full`/`_jac_tqlfre`)、BEND (K1=0 の `lattice._tbendebody0_mat` と K1≠0 の `_tbendecorr/_tbendebody`、エッジ)。lattice の CALC 経路と同じ関数を使う (docs/05 §2)。
- まだ数値微分: `MULT`、ANGLE=0 の BEND、F1 付き BEND、DROTATE≠0 の BEND、空洞 (history/work_emit_report.md §5)。
- 理由 1 (桁落ち): K1=0 の BEND 本体の数値ヤコビアンは rhob≈20 倍の桁落ちで相対 1e-8〜1e-7 の雑音が出る (刻みを小さくすると悪化、1e-5 で 2e-10)。この雑音は量子励起の分割数 `int()` を跨ぐ形で効く。解析式にすると雑音が消え、RADCOD ON も SAD と 2e-11 以内になる。
- 理由 2 (分岐境界): 数値微分の刻み 1e-6 は K1=0 BEND の式の分岐境界に当たり (1,6) に 5.6e-7 の誤差が出るため、残る数値微分の刻みは 1e-7。
- 理由 3 (フリンジ四極): tqlfre を含む四極を数値微分にすると ε・チューンに 1e-5〜2.7e-7 の差が出る (`bench_level26b`)。解析式で 1.3e-9。
- z 行のシンプレクティック性: 数値微分の z 行をそのまま使うと T の誤差 6.6e-5、εy が 2% ずれる (Level 3C)。そのため z 行 (5,1..4) はシンプレクティック条件で組み直す (`_State.piece`)。

## 5. 注意点 (仕様と落とし穴)

1. **Cavity PHI の既定は 0**。`phic = PHI·sign(charge)`、同期位相は trf0 で入る。
2. **4D CALC では空洞は RFSW に依らずドリフト** (bench_level14a_dbg)。空洞が効くのは Emittance[] と追跡。
3. **薄い多極の dv 項**: tthine 内のドリフトにも `z −= L·dv(δ)` が付く (2-1 節)。表示桁 (εx 5 桁、Level 12) の照合ではこの差は見えないので、照合は全桁 (`$FORM="18.14"`、16b/16e) で行う。
4. **閉軌道の δ**: SAD の tcod は周長変化で δ≠0 を選ぶ (3C: δ=−4.85e-8)。close_z がないと U0・εx がずれる。
5. **RADCOD ON の z と trf0 (Equilibrium dz)**: tcod は一周ごとに z を trf0 に移し、収束判定が緩い (z は 1e-7 まで) ので、分け方は反復の履歴で決まる。`_tcod_radcod` が tcod を直訳しているので trf0 は ESRF で 12 桁一致する (trf0 固定の厳密な Newton では dz が 8e-4 mm ずれる)。
   逆に言うと cod と trf0 の「和」だけが物理で、SAD の cod 自体は 1e-9 (dp) 程度しか収束していない (出口と入口がその分違う) — それでも emittance 等は 1e-11 で一致する。
6. **RFSW OFF の RADCOD**: SAD の tcod は RFSW OFF でも RADCOD ON なら放射込みの閉軌道を使う (`_tcod_radcod` の RFSW OFF 枝)。RFSW OFF を RADCOD OFF と同一視すると ESRF で U0 が 7e-5、α_c が 1e-4 ずれる。
7. **RFSW OFF の RF の量**: tturneg/tecalc は RFSW に関係なく空洞の積算量から Vc, wrfeff を作る (Vc, bh, dz, εz, σz は RFSW OFF でも 0 ではない)。
8. **解析ヤコビアンの落とし穴**:
   (a) `lattice._tbendebody0_mat` の z 行 (5,6) (`xsn=sn−phi0`) は 4D の CALC では使われないので、z 行の検証は emit 側の照合 (`tests/test_emit_jac.py`) に依る。
   (b) Newton 法のヤコビアンを使い回す (chord) と軌道が 1.7e-14 変わるだけで εx が 2e-8 変わる — 量子励起の分割数 (`nq`, `ndiv` の `int()`) の境界。
   (c) `tecalc` の `tsolva` (= tsvd、特異値 1e-8 の打ち切り) を numpy の SVD で代用すると y 方向に εy~1e-21 の雑音が出る (SAD は ~1e-38) → `pysad/tsvd.py` (SAD の直訳) を使う。
   emit の結果は軌道の連続関数ではない (1e-8 のレベル) ので、SAD との照合では軌道の差 1e-9 の影響を見分けられない場合がある (解析式で SAD と同じ分割になっている)。

## 6. 照合の要約

| Level | 内容 | 一致 |
|---|---|---|
| 2 | 20 セル FODO, combined BEND, RF なし | εx 表示桁まで一致 (≤4e-10), U0 1e-15 |
| 3C/3D | 結合 (ロール) / 強さ誤差 | εx, εy, σE, U0 ≤7.8e-7 / ≤1.6e-9 |
| 4 | RFSW ON (RADCOD OFF/ON) | εx, εz, σE, σz, Q, 減衰率, U0, trf0, α, bh ≤1e-7 / ≤6e-8 |
| 8 | 四極 F1/FRINGE, 空洞 L=0.5, MULT K2 | εx, σE, σz, α, U0 表示桁 (≤2.3e-6) |
| 11 | MULT K0/K1/K2 + DX/DY/ROTATE | trf0 2.6e-10 |
| 28c/28d/28f | 不安定な光学 (V: QD×1.19、H: QF×0.9) の EMIT 報告・Emittance[] | V は行単位で一致 (雑音の桁のみ)、H は SAD の T の悪条件の範囲 (≤3e-5) |
| 30b | RFSW ON + RADCOD OFF + FSHIFT = 1e-6, 3e-7, −2e-6, 0 | ≤1.2e-11 (tcod の直訳) |
| 25a/25b/25c/25d/25h | ROTATE: 垂直ステアリング (ANGLE=0, ROTATE=−90 DEG, K0)、skew 四極、ROTATE 付き BEND (5/90 DEG, K1 付き)、六極・八極・十極・MULT の ROTATE、pf2 の VD/PV/VS と L=0 の SQ | εx/εy/U0 相対 ≤1e-8 (εy 2e-12〜1e-8)、TUN ≤4e-9 (25h は δ の収束差で 3e-8) |
| 25e | ROTATE/DX/DY/DROTATE/K0 の組み合わせ (癖のない組) | εx/εy 相対 1.6e-6 (δ の収束差)、TUN 4.5e-8 |
| 25f | th≠0 かつずれ≠0 の要素 (tchge の癖) | CALC・追跡は一致、Emittance[] の δ・チューンは SAD の tcod が閉じないため不一致 (2-1 節) |
| 12 | ESRF v03_12 (3563 要素) | 表示桁 (相対 ≤2e-5) |
| 16b/16e | 同、全桁 | εx, εz, σE, σz, νx, νy, νs, α_c, bh, trf0, U0 **≤1.3e-11** |
| item6 (`bench_item6_emit.sad`) | ESRF v03_12, RFSW+**RADCOD**, `Emittance[Matrix->True]` の全桁 | εx 1.1e-12, εz 9e-12, σE 2e-14, σz 9e-12, α_c 2e-11, ν ≤9e-12, 減衰率 ≤8e-13, ビーム行列 ≤4e-11, U0・trf0 (dz) は全桁 |
| item6 (`bench_item6_tiny.sad`) | ESRF 1 セルのリング (ν_s=2e-4): RFSW OFF / ON / RADCOD | 横方向・U0・trf0 ≤1.3e-11、縦方向 (RFSW ON) は 1e-9 (νs が小さく 1e4 倍に拡大)。テスト `test_emit_jac.py` |
| 26a/26c/26d | EMIT 報告の全行 (ESRF の RFSW × RADCOD 4 通り、FSHIFT=1e-6)、RFSW OFF + RADCOD ON の Emittance[] と CALC6D | 報告 148 行 + 37 行すべて一致 (雑音の桁を除く)、26d ≤2.2e-11 |
| 26b | FRINGE=1,2,3 付き四極、ずれ・ロール付き四極・六極 (小リング): CALC と Emittance[] | CALC ≤5e-13、Emittance[] ≤1.3e-9 (数値微分では 1.1e-5) |

対応済み (別モジュール): INTRA (touschek.py: tintrb、TouschekTable/TouschekLifetime)、OneTurnInformation/Matrix 系 (emit_extra.py)、EMITXE/EMITYE/EMITZE の下限、CAVI の DX/DY/ROTATE。
INTRA ON の EMIT 報告 (docs/11 A40): 零電流の報告 → Touschek の表 (運動量アクセプタンス 0.4〜24.2 % の寿命、Nx/Nz・Ny/Nz ごとの寿命/100s) → `Parameters with intrabeam scattering:` (PBUNCH/1e10、MINCOUP %、IBS 込みの EMITXC/EMITYC/EMITZC・SIGE・SIGZ)。PBUNCH の既定は 1e10 (A41)。`ffscmd._print_touschek`、`_print_emit_intra`。
未対応: POL (calpol, spnorm, srequpol)、SOL (tsole)、OneTurnInformation の SpinPrecessionMatrix、Region (部分区間)、InitialBeamMatrix、TRPT 分岐、SynchroBetaEmittance (SAD K64 に実体なし)。
