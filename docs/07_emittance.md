# 07. Emittance[] / emit — SAD の temit と pysad の emit6d.py

> 最終更新: 2026-11-18 (v2.1.0 の後) (現状のスナップショット)。行番号は 2026-10 時点。

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
| SEXT/OCT/… (tthine) | **放射なし** (`enarad=.false.`)。L/6 ドリフト → キック → 2L/3 → キック → L/6。ドリフトは `tdrife0(dvon=.true.)` なので **dv 項 `z −= L·dv(δ)` が付く** (落とすと 実ラティス (16b) で νs 2.8e-6, α_c 5.6e-6 の相対差、累積行列の (5,6) が6 極 1 台ごとに 1.3e-9 ずれる) |
| MULT (tmulte) | キック分割数と itgetqraddiv の max (≤300)、スライスごとに tradke |
| CAVI (tcave) | `ndiv = 1+int(min(\|wL\|, sqrt((v(1/h0+1/h1))²/3/0.01)))`。RFSW: tcavfrie → [ドリフト → 薄いキック]×ndiv → ドリフト → tcavfrie。キック `dh = −va(sin φ + offset1)`, `φ = w t + PHI − trf0·w`, offset1 = sin(trf0 w) (RADCOD OFF) / 0 (ON)。RFSW OFF: 積算のみで長さ L のドリフト |

tturneg が `wrfeff = sqrt|ddvcacc/vcacc|`, `vceff = hypot(vcacc, dvcacc/wrfeff)` を作る。tturneg/tecalc は **RFSW に関係なく**空洞の積算量から Vc, wrfeff を作る (RFSW OFF でも Vc, bh, dz, εz, σz は 0 ではない)。

ROTATE と tchge: BEND は「tchge (入口) → 本体 → tchge (出口)」、薄い多極と L=0 の四極は tthine の先頭・末尾の tchge (pysad `_tbende` → `_tbende_core`、`_tthine`/`_tthine_gen`)。
th≠0 かつずれ≠0 の要素は SAD の tchge の癖 (入口の行列に回転が入らない、docs/11 A25) を k64 モードの `_tchge_jac` が再現する (v2 モードは回転を行列に入れる。docs/05)。この要素を含むリングでは SAD の 6D 閉軌道 (tcod) の z が閉じない/収束が緩く、Emittance[] の δ・チューンは SAD と一致しない (CALC と追跡は一致)。bench_level25f、docs/history/work_rotate.md §3。

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
- Nominal spin tune = `h0·GSPIN`、Polarization time (POL OFF は 0)。POL ON のスピンの行・表は 4b 節 (POL OFF の形は SAD と同じ)。

書式: 実数は `fmt_real` (autofg 直訳)。`vout` は `autofg(x,'10.7')(1:10)` なので 7 桁目が**切り捨てに見える** (dz=20.61866962 mm → `20.618669`、f0=497069.3289 → `497069.32`) — `fmt_real(x,'10.7')`。
`Damping` 系は gfortran の `1p g14.6` (`sadscript.fmt_g(v,14,6)`)、tune は `3(f10.7,10x)` (行末の X 編集は出力されない)。NaN の報告 (閉軌道なし) は幅に右詰めの `NaN`。
実ラティスの RFSW × RADCOD 4 通り (`bench_level26a`) の 148 行が SAD と一致 (0 同然の雑音の桁を除く)。

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
- まだ数値微分: `MULT`、ANGLE=0 の BEND、DROTATE≠0 の BEND、空洞 (history/work_emit_report.md §5)。F1 付き BEND の本体は v2.1.0 で `emit6d._tbende_core` (tblfre・tbendal・断片のループ) として実装した (微分の方式は未確認。F1 の無い BEND は従来と同じ演算)。
- 理由 1 (桁落ち): K1=0 の BEND 本体の数値ヤコビアンは rhob≈20 倍の桁落ちで相対 1e-8〜1e-7 の雑音が出る (刻みを小さくすると悪化、1e-5 で 2e-10)。この雑音は量子励起の分割数 `int()` を跨ぐ形で効く。解析式にすると雑音が消え、RADCOD ON も SAD と 2e-11 以内になる。
- 理由 2 (分岐境界): 数値微分の刻み 1e-6 は K1=0 BEND の式の分岐境界に当たり (1,6) に 5.6e-7 の誤差が出るため、残る数値微分の刻みは 1e-7。
- 理由 3 (フリンジ四極): tqlfre を含む四極を数値微分にすると ε・チューンに 1e-5〜2.7e-7 の差が出る (`bench_level26b`)。解析式で 1.3e-9。
- z 行のシンプレクティック性: 数値微分の z 行をそのまま使うと T の誤差 6.6e-5、εy が 2% ずれる (Level 3C)。そのため z 行 (5,1..4) はシンプレクティック条件で組み直す (`_State.piece`)。

## 4b. 追加した要素・機能の範囲

### POL のスピン (docs/36_spin_csr.md)

`POL;` を付けた `Emittance[]`/`EMIT` は、偏極時間 (`PolarizationTime`)・スピンチューン・偏極ベクトル・スピン歳差行列・1 周の脱偏極ベクトル (表) を出す (`pysad/spin.py` の `tradke_spin`、`spnorm`、`srequpol`)。
放射の 1 ステップ (`tradke`) ごとに `srot` を回して積算し、`tchge` の回転 (ROTATE 付き偏向) でも回す。POL のとき BEND の分割数は `1+int(|φ|·h0·a/0.01)` 以上になる (K64 も同じ)。
- 一致 (gs): 平面リングで偏極時間は全桁、スピンチューン 13 桁、偏極ベクトル 12 桁。縦ずれ・傾き・6 極のあるリングで 13/12/10 桁以上。POL 付きの Emittances・DampingRate は 11〜12 桁。エネルギー 3.0〜3.08 GeV の格子 (aγ=6.8〜6.99) で確認した範囲。
- **MULT**: K1・K2・K3 とずれのある MULT のリングが gs と 1e-11 (`test_mult_ring_matches_gs`、L=0 の薄い MULT は `test_thin_mult_ring_spin_tune_matches_gs`)。`tmulke` の bsir0 の更新も入れた (この例では効かない)。
  K3 が強くずれの大きい MULT リング (K3=300、ずれ 1e-4) は pysad の CALC (qcod) の閉軌道が gs と別の解になり Emittance が NaN (未解決、docs/05 4b)。
- **SOL** (BOUND≠0 の領域): `emit_extra._tturne1` が SOL 領域 (`_solmap`) を通し、偏極時間・スピンチューン・偏極ベクトルが gs と一致 (`test_solenoid_ring_matches_gs`)。領域内の MULT・BEND・CAVI のスピンは未確認 (領域内の要素そのものが未対応)。
- **空洞**: ROTATE+DX+DY 付きの空洞を含む結合リングのスピンは、2026-11-16 に gs と比べて確認済み (偏極時間 13 桁、スピンチューン 1e-14、偏極ベクトル 1e-13)。ただし専用の pytest はまだ無い (gs の出力との手動比較)。
- **MULT ANGLE≠0 (`tmultae`)** の POL スピン: 2026-11-16 に確認済み (偏極時間 14 桁、スピンチューン 1e-14。偏極ベクトルは平面リングなので 0,1,0 で、傾きのあるリングでの偏極ベクトルは見ていない)。**`tbrote`** (BEND の DROTATE+DX): 確認済み (偏極時間 14 桁、スピンチューン 1e-14、偏極ベクトル 1e-14)。
- **Polarization2**: gs の `Emittance[]` は Polarization2 の値を返さず (`Polarization2/.em` が未評価のまま)、pysad も同じ (2026-11-16 確認)。
- 未対応・未確認: `Polarization` (平衡偏極度) は K64 が NaN (B55)。**k64 モードは NaN のまま、v2 モードは有限値** (docs/43、下の「v2 モードで変わった放射・エミッタンス」)。`MAP` のスピン、SOL 内のスピン (`tsconv`) は未対応。
  縦オフセットのある偏向磁石の `bsi2` の項は K64 では捨てられる (B58。v2 モードは入れる。下記)。POL ON は numba の速い経路を使えない (Python 経路)。v2 モードの POL では四極・MULT を細かく分割するので Emittance[] が遅くなる。
- **TrackParticles のスピン追跡**は docs/08 §9 (POL かつ RAD ON のときだけ、`track6.py` の `sprot`、損失の処理なし)。

### SynchroBetaEmittance (`pysad/synchrobeta.py`)

K64 の `tfemits`・`temits` ほか (`tevdif`・`tesolvd`・`tesumb` など) の直訳。`Packages/emit.n` の `SynchroBetaEmittance` が `SynchroBetaEmittance1[nus0, nus1, dnus, mphi]` (`ffscmd.fn_synchrobeta`) を呼ぶ。
運動量のずれ (2·ndims+1 点) ごとに `_tcod_radcod` で閉軌道、`tturne` で T・D・B を求め、(アクション J, 方位角モード m) の連立方程式を反復で解く。
- gs との一致 (12 セルのリング + CAVI VOLT=2e6 HARM=64): nus=0.01・AzimuthalModes=9 で emitx 相対 6e-6 (conv は 2e-4)。nus=0.005〜0.02・AzimuthalModes=5 の 4 点で emitx 4e-7〜1e-5。emity は 1e-18 の雑音で符号も違う。試験は `test_synchrobeta.py` (許容 3e-5)。
- 差の原因は閉軌道の収束の履歴と、条件数 3e6 の連立方程式の解き方 (K64 の tsvd と numpy の SVD) とみているが、確認していない。
- 2026-11-16 の確認: 24 セルの平面リング + CAVI VOLT=2e6、AzimuthalModes->7 の 4 点で emitx は gs と相対 1.4e-7〜1.4e-6 (3.7e-7、4.0e-7、1.4e-6、1.4e-7) で一致、emity は雑音 (gs −2.5e-19 など、pysad 1.3e-20 など)。**結合のあるラティス** (QT・SF・DY のずれ入り 8+4 セルのリング) は gs が SIGSEGV で落ちて比較できず、pysad も `ValueError: cannot convert float NaN to integer` で失敗する (エラーの出し方が粗い)。**v2.1.0 の確認 (C)**: 24 セルの平面 FODO + skew 四極 1 台 (3 GeV、RF 4 MV、6 極なし) の結合のあるリングでは gs は落ちず、k64 の pysad は gs と 7 桁以上一致した (εx 1.3603366e-6 対 1.3603365e-6、εy 4.7516654e-8 対 4.7516654e-8)。Emittance[] と、独立な 6x6 の平衡エミッタンス (自前の 1 周写像 + 放射の減衰 + 量子励起の拡散 + Lyapunov 方程式) とも一致 (ASTRID2 + skew 四極 ks = 0.05 で ε₁ 6e-5、ε₂ 1.3e-4)。SB と Emittance[] の差 (ASTRID2 +0.8%、PF +1.1%・+2.9%) は結合の有無に依らず、シンクロ・ベータ結合と色収差の分。SOL・強い結合で gs が落ちるかは未確認。
- **大きなリング A (3 GeV、周長 603 m、νs = 0.0041、RF 1.4 MV) では使えない** (B128): gs も pysad k64 も 6.18e-9 (AzimuthalModes ≥ 3 で変わらず) を返すが、Emittance[] の 1.55e-10、独立な 6x6 の平衡 1.5173e-10、RAD + FLUC の 6 次元追跡 (60 粒子、4000 周) の 1.6e-10 と 40 倍違う。AzimuthalModes 12・20・30・32 以上で 3.19e-9・1.46e-9・8.50e-10・7.87e-10 と動き続け、νs の入力・エネルギー・6 極の強さにも依らない。原因は未特定 (K64 の式の限界か、gs・pysad が共通に持つ誤りか)。PF・ASTRID2・SP8 (大きなリング B、8 GeV) は Emittance[] と 0.5〜1%。**v2 は AzimuthalModes m と 1.5 倍の m の emitx が 20% 以上違うと警告を出す** (値は変えない。大きなリング A の m = 9 で 5.1e-9 対 2.5e-9)。使うときは AzimuthalModes を変えて収束を確認し、Emittance[] と比べる。`AzimuthalModes->1` は pysad が IndexError (gs は 8.6e-7 を返すが意味がない)、2 以上を使う。
- 未確認: 不安定なチューン、`emiout`、SOL・強い結合のラティスでの gs との一致、RFSW の組み合わせ。検証は最終の 6 つの値 (nus, emitx, emity, emitxp, emityp, conv) だけで、中間配列は比べていない。

### SOL・COORD・TCAVI・UND・MAP・MULT (ANGLE≠0) と Emittance[]

| 要素 | Emittance[]・EMIT の扱い | gs との一致 | 文書 |
|---|---|---|---|
| SOL (BOUND≠0 の領域) | `solemit.py` (領域内を K64 の `tsole1` の順に、写像は追跡カーネルの数値微分、放射は bz 付きの `_State.tradke(al, 0, bzh)`) | ε_x 2e-12〜4e-9、ε_y 4e-10〜2e-8、ε_z 1e-7〜3e-6、チューン 1e-13〜1e-11 (BZ=2 と BZ=5 の 12 セルのリング、RFSW) | docs/29 |
| COORD | `emit6d._tcoord` (`piece`、z 行のシンプレクティック組み直しなし) | ε・チューン・COD・運動量収縮因子 5e-10 以下 | docs/25 |
| TCAVI | `emit6d._ttcave` (tchge → 半ドリフト → キック → 半ドリフト [tradke] → tchge)。RFSW OFF はドリフト | **L=0 のみ gs と一致** (ε 1e-11)。L≠0 は gs が壊れる (B36: EnergyLossU0 が NaN など) ので pysad は NaN を再現せず自然な値を出す | docs/26 |
| UND | k64 モード: 恒等 (K64 の `tturne` は `case default`)。K64 は磁場を設定できず磁場なしのドリフト (B37)。v2 モード: 磁場がなくても長さ L のドリフト、`FBX`・`FBY` つきは追跡のカーネルの数値微分と 1 周期 64 分割以上の `tradke` | k64: 1e-12 台。v2: gs は磁場を設定できないので比較なし (解析式で確認) | docs/27 |
| MAP | `emit6d._tmap` (`ExternalMap["EMIT", l, cod]` の `{cod1, trans[, dtrans, dbeam]}`)。定義がなければ恒等・長さ 0 | ε・チューン・減衰率・軌道 1e-11 | docs/30 |
| MULT (ANGLE≠0) | `emit6d._tmulte` → `multae.emit` | ε 1e-11、チューン 1e-13 (12 セルのリング) | docs/28 |

- **SOL の長さを周長に入れない**: K64 は SOL (BOUND≠0) の長さを周長 (omega0・MomentumCompaction・偏極時間・空洞の位相) に入れない。`emit6d._emit_circ` が同じにしている (PLANAR + SOL で周長 45.6 のときと偏極時間が 15 桁一致)。
  追跡の `circ` も同じ。GEO のない SOL では K64 の S 座標が途切れて周長がさらに狂う (偏極時間が約半分) が、pysad は再現しない (GEO=1 の参照だけ照合)。
- SOL を含む Emittance[] は `emit_extra.emittance_rules` (Orbit/Matrix の規則、BeamMatrices) を使わず `run_emit` に落ちる。EMT・TUN・ALP・COD は合うが、**Stable は 0、HarmonicNumber は NaN、OneTurnTransferMatrix などの規則は出ない** (docs/29)。
- SOL の限界: CALC6D・RADCOD つきの CALC・Touschek (INTRA) との組み合わせは未対応、RADCOD の閉軌道探索中に K64 が内部分割の間で固定するパラメータを pysad は毎スライスの δ で計算し直す (影響は 2 次とみているが未確認)、
  領域内の MULT・BEND・CAVI は k64 モードでは未対応 (v2 モードでは solemit が扱う。docs/44)、行列は数値微分なので解析行列の gs と 1e-10 程度の差。
- v2 モードの SOL 領域: 領域内の BEND・MULT・CAVI (RF つき MULT を含む)・F1 つきの端の放射・ずれて傾いたソレノイド・TRPT の加速・AUTOPHI が Emittance[]・EMIT に入る (docs/29・44)。**POL の量 (SpinTune・Polarization など)**: k64 は領域の中が K64 の式 (縦磁場を z 軸の向きに足す: B130) のままで、ずれて傾いたソレノイドではスピンチューンが追跡と 1.5e-4〜1.1e-3、n0 の向きが 2.5e-3〜1.5e-2 rad 違う。**v2 は v2.1.0 で修正** (B⊥ と B∥ を運動量の向きに分けて足す。`spin.tradke_spin`・`solemit._trk`): 追跡との差は SpinTune で 1e-13〜5e-7、n0 で 0〜1e-5 rad (小さな FODO・大きなリング A・PF・ASTRID2・MPW 入りで確認)。回転した縦磁場の歳差の閉じた式とも 5e-7。条件 1〜3 T・0.3〜3 GeV。bz ≥ 3 1/m の QUAD・MULT の D は 1e-4 台 (docs/44)。
- MULT (ANGLE≠0) の限界: FB1/FB2・CHI2 は Twiss・Emittance[] で未対応、DROTATE≠0 は ε で 5e-6 ずれる、K0 + DX DY + ROTATE の同時指定は gs と閉軌道が一致しない (docs/28)。

### v2 モードで変わった放射・エミッタンス (RADINT・Jx・MomentumCompaction・EMITDIV ほか)

k64 モードは K64 (gs) と同じ。v2 モード (既定) は下の表の「v2」の動作。試験 `test_v2_A_radint.py`・`test_v2_A_damping.py`・`test_v2_A_spin.py`・`test_v2_A_tchge.py`・`test_v2_A_track.py`・`test_v2_A_ffs.py` (計 67 件。詳細は `docs/history/v2_changes_A.md`)。根拠は解析式・独立な数値計算・保存則で、gs は根拠にしていない。条件は 3 GeV の FODO リング (24 セル 105.6 m / 91.2 m、12 セル 45.6 m / 62.4 m)、偏向磁石 L=1.2 m・15 度 (ρ=4.58 m)、RF 500 MHz。**大きなリング・挿入光源・強い結合では確かめていない。**

**(a) RADINT (INTGRL.f)**

| ID | K64 の動作 | v2 の動作 | 根拠 |
|---|---|---|---|
| B06 | 全 BEND に最後の BEND の L/ANGLE を使う | その BEND 自身の L/ANGLE | 曲率の違う 2 種の BEND で自前の 3x3 行列の積分と I1〜I5 が 2e-5 (K64 は 1e-3 以上ずれる) |
| B07 | 6 極が直前の DL を使う。L=0 は 0.2 m。先頭が SEXT だと DL=0 で NaN。DKL2 を代入しない | 6 極は自分の L と K2。L=0 は寄与なし。モジュール変数を毎回初期化 | 軌道 0 なら6 極は I1〜I5 を変えない (K2=3 と 0 で 12 桁同じ)。NaN にならない |
| B12 | 定数が単精度 (1e27、1.E-6、0.2、2.3026) | 倍精度 | ELOSS = Cγ E⁴ I2/(2π) と 1e-13 (K64 は相対 1.15e-8 小さい) |
| B13 | CAVI は最後の空洞の値 | VOLT は総和。FREQ・HARM は電圧のある最初の空洞 (違う周波数が混ざると注意) | 2.5 MV + 1.5 MV → 4000 kV。νs と 2e-3 |
| B72 | `QUANTUM LIFETIME(hours)` が exp(2.3026·τ[秒]) | τ[秒]/3600 | τq = τE/2·exp(ξ)/ξ |
| B73 | `LINEAR CHROMATICITY` が −I9/4π (磁石の中だけの ∮γ ds で色収差にならない。例: −2.36、実際は −7.30) | 運動量を ±1e-4 ずらした光学のチューンの差分 dν/dδ (6 極・偏向の寄与を含む) | `DP0` をずらした CALC の差分と 2e-5 |
| B74 | MULT・SOL を黙って無視する | 場のある MULT・SOL があれば「入っていない」と表示する (積分には入れない) | — |

(番号は docs/21 の B 番号。B72〜B77 と B106〜B110 は v2 の作業で採番した K64 のバグ。)

**(b) Jx (減衰分配) の差と MomentumCompaction** (旧 docs/46 の「RADINT と EMIT の Jx の差」「BUCKET HIGHT の差」の結論)

- **物理として正しいのは RADINT の Jx = 1 − I4/I2。** K64 の Emittance[]・EMIT と追跡は、偏向磁石の端の角 e の蹴りを偏向角に入れるが場の中の行路 al1 を変えないので、放射 ∝ θ²/al1 の I4 の端の項 −η tan(e)/ρ² が 2 倍になる (矩形磁石のリングで Jx = 1.097、正しくは 1.006)。v2 は端の角のある面で場の中の行路を x·tan(e) だけ短くする (B75)。
- 確認 (自前の I4/I2 と、硬い端の厳密な幾何で運動量のずれた閉軌道を求めた Je = 2 + dln l/dδ が 2e-6 で一致): E1=E2=0 で 0.920824、0.1309 で 0.941466、0.5 (矩形) で 1.005751。RADINT はこれと一致。Emittance[] は k64 が 0.92234 / 0.9648 / 1.09735、v2 が 0.92234 / 0.94263 / 1.00586 (`EMITDIV=100` で 0.94154 / 1.00580)。追跡の一周写像の固有値 (E1=E2=0.5) は k64 1.0964、v2 1.0033。gs (2026-10-19 ビルド) の Emittance[] の DampingRate は Jx = 1.0974 (pysad の k64 と 12 桁一致)、同じ gs の RADINT は 1.00575100: K64 の中で EMIT と RADINT が食い違っており、正しいのは RADINT。
- 影響: K64 の Emittance[] は、端の角のある偏向磁石のリングで εx を (矩形磁石の例で) 8% 小さく、σδ を 2% 大きく出す。v2 は RADINT と合う。
- 既定の分割 (EMITDIV=1) では v2 の Emittance[] の Jx が解析値より 1.5e-3 大きい (放射の刻みが粗い: K64 の精度)。`EMITDIV=20;` で 1e-4。Jx+Jy+Je は 3.9990 (4 との差 2.5e-4 = U0/2E の 2 次の効果)。
- B76: 追跡の tbendi (K1≠0 の BEND) は放射で運動量が変わっても曲率半径 rhoe を入口の値のままにするので、減衰の和 Jx+Jy+Je が 4 にならない (K1=0.1 で 4.05、Jx=0.7385)。v2 は rhoe = rhob·(1+δ) を更新し、和が 4.000、Jx = 0.6875 (解析 0.68762、RADINT 0.68762、Emittance[] 0.6883)。
- **MomentumCompaction (B77)**: K64 の RFSW OFF は −R56/C だけで入口の分散の項 Σ R5j η_j が抜ける (FODO の例で 0.022898、正しくは 0.022137: 3.4% 大きい)。RFSW ON はシンクロトロンチューンから μ² で求めるので μ²/12 だけ大きい (0.25%)。v2 は RFSW OFF が −(R56 + Σ R5j η_j)/C、RFSW ON が 4 sin²(μ/2)。自前の I1/C と 2e-4 (FODO の例は 6 桁: 0.0221368)。RFSW OFF の BucketHeight・BunchLength・νs (√α に比例) も直る: BucketHeight 0.009877 → 0.010045 (RADINT 0.010045、RFSW ON 0.010047。残りの 2e-4 は空洞と放射の分)。旧「`BUCKET HIGHT` の 0.010045 と 0.010035」の差の原因はこの α の差で、解決。
- **EMITDIV**: 大域変数 `EMITDIV`・`EMITDIVB`・`EMITDIVQ` (放射の分割の倍率) が Emittance[] に効く。**両モードで効く** (以前の pysad は無視していた。MAIN 部の `EMITDIV = 20;` も効く)。gs: `EMITDIV=20` で DampingRate −.000286831320204 (pysad k64 と 12 桁一致)。

**(c) 追跡・スピン・そのほか**

| ID | K64 の動作 | v2 の動作 | 根拠 |
|---|---|---|---|
| B14 | tsolqum (放射なしの MULT) の z の dv の項が −0.25·(ap−dv)·aln0 | tsolqu と同じ −dv·aln0 | 同じ K1・L の QUAD と MULT の z が 1e-12 (K64 は 1.25·dv·L ずれる。20 MeV・L=0.5 m で 1.9e-5 m) |
| B08 | tsolqur (RAD、SOL 領域内の四極) が最後の半ステップで最後の粒子の bzp を全粒子に使う | 粒子ごとの BZ/(1+δ) | 1 粒子の結果が一緒に追う粒子に依らない (K64 は 2.6e-6) |
| B39 | tradke (SOL 内) の dbeam(8) が B12 − B44 + b·B13 | B12 − B34 + b·B13 | 力学的運動量 (px + b y, py − b x) の共分散 |
| B58 | BEND の端の縦磁場の積分 bsir0: 出口の項が常に捨てられる、K1 の項は積分値 akx を掛ける (次元が合わない)、スピン軌道結合に端の項の y・x 依存が入らない | 入口で +(y/ρ + k x y)、出口で −(y/ρ + k x y)。軌道についての微分もスピン軌道結合に入れる。追跡の tbendi も同じ | 1 要素のスピン回転が Emittance[] と追跡で一致: K1=0.3 の BEND 1.7e-6 (K64 2.0e-3)、K1=0 の BEND 7.5e-5 (K64 1.0e-3) |
| B64 | tmulke の bsir0 が imag(cx)/aln (次元 1/m)。四極・MULT の分割は放射の大きさだけで決まる | ψ = Im Σ k_n z^(n+1)/(n+1)! を入口で足し出口で引く。POL のときは四極・MULT を k·(L/n)² ≤ 2e-3 に分割 | MULT (K2=30) 3.1e-7 (K64 9.4e-5)。四極の分割が 1 個だと歳差の微分が追跡と 50% 違い、細かくすると 2% |
| B55 | `Polarization` が NaN | 有限値 (docs/43) | 平面リング 92.376% (Sokolov–Ternov) |
| B56 | tbrote の srot の更新が行の逐次更新 | (I+rr)·srot (直交性 1e-14) | 入口と出口が互いに逆 |
| B111 | tradke (Emittance[]・閉軌道) の放射の反作用の角度の項の符号が tradkn (追跡) と逆。RADCOD の閉軌道が追跡で閉じない (小リングで x 4e-7 m) | 符号を反転 (`ddpx`・`ddpy` と行列 D・B の微分)。RADTAPER の追跡の BEND の重みも tturne と同じ (3δ_l + δ_l+1)/4 | 閉軌道の 1 周の追跡の動きが 4e-7 → 1e-9 m (小リング)、大きなリング A −2.0e-8 → −4.5e-11 m。**v2 の Emittance[] (RADCOD) は gs と 1e-4〜1e-3 ずれる** (ε_x +1.9e-4、U0 +1e-7。大きなリング A)。符号は導出と閉軌道の閉じ方による (K64 の作者の意図の裏づけは無い) |
| B112 | 放射の古典半径が質量によらない (`re0 = elradi`)。陽子 (30 GeV/c) の U0 が 1.2651 eV | `re0 = elradi·elmass/amass` | 古典式 U0 = (2/3) r_p mc² β³γ⁴ Σθ²/L = 6.88e-4 eV と追跡の 1 周の δ が 1e-6。電子は両モードで同じ |
| B114 | ずれた BEND の CALC (tbrote の出口) の dx1・dz1 の符号 | 出口も入口と同じ回し方 | 厳密な幾何と CALC が 2e-8。矩形磁石 1 台の DX = 2 mm でチューンが動かない (K64 は 2.7e-5) |
| B115 | RADINT が ANGLE=0・K0≠0 の BEND (挿入光源の極) を無視 | K0 を曲げ角として積分に入れる (端の角なし) | ASTRID2 + MPW で U0 6.17 keV → 7.81 keV (Emittance[] 7.80 keV)。εx は 3.5%、Jx は 0.7% Emittance[] とずれる |
| B117 | K1≠0 の BEND の追跡 (tbendi) の端が薄い写像 | `DISFRIN=0` かつ端の角 ≠ 0 のとき CALC と同じ厳密な tbedge | 厳密な写像と 1e-10。CALC との行列の差 1.1e-5 → 4e-8 |
| B107・B126 | RADINT の BMAG の CKL1Y = −cos θ·DKL1。QUADRU が四極自身の DX・DY を引かず K1<0 で勾配ベクトルの符号が逆 | CKL1Y = +sin θ·DKL1。(dx, dy) = 閉軌道 − 要素のずれ | I4U + I4V が物理的な I4 と 1e-7 (K64 は 1.6〜2.7%)。ΔJx = +0.0078 (追跡 +0.00781、K64 +0.0028) |
| B106・B109 | tradke の b = bzh·0.5。tbgrote の trans(5,6) | b = bzh。コメントアウトされた式 | B106 は数値の効果が 1e-13 以下 (式の一貫性)、B109 は有限差分と 4e-12 (K64 は px² + py² 分) |
| B108 | ANGLE=0 の BEND の端の角 AE1・AE2 の放射の項が約 2 倍 (B75 と同じ) | 行路を x tan(e) だけ短くする (`tin`・`tout`) | 厳密な幾何 + ロビンソンの定理で ΔJx が厳密 : v2 : k64 = 1 : 1.03 : 2.08 (傾き 5 mrad の極。24 セル FODO の極 3 個の組) |
| B127 | WSPAC の最後のスライスで `gammab(nlat+1) = 0` を読み、レンズが 8 倍強い (Emittance[] のチューンシフトが 3〜4% 強い) | p1 = (1+δ) p0 | 一次のラスレットの積分と 0.2〜1.8% (大きなリング A・PF・ASTRID2。K64 は 0.9〜7%) |
| B130 | Emittance[] の tradke が縦磁場の積分を z 軸の向きに足す (傾いたソレノイド) | 運動量の向きに分けて足す | 上の POL の節 |
| B132 | PhotonList のエネルギー列に γ でなく βγ を掛ける | γ を掛ける | gs と E 列だけ相対 1/(2γ²) (3 GeV で 1.45e-8) |
| B63・B36・B21 | (バグではない / pysad は有限 / z を閉じる δ) | 変更なし | B63: SOL に L のキーワードは無い (L を書いても軌道・光学に入らない)。pysad は S 座標 (LENG) にだけ入れていたのを直した。B36: L≠0 の TCAVI のリングの Emittance[] は gs が NaN、pysad は有限 (K0=0 で同じ長さのドリフトと 12 桁一致)。B21: RFSW OFF・RADCOD OFF の閉軌道は z が閉じる δ (1e-15) |

- 追跡の FLUC の乱数列は K64 と同じだが、v2 では端の角のある BEND の放射の大きさが変わるので値は K64 と違う。
- SOL (BOUND≠0) の長さを周長に入れないのは K64 の仕様 (B63。バグではない)。

**使える範囲と条件・使えない範囲 (v2)**

| 機能 | 使える範囲 | 使えない・注意 |
|---|---|---|
| RADINT | BEND (K1 つきも)・QUAD・SEXT のリングと輸送路 (INS)。先頭が SEXT でもよい。空洞が複数あってもよい (同じ周波数)。GX・GY は一次の色収差 dν/dδ | MULT・SOL は積分に入らない (表示が出る)。ANGLE≠0 の BEND の K0 は入らない (ANGLE=0 の K0 の極は v2 で入る: B115。εx は Emittance[] と 3.5%)。結合のあるリングの V 面 (I4V) は v2.1.0 で確認 (I4U + I4V が独立な I4 と 1e-7、DU・DV・Je が 1 次の摂動論と 3e-3。K1≠0 の BEND は B107 の修正後)。軌道のずれのある四極は入る (B126)。SEXT・BEND のずれは入らない。空洞は VOLT を単純に足す (PHI を見ない。位相の違う空洞のベクトル和にするかは仕様の判断)。I1〜I5 は独立な積分 (大きなリング A・B・PF・ASTRID2) と 1e-6 以下 |
| Emittance[]・EMIT | 端の角のある偏向磁石 (E1・E2)、勾配つき偏向磁石。Jx は解析値と 1e-3 (既定)、`EMITDIV=20` で 1e-4 | ANGLE=0 の BEND (補正磁石・挿入光源の極) の AE1・AE2 の端の角は K64 が端の項を約 2 倍にする (B108。v2 は直した)。K1≠0・K0≠0 の ANGLE=0 BEND と MULT の端の角 (AE) は pysad が読まず未確認。既定の分割 (EMITDIV=1) では低エミッタンスのリング (大きなリング A) で εx が 2.1% 大きい (`EMITDIV=5` で 8e-5)。v2 の Emittance[] と RADINT が矩形磁石で 1.4e-3 違う (SP8 (大きなリング B)。未調査) |
| MomentumCompaction・BucketHeight・BunchLength (RFSW OFF) | I1/C と 2e-4 | — |
| CALC のチューン・結合 (ロール + ずれのある BEND・6 極・薄い四極・空洞) | 追跡の一周写像と一致 (6 極 1e-9、BEND 4e-6) | ずれた BEND の CALC は v2 で追跡と一致 (B114。矩形磁石 1 台の DX = 2 mm でチューンが動かない。k64 は 2.7e-5 動く)。K1≠0 の BEND は B117 の修正で追跡と 4e-8。CALC の解析行列は y ≠ 0 の軌道で自身の写像の数値微分と y = 2 mm で 1.5e-6 違う |
| Polarization (POL) | docs/43 | docs/43 |

**要判断・直していない**: K1=0 の BEND の追跡のスピン (tbend + tradk) が Emittance[]・tbendi (v2) と 2.2e-5 rad 違う (放射の分割が粗く、スライスの前のエネルギーで歳差を計算する。K64 の精度の範囲)、tthine (長さのある6 極など) の bsir0 (B110。次元が合わない式だが放射なしで呼ばれ入口と出口で打ち消すので影響なし)、RADCOD の 1 周の解析 T と写像の数値微分の 2e-3 の差 (K64 の tradke の微分が近似)、RADINT の SEXT・BEND のずれ、K0 の極の εx (3.5%)、tradke の符号 (B111) と RADTAPER の重みの選択 (「閉じる」を根拠にした判断で、gs の Emittance[] (RADCOD) と 1e-4〜1e-3 違う)。詳細は docs/47 の v2.1.0 の (d) と docs/46。

### v2.1.0 の追加 (F1・FRINGE・K0 の極・WSPAC・SynchroBeta の警告)

試験 `test_v2_B1_optics.py`・`test_v2_E_steer.py`・`test_v2_D1_emit.py`・`test_v2_A2_radint.py`・`test_v2_C3_radint.py` ほか。記録 `history/v2next_{A1,A2,B1,C,D1,E}.md`。

| 機能 | モード | 使い方・使える範囲と条件 | 使えない・未確認 |
|---|---|---|---|
| **BEND の F1・FRINGE・FB1・FB2** の Emittance[]・EMIT・CALC6D | 両 (K64 の動作の追加) | `F1=0.1 FRINGE=3` など (FRINGE = 1・2・3・−1 (入口側だけ)・−2 (出口側だけ)、F1 + FB1 + FB2、K1≠0 (tbendebody の分割)、K0、RADCOD、DX・DY・ROTATE・DROTATE、POL ON の減衰率)。K64 の `tbende` の長さ補正 `alcorr`・入口の `tblfre`・断片 (K1=0 は n = 0 と ndiv+1 に `tbendef1`、K1≠0 は n = −1, 0 と ndiv+1, ndiv+2 の `tbendal`) の直訳。**周長に F1 の長さ補正が入る** (ω0・運動量圧縮・周回周波数)。k64 は gs と全桁一致: 12 セル FODO (3 GeV、F1 = 0.1 FRINGE = 3) で U0 739673.742271 eV (F1 なし 781633)、Tunes {−.2933271700、−.1395704757、−.0179539450}。大きなリング A (全 BEND に `F1=0.03 FRINGE=3`) の U0 285476.566791 eV (F1 なし 296074.86)・運動量圧縮・減衰率・εx が 1e-7、PF・ASTRID2・MPW 入りも 1e-9。放射・RF なしの Emittance[] のチューン = CALC のチューン (1e-9)。F1 の無い BEND は従来とビット単位で同じ | CALC の解析行列と `tbende` は DISFRIN=1 でも Maxwell のフリンジを入れる (追跡の tbendi は DISFRIN に従う。B134)。逆向きの要素 (`direlc<0`) の F1 の入れ替えは未対応。隣り合う BEND の `bradprev` は tbende では使われない。**TrackParticles の F1 は ANGLE=0 の BEND・MULT で未対応** |
| **ANGLE=0 の BEND (補正磁石・挿入光源の極) と MULT の F1** の Emittance[]・EMIT・CALC | 両 | ステアリング (`tsteee`: fb1 = F1 + FB1、fb2 = F1 + FB2、入口の `tblfre`・断片) と MULT (`tmulte`: F1・F2 は L≠0 かつ K1≠0 のときだけ、`tqfrie`・`tbfrme`・`tqlfre`・`tradke(f1in)` の順、`tbfrmle` の線形フリンジ)。k64 は gs と全桁一致 (12 セル FODO の 14 通りの F1・FB1・FB2・AE1・AE2・FRINGE・DISFRIN・K0・K1・K2・DX・DY・ROTATE の組み合わせ。大きなリング A・PF・ASTRID2・MPW 入りにステアリングまたは MULT を 2 台挿入)。F1 の効果は U0 で −4.4 eV (ステアリング)、運動量圧縮で 2e-3 (MULT) | MULT の AE1・AE2 (読まない)、`FB1/FB2` の `tblfre` (`multae.py`)、ステアリングに ROTATE と DX/DY を同時に付けると gs と CALC のチューンが 1e-4 違う (既存・未調査) |
| **ANGLE=0 の BEND の端の角 AE1・AE2 (B108)** | v2 | 行路を x tan(e) だけ短くする。k64 は K64 のまま (端の項が約 2 倍) | 傾き 5 mrad の極で独立な値と 3〜4%、30 mrad で 2〜3% |
| **RADINT の K0 の極 (B115)・四極のずれ (B126)・BMAG (B107)** | v2 | 上の表と使える範囲の表 | |
| **WSPAC の Emittance[]・EMIT・CALC6D・CALC (RFSW ON)** | Emittance[]・EMIT は両、CALC・CALC6D は v2 | `WSPAC` を ON にして `Emittance[]` を呼ぶと、反復して空間電荷つきの一周の写像からチューン・エミッタンス・Twiss を作る (`PBUNCH`・`MINCOUP`。INTRA とも併用可)。K64 の `twspace`・`tintrb`・`tintraconv` の直訳。0.2 GeV の 12 セル FODO (PBUNCH 1e5〜1e8)。**k64 は gs と全桁一致** (K64 の最後のスライスの 8 倍のレンズ: B127 を再現。小リング、大きなリング A・PF・ASTRID2 で PBUNCH 1e9〜1e10)。v2 は一次のラスレットの積分 (Δν = −(1/4π) ∮ β_u 2 N r0 /(β²γ³ √(2π) σz σ_u (σx+σy)) ds) と 0.2〜1.8%、PBUNCH に比例 (1e5 → 1e7 で比 100 ± 0.3%)。CALC・CALC6D は v2 だけ (k64・gs は入口以外の Twiss が NaN: B133)。CELL の CALC6D の Twiss は反復の最終パスの T から。CALC のチューンは Emittance[] のチューンと 2e-6 | **追跡の WSPAC (`LINE["SIG…"]`) のビーム行列は空間電荷なし・MINCOUP なしのまま**。強い空間電荷 (0.2 GeV で PBUNCH ≳ 1e9) は反復が収束しない。RF なしの CALC。TouschekLifetime との併用。要素ごとのレンズの位置は K64 の近似。4D の CALC は K64 の `qwspac` の直訳ではない |
| **SynchroBetaEmittance の収束の警告** | v2 | 本節の SynchroBeta の項。大きなリング A では使えない (B128) | |

## 5. 注意点 (仕様と落とし穴)

1. **Cavity PHI の既定は 0**。`phic = PHI·sign(charge)`、同期位相は trf0 で入る。
2. **4D CALC では空洞は RFSW に依らずドリフト** (bench_level14a_dbg)。空洞が効くのは Emittance[] と追跡。
3. **薄い多極の dv 項**: tthine 内のドリフトにも `z −= L·dv(δ)` が付く (2-1 節)。表示桁 (εx 5 桁、Level 12) の照合ではこの差は見えないので、照合は全桁 (`$FORM="18.14"`、16b/16e) で行う。
4. **閉軌道の δ**: SAD の tcod は周長変化で δ≠0 を選ぶ (3C: δ=−4.85e-8)。close_z がないと U0・εx がずれる。
5. **RADCOD ON の z と trf0 (Equilibrium dz)**: tcod は一周ごとに z を trf0 に移し、収束判定が緩い (z は 1e-7 まで) ので、分け方は反復の履歴で決まる。`_tcod_radcod` が tcod を直訳しているので trf0 は 実ラティス で 12 桁一致する (trf0 固定の厳密な Newton では dz が 8e-4 mm ずれる)。
   逆に言うと cod と trf0 の「和」だけが物理で、SAD の cod 自体は 1e-9 (dp) 程度しか収束していない (出口と入口がその分違う) — それでも emittance 等は 1e-11 で一致する。
6. **RFSW OFF の RADCOD**: SAD の tcod は RFSW OFF でも RADCOD ON なら放射込みの閉軌道を使う (`_tcod_radcod` の RFSW OFF 枝)。RFSW OFF を RADCOD OFF と同一視すると 実ラティス で U0 が 7e-5、α_c が 1e-4 ずれる。
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
| 25a/25b/25c/25d/25h | ROTATE: 垂直ステアリング (ANGLE=0, ROTATE=−90 DEG, K0)、skew 四極、ROTATE 付き BEND (5/90 DEG, K1 付き)、6 極・8 極・十極・MULT の ROTATE、L=0 の SQ | εx/εy/U0 相対 ≤1e-8 (εy 2e-12〜1e-8)、TUN ≤4e-9 (25h は δ の収束差で 3e-8) |
| 25e | ROTATE/DX/DY/DROTATE/K0 の組み合わせ (癖のない組) | εx/εy 相対 1.6e-6 (δ の収束差)、TUN 4.5e-8 |
| 25f | th≠0 かつずれ≠0 の要素 (tchge の癖) | CALC・追跡は一致、Emittance[] の δ・チューンは SAD の tcod が閉じないため不一致 (2-1 節) |
| 12 | 試験用のラティス (3563 要素) | 表示桁 (相対 ≤2e-5) |
| 16b/16e | 同、全桁 | εx, εz, σE, σz, νx, νy, νs, α_c, bh, trf0, U0 **≤1.3e-11** |
| item6 (`bench_item6_emit.sad`) | 試験用のラティス, RFSW+**RADCOD**, `Emittance[Matrix->True]` の全桁 | εx 1.1e-12, εz 9e-12, σE 2e-14, σz 9e-12, α_c 2e-11, ν ≤9e-12, 減衰率 ≤8e-13, ビーム行列 ≤4e-11, U0・trf0 (dz) は全桁 |
| item6 (`bench_item6_tiny.sad`) | 実ラティス 1 セルのリング (ν_s=2e-4): RFSW OFF / ON / RADCOD | 横方向・U0・trf0 ≤1.3e-11、縦方向 (RFSW ON) は 1e-9 (νs が小さく 1e4 倍に拡大)。テスト `test_emit_jac.py` |
| 26a/26c/26d | EMIT 報告の全行 (実ラティスの RFSW × RADCOD 4 通り、FSHIFT=1e-6)、RFSW OFF + RADCOD ON の Emittance[] と CALC6D | 報告 148 行 + 37 行すべて一致 (雑音の桁を除く)、26d ≤2.2e-11 |
| 26b | FRINGE=1,2,3 付き四極、ずれ・ロール付き四極・6 極 (小リング): CALC と Emittance[] | CALC ≤5e-13、Emittance[] ≤1.3e-9 (数値微分では 1.1e-5) |

対応済み (別モジュール): INTRA (touschek.py: tintrb、TouschekTable/TouschekLifetime)、OneTurnInformation/Matrix 系 (emit_extra.py)、EMITXE/EMITYE/EMITZE の下限、CAVI の DX/DY/ROTATE。
INTRA ON の EMIT 報告 (docs/11 A40): 零電流の報告 → Touschek の表 (運動量アクセプタンス 0.4〜24.2 % の寿命、Nx/Nz・Ny/Nz ごとの寿命/100s) → `Parameters with intrabeam scattering:` (PBUNCH/1e10、MINCOUP %、IBS 込みの EMITXC/EMITYC/EMITZC・SIGE・SIGZ)。PBUNCH の既定は 1e10 (A41)。`ffscmd._print_touschek`、`_print_emit_intra`。
対応済み (4b 節): POL (calpol, spnorm, srequpol)、SOL (tsole)、SpinPrecessionMatrix の表、SynchroBetaEmittance (`pysad/synchrobeta.py`、gs と相対 1e-5)。
未対応: Region (部分区間)、InitialBeamMatrix、TRPT 分岐、`Polarization` (K64 が NaN)、SOL 領域内の MULT・BEND・CAVI。
