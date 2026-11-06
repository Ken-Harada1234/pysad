# 05. CELL/INS の光学計算 (CALC) — SAD の qcell/qcod/qtwiss と pysad の lattice.py

> 最終更新: 2026-11-06 (現状のスナップショット)。SAD ソース (`<SAD>/src`) と `pysad/lattice.py`, `pysad/ffscmd.py`, `pysad/match6.py` を照合して記述。行番号は 2026-10 時点。
> 旧版ノート (`<work-olddocs>/sad_cell_calc_sources.md`, `sad_cod_notes.md`, `sad_twiss_storage_notes.md`) の内容は本書に含む。

## 1. 全体の流れ

```
tffsa (CALC)            FFS コマンド解釈。cell なら tfsetupcell、dp(-nfr:nfr) の設定   tffsa.f:544-559, 1079-1207
 └ tffsmatch            反復ループ (CALC では fitflg=.false. なので 1 回)            tffsmatch.f:139
    └ tffscalc          twmov で第1要素の値 → twiss(1,0,:) (前回の解を初期値に)      tffscalc.f:93
       ├ qcell1(ibound,0,…)   中心運動量                                           tffscalc.f:109
       │   ├ (calc6d → qcell6d)                                                   qcell.f:62
       │   ├ cell: qcod → qmdiag → lb の周期 Twiss・分散                           qcell.f:68-279
       │   └ qtwiss → qtwiss1(mat=.false.) で各要素へ伝搬                          qcell.f:282-290
       ├ 安定なら twmov(…,.false.) で twiss(1,0,:) を第1要素に書き戻す              tffscalc.f:118-120
       └ ix=±1..±nfr: qcell1(ibound,1,…) → utwiss に退避 (オフモメンタム)           tffscalc.f:194-215
```

### qcod (`qtwiss.f:810-970`) — 閉軌道の反復

| 項目 | SAD | pysad (`lattice.py`) |
|---|---|---|
| 1 周の評価 | `qtwiss1(mat=.true.)`: 軌道 cod の上で要素ごとの 4×5 行列 (第5列 = ∂/∂δ) を積算 | `_one_turn_cod_matrix`: 要素ごとの 6×6 Jacobian と出口軌道 |
| 閉じる成分 | (x, px, y, py) のみ。δ は固定、z は閉じない | 同じ (`cod0[5]=delta`, `cod0[4]=0`) |
| 収束判定 | `conv = min(1e-10, 1e-18(1+(δ/0.001)²+nlat))`、残差 `r = dx²/βx + βx(dpx+αx/βx·dx)² + (y 同様)` (1 周のずれの Courant–Snyder 不変量)。α, β は一周行列から `cx=(T11+T22)/2`、`|cx|>1 なら 1/cx` で推定 | `_qcod_periodic` に同式を直訳 |
| 更新 | Newton: `(T−I)·cod0_new = T·cod0 − cod_end` を `tsolvg` (SVD, 閾値 1e-8) で解く | `np.linalg.solve` |
| 減衰 | r が悪化: `fact≥1e-3` なら `cod0=(1−fact)·cod00+fact·cod0`, fact を半減。`fact<1e-3` なら fact×16 で逆方向に外挿。改善時は `fact=min(0.5,2fact)` | 同じ |
| 上限 | itmax=10 (最大 11 周) | 同じ |
| 返す行列 | **最後に 1 周させた出発点 cod0 で評価した trans**。収束後に評価し直さない | `cod_in, M_cod` を返す (同じ) |
| 初期値 | `twiss(lb,idp,DX:DDP)` = tffscalc.f:93 の `twmov(..,.true.)` が読む**先頭 MARK の DX/DPX/DY/DPY/DZ/DDP**。直前の CELL の周期解 (後述 twmov) か、スクリプトが `LINE["DX",..]=` で書いた値 | `FFSState.calc()` が先頭 MARK の DX/DPX/DY/DPY を `twiss_table(cod=…)` に渡す (`Lattice._cod_prev` は MARK がないときの代わり。詳細は `history/work_ins_offset.md`) |

この「評価し直さない」仕様のため、オフモメンタムの β は初期値 (= 直前に計算した運動量点) に依存して相対 1e-6 程度揺らぐ。pysad で SAD と一致させるには、SAD と同じ順 (0 → −nfr … → +nfr) で同じ `Lattice` に対して呼ぶ。

### qtwiss (`qtwiss.f:14-516`)

`qtwiss1` のループは `do l=la+1,lb` で要素 `l1=l−1` を通し、結果を `twiss(l)` に書く。したがって **`twiss(l)` は要素 l の入口**、`twiss(nlat)` (`$$$`) は最終要素の出口である。`$$$` 行の COD は最後の一周の出口 COD。`qmat2twiss` が入口の α, β, R, 分散から出口の値を作り、位相 NX/NY は増分 `dpsix = atan2(T12, aa)` を足していく。`ltyp > icMARK` (MONI, APRT など) は前の値をコピーする。

pysad: `Lattice.twiss()` が `_qcod_periodic` → `periodic_optics` → 分散 (後述) → `_transport_twiss` で各要素行列を伝搬。`ffscmd.FFS.calc()` の `rows[l-1]` が SAD の `twiss(l)` に対応する。

### INS (cell=.false.) と twmov (詳細: `history/work_ins_offset.md`)

流れ:
1. **CELL CALC が安定なら周期解を先頭 MARK に書き戻す** (`tffscalc.f:119` `twmov(1,twiss,nlat,ndim,.false.)`): AX BX AY BY EX EPX EY EPY R1–R4 DETR と閉軌道 DX DPX DY DPY DZ DDP (PSIX=PSIY=0)。以後 `LINE["BX","MKS"]` はこの値を返す。
2. **CALC の入口で MARK の値を読む** (`tffscalc.f:93` `twmov(..,.true.)`): ntfun は `orbitcal` (既定 .true.) か calc6d なら 28 成分すべて (DX..DDP を含む)。NX=NY=0。`DP0` は `LINE["DDP",1]` の別名 (`track.f:39`)。
3. CELL は MARK の DX..DDP を qcod の初期値にして周期解を求める。**INS は qcod を呼ばず**、`qtwiss1` が `cod=twiss(ipa,DX:DDP)` (MARK の値) から要素ごとに軌道を更新し、その軌道上の解析行列 (4×5) で `qmat2twiss` を呼ぶ。R1–R4 は MARK の値、`normal = DETR<1`、`|値|<1e-30` は 0。
4. したがって「CELL → INS に切り替える」と、直前の周期解 (βx=4.643 など) が入口になる。先頭 MARK の DX/DPX/DY/DPY を書き換える (`LINE["DX","MKS"]=…`、amp_dep_ts のやり方) と、それは**初期軌道**になる (分散には入らない)。`BXI/AXI/BYI/AYI` も MARK の値を書く。

pysad: `FFSState._twmov_to_mark(rows)` が CELL CALC (4D/CALC6D) の後に MARK の raw 値へ書く。`_ins_rows` が MARK の AX..R4, DETR と DX..DZ を `lattice.ins_twiss` に渡す。`ins_twiss` は `_one_turn_cod_matrix` の要素処理で軌道を追い、その上の行列を `_transport_twiss` (SAD `qmat2twiss` の直訳、結合込み) で伝搬する。
ベンチ: `bench_level22a..g.sad` (SAD と 1e-9〜1e-10)、esrf `amp_dep_ts` の 1 点目 47.5800386898 (SAD 47.5800386899)。

**INS + CALC6D** (docs/history/work_emit_unstable.md B): `qcell` は calc6d なら CELL/INS に関わらず `qcell6d` に入り (qcell.f:62)、INS (cell=.false.) の枝は閉軌道も tecalc も使わず
先頭 MARK の 28 成分 (AX..ZPY、`twmov(right)`) と軌道 DX..DDP から `ri0=etwiss2ri(twiss(1))`、`twiss(1)=tfetwiss(ri0,cod,.true.)`、`ri=etwiss2ri(twiss(1))` (tfinibeam) を作り、`tturne0` で 6D 伝搬する
(各要素で `tfetwiss(ri·inv(T), cod, normali)`)。MARK には書き戻さない。pysad: `calc6d.calc6d_ins` + `ffscmd._mark_tw28`、CELL の CALC の後に MARK へ縦方向の列 (AZ BZ PSIZ ZX..) も書く。
`bench_level30a.sad` (ESRF、5 通り: 4D の CELL → INS、6D の CELL → INS、MARK にオフセット、初期 Twiss・BZ・DDP を変更、RFSW RADCOD) で SAD と ≤1.6e-10。

## 2. 要素ごとの CALC 用の写像 (`qtwiss.f:184-400` の分岐)

要素パラメータは `lp = merge(ld, le, ideal)`: 通常はインスタンス値、IDEAL フラグで要素 (ファミリー) 値。

| 要素 | SAD の経路 | 内容 | pysad |
|---|---|---|---|
| DRIFT | インライン | `pzi = (1+δ)+dpz`, `a12 = L/pzi + px²·L/pzi³` など正準運動量での厳密ドリフト。`cod(5) += dpz·ale − dvemit·L` | `_drift_matrix_at_cod` + 同式で cod 更新 |
| QUAD | `qquad` → `tquade` | `L=0` なら tthine(nord=4, ak=K1 の積分値、K1<0 の −90° 回転なし。pf2 の SQ `L=0 K1 ROTATE=45 deg` が該当)、`K1=0` なら tdrife0。`K1·L<0` は θ−π/2 回転して `−K1`。`tsolrote(…,.true.)` で DX/DY/ROT の座標へ入り、`tqfrie` (ハードエッジのフリンジ、`DISFRIN=0` のとき)、`tqlfre` (線形フリンジ、`FRINGE=1,3` で入口、`2,3` で出口。係数は `tsetfringepe`: f1in = F1+F1K1F 等、`K1≠0, L≠0` のときだけ)、`tsolque` 本体 (`ndiv = 1+int(|L·K1|/0.1)`)、出口で逆順に戻す。z 写像は `xs = x − sin x` | `_quad_offset_map` → `_tquade_cod` (`_tqfrie_cod`, `_tqlfre_cod`, `_tsolque_cod`)。F1/F2/FRINGE/DISFRIN を `_quad_fringe_params` で読む。逆向き要素の mfr 変換は未対応 |
| BEND | `qbend` → `tbende` | 入口角 `pe1 = E1·ANGLE + AE1` (E1 は ANGLE 比)。`fb1 = F1+FB1`, `fb2 = F1+FB2`。`phib = ANGLE+K0`。`tchge` (DX/DY/DZ/ROT/DROTATE/CHI2) → `tblfre` (F1) → `tbedge` → 本体 → `tbedge` → 逆変換。`K1=0`: `tbendebody0` (F1≠0 なら `tbendef1` の 2 区間)。`K1≠0`: `tbendiinit` + `[tbendecorr + tbendebody]×ndiv` | `_bend_cod_update` (`_tbedge_cod`, `_tbendebody0_cod`, `_tblfre_cod`)、K1≠0 は `emit6d._tbende` を流用した `_bend_k1_cod_matrix` (入口・出口エッジは解析行列に差し替え)。**K1=0 でずれ・F1 のないものと K1≠0 でずれのないものは、行列を SAD の解析式 (`_tbedge*_mat`, `_tbendebody0_mat`) で作る**。ロール誤差は `_tbrote_cod` (`tchge`→`tbrote`、LGEO=L/2, ANGLGEO=ANGLE/2 まわり) |
| SEXT/OCT/DECA/DODECA | `qthin` → `tthine` | `ak=K` (積分値), `kord = nord/2−1`, `aki = ak/kord!`。L≠0 なら `L/6` ドリフト → キック (aki/2) → `2L/3` ドリフト → キック → `L/6` ドリフト。キック: `px −= aki·Re((x−iy)^kord)`, `py −= aki·Im((x−iy)^kord)`。DX/DY/ROT は `tchge` (alg=phig=0) | `_tthine_cod`, `_thin_offset_map`。pysad の K は長さあたりなので `K·L` を渡す (L=0 なら K そのまま) |

| MULT | `qmult1` → `tmulte` (seg なら `qmultseg`) | ANGLE=0, BZ=0 の多極写像 (分割数は tpara の式) | `_mult_cod_map` (`_mult_sequence`) |
| CAVI | `qcav` → `tcave` | `qcav` は `rfsw = rfsw .and. trpt` にして呼ぶので、**リングの 4 次元 CALC では RFSW に依らずドリフト** (bench_level14a_dbg で確認) | `lattice.RFSW=False` で `_drift_matrix_at_cod(L)` |
| SOL | `tfbndsol` で区間判定 → `qsol` → `tsole1`。SOL 要素単体は "implementation error" | 区間扱い | **未対応** |
| STEER | abort ("Use BEND with ANGLE=0") | — | — |
| MARK, MONI 等 | 前の値をコピー | — | 単位行列 |

**ROTATE (ROT/TILT) の扱い (docs/history/work_rotate.md)**: BEND は `tchge` (入口: theta で回す→回した座標系でずらす→tbrote、出口は逆) を外側に付けて中身 (tsteee / tbendebody0 / K1≠0 の tbendebody / tquade / tmulte) を誤差なしで計算する
(`lattice._bend_element_cod_matrix`、`_bend_chg_args` が tbende.f の分岐ごとの theta (ANGLE=0 の tsteee は DROTATE を渡さない、K0=0 はドリフトで tchge を呼ばない) を決める)。
薄い多極 (SEXT/OCT/DECA/DODECA と L=0 の四極) は tthine の先頭・末尾の tchge (`_thin_cod_matrix`)。四極と MULT の ROT は tsolrote (ずらす→回す)。
**SAD の癖**: tchge の入口で th≠0 かつ回した座標系のずれ≠0 (または dz≠0, dtheta≠0) のとき、tbrote が行列を上書きするので行列に回転が入らない (COD は回る)。一周行列が R^T·M になり結合が物理とずれる。`_tchge_jac` が再現 (docs/11 A25)。

### BEND の K0 と ANGLE=0 のステアリング (docs/history/work_bend_k0.md)

- ANGLE≠0 の BEND の K0 は phib = ANGLE+K0 として rhob に入り、`drhob = rhob−rho0` が `tbendebody0`/`tbendebody` の `dpza = dpzi+(drhob−xi)/rhob` でキックになる (px が約 −K0 変わる)。pysad: `_bend_cod_matrix`, `_bend_k1_cod_matrix`, emit6d `_tbende`。
- **ANGLE=0 の BEND (pf2 の K1..K5 など)** は `tbende.f` の phi0==0 分岐: K1=0 → `tsteee` (L=0: 薄いキック `px −= K0`、K0=0: ドリフト、他: 入口エッジ (rhob=+L/K0) → 本体 (rhob=−L/K0) → 出口エッジ)、K1≠0 かつ K0=0 → tquade、K1≠0 かつ K0≠0 → tmulte (akm(0)=K0)。pysad: `lattice._tsteee_cod_matrix`, `_tsteee_body`, `emit6d._tsteee`, `track6._tsteer_ens`。
- `LINE["K0",name]` は次の CALC でファミリーに展開される (docs/02_lattice_main.md)。

### 行列の作り方: 解析行列と数値 Jacobian

SAD の各 `t*e` ルーチンは非線形写像と同時に線形化行列 `trans` を解析的に積算する。pysad は写像を直訳し、COD 位置での 6×6 行列を要素の種類に応じて次のどちらかで作る。

- **解析行列 (SAD の式の直訳)**
  - **BEND** (K1=0・K1≠0 とも、ずれなし・F1 なし): `lattice._tbedgemaxwell_mat`, `_tbedge*_mat`, `_tbendebody0_mat` ほか。SAD の解析行列 (`tbedge`, `tbendebody0` の `trans1` など) は写像の厳密な微分ではなく、px, py≠0 では数値微分と 2e-8 程度ずれる。軌道が大きい INS や δ≠0 の CELL ではこれが Twiss の 1e-7〜1e-6 の差になるので、数値微分ではなく SAD の式を使う。これで 14A の DP0 スキャンのチューンが 1e-10、3A のチューンが 1.5e-10、3B の EY が SAD と一致する。
  - **四極** (ずれ DX/DY・ロール ROTATE・K1<0 の −π/2・ハードエッジ tqfrie・線形フリンジ tqlfre (FRINGE=1,2,3)) と**薄い多極** (SEXT〜DODECA。ずれ・ロールは tchge の行列): `emit6d._jac_quad_full`、`_jac_tqlfre`、`_jac_thin`。lattice の CALC 経路と emit の `_USE_ANALYTIC_JAC` が同じ関数を使う (docs/history/work_emit_report.md)。
    tquade の行列は `R_出口 · [tqfrie(+) · tqlfre(入口) · (tqente + tsolque 本体) ×ndiv · tqlfre(出口) · tqfrie(−)] · R_入口` (ずれは微分に効かない)。数値微分との差は z 行を除いて 2e-12。
    SAD との差: 22e の CELL (DP0=0.004、四極ずれ、閉軌道≠0) **1.4e-11**、23a 1.7e-11、25d 4.5e-11。FRINGE=1〜3 の四極を含む `Emittance[]` (bench_level26b) は 1.3e-9 (残りは SAD の tcod 収束の差。数値微分だとチューンが 1.1e-5 違う)。
- **数値 Jacobian** (中心差分 `_num_jacobian(f, cod, eps=1e-7)`): MULT、ANGLE=0 の BEND (tsteee/tmulte)、F1 付き BEND、DROTATE≠0 の BEND、空洞。刻みは 1e-7 (1e-6 では K1=0 の BEND の `phib`/`psi` の分岐境界 (値が 0 かどうかで経路が変わる箇所) に当たり、(1,6) 成分に 5.6e-7 の誤差が出る)。SAD の解析行列との差は 1e-10 程度。
- dp=0、軌道 0、ずれ・フリンジなしの要素では `el.transfer_matrix()` の解析行列を使う (近道)。

## 3. Twiss の格納

### `twiss(nlat, -1:1, 28)` (`tffs.f:739-742`, module `tffitcode` `tffs.f:101-130`)

| 列 | 名前 (FFS / SADScript) | 内容 | pysad `TwissRow` |
|---|---|---|---|
| 1–3 | AX BX NX / PSIX | αx, βx, 位相 ψx [rad] | alfx, betx, mux |
| 4–6 | AY BY NY / PSIY | y | alfy, bety, muy |
| 7–10 | EX EPX EY EPY | **モード座標の分散** (結合があれば物理分散 PEX とは別) | etx, etpx, ety, etpy |
| 11–15 | R1–R4 DETR | 結合行列 R と det R (モード入れ替え時は `1+xyth−detR`, xyth=0.375 のフラグ) | r1–r4, detr |
| 16–21 | DX DPX DY DPY DZ DDP | 閉軌道 (px, py は正準運動量) | x_cod … delta_cod |
| 22–24 | AZ BZ NZ / PSIZ | 縦方向 (CALC6D) | alfz, betz, muz |
| 25–28 | ZX ZPX ZY ZPY | z 結合 (CALC6D) | zx, zpx, zy, zpy |

第 2 次元は idp=0 が中心運動量、+1 がオフモメンタムの作業層 (結果は `utwiss` に退避)。位相 NX は rad で、`tffsa.f:141` の scale=2π で割った値がチューン。`Twiss["NX","$$$"]/(2π)` = Qx。

### 周期解の分散 (`qcell.f:267-279`)

"Dispersion is defined in 2*2 world": R と `c=√(1−detR)` から S を組み、`A, B` ブロックと `(a13,a23) = S·T(:,5)` を作って
`EX = 0.5(a13 + a12·a23 − a22·a13)/dcosμx`, `dcosμx = 2 sin²(μx/2)`。x_cod/δ の差分ではない。pysad `Lattice.twiss()` に同式を直訳 (`lattice.py:1800-1870`)。

### オフモメンタム (DP0 / DP)

`tffscalc` が `twiss(1,1,DDP) = twiss(1,0,DDP)+dp(ix)` として `qcell1(ibound,1,…)` を呼び、初期軌道は中心運動量の閉軌道。pysad `match6._calc_all_dp` は `dp0 + dps[idp]` で同じ順に `_calc_optics` を呼ぶ。

### 端数位置 (`Twiss["BX", {"Q1", 0.3}]`)

SAD: `qtwissfrac` → `qtwissfrac1` → `qfraccomp` で部分要素を作って通す。`r = rx2−rx1` として ANGLE·r, K0·r, K1·r, K(THIN)·r, VOLT·r。BEND の E1 は `rx1=0` のときだけ残し (E1·f1/r)、E2 は `rx2=1` のときだけ (それ以外 0)。FRMD を `−f1−2·f2` にして端の側だけフリンジを付ける。
pysad `ffscmd.row_frac`: `L·fr`, `ANGLE·fr`, `E1` はそのまま, `E2=0` の部分要素の**線形行列**で入口から伝搬する近似 (K1 のスケール、フリンジの切り替え、閉軌道上の非線形写像は未対応)。

## 4. 安定性判定

| 場面 | SAD | pysad |
|---|---|---|
| qcod のオーバーフロー / 非収束 | `over` (一周の途中で行列の和が NaN) → `*****qcod---> Overflow & closed orbit not found n`、twiss を書かずに戻る (前回の値が残る、条件はその値で評価、nstab は両面 +2)。収束しないだけ (`codfnd` 偽) → `*****qcod---> Closed orbit not found n` を出して stab=偽にし、Twiss は続けて計算する。メッセージは nmes<4、反復表を 1 行出すたびに 0 に戻る。qcod は一周行列を `limitnan(·,1e10)` で ±1e10 に抑える | `twiss(sad_stab=)` が一周行列を ±1e10 に抑え (`_TRMAX`)、`_LAST_QCOD_FOUND` で codfnd を返す。非有限・成分が 1e18 超 (SAD は COD の丸め雑音 1e-18 が増幅されて NaN になる) は RuntimeError → `Matcher._calc_optics` が同じメッセージ (nmes 付き) を出し、`_stale_rows` で前回の rows の残差を評価。`tsolvg` は特異に近いときだけ `tsvd` |
| qmdiag | `(TrA−TrB)² < 0` → `***qmdiag---> Sum resonance: …` (lfno=0 で次の表示まで 1 回)、`stab=.false.` (twiss(lb) は以前の値のまま); NaN も不安定。trab2 の絶対値が 1e-14 以下は特殊 6 分岐 | `lattice._qmdiag` (全分岐の直訳)。通常の場合は従来の `periodic_optics` (ビット単位で不変) |
| `|cosμ| > 1` (`qcell.f:209-248`) | 常に `cosμ := 1/cosμ` に置き換えて α, β を計算し続ける。`INTRES` (cosμ>1) / `HALFRES` (cosμ<−1) が立っていれば `stabx=.false.`、立っていなければ trace も `2cosμ` にして安定扱い。分散は stabx で dcosμ≠0 のときだけ書く (不安定な面は以前の値のまま) | `lattice._qcell_periodic_sad` (直訳)。不安定な面の分散の「以前の値」は `Lattice._mark_init` (ffscmd が先頭 MARK の値を入れる)。**FFS のフラグは `NOINTRES;` と書く (`OFF INTRES;` は SAD では無効)** |
| 伝搬 (qtwiss1) | 要素に入る前に βx, βy > 0 でなければ over: 直前の正常な行のコピーで以降を埋め stabx=staby=偽 | `twiss(sad_stab=)` の伝搬 (`_transport_twiss(clamp=False)`) |
| 伝搬後の再判定 (`qcell.f:351-397`) | lb と終点の Twiss から trace=x11+x22 を作り直し、CELL では再び INTRES/HALFRES 判定 (`fam` では判定しない。立っていなければ trace=4/trace)、最後に `stab .and. (codfnd .or. fam)` | `twiss(sad_stab=)` (`Lattice._qstat` の tracex, tracey, stabx, staby) |
| 不安定の表示 | `tffsa.f:1229-1247`: 全 idp の optstat に偽があれば `Horizontal` / `Vertical` / `Horizontal/Vertical` + ` Unstable.`。CALC/GO の条件表に `1/$$$ BMAGX/BMAGY` が加わる (tfshow: cellstab が偽)。nstab は残差が 0 でない idp ごとに偽を数える (条件が何もなければ 0 で `Matched.`) | `FFSState.print_unstable`、`print_fit_table`、`Matcher._tffscalc` (`cellstab`, `last_stab`) |
| 整数共鳴帯をまたぐ (`tffscalc.f:196-213`) | オフモメンタム点で `aint(NX/2π)` (INTRES)、`aint(NX/π)` (HALFRES)、`aint((NX+NY)/2π)` (SUMRES)、`floor((NX−NY)/2π)` (DIFFRES) が中心運動量と違えば unstable | `match6.py:405-432` に直訳 (既定 INTRES/HALFRES/SUMRES=T, DIFFRES=F) |

安定時のみ (全 idp が安定のとき) `twmov(…,.false.)` で第1要素に書き戻すので、次の CALC は最後に安定だった解から始まる。pysad も `_cod_prev` と `_mark_init` を安定時だけ更新する。詳細・ベンチ (Level 28) は docs/history/work_unstable.md。

## 5. ベンチマークの要約 (CALC 関連)

| Level | 内容 | SAD との差 |
|---|---|---|
| 1 | 20 セル FODO: 閉軌道 (全 QD に DX=1e-4)、放射積分 I1–I5 | 閉軌道 入口 x −110.8245 µm 一致、I1–I5 表示桁一致 |
| 2 | オフモメンタム Twiss (dp=0, ±0.01, ±0.02)、複合機能 BEND | βx, αx, ηx ≤3e-8 / βy, Qx, Qy ≤6e-11 |
| 3A | 六極と色収差 (DP0 スキャン 21 点 + 3 次当てはめ) | チューン ≤2.6e-9、ξ 相対 ≤3e-7 |
| 3B | 六極を通る COD (DELX QD.1, DELY SF.1) | 閉軌道 9–10 桁、β/α/η/Q ≤2.4e-8、EY 8e-10、R1–R4 5–6 桁 |
| 3C | x–y 結合 (DTHETA QF.3, B.5) | 閉軌道 8–9 桁、光学 ≤1.6e-8 (EX/EPX は 3e-11)、R1–R4 7 桁 |
| 3D | 強さの誤差 (DK, DDK) | 閉軌道 9–10 桁、光学 ≤9e-8 |
| 5 | INS 輸送ライン、F1/FRINGE/DISFRIN のフリンジ | Twiss ≤2e-8、軌道 ≤1.3e-11、位相 ≤5e-11 |
| 12 | 実ラティス ESRF 20cells/v03_12 (3563 要素、NLAT=3564) の CELL CALC | Twiss・チューン ≤4e-11 |
| 13/14 | マッチング後の Twiss、esrf/calc_chromaticity の DP0 スキャン | K1 ≤1.1e-10、チューン ≤2e-8 |
| 28 | 不安定ライン (ESRF の QD/QF を ×1.19, ×0.9, …): CELL CALC の続行、`Unstable.`、GO の nstab、INTRES/HALFRES | Twiss ≤2.5e-10 (β≲1e5)、GO の反復表・K1 一致 (docs/history/work_unstable.md) |

## 6. 対応状況と既知の未対応

対応済み (本書の該当節を参照):
- **CALC6D** (`qcell6d`、列 22–28): `calc6d.py` (docs/history/work_calc6d.md)。INS + CALC6D は `calc6d_ins` (1 節)。
- **`|cosμ|>1` の 1/cosμ 置き換え**と伝搬後の再判定: `twiss(sad_stab=)` / `lattice._qcell_periodic_sad` (4 節)。
- **結合時の途中要素の位相**: `qmat2twiss` の coup 分岐を `_transport_twiss` に移植済み。

未対応:
- **GEO 系** (`qtwissfracgeo`, `tfgeo`, GX/GY/GZ/CHI1–3): 幾何座標は計算しない。
- **SOL の CALC** (`qsol`/`tsole1`、BZ≠0 の区間): 未移植 (追跡 track6 では対応済み)。
- **モード入れ替え**: 一周の位相の小数部は `periodic_optics` のモード (qcell と同じ) だが、`periodic_optics` はモード入れ替え (`normal=.false.`) を扱わない。
- **端数位置** `row_frac`: 線形行列の近似 (qfraccomp の K1·r、フリンジ切り替え未対応)。
- **INS の fit の一部**: `FIT @ (MARK value)` 構文、逆向き要素の mfr 変換、`BXI` 以外の初期値コマンドの一部。
- **解析行列がない要素** (数値微分): ずれ・ロール付き BEND (DROTATE≠0)、MULT、F1 付き BEND、ANGLE=0 の BEND (四極と薄肉多極は解析行列。2 節)。
- 四極の逆向き (`dir<0`) での `mfr` 変換、MULT の ANGLE≠0 (`tmulta`)、ドリフトなしで続く要素が一緒に動く仕様。
