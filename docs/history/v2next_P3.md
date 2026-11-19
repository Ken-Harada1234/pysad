# v2next の記録: 担当 P3 (TrackParticles の F1 (tsteer・tmulti)、MULT の AE1・AE2、B108 の残り、2026-11-19)

対象 (`docs/46` の「TrackParticles の F1」「MULT の AE1・AE2」、`docs/47` §(d)・`v2next_E.md` §3.1・§3.2 の未解決): (1) TrackParticles の F1 (ANGLE=0 の BEND の `tsteer`、MULT の `tmulti` の F1・F2・FB1・FB2・FRINGE)、
(2) MULT の端の角 AE1・AE2、(3) B108 の残り (K1≠0・K0≠0 の端の角、追跡 (`tsteer`) の端の放射、傾き 30 mrad の極の 0.02〜0.03)。
変えたコード (追跡の核とその周辺): `pysad/track6.py` (`_tsteer_ens`・`_tblfri`・`_tmulk`・`_tmulti_ens`・`_bend_k1_proxy`)、`pysad/track6n.py` (numba の `k_steer` と `Prog` の K_STEER の引数)、
`pysad/track6n_ext.py` (`compile_mult`: F1・FRINGE・FB・AE つきの MULT は Python の track6 に落とす)、`pysad/lattice.py` (`_mult_params`・`_bend_as_mult`・`_mult_sequence`)、`pysad/emit6d.py` (`_tmulte`・`State.tradke`)、
`pysad/perfkern.py` (`tradke_db` の行路の補正に y の係数)。numba のカーネルを編集したので `py -3.12 -m pysad.nbcache` でキャッシュを掃除した。
試験 (新規 124 件): `pysad/tests/test_v2_P3_f1.py` (100 件)・`test_v2_P3_ae.py` (24 件)、補助 `_p3_indep.py` (独立な基準・gs と同じ台本の組み立て)・`_p3_gsdata.py` (gs の値)。
例: `example/v2_next/P3_*` (`P3_README.md`)。gs は 2026-11-19 の 10-19 ビルド。**gs の追跡・Emittance[] は数秒で終わった (timeout 300〜900 s)。終了後に `ps` で sad1.exe の残りは無かった。**

## 1. 項目の状態 (要約)

| 項目 | 状態 | モード | 根拠・試験 |
|---|---|---|---|
| 1. TrackParticles の `tsteer` の F1・FB1・FB2・FRINGE (ANGLE=0 の BEND) | **実装した**。k64 は gs と 10 桁表示の全桁一致 (6 通り + 大きなリング A)。numba の高速追跡と Python は 1e-13 | 両モード (K64 の動作の追加、gate なし) | `test_tracking_matches_gs_k64`・`test_steer_f1_fast_equals_python`・独立な基準 (§2) |
| 1'. `tmulti` の F1・F2・F1K1F/B・FB1・FB2・FRINGE=1,2,3 (MULT) | **実装した** (Python の `_tmulti_ens`。F1 等つきの MULT を含むリングは track6n が track6 に落ちる)。k64 は gs と全桁一致 (9 通り + 大きなリング A) | 同上 | 同上 |
| 2. MULT の AE1・AE2 (CALC・Emittance[]・追跡) | **実装した**。pysad は MULT の AE を読んでいなかった (黙って無視)。k64 の CALC・Emittance[] (`tmulke` の apsi)・追跡 (`tmulk`) は gs と全桁一致 (5 通り + 大きなリング A) | k64・v2 とも (K64 の動作の追加)。v2 だけ出口の蹴りの位置 (§3.2) と行路 (B108) | `test_mult_edge_angle_emittance_matches_gs_k64`、`test_tracking_matches_gs_k64[m_f,m_g,m_i]` |
| 3a. B108: 追跡 (`tsteer`・`tmulti`) の端の放射 | **v2 で修正**: 端の角のある面で場の中の行路が x·tan(e) 短い (放射の `zr0`)。厳密な幾何 + ロビンソンの定理と比 1.03〜1.04 (追跡の一周写像の固有値) | v2 | `test_edge_angle_term_v2_matches_exact_geometry_emit_and_tracking` |
| 3b. B108: MULT (K0 つき) の Emittance[] | **v2 で修正**: `tmulte` の最初・最後の `tradke` に tin・tout (回した座標系・ずれつきは y の係数も)。比 1.03〜1.04 (K1 = 0 と K1 = 0.06 の勾配つきの極) | v2 | 同上、`..._for_combined_function_poles`、`..._rotated_mult_emit_equals_tracking_v2` |
| 3c. K1≠0・K0≠0 の ANGLE=0 BEND の AE | CALC・Emittance[] (tbende → tmulte) は FB・AE を渡す (pysad は渡していなかった)。追跡 (tbend → tmulti) は K64 が psi = 0 を渡すので k64 は AE を無視、v2 は AE を使う。**gs はこの BEND で USE の中で SIGSEGV (`work_rotate.md` の既知) するので K64 の値は確認できない** | k64 (K64 のソースの読み)・v2 | `test_bend_k1_k0_edge_angle_tracking_uses_ae_in_v2_only` |
| 3d. 傾き 30 mrad の極の 0.02〜0.03 | **切り分けた**: 傾き (K0) には依らない。AE に比例する 2.7% の系統差で、極の長さ・間隔に依る (§3.4)。原因は未特定 (**要判断**) | — | `P3_b108_exact.py` |
| 新しい K64 の誤りの候補 (番号なし) | §5 | — | — |

## 2. 項目 1: TrackParticles の F1 (tsteer・tmulti)

### 2.1 実装

| 欄 | 内容 |
|---|---|
| 状況 | v2.1.0: `track6._tsteer_ens`・`track6n` は `F1/FB1/FB2` で NotImplementedError、`_tmulti_ens` は F1・FRINGE で NotImplementedError (Emittance[]・CALC は E で実装済み) |
| K64 の動作 | `tsteer.f`: tpara が fb1 = F1+FB1、fb2 = F1+FB2 (FRINGE ≤ 0 は −1 (入口だけ)・−2 (出口だけ) 以外 0)。入口は dxfr1 = −fb1²/ρb/24、dyfr1 = fb1/ρb²/6、dyfra1 = 4 dyfr1/fb1² (DISFRIN=0)、放射ありは長さ f1r = fb1/2 の 断片 n = −1, 0 (`tbendal`: 偏向 rbl·f1r・rbh·f1r、放射の長さ f1r)、出口は断片 ndiv+1, ndiv+2 と dxfr2 = −fb2²/ρb/24。`phic = phib(1 − (f1r+f2r)/al)`、本体 alc = al·phic/phib。`tmulti.f`: 入口 `ttfrins` (mfring ≠ 2) → `tblfri(akr0, fb1)` (mfring = 1, 3) → 四極の線形フリンジ f1in・f2in (mfring = 1, 3) → `tmulk(psi1)` → 分割 → 出口は逆順 (`tmulk(psi2)` → f1out (`−f1out`)・`tblfri(−akr0, fb2)` (mfring = 2, 3) → `ttfrins` (mfring ≠ 1))。追跡の f1in 等は `tsetfringep` の係数 (f1 = −\|K1/L F1²\|/24、f2 = \|K1/L\| F2。Emittance[] の `tsetfringepe` (F1・F2 そのもの) とは別) |
| 実装 | `_tsteer_ens` (Python) と `elkern` 相当の numba `k_steer` (引数に fb1・fb2・DISFRIN・tan(AE2)・v2 を追加) に tsteer の直訳。`_tblfri`・`_tmulk` を追加して `_tmulti_ens` に F1・F2・F1K1F/B・FB1・FB2・FRINGE・AE1・AE2 を入れた。K1≠0・K0≠0 の ANGLE=0 BEND (tmulti の代理) は FB1・FB2 (F1+FB) と mfring を渡す。`track6n_ext.compile_mult` は F1・FRINGE・F2・F1K1F/B・FB・AE のどれかがある MULT で NotImplementedError にして track6 に落とす (1 周 1 粒子 0.4 s 程度の遅さ。大きなリング A に MULT 2 台を入れた追跡 3 周は数秒) |
| k64 モードの扱い | K64 の動作の追加 (F1 つきの追跡は NotImplementedError だった) なので gate なし。F1 を使わない追跡はビット単位で変わらない |
| gs との照合 | k64 は gs と 10 桁表示の全桁一致: 12 セルの FODO (3 GeV、周長 45.6 m) のステアリング 6 通り (F1・FB1・FB2・FRINGE = 3・AE1・AE2、FRINGE = 1/2、−1/−2、F1 つき FRINGE=0 (効かない)、AE だけ、E1・DX・ROTATE・DISFRIN)、MULT 9 通り (F1、F1+F2、FRINGE = 1/2、K0+FB1+FB2、K0+FB+K2、AE+K0、AE だけ、FRINGE=2/1+K2、F1K1F/B+AE+K0+FB+ROTATE+DX)、NORAD と RAD の 3 粒子 3 周。大きなリング A (`kek-ls_v03_12.lat`、3 GeV、周長 603 m、3563 要素) の先頭に補正磁石 2 台 (K0 = ±5e-4、L = 0.2 m): BEND (F1・FB1・FB2・FRINGE=3・AE)・MULT (K1 = ±0.05、F1・F2・FB・AE)・AE だけの BEND で追跡・CALC・Emittance[] が一致 (差 ≤ 1e-8 相対) |

### 2.2 独立な基準 (v2・k64 とも同じ。F1 の追跡は K64 の動作の追加)

| 基準 | 内容 | 結果 |
|---|---|---|
| 軸上の偏向 | ∫B dl は磁場の立ち上がりの形・δ に依らない: p0 で規格化した px の合計 = −K0 (F1・FB・FRINGE −2..3 の 7 通り × δ = 0, 0.01, −0.02) | 1e-14 |
| 放射の積分 | 磁場が長さ fb で直線に 0 → B0 と立ち上がると ∫B² dl = B0² al (1 − (fb1+fb2)/(6 al))。F1 なしとのエネルギー損失の比 | 解析式と 1e-5 (3 通り、高速追跡と Python) |
| 鉛直の収束 | x 方向に一様な Maxwell 場 (By = B(s) − y²B''/2、Bs = yB'(s)) の端の厳密な運動 (ベクトルポテンシャル Ax = F(s) − y²B'/2 の正準方程式、DOP853 rtol 1e-13)。Δpy/y (F1 あり − なし) = (fb1+fb2)/(6ρ²) (解析式 y∫F B' ds) | 厳密な運動は 1e-3 (積分誤差)、K64 は 1e-4 (3 通り)。MULT の FB1・FB2 (`tblfri`) も同じ |
| 水平の x・px | 厳密な運動 (出口の面に戻す)。fb1 = fb2 のとき δ = 0, ±0.01, −0.02 で x が 1.3e-11、px が一致 | 1e-10 |
| 追跡と CALC | 追跡の一周写像の固有値 (閉軌道のまわりの 4 次元の数値微分) と CALC (tsteee・tmulte。別の経路) のチューン。F1・FB・AE つきのステアリング、MULT (K1・K0・K2・F1・F2・FB) | 2e-6 以内 (k64・v2) |
| 高速追跡と Python | numba の `k_steer` と `_tsteer_ens` (NORAD・RAD・FLUC、乱数の最後の状態も、12 通り) | 1e-13 |

試験: `test_steer_on_axis_deflection_is_k0_for_any_f1`・`test_steer_radiated_energy_equals_ramp_integral`・`test_steer_f1_vertical_focusing_matches_maxwell_ramp`・`test_steer_f1_horizontal_matches_maxwell_ramp_for_equal_fb`・
`test_mult_fb_vertical_focusing_matches_maxwell_ramp`・`test_steer_fringe_zero_ignores_f1`・`test_steer_f1_tracking_tune_equals_calc`・`test_mult_f1_tracking_tune_equals_calc`・`test_steer_f1_fast_equals_python`。例: `P3_maxwell_ramp.py`。

### 2.3 使える範囲・使えない範囲 (項目 1)

| 欄 | 内容 |
|---|---|
| 使える範囲 | ANGLE=0 の BEND (K1 = 0) の F1・FB1・FB2・FRINGE (−2..3)・DISFRIN・AE・DX・DY・ROTATE の追跡 (高速追跡も)。MULT の F1・F2・F1K1F/B・FB1・FB2・FRINGE (1・2・3)・AE・ROTATE・DX・DY・K0・K2.. の追跡 (RF なし)。3 GeV の FODO と大きなリング A (3 GeV、603 m) で確認 |
| 使えない・未確認 | (1) 逆向きの要素 (`direlc < 0`) の F1 の入れ替え (BEND の F1・FB1・FB2、MULT の入口・出口) は未対応 (Emittance[] と同じ)。(2) **K64 の F1 の模型は fb1 ≠ fb2 のとき x が厳密な Maxwell 場の端の運動と (K0/al)(fb1² − fb2²)/24 だけ違う** (K64 の仕様で直さない。fb1 = fb2 では 1e-10 で一致。例: K0 = 0.003、al = 0.2、fb = (0.3, 0.0001) で 5.6e-5 m、`P3_maxwell_ramp.out.txt`)。(3) MULT の F1 つきの追跡は Python の経路 (高速追跡には落ちない)。(4) RF つき MULT (`tmultiacc`) の F1 は NotImplementedError のまま。(5) スピン追跡 (POL) の MULT の F1 の端の項 (K64 の tmulti は fringe の bsi を足さない) は確認していない。(6) 厳密な Maxwell 場の鉛直の運動 py/y = −B0² al は MULT (K0 のみ、tmulti) も BEND (tsteer) も再現する (K0 = 0.003・0.02 で相対 5e-6 以内) |

## 3. 項目 2・3: 端の角 AE1・AE2 (MULT、追跡、B108 の残り)

### 3.1 状況と K64 の動作 (gs で確認)

| 欄 | 内容 |
|---|---|
| 状況 | v2next E: pysad は MULT の AE1・AE2 を `_mult_params` で読まず、CALC・Emittance[]・追跡のどれでも黙って無視していた (`tbende` から tmulte に回る K1≠0・K0≠0 の ANGLE=0 BEND の AE も) |
| K64 | Emittance[]・CALC: `tmulte` が `tmulke(apsi)` で端の蹴りを入れる (入口は tsolque の前、出口は **最後の半スライス (tsolque) の前**)。追跡 `tmulti`: `tmulk(psi)` (入口は分割の前、出口は **最後の半スライスの後**)。k64 の pysad は AE を読むようにして、CALC・Emittance[] 5 通り (K0 = ±0.01 の強い極を含む。減衰が反転するリング)・追跡 (MULT 3 通り)・大きなリング A が gs と全桁一致 (出口の順序も) |
| 追加の確認 (gs) | `P3_steer_ae2_tune.sad`: 追跡の `tsteer` は AE2 を **fb2 ≠ 0 のときだけ**入れる (CALC・Emittance[] の `tsteee` は常に入れる)。AE2 だけのステアリングの追跡の一周写像のチューン (gs) は AE なしと同じ 0.1340952、FB2 = 1e-6 を付けると 0.1338835 (CALC は AE2 を入れた値) |

### 3.2 v2 の動作

| 欄 | 内容 |
|---|---|
| (a) 追跡の AE2 | v2: `tsteer` は AE2 を常に入れる (K64 は fb2≠0 のときだけ)。CALC・Emittance[] と同じにする。根拠: 追跡の一周写像の固有値 (閉軌道のまわり) = CALC のチューン。K64 の動作は AE2 = 0.3・K0 = 0.0015 のステアリングで 1.2e-4、大きなリング A (ra3: AE2 = 0.3) で 1.4e-4 違い、v2 は 6e-6 (K0 の 2 次、§3.5) |
| (b) MULT の出口の蹴りの位置 | v2: Emittance[]・CALC の `tmulte` は出口の端の蹴りを最後の半スライスの後 (放射の刻みの前) に置く (追跡の `tmulti` と同じ。物理的にも出口の面は本体の後)。k64 は K64 のまま (半スライスの前)。根拠: MULT (K1 = 0.7・−0.6、K0 = ±0.01、AE2 = 0.3) の追跡と CALC のチューンの差が v2 で 3e-8、k64 で 3e-5 (AE1 + AE2 で 1.2e-4)。AE1 だけなら両モードで 1e-7 以内 |
| (c) K1≠0・K0≠0 の ANGLE=0 BEND の追跡 | v2: 追跡も AE を使う (K64 は psi = 0 を渡す: ソースの読み)。根拠: 同じ BEND の CALC (tbende → tmulte) のチューンと追跡が一致 (2e-6。K0 の 2 次)。k64 は 5e-4 以上違う。K64 の値は gs が SIGSEGV するので確認不可 |
| (d) B108 の放射の行路 | v2: 端の角 e のある面で場の中の行路が x·tan(e) 短い。(1) 追跡: `tsteer` と `tmulti` の `zr0` (tradk の al1 = al − z + zr0) を入口で −tan(AE1)·x、出口で −tan(AE2)·x (K1 < 0・ROTATE の MULT は回した座標系の x')。(2) Emittance[]: `_tmulte` の端の蹴りの後の最初の `tradke` に tin・出口の蹴りの後の最初に tout (回した座標系・DX・DY は x' = cos θ2 (x − dx) − sin θ2 (y − dy) の係数で渡す: `tradke_db` に y の係数 tiny・touty を追加)。k64 は K64 のまま (行路を変えない) |

### 3.3 独立な基準: 厳密な幾何 + ロビンソンの定理 (v2next E と同じ構成)

24 セルの FODO (3 GeV、周長 100.8 m、QF D2 B D2 QD DS [極 K0 DW 極 −2K0 DW 極 K0] DS D1、極 L = 0.2 m、DW = 0.05 m、K0 = 0.005 rad) の極の AE による ΔJx (AE なしとの差) と、硬い端の磁場の極の厳密な運動 (面は実験室系で z = x tan a。円弧と直線、K1 ≠ 0 は B(x) = (K0 + K1 x)/L の微分方程式の数値積分) の Je = 2 + dln(Σ∫B² dl)/dδ、Jx = 3 − Je。
(`P3_b108_exact.out.txt`。比 = pysad / 厳密)

| 極 | AE | 厳密 ΔJx | v2 Emittance[] | v2 追跡 | k64 Emittance[] | k64 追跡 |
|---|---|---|---|---|---|---|
| MULT K1 = 0 | (0.3, 0.3) | +0.01440 | 1.029 | 1.029 | 2.089 | 2.039 |
| MULT K1 = 0 | (−0.2, 0.4) | +0.00574 | 1.040 | 1.040 | 2.012 | 2.064 |
| MULT K1 = 0 | (0.3, −0.3) | −0.00061 | 1.072 | 1.071 | 1.197 | 2.021 |
| MULT K1 = 0.06 (勾配つき) | (0.3, 0.3) | +0.01296 | 1.027 | 1.027 | 2.080 | 2.030 |
| MULT K1 = 0.06 | (−0.2, 0.4) | +0.00517 | 1.039 | 1.040 | 2.028 | 2.058 |
| BEND ANGLE=0 | (0.3, 0.3) | +0.01440 | 1.027 | 1.030 | 2.077 | 0.975 (AE2 を入れない) |
| BEND ANGLE=0 | (−0.2, 0.4) | +0.00574 | 1.033 | 1.042 | 2.156 | −1.581 |

v2 は MULT・BEND・勾配つきの極のどれでも、Emittance[] と追跡の両方が厳密な幾何と 3〜4% 以内 (K1 ≠ 0 の勾配つきの極は新しく確認)。K64 は端の項が約 2 倍 (B75 と同じ構造)。くさび (AE1 = −AE2) は 1 次が打ち消す。
K1 < 0 の MULT は K64 が多極全体を主軸へ回す (K64 の約束: 端の面の傾きの軸が実験室系の y になる) ので実験室系の x に傾く厳密な幾何とは別の幾何: Emittance[] と追跡の ΔJx が差 3e-4 以内で一致することだけを確認 (`..._rotated_mult_emit_equals_tracking_v2`)。

### 3.4 残る 2.7% の切り分け (3d: 「傾き 30 mrad の極で 0.02〜0.03」)

`P3_b108_exact.out.txt` の末尾 (v2 の Emittance[]、BEND の極):

| 条件 | 比 (v2 / 厳密) |
|---|---|
| AE = (0.05, 0.05)・(0.1, 0.1)・(0.3, 0.3)・(0.5, 0.5) (K0 = 0.005) | 1.0276・1.0275・1.0274・1.0261 (AE の大きさに依らない = tan(e) に比例) |
| EMITDIV = 1・5・20・80・200 | 1.0275・1.0269・1.0274・1.0256・1.0278 (分割に依らない) |
| K0 = 0.0025・0.005・0.03 (AE = 0.1) | 1.0285・1.0275・1.0257 (傾き (K0) に依らない) |
| 極の長さ LP = 0.1・0.2 m、間隔 DW = 0.05・0.3 m (K0 は B を同じにする) | 1.009 (LP = 0.1)・1.027 (0.2)・1.105 (DW = 0.3) |

**v2next E の「傾き 30 mrad で 2〜3% 違う」は 30 mrad に特有ではなく、K0 に依らない同じ 2.7% の系統差だった**。Emittance[] と追跡が同じ値になること (比 1.029 と 1.029)、分割・K0・AE の大きさに依らないこと、極の長さ・間隔に依ることから、K64 の模型 (端の蹴り + 行路の補正) の中の幾何的な項 (隣り合う面の配置や x' の影響) が候補。
追跡の ΔJx + ΔJy + ΔJe は 0 でなく +2.7e-4 (AE (0.3, 0.3)、Jy = 1.00004 → 1.00011): 放射の減衰の和則が AE のある面で 2.7e-4 破れている。原因は特定していない (**要判断**: 厳密な幾何が仮定する端の面の配置と K64 の写像の違い、または K64 の端の蹴りの 2 次の項)。実用上は AE のある極の Jx・εx を 3% 以内の精度で与える。

### 3.5 使える範囲・使えない範囲 (項目 2・3)

| 欄 | 内容 |
|---|---|
| 使える範囲 | MULT・ANGLE=0 の BEND・勾配つきの極の AE1・AE2: CALC・Emittance[]・追跡 (k64 は gs と全桁、v2 は厳密な幾何と 3〜4%)。K0 = 0.005 rad/極 (軌道の傾き 5 mrad)〜0.03 rad/極 (30 mrad) の極、AE ≤ 0.5 rad で確認。大きなリング A (3 GeV、603 m) で、追跡のチューンと CALC のチューンが ra1 で 1.1e-7、ra2 (MULT) で 8e-8、ra3 (AE だけ) で 1.8e-7 以内で一致 (v2。`P3_ring_a.out.txt`) |
| 使えない・未確認 | (1) 追跡と CALC の差は K0 の 2 次で残る (K64 の `tsteer`・`tmulti` と `tsteee`・`tmulte` の端の蹴りの 2 次の違い): AE = 0.3 のステアリングで K0 = 0.0015・0.003・0.006 に対して 4e-6・1.6e-5・6e-5 (K0 の 2 乗)。K64 の仕様で直さない。(2) 残る 2.7% (§3.4)、和則の破れ 2.7e-4。(3) 実験室系で y に傾く端の面 (K1 < 0 の MULT・ROTATE のある MULT) の厳密な幾何との比較は無い (K64 の約束は要素と一緒に面を回す)。(4) BEND ANGLE≠0 の MULT (`multae.py`) は変えていない (ANGLE≠0 は元から AE を読んでいた)。(5) K1≠0・K0≠0 の ANGLE=0 BEND の K64 の値 (gs が SIGSEGV)。(6) 逆向きの要素の AE1・AE2 の入れ替え (`direlc < 0`)。(7) POL (スピン) つきの MULT の端の角の bsi・スピンの端の項 (K64 の `tmulke` の bsir0 は AE を含まない。未確認) |

## 4. 大きなリングでの確認 (大きなリング A、3 GeV、周長 603 m、3563 要素)

`example/v2_next/P3_ring_a.py` (`.out.txt`)。先頭に補正磁石 2 台 (K0 = ±5e-4、L = 0.2 m、周長 +0.4 m) を入れた:

- k64: 追跡 (NORAD・RAD、3 粒子 3 周) が gs と ra1・ra3 で 3e-26 以内 (10 桁表示で全桁一致)、ra2 (MULT) で差 0。CALC のチューン・Emittance[] のチューン・U0・運動量圧縮・減衰率・エミッタンスが 1e-7 以内 (表示桁の丸め)。
- v2: 追跡の一周写像のチューン (4 次元) と CALC のチューンの差 (追跡 − CALC、x・y): ra1 (BEND、F1・FB・AE): −9e-9・−1.1e-7 (k64 も同じ。fb2≠0 なので K64 の追跡も AE2 を入れる)、ra2 (MULT): −4e-8・−8e-8 (k64 は −4e-7・+1.2e-6)、ra3 (AE だけ): +8e-8・−1.8e-7 (**k64 は +1.4e-4・−1.1e-4**)。
- 追跡の時間: ra1・ra3 は numba。ra2 (MULT の F1 つき) は Python に落ちる (この例の全体 (追跡 3 周と閉軌道の数値微分) で数秒)。

## 5. K64 の新しい誤りの候補 (B 番号は付けない。親が番号を付ける)

| 候補 | 内容 | 根拠 | 状態 |
|---|---|---|---|
| (a) `tsteer` の出口の端の角 | 追跡の `tsteer` は AE2 を fb2 ≠ 0 のときだけ入れる (`if(n .eq. n2 .and. fb2 .ne. 0.d0)` の中で tanp2 を設定)。CALC・Emittance[] の `tsteee` は常に入れる | gs: AE2 = 0.3 のステアリングの追跡のチューン (0.1340952) が AE なしと同じ。大きなリング A (AE2 = 0.3) で追跡と CALC のチューンが 1.4e-4 違う | v2 で修正 (AE2 を常に入れる) |
| (b) `tmulte` の出口の端の蹴りの位置 | Emittance[]・CALC の `tmulte` は出口の端の蹴り (`tmulke` の apsi2) を最後の半スライス (`tsolque`) の前に置く。追跡の `tmulti` は後 (`tmulk`) | gs の CALC (pysad k64 が全桁一致) と追跡 (全桁一致) のチューンの差 3e-5〜1.2e-4 (K0 = ±0.01・AE = 0.3 の MULT) | v2 で修正 (CALC を追跡と同じ後ろに) |
| (c) K1≠0・K0≠0 の ANGLE=0 BEND の追跡の AE | `tbend` が `tmulti` に psi = 0 を渡す (`0.d0,0.d0,0.d0,0.d0` = bz・phia・psi1・psi2)。`tbende` は `tmulte` に AE (apsi) を渡す | K64 のソースの読み。gs はこの BEND で SIGSEGV するので実行で確認できない | v2 で修正 (追跡も AE を使う)。要確認: ソースの読みが正しいか |
| (d) B108 の MULT への拡張 | `tmulte`・`tmulti` の端の蹴りは放射の行路を変えない (B75 と同じ)。端の項が約 2 倍 | 厳密な幾何 + ロビンソンの定理 (表 3.3) | v2 で修正 |
| (e) 補正磁石の ROTATE + DX/DY の CALC | gs の CALC の νy (ROTATE = 0.1・DX = 1e-4 のステアリング、12 セルの FODO) が 0.1258111 で、**gs の追跡 (pysad k64 と全桁一致) の一周写像のチューン 0.13371 と違う** (νx は 0.2668346 対 0.2668366 で合う)。pysad の v2 の CALC は追跡と 5e-7 で一致 (v2 は B09 の修正が入る)。k64 の pysad の CALC は gs と 1e-4〜7e-4 違う (docs/46 の「ROTATE と DX/DY を同時に付けると gs と CALC のチューンが 1e-4 違う」の未調査分。pysad k64 の CALC が K64 を再現していない) | gs の追跡と CALC の不整合 (B09 と同じ型: ロール + ずれの tchge の行列) | 追跡の核の外 (CALC の k64 の再現は未調査のまま)。v2 は修正済み |

## 6. 要判断・残した項目

1. K64 の F1 の模型は fb1 ≠ fb2 で x が厳密な Maxwell 場の端の運動と (K0/al)(fb1² − fb2²)/24 違う (§2.3)。K64 の仕様 (模型の近似) で直さない。直すなら入口の面の取り方 (磁場の立ち上がりの中心を要素の端に置くか) の判断が要る。
2. B108 の 2.7% の系統差と、AE のある面の放射の和則の破れ 2.7e-4 (§3.4)。原因未特定。
3. K1≠0・K0≠0 の ANGLE=0 BEND (gs が SIGSEGV) の追跡・CALC の AE の K64 の値の確認。
4. 逆向きの要素 (`direlc < 0`) の F1・AE の入れ替え。
5. 実験室系で y に傾く端の面の MULT (K1 < 0・ROTATE)、K64 の約束 (面が要素と一緒に回る) が意図かの判断。
6. MULT の F1 つきの追跡の高速化 (numba の MULT カーネルは F1・FRINGE・FB・AE に対応していない)。

## 7. 全試験・回帰 (2026-11-19)
全試験 (`py -3.12 -m pytest pysad/tests -q --ignore=test_bend_feeddown.py --ignore=test_coord.py --ignore=test_sext_radiation.py`、`PYSAD_SAD_PACKAGES` つき、k64 が既定): **1512 通過・8 skip・失敗 0** (28 分。基準 1345 通過・8 skip に P3 の新規 124 件と、同時に走った他の担当 (P1・P2・P4・P5) の新規試験が加わる)。
回帰 (`QUICK=1`): **50 ベンチ FAIL 0** (132 s)。変えた既存の試験は無い (期待値の追従なし)。MULT の AE を読むようにした (k64 の動作の追加) が、既存の試験・ベンチに AE つきの MULT・K1≠0・K0≠0 の ANGLE=0 BEND は無かった。

## 8. 付けた B 番号 (v2.2.0 の文書化で。docs/21)

- §5 の (a) `tsteer` の出口の AE2 = **B140**、(b) `tmulte` の出口の蹴りの位置 = **B141**、(c) K1≠0・K0≠0 の ANGLE=0 BEND の追跡の AE = **B142**、(d) B108 の MULT への拡張 = **B143**、(e) 補正磁石の ROTATE + DX/DY の CALC = **B144**。
- §6 の 1 (F1 の模型が fb1 ≠ fb2 で厳密な Maxwell 場と違う) = **B145**。B108 の更新 (2.7% の系統差) は B108 の行。
- 要判断は docs/47 の v2.2.0 の (d)。
