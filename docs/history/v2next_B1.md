# v2next 担当 B1 の記録 (2026-11-18): 光学・Emittance・SOL まわりの未対応の実装と確認 (6 項目)

試験: `pysad/tests/test_v2_B1_optics.py` (42 件)・`pysad/tests/test_v2_B1_sol.py` (45 件)。A1 の試験 `test_emittance_warns_about_ignored_bend_f1` は動作が変わったので `test_emittance_no_longer_warns_about_bend_f1` に置き換えた。
例: `example/v2_next/B1_*` (説明は `B1_README.md`。台本と pysad の v2・k64 の出力、gs の出力 `.gs.txt`)。gs の確認は `local/wake/gsrun.sh` (10-19 ビルド、2026-11-18)。互換モードの規則は `local/wake/v2_brief.md`。
確認に使ったリング: (小) 3 GeV・周長 45.6 m の 12 セル FODO、0.2 GeV の同形のリング (空間電荷)。(大) 大きなリング A (`pysad/tests/kek-ls_v03_12.lat`、3 GeV、603 m、3563 要素)、PF リング (`example/pf/lattice/v9_16_asym.sad`、2.5 GeV、187 m)、ASTRID2 (`example/astrid/lattice/ASTRID2_v01_01.sad`、`ASTRID2_v01_02_id.sad` の `RINGMPW`、0.58 GeV)。大きなリング B (8 GeV) は使っていない (項目 1 の gs 比較は A・PF・ASTRID2 で足りた)。
**numba のキャッシュ**: `track6n.py`・`track6n_ext.py` を編集したので `pysad/__pycache__` の `track6n*`・`multa*` の `.nbi`・`.nbc` を消した (消さないと古い結果が出る)。

## 一覧

| 項目 | 結果 | モード |
|---|---|---|
| 1 Emittance[]・EMIT・CALC6D の BEND の F1・FRINGE | **実装した**。k64 は gs と全桁一致 (小リング 9 通り、大きなリング A・PF・ASTRID2・ASTRID2+MPW の全偏向磁石に F1 を付けて) | 両モード (gate なし。K64 の動作の追加) |
| 2 K1≠0 の BEND の CALC と追跡の 1e-5 の差 | **誤りは追跡 (tbendi) の薄い端の写像**。CALC は厳密な写像と 1e-11 で一致。追跡の端を CALC と同じ厳密な tbedge にして v2 で直した | v2 だけ。k64 は K64 のまま |
| 3 光学・Emittance への空間電荷 (WSPAC) | **実装した** (Emittance[]・EMIT・CALC6D・RFSW ON の CALC)。一次のラスレットの式と 0.15%。gs の Emittance[] より 4% (x)・3% (y) 弱い (原因未特定) | Emittance[]・EMIT は両モード、CALC・CALC6D は v2 だけ (gs は NaN) |
| 4 CALC6D・RADCOD つきの CALC と SOL、INTRA と SOL | **実装した**。k64 は gs と 1e-9。RADCOD の CALC6D の後の 4D の CALC (DP0 の別名) も直した | 両モード (gate なし) |
| 5 CALC (4 次元) の RADTAPER | **実装しない (理由を調べた)**: K64 の 4D の CALC は RADTAPER を見ない (gs で確認)。CALC6D は k64 でも RADTAPER を入れるようにした (gs と一致) | 4D は仕様、CALC6D は両モード |
| 6 SOL 領域の POL (MULT・CAVI)・加速 + RAD・bz ≥ 3 | 確認した (bz = 20 まで、SpinTune の解析式 1e-9、放射損失 1e-3)。**足りなかった DRIFT の刻みを足した** (スピンが 1e-2 ずれていた)。POL + TRPT の加速は明示のエラーにした | v2 だけ (DRIFT の刻み)。ほかは確認の試験 |

## 項目 1: Emittance[]・EMIT・CALC6D の BEND の F1・FRINGE

| 欄 | 内容 |
|---|---|
| 状況 | A1 の発見: `emit6d._tbende` は F1 = 0 を前提に黙って F1 なしの値を返す (gs と U0 で 5%、チューンで 0.006 違う。CALC は F1 を入れるので CALC と Emittance[] のチューンが食い違う) |
| K64 の動作 | `tbende.f`: `alcorr` の長さ補正 `al = al0 − ((φ fb1)² + (φ fb2)²)/48/al0 · sin(½(φ − pe1 − pe2))/sin(½φ)` (`mfring≠0`)、入口の `tblfre` (FRINGE>0 か −1)、`f1r = fb1/2`・`f2r = fb2/2` (FRINGE>0 か −2)、本体の長さを `rbc = 1 − (f1r+f2r)/al0` 倍。K1=0 は n = 0 に入口の断片 2 つ (`tbendef1`: 長さ f1r·rbl、f1r·rbh の `tbendebody0`、それぞれ放射の長さ alr = f1r)、本体 n = 1..ndiv、n = ndiv+1 に出口の断片 2 つ (放射なし。放射は出口の後の 1 回 alr = f2r)。K1≠0 は n = −1, 0 と ndiv+1, ndiv+2 の `tbendebody` (`tbendal`)。出口の後に `tblfre` (dxfr2 = −fb2²/ρb/24)。**周長** `circ = pos(nlat) − pos(1)` (tfgeo の pos) は F1 の長さ補正を含む (Emittance[] の ω0・運動量圧縮・周回周波数) |
| v2 の動作 | k64 と同じ (gate なし。K64 の動作の追加)。`emit6d._tbende_core` に tblfre (`_tblfre`、K64 の trans1 のとおり trans1(4,6) は dyfra の項を含まない近似)・`_tbendal`・断片のループを入れた。F1 の無い BEND は今までと同じ演算 (ビット単位で同じ)。`_emit_circ` は F1・FRINGE のある BEND に `geo._geo_length` を使う。F1 (FRINGE≠0) のある BEND は numba の `k_bend` を使わず Python の経路。`ffscmd` の A1 の警告を削除。`touschek_table` の周長も `_emit_circ` に |
| 根拠 | gs (Emittance[]、RFSW RAD NOFLUC): 12 セル FODO (3 GeV、F1=0.1 FRINGE=3): Tunes {−.2933271700、−.1395704757、−.0179539450}、U0 739673.742271 eV (F1 なし 781633)、EquilibriumPosition .0282503369、MomentumCompaction .0326541595、DampingRate・Emittances・BunchLength も全桁。FRINGE = 1・2・3 は同じ値、−1 (入口側だけ)・−2 (出口側だけ)、F1+FB1+FB2、K1=−0.05 (tbendebody の分割)、K1=0.04、K0=1e-3、RADCOD、DX・DY・ROTATE・DROTATE つき (表示の丸めまで)、POL ON (減衰率)。検算 (gs を使わない): 放射・RF なしの Emittance[] のチューン = CALC のチューン (1e-9。F1・K1≠0・FB1/FB2 つき) |
| 大きなリング | 大きなリング A (全 BEND に `F1=0.03 FRINGE=3`: U0 285476.566791 eV (F1 なし 296074.86)、Tunes {−.4199675211、−.3846720252、−.0040789829}、運動量圧縮・減衰率・εx が gs と 1e-7 以内)。PF リング (`F1=0.05 FRINGE=3`: U0 392531.84006)・ASTRID2 (`F1=0.04 FB2=0.02 FRINGE=3`: U0 6045.54954484)・ASTRID2+MPW (`F1=0.04 FRINGE=2`: U0 7704.31936504) も Tunes・U0 が gs と 1e-9。RADCOD (PF・ASTRID2) も 1e-7 |
| 試験 | `test_emittance_bend_f1_matches_gs_k64[8 通り]`、`test_emittance_bend_f1_radcod_and_pol_match_gs_k64`、`test_emittance_tunes_equal_calc_tunes_with_bend_f1[k64・v2 × 5]`、`test_bend_f1_changes_u0_in_v2`、`test_large_ring_a_bend_f1_emittance_matches_gs_k64`、`test_large_ring_a_bend_f1_v2_consistent`、`test_other_rings_bend_f1_emittance_matches_gs_k64[3]`、`test_emittance_no_longer_warns_about_bend_f1` (A1 の試験を置き換え) |
| k64 モードの扱い | gs と全桁一致 (表示 10 桁) |
| 使える範囲 | ANGLE≠0 の BEND (K1 の有無、FRINGE = 1・2・3・−1・−2、F1・FB1・FB2、K0、DX・DY・ROTATE)。Emittance[]・EMIT・CALC6D・INTRA (touschek)・Touschek。CALC (4D) は元から対応 |
| 使えない・未確認 | (1) **ANGLE=0 の BEND (tsteee: ステアリング・挿入光源の極) の F1・FB1・FB2 は `NotImplementedError` のまま** (K64 の tsteee.f の tbendal の断片の直訳が要る。CALC も同じ)。(2) MULT の FB1・FB2 の tblfre (`multae.py`) は NotImplementedError のまま。(3) 逆向きの要素 (`direlc<0`) の F1 の入れ替え。(4) v2 の絶対値は A の項目 (端の角の行路の補正) で gs と 1e-3 違う (k64 は全桁一致)。(5) 隣り合う BEND の `bradprev` は tbende では使われない (tradke に渡らない) ので扱っていない。(6) CALC の解析行列と tbende は DISFRIN=1 でも Maxwell のフリンジを入れる (追跡の tbendi は DISFRIN に従う。K64 の不整合、直していない) |

## 項目 2: K1≠0 の BEND の CALC (tbende) と追跡 (tbendi) の 1e-5 の差

| 欄 | 内容 |
|---|---|
| 状況 | A2: 軌道が 1 mm ずれた K1≠0 の BEND の x-px の行列が CALC と追跡で 1e-5 違う (DX 無しでも、K1=0 は 1e-13 で一致)。再現: L=1.2 m・15 度・K1=−0.05・E1=E2=0.5·ANGLE・x0=1 mm・δ=3e-3 で 1.1e-5 (E=0.25: 2.7e-6、E=0: 4e-8) |
| 調べたこと (独立な厳密な写像) | **厳密な写像**を、ローレンツ力 dp/dl = u × B の数値積分 (DOP853、rtol 1e-13) で作った (`test_v2_B1_optics._exact_map`)。磁場の形は SAD の BEND の模型: 双極 B_y = h は端面 (半径方向から E1・E2 回った平面) の間、勾配 B_y = k x + c k h x²、B_r = k y は半径方向の入口・出口の面の間 (**K64 の tbende は端面で四極の磁場を入れず、本体の扇形の長さ (半径方向の面の間) だけに四極が働く模型**)。c = −1/2 は `tbendicorr` の蹴り −K1 x² h/2 から決めた (e=0 で CALC の行列と 1e-8 で合うことで確認)。端面は硬い端でドリフトで繋ぐ |
| 結果 | **CALC (tbende) の写像 (x, px) は厳密な写像と 1e-11 で一致** (軌道 x0 = 1・3 mm、px0 = 5e-4、δ = 5e-3、E = 0.5・0.25)。K64 の追跡 (tbendi) は 1e-10〜9e-9 ずれる (x0=3 mm で 1.6e-9、x0=1 mm・δ=5e-3 で 9.2e-9、px0 = 5e-4 で 1.1e-9。端の角に比例して増える。e=0 では CALC と 1e-10 で合う)。行列 (x0 = 1 mm の 4x4): CALC 対 追跡の数値微分の差は K64 で 1.1e-5 (M12 が 9.2e-6)。**誤りは追跡の側** |
| 原因 | `tbendi` の端は薄い写像 (`pxi = px + tanp1·x/ρe`、`pyi = py − (tanp1 + fpx) y/ρe`、`xi = x + ff`) で、端面までのドリフトと端面から扇形の面までの弧 (x tan e の項) を入れない。CALC の `tbedge` (tbedgedrift → tbedgemaxwell → tbedgebody) は厳密。K1=0 の追跡 `tbend` (tbendcore) は厳密な閉じた式なので CALC と合う |
| v2 の動作 | `track6n.k_bendi` (numba) と `track6._tbendi` (Python) の入口・出口の端を、`DISFRIN=0` かつ端の角 ≠ 0 のとき `lattice._tbedge_cod` と同じ厳密な tbedge にする (`track6n.edge_exact1`: 追跡の px・py は p で割った量なので p を掛けて tbedge に渡し、戻す)。DISFRIN=1 か端の角 = 0 は K64 の薄い写像のまま。MULT の曲がった要素 (multa の k_bendi の呼び出し) にも入る |
| 根拠 (数値) | v2 の追跡: 厳密な写像と 1e-10 以下 (x0 = 1 mm で 6e-13、3 mm で 1.8e-11、δ = 5e-3 で 4.6e-11)。CALC との行列の差 4e-8 (x0 = 1 mm、δ = 3e-3、E=0.5・0.25・0)。k64: 1.1e-5・2.7e-6 (記録) |
| 試験 | `test_k1_bend_calc_map_equals_exact_hard_edge_map[E=0.5・0.25]` (CALC)、`test_k1_bend_tracking_map_equals_exact_in_v2_not_in_k64[2]`、`test_k1_bend_calc_matrix_vs_tracking_derivative` (v2)・`..._k64` |
| k64 モードの扱い | K64 のまま (1e-9 台の差。gs と同じ) |
| 使える範囲 | K1≠0 の BEND (複合機能磁石) の追跡と CALC・Emittance[] の整合 (E1・E2 ≠ 0 の矩形・中間の磁石、δ・軌道が大きいとき)。大きなリング A の BPI6・BPI7 など |
| 使えない・未確認・別の差 | (1) **四極の磁場が端面でなく半径方向の面で始まる模型**は CALC・追跡とも同じ K64 の仕様 (直さない)。端面で四極を入れる厳密な硬い端の磁石とは、行列で x0 tan e に比例する差 (x0 = 1 mm、E=0.5 で 2.3e-5) が残る。(2) **CALC の解析行列は軌道が y ≠ 0 のとき自身の写像の数値微分と y = 2 mm で 1.5e-6 違う** (K1=0・K1≠0 とも。tbfrie の解析行列の近似)。K1=0 の追跡 (tbendcore) と CALC の写像の微分も y = 2 mm・E=0.5 で 1.5e-6 違う (どちらが正しいかは調べていない: Maxwell のフリンジの y² 項)。v2 の K1≠0 の追跡は CALC の写像と 1.3e-7 で合う。(3) DISFRIN=1 の K1≠0 の端は K64 の薄い写像のまま |

## 項目 3: 光学・Emittance への空間電荷 (WSPAC)

| 欄 | 内容 |
|---|---|
| 状況 | 旧: v2 の WSPAC は追跡だけで、ビーム行列は空間電荷なし (ゼロ電流) の平衡値。Emittance[]・CALC のチューンに空間電荷のシフトが出なかった。K64 のソース: `tintrb.f` の `twspace` (ビーム行列 σ から条件付きの横の分散 xx1 = σxx − σxz²/σzz、yy1、xy1 を作り、主軸 σ1・σ2 と角度から f = 2 N r0 L/(p1² h1 (σ1+σ2) √(2π σzz))、kx = f/σ1、ky = f/σ2 の thin レンズ)。`tintrb` が スライス中心 (`trans1` のドリフトで送る) にレンズを入れる (`trans2 = trans1⁻¹ · transsp · trans1`、`trans`・`beam(1:21)` だけに掛ける。`tmulbs(..,.false.)`)。temit の反復 `tintraconv` (INTRA と同じ。`econv = .not.(intra .or. wspac)`) で B2 = 前回の平衡ビーム行列を自己無撞着にする。`qtwiss` の CALC は `qwspac` で要素ごとのレンズ (長さ = 中点間) を入れる |
| K64 の動作 | Emittance[]・EMIT: 反復して空間電荷つきの一周の写像から Tunes・Emittances・Twiss を作る。MINCOUP の下限 εy ≥ MINCOUP (εx+εy) が入る。**CALC・CALC6D (WSPAC ON) は入口以外の Twiss が NaN になる** (gs で確認: 0.2 GeV のリングで PBUNCH = 1e4 でも「Horizontal/Vertical Unstable」と NaN)。ビームサイズ (`beamsize`) は EMIT の最後の plot の一周から |
| v2 の動作 | `touschek.touschek_table(.., intra=False, wspac=True)` に `_twspace` と `_tintrb` の wspac の分岐を入れた (INTRA とも併用可)。ffscmd: `Emittance[]`・EMIT は WSPAC ON なら反復 (`run_emit`、`emittance_rules(sc=…)`)。空間電荷つきのビーム行列は `touschek_table` の `BeamMatricesSC` に残すが、**`LINE["SIG…"]`・追跡の WSPAC が使うビーム行列は今までどおりゼロ電流の平衡値** (MINCOUP の下限なし。`test_v2_C2_spch` の WSPAC の試験の前提。K64 は空間電荷つきの値と MINCOUP の下限)。CELL の CALC6D と、RFSW ON の 4D の CALC は 6D の経路 (`calc6d(sc=…)`: 反復の最終パスの T_{1→l} と正規形から Twiss。縦方向の列は 4D の既定値)。**CALC・CALC6D は v2 だけ** (k64 は今までどおり入れない: gs が NaN) |
| 根拠 | (a) 一次のラスレット: 0.2 GeV、PBUNCH = 1e5、MINCOUP = 0.5: Δν = −(1/4π) ∮ β_u 2 N r0 /(β²γ³ √(2π) σz σ_u (σx+σy)) ds (σx² = β εx + η² σδ²、σy² = β εy。円形で分散なしなら公称の ラスレット Δν = −N r0 C/(4π β²γ³ ε √(2π) σz))。要素を 30 分割したリングで pysad の Δν (−4.51e-6、−9.48e-6) / 式 = 1.0008・1.0006 (**0.1%**)。PBUNCH に比例 (1e5 → 1e7 で比 100 ± 0.3%)。(b) CALC6D・CALC の Twiss: BX($$$) = BX(1) (周期的)、チューンが Emittance[] のチューンと 2e-6 |
| gs との比較 | Emittance[] (PBUNCH = 1e8、MINCOUP = 0.5): ZERO (WSPAC なし) は全桁一致。**SC のチューンシフトは gs の方が 4.3% (x)・3.0% (y) 強い** (gs −4.355e-3・−8.848e-3、pysad −4.176e-3・−8.594e-3)。比は PBUNCH (1e7・1e8)・RF 電圧 (0.2・0.8 MV: σz が 2 倍)・MINCOUP (0.5・0.2) に依らず一定 (1.043・1.030)。gs のプロット後のビームサイズ (SIGX・SIGY) に K64 の twspace の式を当てても gs のシフトは 2% 強い。**原因は未特定** (tintrb の式・`tmulbs(,.false.)`・B2 の作り方・al の分け方は K64 のソースを読み直して一致するが、gs の動作が 2〜4% 強い)。ラスレットの式 (教科書の式) とは pysad が 0.1%、gs が 2〜4% |
| 試験 | `test_wspac_emittance_matches_gs_within_known_difference_k64`、`test_wspac_tune_shift_equals_first_order_laslett_formula` (v2)、`test_wspac_calc6d_and_calc_rows_have_the_space_charge_tune_shift_v2`、`test_wspac_calc6d_k64_keeps_ignoring_space_charge` |
| k64 モードの扱い | Emittance[]・EMIT は空間電荷つき (gs と 4% の差、上)。CALC・CALC6D は今までどおり (空間電荷なし) |
| 使える範囲 | 電子 0.2 GeV の 12 セルのリング (RF あり、MINCOUP の下限つき)。低電流 (Δν ≲ 1e-2 で反復が収束)。TouschekLifetime との併用は未確認 |
| 使えない・未確認 | (0) 追跡の WSPAC のビーム行列は空間電荷なし・MINCOUP なしのまま (上)。(1) gs と 4% の差 (上)。(2) **RF なし (σz が決まらない)・垂直エミッタンス 0 の平らなリング**は MINCOUP の下限で動くが、CALC (4D) への経路は RFSW ON が要る。(3) 強い空間電荷 (0.2 GeV で PBUNCH ≳ 1e9: 反復が収束しない・不安定)。(4) 要素ごとのレンズの位置は K64 の tintrb (スライスの中心を drift で送る近似。磁石の中では厳密でない)。(5) 追跡の WSPAC は強-弱のまま (ビームサイズ固定)。(6) 大きなリング・実機のリングでは確かめていない。(7) 4D の CALC (qwspac の中点間のレンズ) は K64 の式の直訳ではなく CALC6D と同じ 6D の経路 |

## 項目 4: CALC6D・RADCOD つきの CALC と SOL、Touschek (INTRA) と SOL

| 欄 | 内容 |
|---|---|
| 状況 | SOL 領域があると SADError (CALC6D・RADCOD つきの CALC)・NotImplementedError (INTRA)。領域の写像は `solmap` (CALC)・`solemit` (Emittance[]) にあった |
| K64 の動作 | `qtwiss` の 4D は SOL 領域を `qsol` (`tsole1(.. enarad=.false.)`) で、CALC6D は `qtwissfrac` の `tsole1(…, .true., .true.)`。INTRA の tturne は領域内の要素を tsole1 で通す |
| v2 の動作 | 両モード (gate なし。K64 の動作の追加): `_check_sol_optics` は CALC6D・RADCOD でも領域を annotate して通す。`calc6d._turn`・`touschek._tturne_intra` に `_solmap` の要素を `solemit.emit` で通す分岐を足した (`emit_extra` と同じ)。**CELL CALC6D のあと `DP0` を先頭の `DDP` (閉軌道の δ) にする** (track.f: DP0 は LINE["DDP",1] の別名。RADCOD の CALC6D のあとの 4D の CALC は、そのエネルギーで光学を計算する。SOL なしの RING でも gs と一致: チューン 2.70642542976106 (DP0 なしの値 2.70675785)) |
| 根拠 | gs (SOL 領域のあるリング、3 GeV、BZ = 2 T、回転四極、QFE の DX・DY): CALC 4D 2.73475252978800・1.87531224259842 (pysad k64 2.73475252978794)、CALC6D 2.73495437251036・1.87554377922481・NZ −.01173531080122 (2.73495437251039 など。1e-9。NZ は 1e-9 絶対)、RADCOD つき CALC6D 2.73497629185364 (1e-10)、RADCOD の CALC6D の後の 4D の CALC 2.73465818826401 (pysad 2.73465818826743)、β (1・30・70 番) も 1e-9。INTRA + SOL: Emittances {2.706114556595e-7、2.245979910482e-8、2.778701736555e-5} (pysad の εx・εy は 2e-8 以内、εz は 9e-8 の相対差)。TouschekLifetime は gs・pysad とも NaN (この条件) |
| 試験 | `test_sol_region_calc6d_radcod_calc_match_gs_k64`、`test_calc6d_radcod_sets_dp0_for_the_following_4d_calc_k64`、`test_sol_region_intra_emittance_matches_gs_k64` |
| k64 モードの扱い | gs と 1e-9 |
| 使える範囲 | SOL 領域 (領域内の DRIFT・QUAD・MULT・SOL。`docs/44` §3) のあるリングの CALC (4D・6D・RADCOD)・Emittance[]・EMIT・INTRA。v2 では領域内の BEND・CAVI も (solmap の対応範囲) |
| 使えない・未確認 | Touschek の寿命 (TouschekLifetime) の値は SOL のあるリングで確認していない (gs が NaN)。BOUND=0 の SOL (BZ≠0) の CALC は NotImplementedError のまま。RADTAPER + SOL は未確認。WSPAC・INTRA + SOL の併用は走るが値の確認なし |

## 項目 5: CALC (4 次元) の RADTAPER

| 欄 | 内容 |
|---|---|
| 状況 | docs/46: 「RADTAPER: CALC (4 次元) のテーパーは未対応」 |
| K64 の調査 | `qtwiss.f` の `rtaper` は `qtwissfrac` (CALC6D の `calc6d` の分岐、行 1011: `rt=radtaper`、`tsole1` に `rt1` を渡す) にだけ。4 次元の qtwiss の本体 (要素ごとの行列) に rtaper は無い (`grep rtaper qtwiss.f`)。`tcod.f` (閉軌道。radcod .and. radtaper) と `temit.f`・`tfemit.f`・`tturne.f` (6D の追跡) で使う。**gs で確認**: RFSW・RAD・RADCOD・RADTAPER ON の CELL CALC (4D) は、RADTAPER なしの CALC とチューン 2.70675785189057・1.86687229750657、β、DDP (0) が全桁同じ。pysad も同じ値 (両モード) |
| 結論 | **実装しない** (K64 の仕様で、実装すると K64 と別の値になる)。理由: 4 次元の光学には放射によるエネルギー変化 (閉軌道の δ の s 依存) が無く、テーパー (局所の運動量に磁石の強さを合わせる) は何も変えない。テーパーの光学は CALC6D (6D の閉軌道の δ(s) を使う) で出る。`CALC6D; CALC;` の動作が使える (下) |
| 見つけた差 (直した) | gs の **CALC6D + RADTAPER はテーパーを入れる** (チューンが Emittance[] の値と同じ 2.70673352、DDP {0、−1.0854e-4、−2.3872e-4})。pysad の k64 は CALC6D にテーパーを入れていなかった (v2 だけが入れていた。コメントの「gs は Unstable optics」は今のビルドでは当たらない)。k64 も `_radtaper_lattice()` を使うようにした。gs との差は Emittance[] と同じ 8e-8 (k64 の RADTAPER の Emittance[] の gs との差)。**`DP0` の別名**(項目 4) で RADCOD の CALC6D のあとの RADTAPER の CALC6D は DP0 = 1.28e-4 を閉軌道の探索に使う (k64、gs のチューン 6.696366849 と 1.2e-8。使わないと 1.3e-7)。v2 は CALC6D が書いた DP0 を RADTAPER の探索に使わない (Emittance[] の軌道と 1e-9 で整合する C2 の試験を保つ) |
| 試験 | `test_calc_4d_ignores_radtaper_like_k64[k64・v2]`、`test_calc6d_radtaper_matches_gs_k64` |
| 使えない・未確認 | SOL 領域 + RADTAPER (未確認)、RFSW OFF の RADTAPER (従来どおり「テーパーなしで続ける」) |

## 項目 6: SOL 領域内の POL (MULT・CAVI)・加速 + RAD・bz ≥ 3 1/m

| 欄 | 内容 |
|---|---|
| 状況 | docs/46 §3・docs/44 §7a の「領域内の未確認」: bz ≥ 3 1/m、領域内の MULT・CAVI の POL (Emittance[])、加速 + RAD、RF つき MULT の RAD |
| bz ≥ 3 | 領域内の DRIFT・BEND・MULT (K0・K1・K2) の写像 (力学的運動量) がローレンツ力の積分と bz = 6・10・20 1/m で 1.2e-8 以内 (多極の分割 EPS の分)。放射損失 −(2/3) r_e γ³ (K0/L)² L が bz = 6・20 で 4e-5。bz = 3 は C の試験 |
| **足りなかった実装 (v2 で直した)** | **領域内の DRIFT のスピン・放射の刻み**。領域内の BEND・MULT は bz·al/ndiv ≤ 0.02 に刻んでいた (K64 の tradk は 1 刻みの前後の向きの変化 (弦) から損失とスピンの歳差を出すので、らせんの回転角が大きいと不正確) が、DRIFT は 1 刻みだった。0.4 m の DRIFT・軸を外れた粒子 (δ = 1e-3・−5e-3) のスピンが Thomas-BMT の積分と、bz = 4.5 (0.1 GeV) で 1.5e-3、9 で 8.7e-3 (k64 の記録では 1.6e-2)、0.3 GeV の bz = 3 で 1.3e-4 ずれた (軸上の粒子は厳密)。`track6n_ext.k_sol_drift` (numba) と `track6._sol_drift_v2` (Python のスピン) に v2 だけ ndiv = 1 + int(|bz·al|/0.02) の刻みを足した。結果: 2.3e-7・8.9e-7・4.5e-8。**放射も**: 横運動量 θ0 のらせんの放射損失 −(2/3) r_e γ³ (bz θ0)² L に 1e-4 (3 GeV、bz = 3・9)。K64 の 1 刻みは bz = 9 で 1.7 倍外れる (k64 の試験で記録) |
| 領域内の QUAD・MULT のスピン (bz ≥ 3) | 大きな bz で誤差が増える: QUAD (K1 = 0.3、0.4 m) は bz = 4.5・9 で 2.0e-5・9.0e-5、MULT (K0 + K1) は 1.3e-4〜1.6e-4、BEND は 6e-7・9e-7。放射の分割が EPS 由来で bz·al/ndiv ≤ 0.2 のため (QUAD)。DRIFT と同じ刻みを足せば小さくなるが、v2 の既存の試験値 (1 GeV・bz ≤ 0.6 で 5e-5) を変えないため足していない。bz = 9 でも誤差 1e-4 台 |
| POL の Emittance[] (MULT・CAVI) | 平面リング (24 台の偏向磁石、0.5 GeV: aγ = 1.135、0.1 GeV: aγ = 0.227) の SOL 領域 (DRIFT・MULT (K1 = 0・0.05)・VOLT=0 の CAVI) の SpinTune が cos(μ/2) = cos(π ν0) cos(φ/2)、φ = (1+a) bz L (L = 領域の長さの合計) と **1e-9 以内で一致** (bz = 0.6・4.5・9 1/m)。ν0 は同じリングの BZ=0 の SpinTune (aγ との差 2.7e-5 (0.5 GeV)・3.9e-4 (0.1 GeV) はこのリングの周長・幾何の値で、Emittance[] 自身の BZ=0 の値を基準にして 1e-15 で一致)。例: 0.1 GeV・bz = 9 の drift 0.341980851748 (解析 0.341980851748) |
| 加速 + RAD | 領域内の RF つき MULT (TRPT、VOLT = 0〜100 MV、基準運動量 +0〜3.1%) の放射損失 (RAD ON − OFF の δ の差) が −(2/3) r_e mc² (bz γ0² θ0)² L/E2 (横運動量 θ0 = 1e-2、3 GeV、L = 0.5 m、bz = 3・9 1/m) と 1e-3 以内で一致 (比 1.00002〜1.00017、100 MV の bz = 3 で 0.9990)。**加速しても損失は変わらない** (横運動量はソレノイドの縦磁場で変わらず曲率 κ ∝ p⊥/p² なので γ⁴κ² は一定) の解析式。TRPT の加速そのもの (運動方程式の積分) は C2 の試験 (bz = 3 まで) |
| **POL + TRPT の加速** | track6 の経路 (POL + RAD) は TRPT の基準運動量の変化を扱わず、黙って δ = 0.27 の値を返していた (領域内の CAVI 2 台の輸送路)。**明示のエラーにした**: `SADError: POL (スピン追跡) + TRPT + RFSW の加速 (VOLT≠0 の CAVI・MULT・TCAVI) は未対応`。加速のある輸送路の追跡は POL なしの `TrackParticles` (領域内も対応)。track6 に基準運動量の更新を入れる実装は未着手 (領域の外の加速も同じ。スピンの歳差に energy-dependent な a·γ(s) と bz·p0/p(s) が要る) |
| 試験 | `test_sol_region_element_matches_lorentz_force_at_large_bz[3 × 5]`、`test_sol_region_radiation_loss_at_large_bz`、`test_sol_region_drift_radiation_loss_helix`、`test_sol_region_drift_radiation_k64_is_a_single_step`、`test_spin_in_sol_drift_at_large_bz_matches_thomas_bmt[3]`、`test_spin_in_sol_drift_at_large_bz_k64_is_coarse`、`test_emittance_pol_spin_tune_with_sol_region_equals_analytic[3 × 4]`、`test_rf_mult_acceleration_does_not_change_radiation_loss_in_sol_region[2 × 4]`、`test_pol_with_trpt_acceleration_is_an_explicit_error` |
| k64 モードの扱い | DRIFT の刻みは v2 だけ (k64 は K64 の 1 刻み) |
| 使えない・未確認 | (1) **ずれて傾いたソレノイドの Emittance[] の POL の量**は K64 の式のまま (SpinTune で 2.3e-4、n0 で 4e-3 rad の差。`spin.tradke_spin` の `xpzb`・dxpzb の修正が要る)。(2) POL + TRPT + 加速 (track6) の実装。(3) 領域内の CAVI の POL (Emittance[]) は VOLT=0 の CAVI だけ確かめた (加速する CAVI の POL の Emittance[] は RF あり領域の RFSW で、一周でエネルギーが元に戻るので解析式は縦の運動を含む)。(4) bz ≥ 3 の QUAD・MULT のスピンは 1e-4 台 (上)。(5) 大きなリングでの SOL 領域の確認 (大きなリング A に SOL を入れた確認は項目 4 の小リングのみ) |

## 要判断・残した項目

1. **WSPAC の Emittance[] が gs より 4% (x)・3% (y) 弱い (項目 3)**: ラスレットの式とは 0.1%。gs の 2〜4% の強さの原因が K64 のソースから見つからなかった (tintrb の式は同じ)。gs が正しいか pysad が正しいかの独立な根拠 (K64 の beamsize の定義など) が要る。
2. K1≠0 の BEND の端の模型 (四極が半径方向の面で始まる): CALC・追跡の共通の K64 の仕様。端面で始める模型に直すか (項目 2 の (1))。
3. ANGLE=0 の BEND・MULT の F1 (tsteee・tmulte の tblfre) の実装 (項目 1)。
4. POL + TRPT + 加速の track6 の実装 (項目 6)。
5. 4D の CALC の WSPAC は CALC6D と同じ 6D の経路 (K64 の qwspac の直訳ではない): K64 が NaN を返す機能なので、直訳の検証ができない。

## 変更したファイル

`pysad/emit6d.py` (F1: `_tblfre`・`_tbendal`・`_tbende_core`・`_bende_fast`・`_emit_circ`)、`pysad/touschek.py` (WSPAC、SOL)、`pysad/emit_extra.py` (`sc`)、`pysad/calc6d.py` (`sc`・SOL)、`pysad/ffscmd.py` (WSPAC の EMIT・Emittance[]・CALC・`_sc_opts`・DP0・CALC6D の RADTAPER・SOL の annotate・POL + TRPT のエラー・A1 の警告の削除)、`pysad/track6n.py` (`edge_exact1`・`k_bendi`)、`pysad/track6.py` (`_tbendi`・`_sol_drift_v2` の刻み)、`pysad/track6n_ext.py` (`k_sol_drift`)、試験 2 つ。**既存の試験の変更 3 件** (動作が変わったため): A1 の `test_emittance_warns_about_ignored_bend_f1` → `test_emittance_no_longer_warns_about_bend_f1`、`test_sol_optics.py::test_sol_calc6d_is_refused` → `test_sol_calc6d_runs`、`test_v2_C2_optics.py::test_radtaper_in_calc6d` の k64 の期待値 (CALC6D にテーパーが入る: −7.887e-5 → −2.1736e-7。gs は今のビルドで RADTAPER の CALC6D が動く)。

## 全試験・回帰

全試験 (`--ignore` の 3 本を除く、`PYSAD_SAD_PACKAGES` つき、他の担当の追加分を含む): **1131 通過・5 skip・失敗 0** (B1 の追加 87 件を含む)。回帰 (`QUICK=1`、k64): **50 ベンチ FAIL 0** (132 s)。最初の全試験の実行では 5 件が失敗した (動作が変わった既存の試験 3 件と、`beam_matrix_at` を空間電荷つきにしたことによる WSPAC の追跡の試験 3 件)。前者は期待値を更新し (上の「変更したファイル」)、後者は `beam_matrix_at` を元に戻して直した。
