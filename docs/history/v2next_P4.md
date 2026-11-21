# v2next の記録: 担当 P4 (SOL・空間電荷・DAPWIDTH・FREE の残り、2026-11-19)

対象: `pysad/track6.py` (領域内の QUAD の刻み・ψ・光子)、`pysad/solemit.py` (ψ の規則)、`pysad/spch.py`・新規 `pysad/spch_k64.py` (k64 の PSPAC)、`pysad/dapert.py` (DAPWIDTH の窓)、
`pysad/ffscmd.py`・`ffsextra.py`・`match6d.py` (成分ごとの FREE・DAPWIDTH の既定・WSPAC のビーム行列)。numba のカーネルは編集していない (`py -3.12 -m pysad.nbcache` は走らせた)。
P3 が編集中の `track6*.py` は、領域内の QUAD・MULT・光子の数行の追加だけ (`track6n_ext.py` は触っていない。実験で一時的に書き換え、`.bak` から戻した。P3 の差分はそのまま残っている)。
試験 (新規 6 ファイル 57 件): `test_v2_P4_free.py` (6)、`test_v2_P4_dapwidth.py` (11)、`test_v2_P4_pspac_k64.py` (14)、`test_v2_P4_wspac_beam.py` (4)、`test_v2_P4_sol_spin.py` (18)、`test_v2_P4_sol_photons.py` (4)。
既存の試験の変更 (動作が変わったための期待値の追従): `test_v2_B2_free.py::test_mixed_variables_fit_and_table_order_like_gs` の VAR の並び (何もしない変数 QF.3 が先頭に出る。gs と同じ順)、
`test_v2_C2_spin.py::test_spin_path_follows_the_same_orbit_as_fast_tracking` の許容 1e-10 → 5e-9 (領域内の QUAD は POL のとき細かく刻むので座標が 2e-9 だけ精密になる)。
例 (`example/v2_next/`): `P4_free_noop_variable.sad`・`P4_dapwidth_window.sad`・`P4_pspac_k64_shape.sad`・`P4_pspac_astrid2.sad`・`P4_wspac_tracking_beam.sad`・`P4_sol_photons.sad`
(各 `.gs.txt`・`.pysad.txt` (v2)・`.pysad_k64.txt`)、`P4_sol_spin_quad_mult.py` (`.out.txt`・`.out_k64.txt`)、説明 `P4_README.md`。
使ったリング: 12 セル FODO (0.2・0.5・3 GeV)、24 セル FODO (3 GeV、6 極つき)、大きなリング A (3 GeV、603 m。DAPWIDTH のみ)、ASTRID2 (`example/astrid/lattice/ASTRID2_v01_01.sad`、0.58 GeV、45.7 m)、5 MeV の輸送路。
大きなリング A (`kek-ls_v03_12.lat`、603 m) は DAPWIDTH だけに使い、PF リングは使っていない (下の「未確認」)。gs の追跡は NPARA = 1 で取った (gs の既定の NPARA は 1 より大きく、`TrackParticles`・`DynamicApertureSurvey` は
粒子を子プロセスに分け、PSPAC は各プロセスが自分の粒子だけで電位を作る。この記録の gs の値はすべて NPARA = 1)。

## 1. 一覧

| 項目 | 結論 | モード | 試験 |
|---|---|---|---|
| 1 SOL 領域内の QUAD・MULT・CAVI のスピン | QUAD の刻みを bz·al/ndiv ≤ 0.02 にした (9e-5 → 1e-6)。**端の縦磁場の積分 ψ は DISFRIN=0 のときだけ**にした (MULT の K0 つきで 1.6e-4 → 1e-6)。CAVI は VOLT = 0.2 MV の領域内の CAVI で Emittance[] と追跡の SpinTune が 1.6e-7 (弱い試験) | v2 だけ (k64 は K64 のまま) | `test_v2_P4_sol_spin.py` |
| 2 SOL 領域内の光子・加速 | **光子: v2 に DRIFT・QUAD・MULT を実装** (位置が軌道上、エネルギーの合計が古典値)。**加速 (POL + TRPT + 加速) は実装しない** (理由と条件は §3) | v2 だけ | `test_v2_P4_sol_photons.py` |
| 3 PSPAC の K64 の形状関数 | **k64 モードに実装** (`spch_k64.py`)。gs と rms の 12 桁以上一致 (5 条件、多周回・凍結の 4 条件、ASTRID2 の 3 SPCH 10 周) | k64 (v2 は変えない) | `test_v2_P4_pspac_k64.py` |
| 4 WSPAC のビーム行列を追跡に | **v2 のオプション `WSPACSC = 1;`: 空間電荷つきの自己無撞着な行列 (MINCOUP の下限つき) を追跡の WSPAC が使う**。チューンシフトが同じ行列の一次の式と 2e-4。**既定は変えない** (gs の追跡のビーム行列には MINCOUP の下限が入らず、ゼロ電流の行列と 1〜2% で合う) | v2 だけ | `test_v2_P4_wspac_beam.py` |
| 5 DAPWIDTH の窓 | **trackd.f のスロット管理を直訳** (`dapert.py`)、DAPWIDTH の既定を 7 に (gs と同じ)。NORAD の地図・Score と、RAD + FLUC の最終 seed が gs と一致 | 両モード (K64 の動作の再現) | `test_v2_P4_dapwidth.py` |
| 6 要素ごとの FREE の何もしない変数 | **gs と全桁一致** (QF.3 の変数が結合先 QF.1 を持つ何もしない変数、変数の数・VAR の行に数える、エラーの規則、Saved の欄、RESET) | 両モード | `test_v2_P4_free.py` |

## 2. 項目ごとの記録

### 2.1 項目 6: 要素ごとの FREE の何もしない変数

| 欄 | 内容 |
|---|---|
| 状況 | B2 は「先頭でない成分 (`FREE QF.3 K1`) の変数は pysad は作らず Info だけ」で、gs と反復の行の最後の数 (変数の数) と VAR の行が 1 つ違った。gs で規則を調べ直した |
| K64 の動作 (gs で確認) | (a) `FREE QF.3 K1` は VAR に `QF.3 K1 .7 ! … QF.1 1.0` の行を作る (結合先 QF.1)。その変数だけの `GO` は `?No variable.` を出さず、反復なしで残差の行だけ (0.7344)。(b) `FREE QF.1 K1` の変数のラベルは QF.1 (`FREE QF K1` は QF)。(c) `FREE QF K1` のあとの `FREE QF.3 K1` は `Element already used as variable: QF`、`FREE QF.3 K1` のあとの `FREE QF K1` は `A component has been already used as variable: QF.3`。Info (`QF.3 is coupled to QF.1`) は変数を作るときだけ。(d) `FIX QF.1 K1`・`FIX QF.3 K1` はその成分だけ、`FIX QF` (要素全体の語) は成分の変数に触れず、`FIX QF.*` は成分も消す。(e) VAR の Saved の欄は設計値 (定義の値、`SAVE` で更新)。GO では変わらず、`RESET` で Now が設計値に戻る。(f) `FREE QF.3 K1; FREE QF.1 K1; GO` は変数 2 つ (反復の行の最後の数 2)。QF.3 の Now は QF.1 に追従 |
| v2 の動作 | k64 と同じ (どちらのモードも K64 の動作の再現)。`ffscmd._free_components` が既定のキーの成分の変数を `comp`・`compdef` つきで作る (先頭は要素全体の値、先頭でないものは `noop`: `set_var` は何もしない、微分の列は 0)。`_tfvcomp` は `compdef` を数えない (既定でないキーの成分の変数だけが NOKEEPEXP の展開を止める)。`ffsextra.cmd_variables`: 結合先の欄・Saved の欄 (`_saved_value`)。`reset_cmd` は変数の Now も戻す |
| 根拠 | gs の全桁一致: `FREE QF.3 K1; FREE QF.1 K1` で NX を 1.30 にする GO が残差 5.3746E-02・5.3438E-02、QF.1 = QF.3 = .743599255507 (12 桁)、反復の行の最後の数 2。VAR の表 (並び・Saved・Previous) が gs と一致 (`P4_free_noop_variable.gs.txt` と `.pysad.txt`) |
| 試験 | `test_v2_P4_free.py` (6 件: 何もしない変数・先頭と先頭でない変数のフィット・競合のエラー・FIX・Saved と RESET・NOKEEPEXP の規則) |
| k64 モードの扱い | 同じ (ゲートなし。K64 の動作) |
| 使える範囲と条件 | 3 GeV の周長 18.4 m の FODO (QF・QD 各 4 台、BEND 8 台)。大きなリングでは B2 が既定でないキーの成分の FREE を大きなリング A で gs と確かめたが、今回の何もしない変数は小さなリングだけ |
| 使えない・未確認 | `COUPLE` で結合した要素の成分、`ElementValues` との連動、`BYE ALL`。反復の途中の `qcod`・`qtwiss` のオーバーフロー行の数が gs と違う (カオス的なフィットの途中。値には影響しない) |

### 2.2 項目 5: DAPWIDTH の窓

| 欄 | 内容 |
|---|---|
| 状況 | pysad の `DynamicApertureSurvey` は常に全 51 点を追い (DAPWIDTH の既定 51)、`Maximum number of particles` も違った (gs の既定は 7: 3 行で 42、pysad 306)。窓の動きで追跡する点が変わると RAD + FLUC の乱数の列が変わる |
| K64 の動作 (trackd.f) | 行 i ごとに損失が見つかった最小の j (nxm、初期 52) から下へ ncons = max(1, nint(DAPWIDTH)) 点の窓の中だけを追跡 (`j < nxm − ncons` で打ち切り、損失点が見つかると nxm = j)。同時に追うスロットは npmax = n1p·ncons (NPARA = 1)。損失・完了でスロットが空くと配置ループに戻り、空きがあるとき (np1 < npmax) は 16 周ごとに戻る。空きスロットは (i, j) の順で点を探し、なければ末尾の粒子と入れ替えて np1 を減らす。窓の外の点は ntloss = 2^29 のまま (地図は `*`) |
| v2 の動作 | `dapert.dynamic_aperture` に直訳 (配列位置 ↔ スロットの ID の入れ替え、配置ループ、nw = 16、`track6n.turn` の `kid`)。`NP = min(4000, n1·min(DAPWIDTH, 51)·2)` の表示。DAPWIDTH の既定を 7 に。**どちらのモードも同じ** (K64 の動作の再現) |
| 根拠 (gs) | 24 セルの FODO (3 GeV、SF = 4・SD = −7 の6 極、RF 4 MV)、`{{1,300},{1,1},{-0.5,0,0.5}}`、100 周、NORAD: DAPWIDTH = 既定 (7)・3・12・51 の 4 通りで Score 24 と地図 (窓の外 `*` を含む) が一致 (損失の境界の 1 点 j = 14 の 'H' と 'h'、46 周と 50 周は 1 周ほどの差: カオス的な境界の粒子)。RAD + FLUC (NPARA = 1) の最終 seed: DAPWIDTH = 2、z = {−0.5, 0, 0.5, 1}、x = {1, 150}、5・10・20・25 周の 4 例が一致 (15 周だけ境界の粒子で違う)、DAPWIDTH = 1 (スロット 1 つ) の xm = 20〜130 が一致 (xm = 140 だけ境界の粒子で違う)、DAPWIDTH = 5・6 の 1〜4 周が一致 |
| 試験 | `test_v2_P4_dapwidth.py` (11 件: 地図 4・窓の動きの seed 4・スロット 1 つの seed・ResultOfDAPERT の 2^29・大きなリング A の地図) |
| k64 モードの扱い | 同じ。k64 は gs と RAD + FLUC の seed まで一致する (上)。**v2 は放射の式 (B111・B112 など) が K64 と違うので RAD + FLUC の seed は gs と合わない** (NORAD の地図は合う) |
| 大きなリング A | `kek-ls_v03_12.lat` (3 GeV、周長 603 m、3563 要素、RFSW + NORAD)、`{{1,400},{1,1},{-0.5,0,0.5}}`、60 周: DAPWIDTH = 7 (j = 1〜14 が窓の外) と 51 の地図・Score 63 が gs と一致 (境界の 1 点 (z = 0.5 の j = 34) だけ 1〜2 周違う)。gs 11.7 s、pysad 4.9 s。`P4_dapwidth_ringA.out.txt` |
| 使える範囲と条件 | NPARA = 1、軸 `Axes->"ZX"`、`Output->6` か 0、損失の判定は LOSSAMPL・LOSSDZ。最大 n1p = 201、DAPWIDTH 任意 |
| 使えない・未確認 | **NPARA > 1** (gs は子プロセスに分け、子ごとに乱数の種をずらす。再現しない)、`Damping->True`・`DampEnough`。損失の境界の粒子は gs と 1〜数周違うことがある (追跡の丸め。NORAD で 1 点、RAD + FLUC で境界の違いが乱数の列を変えて最終 seed が違う例がある)。大きなリングでの時間: 窓のスロット管理は Python のループ (点の数 × 行の数の走査)。1 行 51 点・100 行で数秒 |

### 2.3 項目 3: PSPAC の K64 の形状関数 (k64 モード)

| 欄 | 内容 |
|---|---|
| 状況 | v2 の PSPAC は K64 と同じ物理モデルを別の形状関数 (2 次の B スプライン) で書き直したもの。粗いメッシュの gs との差が残っていた (A1)。k64 モードは追跡に空間電荷を入れていなかった |
| K64 の動作 (spch.c・psn.c・tspch_.c の読みと gs) | 電荷: 横は CIC (双線形、|x| ≤ dx(nx/2 − 1) の外は入れない)、縦は「ほぼ最近接」(`rhocal`: セルの前半は 1 を下のノードへ、後半は 2(1 − dzic)² を下に)。電位: 5 点差分 φ(i±1, j) + (dx/dy)² φ(i, j±1) − 2(1 + (dx/dy)²) φ = −2π(dx/dy)ρ (∇²φ = −2π ρ/(dx dy)、境界 φ = 0。FACR は離散方程式の厳密解) を双 3 次補間 (`makebcu`: 中心差分、NR の bcucof)。線密度は 3 次スプライン (`splined`、端の傾き 0、ノード数 nz)。キック: Kfac = K L/(np² dz)、K = 2 N r_0/(β²γ³)、`px_canon −= Kfac λ ∂φ/∂x`、`py` も、**δ −= Kfac λ' φ** (δ を直接変える。z は変えない)、dv = −δ(2 + δ)/h1/(h1 + (1 + δ) h0)。K は最初の呼び出しの h0 で固定、SPCH ごとの電位 (SP[isp]) は PSPACNTURNCALC ターンを超えると凍結、λ (rhoz) は大域変数で最後に計算した SPCH のもの |
| v2 の動作 | k64 モードに `spch_k64.K64Pspac` (`SpaceCharge._kick2d` から呼ぶ)。k64 は PSPAC だけ入れる (SPAC は K64 の 1/γ の誤り (B89)・WSPAC の追跡は未実装で今までどおり入れない)。v2 は変えない |
| 根拠 (gs と全桁) | 5 MeV・1e11 個・半径 1 mm の一様な円 (準乱数 4096 個)、`D(0.5) SPCH D(0.5) SPCH`、NPARA = 1。rms x・y・px・py・平均 δ・rms δ が 5 条件 (縦 20・4・2 セルの 3 つ、DZ = 1 mm、64×64・16 ノード) で相対 3e-12 以内 (14 桁)。**多周回**: PSPACNTURN = 3 (毎周再計算)・PSPACNTURNCALC = 1 (凍結)・PSPACNTURN = 4 かつ PSPACNTURNCALC = 2・PSPACNTURN = 1 で 2 周 (gs は解放済みのメモリの中身 = 凍結した電位を読む) で全桁一致。**大きなリング (ASTRID2 に SPCH を 3 か所、PBUNCH = 1e14、10 周)**: rms x = .00229757974981 mm が gs と 14 桁、rms δ・平均 δ も 11 桁以上一致 (v2 は x が 8% 違う)。独立な確認: spch.c をそのままコンパイルした C の関数と、同じ入力で電位 1e-9 (相対 1.5e-13)・δ のキック 1e-16。内部の部品 (splined・ポアソン・双 3 次の係数) は 2 次式・差分方程式で確認 |
| v2 との差 (同じ条件) | rms x の gs との差: 縦 20 セル 0.6%、4 セル 4.5%、2 セル 12% (v2 は粗いメッシュで gs・解析式から離れる。A1 の記録)。平均 δ は gs が 1.8e-5、v2 が 1.9e-7 (下の K64 の候補) |
| 試験 | `test_v2_P4_pspac_k64.py` (14 件: 5 条件・多周回 4 条件・v2 との差の記録 3・PSPACNZ < 3 のメッセージ・部品の確認) |
| k64 モードの扱い | 新しい機能 (k64 だけ。gs と全桁一致)。既存の k64 の試験は PSPAC を使わないので変わらない |
| 使える範囲と条件 | 電子 (r_0 = r_e m_e/mass の係数は K64 のまま。**電荷の 2 乗は入れない**: K64 は q を使わない)。RFSW の ON/OFF は問わない。SPCH の位置の間隔 (前の SPCH の中点から) は v2 と同じ。PSPACNZ ≥ 3 (ノード nz + 1)、nx・ny は 4 以上 (FACR は 2 のべきだが厳密解なので一般の値も同じ解) |
| 使えない・未確認 | **gs の NPARA > 1 の再現** (子プロセスごとの電位。NPARA = 1 で比べる)。gs は 1 セッションの 2 回目の `TrackParticles` で PSPAC が働かない (SPC8) ので、k64 は TrackParticles ごとに新しく始める (gs と違う)。PSPACNZ = 1 (既定) は gs が SIGSEGV (B93) で、k64 は PSPACNZ < 3 をエラーにする。メッシュの端のセル (ii = nx − 1 または jj = ny − 1) の bcucof は配列の外を読み、縦の kk = nz − 1 のスプライン係数は未初期化 (k64 は外側 0)。SPAC・WSPAC との併用、ウェーク・放射との組み合わせ、PF・大きなリング A |

### 2.4 項目 4: 追跡の WSPAC のビーム行列

| 欄 | 内容 |
|---|---|
| 状況 | B1: 追跡の WSPAC のビーム行列は空間電荷なし (ゼロ電流) の平衡値・MINCOUP なし。垂直のエミッタンスが 0 の平らなリングは「WSPACBEAM を与えること」のエラー。「追跡に空間電荷つきの行列を入れる」を調べた |
| K64 の動作 (gs で確認、`P4_wspac_tracking_beam.gs.txt`) | tturn の `twspac` は `beamsize(:,l)` (**`CODPLOT; WSPAC; EMIT;`** の最後の一周のビーム行列) を使う (`Emittance[]` のあとの追跡は入力を返すだけで動かない)。結合のある 0.5 GeV のリング (C2s_wspac と同じ、2e10 個) の試験粒子のチューンシフトは gs で (−2.938e-5、−1.210e-3)。**MINCOUP の下限を入れた空間電荷つきの行列 (−2.776e-5、−3.807e-4) とは合わず**、ゼロ電流の行列 (v2 の既定: −2.982e-5、−1.201e-3: x 1.5%・y 0.7%) と下限なしの空間電荷つき行列 (−2.982e-5、−1.200e-3) と合う。垂直のエミッタンスが 0 の平らな 12 セルのリングは gs が NaN (ビームの高さ 0)。すなわち K64 の追跡のビーム行列には MINCOUP の下限が入らない |
| v2 の動作 | **既定は変えない** (ゼロ電流の行列。gs と 1〜2%)。**オプション `WSPACSC = 1;`** (v2、K64 にはない): `ffs.beam_matrix_at` が `emittance_rules(sc=…)` (`touschek_table(wspac=True)` の最終パス = `BeamMatricesSC`、MINCOUP の下限つき) を返す。強い空間電荷でビームが膨らむ場合と、垂直のエミッタンスが 0 のリングで WSPACBEAM なしに使う場合のため。キャッシュのキーに PBUNCH・MINCOUP を足した。k64 は変えない |
| 根拠 | (1) WSPACSC = 1 の行列の εy の下限: EMITYC/EMITXC = 0.5 (MINCOUP = 0.5)、Σyy(1)/βy(1) = EMITYC (PBUNCH = 1e8 で 0.5%、1e6 で 2e-5)。(2) 試験粒子の追跡のチューンシフト (PBUNCH = 1e6、0.2 GeV の 12 セル FODO): 追跡 (+4.1102e-5、+8.4590e-5) が、同じ行列を各要素の入口に置いた薄いレンズ (長さは前の要素の中点から) の一次の式 Δν = −(1/4π) Σ β k と比 0.99980・0.99976 |
| 試験 | `test_v2_P4_wspac_beam.py` (4 件: WSPACSC の行列の εy の下限・チューンシフトの一次の式・既定はゼロ電流の行列・k64 は入れない) |
| k64 モードの扱い | 変えない (追跡に WSPAC を入れない) |
| 使える範囲と条件 | WSPACSC = 1: 0.2 GeV の 12 セル FODO、PBUNCH = 1e6〜1e8、MINCOUP = 0.5、RF あり。強-弱 (ビームの大きさは固定) |
| 使えない・未確認 | WSPACSC = 1 の行列は gs の追跡と合わない (MINCOUP の下限が入るため)。強い空間電荷 (0.2 GeV で PBUNCH ≳ 1e9) は反復が収束しない、RF なし、大きなリング (反復の収束・時間)。`LINE["SIG…"]` も WSPACSC = 1 のとき同じ行列 |

### 2.5 項目 1: SOL 領域内の QUAD・MULT・CAVI のスピン

| 欄 | 内容 |
|---|---|
| 状況 | B1・E: bz ≥ 3 1/m の QUAD・MULT のスピンは 1e-4 台。CAVI の POL (Emittance[]) は VOLT = 0 だけ確認 |
| 調べたこと | Thomas-BMT の積分 (`test_v2_B1_sol._tbmt`) と比べると、(a) QUAD は bz = 4.5・9 で 2e-5・9e-5 (刻み bz·al/ndiv ≈ 0.2 の弦の誤差)、(b) **MULT (K0 + K1) は bz に依らず 2e-5〜1.6e-4**: 刻みを 16 倍にしても変わらず、端の縦磁場の積分 ψ の双極成分 k0 y (5e-5) を落とすと 1e-6 台になる。K64 の式も v2 (E) も ψ を DISFRIN に関係なく足していた |
| K64 の動作 | 領域内の QUAD のスピンは `bsi` の式で最初・最後のスライスを半分にして全体で 1 スライス分足りず、Thomas-BMT の積分と 6e-2〜0.19 ずれる (bz = 0.3〜9 1/m。C2 が v2 を `bpar` に直した理由)。MULT は k64 では領域内が未対応 |
| v2 の動作 | (a) `track6._tsolqu_bz_ens`: POL のとき bz·al/ndiv ≤ 0.02 に刻む (Emittance[] の `solemit._tquade` と同じ)。(b) ψ (QUAD の k x y、MULT の k0 y + k1 x y + …) は端のフリンジを入れるとき (DISFRIN=0) だけ (追跡の `_RadHook`・`_tsolqu_bz_ens` と Emittance[] の `solemit._mult_psi`・`_tquade`)。DISFRIN=0 の動作は変えない |
| 根拠 | Thomas-BMT の積分 (0.4 m、DISFRIN=1、軸上と軸を外れた粒子・3 つのスピンの向き): QUAD・MULT とも bz = 0.3・3・4.5・9 1/m で **9.3e-7 以下** (以前: QUAD 9.0e-5、MULT 1.6e-4)。スピン軌道結合 D (x・px・y・py の微分): 追跡は 1e-8 以内 (以前 1e-6)、Emittance[] は 2e-6 以内、軸上の回転 R は 1e-6 以内 (以前の許容 3e-5)。リングの確認 (ASTRID2、RF の前に BZ = 0.5 T・0.9 m の領域に QUAD・MULT): Emittance[] の SpinTune が追跡の 1 周のスピン行列と 1.2e-7 以内。CAVI: 12 セルのリング (0.5 GeV) の領域内の CAVI (VOLT = 0 と 0.2 MV) で 1.6e-7 (弱い試験: 位相 0 なので SpinTune がエネルギー変化にほぼ鈍い) |
| 試験 | `test_v2_P4_sol_spin.py` (QUAD・MULT × 4 条件の Thomas-BMT 8、k64 の QUAD の記録 1、D・R の 4、CAVI 2、ASTRID2 の 3 = 18 件) |
| k64 モードの扱い | 変えない (QUAD は K64 の式のまま) |
| 使える範囲と条件 | 0.1〜1 GeV、bz = 0.3〜9 1/m、長さ 0.4〜0.5 m の QUAD・MULT・DRIFT・BEND の領域 (DISFRIN=1)、ASTRID2 までのリング |
| 使えない・未確認 | **DISFRIN=0 (端のフリンジあり) の ψ を、硬い端の B_s と軌道のフリンジ写像を合わせた厳密な基準で確かめていない** (DISFRIN=1 だけ)。ずれ・ROTATE を持つ MULT の ψ、CAVI の加速・位相が 0 でないときの SpinTune の確認 (弱い)、大きなリング A・PF リング・8 GeV (使っていない) |

### 2.6 項目 2: SOL 領域内の光子・加速

| 欄 | 内容 |
|---|---|
| 状況 | docs/46: 光子は SOL 領域内は記録しない。領域内の加速 (POL + TRPT) は SADError |
| K64 の動作 (gs で確認、`P4_sol_photons.gs.txt`) | 3 GeV、BZ = 1 T の領域の中の MULT (K0 = 0.01、L = 0.4)・DRIFT の光子を**記録する** (MULT: 955 個、gx は 0.5006〜0.8999 m = 要素の範囲、エネルギーの合計 312.5 eV/粒子 = 古典値 285 eV の 1.10 倍 (統計の範囲))。**領域内の QUAD の光子は記録しない** (0 個。B68 の疑いと同じ: RAD ON・bz≠0 の領域内の四極の放射がない)。領域内の BEND・MULT の偏向はほぼ無視 (B69) |
| v2 の動作 | `track6._tsol_region_v2` が領域内の要素ごとに `l_track` を進め、DRIFT・QUAD・MULT で `_ph_set` (領域の外の同じ要素と同じ位置の情報)・刻みごとに `_ph_adv`。BEND は未対応 (記録しない)。k64 は変えない (MULT は未対応のまま、DRIFT の光子は領域の始点 gx = 0 で記録する今までの動作) |
| 根拠 | (1) 光子の放射点が粒子の軌道上: MULT の K0 = 0.01 の放物線 x = (K0/2L) s² と 5e-7 以内 (gx = 進行方向、gy = x、gz = y の SAD の大域座標)、DRIFT・QUAD の光子は要素の範囲。(2) エネルギーの合計: MULT で 309 eV/粒子 (1500 粒子、古典値 285 eV。統計の誤差 8%)、領域の外の同じ並び (bz = 0) の MULT と同じ精度。(3) RAD NOFLUC の損失は領域内と領域外 (bz = 1e-9) で 1e-17 まで一致 (光子の数 (FLUC) は乱数の使い方が違い一致しない) |
| 試験 | `test_v2_P4_sol_photons.py` (4 件) |
| k64 モードの扱い | 変えない |
| 使える範囲と条件 | 3 GeV、bz = 0.1 1/m、傾き・ずれのない SOL の領域の DRIFT・QUAD・MULT (ANGLE = 0) |
| 使えない・未確認 | 領域内の BEND・SEXT・OCT の光子、傾いた・ずれた SOL の領域の光子の位置 (`geo_frac` は要素の幾何を使うので、領域の軸が傾くと別の系になりうる)、偏光 (ξ) の値、gs との位置・数の比較 (gs は QUAD を記録せず MULT の偏向も無視するので比べられない)、大きなリング |
| **加速 (POL + TRPT + 加速) を領域内で実装しない理由と条件** | track6 のスピンの経路は TRPT の基準運動量の変化を扱わない (領域の外の CAVI は D1 が `tcav_acc` で扱う)。領域内では (i) 刻みごとに基準運動量 p0(s) を進め、(ii) ソレノイドの bz は入口の p0 で規格化した値なので bz p0/p(s) に直し、(iii) スピンの歳差の aγ(s)・運動量の向きの磁場の積分を p(s) で作り直す必要がある。独立な基準 (加速する空洞の中のソレノイド + RF の Thomas-BMT の数値積分: 縦の運動も含む) は作っていない。**実装するなら**: `sol_mult_full` (加速の重み ws) の刻みで `_RadHook` が基準運動量を更新し、`_tradk_sol` に p(s) を渡す。独立な基準は RF の電場を含む運動方程式の積分 (E_z(t) cos(ωt))。現状は明示のエラーのまま (`test_v2_B1_sol.py::test_pol_with_trpt_acceleration_is_an_explicit_error`)。リングの RF (TRPT でない) の領域内の CAVI は上の通り Emittance[] と追跡が合う |

## 3. K64 の新しい誤りの候補 (B 番号は親が付ける)

| 候補 | 内容 | 根拠 |
|---|---|---|
| K64-P4-a | **SOL 領域内の端の縦磁場の積分 ψ が DISFRIN に関係なく入る** (MULT の双極成分 k0 y、QUAD の k x y)。フリンジなしの模型でスピンが 1e-4 ずれる (k64 の領域内の MULT は未対応なので、v2 の式 (E) が受け継いだもの。v2 は DISFRIN=1 で入れない) | Thomas-BMT の積分 (§2.5) |
| K64-P4-b | **PSPAC のキックは平均のエネルギーを保存しない**: 縦の 3 次スプライン λ' と形状関数の不整合で、一様な円・長さ 0.1 m のバンチ (4096 個) の平均 δ が +1.8e-5 (rms δ 3.6e-4 の 5%)。v2 は 1.9e-7。粗い縦メッシュ (4 セル) では平均 δ が rms の 12% | `test_v2_P4_pspac_k64.py` の表 |
| K64-P4-c | PSPAC: (1) 電位・λ のメモリを PSPACNTURN 周で解放するのに、その後の呼び出しが解放済みのメモリを読む (値は凍結した電位と同じで壊れないが未定義の動作)、(2) メッシュの端のセルの `makebcu` が配列の外を読む、(3) 縦のスプライン係数 `rhoz_sp[nz]` が未初期化、(4) 1 セッションの 2 回目の `TrackParticles` で PSPAC が働かない (SPC8、既出)、(5) `PSPACNZ = 1` で SIGSEGV (B93、既出)、(6) **NPARA > 1 の `TrackParticles` は粒子を子プロセスに分け、PSPAC は各プロセスが自分の粒子だけで電位を作る** (粒子数が 1/NPARA で、キックの係数 Kfac も子ごとの np で変わる。統計の違いだけで物理的な誤りではないが、NPARA に依存して PSPAC の結果が変わる。A1 の gs の 16% の x・y の非対称はこれに由来する可能性がある: NPARA = 1 に替えると gs の値が変わる (rms x 1.23379 → 1.23401 mm)。再現は未確認) | spch.c・tspch_.c の読み、gs |
| K64-P4-d | `DynamicApertureSurvey` の DAPWIDTH の既定 (7) と trackd の窓は文書化されているが、pysad が全点を追っていた (これは pysad の不具合で、直した) | 項目 5 |
| K64-P4-e | gs の `TrackParticles` は WSPAC ON でも `Emittance[]` のあとでは入力をそのまま返す (`CODPLOT; WSPAC; EMIT;` が要る。beamsize が EMIT の plot でだけ作られるため)。垂直のエミッタンスが 0 の平らなリングは NaN (ビームの高さ 0 のまま追跡する) | `P4_wspac_tracking_beam.gs.txt` |
| K64-P4-f | gs: RAD ON・bz≠0 の領域内の QUAD は光子を記録しない (B68 の延長の疑い) | `P4_sol_photons.gs.txt` |

## 4. 要判断・残した項目

1. **DAPWIDTH の既定を 7 にした** (gs と同じ。以前は 51)。51 を前提にした台本 (DAPWIDTH を書かないもの) は窓の外が `*` になり、Score が変わる (gs の動作に合う)。全点を追うには `DAPWIDTH = 51;`。
2. **ψ を DISFRIN=1 で入れない** (v2): 1e-4 の変化。DISFRIN=0 (フリンジあり) の ψ が物理的に正しいかは、硬い端の B_s と軌道のフリンジ写像を合わせた厳密な基準が要る (未着手)。
3. **`WSPACSC` を既定 0 にした** (gs の追跡のビーム行列に MINCOUP の下限が入らないため、ゼロ電流の行列が gs と 1〜2% で合う)。WSPACSC = 1 は K64 にないオプションで、gs と合わない (下限つき)。既定を空間電荷つきの行列に変える (K64 が `EMIT` の beamsize を使う意図) かは判断が要る。
4. 領域内の加速 (POL + TRPT)・領域内の BEND の光子・DISFRIN=0 の ψ の独立な確認は実装しない (§2.5・§2.6、理由と条件を記載)。
5. **大きなリング A での確認は DAPWIDTH だけ** (地図・Score が gs と一致)。PSPAC (3 SPCH・10 周)・SOL 領域内の QUAD・MULT のスピンは ASTRID2 (0.58 GeV、45.7 m) で代えた。603 m のリングでの PSPAC・領域内の光子・WSPAC の追跡は未確認 (DAPWIDTH のスロット管理は Python のループで、100 行 × 51 点・100 周で 5 秒)。PF リングは使っていない。

## 5. 変更したファイル

`pysad/ffscmd.py` (成分ごとの FREE の何もしない変数・`reset_cmd`・DAPWIDTH の既定 7 と `dynamic_aperture` への引数・`beam_matrix_at` の WSPACSC)、`pysad/ffsextra.py` (VAR の Saved・結合先)、`pysad/match6d.py` (何もしない変数の微分)、
`pysad/dapert.py` (trackd のスロット管理)、`pysad/spch.py`・新規 `pysad/spch_k64.py` (k64 の PSPAC)、`pysad/track6.py` (領域内の QUAD の刻み・ψ・光子)、`pysad/solemit.py` (ψ)、
試験 (新規 6 ファイル 57 件、変更 2 件)、例 (`example/v2_next/P4_*`)。コミットしていない。

## 6. 全試験・回帰 (2026-11-19、`PYSAD_SAD_PACKAGES` つき、k64 が既定)

- 全試験 (`--ignore` の旧式 3 本を除く): **1568 通過・8 skip・1 失敗**。失敗は `test_v2_C_measure.py::test_k64_mode_measure_and_trcod_are_ignored_as_before` (`measure.py` の TRCOD のメッセージの文言が英語に変わったため。別の担当の変更で P4 とは無関係。P4 の変更の前後は未比較だが、メッセージの出力は P4 が触っていない)。
- 回帰 `QUICK=1 run_regression.py`: **50 ベンチ FAIL 0** (128 s)。
- P4 の新規 57 件は全部通過。変えた既存の試験は `test_v2_B2_free.py` 1 件 (VAR の並び)、`test_v2_C2_spin.py` 1 件 (許容)。WSPACSC を既定 0 にしたので `test_v2_C2_spch.py` の WSPAC の試験は変えていない。

## 7. 付けた B 番号 (v2.2.0 の文書化で。docs/21)

- §3 の K64-P4-a (領域内の ψ が DISFRIN に関係なく入る) = **B146**、K64-P4-b (PSPAC の平均エネルギー) = **B147**、K64-P4-c (PSPAC の C コードの不整合 (1)〜(6)) = **B148**、K64-P4-e (TrackParticles の WSPAC の入力) = **B149**、K64-P4-f (領域内の QUAD の光子) = **B150**。
- K64-P4-d (DAPWIDTH の窓を pysad が追っていなかった) は pysad 自身の不具合で B 番号を付けない (docs/47 の (c))。
- 要判断 (DAPWIDTH の既定、`WSPACSC`、ψ) は docs/47 の v2.2.0 の (d)。
