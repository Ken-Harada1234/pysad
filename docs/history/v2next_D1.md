# v2next の記録: 担当 D1 (物理の残り: NORFSW の軌道、WSPAC、SynchroBeta、POL の加速、光子、spdepol、不安定なリングの POL、2026-11-18)

対象: `pysad/emit6d.py` (`closed_orbit6d`・`emittance6d`・`_emittance6d`)、`pysad/touschek.py` (`_twspace`・`_tintrb`)、`pysad/track6.py` (光子の記録、BEAMBEAM、TRPT の加速・分割)、`pysad/ffscmd.py` (光子の経路、SB の警告)、`pysad/spin.py`。
試験は `pysad/tests/test_v2_D1_emit.py` (22 件、うち `PYSAD_SLOW=1` が 1 件)。既存の試験の変更は 3 件 (§10)。
走る例 (`example/v2_next/`、gs の出力 `.gs.txt`、pysad の v2・k64 の出力つき): `D1_norfsw_orbit.sad` (項目 1)、`D1_photons_mult.sad`・`D1_photons_bb_bend.sad` (項目 9)、`D1_spin_accel.py` (項目 5、`.out.txt`)。numba のカーネルは編集していない (キャッシュの掃除は不要)。
gs の確認は WSL の gs (10-19 ビルド)。**gs の光学の中身を見るために、ソースの複製 (`~/sadd1`、元のフォルダは書き換えていない) の `tintrb.f` に出力文を足して再ビルドしたもの**を使った (項目 2。ここの数値は元の gs と同じ)。
使ったリング: 12 セル FODO (3 GeV・0.2 GeV、周長 45.6 m)、大きなリング A (`pysad/tests/kek-ls_v03_12.lat`、3 GeV、周長 603 m、3563 要素)、SP8 (大きなリング B、8 GeV、周長 1436 m、1998 要素)、PF リング (`example/pf/lattice/v9_16_asym.sad`)、ASTRID2 (`example/astrid/lattice/ASTRID2_v01_01.sad`)。

## 1. 一覧

| 項目 | 結論 | モード | 試験 |
|---|---|---|---|
| (1) NORFSW の Emittance[] の閉軌道のずれ | **pysad の不具合を直した**: z の戻りを δ で閉じていた (K64 は閉じない)。gs と全桁一致 | 両モード (gate なし) | `test_norfsw_*` 3 件 |
| (2) WSPAC の 4% (x)・3% (y) | **K64 の誤り**: 最後のスライスで `gammab(nlat+1) = 0` を読み、そのレンズが 8 倍強い。k64 は再現、v2 は物理的な値 (ラスレットと一致) | v2 だけ (k64 は K64 のまま) | `test_wspac_*` 4 件 |
| (3) SynchroBetaEmittance (大きなリング A) | **K64 の SB は大きなリング A で信頼できない** (独立な 6x6 の平衡・Emittance[]・追跡と 40 倍違い、モード数で動き続ける)。原因は未特定。v2 は収束しない値を警告 | v2 だけ (警告) | 警告 2 件 (1 件は SLOW) |
| (4) 傾いてずれたソレノイドの POL の量 | **未着手** (§7) | — | — |
| (5) POL + TRPT + 加速 | **実装した** (CAVI、SOL 領域の外)。断熱加速の解析式と 5e-7 | v2 だけ | `test_spin_tracking_with_trpt_*` 4 件 |
| (6) ANGLE=0 の BEND・MULT の F1 | **未着手** (K64 の tsteee を読んだ内容だけ §7 に) | — | — |
| (7) 側帯の項と追跡の 2.3 倍 | **未着手** (追跡の分解は CPU を使うので実行していない) | — | — |
| (8) SP8 (大きなリング B) の POL | **原因を特定して直した**: 光学が不安定 (閉軌道が発散して NaN) なのに、NaN の 1 周 (POL では 21〜65 秒) を 30 回以上繰り返していた | 両モード (不具合) | `test_pol_emittance_of_unstable_*` |
| (9) 光子 (MULT・SEXT・OCT・ドリフト) と BEAMBEAM + BEND | **実装した**。gs と全列一致 (MULT・QUAD・BEND の 325 個、BEAMBEAM + BEND + TRPT の 71 個)、SEXT・OCT は解析式。SOL 領域内は未対応 | 両モード (追加) | `test_photons_*` 5 件 |
| (10) `spdepol` の誤りの位置 | **位置を特定**: K64 の脱偏極率は emittance/減衰率に比例 (1/(1−e^{−2λ}) の二重数え)。正しい率は拡散で決まり減衰率に依らない | v2 は独立な式 (B3) | `test_k64_spdepol_*` |
| (11) `test_b107_bend_k1_with_vertical_dispersion_v2` の失敗 | **再現しなかった** (§9) | — | — |

## 2. 項目 1: NORFSW; RAD の Emittance[] と閉軌道のずれ

| 欄 | 内容 |
|---|---|
| 状況 | 12 セル FODO (3 GeV、周長 45.6 m、RF 3 MV)・QF の DX = 1e-4: gs の εx 2.70000647e-7 に対し k64 は 2.65085139e-7 (−1.8%)、チューンも −0.2932001 対 −0.2921601 (1e-3)。軌道のずれがなければ一致。RFSW ON でも一致。v2.0.0 から |
| 調べたこと | 要素ごとの比較はせず、gs の `OrbitAtExit` を見た: gs は x = 2.156e-4、px = 5.8e-5、**δ = 0、z = −5.98e-4** (z は閉じない)。pysad (k64) は x = −6.6e-5、**δ = −4.0e-4** (z を 0 に閉じる δ を探していた)。`tcod.f` (K64): RFSW OFF は `im = 5`、`dcod1(5) = 0`、`trs(:,5:6) = 0` で **x・px・y・py だけを閉じ、δ も z も動かさない** (`codw` も 0.01 倍)。RFSW ON でも空洞がない (`vceff = 0`) と `cod(5)` は決まらず z は閉じない |
| 原因 | `emit6d.closed_orbit6d` の最後の分岐 (RFSW OFF か、空洞のない RFSW ON) が、「SAD の Emittance[] は z も閉じる」というコメントのもと、一周で z が戻る δ を割線法で探していた。**空洞があって RFSW OFF のとき**、K64 の tcod は z の戻りを δ でなく tturne の update (RF の位相の基準 trf0) に吸収し、δ は 0 のまま。軌道のずれ (四極の DX は経路長を変える) がなければ z の戻りが 0 で δ = 0 のままなので、これまで見えなかった。**空洞のないリング** (回帰 bench_level3c・3d・26b: 強さ・ロールの誤差つき、保存した gs の値と 1e-13 で一致) では、gs は従来の δ による閉じ方に合う (今の gs で 3d を再確認: EMITX 1.398193521e-7、SIGE 9.759769841e-4、U0 678485.9277717 が k64 と一致) |
| 直し | 「空洞があって RFSW OFF」のときだけ z を δ で閉じない (`not (has_cav and not rfsw)`)。空洞のないリングと RFSW ON の探索は従来どおり。両モード (pysad 自身の不具合: K64 と違っていた)。**最初の直し (空洞があれば RFSW ON のときだけ閉じる) は、空洞のないリングの回帰 3 本 (3c・3d・26b) を壊したので改めた** |
| 根拠 (gs) | FODO: εx 2.70000645e-7 (gs 2.70000647e-7)、RADCOD つき 2.68403700e-7 (一致)、チューン 9 桁。大きなリング A (QF に DX = 1e-4、QD に DY = −8e-5): εx 1.5803933e-10・εy 3.1222687e-11 (gs と同じ)。PF リング: 3.62598585e-8・6.7367214e-12 (一致)。**SP8 (大きなリング B、8 GeV、D2 の例の QF5 に DX = 5e-5): gs 7.1529048e-9、k64 7.1529045e-9** (以前 6.8685e-9)。v2 は 7.1627e-9 (+0.14%。A1・A2・C の矩形の BEND の修正) |
| 根拠 (v2、独立な値) | NORFSW (δ = 0 の軌道。RF の周波数と無関係) の Jx = 1 − I4/I2 の、QF の DX = 1e-4 による変化: Emittance[] の ΔJx = −0.0054、平面リングの放射積分 (`_c3_indep.planar_quad_offset_jx`: 3x3 行列) −0.0052 (4%)。RFSW ON では経路長の変化で δ が動くので比べられない (C の試験は経路長が変わらない組) |
| 使える範囲 | 空洞のあるリングの RFSW OFF・RADCOD OFF/ON の Emittance[] (FODO・大きなリング A・PF・SP8 (大きなリング B))。空洞のないリングは従来どおり (回帰 3c・3d・26b で gs と一致)。RFSW ON で空洞のないリングの z の閉じ方は gs と比べていない。**未確認**: v2 の Emittance[] と RADINT の差 (矩形磁石で 1.4e-3、D2 が指摘した件) には触れていない |

## 3. 項目 2: 空間電荷 (WSPAC) の Emittance[] のチューンシフト

| 欄 | 内容 |
|---|---|
| 状況 | 0.2 GeV の 12 セル FODO (PBUNCH = 1e8、MINCOUP = 0.5): gs のシフト (−4.355e-3、−8.848e-3) は pysad (−4.176e-3、−8.594e-3) より 4.3%・3.0% 強い。PBUNCH・σz・εy・RF 電圧・エネルギー (0.1〜1 GeV)・BEND の分割に依らず比が一定 |
| 調べたこと | gs の `tintrb.f` の複製に、`twspace` の各呼び出しで (l, al, σx, σy, σzz, kx, ky, ⟨xz⟩, xx1, yy1) と一周の行列の跡、`gammab(l)`・`gammab(l+1)`・`p1` を出力させ、**最後の反復の 109 個のレンズ全部**を pysad (同じ量を出力) と比べた: 一周のうち 108 個は 1e-3 以内 (残りの差は反復の履歴) で一致し、**最後の 1 個 (l = nlat、ドリフトの最後の 1/4 = 0.25 m) だけ 8.0 倍** (gs kx = 5.778e-4、pysad 7.22e-5)。そのときの出力: `GAMB 75  391.3902  0.000000  195.6951  391.3902` (l、gammab(l)、**gammab(l+1) = 0**、p1、p0) |
| K64 の誤り | `tturne.f` の最後の `tintrb(…, iend+1)` (一周の最後のスライス、l = nlat) の `twspace` が `p1 = (1+δ)·(gammab(l) + gammab(l+1))/2` で、**gammab(nlat+1) は配列の外で 0**。p1 が半分になり、レンズ ∝ 1/(p1² h1) が 8 倍になる。最後のスライスはリングの先頭の前の位置に 1 個あるので、一周のチューンシフトに効く (3〜4%)。energy・電流・光学に依らない一定の割合なのはこのため (最後の要素の β と長さで決まる) |
| K64 の動作 | gs の値。k64 は `_twspace(closing=True)` で p1 を半分にして**再現** (チューンシフトが gs と全桁一致) |
| v2 の動作 | p1 = (1+δ) p0 (最後のスライスも物理的な値) |
| 根拠 | 一次のラスレットの積分 Δν = −(1/4π) ∮ β_u 2 N r0 /(β²γ³ √(2π) σz σ_u (σx+σy)) ds (ガウス分布、σx² = β εx + η² σδ²): 12 セル FODO (30 分割) で v2 は 0.15% (B1 の試験)。大きなリング A (要素ごとに入口・中点・出口の Simpson、`row_frac`): **v2/ラスレット = 0.998 (x)・1.004 (y)、k64 = 1.0085・1.0109**。PF リング: v2 1.003・1.008、k64 1.067・1.032。ASTRID2: v2 0.994・1.018、k64 1.076・1.050。**v2 (物理的な p1) がラスレットと合う** |
| k64 が gs と一致する確認 | FODO (0.2 GeV): 全桁。ASTRID2 (PBUNCH = 1e9): チューン −0.1799774472・0.1401141997 が全桁、εx 1.26803995e-8。PF リング (1e10): −0.4024743242・0.2789345926、3.45421908e-8。大きなリング A (1e10): −0.4225238563・−0.401707958、1.5497111e-10 |
| 試験 | `test_wspac_last_slice_uses_gammab_zero_in_k64_and_physical_p1_in_v2`、`test_wspac_emittance_matches_gs_k64_other_rings[3]`、`test_wspac_v2_tune_shift_large_ring_a_matches_laslett_integral_v2`、B1 の `test_wspac_emittance_matches_gs_k64` (gs と全桁一致に更新) |
| 使える範囲 | WSPAC の Emittance[]・EMIT・CALC6D・RFSW ON の CALC (B1 の範囲)。大きなリング A・PF・ASTRID2 で k64 = gs、A・PF・ASTRID2 で v2 とラスレット 0.2〜1.8% (Simpson の積分誤差が主) |
| 備考 | gs の CALC・CALC6D (WSPAC ON) が NaN になるのは別の K64 の誤りで、そのまま (B1) |

## 4. 項目 3: SynchroBetaEmittance (大きなリング A: 6.18e-9 対 1.55e-10)

| 量 (大きなリング A、3 GeV、603 m、νs = 0.004075、RF 1.4 MV) | 値 |
|---|---|
| Emittance[] の εx | 1.5497e-10 (σδ = 6.9469e-4、U0 = 296074.88 eV) |
| **独立な 6x6 の平衡** (`_c3_indep.coupled_equilibrium`: 自前の 6x6 一周写像 + 偏向磁石の刻みごとの放射の減衰 + 量子励起の拡散 + Lyapunov 方程式。RF の位置の分散・シンクロ・ベータ結合を含む) | **εx = 1.5173e-10** (σδ = 6.9458e-4、U0 = 296074.84 eV)。Emittance[] とは 2.1% (Emittance[] 側の差は A の範囲) |
| 6 次元追跡 (C の記録: RAD + FLUC、60 粒子、4000 周 = 減衰時間の 0.3) | 1.6e-10 のまま (SB の平衡なら 2.7e-9 まで上がるはず) |
| SynchroBetaEmittance (AzimuthalModes = 3 / 5 / 8) | 6.179e-9 (**3〜8 で変わらない**) |
| 同 12 / 20 / 30 / 32 以上 | 3.19e-9 / 1.46e-9 / 8.50e-10 / 7.87e-10 (**モード数で動き続ける**。pysad は 32 で打ち切り) |
| SB の νs の入力 (0.002〜0.010) | 6.1792e-9〜6.1806e-9 (**νs に依らない**: シンクロ・ベータの共鳴の効果なら強く依存するはず) |
| SB のエネルギー依存 (E/2、E、2E、4E) | 1.20e-8、6.18e-9、4.04e-9、5.72e-9 (**εx ∝ γ² に比例せず一定**: 物理的な励起ではない) |
| SB の6 極の強さ依存 (×0.5〜1.5) | 6.17e-9〜6.37e-9 (ほぼ一定)。×0.25: 6.2e-8、×0: 1.35e-7 (自然クロマティシティ) |
| ほかのリング | SP8 (大きなリング B、8 GeV): SB 7.00e-9 対 Emittance[] 6.965e-9 (+0.5%)。ASTRID2・PF: 1% (C)。24 セル FODO と一致 |

- 結論: **K64 の SB は大きなリング A で独立な値 (6x6 の平衡と追跡) と 40 倍違い、モード数に収束しない** ので使えない。pysad と gs の一致 (6.18e-9) は両方が同じ K64 の式を計算している証拠にしかならない (K64 の誤りか、適用範囲の限界かは未特定)。
- 確かめたが原因ではないもの: 閉軌道の探索の失敗 (δ の 11 点で全部収束)、二次の分散 η2 = 0.05 m・クロマティシティ ξ = −0.062 (ξσδ/νs = 0.01)、RF の位置の分散 (η = 0.0445 m、H σδ² = 2.1e-10)。
- v2 の動作: **値は変えず、収束していないときに警告する**。AzimuthalModes m の SB と 1.5 倍の m の SB の emitx が 20% 以上違えば `(pysad) SynchroBetaEmittance: AzimuthalModes->9 の emitx … は AzimuthalModes->14 の値 … と 51% 違う (モード数について収束していない。Emittance[] と追跡で確かめること)` を出す (m = 9 の大きなリング A: 5.1e-9 対 2.5e-9)。SP8 (大きなリング B)・ASTRID2 は警告しない。時間は v2 のとき 2 倍 (A は 19 秒)。
- 使える範囲: 警告が出ないリング。使えない: 警告が出るリング (大きなリング A)。

## 5. 項目 8: SP8 (大きなリング B、8 GeV) の POL の Emittance[]

| 欄 | 内容 |
|---|---|
| 再現 | SP8 (大きなリング B、1998 要素) の四極 60 台に DY (または DX) の rms 0.1 mm を付けた POL の Emittance[]: 20 分以上。誤差なし・DX = 5e-5 の QF 1 台・DY の少数台は 12〜13 秒で終わる (P = 0.9238・0.9177) |
| 原因 | **光学が不安定** (四極の誤差と強い6 極のフィードダウン。CALC も `Horizontal/Vertical Unstable.` と NaN) で、閉軌道が発散して δ = −0.999999・座標 NaN。`closed_orbit6d` の RFSW + RADCOD の分岐が (a) NaN の 1 周 (POL では v2 が四極・MULT を細かく分割するので 21〜65 秒) を、NaN でも 30 回の Newton 法で繰り返し (`np.max(abs(res)) < 1e-15` が NaN で偽のまま)、(b) 最後に 2000 要素ごとの tturne (スピンの線形の脱偏極率) を NaN で回していた。**性能の問題と不安定なスピン行列の問題が重なったもの**で、安定なリングの POL (12 秒) は遅くない |
| 直し (両モード) | (1) 閉軌道の反復で NaN になったら打ち切る。(2) POL ON・RFSW・RADCOD の `emittance6d` は先にスピンなしで軌道を確かめ、NaN なら POL の計算をせず NaN の結果を返す (`depol_error` つき)。(3) v2 の `spin_results` は軌道が NaN なら要素ごとの計算をしない |
| 結果 | 同じ 60 台 DY の例 (B3 の条件): **20 分以上 → 2 秒** (NaN と `depol_error = closed orbit not found (unstable optics)`)。DX 0.1 mm × 60 台 (CALC も NaN): 同様。安定なリングの結果は変わらない (SP8 (大きなリング B) の P = 0.9238・0.9177、12〜13 秒) |
| 試験 | `test_pol_emittance_of_unstable_ring_returns_nan_without_exception_v2` (小さな不安定なリング。SP8 (大きなリング B) の時間は試験にしていない) |
| 使えない | 不安定なリングで POL の値は出ない (NaN)。誤差の大きさ (SP8 (大きなリング B) で四極 60 台に 0.1 mm) が不安定の原因なので、補正してから使う |

## 6. 項目 10: K64 の `spdepol` (P = 1.27e-4 対 v2 の 0.119) の誤りの位置

方法: v2 の `spin.spdepol` (K64 の式のまま。`cerfc` の桁あふれだけ直したもの = `equpol_k64`) の引数を変えて、率 r = sdamp (P_ST/P − 1) の振る舞いを見た (4.2 GeV・12 MV の 12 セル FODO、`example/v2_next/B3_pol_sideband.sad`)。

| 変えたもの | K64 の式の率 [/周] | 注 |
|---|---|---|
| そのまま | 4.30e-5 | 線形の脱偏極率 (独立な導出) は 2.318e-8、側帯を入れて 4.5e-8、追跡 (10.4 ± 1.6)e-8 |
| 光子の跳びの大きさ e1・e2 ×0.1 / ×10 | 4.30e-6 / 4.30e-4 | **emittance に厳密に比例** (e は平衡のエミッタンス) |
| 減衰率 ×0.1 / ×0.5 / ×2 / ×10 | 4.09e-4 / 8.4e-5 / 2.27e-5 / 6.4e-6 | **減衰率にほぼ反比例** (率 × 倍率 = 4.1e-5・4.2e-5・4.5e-5・6.4e-5) |
| `ssprd` ×1e-3 〜 ×1e3 (スピンの位相の広がり) | 4.294e-5〜4.300e-5 (×1e-3 は NaN) | 広がりに依らない |
| `M_2PI` を 1 にする (`cm` の 2π を外す) | P = 1.265e-4 → 1.359e-4 | 影響は小さい (率は 2π ほどは動かない) |
| エネルギー 1・2・3 GeV (桁あふれを直したあとの K64 の式) | 1.8e-7・2.4e-5・2.6e-5 (線形は 1.3e-14・3.3e-10・5.5e-11) | 比は 1e7・7e4・5e5 と共鳴の分母で大きく変わる |

- **誤りの位置**: K64 の率が (emittance)/(減衰率) に比例する。独立な導出 (B3 §2.1 と DK の式): 光子の跳び (共分散 B_i、減衰に依らない) ごとの n の変化から、率 = ½ Σ tr(N B Nᵀ) は **拡散 D だけで決まり、減衰率に依らない**。平衡のエミッタンス e = D/(2λ) を使うなら率 = λ·|N|² e (減衰率に比例)。K64 は e (すでに 1/(2λ) を含む) に、2 次の摂動の**零位相の項** (`spintrm` の `cm0 = −mat/(e^{cphi} − 1)`、位相が打ち消し合う 2 次の積 `indn(i,2,0,0,0,0)`・`indn(i,1,1,0,0,0)` で cphi = −2λ のとき分母 ≈ 2λ) をもう一度掛けるので、**減衰時間の逆数の分 (1/(2λ) ≈ 1000〜1700) だけ大きい**。数値の効果は光学で決まる共鳴の分母 (1/(1−e^{iθ})) に依って 1e3〜1e7 倍 (上の表)。
- 2π について: `spintrm` の Gaussian 平均の項 `cm = −rm·2π√π e^{cb²−cφ} cerfc(cb)/a` は、`a → 0` (広がり 0) の極限で `2π·cm0` に行く (`e^{cb²}cerfc(cb) → 1/(cb√π)` と `cb = −(e^{−cφ} − 1)/a` から)。正しい規格化 (幅 0 で `cm0` に一致) は `√π/a` で、2π が余分。ただし 2π を外しても P は 1.265e-4 → 1.359e-4 にしか変わらず、1000 倍の差の原因は上の二重数えのほう。
- **docs/21 に載せる形 (番号は親が付ける)**: 「`photons.f` の `spdepol`: 2 次の零位相の項が、平衡エミッタンス (すでに 1/(2λ) を含む) に更に 1/(1 − e^{−2λ}) の分母を掛けるため、脱偏極率が emittance/減衰率に比例する (正しくは拡散 D だけで決まる)。例 (4.2 GeV の 12 セル FODO、RF 12 MV): K64 の P = 1.27e-4 (率 4.3e-5 /周)、独立な式 P = 0.119〜0.205 (率 2.3e-8〜4.5e-8)、追跡 (10.4 ± 1.6)e-8。`spintrm` の Gaussian の項の 2π も余分 (幅 0 の極限で 2π 倍)。`cerfc` の桁あふれ (B55) も合わせ、K64 の Polarization は使えない」。確度: 位置の特定は数値の依存性 (e に比例、λ に反比例、ssprd に依らない) による。K64 の式の導出そのもの (Mathematica 由来の 6 次までの展開) には戻っていない。
- 試験 `test_k64_spdepol_rate_scales_as_emittance_over_damping_v2` (率 4.0〜4.6e-5、減衰率 ×0.5・×2 で率 ×1.95・×0.53、線形の 1000 倍以上)。

## 7. 未着手の項目 (理由と、調べたこと)

| 項目 | 調べたこと・計画 | 理由 |
|---|---|---|
| (4) 傾いてずれたソレノイドの POL の量 | `spin.tradke_spin` の `xpzb = xpz + (bsir0 + bzh 2 al) pr²` と `dxpzb` (6 次元の微分) を、傾いたソレノイドの端 (trots) の前後でスピンの枠の回転を含めて直す必要がある。track6 の追跡 (スピンの枠の回転は幾何の回転そのもの: C2 の試験) が基準。小さなリングで SpinTune の差 (2.3e-4) を、傾き・ずれ・bz ごとに分けて測り、ずれたゲージの正準運動量 (px + (bz/2)dy) の補正が `bsi` に入るかを調べる | 残りの時間・文脈の制約 |
| (6) ANGLE=0 の BEND・MULT の F1 | `tsteee.f` を読んだ: `phib = −K0` (pysad の `rhob`)、入口の `tblfre(dxfr1 = −fb1²/rhob/24、dyfr1 = fb1/rhob²/6、dyfra1)` (mfring > 0 か −1 のとき、f1r = fb1/2、rbc −= f1r/al、n1 = −1)、出口は mfring > 0 か −2 (f2r = fb2/2、n2 = 2、dxfr2 = +fb2²/rhob/24)、本体の長さ `alc = al·rbc`、スライス alx・alr は `tbendal(n, ndiv, f1r, f2r, aln)` (B1 の `_tbendal` と同じ)、`tbedge(… −phib, apsi1)` は mfring ≠ −2、出口は ≠ −1。B1 の `_tblfre`・`_tbendal` が使える。**CALC (`lattice._tsteee_cod_matrix`・`_tsteee_body_mat`) も NotImplementedError** なので、Emittance[] だけでなく CALC・追跡の `tsteee` 相当も要る (MULT の `tblfri`・`tmulte` も同様)。B108 (AE1・AE2 つきの極の独立な値) は、独立な減衰分配の自前計算 (`coupled_equilibrium` に極と端を足す) が要る | 同上 |
| (7) 側帯の項と追跡の 2.3 倍 | 計画: 4.2 GeV・12 MV の FODO で、(a) ベータトロン振幅だけ (δ = 0、z = 0 の初期値)、(b) シンクロトロン振幅だけ、(c) 両方 を FLUC つきの追跡 (粒子数 480 以上、1 粒子 1 周 4.5 ms) で分けて、(a) 線形の軌道 + 厳密なスピンの MC (B3 の `B3_spin_mc.py`) との差、(b) 振幅依存のスピンチューン (側帯の式は線形の軌道) の寄与を切り出す。表 (リング・エネルギー・RF 電圧に対する追跡との比) は未取得 | CPU を長く使う (1 条件 数時間) |
| SOL 領域内の光子・RF つき MULT・SOL 領域内の加速 (POL) | `track6` の `_tsol_region*` の中の要素に `_ph_set` を足せば光子は出る (幾何は領域の座標系: `geo_frac` の使い方は未確認)。加速の SOL 領域は `_tcav_acc_ens` の領域版が要る | 同上 |

## 8. 項目 5 (POL + TRPT + 加速) と項目 9 (光子) の詳細

### 8.1 加速 (`track6._tcav_acc_ens`)
- TRPT のとき VOLT≠0 の CAVI (SOL 領域の外) を、高速追跡と同じ `track6n.cav_full` (tcav.f の全体: ndiv 分割と空洞の後の p0・h0 の更新) で通し、全粒子の `ph` を空洞の後の基準運動量にする (`ph.set_pgev`)。以降の放射・スピンの歳差 (aγ) はこの新しい p0・h0 を使う。DX・DY・ROTATE は `tent_all`・`texit_all`。k64 モードは従来どおり明示のエラー (SOL 領域と RF つき MULT・TCAVI も v2 で従来のエラー)。
- 確かめ: 1 GeV、扇形の偏向磁石 (θ1 = 0.3、θ2 = 0.2) の間に空洞 (VOLT = 0・0.2・0.5 GV、PHI = −π/2) を置いた輸送路で、縦向きのスピンの (sx, sz) が断熱加速の解析式 Φ = aγ1θ1 + aγ2θ2 (γ2 = γ1 + VOLT/mc²) の (sin Φ, cos Φ) と 5e-7 以内 (VOLT = 0.5 GV: Φ = 1.36163 rad、加速なしなら 1.13469 rad)。VOLT = 0 は加速なしの値 0.906404 と一致。
- 使えない: SOL 領域内の加速 (従来どおり SADError)、RF つき MULT・TCAVI の加速、リングの何周もの TRPT (K64 の TRPT は輸送路用: 周の終わりで基準運動量を元に戻す)、スピン追跡の損失の処理。

### 8.2 光子
- `track6` で記録するようにした: `_tthin_sext` (SEXT・OCT・…: `tsetpcvt(dx, dy, theta, 0, 0, 0, L)`)、`_tmulti_ens` (MULT: `tsetpcvt(dx, dy, theta2, …, al)` と各スライスの `fr0 += al1/al`)、磁場つきドリフト (`tdrift` の `tsetpcvt` 0)。`ffscmd` は BEND・QUAD・MULT・SEXT… (L≠0) のいずれかがあれば Python の経路 (`photons_py`)。COORD・TCAVI・MAP があるラティスは従来の高速追跡 (BEND がなければ光子を記録しない)。
- **BEAMBEAM + BEND**: `track6._step_beambeam` が `beambeam.apply` を 1 粒子ずつ呼ぶ (乱数は SadRandom の種を状態配列に出し入れ、PhotonList の粒子番号を p.idx + 1 に直す)。SADError を削除。
- **TRPT の放射の分割**: track6 の `_ndivrad` が TRPT の 10 倍を見ていなかった (高速追跡は見ていた) ので、`track_particles(trpt=…)` を足した。BEAMBEAM + BEND + TRPT の gs との乱数の状態の一致に必要。**これは POL・RADPOL・光子の追跡 (track6 の経路) の TRPT の結果を変える** (gs に近づく)。
- 照合 (k64): 3 GeV 12 セルのリング (BEND、QUAD、K0・K1・K2・DX・DY・ROTATE つきの MULT、DX・ROTATE つきの SEXT・OCT、`D1_photons_mult.sad`): 光子 325 個が gs と全列一致 (E の最大 2.7e-14)。**E の列は K64 が γ でなく gammab = p/(mc) = βγ を掛ける** (相対 1/(2γ²): 3 GeV で 1.45e-8。B3 が「amass の定数」と推定した差の正体) ので、k64 は p0、v2 は γ (h0) を掛ける。`D1_photons_bb_bend.sad` (1.7 GeV、TRPT、BEAMBEAM + BEND + QUAD): 光子 71 個 (BEAMBEAM 7、BEND 64)、追跡のあとの乱数の状態と座標の和まで gs と一致 (BEAMBEAM の行の E は 2e-9)。
- 解析式 (v2): 3 GeV、偏向角 θ の MULT (K0 = 0.02)・SEXT (K2 = 40、DX = 1 cm: θ = ½K2 DX² = 2e-3、SAD の K2 は積分値)・OCT (K3 = 1.2e4: θ = K3 DX³/6): 光子数 (5/(2√3)) α γ θ per 粒子が 4σ 以内、平均 ⟨u⟩ = 8/(15√3) u_c (u_c = (3/2) ħcγ³/ρ、ρ = L/θ) が統計の範囲 (6000 粒子)。
- 使えない: SOL 領域内の光子 (記録しない: 乱数だけ引く)、強い6 極・八極の大きな軌道のずれでリングが不安定になる条件 (pysad の CALC が NaN で gs と違う動きをした例があった)。

## 9. 項目 11: `test_b107_bend_k1_with_vertical_dispersion_v2` の失敗

現在のツリー (B3・D2 の変更後) で、単独 (2 件とも通る)・`test_v2_C3_radint.py` 全体 (21 件)・B3 の試験 5 ファイルと C3 の 2 ファイルを続けて (69 件) 通る。再現しなかった。numba のカーネルの古いキャッシュ (C の B109 の `elkern` の編集とキャッシュの削除、B3 の `track6n_ext` の編集が同時に進んでいた) か、他の担当が編集中の共有ファイルを読んだことによる一時的な失敗と考えられるが、**失敗した実行の記録がないので確かめていない**。`py -3.12 -m pysad.nbcache` は変更なし (0 ファイル)。最後の全試験でも通る (§10)。

## 10. 全試験・回帰

2026-11-18 に実行 (`PYSAD_SAD_PACKAGES` つき、k64 が既定、Python 3.12。全体を 1 本で回すと他の担当の gs・ComfyUI と競合して遅いので、ファイルを 6 組に分けて並行に実行した):

- 全試験 (`--ignore` の 3 本を除く): **1303 通過・7 skip・失敗 0** (D1 の試験 22 件を含む。`PYSAD_SLOW=1` の 1 件は skip)。最初の実行で 1 件失敗し、直した: B3 の `test_photons_with_beambeam_and_bend_not_supported` は BEAMBEAM + BEND が未対応だったことを確かめる試験だったので、対応したあとの `test_photons_with_beambeam_and_bend_are_recorded` に置き換えた。
- 回帰 (`QUICK=1`): **50 ベンチ、FAIL 0** (128 秒)。最初の実行で 3 本 (3c・3d・26b) が失敗した: 項目 1 の最初の直しが、空洞のないリング (これらの回帰の参照は gs の値) の z の閉じ方まで変えていたため。「空洞があって RFSW OFF」のときだけに限り (§2)、3 本とも元に戻った。
- 既存の試験の変更 (動作が変わった・脆かったため、3 件): `test_v2_B1_optics.py::test_wspac_emittance_matches_gs_*` (gs と全桁一致に更新)、`test_v2_B3_photons.py` (上)、`test_geo_fit_sig.py::test_beam_size_requires_emit_and_uses_twiss` (σx² = 6780 の周期解 (νx = 3.75、marginal な光学) の行列の成分との比較が絶対許容 1e-25 で、1 ulp の違いで落ちた。相対 1e-14 に)。

## 11. 直したファイル

`pysad/emit6d.py` (`closed_orbit6d` の z の閉じ方、NaN の打ち切り、POL の事前確認)、`pysad/touschek.py` (`_twspace(closing)`)、`pysad/track6.py` (光子: `_ph_adv`・`_tthin_sext`・`_tthin_gen`・`_tdrift_kick`・`_tdrift_solenoid`・`_tmulti_ens`、E の列、`_step_beambeam`、`_TRPT`、`_tcav_acc_ens`、`_cav_vnom`)、`pysad/ffscmd.py` (光子の経路、TRPT を track6 へ、POL + TRPT の誤りの条件、SB の警告)、`pysad/spin.py` (不安定なリングのガード)。
試験: `pysad/tests/test_v2_D1_emit.py` (新規)、`test_v2_B1_optics.py`・`test_v2_B3_photons.py`・`test_geo_fit_sig.py` (各 1 件、§10)。例: 冒頭の 4 組。

整理: 作業用の一時ファイルは `local/scratch_D1/` (git の追跡外。削除してよい)。WSL の `~/sadd1` (gs のソースの複製と `tintrb.f` に出力を足した再ビルド。約 380 MB) も作業用で削除してよい (元の SAD のフォルダは書き換えていない)。
