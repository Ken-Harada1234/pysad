# v2 の続き (担当 C): gs と比べられない・根拠が弱いものを、独立な物理・数値の根拠で確かめる (2026-11-18)

対象: MEASURE・TRCOD (B42)、`RadiationField` (B49)、`SynchroBetaEmittance`、DELCOR の規格化 (B99)・CAVI の符号 (B40)、A が「要判断」にした K64 の疑い (B106〜B109)、
`FITP>1` の最終解・`tffamsetup`・`CSROYMatrix`、FLUC つきの追跡の乱数列。
試験は `pysad/tests/test_v2_C3_*.py` (4 ファイルと独立な計算の 1 ファイル、計 53 件)。自前の計算は `pysad/tests/_c3_indep.py` (Liénard-Wiechert の遅延時間の求解、平面リングの放射積分 + 軌道のずれ、
結合リングの 6x6 の平衡エミッタンス (Chao の方法))。走る例は `example/v2_next/C3_*.sad` (pysad v2・k64 の出力、gs の出力つき)。
**この記録の名前の対応**: 担当名 `C` の記録 (ファイル名の `C3` は、既に `test_v2_C_*`・`test_v2_C2_*` があるための区別)。
互換モード: 修正は `compat.v2()` で分け、k64 は K64 (gs) のまま。k64 の結果を変える変更は無い (下の「k64 モードの扱い」)。

## 使ったリング (「大きなリング」)

| 名前 | 条件 |
|---|---|
| 小さなリング | 24 セルの平面 FODO (3 GeV、周長 105.6 m、RF 1.8〜4 MV)、ASTRID2 (0.58 GeV、周長 45.7 m。`example/astrid/lattice/ASTRID2_v01_01.sad`) |
| 大きなリング A | 3 GeV、周長 603 m、3563 要素、RF 1.4 MV (`pysad/tests/kek-ls_v03_12.lat`) |
| PF リング | 2.5 GeV、周長 187 m、760 要素、六極つき (`example/pf/lattice/v9_16_asym.sad`) |
| SP8 (大きなリング B、8 GeV) | 今回の項目では使っていない (SB・RADINT の確認は A・PF・ASTRID2 で足りた) |

## 1. 項目ごとの状態 (まとめ)

| 項目 | 結論 | 状態 |
|---|---|---|
| (1) MEASURE・TRCOD・ALIGN・FIXCOD | 独立な基準 (閉軌道の式、統計、固有エミッタンス、Twiss のサイズ) で確かめた。誤りは見つからず、変更なし | 確認済み (gs との比較は依然として無い) |
| (2) RadiationField | Liénard-Wiechert (遅延時間を解く) と、振幅・位相・スペクトル・全エネルギー・複数セグメントで一致。複数周の軌道で失敗する不具合を直した (重なった点を除く。両モード) | 確認済み + 直した |
| (3) SynchroBetaEmittance | 結合のあるリングで走り、gs とも (k64) 7 桁一致 (gs は落ちなかった)。Emittance[] と独立な 6x6 の平衡エミッタンスとも一致 | 確認済み (大きなリング A の値は要判断) |
| (4) DELCOR | B99 の規格化は K64 のソースから「複素数の <\|c\|²> = 1」と確定。作者の意図 (rms = delx) は推定のまま。B40 の CAVI の符号は観測不能 (追跡に効かない) | 結論 (要判断は 1 点) |
| (5) B106 | v2 は b = bzh (力学的運動量 px + bzh y)。数値の効果は 1e-13 以下 | 直した (根拠は式の一貫性) |
| (5) B107 | **誤りを 2 つ確認して v2 で直した**: (a) K1≠0 の BEND の勾配ベクトルの y 成分、(b) 軌道がずれた四極 (QUADRU) の要素の DX・DY の無視と K1<0 の符号 | 直した |
| (5) B108 | 独立な値が取れず、直していない (どう調べれば決まるかを §6 に書いた) | 要判断 |
| (5) B109 | 写像の有限差分と比べて確定: tbgrote の trans(5,6) は誤り (px² + py² のオーダー)。v2 で直した | 直した |
| (6) FITP>1・tffamsetup・CSROYMatrix | FITP の最終解は独立な最小化と一致 (pysad が最小値、gs は手前で止まる)。tffamsetup は独立な基準を作れず。CSROYMatrix は原因を特定できず | 一部確認 |
| (7) FLUC の乱数列 | **差の原因は gs の既定 NPARA=4 (子プロセスに分ける)**。NPARA=1 では pysad (k64) は gs と全桁一致、`NPARA=4;` で分け方を再現 (新機能) | 解決 |

## 2. 直した項目 (詳細)

| ID | K64 の動作 | v2 の動作 | 根拠 (式・数値) | 試験 | k64 モードの扱い |
|---|---|---|---|---|---|
| B107a (RADINT の BMAG) | `CKL1Y = −cos(θ)·DKL1`。θ=0 (回転なし) の K1≠0 の BEND でも y 成分が 0 にならず、Y 方向の分散 η_y に 2 k h η_y の項が余計に入る | `CKL1Y = +sin(θ)·DKL1` (磁場の勾配ベクトル = k·(cos θ, −sin θ)、F4 の式は CKL0Y − 2 CKL1Y/CURV)。θ=0 で y 成分は 0 | 24 セル FODO + skew 四極 + K1 = ±0.1 の BEND: I4U + I4V = 物理的な I4 (結合のある 5x5 行列、η_x h(h²+2k) の積分) は v2 で 1e-7 (K64: K1 = +0.1 で +1.6%、−0.1 で +2.7% ずれる)。DU・DV・Je は 1 次の摂動論 (6x6 の固有ベクトル) と 3e-3 以内 (K64 は 0.3 ずれる)。K1=0 または結合なしでは両モード同じ | `test_v2_C3_radint.py::test_b107_bend_k1_*` | 従来どおり (gs と全桁一致: `C3_radint_k1_bend_vdisp.gs.txt`) |
| B107b (RADINT の QUADRU) | 四極の中の軌道の位置 (dx, dy) に Twiss の閉軌道 DX・DY だけを渡す (四極自身の据え付けのずれ DX・DY を引かない)。K1<0 の四極では勾配ベクトル CKL1 の符号が逆 (DKL1 の符号つき) | (dx, dy) = 閉軌道 − 要素の DX・DY。勾配ベクトル |K1|·(dx, dy)/r (|B| = |k| r の勾配) | 平面リング (24 セル FODO、QF・QD の DX = 1e-4・2.2974e-4 (経路長が変わらない組)): 3x3 行列の放射積分 (h = k x_o、I4 の 2k²x_o η) Jx = 0.949255、RAD つき 6 次元追跡の ΔJx = +0.00781、Emittance[] +0.00786、v2 の RADINT +0.0078 (K64 +0.0028)。QF だけ・QD だけのずれでも独立な積分と 3e-4 以内。Y 方向のずれ + skew 四極: RADINT の ΔJy = +0.00042・−0.00042・+0.00121 が追跡と 1e-4 以内 (K64 は符号が逆) | `test_b107b_*` (10 件) | 従来どおり (gs と全桁一致: `C3_radint_quad_offset.gs.txt`) |
| B109 (tbgrote の trans(5,6)) | `r33(dpzf pr/pzf + dpzi/pzi) − dr33` (dr=0 でも dpzi(pr+1)/pzi = −(px²+py²) が残る)。コメントアウトされた式 `zf/pzf (1 − pr/pzf r33 dpzidp)` が写像 z += pr zf/pzf の厳密な微分 | コメントアウトされた式 (dr=0 では 0) | tchge の解析的な 6x6 行列と写像 (`_tchge_cod`) の有限差分 (3 つの条件、入口・出口): v2 は全成分 4e-12 以内、K64 は (5,6) が 1.3e-5 (= px² + py²) ずれる。大きな差ではない (2 次の量) が、写像と行列の整合は厳密であるべき | `test_b109_*` (8 件) | 従来どおり |
| B106 (tradke のビームの広がりの項) | `b = bzh·0.5`。同じ関数の 10 行上で力学的運動量を `px = cod(2) + bzh·cod(3)` と定義している | `b = bzh` (⟨(px + b y)²⟩ の b。B39 と同じ `fixb39` の v2 だけ) | **数値の効果は出ない**: 結合のある SOL つきリングの Emittance[] は 13 桁変わらない。拡散の項は ⟨u²⟩⟨P_⊥²⟩ (u = 光子のエネルギー比) で、横方向の減衰 ⟨u⟩ に対し ⟨u²⟩/⟨u⟩ ≈ 1e-6 だけ小さい。式の一貫性による変更 | `test_v2_A_spin.py::test_tradke_px_py_diffusion_uses_b_y_py_v2` (b = bzh の式に更新) | 従来どおり (b = bzh/2) |
| RadiationField の複数区間 (pysad の不具合) | 区間の始めの点と直前の要素の後の点が同じ (同じ時刻・位置) で、`x must be strictly increasing sequence` で失敗 (gs は RadiationField で SIGSEGV) | 完全に同じ点を除く。時刻が戻る・同じ時刻で位置が違う点は `SADError` (軌道の誤り) | 戻りの点数が重なりなしの点数、場は重なりなしの軌道と同じ | `test_v2_C3_radfield.py` の 3 件 | 両モードで動く (pysad 自身の不具合: gate なし) |
| `NPARA` の fork の分け方 (新機能) | gs の既定は NPARA=4: `TrackParticles` は np ≥ 13 で子プロセスに粒子を分け、子の乱数の種は親の種 + 2k、親が残りを元の種で追う | `NPARA=n;` を FFS で設定すると同じ分け方で追う (有効な並列数 = min(NPARA, np/9 + 1, ne·np/3000 + 1)、RAD + FLUC の追跡だけ。POL・RADLIGHT・PHOTONS・WAKE は分けない) | gs (4 並列の既定と NPARA=1) と全桁一致 (§3.7) | `test_v2_C3_misc.py::test_fluc_*` (5 件) | 既定は NPARA=1 (従来どおり。gate なし: 純粋な追加) |

## 3. 項目ごとの確認の詳細

### 3.1 (1) MEASURE・TRCOD・ALIGN・FIXCOD (すべて v2。`test_v2_C3_measure.py` 8 件)

独立な基準と結果 (大きなリング): **変更なし (誤りは見つからなかった)**。

| 確かめたこと | 独立な基準 | 結果 |
|---|---|---|
| TRCOD が測った閉軌道 (薄いキック θ = 1e-5 が 1 か所) | (a) CALC の閉軌道、(b) 解析式 x(s) = θ√(β_s β_k)/(2 sin πν) cos(πν − \|φ_s − φ_k\|) (Twiss の β・位相) | TRCOD と CALC の差 / 最大: ASTRID2 1.8e-6、PF 6e-7、**大きなリング A 1.3e-5**。解析式との差: ASTRID2 4.3e-4、PF 6.4e-4、A 1.9e-3 (六極のフィードダウン: 軌道 8e-5 m の六極の中心からのずれで、線形の式の外) |
| MEASURE の初期分布 (結合のある PF の周期解、R1..R4 ≠ 0、DP=0) | 4x4 の固有エミッタンス (Σ·J の固有値) = EMITX・EMITY | N = 2e4: 3% 以内、N = 2e5: 1% 以内 (誤差 ∝ 1/√N) |
| MEASURE のビームサイズ (PF、EMITX = 1e-9・EMITY = 3e-10、DP = 1e-4、N = 20000、3 か所の四極の入口) | σ = √(ε β + (η δ)²) (Twiss) | 4% 以内 (統計誤差 0.5%) |
| FIXCOD・補正後の残差 (PF) | 基準 (キック 1e-5) とキック 3e-5 の CALC の軌道の差 = `Twiss["DDX"]`。キックを戻したあとの DDX | DDX = CALC の差と 1e-12。戻した後の残差 5e-8 (CALC の閉軌道の収束の精度、軌道の最大値 8e-5 の 2e-4) |
| ALIGN (PF、QUAD・SEXT・BEND の区間) | 要素1 の DX・DY の値 | 区間の磁石の DX・DY が 0 差で一致 (注: 同じ名前の要素は `NAME.n` で指定。名前だけだと最初の 1 個) |

一度に確かめた周回の整合 (matched distribution を 1 周追った分散行列が元と一致するか) は、色収差が効くので小さな δ (DP = 1e-5) でも 5% 程度ずれた。
(クロマティシティ ξ δ の位相のずれと RF の項。DP = 1e-3 では 40% ずれる: 4D の周期解の分布は δ のあるビームの不変分布ではない。MEASURE の問題ではなく、「周期解 = 不変な分布」の前提が δ で壊れる。試験には入れていない。)

### 3.2 (2) RadiationField (`test_v2_C3_radfield.py` 8 件、例 `C3_radfield_bend.sad`)

独立な計算 (`_c3_indep.lw_field`): 円弧・直線の軌道 (解析的な位置・速度・加速度) で、観測時刻 τ ごとに遅延時間 t' を brentq で解き、
E = q/(4πε0)[(n−β)/(γ²(1−n·β)³R²) + n×((n−β)×β̇)/(c(1−n·β)³R)] (Liénard-Wiechert) を評価する。pysad の RadiationField は K64 の Feynman-Heaviside の式 (軌道の点列と自然スプラインの微分)。

| 確かめたこと | 条件 | 結果 |
|---|---|---|
| 振幅・位相 (τ ごとの場の一致) | 電子 γ = 20、ρ = 1 m、弧 ±0.5 rad の接線 50 m。点 2001・4001・8001 | 主成分 (水平) の最大の差 / ピーク: 6.7e-5 → 1.7e-5 → **4.2e-6** (2 次収束)。全成分の差を全体の最大値で規格化して 2e-5 以内 (8001 点) |
| 振幅 (v2 の RADLIGHT の軌道) | 10.2 MeV、ρ = 1 m、BEND 0.4 rad、RADLIGHTDIV = 16・64、観測点 50 m | \|E\| の最大: 16 点で 2.6e-3、64 点で 1.6e-4 (docs/45 の 0.25%・5e-5 と整合) |
| スペクトル | 同、観測点 500 m。RadiationSpectrum の \|c + is\|² → d²I/dωdΩ | Jackson の式 (e²/4πε0)/(3π²c)(ωρ/c)²γ⁻⁴K²₂/₃(ωρ/3cγ³) と ω = 0.1〜2.6 ωc で **1% 以内** (実測 0.3%)。ピークの位置 (臨界周波数 ωc = 3γ³c/2ρ に対する) も一致 |
| 全エネルギー | 観測点を水平 ±0.5 rad の弧の接線方向に置き、垂直 ±0.3 rad (31 点) で ε0 c R² ∫\|E\|² dτ cosθ dθ | 偏向 1 rad あたり e²γ⁴β³/(6πε0ρ) (Larmor の式 e²cγ⁴β⁴/6πε0ρ² × 弧の時間) と **1e-3 以内** (4e-4〜1e-3) |
| 複数セグメント | 直線・円弧 (+)・直線・円弧 (−) のシケイン、区間の境目の点を重ねた 6001 点 | 厳密な LW (2 つのパルス) と最大値の 1.5e-4。重なった点が除かれる |
| 1 GeV (v2 の RADLIGHT、BEND L = 1.2 m・ANGLE = 0.2618、8200 点、RADLIGHTDIV = 16) | 軸上のスペクトル、観測点 **5 km** (ω = 0.2〜5 ωc、λc = 2.56 nm) | Jackson の式と 4e-4 (ω によらず一定)。観測点が 50 m だと光源の大きさ ρ/γ に対し Fresnel 領域で 4.6% 小さい (物理) |
| 複数周の RADLIGHT の軌道 | k64・v2。`TrackParticles[d, 9, 1, 2]` の軌道 (22 点 / 32814 点、重なった点が 1) | RadiationField が場を返す (以前は strictly increasing のエラー) |

使える範囲: 電子の円運動・直線の軌道 (RADLIGHT の v2 の軌道)、点の間隔が 1/γ の角度あたり 16 点以上 (振幅 3e-3)。64 点で 2e-4。観測点は光源の大きさ ρ/γ に対しスペクトルの λ で Fraunhofer の距離 (≫ (ρ/γ)²/λ)。
使えない・注意: k64 モードの軌道は BEND の中を記録しないので放射場が意味を持たない (docs/45)。軌道の点に位置・時刻のずれがあると、パルスの時間幅 ρ/(2γ³c) (1 GeV で 1e-18 s、位置 3e-10 m) の精度が要る: v2 の RADLIGHT の位置・時刻は 1e-13 m だが、RAD が ON のため軌道が放射で変わる
(1 GeV の 1 台の偏向磁石で横方向に 8e-8 m): 同じ軌道を別の計算に渡すときは放射の効果が同じかを確認する。gs は RadiationField で SIGSEGV (B49。`C3_radfield_bend.gs.txt`)。

### 3.3 (3) SynchroBetaEmittance (`test_v2_C3_misc.py` 5 件、例 `C3_synchrobeta_coupled.sad`)

- **結合のあるリングで動く**。24 セル FODO + skew 四極 1 台 (3 GeV、RF 4 MV、六極なし): **gs も落ちず、k64 の pysad は gs と 7 桁以上一致** (εx 1.3603366e-6 対 1.3603365e-6、εy 4.7516654e-8 対 4.7516654e-8、`C3_synchrobeta_coupled.gs.txt`)。docs/46 の「結合のあるラティスで gs が落ちる」は再現しなかった。
- 独立な計算 (`_c3_indep.coupled_equilibrium`: 自前の 6x6 の 1 周写像 + 偏向磁石の刻みごとの放射の減衰 (摂動) + 量子励起の拡散 D = Σ⟨u²⟩ S eδ eδᵀ Sᵀ + Lyapunov 方程式、固有エミッタンス): ASTRID2 + skew 四極 (ks = 0.05、結合が強い) で
  ε₁ 1.21049e-8 (Emittance[] 1.21042e-8: 6e-5)、ε₂ 1.14506e-9 (1.14491e-9: 1.3e-4)、σδ 4e-5、U0 2e-6。**Emittance[] の結合のあるリングの値は独立に確かめられた**。
- SB は Emittance[] と、ASTRID2 で +0.8% (x)・+0.8% (y)、PF で +1.1% (x)・+2.9% (y、εy が 1e-10 と小さい)。この差は結合の有無で変わらない (比が 0.5% 以内): SB の差はシンクロ・ベータ結合と色収差の分 (結合の効果ではない)。
- 和則: 弱い結合 (ASTRID2、ks = 0.005) では ε₁ + ε₂ が結合なしの値 (EMIT・SB とも 2e-3 以内。+0.05%)。強い結合 (ks = 0.05) では +4.2% (ν_x − ν_y ≈ 0.04 の近く、Y 方向の分散ができる) で和は保存しない。
- **要判断 (大きなリング A)**: SB は 3 GeV、周長 603 m、νs = 0.0041 のリングで、gs も pysad k64 も **6.18e-9 (AzimuthalModes ≥ 2 で収束)** を返すが、Emittance[] の 1.55e-10 の 40 倍で、RAD + FLUC の 6 次元追跡 (60 粒子、Emittance[] の分布から開始) の ε は 4000 周 (減衰時間 13000 周の 0.3) で 1.6e-10 のまま
  (SB の平衡に向かうなら 4000 周で 2.7e-9 になるはず)。追跡は ring A + 0.1 m のドリフトで行った (SB はこの変更で不安定になり、モード数 1・2・3・5 で 3.6e-9・1.3e-9・6.1e-9・7.8e-10: 収束しない)。
  つまり K64 の SB は大きなリング A では追跡と合わない (K64 の式の限界か、pysad・gs が共通に持つ未特定の誤りか)。PF (2.5 GeV) と ASTRID2 では Emittance[] と 1% で合う。使うときは AzimuthalModes を変えて収束を確認し、Emittance[] と比べること。
- pysad の `AzimuthalModes->1` は直接の呼び出しで IndexError (K64 の配列の外を読む。gs は 8.6e-7 という値を返すが意味がない)。2 以上を使う。

### 3.4 (4) DELCOR (`test_v2_C3_misc.py` 3 件)

K64 の `tcorr.f`・`tvcorr.f`・SADHelp の説明 (SADHelp.html に DELCOR の記述は無い) から:
- **B99 (corlength ≠ 0 の規格化)**: スペクトルの重み w_k = exp(−(κ x0)²/4) は相関関数 ⟨δ(s)δ(s')⟩ ∝ exp(−Δs²/(2 corlength²)) を与える (既存の試験で確認済み)。K64 の規格化 s = 1/√(4Σw² + 2) は **複素数 (Re + i Im) の ⟨|c|²⟩ = 1** にするもので、実部の分散は 1/2、DX = Re·delx の rms は delx/√2 (試験 `test_delcor_tvcorr_normalisation_of_the_correlated_branch`: k64 は ⟨|c|²⟩ = 1.0、v2 は 2.0)。
  一方 corlength = 0 の枝 (ランダムウォークの橋) は実部・虚部それぞれ分散 al/n の刻み (1/√2 の補正なし。中央の分散 al/4: `test_delcor_tvcorr_zero_corlength_is_brownian_bridge_per_component`)。K64 の 2 つの枝で「実部・虚部の分散の約束」が揃っていない。
  実部と虚部を別の自由度 (DX と DY) に使う使い方では、各成分の rms が delx・dely になる約束 (v2) が自然だが、**作者の文書による裏づけは無い** (推定のまま)。K64 の値が必要なら k64 モードか delx に 1/√2 を掛ける。
- **B40 (CAVI の符号)**: `tcorr.f` の `id=20` (SOL) と `id=31` (CAVI) の枝は `dx1 = (設計の値) − delx` (add のときは −delx) で、ほかの磁石 (+delx) と符号が逆である。`id` は tfstk.f で 20 = icSOL、31 = icCAVI、41 = icMARK。
  SOL の枝は境界の SOL の組 (BOUND の入口と出口: `istep = ie − i0`) に同じ値を書き、tsol は入口で x −= dx、出口で x += dx なので、−delx は領域全体をずらす符号と読める (SOL の領域の DX の約束)。CAVI の枝はこの SOL の枝の写し (書く先の添字 +13・+14 は今の CAVI の DX・DY)。
  しかし K64 の CAVI は追跡で DX・DY を使わない (wake の `txwake` に dx = dy = 0 を渡す: docs/38)。**符号は観測できない** (`LINE["DX","RFC"]` の読みだけ)。v2 は +delx のまま (変更なし、要判断としても実害なし)。

### 3.5 (5) B106〜B109

(§2 の表のとおり。B108 は未解決で §6。) B107 の確かめ方の要点:

1. A2 の記録 (§4.3) は K1=0 の BEND だけで I4V を確かめ、K1≠0 の BEND に Y 方向の分散がある場合は「未確認」だった。この場合 (24 セル FODO + skew 四極 + K1 = ±0.1 の BEND) で、RADINT の I4U + I4V が物理的な I4 と 1.6〜2.7% 合わないことを確認 (K64)。
2. y 成分の勾配ベクトルを +sin(θ)·K1 にすると 1e-7 で合う。物理: 偏向磁石の磁場 B = (B0 + k u) v̂ + k v û (û = 偏向の向き = (cos θ, −sin θ))、|B|² の勾配 = k û 方向。
3. 四極の中の軌道: 平面リングの 3x3 行列の放射積分 (h = k x_o) と、RAD つき 6 次元追跡の 1 周写像 (RADINT・Emittance[] のどちらも使わない) の両方が、RADINT に「ずれのある四極」の項が足りないことを示す。
   **注意**: 軌道のずれが経路長を変える (δ の閉軌道がずれる、RF 周波数固定) とき、Emittance[]・追跡は δ ≠ 0 の軌道 (RF が δ を決める) の放射を含み、RADINT は δ = 0 の軌道だけ。独立な計算は経路長の変化が 0 になる QF : QD のずれの比で比べた。
   一般のずれでは RADINT と Emittance[] の J がずれの分の δ の寄与だけ違う (QF の DX = 1e-4 で ΔJx = −0.005 対 +0.007)。RADINT は「RF の周波数を軌道長に合わせて δ = 0」の約束で読む。
4. SEXT・BEND の DX・DY・軌道のずれ (K1≠0 の BEND の勾配の項、六極のフィードダウン) は v2 でも RADINT に入れていない (独立な値を作っていない)。傾いた (ROTATE ≠ 0) K1 の BEND の y 成分 (sin θ) は θ = 0 の場合の確認から物理的に導いた式で、θ ≠ 0 の数値の確認はしていない (傾いた BEND は閉軌道が大きく動いて線形の確認ができない)。

### 3.6 (6) FITP>1・tffamsetup・CSROYMatrix

- **FITP>1 の最終解**: `example/03_matching/FFS_FITPOINTS.sad` (24 セル FODO、QF QD SF SD、δ = −1%・0・+1% の 3 点で NX = 6.7・NY = 4.9)。K64 の残差の定義 (Σ w (f − goal)² scale²、w = 1 (δ=0)・(OffMomentumWeight/2/√npoints)² = (1/4)² (δ=±1%)) の独立な最小化 (CALC のチューンに scipy の最小二乗) の最小値:
  K1(QF, QD) = (0.70021085, −0.59962825)、K2(SF, SD) = (1.69527587, −2.43892813)、残差 6.33178e-6。**pysad の GO の最終解 (0.70021085、−0.59962825、1.69527582、−2.43892820) は独立な最小値と 5e-8 以内で一致**。
  gs は (0.700211, −0.599628、1.695226、−2.439058) で止まり、残差 6.3321e-6 は最小値より 5e-5 (相対) 大きい (K2 が 5e-5・1.3e-4 違う)。残差の谷が平らで、`*****qcod---> Overflow & closed orbit not found` (gs 2 回、pysad 4 回) の試行点での打ち切りで gs は最小値の手前で止まる。再起動 (docs/46 の「GO で再起動」) は今回の gs の実行では起きなかった。
  試験: `test_fitp_final_solution_is_the_independent_minimum` (k64 の GO と独立な最小化: K 値 5e-6、残差 2e-4)。
- **tffamsetup (fam 点、MatchingAmplitude・NFAMP)**: 独立な基準を作れなかった。有限振幅の周期解の条件 (振幅のある粒子のチューン) の独立な計算 (追跡の NAFF、解析の振幅依存チューンシフト) が要る。既存の `test_fam_points.py` (自己検証) のまま。
- **CSROYMatrix**: gs が wrongtype を返す原因を特定できなかった。`tfcsroymat` (csroy.f) の型検査は「6 要素のリスト (全部が実数のリスト) と 5 個の実数のリスト」で、`tfcsroysetup` の戻りは `ktadaloc(-1,6)` の 6 要素のリストに見え、ソースの読みでは通るはずだが gs が拒否する (登録名の綴り・`tflistq`/`tfreallistq` の判定のどちらか)。
  gfortran でソースを再ビルドしないと確かめられない。pysad の `CSRSetupOY` の zj・dpj の修正は蛙跳びの関係で確認済み (B59)。

### 3.7 (7) FLUC の乱数列

調べ方: 同じ台本を gs (WSL) と pysad (k64) で走らせ、粒子ごとの dp・最終の種を比べた。
1. np ≤ 12 では gs と pysad が全桁一致、np ≥ 13 で全粒子が食い違う → 粒子数で動作が変わる。
2. `tftrack.f:124-195`: `npara = min(NPARA, 256, np/9 + 1, ne·np/3000 + 1)` (npparamin = 9、npnlatmin = 3000)。npara > 1 なら `np1 = np/npara` ずつ粒子を分け、子プロセス k (k = 1, 2, …、条件 `ipn + np1 < np` かつ `npr < npara`) は **乱数の種 = 親の種 + 2k (`tfaddseed`: `SeedRandom[vn + kseed]`)** で `[ipn, ipn + np1)` を追い、親は残りを元の種で追う。親の最終の種だけが残る。
3. gs の NPARA の既定は 4 (計算機で決まる)。pysad の NPARA の既定は 1。そのため「gs 4.21e-4、pysad 4.52e-4」の差は NPARA=4 と NPARA=1 の比較だった。
4. pysad に同じ分け方を足した (`ffscmd.FFS._npara_chunks`・`_fn_track_npara`: 子に種 + 2k を `rng.restore` し、子から順に追い、最後に親を元の種で追う)。

確認 (k64、gs との一致): FODO リング (np = 200、1・10・100 周): NPARA=1 の rms dp 5.252739925347e-4・3.14056869362e-3・2.36045313644e-3 と最終の種 (4063185639・3350258407・4157604543)、NPARA=4 の 5.205643226207e-4・3.14042063085e-3・2.37021386166e-3 と種 (3890937727・2381179799・4001722967) がすべて gs と全桁一致。
np = 16 (有効 npara = 2) の粒子ごとの dp も 16 個すべて 12 桁一致。**PF リング (760 要素、六極つき)、np = 40、1・10・40 周: NPARA=1 と 4 の両方で rms dp が 13 桁、最終の種が gs と同じ**。v2 は A1・A2 の追跡の放射の式の修正 (ddpx の符号など) で rms が 6e-6 (1 周)〜1e-3 (100 周) 違う (乱数の引き方・最終の種は同じ)。
このため「統計的な一致の範囲」を調べる必要は無くなった (乱数列が一致する)。

## 4. docs/46・docs/47 に反映してほしい訂正 (親が更新)

| 文書の記述 | 今回の確認 |
|---|---|
| docs/46 §5「`SynchroBeta` の結合ラティスで gs が落ちる」、§3「結合のあるラティスは gs が落ちる」 | 24 セル FODO + skew 四極 1 台の結合のあるリングで、gs は落ちず k64 の pysad と 7 桁一致 (`C3_synchrobeta_coupled`)。ほかの結合の条件 (SOL・強い結合) で落ちるかは未確認 |
| docs/46 §3「乱数列 (FLUC つきの追跡) δ の rms が gs 4.21e-4、pysad 4.52e-4」 | 原因は gs の既定 NPARA=4。NPARA=1 で全桁一致、pysad に `NPARA=4;` で同じ分け方を足した (§3.7) |
| docs/46「`FITP>1` の最終解 … gs が `GO` で再起動する」 | 今回の gs は再起動せず完走し、最小値の手前で止まった。pysad は独立な最小化と 5e-8 で一致 |
| docs/46 §3「`RadiationField` … 速度場との比較のみ」 | Liénard-Wiechert の厳密な数値計算と振幅・位相・スペクトル・全エネルギーで一致 (§3.2) |
| docs/46 §4「複数セグメントの RadiationField の失敗」 | 直した (§2) |
| docs/46 §5「gs が落ちる: `RadiationField` の SIGSEGV」 | 再現 (`C3_radfield_bend.gs.txt`)。ほかの引数の渡し方でも落ちる |
| docs/47 (e)「DELCOR の corlength≠0 の規格化 (B99): v2 は rms = delx (推定)」 | 規格化の根拠を確定 (K64 は複素数の ⟨\|c\|²⟩ = 1、corlength = 0 の枝は成分ごと)。意図は推定のまま |
| docs/47 (e)「DELCOR の CAVI の符号 (B40)」 | 追跡で使われない値なので観測できない (§3.4) |
| docs/47 (e)「RADINT の CKL1Y (B107)」 | 誤りと確定して v2 で修正 (§2)。さらに QUADRU の 2 つの誤りを追加 |
| docs/47 (e)「tbgrote の trans(5,6) (B109)」 | 誤りと確定して v2 で修正 |
| docs/47 (e)「tradke の b = bzh·0.5 (B106)」 | v2 で b = bzh に直した (効果は 1e-13 以下) |
| docs/21 B107 | 「K1≠0 の BEND の CKL1Y (y 成分が θ=0 でも 0 にならない)」と「QUADRU が四極自身の DX・DY を引かず、K1<0 の勾配ベクトルの符号が逆」の 2 つに分けて確度「確定」に。B106・B109 も確度「確定 (B106 は式の一貫性)」 |

## 5. 試験の数と実行結果

`test_v2_C3_radint.py` 21 件、`test_v2_C3_radfield.py` 8 件、`test_v2_C3_measure.py` 8 件 (うち 1 件は `PYSAD_SLOW=1` の大きなリング A)、`test_v2_C3_misc.py` 16 件 (SB の gs 参照 1 件を含む)、計 53 件。
既存の `test_v2_A_spin.py::test_tradke_px_py_diffusion_uses_b_y_py_v2` は B106 の b を変えたので、期待値を b = bzh の式に直した (他の試験の変更は無い)。全試験・回帰の結果は §7。

## 6. 要判断・残した項目

| 項目 | 調べたこと | 足りないもの |
|---|---|---|
| B108 (ANGLE=0 の BEND・K0 つき MULT の端の角) | `tsteee.f` は端の角 apsi1・apsi2 (SAD のキーワード AE1・AE2。E1・E2 は ANGLE との比なので ANGLE=0 では 0) を `tbedge` で蹴りに入れるが、`tradke` に行路の補正 (B75 の tin・tout) を渡さない。ANGLE≠0 の BEND の B75 と同じ構造 | 独立な値: ANGLE=0・K0≠0 の極 (AE1・AE2 ≠ 0) のあるリングの減衰分配の独立な計算 (`coupled_equilibrium` に極と端を足す。軌道の曲がり・端のフォーカスを自前で入れる必要がある)。K0 の極の磁場の端の幾何も独立に決める必要がある。リポジトリの ASTRID2 の MPW の極は AE1 を使っていない (E1 = 0.5 だが ANGLE=0 なので 0) ので、効かない |
| tffamsetup (fam 点) | 有限振幅の条件の独立な基準を作れなかった | 振幅のある粒子の追跡のチューン (NAFF) か、振幅依存チューンシフトの解析式 |
| CSROYMatrix | gs の wrongtype の原因が特定できない (§3.6) | gfortran でのソースの再ビルドとデバッグ |
| SB (大きなリング A) | §3.3 | SB の式の導出 (K64 の temits の `amus`・`dispp` の使い方)、または 6 次元追跡の長時間 (4τ ≈ 50000 周、60 粒子で 30 分) の平衡との比較 |
| RADINT の軌道のずれ | 四極の DX・DY を v2 で入れた。SEXT (K2 のフィードダウン)・BEND (DX・K1 の勾配) は K64 の RADINT が軌道を無視する作りのまま | 独立な基準 (追跡の J) で ΔJ を測る (四極より桁が小さく、試験の分解能 1e-4 未満だった) |
| DELCOR の意図 (B99) | §3.4 | 作者の文書 (SAD 公式の説明) |
| 周期解を matched 分布として 1 周追う試験 | δ が大きいと色収差で崩れる (DP = 1e-3 で 40%) | 周期解の δ 依存の分布 (4D の ε でなく 6D の不変量) |

## 7. 変更したファイルと実行結果

- `pysad/radint.py`: B107a (`bmag` の CKL1Y)、B107b (`quadru` の |DKL1|、`compute` の四極の DX・DY)。いずれも `_compat.v2()`。
- `pysad/lattice.py`: B109 (`_tbgrote_mat` の (5,6)、`_tchge_jac` の simple の分岐)。`pysad/elkern.py`: B109 (`_quirk56`・`_tbgrote` に `v2` を足して `_tbrote` から渡す)。numba のキャッシュ (`elkern.*`・`turnkern.*`・`perfkern*`) を削除した。
- `pysad/perfkern.py`: B106 (`tradke_db` の b。B39 と同じ `fixb39`)。
- `pysad/radfield.py`: 重複点の除去 (`_drop_repeated_points`)。`pysad/ffscmd.py`: NPARA の fork の分け方 (`_npara_chunks`・`_fn_track_npara`、`fn_track` に 2 行)。
- 試験: `pysad/tests/test_v2_C3_{radint,radfield,measure,misc}.py`、`_c3_indep.py`。`test_v2_A_spin.py` の 1 件を更新。
- 例: `example/v2_next/C3_*.sad` (6 本) と `.gs.txt`・`.pysad.txt`・`.pysad_k64.txt`、`C3_README.md`。
- 作業用 (削除してよい): `local/scratch_C/` (git の追跡外)。
