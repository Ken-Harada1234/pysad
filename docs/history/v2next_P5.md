# v2next の記録: 担当 P5 (検証の拡充、2026-11-19)

v2.1.0 の `docs/46`・`docs/47` §(d) と `docs/history/v2next_P1〜P3.md` で「未確認」とされた 9 項目を、gs (SAD K64) との比較で確かめた。
gs は WSL の 10-19 ビルド (2026-11-19、`timeout` 120〜900 s、各回の後に `ps` で sad1.exe の残りが無いことを確認。1 回だけ `CorrectOrbitBump` を直接呼ぶ台本が 300 s で終わらず timeout で止まった)。
`docs/46`・`docs/47`・`docs/21`・README・CHANGELOG は触っていない (新しい K64 の誤りの候補は §11 に番号なしで書いた)。コミットしていない。
例: `example/v2_next/P5_*` (説明 `P5_README.md`)。試験: `pysad/tests/test_v2_P5_{measure,detr,micado_fitp_oy,specfn_cx,proton_beta,syntax,tk90}.py`、補助 `_p5_util.py`・`_p5_gscx.py`。
使ったリング: 小さな FODO (3 GeV、周長 48〜52.8 m)、大きなリング A (`kek-ls_v03_12.lat`、3 GeV、603 m)、PF (2.5 GeV、187 m)、ASTRID2 (0.58 GeV、周長 45.7 m)。

## 1. 項目の状態 (要約)

| 項目 | 確認できたこと | できなかったこと | コードの変更 | モード |
|---|---|---|---|---|
| 1. MEASURE・TRCOD | **gs は今は MEASURE を出力する** (B42 の更新): 先頭の要素 (入口の MARK) の位置だけ。k64 を gs と同じにし、4 つのリングで統計ブロックが 4 桁一致 | 先頭以外の位置は gs が追跡で全粒子を落とし (無出力、TRCOD は SIGSEGV)、v2 の追跡経由の値は gs と比べられない | `measure.py` に k64 の MEASURE・TRCOD・LOGMEASURE | k64 (K64 の動作の追加) |
| 2. DETR >= 1 の INS | 90 度回転 (COORD) の対で DETR = 1.375 の区間を作り、CELL CALC・INS・Emittance[] を 3 つのリングで gs と比較: β・DETR・Emittance[]・U0 は一致。**一周の位相の整数部が ±1 ずれる** (K64 の丸め誤差の符号) | MARK の生の DETR >= 1 を入口にする INS (gs の挙動が説明できない) | `lattice._exch_phase`・`perfkern` の kernel (v2) | v2 で修正、k64 は K64 のまま |
| 3. 複素の特殊関数 | 732 式を走査: 追加前 436 式 → **641 式が gs と 1e-9**。15 関数を追加。残り 91 式は gs が誤る・不正確 | gs の誤りの再現 (k64) は K64 の gamma.f・D.n の直訳が要るので要判断のまま | `specfn_cx2.py` (新)・`specfn.py`・`sadscript.py` | 両モード (追加) |
| 4. OYStability・OYScan | 後半 (CSROYMatrix → Eigensystem → 成長率) を gs と 1e-12 で比較 (安定・不安定の 2 通り) | 通しの `OYStability`・`OYScan` は gs が動かず (B60)、pysad でも CSR の構築を含む通しは未確認 | なし | — |
| 5. Micado の変種 | 水平・2 平面同時・`SVDmode->1`: k64 が gs と 5e-6。v2 (B129) は単調に下がる。`MakeBump` は gs と 1e-12 | `CorrectOrbit` の `Bump->` は K64 の誤りで gs も pysad も動かない。`CorrectOrbitBump` の直接呼びは gs が 300 s で終わらず | なし | — |
| 6. FITP>1 | 5 点・3 点 + β・FITP 1: 最終の残差が gs と 2e-4 以内、K1 が 6 桁 | 反復の経路 (gs は早く打ち切る) の再現 | なし | — |
| 7. Tk 9.0 | Python 3.14.9 + Tk 9.0 で GUI が動き、窓の大きさが Tk 8.6 と同じ。見た目を PNG に撮った | Tk 9.0 の `cget` の単位つきの値 (`4p`) を使う SAD の `.n` (見つからず) | なし | — |
| 8. 陽子の β の因子 (B112) | k64 は gs と全桁。v2 の U0 は古典式のちょうど 1/β³ 倍 (1.5・3・30 GeV/c で 64%・15%・0.15%) | β の因子の修正 (要判断) | なし | — |
| 9. 文法の細部 D1〜D3 | `SIN[0.1]` を構文エラー (両モード)、`45DEG` を k64 で構文エラー、MAIN の `E`・`Pi`・`Degree` の代入を使う (両モード)、末尾の演算子 `1+` を k64 で値 1 | 未定義の名前・未代入の `E` のエラーの文言 | `parser.py` | 混在 (下の表) |

## 2. 項目 1: MEASURE・TRCOD

| 欄 | 内容 |
|---|---|
| 状況 | v2.0.0 から「gs は MEASURE が無出力、TRCOD が SIGSEGV (B42) で、gs との比較は 1 件も無い」だった。今回 `MEASURE ORIGIN` (先頭の MARK) を gs で走らせると**統計が出た** |
| K64 (gs) の動作 | (a) **出力が出るのは位置が先頭の要素 (入口の MARK) のときだけ**。`MEASURE`・`MEASURE $$$`・`MEASURE QF.1`・`MEASURE M1` (MARK だけを通る位置も) は無出力で、`FFS["MEASURE"]` は `{1,{0,0,0,0,0,0},{{0,…}}}` (np = 1)。粒子の追跡 (`tracker` の `tturn`) で全粒子が落ち np ≤ 1 になる (原因は未特定: 条件・フラグ・エネルギー・EMITX を変えても同じ、大きなリング A・PF・ASTRID2 でも同じ)。(b) TRCOD も先頭の位置だけ動き、それ以外は SIGSEGV。(c) 初期分布は `ttinit` の B41 のとおり (軌道の揺らぎを `DP`·gauss から読む): py の重心が 4e-4 ほどずれる (x・px・y は 1e-5 以下)。JITTER OFF の TRCOD の重心も py が −8.25e-4 (B87 と同じ重心合わせの誤り: 他の成分は 1e-21)。(d) **MEASURE を続けて呼ぶと分布が変わる**: `NISTACK$FIXSEED@Push/Pop` が効かず種が戻らない (FIXSEED でも NOFIXSEED でも)。位置が先頭でない MEASURE (無出力) も種を進める。(e) TRCOD の粒子数は TRCOD 時の `NP` でなく、最後の `tfsetparam` (CALC・EMIT・MEASURE など) の時の値 (なければ 1000)。(f) `NP` の既定は 1000 (`P5_measure_default_np.sad`。k64 も同じ。v2 は NP の指定を要求する) |
| v2 の動作 | 追跡を通せる (v2.1.0 の C・C3)。種は戻す (K64 の設計の意図)、py の重心は 0 のまわり。変えていない |
| k64 モードの扱い | **K64 の動作の追加** (これまで k64 は MEASURE・TRCOD を無視): `measure._k64_measure`・`_k64_trcod`、`ffscmd` の `is_word`・分派・`_k64_np0` (CALC・EMIT などの時点の NP を覚える)。先頭以外の位置は無出力 (`FFS["MEASURE"]` は np = 1)、TRCOD は警告 (`SIGSEGV (B42)` と書く) を出して何もしない。種を戻さない。ORBITJITTER・PUTMEASURE・GETMEASURE・TRCOD FIXCOD/MVCOD/STNCOD/ALIGN は k64 で従来どおり無視 (FIXCOD は gs が SIGSEGV、ほかは gs でも効果が見えない) |
| 根拠 | 統計ブロック (重心・6×6 の分散行列・射影/非結合のエミッタンス・β・R1..R4・σ) が **gs と 4 桁で全一致**: FODO + スキュー四極 (結合あり) の `MEASURE ORIGIN` を GAUSS・NOTRGAUSS・DP0・NP・EMITX 等の大域変数を変えて 7 通り、LOGMEASURE、TRCOD ORIGIN 2 回 (NP の古い値の再現を含む)、大きなリング A・PF・ASTRID2 の先頭に MARK を足したもの |
| 試験 | `test_v2_P5_measure.py` 10 件: `test_measure_and_trcod_statistics_match_gs_k64` (5: NP の既定を含む)・`test_k64_py_centroid_is_the_b41_dp_jitter_not_zero`・`test_k64_measure_at_other_positions_prints_nothing_like_gs`・`test_k64_measure_consumes_the_random_seed_in_every_call_like_gs`・`test_v2_initial_distribution_is_centered_on_the_design_orbit`・`test_k64_trcod_beyond_first_element_is_refused_like_gs_sigsegv` |
| 使える範囲 | k64: 先頭の MARK の位置の MEASURE・TRCOD (初期分布と統計の検証)。v2: 追跡を通した MEASURE・TRCOD (小さなリング、PF、ASTRID2、大きなリング A は v2.1.0 の C・C3) |
| 使えない・未確認 | **v2 の追跡経由の MEASURE・TRCOD の値は gs と比べられない** (gs が追跡で落ちる)。k64 の「先頭以外は何も出ない」は先頭の位置の実験 (MARK 1〜2 個の線、4 つのリング) から決めた規則で、gs の落ちる原因は未特定。TRCOD FIXCOD などの k64 の実装は無い。`GETMEASURE`・FIFO は未確認のまま |

## 3. 項目 2: DETR >= 1 の INS

| 欄 | 内容 |
|---|---|
| 状況 | v2.1.0 の A1 (N12) は「DETR ≥ 1 (モードが入れ替わる) の INS と MEASURE は未確認」とした。強い結合のリングは DETR < 1 のままなので、**座標の 90 度回転 (COORD CHI3 = ±π/2) の対**をリングの両端に入れて、本体を入れ替わった状態にした (DETR = 1 + xyth − R1R4 + R2R3 = 1.375、R ≈ 0) |
| K64 (gs) の動作 | 回転の間は x モードが鉛直、y モードが水平。CELL CALC (周期解) も INS CALC も DETR = 1.375 の区間を正しく伝搬する。**一周の位相の整数部**: 純粋な 90 度回転の要素は `atan2(y, x)` の x < 0 で y が丸め誤差 (0 か ±1e-31) になり、位相の跳びは丸め誤差の符号で +π か −π になる (`qtwiss.f` の `dpsix=atan2(u32,aa)` など)。gs は 3 つのリング (大きなリング A・PF・ASTRID2) とも 2 つの跳びが打ち消し、NX・NY の整数部は本体の物理的なチューン (ASTRID2: NY = 5.18、PF: 9.5975、大きなリング A: 47.58)。丸め誤差の符号は計算の細部で決まるので、同じ式でも他の実装では変わる |
| 確認 (k64 と gs) | β・DETR・R・Emittance[] の Tunes・減衰率・エミッタンス・U0 が gs と一致 (3 つのリング、CELL CALC と INS CALC)。**pysad の k64 は一周の NX・NY の整数部が ASTRID2 で +1、PF の INS で +1、大きなリング A で −1 ずれる** (小数部は 8 桁で一致) |
| v2 の動作 | `lattice._exch_phase` (と numba の `perfkern._step` の `v2` 引数): 入れ替わりに**入る**要素で x < 0 かつ \|y\| < 1e-9 \|x\| の位相の跳びを +π、**出る**要素を −π と決め、対で打ち消す (丸め誤差に依らない)。3 つのリングで gs と同じ整数部 |
| k64 モードの扱い | K64 のまま (丸め誤差の符号で決まる跳び。kernel は `v2=False` で今までと同じ演算)。`perfkern.py` を編集したので `py -3.12 -m pysad.nbcache` で掃除した |
| 根拠 | gs の整数部 (本体の物理的なチューン)、v2 の Emittance[] は回転なしのリングの Emittance[] と、モードの並びを除いて 1e-6 で同じ (PF・ASTRID2。エミッタンス・減衰率・U0) |
| 試験 | `test_v2_P5_detr.py` 11 件: `test_detr_rows_and_beta_match_gs_in_both_modes` (3)・`test_phase_advance_fraction_matches_gs_k64_and_v2_integer_too` (3)・`test_emittance_with_swapped_modes_matches_gs_k64` (3)・`test_v2_emittance_is_independent_of_the_mode_exchange` (2) |
| 使える範囲 | 90 度回転で入れ替わるリングの CALC・INS・Emittance[]。v2 の整数部は純粋な回転の跳び (\|y\| < 1e-9 \|x\|) に限る |
| 使えない・未確認 | **MARK に生の DETR ≥ 1 を書いて入れ替わった入口にする INS は gs と合わない** (調べたが原因不明): gs は `DETR=1` を書いても入口の DETR を R から計算し直した値 (R1 の値が出る) にして、伝搬は DETR を書かないときと途中まで同じ (QS の後で違う)。pysad は生の DETR ≥ 1 を入れ替わった入口として使う (A1 の「入口の規格化が違う」の原因)。これを k64 で gs に合わせる修正は途中で止めた (gs の挙動を説明できないため。`DX` を書くと K64 の B40 の取り違えで 1e200 になる)。**MEASURE の入口が DETR ≥ 1 の場合**は未確認 (gs は先頭の MARK の位置だけ動く + MARK の生の DETR が gs で効かない)。非常に強い結合・回転以外のモード入れ替わり (例: ソレノイドを含む位相のゆっくりした入れ替わり) の整数部は未確認 |

## 4. 項目 3: 複素の特殊関数の k64 の挙動

| 欄 | 内容 |
|---|---|
| 状況 | 「k64 は gs が誤る点でも真値を返す (B118〜B120)」だった。どこまでが gs の誤りで、pysad に何が足りないかを、12 個の複素数 (|z| 0.05〜22、各象限) × 60 の関数・形 (732 式、`P5_specfn_cx_scan.sad`) の走査で調べた |
| 追加前の結果 | gs と 1e-9 で一致: 436 式。**pysad に足りない・誤る**: `ExpIntegralE`・`PolyGamma`・`WhittakerM/W`・`ProductLog`・`CosIntegral`・`GammaRegularized`・`Beta`・`Factorial`・`Gamma[a,z]`・複素の次数の `BesselJ` は何も返さない。`ExpIntegralEi`・`SinIntegral`・`FresnelS/C`・`HarmonicNumber`・`LegendreP` は**虚部を黙って捨てた実数**を返した (誤った値) |
| 追加 (両モード) | `pysad/specfn_cx2.py`: scipy の複素対応 (`exp1`・`expi`・`digamma`・`sici`・`fresnel`・`lambertw`・`hyp1f1`・`hyp2f1`) と、整数次の ExpIntegralE の漸化式・級数、`polygamma` (Hurwitz ゼータ経由)、複素の次数の BesselJ の級数 (\|x\| ≤ 25)、Whittaker W (`specfn_cx.hyperu_c`)。`specfn.install` が、複素数が 1 つでもあれば複素の実装に切り替える (実数は従来のまま) |
| 追加後の結果 | **641 / 732 式が gs と 1e-9**。残り 91 式は独立な基準 (scipy、級数) と比べて pysad が正しく、gs が誤る・不正確: (a) 大きな誤り: `LogGamma` (Re z < 0 で虚部の分枝が 4π 違う、4/12 点)、`Erf`・`Erfc` (Re z < 0: 例 Erfc[−1.2+.7I] gs −.0664+.1202I、真値 2.0664−.1202I。5 点)、`ExpIntegralEi` (複素引数 12 点すべて: 例 Ei[.3+.4I] gs −.629+2.409I、級数 .1595+1.389I。D.n の式で、実引数の B44 と同じ系統)、`ArcCosh` (6 点で符号が逆)、`FresnelS`・`FresnelC` (\|z\| ≥ 12 で 1e38 や符号違い)、`BesselY` (複素引数は gs が NaN、24 点)、`HypergeometricU` (整数 b で NaN、14 点: B120)、`HurwitzZeta[−20−10I,.3]` (gs 2.7e7、pysad 5e16: B118 の周辺)、`PolyLog[1.5,…]` (B119)。(b) 精度が低い: `BesselK` (\|z\| 大で 1e-9〜1e-1)、`AiryAi/Bi` (\|z\| 大)、`ExpIntegralE[1,−20−10I]` (4e-6)、`Factorial` (1e-9〜3e-8) |
| k64 モードの扱い | 追加した関数は K64 の動作の追加 (gate なし)。gs が誤る点は k64 も真値 (B118〜B120 と同じ扱い。gs の誤りの再現は K64 の `gamma.f`・`D.n` の直訳が要るので**要判断**: 上の (a) のうち `LogGamma` の分枝・`ArcCosh` の符号・`Erf/Erfc` の Re z < 0・`ExpIntegralEi` は規則が単純なので k64 で再現できる可能性がある) |
| 根拠 (独立) | `ExpIntegralEi` の級数、`E_n` の漸化式と dE_n/dz = −E_{n−1}、ψ(z+1) = ψ(z) + 1/z と ψ'(z) = Σ 1/(z+k)²、H(z) = ψ(z+1) + γ、Γ の漸化式、B(z,1) = 1/z、W e^W = z、Si・Ci の級数、Fresnel と不完全ガンマの積分 (原点から z への直線の複素の線積分)、Legendre・Whittaker M/W・複素の次数の Bessel の微分方程式の残差 (中心差分、3e-6)。使った点 \|z\| ≤ 5.1 |
| 試験 | `test_v2_P5_specfn_cx.py` 23 件 (gs と一致 2、gs の誤りの点で真値 1、恒等式・微分方程式 20)、補助データ `_p5_gscx.py` (gs の値: 一致 199 点、誤り 55 点) |
| 使える範囲 | 追加した関数: \|z\| ≲ 30 (Legendre・WhittakerW は ≲ 20、複素の次数の BesselJ は \|x\| ≤ 25)。ExpIntegralE は整数の次数 |
| 使えない・未確認 | 複素の次数の `ExpIntegralE`・`HarmonicNumber[z,r]`・`LegendreP[ν,m,z]` (m ≠ 0)、複素の a の `GammaRegularized`、`BesselY/I/K` の複素の次数、\|z\| > 30 の W・U は SADError か未確認 |

## 5. 項目 4: OYStability・OYScan

| 欄 | 内容 |
|---|---|
| 状況 | gs は `CSRSetupOY` の出力を `CSROYMatrix` に渡せない (B60: zj の head の上書き、D2) ので `OYStability`・`OYScan` が完走せず、比較が無かった |
| 方法 | `OYStability` の後半 (`oym=CSROYMatrix[…]` → `{ev,u}=Eigensystem[oym]` → `ev=Sqrt[ev]` → `gr=Max[Abs[Im[ev]]−damp Restrict[Abs[Re[ev]/(2Pi nus)],1,Infinity]]`) を、pysad v2 の `CSRSetupOY` の出力 (nr = 256、nj = 6、nphi = 3、dz = 2e-5。zj が正しい実数リスト) を**台本にそのまま書き込んで** gs に通した。`CSRHaissin`・`CSRSetupOY` の値は gs と一致が確認済み (`test_csr.py`) |
| 結果 | 行列 (18 × 18) の要素の和・絶対値の和、18 個の固有値 (実部・虚部)、成長率が gs と 1e-12 以内。安定 (aw = 1e-8): gr = −damp = −1e-4、不安定 (aw = 1e-3): gr = +0.027279622602 |
| K64 の動作 | `OYStability` は gs で wrongtype (`CSROYMatrix` の引数)。k64 の pysad は `CSRSetupOY` の出力を受けて走る (gs の wrongtype を再現しない: 変えていない。**要判断**: 再現すると k64 で OYStability・OYScan が使えなくなる) |
| 試験 | `test_v2_P5_micado_fitp_oy.py` の `test_oy_pipeline_matches_gs` (2)・`test_oy_growth_rate_values` |
| 未確認 | `OYStability`・`OYScan` の通し (CSRStability の構築・ZL の用意を含む) は gs が動かないので gs とは比べていない (MakeZL は 7 分かかる)。実用の物理 (閾値の強度) の独立な基準 (例: Vlasov 方程式の既知の解) は作っていない |

## 6. 項目 5: Micado の変種

12 セルの FODO (3 GeV、周長 52.8 m、BPM 24 台、水平・垂直ステアラ各 12 台、四極の据え付け誤差。`D2_micado.sad` のラティス)。n = 0 (最小二乗)・1・2・3・6・12。

| 変種 | gs と k64 | v2 |
|---|---|---|
| 水平 (`"X"`) | rms・最大のキック・K0 の和が 5e-6 以内。**n = 1 で補正前 6.85e-4 より悪化 (5.86e-4 は n = 0 の 7.3e-5 の 8 倍)** = B129 (K64 が選んだ n 本でなく最小二乗解の先頭 n 個を使う) を水平でも再現 | n = 1・2・3・6・12 で 4.3e-4・3.2e-4・2.3e-4・1.2e-4・7.26e-5 と単調に下がり、n = 12 は最小二乗 (7.26e-5) |
| 2 平面同時 (`{"X","Y"}`、水平ステアラだけ) | 同上 (5e-6)。垂直の rms は変わらない | 水平は単調 |
| `SVDmode->1` (垂直) | 5e-8 以内。**n = 1・2・3・6・12 の結果が互いに同じ (n = 0 とも 1e-8 の差だけ)**: `SVDmode->1` は Micado の n を無視する (K64 の `LinearSolveConditionedS` が `nx` を使わない?)。pysad も同じ | 同じ (v2 は `SVDmode->1` を直していない) |
| `MakeBump` (`Bump[sv,Condition->{"BA5","Y",1e-3},Range->{"BA2","BA8"}]`) | キック 6 本が gs と 1e-12 | 同じ |
| `CorrectOrbit[…, Bump->b]` | K64 の誤り: `correction.n` が `CorrectOrbitBump[kind,mon,str,opt,op]` と 5 引数で呼ぶが、定義は `CorrectOrbitBump[orbitkind_,bumpset_,bpm_,str_,opt_,options___]` (bumpset が抜けている) → gs は `Index 15 out of range` で `r` は評価されず、補正は何も起きない。pysad も同じエラーで同じ | 同じ |

試験: `test_v2_P5_micado_fitp_oy.py` の `test_micado_variants_match_gs_k64` (3)・`test_micado_horizontal_v2_decreases_monotonically_k64_does_not`・`test_make_bump_matches_gs_and_corrector_option_is_broken_in_both`。
未確認: `CorrectOrbitBump` を直接呼ぶ使い方 (`Bump` の組を渡す) は gs が 300 s で終わらず (pysad も 2 分以上) 比べていない。実機の BPM の誤差 (`Weight`・`Orbit->` の乱数) つきの補正は未確認。

## 7. 項目 6: FITP > 1

12 セルの FODO、`FREE QF QD SF SD`。「5 点」: DP = 0.02、`NX 3.35 5; NY 2.45 5; FITP 5`。「beta」: DP = 0.01、`NX 3.4 3; NY 2.5 3; BX 5.2 3; BY 2.5 3; FITP 3`。「one」: FITP 1 (基準)。

| 例 | gs の反復 → 最終の残差 | pysad (k64 = v2) | K1 (QF・QD) | K2 (SF・SD) |
|---|---|---|---|---|
| 5 点 | 2.0073e-5 → 1.8081e-5 (3 反復。`Overflow & closed orbit not found` 2 回で打ち切り) | 1.8557e-5 → 1.8076e-5 (3 反復 + 4 回の overflow) | .700189・−.599611 (6 桁一致) | gs 1.697215・−2.441658、pysad 1.697386・−2.441213 (1e-4) |
| beta | 9.2310e-4 → 9.0772e-4 (2 反復) | 4.1621e-2 → 9.0816e-4 → 9.0772e-4 (4 反復) | .709543・−.608680 (一致) | 1.748287・−2.529946 / 1.748324・−2.529775 (2e-5〜7e-5) |
| FITP 1 | 1.0454e-5 → 3.5e-14 | 同じ | .709545・−.608677 (一致) | −3.6e-16 / 0 (雑音) |

- 最終の残差は FITP>1 で gs と 2e-4 以内 (beta は 9.0772e-4 で全桁一致)。K2 は sextupole の谷が浅く、反復の経路の違いで 1e-4 ずれる (残差には 1e-4 しか効かない)。経路が違う原因は、gs が試行点の「Overflow & closed orbit not found」で 1〜2 反復早く止まる (pysad の `qcod` の失敗の扱いと K64 の違い)。v2.1.0 の C は FODO 24 セルで同じ傾向 (gs は最小値の手前で止まる) を記録している。5 点・beta の独立な最小値は未計算 (24 セルでは C が独立な最小化と 5e-8 一致を確認)。
- k64 と v2 は FITP で同じ結果 (最終の残差・K・反復)。
- 試験: `test_fitp_final_solution_matches_gs_k64` (3)・`test_fitp_first_iteration_residual_matches_gs_for_fitp_1`。
- 未確認: 結合・六極のあるリングの fam 点 (`MatchingAmplitude`・`NFAMP`)、大きなリング・PF・ASTRID2 の FITP>1 (オフモメンタムの closed orbit が不安定で gs が長くかかる/落ちる可能性があり、今回は時間の都合で走らせていない)。

## 8. 項目 7: Tk 9.0 の見た目

- Python 3.14.9 (専用の venv、numpy 2.5.5、scipy 1.18.1、numba 0.68.1) の Tk 9.0 で、pysad が SAD の `Tkinter.n`・`KBFrame.n` を解釈して作る窓 (Window・Button・CheckButton・TextLabel・Entry・ListBox・Scale・Canvas・Menu・KBMainFrame・KBFComponentFrame) が動く。`example/v2_next/P5_tk90_shot.py` (PIL なしの PrintWindow + 自前の PNG) で撮った 4 枚: `P5_tk90_window_button.png`・`P5_tk90_canvas.png`・`P5_tk90_entry_listbox_scale.png`・`P5_tk90_kbmainframe.png`。
- Python 3.12 (Tk 8.6) と窓の大きさが同じで、白以外の画素数の差は 0.1% 以内 (チェックボックスの描画など Tk 9.0 のテーマだけが違う。画像を目で確認した)。
- 試験: `test_v2_P5_tk90.py` (3 件、デスクトップ・venv が無い環境は skip)。D2 の記録の 25 + 90 件 (Python 3.14 の GUI・パーサ・補正など) は全試験の中で再確認していない。全試験を 3.14 で走らせることは今回もしていない。
- 未確認: gs の窓 (WSLg の X11) の撮影、Tk 9.0 の `cget` の単位つきの値 (`4p`) を数として使う SAD の `.n`、実機の高 DPI での見た目。

## 9. 項目 8: 陽子など電子以外の β の因子 (B112)

| 欄 | 内容 |
|---|---|
| 状況 | A1 は「放射の古典半径が質量によらない (B112) を v2 で直した (`re0 = elradi·elmass/amass`)。β の因子は直していない」とした。今回は低い γ (β = 0.85〜0.9995) で量を測り、k64 が gs と一致するかを確かめた |
| 条件 | 陽子 (MASS = 938.272088 MeV)、運動量 1.5・3・30 GeV/c (γ = 1.89・3.35・32.0)、12 セルの FODO (周長 45.6 m、偏向 15 度 × 12)、RF 2 MV・h = 24、`Emittance[]` の `EnergyLossU0` |
| 結果 | k64: gs と全桁 (U0 = 1.52749024e-5・1.52172596e-4・1.2651 eV)。mp/me = 1836.15 倍 (v2 との比が 1e-8 で mp/me)。**v2 の U0 は古典式 U0 = (2/3) r_p mc² β³γ⁴ Σθ²/L のちょうど 1/β³ 倍 (6e-10 で一致)**: 1.5 GeV/c で 1.641、3 GeV/c で 1.150、30 GeV/c で 1.0015 |
| 判断 | **直していない (要判断)**。電子以外で放射が意味を持つのは γ ≫ 1 (LHC の陽子など) で、そこでは 1/β³ − 1 ≈ 3/(2γ²) < 1e-3。低い γ の陽子では U0・減衰率が大きく出る (1.5 GeV/c で 64%) ので、その範囲で放射を使わない。直すには、`tradkn` の光子数 (5/(2√3) α γ θ) と光子のエネルギー分布の β 依存を、厳密な放射のスペクトル (β < 1 のシンクロトロン放射) で入れ、独立な基準 (厳密な放射の式の積分) と比べる必要がある。`crad` に β³ を掛ける近道は U0 を合わせるが、量子励起 (光子数と ⟨u²⟩) の β 依存を無視するので入れなかった |
| 試験 | `test_v2_P5_proton_beta.py` 6 件 (k64 = gs、k64/v2 = mp/me、v2 = 古典式/β³) |
| 使えない・未確認 | `TouschekLifetime`・`Intrabeam`・RADINT の `rclassic` (陽子での確認なし。v2.1.0 のまま)。重イオン (電荷 > 1) の放射の β の因子、β < 0.8 の粒子 |

## 10. 項目 9: 文法の細部 (D1〜D3)

docs/47 §(d) の 3 件 (`45DEG`・`SIN[0.1]`・未定義の名前 / MAIN の `E`・`Pi` の代入 / 末尾の `1+`) を、gs と 1 式ずつ比べた (gs は最初の構文エラーで止まるので別の台本)。`BEND B1=(L=1 ANGLE=式)` のあと `Element["ANGLE","B1"]`。gs の値は `P5_README.md`。

| 形 | gs | k64 (今回) | v2 (今回) | 変更 |
|---|---|---|---|---|
| `SIN[0.1]` | 構文エラー | 構文エラー | 構文エラー | pysad は値 0 を入れていた (誤り): **両モードで構文エラー**にした |
| `45DEG` (数と単位の間に空白なし) | 構文エラー (`Undefined keyword: DEG…`) | 構文エラー | 受ける (45 度) | k64 だけ gs と同じにした。v2 は寛容さを残した (空白の位置を `Token.pos/end` で判定) |
| MAIN の `E=0.5;` `Pi=3;` `Degree=0.1;` のあとの `ANGLE=E`・`Pi/10`・`Degree` | 利用者の値 (.5・.3・.1) | .5・.3・.1 | .5・.3・.1 | **両モードで直した** (pysad は SADScript の定数 e・π を使っていた)。代入は別の辞書 (`env._main_user_vars`) に入れ、要素定義の式だけが見る (FFS の SADScript の `E` は変わらない: gs の実測 `EE 2.718…`) |
| 末尾の演算子 `1+`・`2*3+` | 構文エラーにせず**値 1** (`2*3+` も 6 でなく 1) | 値 1 | 構文エラー | k64 だけ gs と同じ (K64 の誤りの候補)。v2 は入力ミスとして受けない |
| 未定義の名前 `ZZZ`・`ZZZ*2`・`pi/10` | 構文エラー | エラー (文言が違う) | 同じ | 変更なし |
| MAIN で代入していない `E` (`ANGLE=E`) | 構文エラー | 定数 e を使う (未解決) | 同じ | 変更なし |

試験: `test_v2_P5_syntax.py` 13 件。`parser.py` の `Token` に `pos`・`end` を足した (他の動作は変わらない)。

## 11. K64 の新しい誤りの候補 (B 番号は付けない。確かめたこと・足りないこと)

| 仮の名前 | 内容 | 確認 | 足りないこと |
|---|---|---|---|
| P5-1 | `MEASURE`・`TRCOD` の追跡 (`tracker` の `tturn`) が、先頭の要素より後ろを通ると全粒子を落とす。MEASURE は無出力、TRCOD は SIGSEGV (B42 の更新: 先頭の MARK の位置だけ動く) | gs の 10-19 ビルドで、FODO・大きなリング A・PF・ASTRID2、エネルギー・エミッタンス・フラグ (NORFSW・NOCOD・RAD・CALC6D) を変えて再現 | 落ちる原因 (`tturn` の入力の設定のどこか) はソースの読みで未特定 |
| P5-2 | `NISTACK$FIXSEED@Push/Pop` が MEASURE で効かず、FIXSEED でも種が戻らない (MEASURE を続けると分布が変わる) | gs (FIXSEED・NOFIXSEED、2 回・3 回の MEASURE) | `NISTACK.n` の `GetSeed[] := SeedRandom[]` の戻り値が種になっていない可能性 |
| P5-3 | 入れ替わり (DETR ≥ 1) をまたぐ純粋な 90 度回転の位相の跳び ±π が丸め誤差の符号で決まり、一周の NX・NY の整数部が不定 (gs と pysad k64 で ±1) | gs と pysad の 3 リング。v2 が対で打ち消す規則で gs と一致 | gs が 3 つともたまたま打ち消した可能性 (丸め誤差の符号は gs の計算順で決まる) |
| P5-4 | TRCOD の粒子数が最後の `tfsetparam` の `NP` (CALC・EMIT・MEASURE の時点) で、TRCOD 時の NP でない (なければ 1000) | gs (`NP=200; CALC; TRCOD ORIGIN` で 200、`CALC; NP=200; TRCOD` で 1000) | 他の `tfsetparam` を呼ぶ命令の全リスト |
| P5-5 | `CorrectOrbit[…,Bump->b]` が `CorrectOrbitBump` を引数 1 つ足りない形で呼ぶ (bumpset の抜け)。gs は Index エラー | gs と pysad | — |
| P5-6 | MAIN の要素定義の式の末尾の演算子 (`1+`・`2*3+`) が構文エラーにならず値 1 になる | gs | — |
| P5-7 | `CorrectOrbit` の `SVDmode->1` が Micado の n を無視する (n ≥ 1 の結果が互いに同じ) | gs と pysad (k64) | `LinearSolveConditionedS` のソース |
| P5-8 | 複素の特殊関数: LogGamma の分枝 (Re z < 0)・Erf/Erfc (Re z < 0)・ExpIntegralEi (複素引数)・ArcCosh (複素の符号)・FresnelS/C (\|z\| ≥ 12)・BesselY (複素引数、NaN)・BesselK/AiryAi/Bi/ExpIntegralE/Factorial の精度 | gs と scipy・級数・恒等式 (§4) | 原因 (D.n・gamma.f の式) の個別の特定 |
| P5-9 | MARK の `DETR` キーワードが INS の入口で期待どおりに効かない (入口の DETR が R から計算し直され、R1 の値が出る。伝搬も途中から違う) | gs (R と DETR を変えた 6 通り) | 原因 (MARK の値の添字の取り違え: B40 の周辺か) |

## 12. 追加した試験・変えた既存の試験

- 追加: `test_v2_P5_measure.py` 10、`test_v2_P5_detr.py` 11、`test_v2_P5_micado_fitp_oy.py` 12、`test_v2_P5_specfn_cx.py` 23、`test_v2_P5_proton_beta.py` 6、`test_v2_P5_syntax.py` 13、`test_v2_P5_tk90.py` 3 = **78 件**。補助 `_p5_util.py`・`_p5_gscx.py`。
- 変えた既存の試験: `test_v2_C_measure.py::test_k64_mode_measure_and_trcod_are_ignored_as_before` → `test_k64_mode_measure_and_trcod_at_the_end_do_nothing_like_gs` (k64 の MEASURE・TRCOD が gs の実動作になったため: 終端の位置は統計なし・SIGSEGV の警告・Twiss 不変・`FFS["MEASURE"]` は np = 1)。ほかは変えていない。
- 変えたコード: `pysad/measure.py` (k64 の MEASURE・TRCOD、`Tracker.reset` の v2 だけ種を戻す)、`pysad/ffscmd.py` (MEASURE・TRCOD・LOGMEASURE の k64 の分派、`_k64_np0`)、`pysad/lattice.py`・`pysad/perfkern.py` (`_exch_phase`・kernel の `v2` 引数。`pysad.nbcache` で掃除)、`pysad/parser.py` (`Token.pos/end`、D1〜D3)、`pysad/specfn.py`・`pysad/specfn_cx2.py` (新)・`pysad/sadscript.py` (複素の特殊関数、Factorial の複素)。改行は元のまま (`ffscmd.py`・`perfkern.py`・`parser.py`・`specfn.py`・`sadscript.py` は CRLF、`measure.py`・`lattice.py` は LF、新規の `specfn_cx2.py` は CRLF)。

## 13. 全試験・回帰

全試験 (`PYSAD_SAD_PACKAGES` つき、旧式 3 本を除く): 1644 通過・9 skip・1 失敗 (上の `test_v2_C_measure.py` の 1 件。直して 19 件通過を確認。P5 の 78 件は最終に 56 s で全通過)。基準の 1512 通過・8 skip 前後に対し、他の担当の追加分と P5 の 78 件が加わった数。回帰 `QUICK=1`: 50 ベンチ FAIL 0 (129 s)。他の担当が同時に走っていたので全試験は 28 分かかった。

## 14. 付けた B 番号 (v2.2.0 の文書化で。docs/21)

- §11 の P5-1 = **B151**、P5-2 = **B152**、P5-3 = **B153**、P5-4 = **B154**、P5-5 = **B155**、P5-6 = **B156**、P5-7 = **B157**、P5-8 = **B158**、P5-9 = **B159**。
- B42 (MEASURE・TRCOD)・B112 (β の因子)・B118〜B120・B129 の更新は各行に追記した。
- 要判断 (複素の特殊関数の k64 の値、`CSROYMatrix` の wrongtype、β の因子、文法) は docs/47 の v2.2.0 の (d)。
