# 48. 周波数マップと局所運動量アクセプタンス (LMA) の解説

最終更新: 2026-11-20 (v2.3.0、S1 の検証を §9 に追加)

関連: [docs/49 (共鳴駆動項・ADTS・クロマティシティ。v2.4.0 の J²・δ 依存は §13)](49_resonance_terms_adts.md)。**入門 (何をどう求めて何と比べたか) は [docs/50](50_adts_rdt_primer.md)**。

K64 にない pysad の追加機能。1 つの位置の初期条件の格子から **周波数マップ** (チューンの拡散の地図) を作る機能と、位置 s ごとの **局所運動量アクセプタンス**
(δ+(s)・δ−(s)) を求める機能。どちらも SAD の追跡 (`TrackParticles` と同じ損失判定の高速追跡、numba) の上に作り、SAD の台本から使える関数
(`FrequencyMap[]`、`MomentumAcceptance[]`、`TouschekLifetimeLMA[]`) と Python の API の両方で使える。K64 の MAIN の TRACK には Fourier 解析 (FOURIER フラグ) があるが、
FFS からは使えない (docs/history/work_main_track.md)。チューンの取り出しは [`pysad/naff.py`](../pysad/naff.py) (NAFF、v2.1.0) を使う。

## 1. 何ができるか

| 機能 | 状況 | 実装 |
| --- | --- | --- |
| 周波数マップ (x–y、x–δ、y–δ の格子、前半・後半のチューンの差 D = log10 √(Δνx² + Δνy²)、失われた点の区別、表・図) | 対応 (gs の追跡と座標が 3e-7 以内、D が 4 桁一致) | `freqmap.frequency_map` |
| 共鳴線 (m νx + n νy = p) の描画、作動点に近い共鳴線の表 (整数・半整数を必ず含む) | 対応 | `freqmap.resonance_lines`・`resonance_check`・`line_distance` |
| LMA (位置ごとの δ+・δ−、粗い走査 + 二分法、放射・RF つき、RF バケットから出た粒子も損失) | 対応 (gs と損失判定が一致) | `lma.momentum_acceptance` |
| 周期的なリングの 1 周期 (セル) だけの計算と、別のセルでの確認 | 対応 | `trksetup.find_periodicity`・`select_locations` |
| Touschek 寿命への受け渡し (s, δ+, δ−)、位置ごとの受容での寿命 | 対応 (既存の TouschekTable を使う) | `LMAResult.touschek_input`・`lma.touschek_lifetime_local` |
| SAD の台本の関数 | 対応 (pysad のみ。gs では走らない) | `q1cmd.py` |
| FLUC (量子励起) つきの LMA | 対応 (`TrackSetup.fluc = True`・`seed`。複数の種の平均・分散で使う。2 周期では効果は二分の幅以下、§4.5 (f))。周波数マップは S2 §10.3 | `lma.momentum_acceptance` |
| 誤差つきのリング (位置のずれ・ロール、全周の追跡) | 周波数マップは走る (S1、ASTRID2、§9.4)。LMA も走る (L2、ASTRID2 level 1・3 種。v2.4.0 A4: level 5 = 結合の見積もり |c−| 0.107 が差の共鳴までの距離 0.04 を超える強い結合・閉軌道 4.8 mm でも走る。大きなリング A は `docs/history/v24_A.md` §A4)。補正つき (COD 補正後) は閉軌道を Newton 法で求めるので補正の状態をそのまま使えるが、台本から取り込む経路は未検証 | `fmapres.apply_alignment_errors` |
| δ を含む周波数マップ (x–δ、y–δ) | 対応 (既存の `frequency_map` の dp の格子。S1 で ASTRID2・KEK-LS を確認、§9.3) | `freqmap.frequency_map` |
| 共鳴駆動項 (RDT)・ADTS からの共鳴の領域の予想 | 対応 (1 次の RDT が駆動する共鳴だけ。S1、§3.5・§9.1) | `fmapres.ResonancePredictor` |

## 2. 共通の部分 (`pysad/trksetup.py`)

`TrackSetup` が FFS の状態 (要素の並び、運動量、RF、RAD・RFSW フラグ、`Emittance[]` の閉軌道と同期位相、LOSSAMPL・LOSSDZ) を 1 つにまとめる。

```python
from pysad.trksetup import TrackSetup
st, ffs = TrackSetup.from_source(lattice_text, 'RING', rad=True, rfsw=True, alost=0.01)   # USE RING; CELL; CALC; Emittance
st = TrackSetup.from_ffs(ffs)                                                              # 台本で USE・CALC したあとの FFS オブジェクトから
st2 = st.copy(rad=False, rfsw=False)                                                       # 周波数マップ用
```

- **閉軌道**: 位置 1 の閉軌道は `Emittance[]` (RADCOD) の値から Newton 法で求め直す (RFSW ON・放射つきは 6 次元、RFSW OFF は δ を決めた 4 次元)。ほかの位置へは 1 粒子を運ぶ
  (`closed_orbit(locs)`)。1 周後との差を `conditions()['orbit_residual']` に記録する (例: ASTRID2 で 2e-18、大きなリング A で 4e-16)。
- **損失の判定**: SAD の `TrackParticles` (`track6n.TrackRun`) と同じ。tapert は最初の要素・15 要素ごと・区間の終わりで、|x|,|y| < LOSSAMPL (正方形)、|z| ≤ LOSSDZ、APERT 要素。
  - `setup.track(loc, coords, nturn)` の区間の分け方は、`TrackParticles[{loc, beam}, loc, 1, nturn]` の 1 回の呼び出しと同じ (loc→nlat、(nturn−1) 回の 1→nlat、1→loc)。
  - 周ごとの座標を記録するとき (周波数マップ) は 1 周ずつ loc から回す (1 周ずつ `TrackParticles[{loc, beam}, loc]` を呼ぶのと同じ)。
  - **`check_every = k`** (`TrackSetup` の既定は 0 = SAD と同じ。**LMA (`momentum_acceptance`・`MomentumAcceptance[]`) の既定は 1**、§4.5): k > 0 なら k 要素ごとに区間を分けて判定する。SAD の判定は 15 要素に 1 回なので、
    アパーチャが限界のとき、判定点のあいだの最大のずれを見逃し、δ の限界を大きく見積もる (ASTRID2、±10 mm、位置 17 で δ+ = 1.05% (k = 0) に対し 0.83% (k = 1)。線形の見積もり 0.865%)。
- **周期性**: `find_periodicity(els, rows)` は、マーカー・モニター・空洞を除いた要素の並び (種類と数値の属性) が同じ列の n 回の繰り返しになる最大の n を返す
  (ASTRID2 で 6、`kek-ls_v03_12.lat` で 20、24 セルの FODO で 24。セルの長さが全部同じかも確かめる)。誤差・補正は見ない。セルの中の鏡映対称は見ない。
- **位置の選び方**: `select_locations(setup, ls, le, step, ds)`: 要素の境界 (要素の入口、1 始まり) を、範囲・k 要素おき・間隔 ds [m] 以上で選ぶ。同じ s の位置 (長さ 0 の要素の列) は 1 つだけ。
- 並列の子プロセスは `lower_priority()` で優先度を下げる (Windows: BELOW_NORMAL)。

## 3. 周波数マップ (`pysad/freqmap.py`)

### 3.1 アルゴリズム

位置 loc の閉軌道に、x・y (px = py = 0、z = 0) と δ の格子点を足した粒子を `nturn` 周追跡し、位置 loc の (x, px, y, py) を周ごとに記録する。
前半 (nturn/2 周) と後半の各面の信号 z = x̂ − i p̂ (x̂ = x/√β、p̂ = (αx + βpx)/√β) に NAFF (Hann 窓、`naff.tune_xp`) をかけてチューン ν1・ν2 を求め、

  D = log10 √( (νx2 − νx1)² + (νy2 − νy1)² )

を出す (チューンは端数 0〜1、差は −0.5〜0.5 に畳む)。粒子が `nturn` 周の前に失われた点は、`lost_turn > 0` (落ちた周)、D = nan とし、図では黒い × で区別する。

- 規格化 (α, β): 既定 (`norm='data'`) は、その半分の追跡の標本の 2 次のモーメントから (x, px) の楕円の α・β を求める。δ ≠ 0 の軌道の光学関数は位置 loc の Twiss (δ = 0) と違い、
  規格化がずれると NAFF の周波数に 1e-6 の程度の系統誤差が入るため (δ = 0.5% で D の下限が −5.8 → −9.6 に下がった)。`norm='twiss'` は位置 loc の Twiss の値を使う。
- δ ≠ 0 の点: RFSW OFF は δ の閉軌道 (4 次元、Newton 法) の上に置く。RFSW ON は線形の分散 (η, η') の上に置く。x・y は δ の軌道からのずれ。
- D の下限: 線形の運動 (振幅 1 μm) で、δ = 0 のとき D ≈ −11 (FODO)。NAFF の精度で決まる。実用上は D < −6 を規則的とみなす。

### 3.2 Python の API

```python
import numpy as np
from pysad.freqmap import frequency_map, resonance_check
st = TrackSetup.from_source(src, 'RING', rad=False, rfsw=False, alost=0.1)[0]
r = frequency_map(st, 2, x=np.linspace(5e-4, 18e-3, 41), y=np.linspace(5e-4, 14e-3, 31), nturn=1024, workers=1)
r.save('fmap.table.txt'); r.plot('fmap.png', order=4)          # 表、図 (チューン平面 + 初期条件平面)
r.D, r.nux1, r.nuy1, r.nux2, r.nuy2, r.lost_turn, r.x, r.y, r.dp
resonance_check(*st.tunes_full(), order=4)                      # 作動点に近い共鳴線の表 (整数・半整数を必ず含む)
```

- x, y, dp は数または配列。配列が 2 つなら 2 次元の格子。`nturn` は偶数にする (奇数なら 1 引く)。`workers` は子プロセスの数 (既定 1)。
- 表 (`save`): 条件のコメント行のあと、列 `x y dp nux1 nuy1 nux2 nuy2 D lost_turn`。`FMapResult.load_table(path)` で読み直せる。
- 図 (`plot`): 左 = チューン平面 (νx1, νy1) を D で色づけして order 次以下の共鳴線を重ねる (点の 2%〜98% を窓にする。外れた点は見えない)、右 = 初期条件の平面を D で色づけ、失われた点は黒い ×。
  カラーバーは既定で D = −8 〜 −2。
- v2.3.0 S2 の追加 (§3.6、§10): `frequency_map(..., scaling=True, warn=True)`。`r.D_half` (最初の nturn/2 周の D)・`r.slope` (周回数への依存の傾き)、`r.D_norm` (1024 周に換算した D = D + log2(nturn/1024))、
  `r.chaotic(thr=−5)` (周回数に依らないカオス的な点の判定 = D_norm > thr)、`r.warnings` ([(コード, 文)])。信頼できない条件では `FMapWarning` が出る。表の末尾に列 `D_half` と、警告の行 `# warning = ` が入る (古い表も読める)。
  閉軌道が求まらない (誤差つきで線形に不安定) ときは `RuntimeError`。誤差つきの TrackSetup は `fmapres.apply_alignment_errors(st, seed, quad_dxy=…, streams='type')`、統計・フィードダウンの見積もりは [`pysad/fmapcheck.py`](../pysad/fmapcheck.py)。

### 3.3 SAD の台本から

```
FFS;
USE RING; CELL; CALC;
NORFSW; NORAD;
LOSSAMPL = 0.1;
f = FrequencyMap[2, 512, X->{0.0005, 0.018, 11}, Y->{0.0005, 0.014, 8}, Output->"fmap.txt", Plot->"fmap.png"];
Print[Diffusion /. f];
```

`FrequencyMap[位置, nturn, X->…, Y->…, DP->…, Window->1, Workers->1, Output->"file", Plot->"file.png", Order->4]`。X・Y・DP は数 (固定)、`{min, max, n}`、または値のリスト。
戻りは `{X->…, Y->…, DP->…, NuX->…, NuY->…, NuX2->…, NuY2->…, Diffusion->…, LostTurn->…}` (点ごとのリスト。失われた点の Diffusion は −99)。追跡の条件: **既定は RFSW OFF・RAD OFF・FLUC OFF** (v2.3.0 S2。現在のフラグは使わない。FFS に入ると RFSW は ON なので、`NORFSW;` を書き忘れると D > −5 の割合が 0.03 から 0.57 になるため。現在のフラグを使うときは `UseFlags->True`)、LOSSAMPL・LOSSDZ は現在の値。警告は ` FrequencyMap: WARNING [コード] …` で出力される (§3.6)。
例: [`example/v2_3/Q1_script_commands.sad`](../example/v2_3/Q1_script_commands.sad)。

### 3.4 共鳴線との一致

作動点の近くの共鳴線 (距離はチューン平面の距離 |m νx + n νy − p| / √(m² + n²)):

| 例 | 整数 νx = p | 整数 νy = p | 半整数 2νx = p | 半整数 2νy = p | 最も近い次数 4 以下の線 |
| --- | --- | --- | --- | --- | --- |
| ASTRID2 (νx = 5.18、νy = 2.14) | 0.180 | 0.140 | 0.180 | 0.140 | νx − νy = 3 (差の共鳴、距離 0.028) |
| 大きなリング A (νx = 47.58、νy = 17.62) | 0.42 | 0.38 | 0.08 | 0.12 | νx − νy = 30 (距離 0.028)。3 次は 3νy = 53 (0.047)・νx + 2νy = 83 (0.080)・3νx = 143 (0.087) |

大きなリング A は KEK-LS (ESRF 型の HMBA。ESRF の Pantaleo Raimondi 氏と Simone Liuzzo 氏の協力で設計した)。整数・半整数の共鳴は、どちらの作動点でも距離 0.08 以上離れている (大きなリング A の 2νx = 95 の半整数が最も近く 0.08)。

作動点の整数・半整数の共鳴線は、どちらの例でも拡散の大きい領域 (格子の範囲のチューン) の外にある (ASTRID2: 格子のチューンは νx = 0.18〜0.25、νy = 0.14〜0.23 で、整数 0 と半整数 0.5 から遠い。
大きなリング A: 大部分が νx = 0.56〜0.79、νy = 0.60〜0.76 (数点は外れる) で、半整数の νx = 0.5 は格子の左端の外)。**振幅が大きくなるとチューンが動いて、次の共鳴線に近づく**:

- ASTRID2 (位置 2、41 × 31、x = 0.5〜18 mm、y = 0.5〜14 mm、1024 周): 失われた点 0。D の範囲 −11.45〜−1.70。D > −4 の 21 点のうち 10 点は νx − νy = 3 (差の共鳴。次数 2、同じ線が 2νx − 2νy = 6 の次数 4)、8 点は 5νy = 11 (次数 5)、
  2 点は 4νx = 21 (次数 4)。νx − νy = 3 の線の上では、x の信号に y のチューンの成分が現れる (νx = νy にロックされた点が 377 点。D は −7.5 の程度の規則的な運動で、島に捕らえられている)。
  次数 4 以下の共鳴線から距離 0.003 以内の点 (379 点) とそれより遠い点 (892 点) で、D > −5 の割合は 0.061 / 0.028 (2.2 倍)。
- 大きなリング A (位置 3、41 × 31、x = 0.2〜6 mm、y = 0.2〜8 mm、1024 周): 失われた点 460 (x ≳ 4.2 mm)。D の範囲 −10.80〜−1.09。D > −4 の 62 点のうち 32 点は νx − νy = 30 (差の共鳴)、
  10 点は νx + 2νy = 83、7 点は 3νy = 53、3 点は 3νx = 143、3 点は νx − 4νy = −5 (最も近い次数 5 以下の線)。図 (`Q1_fmap_kekls.png`) では (νx, νy) = (2/3, 2/3) に集まる次数 3 の線 (3νx = 2、3νy = 2、νx + 2νy = 2 の端数の線)
  の上に D の高い点が並ぶ。次数 4 以下の共鳴線から 0.003 以内の点 (507 点) とそれより遠い点 (304 点) で、D > −5 の割合は 0.195 / 0.082 (2.4 倍)。

この統計 (近い/遠いの比 2 倍) は図の印象より弱い。格子のチューンの広がり (νx で 0.2〜0.3) のなかに次数 4 以下の線が密に通っていて、点の約 6 割が「近い」側に入るため。結論は「D の大きい領域が、
ラティスの対称性が許す共鳴線 (差の共鳴、3 次、5 次) の上に載る」であり、個々の線が D を決める強さ (駆動項) の定量的な確認は [`docs/49_resonance_terms_adts.md`](49_resonance_terms_adts.md) (共鳴駆動項) の仕事。

### 3.5 RDT・ADTS から予想する共鳴の領域 (`pysad/fmapres.py`、v2.3.0 S1)

```python
from pysad.fmapres import ResonancePredictor, coverage_stats, apply_alignment_errors
pr = ResonancePredictor.from_lattice(lat, loc, rows)      # lat = Lattice、rows = lat.twiss_table(delta=δ)。RDT と ADTS を計算する
m = pr.predict(r.x, r.y, nu=None)                          # {'any': 位相つきで島に捕らえられる点, 'envelope': 位相を問わない包絡, 'by_harmonic': {調和: ...}}
s = coverage_stats(r.D, m['any'], r.grid, thr=-5.0)        # D > thr の点が領域に入る割合 (cov)、面積 (area)、領域を格子の上でずらした対照 (shift_*)、p 値
st2, cnt = apply_alignment_errors(st, seed=1, quad_dxy=5e-5, quad_rot=2e-4, sext_dxy=5e-5)   # 要素ごとの位置のずれ・ロールの誤差つきの TrackSetup
```

式 (docs/49 の RDT h_jklm と ADTS の行列 A = ∂νi/∂Jj): 調和 a = (j−k, l−m) = g(m0, n0) の 1 次の項 V = 2|F|cos(a·φ + arg F)、F = Σ h_q (2Jx)^((j+k)/2) (2Jy)^((l+m)/2)。
共鳴 m0 νx + n0 νy = p/g の方向の作用 I を取ると、Δ0 = m0 νx(J) + n0 νy(J) − p/g は I の一次式 (傾き G0 = (m0,n0)·A·(m0,n0)ᵀ) で、振り子の式から
島の幅 W = √(2|G0|V/π) (チューン平面)、捕らえられる条件 πΔ0²/|G0| < V(1 − sgn(G0) cos θ)。V は初期条件の J で評価する。ν(J) は線形の ADTS ν0 + A·J か、追跡のチューン。
normal 6 極では調和は 3νx、νx、νx ± 2νy の 4 つだけ。**νx − νy などの 2 次の共鳴 (6 極の 2 次の項)、4 次以上の共鳴の強さの解析式はない**ので、予想できない (S1、§9.1)。
共鳴線の判定 (`lines_distance`・`other_mn`)、図 (`plot_with_prediction`) も同じファイル。

**S3 の作り直し (v2.3.0、[`pysad/fmappred.py`](../pysad/fmappred.py)、`ResonancePredictor.from_lattice(..., model='s3')`、導出・限界は docs/49 §12)**: 駆動項に 4 次の h (6 極の 2 次ほか、`rdt2nd.rdt4`) と非共鳴の 3 次の補正 (`rdt4_normal_form` の F4')、
ν(J) は解析の線形の ADTS (R6)、幅は上の振り子の式 (J ごとの G0)。点ごとに連続のスコア log10(W/|Δ0|)、包絡に入る調和の数、信頼度 (0/1/2) と全体の警告を返す。**追跡に合わせた当てはめの項は足していない。**
結果 (D > −5 を正例、保存済みの 41×31): FODO の 3νx は AUC 0.94 (振幅だけ 0.42)、ASTRID2 MPW は 0.79 (0.73)、**ASTRID2 は 0.68・KEK-LS は 0.59 (振幅だけ 0.78 より低い) で判別力がない**。
予想した 2νx−2νy の島の中は測ったチューンが線に固着している (ASTRID2 0.93、外側の輪 0.21) が、島の中は規則的で D は小さく、D の大きい点は島の縁・重なり・次数 5 以上の共鳴にある。
予想できる条件は小振幅・孤立した低次 (3・4 次) の共鳴。できない条件は 5 次以上、共鳴の重なり (島の幅が線の間隔と同程度)、大振幅 (J² の項が未実装)、δ ≠ 0。詳細は [`example/v2_3/FMAP_summary.md`](../example/v2_3/FMAP_summary.md) §13。

### 3.6 信頼して使うための条件と警告 (`pysad/fmapcheck.py`、v2.3.0 S2)

周波数マップの D は、条件によっては拡散の指標でなくなる。`frequency_map(warn=True)` (既定) は、次の条件で `FMapWarning` を出し、`FMapResult.warnings` ([(コード, 文)])・表の `# warning = `・条件の `warnings` にも入れる
(`warnings.simplefilter('ignore', FMapWarning)` か `warn=False` で消す。台本 `FrequencyMap[]` は ` FrequencyMap: WARNING [コード] …` を出力)。根拠の測定は §10。

| コード | 条件 | 理由 |
| --- | --- | --- |
| `nturn_short` | nturn < 512 | NAFF の精度。D の周回数への換算 (下) が 256 周では成り立たない |
| `window_none` | window = 0 | D の下限 −7.95、D > −5 の割合が Hann の 5 倍 |
| `rfsw` | RFSW ON | z が戻らず δ が変調される。ASTRID2 で D > −5 が 0.03 → 0.57 |
| `rad`・`fluc` | RAD ON・FLUC ON | D が減衰・乱数の指標になる (0.03 → 0.99)。RF なしの RAD ON は δ が下がり続ける |
| `tune_integer`・`tune_lowres` | 作動点が整数・半整数の 0.03 以内、2 次の和・差の共鳴の 0.02 以内、3 次の 0.01 以内 | 線形の運動が不安定に近い、結合・低次の共鳴で D が大きい領域 |
| `errors`・`coupling`・`unstable`・`orbit_large` | 要素にずれ・ロールがある、\|c−\| の見積もりが差の共鳴までの距離の半分超、閉軌道が求まらない・不安定、閉軌道のずれが LOSSAMPL の 20% 超 | 周期性なし・位置固有、D > −5 が誤差なしの 7〜12 倍 (§10.1) |
| `orbit` | 閉軌道の Newton 法が収束しない (残差 > 1e-9) | |
| `lost_many`・`lost_all` | 失われた点が 20% 超・全部 | D の統計は残った点だけ |
| `naff_fail`・`tune_edge` | 残った点で D が nan、チューンが 0・0.5 の 0.01 以内の点が 2% 超 | NAFF の失敗、分解能と同程度 |

**周回数に依らない指標**: 規則的な点の D は周回数が 2 倍になると約 1 下がる (NAFF の誤差が N^−3.9)。カオス的な点の D は周回数にほとんど依らない。D > −5 は 1024 周の目安なので、別の周回数の結果は
`r.D_norm` (= D + log2(nturn/1024)、1024 周に換算した D) で比べ、`r.chaotic()` (D_norm > −5) で判定する。512〜4096 周で使える (256 周は不可)。点ごとの傾き `r.slope` は診断用 (判定には揺れが大きい)。

**設定の既定** (安全な値): NORFSW・NORAD・NOFLUC、1024 周、Hann 窓、規格化は標本の楕円、台本 `FrequencyMap[]` もこれ。RFSW ON・RAD OFF の同期位相は 0 (放射損失がないので正味の加速なし。gs の CALC のあとの TrackParticles と同じ)。
並列でない (親プロセスで走る) ときも、TrackSetup が持つ `bend_k1_edge`・`bend_k1_def` (`compat.set_bend_k1_edge`・`set_bend_k1_def`) に合わせて追跡する。

## 4. 局所運動量アクセプタンス (`pysad/lma.py`)

### 4.1 アルゴリズム

位置 loc (要素の入口) の閉軌道 (x, px, y, py, z, δ) に δ = ±Δ を足した粒子を、放射 (FLUC なし) と RF つきで `nturn` 周追跡し、失われないΔの限界を求める。位置ごと・符号ごとに:

1. Δ = 0 と Δ = dpmax·k/ncoarse (k = 1..ncoarse、既定 12) を同時に追跡する。外向きに最初に失われた点を挟む区間 (lo: 残った最大の Δ、hi: 失われた最小の Δ) を作る。全部残れば Δ = dpmax ('sat')。
2. 区間を nsub + 1 等分した nsub 点 (既定 3。**nsub = 1 なら普通の二分法**) を追跡し、最初に失われた点とその手前の残った点で区間を狭める。hi − lo ≤ tol (既定 1e-4) まで繰り返す。
3. 結果は lo (確かに残った δ)。真の限界は (lo, hi) の中 (幅 ≤ tol)。

- 外向きの最初の損失を限界とする。走査の点のあいだの狭い損失の窓は `verify=n` で探せる (§4.5)。その外側の残る島は見ない (保守的)。
- **周回数**: 既定はシンクロトロン振動 `nsynch` (既定 2) 周期分、`ceil(nsynch/|νs|)` (νs は `Emittance[]` のシンクロトロンチューン)。例: ASTRID2 は νs = 0.00152 で 1312 周、大きなリング A は νs = 0.00408 で 491 周。
  周回数を直接 `nturn` で指定して 1 周期分より短いと警告する。RF なし (RFSW OFF) は既定 512 周。バケットで決まる限界の周回数への依存は §4.5。
- **RF バケット**: `bucket=True` (既定、RFSW ON のとき) は、|z| がバケットの長さの半分 (π + 2|φs|)/k_rf (φs = trf0·k_rf = 同期位相。安定固定点の周り ψ = −π−2φs … π−2φs の外側の端) を超えた粒子も失われたとする
  (内部では LOSSDZ をこの値にする)。損失の理由は `reason_plus`・`reason_minus` に 'rf' (|z|)・'aper' (横の |x|,|y| または APERT)・'sat' (走査の上限まで残った) で記録する。
  バケットの外へ出る粒子は、位相の滑りが始まるまでに時間がかかるので、バケット高さの付近の限界は周回数に依存し、`nsynch` を増やすとバケット高さに近づく。
- **アパーチャ**: `alost` (既定は LOSSAMPL。±alost の正方形)・APERT 要素。判定の間隔は `check_every` (**既定 1 = 要素ごと**。0 は SAD の tapert のまま = gs と同じ判定。§2、§4.5)。
- **dpmax** の既定: RF あり = バケット高さの 1.25 倍、RF なし = 0.05。
- **並列**: 位置ごとに独立なので子プロセス (`workers`、既定 2、優先度を下げる) に分ける。`workers=1` は親プロセスで順に。コア数 − 1 を超えない。

### 4.2 Python の API

```python
from pysad.trksetup import TrackSetup, find_periodicity, select_locations
from pysad.lma import momentum_acceptance, touschek_lifetime_local
st = TrackSetup.from_source(src, 'RING', rad=True, rfsw=True, alost=0.01)[0]
per = find_periodicity(st.els, st.rows)                                    # セルの数、各セルの最初の要素の位置、セルの長さ
locs = select_locations(st, per['starts'][0], per['starts'][1] - 1)        # 最初のセル (1 周期) の要素の境界
r = momentum_acceptance(st, locs, workers=2)                               # nturn の既定 = シンクロトロン振動 2 周期分、check_every の既定 = 1
r.save('lma.table.txt'); r.plot('lma.png')
s, dpp, dpm = r.touschek_input()                                           # (s, δ+, δ−) (δ− は負)
tt = touschek_lifetime_local(touschek_table, r, setup=st, period=per['cell_length'], s0=per['s_starts'][0])['lifetime']
```

`LMAResult`: `loc`、`s`、`name`、`dp_plus`、`dp_minus` (負)、`hi_plus`・`hi_minus` (最初に失われた δ)、`reason_plus`・`reason_minus`、`cond` (使った条件の辞書)、
`nwin_*`・`win_*` (verify の窓)、`tdep_*` (周回数の依存)、`warnings` (警告の文字列のリスト。`cond['warnings']` にも入る)、`flags(i)` (位置 i の旗)。
表 (`save`) は条件のコメント行 + 警告の行 (`# WARNING:`) + 列 `loc s name dp_plus dp_minus hi_plus hi_minus reason_plus reason_minus flags` (flags: `W±` 窓、`T±` 周回数の依存、`S±` 走査の上限、`0±` δ = 0 で失われた、なければ `-`)。

### 4.3 SAD の台本から

```
FFS;
USE RING; CELL; CALC;
RFSW; RAD; RADCOD; NOFLUC;
Emittance[];
INTRA; PBUNCH = 1e9; Emittance[];             ! TouschekTable を作る
LOSSAMPL = 0.01;
m = MomentumAcceptance[Cell->True, Distance->0.9, Workers->2, CheckEvery->1, Output->"lma.txt", Plot->"lma.png"];
Print[S /. m, DPPlus /. m, DPMinus /. m];
Print[TouschekLifetimeLMA[]];
```

`MomentumAcceptance[位置のリスト (番号・名前・パターン), opts]`: 位置を省略すると `Range->{ls, le}`・`Step->k`・`Distance->ds` で選ぶ。`Cell->True` は、並びが n 回の繰り返しのリングの最初のセルだけ (周期的でない (誤差・挿入光源など) か閉軌道がセルで違うときは警告して全体)。
オプション: `Turns`、`NSynch` (2)、`DPMax`、`Tolerance` (1e-4)、`Coarse` (12)、`Sub` (3)、`Bucket` (True)、`Workers` (2)、`Aperture`、`CheckEvery` (**1**。0 で SAD の tapert のまま)、`Verify` (0)、`TurnCheck` (2。0 で無効)、`Output`、`Plot`。
戻りは `{S->…, Loc->…, DPPlus->…, DPMinus->…}`。結果はセッションに残り、`TouschekLifetimeLMA[]` が使う (Cell->True の結果は周期的に延ばす)。

### 4.4 Touschek への受け渡し

既存の Touschek (`touschek.py`: `Emittance[]` を INTRA ON で走らせると作る `TouschekTable`) は、一様な運動量受容 a (0.004〜0.242、0.002 刻み) の関数として、要素 l ごとの損失率 R_l(a)
(`table[3]`、要素の数 × 120) を持つ。`touschek_lifetime_local` は、要素ごとの受容 δ+(l)・|δ−(l)| (要素の中央の s での線形補間。周期なら折り返す) から

  1/τ = Σ_l ½ [ R_l(δ+_l) + R_l(|δ−_l|) ]

を求める (対称な受容 a の損失率は +Δ と −Δ を受け取る両方の粒子の損失を含むので、片側ずつ半分)。R_l(a) は dp の表の上の log-log 補間 (表の外は a⁻² で外挿し、外挿した要素の数を `extrapolated` で返す)。
**一様な受容 δ+ = |δ−| = a (表の点) を渡すと、`TouschekTable` の損失率と一致する** (試験)。既存の `TouschekLifetime[INF, INF, nz]` は一様な受容の寿命。
既存コードは変えていない (表の構造から受け取るだけ)。

**Piwinski 型の式 (v2.3.0 T2、v2 の既定)**: K64 の表の損失率 (`touckf`) は、平坦ビームの関数に ptrans² = tr(P) を入れた近似で、結合が小さくないリングで寿命が ±20% ずれる (B162)。
`TouschekLifetimeLMA[]` と `TouschekLifetime[INF, INF, nz]` の v2 の既定は、2 粒子の Møller 散乱の積分 (`touschek.touschek_piwinski`、`lma.touschek_lifetime_local_piwinski`) を使う。
式: 静止系 (z* = γz、運動量 p* [mc]) で、同じ位置の 2 粒子の相対運動量 q は N(0, 2P) (P = 位置を固定したときの運動量の共分散。η'・結合・z-δ の相関を含む)。質量中心系の非相対論的な Møller 断面積
dσ/dΩ = (4 r0²/k⁴) (4/sin⁴θ − 3/sin²θ) で、散乱で δ が a+ を超えて増える (または −a− を超えて減る) 率は、ps* の広がりが x-y 面の広がりより十分小さいとき、P の x-y ブロックの固有値 λ1 ≥ λ2 (mc²)、α_i = 1/(4λ_i) として
1/τ_inst = r0² c N / (γ V* √(λ1 λ2)) Σ± J(2β a±)、 J(c) = ∫_c^∞ dk k⁻² W(k/c) e^{−k²α1} e^{−z} I0(z) (z = k²(α2 − α1)/2)、 W(κ) = 2π (κ² − 1 − ln κ)、 V* = (4π)^{3/2} √det(xx*)
で、リングの寿命は 1/τ = Σ_slice (al/C) /τ_inst (スライスは K64 の tintrb の分割)。平坦ビームの極限 (λ2 → 0) は Bruck・Piwinski の関数 D(ε) (記憶に基づく式) に一致する。
`Emittance[]` (INTRA ON) は各スライスの静止系のビームを `ffs.touschek['rest']` に記録する。Python: `touschek_piwinski(ffs.touschek, a=0.006)`、`touschek_piwinski(ffs.touschek, ap=δ+ の配列, am=|δ−| の配列)`。
台本: `TouschekLifetime[INF, INF, nz, Formula->"K64"|"Piwinski"]`、`TouschekLifetimeLMA[Formula->...]` (`PYSAD_COMPAT=k64` と `Formula->"K64"` は K64 の表。nx か ny が有限の `TouschekLifetime` と EMIT の報告は K64 の表だけ)。
検証 (`example/v2_3/T2_*`、[`docs/history/v23_T2.md`](history/v23_T2.md)): 3 次元の一般の積分と実験室系の対のモンテカルロ (厳密な相対論の運動学) に 0.3% 以内、elegant (全要素を 8 分割した標本化) と −1.3〜+0.1%、K64 の表は −14〜+21%。
elegant の `touschekLifetime` は要素の出口だけを標本化するので、粗いリングでは寿命が +1〜+5% 長く出る。

**バンチ内の δ の分布 (v2.4.0 A3、オプション)**: 台本 `TouschekLifetimeLMA[BunchDelta->True]`・`lma.touschek_lifetime_local_piwinski(bunch_delta=True)`・`touschek_bunch.touschek_piwinski_bunch`。散乱後の粒子の δ は δ_c + Δ (δ_c = 対の重心の δ ~ N(0, σδ²/2)、Δ は相対運動量の散乱だけで決まる) なので、受容を a ∓ δ_c にして平均する。受容が σδ の 8〜14 倍のとき損失率は ASTRID2 で +1.3〜+2.7%、KEK-LS で +1.8〜+4.3% 増える (寿命が短くなる)。実験室系の 2 粒子の対のモンテカルロ (厳密な相対論の運動学) で δ_c の効果が 2% 以内で再現される。**既定は入れない**: 相対運動量の縦成分 q_s の効果 (MC で通常の kick の判定より損失率が ASTRID2 で 2%・KEK-LS で 0.5% 小さい、つまり絶対値の判定との差がその分大きい) は `BunchDelta->True` だけでは入らない。**v2.4.0 A10: `BunchQs->True` を足すと q_s の効果も入る** (3 次元の積分 `touschek_bunch.qs_ratio`。絶対値 / kick の損失率の比が MC と ASTRID2 +4.81・+2.22% (MC +4.84・+2.22%)、KEK-LS +4.5・+1.7% (MC +4.5・+1.8%) で 0.3% 以内。12 スライスで求めて補間するので全周で 70 秒。ASTRID2 全周で a/σδ = 9・14・23 の寿命が δ_c = 0 に対し δ_c だけ −2.5・−1.2・−0.5%、q_s も入れて −4.1〜−4.3・−1.9・−0.7% (補間の不確かさ ±0.3%))。elegant・K64・Piwinski の式は kick の基準を δ_old = 0 とし、elegant の `touschekLifetime` には δ の分布を入れる機能が無い。a ≲ 3σδ では使えない (警告)。詳細・限界は `pysad/touschek_bunch.py` の冒頭と `docs/history/v24_A.md` §A3、検証は `example/v2_4/A3_mc_bunch.py`。

### 4.5 信頼して使うための条件と警告 (`pysad/lma.py`・`pysad/lmabucket.py`、v2.3.0 L2)

LMA の信頼性の基準: (1) 解析・追跡・elegant の比較が条件つきで揃っている、(2) 既知の不一致に原因がある、(3) 既定で安全な値を返す、(4) 使えない条件で警告を出す。条件は ASTRID2 (0.58 GeV、45.7 m、δ_rf = 1.618%、νs = 0.00152) と
大きなリング A (`kek-ls_v03_12.lat`、3 GeV、603 m、δ_rf = 2.994%、νs = 0.00408) の ±10 mm (例の値)、RF・放射つき。詳細・数値は [`docs/history/v23_L2.md`](history/v23_L2.md)、例は `example/v2_3/L2_*`。

**(a) check_every の既定は 1 (要素ごとの判定)**。SAD の tapert (15 要素ごと = 0) は、アパーチャが限界の位置で δ を大きく見積もる (ASTRID2 の 25 位置で δ+ を最大 1.46 倍、δ− を 1.52 倍。位置 17: 1.048% 対 0.829%)。
要素ごとの判定は elegant (要素ごとの MAXAMP) と δ+ 3.5e-5、δ− 7.0e-5 以内。5 要素ごとは 1 との差 2.6e-5 以下。**0 を明示するとアパーチャで決まる位置があれば警告する** (gs の `TrackParticles` と同じ判定で比べるときに使う)。
時間は ASTRID2 (199 要素) で +7%、大きなリング A (3563 要素、全位置が RF バケットで 0・1・5 が同じ値) で ×2.4 (156 s → 383 s、11 位置)。時間を節約するときは 5。

**(b) RF バケットの限界とアパーチャの限界** (`lmabucket.py`)。LMA の値 = min(RF バケットの限界、実効のアパーチャの限界):
- RF バケットの限界は、解析式 `rf_bucket_analytic` (δ_rf² = (2A/(kCα))[2 cos φs − (π − 2φs) sin φs]、`Emittance[]` の bucket_height と 1e-6 で一致)、縦方向の 1 次元の写像 `bucket_limits_1d` (解析式の RF の力 + 追跡で測った非線形の経路長
  Δz = −C (a1 δ + a2 δ²) と放射の δ 依存)、6 次元の追跡 `bucket_limits_6d` (アパーチャ無限大) の 3 通りで求まる。1 次元の写像と 6 次元の追跡は 1e-4 以内で一致 (ASTRID2: 1.643% と 1.635〜1.6375%、大きなリング A: 3.012% と 3.021%)。
  理想化した解析式 (減衰なし・α 一定) との差 (ASTRID2 で δ+ が +1.5%、δ− が +0.4%) は、放射の減衰・α2 (a2/a1 = −0.76) と離散写像で説明がつき、**LMA の δ+ が δ_rf を 1% 超える (Q1・T1) のは有限の周回数と、この効果**。
  出発点のベータトロン振幅 (δ = 0 の閉軌道から出発すると経路長が変わる) で位置により 1.4% 変わる (分散の大きい位置の 6 次元の限界 1.660%、δ の閉軌道から出発すると 1.6375%)。
- アパーチャの限界は、RF・放射なしの 4 次元の追跡 (δ を固定)。シンクロトロン振動で半周期のうちに δ の符号が反転するので、実効の限界は min(δ+_aper, |δ−_aper|) になる (追跡で確認: δ0 = −0.86% の粒子が半周期後に δ = +0.86% で x = +10 mm に当たって失われた)。
- `limit_types(setup, result)` が位置ごとに rf・aper・both を返す。**ASTRID2 (±10 mm) の 25 位置 (δ+・δ−) はすべてアパーチャで決まる** (LMA − min の最大: δ+ 2.6e-5、δ− 2.6e-4 (分散の大きい 12 位置の δ− で LMA が実効の限界を 3% 超える。原因は切り分けていない))。
  大きなリング A (±10 mm) の 11 位置は RF バケット 7 点 + バケットとアパーチャ (3.03%) が 2 tol 以内 4 点で、LMA = min(rf, aper) が一致 (差 0)。
  **以前の「位置 2 は RF」は check_every = 0 の値 (δ+ = 1.635%、理由 'rf') の取り違え**。`reason` は最初に失われた粒子を落とした直接の原因 ('rf' は |z| の判定が先に効いたとき) で、限界の種類とは別。バケットの限界でも、バケットから出てすべった粒子はアパーチャに当たって 'aper' と記録されうる。

**(c) 周回数**。バケットで決まる限界は、シンクロトロン振動 1.5 周期以上で 1% 以内、2 周期 (既定) で二分の幅以内、3 周期で収束する。1 周期以下は限界が 40〜60% 大きく出る (ASTRID2 の位置 2 の 6 次元: 0.25 周期 2.52%、1.0 周期 1.98%、1.5 周期 1.64%、2.0 周期 1.634%、3.0 周期 1.634%)。
大きなリング A (11 位置) は周回数 491・1000・1500・2000 で δ± の変化が 5e-5 以内。elegant (1500 周) と同じ判定 (`bucket=False`、アパーチャまで追跡) で同じ 1500 周にすると、δ+ の差は +5.2e-5、δ− の差は 1.3e-5〜7.8e-5 (elegant の最終刻み 7.8e-5 以内)。
`turncheck` (既定 2): 限界の tol 手前の粒子を 2 倍の周回数で追跡し、失われれば `T±` (tol より大きく下がる) と警告する。

**(d) 狭い損失の窓**。`verify=n`: 限界の手前 [0, δ] を n 点で再追跡し、窓があれば限界をその手前に下げる (`nwin_*`・`win_*`・旗 `W±`・警告)。ASTRID2 (25 位置、n = 100、刻み 1.6e-4)・大きなリング A (11 位置、n = 150、刻み 2e-4) で窓は **0 個**
(限界は verify なしと同じ)。既定は 0 (ASTRID2 で ×5.9、大きなリング A で ×2.8)。刻みより狭い窓・限界の外側に残る島は見ない。

**(e) 1 セルから全周**。RF 空洞は 1 か所なので周期性は崩れ、閉軌道の δ がセルごとに U0/E0 ずつ変わる (ASTRID2 1.06e-5、大きなリング A 9.9e-5)。全セル (ASTRID2 6、大きなリング A 20) の同じ位置の LMA の差は、ASTRID2 で 2.6e-5 (tol 1e-4)・
4.9e-6 (tol 5e-6)、大きなリング A で 4.9e-5 (二分の刻み 1 つ) で、**空洞からの距離への系統的な依存はない**。周期性の許容 tol = 1e-4 の根拠は、二分の幅と U0/E0 (どちらも 1e-4 以下)。U0/E0 が大きいリング (24 セルの FODO で 5e-4) は測っていない。
`periodicity_report(setup)`・`cell_locations(setup)` が、要素の並びが周期的でない (誤差・挿入光源・セルの長さの違い)・閉軌道がセルで違う (許容 1e-7 + 2 (U0/E0) max(|η|, |η'|)) ときに警告して全周の位置を返す。台本の `Cell->True` も同じ警告。

**(f) FLUC (量子励起)**。`TrackSetup.copy(fluc=True)`・`seed` で走る。結果は 1 つの標本なので、複数の種の平均・分散を使う (警告が出る)。ASTRID2 (2 周期、4 位置、6 種): FLUC あり − なし = −1.8e-5〜+1.3e-5 (平均)、種の分散 1.4e-5〜3.2e-5。
elegant の ISR あり (6 種) は ISR なしと 9e-5 以内。13120 周 (20 周期、σδ の成長は定常の 45%) でも 5 種の平均は FLUC なしと 2.1e-5 以内。
振幅 (0.8〜1.6%) が σδ (4.3e-4) の 20〜40 倍なので、量子励起による振幅の変化 ≈ σδ²/δ が 2e-6〜2e-5 と二分の幅より小さいため。大きなリング A (σδ = 6.9e-4、位置 16・89、6 種、tol 1e-4): pysad の FLUC あり (491 周と 1500 周で同じ) は δ+ = 3.027%・3.029% (σ 1.1〜1.3e-4)、δ− = −3.017%・−3.005% (σ 4e-5〜5e-5)、FLUC なしとの差 +0.8e-5〜+8.1e-5 (σ 以内)。elegant の ISR あり 1500 周 (6 種) は δ+ = 3.008%・3.010% (σ 1.0e-4・1.2e-4)、δ− = −3.022%・−3.009% (σ 5.9e-5・7.7e-5)、pysad との差は δ+ で +1.9e-4 (σ の 1.5〜2 倍)、δ− で +0.5e-4〜+0.4e-4 (T1 の ISR なしの差 +5.7e-5〜+1.05e-4 に近い)。elegant の ISR あり 491 周は位置 89 で δ± = +3.19%・−3.25% (σ 7e-4・1.4e-3) と収束しない (バケットから出た粒子がアパーチャに届かないため。T1 と同じ理由。pysad はバケットの判定があるので 491 周でよい)。

**(g) 誤差つき・挿入光源**。誤差 (`fmapres.apply_alignment_errors`、四極 50 µm・0.2 mrad、6 極 50 µm、偏向 50 µm・0.1 mrad の rms) つきの ASTRID2 (3 種) は全周の位置で走る (25 位置、1 つ約 5.5 分)。`find_periodicity` が `nper = 1` を返し、
`periodicity_report`・`cell_locations`・`Cell->True` が警告する。誤差で閉軌道が歪む (cod1 の x・y 0.3〜0.9 mm) ため、δ± は誤差なしより小さくなる (δ+ の平均で −0.6〜−1.5e-3、閉軌道の歪みの 0.5 mm / アパーチャ 10 mm × 0.8% ≈ 4e-4 と同じ桁)。
1 セルの結果を誤差つきの全周に使うと、同じ位置の別のセルとの差が 0.8〜2.6e-4 になる。MPW 版 ASTRID2 (`ASTRID2_v01_02_id.sad` の `RINGMPW`) は周期性が崩れる (`nper = 1`) ので警告し、全周 32 位置で走る (δ+ = 0.836〜1.590%、δ− = −0.856〜−1.572%、δ_rf = 1.571%。LMA = min(rf, aper) が 2.8e-4 以内)。
結合の強いリング・補正つき・大きなリング A の誤差つきは未確認。

**(h) 警告** (`LMAResult.warnings`・`cond['warnings']`・`warnings.warn`・表の `# WARNING:`): LOSSAMPL ≥ 0.1 m (アパーチャ未設定)、check_every = 0 でアパーチャで決まる位置がある、FLUC ON、RFSW OFF、周回数がシンクロトロン振動 1 周期未満、
走査の上限 (dpmax) まで残った限界、δ = 0 で失われた、周回数の依存 (turncheck)、verify で窓、閉軌道の残差 > 1e-9。位置ごとの旗は `LMAResult.flags(i)` と表の最後の列 (`W±`・`T±`・`S±`・`0±`)。

## 5. 検証

[`example/v2_3/`](../example/v2_3/) の例と [`pysad/tests/test_v23_Q1_fmap_lma.py`](../pysad/tests/test_v23_Q1_fmap_lma.py) (21 件)。使った gs は K64 の 2026-10-19 ビルド (WSL、timeout 1800 s、終了後に `ps` で sad1.exe の残りがないことを確かめた)。

| 項目 | 内容 | 結果 |
| --- | --- | --- |
| 損失判定 (gs) | ASTRID2 (0.58 GeV、RFSW・RAD ON、FLUC なし、LOSSAMPL = ±10 mm、LOSSDZ = バケットの半分 1.5407 m、1312 周)。位置 2 (長直線部。δ+ はアパーチャ、δ− はバケットとアパーチャが近い) と 17 (分散の大きい位置、アパーチャが限界) で、閉軌道に δ = ±(0〜2%) を足した 30 粒子を、1 回の `TrackParticles[{loc, beam}, loc, 1, 1312]` で追跡 | 生き残りの旗は、gs・pysad の台本・`TrackSetup.track` の 3 者で 30 / 30 一致 (落ちた 11 個、残った 19 個)。最後の x の最大の差 1.55e-3 m (限界ぎりぎりの粒子は、gs の `Emittance[]` の閉軌道・同期位相の 1e-7 の違いで軌道が変わる) |
| 二分法 (gs) | gs の台本で `TrackParticles` を 1 粒子ずつ呼ぶ二分法 (9 回、区間 [0, 0.02]) で、位置 2・17 の δ±の区間 (lo, hi) を求め、pysad の台本の同じ二分法と、`momentum_acceptance` (nsub = 3、tol = 2e-5) と比べる | 4 / 4 で gs と pysad の台本の区間が完全に同じ (例: 位置 17、δ+: 0.010469〜0.010508)。LMA の区間 (0.010479〜0.010488) はその内側 |
| 周波数マップの追跡 (gs) | ASTRID2、位置 2、NORFSW・NORAD、6 点 (1 mm から 18 mm)、gs で 1 周ずつ `TrackParticles` を回し、512 周の (x, px, y, py) を印刷。その座標列に NAFF をかけて D を出し、`frequency_map` と比べる | 座標の最大の差 2.8e-7 m (D = −2.9 の最も混沌とした点)、D は 4 桁一致 (最大の差 5.2e-7)、νx は 1e-13 の桁 |
| 周期性 | ASTRID2 (6 セル、セルの長さ 7.6174 m) と大きなリング A (20 セル、30.156 m) で、1 セルの位置のなかの 3 点を、別の 3 つのセル (ASTRID2: 2・3・5 番目、大きなリング A: 3・9・15 番目) の同じ位置でも計算 | 最大の差 |Δδ| = 2.6e-5 (ASTRID2)、4.9e-5 (大きなリング A)。二分の幅 tol = 1e-4 以下 |
| RF バケット | ASTRID2 (δ_rf = 1.6176%、νs = 0.00152)、大きなリング A (δ_rf = 2.994%、νs = 0.00408) | ASTRID2 (±10 mm): δ+ の最大 1.598% (δ_rf の 0.988 倍)、δ− の最小 −1.609% (0.995 倍)。大きなリング A (±10 mm): δ+ の最大 3.026% (1.011 倍)、δ− の最小 −3.021% (1.009 倍)。**δ+ ≲ δ_rf** (1.1% の超過は二分の幅・バケットの判定の余裕の範囲)。SAD の既定の判定 (check_every = 0)・アパーチャ 20 mm の別の計算では ASTRID2 の δ+ が δ_rf の 1.024 倍まで出た (バケット判定の余裕が (π + 2φs)/π で、バケットの端の粒子の滑りが 2 周期より遅いため) |
| L2 (§4.5) | check_every・バケットとアパーチャの区別・窓・周回数・1 セルから全周・FLUC・誤差・MPW | [`example/v2_3/L2_README.md`](../example/v2_3/L2_README.md)、[`docs/history/v23_L2.md`](history/v23_L2.md)。「RF バケット」の行の δ+ の超過の原因と、位置 2 の分類は §4.5 (b) |
| 線形の限界 | ASTRID2 ±10 mm の分散の大きい位置 (10 位置): δ+ = 0.829%、δ− = −0.856%。水平の線形の見積もり (§4 の `linear_momentum_limit`) 0.865% | δ+ / δ_lin = 0.96 (check_every = 1)、1.21 (check_every = 0、SAD の既定の判定)。δ+ ≤ min(δ_rf, δ_lin) × 1.00 (ε の最大 −0.012) |
| Touschek | INTRA ON、PBUNCH = 1e9 粒子/バンチ。ASTRID2: 局所の δ±(s) (1 セル、周期的に延ばす)、一様な δ = δ_rf、一様な δ = LMA の最小値 | 寿命 9502 s、19850 s、2918 s。大きなリング A: 26876 s、26045 s、26473 s (SIGE = 7.53e-4)。一様な受容の値は `TouschekLifetime[INF, INF, nz]` と同じ表。一様な a を渡すと表の損失率と 1e-10 で一致 (試験) |
| 台本の関数 | `FrequencyMap[]`・`MomentumAcceptance[]`・`TouschekLifetimeLMA[]` ([`example/v2_3/Q1_script_commands.sad`](../example/v2_3/Q1_script_commands.sad)) | Python の API と同じ値 (試験 `test_script_commands`)。ASTRID2 の最初のセルの 6 位置: δ+ = 1.598%・0.829%、寿命 9926 s |
| elegant との比較 (T1) | ASTRID2 (全周 199 位置の elegant `momentum_aperture`、pysad は 1 セル 25 位置、check_every = 1、±10 mm、1312 周)、kek-ls_v03_12.lat (elegant 14 位置・1500 周、pysad 52 位置・491 周)。RF・放射つき (ISR なし)、同じアパーチャ・同じ周回数 (ASTRID2)。Touschek は elegant の `touschekLifetime` (Piwinski の式) に pysad の IBS 込みのビームを渡す | **LMA: δ± の相対差 0.5% 以内** (ASTRID2: アパーチャが限界の δ = 0.83% も RF バケットで決まる 1.60% も、最大の差 7e-5。kek-ls: 全位置が RF バケット、0.35% 以内)。**Touschek (一様な δ、表の点): elegant / pysad = −4.3〜+2.4% (ASTRID2)、−3.4〜+1.2% (kek-ls)**。局所の LMA の寿命: ASTRID2 pysad 9502 s・elegant 9069 s (−4.6%)、kek-ls 26876 s・27117 s (+0.9%)。差は LMA でなく寿命の式 (pysad は SAD の `touschek_f` の当てはめと等方な運動量分散の近似) による。elegant の `momentum_aperture` はバケットから出ただけでは損失にしないので、周回数が足りないと δ+ が見つからない (kek-ls は 1500 周)。[`example/v2_3/T1_README.md`](../example/v2_3/T1_README.md)、[`docs/history/v23_T1.md`](history/v23_T1.md)。**T2: 独立な式 (Piwinski 型、3 次元積分・厳密な運動学のモンテカルロで確認) では elegant (標本化を細かくしたとき) と −1.3〜+0.1%、K64 の式が不正確 (B162)。局所の LMA の寿命 (Piwinski 型): ASTRID2 9327 s、大きなリング A 26862 s** |

計算時間 (2 並列 + 周波数マップ 1 本を同時に走らせた状態。CPU を使い切らない): ASTRID2 の LMA 25 位置 364 s (追跡した粒子 1075、1312 周)、大きなリング A の LMA 43 位置 580 s (1849 粒子、491 周、
check_every = 5)、周波数マップ 1271 点 × 1024 周は 305〜315 s (1 プロセス)。追跡の速さは 1 粒子・1 周あたり ASTRID2 (200 要素) で 0.25 ms (放射なし)〜0.5 ms (放射あり)、大きなリング A (3563 要素) で 0.6 ms (放射あり)。

## 6. 使った条件 (例)

| 項目 | ASTRID2 ([`example/astrid/lattice/ASTRID2_v01_01.sad`](../example/astrid/lattice/ASTRID2_v01_01.sad)) | 大きなリング A ([`pysad/tests/kek-ls_v03_12.lat`](../pysad/tests/kek-ls_v03_12.lat)) |
| --- | --- | --- |
| エネルギー、周長、要素数、セル | 0.58 GeV、45.704 m、199 要素 (+1)、6 セル (セルの長さ 7.617 m) | 3 GeV、603.12 m、3563 要素 (+1)、20 セル (30.156 m) |
| RF | 50 kV、104.95 MHz、U0 = 6.17 keV、νs = 0.00152、δ_rf = 1.618%、同期位相 φs = 0.124 rad | 1.4 MV、HARM 992、U0 = 296 keV、νs = 0.00408、δ_rf = 2.994%、φs = 0.213 rad |
| 誤差 | なし | なし |
| LMA | RFSW ON、RAD ON、FLUC なし、RADCOD の閉軌道 (残差 2e-18)、LOSSAMPL = ±10 mm (例の値で実機の値ではない)、check_every = 1、周回数 1312 (= 2/νs)、バケット判定 ON (|z| < 1.5407 m)、粗い走査 12 点 (dpmax = 1.25 δ_rf)、nsub = 3、tol = 1e-4、位置は最初のセルの要素の境界 25 点、粒子 1075 個、2 並列 | 同じ。LOSSAMPL = ±10 mm、check_every = 5、周回数 491、バケット判定 |z| < 0.3452 m、位置は最初のセルを 0.5 m 以上の間隔で 43 点、粒子 1849 個、残差 4e-16 |
| 周波数マップ | RFSW OFF、RAD OFF、δ = 0、位置 2 (長直線部の入口、βx = 3.89 m、βy = 1.62 m)、x = 0.5〜18 mm・y = 0.5〜14 mm の 41 × 31 (1271 点)、1024 周 (前半・後半 512 周)、NAFF の窓 Hann、規格化 = 標本の楕円、LOSSAMPL = 0.1 m、px = py = z = 0 から出発 | 位置 3 (長直線部の入口、βx = 4.64 m、βy = 6.23 m、ηx = 4.5 cm)、x = 0.2〜6 mm・y = 0.2〜8 mm の 41 × 31 (1271 点)、同じ |
| Touschek | INTRA ON、PBUNCH = 1e9 粒子/バンチ (例の値)、SIGE = 4.39e-4、表の dp = 0.004〜0.242 | 同じ。SIGE = 7.53e-4 |
| 計算の互換モード | v2 (`compat.set_mode('v2')`) | v2 |
| gs との比較 | K64 の 2026-10-19 ビルド、NPARA = 1 | — (gs の比較は ASTRID2 だけ) |

## 7. 未確認・限界

- 周波数マップは、誤差なしのほかに誤差つき (全周の追跡、§9.4・§10.1)・MPW つき (§9.5・§10.2)・k64 (§9.6)・RFSW ON・RAD ON・FLUC ON (§10.3) でも走ることを確かめた (S2 で複数の種・強さ・種類を統計、KEK-LS も)。**LMA は、誤差つき ASTRID2 (3 種) と MPW 版 ASTRID2 でも全周で走ることを確かめた (§4.5 (g))。結合の強いリング・補正つき・大きなリング A の誤差つきは未確認**。
  誤差つき・結合のあるリングの D は、誤差なしと別の量 (D > −5 が 7〜12 倍) で、結合が強いと νx・νy の対応が崩れる (警告 `coupling`)。FLUC つきは D が乱数の種に依る。
  FLUC ON の LMA は 1 つの標本 (複数の種の平均・分散で使う。警告)。2 周期では効果は二分の幅以下 (§4.5 (f))。
- 周波数マップの D は、粒子が `nturn` の途中で失われた点では nan。nturn が短いと D の下限 (§3.1) が上がる。共鳴線の上の粒子 (島の中に捕らえられた粒子) は、νx = νy のように
  2 面のチューンが同じ値になり、D は小さい (規則的な運動) ことがある (ASTRID2 の例の νx − νy = 3 の線)。
- LMA の損失の判定は、既定 (`check_every = 1`) では要素ごと。SAD の tapert のまま (`check_every = 0`) だと、アパーチャが限界のとき δ の限界を最大 1.5 倍大きく見積もる (§2、§4.5 (a)。警告が出る)。
  `check_every ≥ 1` は SAD の `TrackParticles` と判定が違う (gs と判定を比べるときは 0 にする)。`check_every = 1` と 5 の差は小さい (ASTRID2 で 2.6e-5 以下)。
- **Touschek 寿命の式 (T1・T2)**: K64 の表の式 (`touckf`) は、同じビーム・同じ一様な δ で独立な式 (Piwinski 型) と −14〜+21% 違う (結合 1% で +0〜+6%、丸いビームで ±20%。原因は x-y 面の運動量の分布を平坦ビームの 1 次元の関数と tr(P) で近似していること。B162)。elegant の式は独立な式と −1.3〜+0.1% (elegant は要素の出口だけの標本化で、粗いリングで +1〜+5% 長く出る)。v2 の既定は Piwinski 型。限界: 受容 a が σδ の 10 倍以下では、バンチ内の δ の分布 (a − δ_old) を入れていない (4 つの式とも共通) ことが 5〜20% 効きうる (概算。[`docs/history/v23_T2.md`](history/v23_T2.md) §7)。K1≠0 の BEND の K1 は B' 一様ではなく (1 + hx)∂b/∂x 一様の定義 (B160・docs/11 A59。誤りではない) で、elegant の CSBEND には K2 = −hK1 を足して合わせた。ASTRID2 の LMA には影響が無かった。
- バケット高さの付近の限界は、周回数に依存する。シンクロトロン振動 1.5 周期以上で 1% 以内、2 周期 (既定) で二分の幅以内、1 周期以下は 40〜60% 大きく出る (§4.5 (c)。`turncheck` と警告)。RF バケットの高さ δ_rf は `Emittance[]` の値 (解析式と 1e-6 で一致)。
  限界は δ_rf より 1〜2% 大きい (放射の減衰・α2・ベータトロン振幅。1 次元の写像と 6 次元の追跡が 1e-4 以内で一致、§4.5 (b))。RF の高調波・バケットの外の運動量の分布は見ていない。
- 1 セルの結果を全周に延ばす (`touschek_lifetime_local` の period) のは、RF の位置の違い (空洞が 1 か所) と閉軌道の δ の違いを無視する。別のセルの同じ位置との差は、ASTRID2 (U0/E0 = 1e-5) で 5e-6〜2.6e-5、大きなリング A (U0/E0 = 1e-4) で 4.9e-5 以下 (§4.5 (e))。U0/E0 が tol (1e-4) より大きいリングは未確認。周期でない (誤差・挿入光源) ときは警告する。
- 線形の見積もり `lma.linear_momentum_limit` は水平の分散・ベータ関数だけを使う (垂直の分散・結合・非線形を入れない)。
- 損失の判定の位置は要素の境界。要素の中の損失は見ない。
- 互換モードは v2 で確かめた。**k64 の周波数マップは S1 で確かめた** (ASTRID2・KEK-LS で v2 と同じ、チューンの差 5e-11、§9.6)。**LMA の k64 は v2.4.0 A4 で確かめた** (ASTRID2 の 8 位置で v2 と δ± が 6e-8 以内。check_every = 1。k64 の最初の tapert の扱い (B33) の違いは check_every = 0 のときだけで、見ていない)。
- 周波数マップの D の絶対値の閾値 (D > −5 など) は周回数・窓で動く (周回数 2 倍で約 −1、§9.2)。S2 で 1024 周への換算 `D_norm` を足した (§3.6、§10.5。512〜4096 周)。RDT・ADTS の予想は 1 次の共鳴だけで、ASTRID2 の D の大きい領域 (νx − νy = 3 など 2 次の共鳴) は予想できない (§9.1)。
- 並列 (`workers` ≥ 2) は子プロセス (Windows の spawn) を使うので、Python のスクリプトは `if __name__ == '__main__':` の中から呼ぶ。台本の関数 `MomentumAcceptance[]` の `Workers` も同じ。
- numba のカーネル (`track6n.py` など) は編集していない。最初の呼び出しはコンパイルのキャッシュで数秒かかる。

## 8. 実装の構成

| ファイル | 内容 |
| --- | --- |
| [`pysad/trksetup.py`](../pysad/trksetup.py) | `TrackSetup` (追跡の設定、閉軌道、`track`、バケット)、`find_periodicity`、`select_locations`、`lower_priority` |
| [`pysad/freqmap.py`](../pysad/freqmap.py) | `frequency_map`、`FMapResult` (表・図・読み込み)、`resonance_lines`・`resonance_check`・`line_distance` |
| [`pysad/fmapcheck.py`](../pysad/fmapcheck.py) | (S2) 誤差の統計・フィードダウンの見積もり (`errors_summary`・`feeddown_estimate`・`center_sextupoles`)、`map_stats`、周回数の換算 (`threshold_for_nturn`)、警告 (`check_setup`・`check_result`) |
| [`pysad/fmapres.py`](../pysad/fmapres.py) | (S1) `ResonancePredictor` (RDT・ADTS から共鳴の領域を予想)、`coverage_stats`、`lines_distance`・`other_mn`、`apply_alignment_errors`、`plot_with_prediction` |
| [`pysad/lmabucket.py`](../pysad/lmabucket.py) | L2: `rf_bucket_analytic`・`longitudinal_params`・`bucket_limits_1d`・`bucket_limits_6d`・`aperture_limits_4d`・`limit_types` (RF バケットの限界とアパーチャの限界の区別、§4.5) |
| [`pysad/lma.py`](../pysad/lma.py) | `momentum_acceptance`、`LMAResult`、`touschek_lifetime_local`、`acceptance_at`、`linear_momentum_limit` |
| [`pysad/q1cmd.py`](../pysad/q1cmd.py) | SAD の台本の関数 `FrequencyMap[]`・`MomentumAcceptance[]`・`TouschekLifetimeLMA[]` (`ffscmd.py` の入口 2 行) |
| [`pysad/tests/test_v23_Q1_fmap_lma.py`](../pysad/tests/test_v23_Q1_fmap_lma.py)・[`test_v23_S1_fmapres.py`](../pysad/tests/test_v23_S1_fmapres.py) | 試験 |
| `example/v2_3/Q1_*`・`S1_*`・`FMAP_summary.md` | 例と出力 (索引は [`example/v2_3/Q1_README.md`](../example/v2_3/Q1_README.md)・`S1_README.md`)。台本ごとのまとめは [`example/v2_3/FMAP_summary.md`](../example/v2_3/FMAP_summary.md)。記録は [`docs/history/v23_Q1.md`](history/v23_Q1.md)・`v23_S1.md` |

## 9. S1 の検証 (v2.3.0、2026-11-20)

`example/v2_3/S1_*`、試験 [`pysad/tests/test_v23_S1_fmapres.py`](../pysad/tests/test_v23_S1_fmapres.py)、台本ごとのまとめと全数値は [`example/v2_3/FMAP_summary.md`](../example/v2_3/FMAP_summary.md)、記録は [`docs/history/v23_S1.md`](history/v23_S1.md)。
条件は §6 と同じ (v2、RFSW・RAD OFF、LOSSAMPL = 0.1 m、41 × 31、1024 周、Hann 窓)。FODO は 24 セル、3 GeV、105.6 m、QF の K1 = 0.697 (νx = 6.6578)、SF・SD の K2 = 4, −6 (玩具の強さ)。

### 9.1 共鳴線と RDT・ADTS の定量的な接続 (a)

| 例 | D > −5 の点 | 1 次の RDT の領域に入る割合 cov (面積 area、対照 ± σ、p) | D > −5 の点の最も近い線が 1 次の調和 / 2 次以上の線 |
|---|---|---|---|
| FODO (3νx の共鳴が格子の中) | 44 | 1.00 (0.26、0.27 ± 0.34、p = 0.010) | 44 / 0 |
| ASTRID2 | 48 | 領域が格子に入らない (面積 0) | 0 / 23 |
| KEK-LS (追跡のチューンを使う。線形の ADTS では対照より悪い) | 124 | 0.73 (0.44、0.28 ± 0.13、p = 0.003) | 11 / 88 |

- FODO では共鳴の中心 x(Δ0 = 0) の予想が D の尾根の x と平均 0.52 ± 0.28 mm (線形の ADTS、追跡を使わない)、0.02 ± 0.19 mm (追跡のチューン) で一致する。予想の「捕らえられる」領域は D の大きい範囲の約 2 倍の太さ (島の中は規則的、D の大きい層は島の縁)。
- 振り子の式 (位相つき) は、同じ生成関数の 1 周の写像の追跡と、捕らえられる J0 の範囲の端が 0.15 J_res 以内で一致する (試験)。
- **限界**: ADTS が大きい (KEK-LS: xx = 5.1 万) と、線形の ADTS の ν(J) は大振幅で外れ (全点の中央値 0.068)、領域が格子の大部分を覆って判別力が小さい。ASTRID2 の D の大きい点は νx − νy = 3 (2 次の共鳴)・5νy などで、1 次の RDT の共鳴は格子のチューンの動く範囲に入らない。
  原因は、6 極の強さを λ 倍した周波数マップ (D > −5 の点: λ = 0、0.5 で 0 点、1 で 10 点、1.5 で 37 点) の通り、2 次以上の共鳴。**6 極の 2 次の共鳴駆動項の解析式が `resdrive.py` に無い** (R4 の担当範囲)。

### 9.2 周回数・窓・格子 (f)

周回数を 2 倍にすると D が約 1 下がる (中央値 −6.79 → −7.86 → −8.97 → −10.17、512〜4096 周)。D > −5 の割合は 0.065 → 0.030 → 0.021 → 0.009。順位相関 0.80〜0.88 (規則的/カオス的の順序は残る)。窓なしは D の下限が上がり (−5.7)、D > −5 の割合が 5 倍になるので使わない。
Hann p = 1 (既定) と p = 2、3 の D の差は小さい (中央値 0.1〜0.3)。格子 11 × 8 から 81 × 61 で D の分布と D > −5 の割合 (0.023〜0.038) はほぼ変わらない。**D > −5 は周回数 1024 での目安で、別の周回数とは比べない。**

### 9.3 δ を含む周波数マップ (b)

x–δ・y–δ の格子は `frequency_map` で動く (RFSW OFF は δ の閉軌道を 4 次元の Newton 法で、x・y は δ の軌道からのずれ)。新しい実装は不要。ASTRID2 x–δ・y–δ (δ = ±2%、31 点) と KEK-LS x–δ で確認。
最小振幅の点の追跡のチューンと、δ ごとの線形のチューン + ADTS (その点の J) の差は 0.0007〜0.0019 以内 (ASTRID2)。D > −4 の点の 70〜75% が、予想した作動点が次数 4 以下の線の 0.003 以内の領域 (面積 7〜10%) に入る (ASTRID2、対照 0.06〜0.11)。KEK-LS は ADTS が大きく予想が外れる (0.013 の中央値)。δ ≠ 0 では偏向磁石の 3 次の項を入れていない。

### 9.4 誤差つき (d)

要素ごとの乱数の位置のずれ・ロール (四極・6 極・偏向磁石 DX・DY 50 μm、四極のロール 0.2 mrad、偏向磁石のロール 0.1 mrad、種 1) つきの ASTRID2 で、周期性を使わない全周の追跡の周波数マップが走る。
閉軌道は 4 次元の Newton 法 (残差 3e-16)。位置 2 の閉軌道は (−0.94, 0.07, 0.06, −0.16) [mm, mrad]。D > −5 の割合が 0.038 → 0.50、失われた点 0 → 7、最小振幅の D の下限が −10.7 → −7.9 (結合)。誤差 ×3 (coarse) は閉軌道 x = −2.8 mm、失われた点 79 / 336、D > −5 が 0.62。誤差のレベルは例。

### 9.5 挿入光源つき (e)

`ASTRID2_v01_02_id.sad` の `RINGMPW` (223 要素) の周波数マップが走る。νy = 2.1402 → 2.2492 (4νy = 9 の 0.0008 以内)、1271 点のうち 950 点が失われる (y の大きい列で 83〜100%)。MPW の内部 (K0 ≠ 0 の偏向磁石) の RDT・ADTS は `resdrive` の対象外。

### 9.6 k64 モード (c)

`PYSAD_COMPAT=k64` で ASTRID2 (1271 点)・KEK-LS (336 点) の周波数マップが走り、v2 と同じ (チューンの差 5e-11、失われた点の食い違い 0)。新機能なので互換モードの分岐はない。

### 9.7 gs との比較 (g)

KEK-LS (3563 要素、NORFSW・NORAD、512 周、位置 3、8 点): 規則的な 5 点 (D < −3.5) の座標が gs と 1.2e-11 m 以内、D の差 2.9e-7。カオス的な点 (D = −2.3、限界の点) の食い違い (座標 9.8e-3 m、落ちた周 369 対 0) は、
pysad の初期値を 1e-13 m ずらしただけの食い違い (座標 1e-2 m、落ちた周 114〜187) と同じ大きさで、丸め誤差の増幅。gs は 2 s、終了後に sad1.exe の残りなし。

### 9.8 S1 の未確認

誤差つきは ASTRID2 の 1 つの種だけ (KEK-LS・補正つきは未確認)。RFSW ON・FLUC つきの周波数マップ、PF リング。2 次以上の共鳴の強さの予想。δ ≠ 0 の偏向磁石の 3 次の項。MPW の RDT。

## 10. S2 の検証: 周波数マップを信頼して使う (v2.3.0、2026-11-21)

`example/v2_3/S2_*` (索引 `S2_README.md`)、記録 [`docs/history/v23_S2.md`](history/v23_S2.md)、試験 [`pysad/tests/test_v23_S2_fmapcheck.py`](../pysad/tests/test_v23_S2_fmapcheck.py)。条件: v2、RFSW・RAD OFF、LOSSAMPL = 0.1 m、Hann 窓、小さな格子 21 × 16 = 336 点、1024 周
(ASTRID2: 位置 2、x = 0.5〜18 mm、y = 0.5〜14 mm。KEK-LS: 位置 3、x = 0.2〜6 mm、y = 0.2〜8 mm)。誤差の標準 (k = 1): 四極・6 極・偏向磁石 (KEK-LS は 8 極も) の DX・DY 50 μm、四極のロール 0.2 mrad、偏向磁石のロール 0.1 mrad (補正なし。実機の許容値ではなく例)。
誤差なしの D > −5 の割合は ASTRID2 で 0.030 (41 × 31 で 0.038)、KEK-LS で 0.142 (失われた点 125 / 336)。

### 10.1 誤差つき

ASTRID2、種 8 つ (all) または 5 つ (種類ごと)、平均 ± 標準偏差:

| 誤差 | D > −5 の割合 (k = 0.3 / 1 / 3) | 失われた点 / 336 (k = 0.3 / 1 / 3) |
| --- | --- | --- |
| all (四極・6 極・偏向磁石の全部) | 0.288 ± 0.029 / 0.470 ± 0.047 / 0.644 ± 0.039 | 0.1 / 3.5 ± 3.4 / 113 ± 39 |
| 四極のずれ | — / 0.367 ± 0.042 / 0.555 ± 0.041 | — / 0.4 / 33 |
| 6 極のずれ | — / 0.256 ± 0.036 / 0.361 ± 0.029 | — / 0.2 / 1.0 |
| 偏向磁石のずれ | — / 0.399 ± 0.013 / 0.586 ± 0.057 | — / 0.8 / 62 |
| 偏向磁石のロール | — / 0.287 ± 0.014 / 0.495 ± 0.034 | — / 0.0 / 3.2 |
| 四極のロール | — / 0.219 ± 0.030 / 0.327 ± 0.047 | — / 0.2 / 0.4 |

- どの 1 種類でも、×1 で D > −5 は誤差なしの 7〜12 倍。足し算にはならない (飽和)。
- **主因は結合 (差の共鳴 νx − νy = 3 を直接駆動する)**: ASTRID2 の作動点は差の共鳴の 0.04 以内で、誤差なしでも 31% の点が νx = νy にロックされた島の中 (規則的)。1 台の四極のロールだけ (閉軌道なし) で |c−| = 0.0003・0.001・0.003・0.01 のとき D > −5 は 0.134・0.193・0.304・0.424。
  D > −5 の割合は 0.871 + 0.190 log10 |c−| (84 組、r = 0.75)。島の縁・外側にもカオスが広がる (島の中 / 縁 / 遠い: 誤差なし 0.07 / 0.05 / 0.01、四極のロール ×1 は 0.28 / 0.58 / 0.15)。
- **6 極のフィードダウン**: `fmapcheck.feeddown_estimate` (ΔK1L = K2L x_rel、ΔK1sL = K2L y_rel、四極のロール 2K1Lθ) が差の共鳴までの距離を ×1 で 2〜8% (×3 で 20%) で再現する。6 極のずれを閉軌道に合わせる (`center_sextupoles`) と、all ×1 の D > −5 が 0.470 → 0.330、|c−| が 0.0092 → 0.0011。
  閉軌道の寄与 (垂直の閉軌道 y の rms 1〜2 mm) は、同じ |c−| の結合だけより D > −5 を約 +0.1 増やす。
- D > −5 の点は NAFF の見かけではなくカオス (4096 周の追跡で D_1024 ≥ −4 の点の傾き +0.1)。最小振幅の D の下限は結合で上がる (−10.7 → ×1 の平均 −8.2)。
- KEK-LS (補正なし): ×0.03 (1.5 μm) で D > −5 が 0.142 → 0.556 ± 0.036、失われた点 125 → 171。×0.1: 0.557、227。×0.3: 0.450、285。×1: 線形に不安定 (5 種のうち 4 つ)、333 失われる。×3: 閉軌道が求まらない。
  四極のロール 20 μrad だけ (|c−| = 0.0002) で 0.427。6 極・8 極を閉軌道に合わせると失われる点が 227 → 128 に戻る。
- **使い方**: 誤差つきの周波数マップは複数の種の統計で見る。警告 `errors`・`coupling`・`unstable`・`orbit_large` (§3.6)。

### 10.2 挿入光源 (MPW)

`ASTRID2_v01_02_id.sad` の `RINGMPW` (336 点): **失われる点の原因は 4νy = 9 ではない**。νy を 2.1992〜2.2992 (±0.05) に四極で動かしても、失われる点は 207〜252 で変わらない (|4νy − 9|/4 = 0.0008 のときの 245 と同じ)。
MPW の強さ s で y の限界が決まる: y の 50% 以上が失われる最小値は s = 0.5: 10.4 mm、0.75: 6.8 mm、1: 5.0〜5.9 mm。s = 1 で νx・νy を MPW なしの値に四極で戻しても 216 / 336 失われる (s = 0.5 は 94)。
原因は MPW の内部の非線形 (K0 ≠ 0 の偏向磁石 (ANGLE = 0) の本体) による y の力学口径の縮小。MPW つきの y の範囲は 5 mm 以下に取る。

### 10.3 RFSW ON・RAD ON・FLUC ON

RF ON で動くが、D は拡散の指標でなくなる: D > −5 の割合 0.03 (既定) → 0.565 (RF ON・RAD OFF) → 0.985 (RF + RAD) → 0.98〜0.99 (+ FLUC、種 2 つ)。RF なしの RAD ON は δ が 1024 周で −1% 動く。
z が振幅に比例して動き (4 mm、δ 1.6e-4)、δ = ξ × z の変調 (1/νs = 660 周) がチューンを動かす。gs (RFSW ON・RAD OFF、ASTRID2、6 点、512 周) と座標が 4e-16〜2e-11 m で一致。
直したもの: RAD OFF のときの同期位相を 0 に (放射損失がないので。`Emittance[]` の φs では正味の加速が入り、gs と 1% ずれた。0 にして gs とビット一致)。台本 `FrequencyMap[]` の既定を RFSW・RAD OFF に。

### 10.4 K1 の設定

ASTRID2 の偏向磁石は K1 ≠ 0 だが端の角が 0 なので、`set_bend_k1_edge('face')` は効かず D が同じ。**`set_bend_k1_def('lin')` (B' 一様) は D の相関 0.47、点ごとのチューンの差の中央値 4.8e-3** (作動点・最小振幅は同じ)。
端の角を足した変形 (E1 = E2 = 0.02): face が D の相関 0.924、lin が 0.634、k64 が 0.983。K1 ≠ 0 の偏向磁石のあるリングの周波数マップは K1 の定義を記録して使う。

### 10.5 D の閾値・周回数

規則的な点の D は周回数 2 倍で −1.00 ± 0.05 (ASTRID2・KEK-LS、256〜4096 周)。1024 周への換算 `D_norm = D + log2(nturn/1024)` で判定すると、D > −5 の割合の周回数による変化が 7 倍 (固定) から 2 倍に減る
(ASTRID2: 512 → 4096 周で固定 0.065 → 0.009、換算 0.021 → 0.045。KEK-LS: 固定 0.269 → 0.095、換算 0.127 → 0.181)。512〜4096 周で使う (256 周は不可)。点ごとの傾きは診断用。窓なしは不可、Hann p = 1〜3 は同等、格子は 88 点以上で分布が同じ。

### 10.6 gs との比較 (誤差つき)

gs の乱数の種は pysad に渡せない (SAD は自分の乱数列)。pysad の誤差の値を、誤差のある実体ごとの固有の要素定義として gs の台本に書いた (`LINE["DX", n] = 値` は同じ名前の全部の実体 (QF3 は 12 個) に効くので使えない)。
ASTRID2 (all ×1、種 1、6 点、512 周): pysad の追跡と台本はビット一致、gs と規則的な 4 点で 1.6e-15〜1.4e-12 m、D は 4 桁一致、カオス的な 2 点は 5e-5・2e-2 m (1e-13 m の摂動の増幅と同じ)。

### 10.7 S2 の未確認

誤差による島・カオス層の幅の理論との比較、MPW の y の口径を決める高次の項の解析、偏向磁石のフィードダウンの見積もり、誤差つき KEK-LS の補正 (軌道・結合補正) つき、PF リング・SP8 (大きなリング B)、FLUC ON の gs との統計の比較。
