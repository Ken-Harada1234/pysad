# 作業ログ: RAD+FLUC の損失周の乱数消費 (その 2) — DAsv small の粒子ごとの突き合わせ

担当ファイル: `pysad/track6n.py`、`pysad/track6.py` (`pysad/track6n_ext.py`、`pysad/dapert.py` は変更なし)。作業用ファイルは `<work-scratch>/`。
日付: 2026-10-16。前担当のログ: `work_loss_rng.md` (単粒子 12 ケース一致)、`work_dasv_timing.md` (実寸 DAsv の計測)。

## 0. 結論 (先に)

1. **小規模 DAsv (153 粒子 × 100 周) で SAD と pysad の乱数列が最初にずれるのは 1 周目の要素 600 (BPI3H の 2 つ目) で、落ちる粒子 1 個 (ID 98: z=0.0 行の 48 番目、x=376。もう 1 個 ID 46: z=−0.5 行の 47 番目、x=368 も同じ要素で光子を引く) の通り道が違うため。**
   その粒子は 592 番の tapert を通り抜けたあと 596〜600 で |px|>1, 座標が 1e6→1e10→1e80 と発散し、要素 600 の
   `tbendcore` の `da = asinz(...)` が −phi0 ちょうどになるか ±1 ulp ずれるかで、その後の `y3, py4 (7.7e110 対 3.6e151)` が 40 桁変わり、
   `tradkfn` が光子を出す (θ=π/2, 乱数 225 個) か出さない (NaN, 乱数 1 個) かが決まる。
2. **この分かれ目は pysad のコードの論理の違いではなく、丸め誤差の 1 ulp の差が増幅された結果である。** 証拠:
   - 要素 600 の入力を SAD のデバッグ版が印字した値 (17 桁) にすると、pysad のコアは SAD の出力と 12 桁一致し、seed も一致する (3140490377)。
   - 各要素の入力を SAD の状態に毎回合わせ直して 1 要素ずつ比べる再同期テスト (`<work-scratch>/r2_resync.py`) では、
     **30 周 (93 030 要素・粒子集合の記録) のすべてで要素後の seed が SAD と完全一致**し、座標が合わないのは 4 記録だけ (どれも落ちた粒子の NaN と inf の違いか、丸めの 1 ulp が増幅されたもの)。
   - SAD 自身が、1 周目の初期座標に相対 1e-16 の摂動を入れただけで最終 seed を変える (12 回中 2 回)。1e-12 で 8/12、1e-10 で 12/12。
     修正後の pysad の 1 周目の最終 seed (3224547513) は、SAD が 1e-16 の摂動で実際に返す値の 1 つと同じである (下記 §3.5)。
3. **自由に走らせた (再同期なしの) 乱数列が SAD と最終 seed まで一致することは、この判定基準のままでは達成できない。**
   pysad (numba) と SAD (gfortran -Ofast + glibc) では演算順序 (-Ofast の再結合・逆数化) と libm (sinh は 25% の引数で最下位ビットが違う、sin/cos/asin/exp が 0.2〜1%) が違い、
   1 要素ごとの差は 1〜数 ulp (規格化した誤差で 1e-16〜7e-15) あり、これが 600 要素の間に 1e-11 まで増える。落ちる粒子は発散するまで乱数を引き続けるので、
   1 個の粒子で分かれ目が逆になると以後の seed がすべて変わる。
4. 再同期テストで **論理の違いが 3 件**見つかり、直した (§4)。直したあとの小規模 DAsv の地図・Score は SAD と一致のまま (Score 63、nsc=21×3)。
   1 周目・2 周目の損失周のリストは SAD と完全一致するようになった (修正前は 2 周目で粒子 84 が 1 周違った)。最終 seed は一致しない (カオス)。
5. 既存の一致 (bench_level12_track_fast、16g_np1、16k_32_np1、16k_n12t20_np1、16j の最終 seed、16l の 12 件、`run_regression.py` QUICK) は壊れていない (§5)。

したがって「損失周の乱数消費のずれ」の論理的な原因はこれで出尽くしたと見てよい (再同期テストで 30 周分が一致)。残りは SAD 自身も 1e-16 の摂動で変わる種類の結果で、
判定基準は「再同期テストの seed 一致」に置くのが妥当 (§6)。

## 1. 作った道具 (`<work-scratch>/`)

| ファイル | 内容 |
|---|---|
| `turns_a.sad` | small の DAsv を `SeedRandom[17]; DynamicApertureSurvey[range, nt]` で nt=1..30 まで回し、周ごとの最終 seed と全粒子の損失周を印字 (SAD 用、pysad でも同じファイルを `run_sad.py` で実行) |
| `cmp_loss.py` | 上の出力を SAD と pysad で突き合わせ |
| `r2_elem.py` | pysad を `track6n.TRACE` フックで 1 要素ごとに追って seed・生存数・tapert を印字、初期座標を `r2_beam.txt` (SAD の beam 行列の形) に書き出す |
| `dest_t1.sad` | その beam を `TrackParticles[beam, dest, 1, nturn]` で dest=2..3564 まで 1 つずつ走らせ、要素ごとの seed と生存数を印字 (SAD の通常版) |
| `cmp_dest.py`, `cmp_coord.py`, `cmp_coord2.py` | 要素ごとの seed と座標 (tapert のクランプを pysad 側にも掛けて) の突き合わせ |
| `<work-debug-patch>/*.diff` | **SAD のデバッグ版** (WSL の (WSL のホームの作業フォルダ) に SAD 一式をコピーして `make exe`。リポジトリの SAD には触っていない) の差分。`tturn.f` に各要素の後の全粒子 (id, x, px, y, py, z, g, dv) と seed (`DBGA`/`DBGN`)、`tbend.f` に |px|>1 の粒子の tbendcore の入力と出力 (`DBGI`/`DBGC`)、`photons.f` に tradkfn の ppx, ppy, ppz, ppa, θ, anp (`DBGR`)、`random_plugin_sad.c` に seed を返す関数 `dbgseed_` |
| `dbgA.sad` / `dbgA3.sad` / `dbgA10.sad` / `dbgA30.sad` | デバッグ版で 1/3/10/30 周を走らせる SAD スクリプト (出力 `dbgA*.result` は最大 1.4 GB: 30 周) |
| `r2_resync.py` | **再同期テスト**。SAD の各要素後の状態と seed (DBGA/DBGN) を読み、pysad の次の要素 (と間の MARK、周末の tapert) を 1 つだけ走らせて SAD の次の記録と比べる。記録ごとに SAD の状態に合わせ直すので、丸め誤差の蓄積を除いた論理の違いだけが残る |
| `r2_nf.py`, `r2_core_locals.py`, `r2_core_sadin.py`, `r2_tradk54.py` | 1 粒子・1 要素を純 Python (numba 無効) で走らせ、最初に非有限になる変数や tbendcore の中間変数を印字 |
| `r2_elem_test.py`, `r2_exact.py`, `r2_noise.py`, `r2_quad2.py` | NORAD/RAD(NOFLUC) の 1 要素ごとのビット一致率・誤差の統計、四極の演算順序の探索 |
| `pert1.sad` | SAD に与える初期座標へ相対 δ の摂動を入れ、1 周後の seed がどう変わるかを測る (§3.5) |
| `ftest/libm.c`, `ftest/cpow.f` | glibc と Python の libm の差、gfortran の複素数べきの NaN の出方の確認 |

## 2. 切り分けの過程

### 2.1 周ごとの DAsv (手順 1)

`SeedRandom[17]; DynamicApertureSurvey[range, nt, Output->0]` を nt=1..30 で SAD と pysad に走らせて (`turns_a.sad`) 最終 seed と 153 粒子の損失周を比べた。

| nt | SAD の最終 seed | pysad の最終 seed (修正前) | 損失周のリスト |
|---|---|---|---|
| 1 | 2957376101 | 3307310777 | 一致 (42 個が 1 周目で落ちる) |
| 2 | 1087225709 | 3177249529 | 1 個違う (粒子 84: SAD 2 周、pysad 1 周) |
| 3 | 2450521345 | 902067809 | 3 個違う |

最初にずれるのは **1 周目** (損失周のリストは同じでも seed が違う)。nt=1 の SAD の seed は `TrackParticles[beam, 3564, 1, 1]` (同じ 153 粒子を pysad の座標で) の最終 seed と同じで、
粒子の初期座標が pysad の dapert の式と SAD の trackd でビットまで同じでなくても 1 周目の結果は変わらなかった (§3.5 のとおり、これは SAD 側の運)。

### 2.2 要素ごとの seed (手順 2)

`dest_t1.sad` で 1 周目の要素ごとの seed を SAD から取り、pysad (`r2_elem.py`) と突き合わせると **DEST 601 (= 要素 600 BPI3H の直後) が最初の不一致**
(SAD 3140490377 対 pysad 623333169)。DEST 600 までは seed も、tapert を考慮した全粒子の座標 (alive な粒子は 6 桁以上) も一致している。

要素 600 で消費する乱数は SAD 375 個、pysad 577 個 (差 202 個)。pysad で粒子ごとの乱数の消費数を数えると (`r2_slotdraws.py`) ID 98 が 225 個、ID 46 が 211 個で、
どちらも光子を約 100 個ずつ出している (tdusrn: θ=π/2 → anp=94)。SAD ではこの 2 個の粒子は光子を出していない。
座標を比べると、他の全粒子は ID 98 の次から (スロット順に) ずれていた (乱数が 202 個分ずれるので当然)。

### 2.3 LOSSAMPL=1e100 にすると全要素の seed が一致する

`LOSSAMPL=1e100` (tapert のクランプ ±3·LOSSAMPL とアパーチャ判定を実質無効にする) で SAD の `TrackParticles` と pysad を走らせると、要素 592〜602 のすべての seed が一致した。
つまり差が出るのは「落ちる途中の粒子 (px>1 を tapert でクランプされた粒子)」の通り道だけで、通常の追跡には問題がない。

### 2.4 SAD のデバッグ版で落ちる粒子の生の値を見る

`TrackParticles` の出力は最後の tapert でクランプされるので、発散中の粒子の本当の値が見えない。そこで SAD を WSL にコピーして ((WSL のホームの作業フォルダ))
`tbend.f`/`photons.f`/`tturn.f` に write 文を足し (`<work-debug-patch>/`)、要素 600 の tbendcore の入力 (`DBGI`)・出力 (`DBGC`)・tradkfn (`DBGR`) を見た。
ID 98 (デバッグ版のスロット 55) は次のとおり:

| | x | px | y | py | z |
|---|---|---|---|---|---|
| 入力 (SAD) | −5.2975135803127869e10 | −7.4422428241049072e10 | 3.1511578250201099e11 | 3.7214853648150600e15 | 7.0959051315318389e21 |
| 入力 (pysad) | −5.2975135802725601e10 | −7.4422428239905746e10 | 3.1511578249556561e11 | 3.7214853647380695e15 | 7.0959051313721031e21 |
| 出力 (SAD) | +4.39e99 | −2.31e80 | 1.13e25 | **3.5754783338584036e151** | 1.01e180 |
| 出力 (pysad) | −2.32e79 | −2.31e80 | (y のまま) | **7.67e110** | 6.04e125 |

入力の相対差は 1e-11 (7.6e-12〜2.7e-11)。SAD の出力側では ppy = dpz·px0 − … が 1e300 を超えて inf になり (DBGR の ppx, ppy = Infinity, ppa = NaN → θ = −π/2)、光子が出ない。
pysad の出力側では ppy = 3e250 の有限値で ppa = 3e250 → θ = π/2 → anp = 94.4 → 約 100 個の光子 (225 個の乱数)。

### 2.5 入力を合わせると pysad も同じになる

SAD の入力 (17 桁) をそのまま pysad の `tbendcore1` に入れると、ID 98 の出力は x, px, y, py, z の 5 つとも SAD の値と **12 桁一致** (`r2_core_sadin.py`)。
その要素を `tradk` まで含めて走らせた再同期テストでも seed は 3140490377 で SAD と一致する。論理は同じで、**入力の 1e-11 の差だけが原因**。

差の効き方: `tbendcore` の `da = asinz((…)/ph2²)` は巨大な項の打ち消しで ±(phi0 ちょうど) になる。pysad の自由走行では da = −phi0 (phi0a = 0, y3 = y1) だが、
SAD では ±1 ulp ずれて phi0a ≈ 1e-18、y3 = y1 + py2·rhoe·phi0a = −1.3e38 となり、f ∝ y3³ を通して py4 が 1e41 倍変わる。
SAD の入力を 1e-16〜1e-10 の相対ノイズで揺らして pysad のコアに入れると (`r2_branch.py`)、phi0a=0 の側に落ちる確率は 9% (1e-16)、40% (1e-14 以上) で、運で決まる。

### 2.6 再同期テストで論理の違いを全部洗い出す

2.5 の考え方を全要素・全粒子に広げたのが `r2_resync.py`: SAD のデバッグ版が各要素の後に印字する全粒子の状態と seed を読み、
その状態から pysad の次の要素を 1 つだけ走らせて SAD の次の記録 (個数・seed・座標) と比べる。要素ごとに SAD に合わせ直すので誤差が溜まらない。

修正前の 1 周目 (3101 記録) は 4 要素で論理の差があった (すべて発散した粒子):

| 要素 | 粒子の入力 | pysad | SAD | 原因 |
|---|---|---|---|---|
| 531 QFA1H | x=1.3e16, y=9.8e11 | x=y=NaN | x=y=inf | `rot1` が無条件に `st*y` を掛ける (0·inf=NaN)。SAD の tsolrot は chi3≠0 のときだけ回す |
| 668 OCFA | x=4.4e3, y=−1.3e5 | g=NaN | g=gmin=−0.999999 | tradkfn の `max(gmin, g+dg)` が NaN で gmin を返す (pysad の fmax は NaN を返した)。g が NaN だと**次の要素の tdusrn が anp=NaN で 10000 回ループして乱数を 20001 個引く** (SAD は anp<0 で 1 個) |
| 3527 SFA | x=−6.5e9 | 同上 | 同上 | 同上 |
| 1558 OCFA | x=−2e6 | py=NaN | py=−inf | 複素数の積・割り算 (ttfrin 出口) の inf/NaN の出方 (gfortran -fcx-limited-range と Python/numba の複素数の違い)。px=NaN なので以後の結果は同じで、tapert でも同じく落ちる |

修正後は **3 周 (9 303 記録) で 4 記録、10 周 (31 010 記録) で 4 記録、30 周 (93 030 記録) で 4 記録**しか差がない (その 4 つは同じもので、後述)。seed が違う記録は 0:

| 周・要素 | 内容 |
|---|---|
| 1 周目 600 BPI3H | §2.4 の粒子。SAD の入力を使っても座標の相対差 4e-9 (入力 17 桁の最後の桁の丸めが増幅されたもの) |
| 1 周目 1558 OCFA | py=NaN 対 −inf (上の表) |
| 2 周目 134 OCFA | py=NaN 対 inf (同種) |
| 3 周目 1769 BID3H | 発散した粒子で asinz の ±1 ulp (§2.5 と同じ) の差が 1e54 倍に増幅 |

要素の種類は QUAD (tsolqur/tsolqu、フリンジ)、薄い SEXT/OCT (tthin、tradk)、BEND (tbendcore, tbshift)、CAV はこのラティスの 1 台なので 1 周に 1 回だけ。

### 2.7 pysad と SAD の 1 要素ごとの誤差 (カオスの種)

`r2_noise.py` (RAD FLUC・再同期、落ちない粒子 |x|,|px|,|y|,|py|<0.05) での 1 要素ごとの誤差 (座標ごとにその要素の最大値で規格化):
ほぼすべて 1e-16〜4e-15 (z の BEND と四極の px, py)、OCFA の py だけ 1.3e-13 (py が 1e-19 程度の小さい値で、もとの丸めが相対的に大きく見えるもの)。
NORAD の `r2_exact.py` でビット一致する割合: ドリフト・MARK・HD は 100% (VD は 99%)、四極 4〜13% (px, py, z)、BEND 9〜47% (x, z。px, py は 99% 以上一致)、六極 19〜20%、八極 16%。

原因は 2 つ (どちらも「pysad の式の間違い」ではない):
1. **演算順序**: gfortran -Ofast (`-freciprocal-math`, `-fassociative-math`) は `a/定数` を `a*(1/定数)`、同じ除数が複数回なら逆数を共有、ループ不変の積 (`w1*s1`, `s1/w1`) を巻き上げ、`pxi + pxi*dc1` を `pxi*c1` にまとめる。
   四極の tsolqu (ndiv=1) で SAD の出力を再現する形を探すと、x は 1230 サンプル中 3 個 (NaN) を除いて完全一致 (`sw = s1*(1/w1)` を巻き上げた `x + (xi*dc1 + pxi*sw)`)、
   px は 26 個 (`pxi*c1 − xi*(w1*s1)`)、y は 53 個、py は 284 個を除いて一致 (y, py は下の libm の差が残る)。これを全要素に実装すれば一致率は上がるが、下の 2. が残る。
2. **libm**: Windows の CPython (numba も同じ CRT) と glibc の最下位ビットが違う割合 (20000 点): sin 0.7%、cos 0.9%、**sinh 25%**、cosh 2.1%、exp 0.2%、log 0.04%、asin 1.1%、atan 0.7%。
   SAD の `xsincos`/`xsincosh` は |x|<0.0235 以外で sin/cos/sinh/cosh を呼ぶ (四極の phi=0.05〜0.2 がこれ)。ビット一致を目指すには glibc の expm1/exp/IBM sin-cos を移植する必要がある。

## 3. SAD 自身の感度 (カオスであることの確認)

### 3.5 初期座標に相対 δ の摂動を入れたときの 1 周後の seed (SAD、`pert1.sad`、各 δ で 12 回)

| δ | 摂動なしの seed (2957376101) と同じ | 異なる seed の種類 |
|---|---|---|
| 0 | 12/12 | 1 |
| 1e-16 | 10/12 | 2 (3224547513 が 2 回) |
| 1e-14 | 10/12 | 2 (3224547513 が 2 回) |
| 1e-12 | 4/12 | 4 |
| 1e-10 | 0/12 | 7 |

修正後の pysad の 1 周目の最終 seed (`turns_a2.py.result`, nt=1) は **3224547513**。SAD が 1e-16 の摂動で返す 2 種類のうち 1 つとちょうど同じ (修正前は 3307310777)。
pysad が SAD の「もう一つの結果」に着地しているので、これは論理の違いではなく丸めの運と見てよい。

## 4. 修正 (`pysad/track6n.py`, `pysad/track6.py`)

1. **`rot1` (track6n) と `_rot` (track6)**: `st == 0 and ct == 1` (th2=0) のときは何もしない。SAD の `tsolrot` は `chi3 /= 0` のときだけ回す。
   無条件だと x=inf, y=inf の粒子で `0·inf = NaN` になり SAD (x=inf のまま) と変わる (要素 531)。track6n_ext の `rot1` 呼び出し (MULT, ソレノイド領域) もこれを使うので同様。
2. **`gfloor` (track6n)**: tradkfn/tradkn の `g = max(gmin, g+dg)` は、SAD (gfortran -Ofast) では `g+dg` が NaN のとき gmin を返す (再同期テストの要素 668, 3527 で SAD の g が −0.999999 になった)。
   `fmax(GMIN, v)` (NaN を返す) をやめて `v if v > GMIN else GMIN`。g が NaN のままだと次の要素で tdusrn が 10000 回ループして乱数を 20001 個引き、seed が大きくずれる。
   前担当の `fmax/fmin` の `NAN_MODE` (a が NaN のとき b を返す) は他の箇所 (sqrt1/pxy2dpz の `max(1e-200,·)` など) では SAD と合っているので残した。
3. **`xsincos`/`xsincosh` (track6n)**: `mathfun` の直訳 (|x| < xthcs1, xthcs2 の 2 段の級数と、それ以外の sin/cos/sinh/cosh、c>0 のとき dc = −s²/(1+c)) に置き換え。
   以前は `sin`, `-2sin²(x/2)` と `sinh`, `cosh` だけで、系統的に 1e-16〜1e-15 SAD と違っていた。級数の `/5040`, `/720`, `/39916800`, `/3628800` は gfortran -Ofast と同じく逆数の掛け算にした (`R5040` など)。
4. **`TRACE` フック (track6n)**: `turn()` の要素ごとに `TRACE(l, prog, S, n, la, rs)` を呼べる (既定 None、デバッグ用)。切り分けスクリプトが使う。

## 5. 照合結果 (手順 3)

| 項目 | 結果 |
|---|---|
| small の地図・Score | SAD と一致 (Score 63、nsc=21×3、地図 3 行は SAD と同じ傾向で、損失周の細部は粒子ごとに違う) |
| small の最終 seed | **不一致** (SAD 966470477、pysad 620061457)。§2〜3 のとおりカオスで、再同期テストでは 30 周分の seed が全記録で一致 |
| small の 1・2 周目の損失周のリスト | 修正後は SAD と完全一致 (修正前は 2 周目で粒子 84 が違った) |
| 再同期テスト (SAD デバッグ版との全要素比較) | 1 周 3101 記録、3 周 9303、10 周 31010、30 周 93030: seed が違う記録は 0、座標・個数が違う記録は 4 (§2.6) |
| 実寸 DAsv (`bench_level20_dasv_full.sad`、手順 4、`work_dasv_timing.md` §6) | Score 273 (SAD 260、10/15 は 275)。nsc は −43.18 行を除く 16 行で SAD と一致 (10/15 は 13 行)、ntloss の完全一致 618/867 (614)、生死の反転 1 粒子 (4)、最終 seed は不一致 (SAD 1343530601、pysad 2334954089)。pysad 1214 s (SAD 785 s の 1.55 倍) |
| bench_level16f_np1 | Score 63 は一致、最終 seed は不一致 (SAD 3824086097、pysad 4221456353)。同じ理由 |
| `bench_level12_track_fast.py` | 変化なし (≤4.7e-11、F T1/T20 の乱数一致、track6 対 track6n 5e-18) |
| bench_level16g_np1 / 16k_32_np1 / 16k_n12t20_np1 / 16j | 最終 seed まで SAD と完全一致のまま (221051865 / 2211729449 / 2926385489 / 17) |
| 16l の単粒子 12 件 (`<work-scratch>/l16l_check.py`) | 10/12 一致 + 2 件は SAD の seed 表示の末尾の丸めだけ (前担当と同じ) |
| `run_regression.py` (QUICK=1) | 28 ベンチ中 FAIL 1 (3C の Emittance[] EMITY の絶対差 5.3e-16 対 許容 3e-16; 追跡を使わない emit6d のベンチで、今回の変更とは無関係) |
| pytest | test_track6n_ext 8 件、test_pysad 49 件、test_item4 6 件 pass。test_new_features/test_error_spreading はスクリプト形式で ALL TESTS PASSED。test_coord (ImportError: `Coord`)、test_sext_radiation (`track_tthin` の引数)、test_bend_feeddown (6 FAIL) は tracking.py/elements.py 側の既存の不具合で今回の変更とは無関係 |

## 6. 今後の方針

- **判定基準**: 「小規模 DAsv の最終 seed 一致」は SAD 自身も 1e-16 の摂動で変わるので基準にならない。再同期テスト (`<work-scratch>/r2_resync.py`、SAD デバッグ版との要素ごとの突き合わせ) の seed 全一致を基準にするのがよい。
  新しい要素や分岐を足したときは、同じ手順 (§1 のデバッグ版 + dbgA*.sad + r2_resync.py) で確認できる。
- ビット一致まで詰める場合の筋道: ① gfortran -Ofast の演算順序を要素ごとに再現 (四極は §2.7 の 1. で探索済み、BEND/六極/八極/tradk も同じ方法で探せる)、
  ② glibc の expm1, exp, sin, cos, asin, atan を numba に移植 (sinh/cosh は expm1 だけで済む範囲 |x|<0.35 が四極では大半)。
  しかし 1e-16 の摂動で SAD 自身が 17% の確率で別の結果になることを考えると、完全一致の保証にはならず、工数に見合わない。
- 実寸 DAsv の 4/867 粒子の生死の違い・最終 seed の不一致 (`work_dasv_timing.md`) も同じ性質のもので、(a) の「別の要素での発散」は再同期テストで見つからず、(b) の「丸めの差がカオス的な粒子で増幅」が主因。
