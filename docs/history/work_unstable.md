# 不安定ラインでの CELL CALC の続行と nstab の判定 (2026-10-16)

担当: pysad/lattice.py, ffscmd.py, match6.py (run_regression.py にベンチを追加)。課題は docs/history/work_descend.md の「残った点」:

- pf2 の tune_sv の `USE RING; cell calc emit` で SAD は νy=7.0293 を返して `ny not matched` に進むが、pysad は `periodic_optics` が sin μy≈0 で ValueError、qcod の閉軌道が見つからずスクリプトが止まる。
- matching02_r (解なし) で SAD は r=3.74 < rbest=4.44 の試行点を不採用にする (nstab を先に比べるので試行点が不安定扱い)。pysad は試行点を採用していた。

## 1. 結論

**両方とも解消した。** 不安定なラインの CELL CALC は SAD の qcell1 と同じ値・同じ表示で続行し、GO の nstab 列・反復表も SAD と全行一致する。
原因は 1 つではなく、SAD のソースを読み直して見つけた次の 8 つの不一致の合算だった (どれも不安定ラインで効く。安定なラインでは SAD と pysad の出力は変わらない)。

| # | 不一致 | 場所 | 修正 |
|---|---|---|---|
| 1 | `|cosμ|>1` で cosμ → 1/cosμ に置き換えてチューンらしい値を作る処理と、INTRES/HALFRES による stabx/staby の判定、伝搬後の再判定 (tracex)、`stab .and. codfnd` がなかった (`periodic_optics` は ±1 に丸めて例外) | lattice.py `twiss(sad_stab=)`、`_qmdiag`、`_qcell_periodic_sad` | qcell.f の直訳。判定は `Lattice._qstat` に残す |
| 2 | qtwiss1 の over (伝搬中に βx, βy ≤ 0 になる) で以降の行を直前の正常な行のコピーにして stabx=staby=偽 | lattice.py `twiss()` | `_transport_twiss(clamp=False)` で β の下限 1e-12 の丸めを外して判定 |
| 3 | tfsetupcell が CELL の GO/CALC で条件の最後に `AX(1→$$$)=0, AY(1→$$$)=0` の 2 条件 (周期性) を足す。周期解が安定なら |df|<1e-8 で捨てられるが、1/cosμ に置き換えた「周期解」では 1 周で閉じず、残差と微分に効く | match6.py `go()` | 条件を一時的に追加 (`cellcond`)。`rl=CONVERGENCE·nfcol` も +2 |
| 4 | tdfun は `|mfitp|` が奇数 (mp·2+1=|mfitp|) の条件を idp=0 では評価しない (中心を除く)。tfsetupcell の条件は mfitp=max|mfitp| なので、ユーザー条件がないとき (mfitp=1) は評価されず `Matched.` になる | match6.py `_tffscalc` | 同じ条件で idp=0 を飛ばす |
| 5 | qcod が overflow (一周行列が NaN) で戻ったとき、SAD は twiss を書かず前回の値が残り、残差はその値で評価される (先頭行は MARK = 最後に安定だった周期解)。オフモメンタム点 (idp≠0) も、前の点が over ならその idp の utwiss の前回の値で評価される。pysad は条件を落として r=0 にしていた (r=0 は crate=1 の「大改善」に見えて、FITP の GO で反復の経路が変わる: pf2 の matching_v31_23 の最後の GO) | match6.py `_stale_rows`、`_rows_idp` | 前回の rows で評価。nstab は残差が非 0 の idp ごとに +2 |
| 6 | `log10(0)` を NaN にしていた (Fortran は -Infinity)。勾配が 0 になる点で crate=0 → badcnv=1 → モード切替になるはずが、NaN の fuzz=1 → badcnv=0 で**同じ点を無限に反復**した | match6.py `_log10` | 0 は -inf |
| 7 | `fmincube` の上限 `min(m_euler·f1, …)` の `m_euler` を**ネイピア数 e (2.718)** にしていた。tfstk.f の `m_euler` は**オイラー定数 γ (0.5772)**。直線探索で fmincube の結果が 0.577·f1 を超えるとき (試行点が不安定で r が小さくなりすぎ、探索が伸びるとき) に違いが出る。これが matching02_r の `fact` の列のずれの原因だった (不安定判定だけの問題ではなかった) | match6.py `M_EULER` | γ に修正 |

| 8 | qcod は一周行列を `limitnan(·,1e10)` で ±1e10 に抑えてから Twiss を作る (NaN は +1e10)。一周行列の成分が 1e10 を超えるほど不安定なとき (QF ×1.16: 成分 7e15) の Twiss は抑えた行列から作られる。さらに成分が 1e18 を超えると SAD では COD の丸め雑音 (1e-18) が一周で非線形領域まで増幅され NaN → `*****qcod---> Overflow` になる (pysad の COD は厳密に 0 で収束するので再現できない → 1e18 超を overflow とする目安) | lattice.py `_qcod_periodic` (`_TRMAX`, `_LAST_MMAX`, `_MOVERFLOW`) | QF ×1.16 は NX 62.1509174812、NY .5、BX 1.4605585589、BY 4.141200142E41 まで SAD と一致 |

加えて SAD K64 の癖を 3 つ確認して合わせた (§3)。

## 2. SAD のソースの読み (照合した箇所)

### qcod (qtwiss.f:810)
- 初期値 cod0 = 先頭 MARK の DX DPX DY DPY (twmov(1,..,.true.) = 最後に安定だった周期解)。`itmax=10`、`conv=min(1e-10, 1e-18(1+(dp/1e-3)²+nlat))`。
  `r=(Courant-Snyder 不変量で測った |出口-入口|²) ≤ conv` で `codfnd=.true.`。使い切ると `cod=cod0; codfnd=.false.`。
- Newton の線形解は `tsolvg` = `tsvd(a,b,x,1d-8)`。整数・半整数共鳴の近くで `M−I` が特異に近いとき、numpy の `solve` と暗黙 QR の滑らかな減衰は違う。
  → `lattice._solve_cod`: 条件の良い行列は従来どおり solve (出力はビット単位で同じ)、特異に近い (最小特異値 ≤ 1e-8·最大) ときだけ `tsvd`。
- `over` は qtwiss1(mat) が一周の途中で `ktfenanq(trtr)` (NaN/Inf) を見つけたとき。`itgetfpe()` は SIGFPE を数える仕組みだが、この K64 ビルドは浮動小数例外をトラップしない (`sad_signal.c` は `feenableexcept` を呼ばない) ので常に 0。
- 閉軌道が見つからないときの挙動: `over` なら `*****qcod---> Overflow & closed orbit not found n`、`over` でなく `codfnd` が偽なら `*****qcod---> Closed orbit not found n` を出し (nmes < 4)、
  後者は stab=偽にして **Twiss の計算は続ける**。前者は twiss を書かずに戻る (前回の値が残る)。

### qmdiag (qmdiag.f)
- `(TrA−TrB)² = trab2`。`≥1e-14` の通常の場合は `periodic_r` と同じ式。`|trab2|≤1e-14` は特殊な 6 分岐、`<0` は `***qmdiag---> Sum resonance: (TrA-TrB)^2 = x` を出して stab=偽 (twiss(lb) は上書きされない)。NaN は `Matrix not a number`。
  → `lattice._qmdiag` に全分岐を移植。通常の場合は従来どおり `periodic_optics` を使い (ビット単位で不変)、`|cos|>1` か stab=偽のときだけ新しい経路 (`_qcell_periodic_sad`) を通る。

### qcell1 の cos μ の処理 (qcell.f:209-248, 351-397)
- A, B ブロック (s 行列で一周行列を対角化) の `cosμ=tr/2`。`>1` なら `cosμ:=1/cosμ`、INTRES なら stabx=偽、そうでなければ `tracex=2cosμ`。`<-1` も同様に HALFRES。
  `sinμ=sign(sqrt(1−cos²), a12)`、`ax=(a11−cosμ)/sinμ`、`bx=a12/sinμ`、`μ=atan2(sinμ,cosμ)`。
- 分散は `stabx .and. dcosμ≠0` のときだけ書かれる。不安定な面の EX/EPX (EY/EPY) は**書かれず、以前の値 (先頭 MARK の値) のまま**伝搬される。
  → `lattice._mark_init` (ffscmd が MARK の生の値を入れる) を使う。
- qtwiss で各要素へ伝搬 (位相は伝搬の積算そのまま: 安定のときだけ pysad は位相の小数部を一周行列のモードに揃えていたが、不安定では行わない) → `tracex=x11+x22` (入口と出口の Twiss から)。
  `>2` で `intres .and. .not. fam` なら stabx=偽、そうでなければ `tracex=4/tracex`。`<-2` は HALFRES。最後に `stabx = stabx .and. stab .and. (codfnd .or. fam)`。
- 先頭 MARK に周期解を書き戻す (twmov) のは**全 idp が安定のときだけ**。不安定なら次の初期値 (閉軌道、MARK の Twiss) は前の安定解のまま。pysad も `_cod_prev` と `_mark_init` を安定時だけ更新する。

### qtwiss1 (qtwiss.f:150-180)
- 要素に入る前に `bx0>0 .and. by0>0` を調べ、だめなら (その行より前の最後の正常な行、なければ入口) のコピーで以降を埋めて `over=.true.`。qcell1 は over のとき `stabx=staby=.false.` で即 return。
  伝搬の最後の行が悪くても (その後に要素がなければ) 検査されない。

### tffscalc / tffsa
- nstab (CELL): `residual1(i)%r ≠ 0` の idp i だけ、`.not. stabx` と `.not. staby` を 1 ずつ数える。どの idp も残差が 0 (条件なし) なら nstab=0 で、`cellstab` だけ全 idp の stabx・staby の積。
  INS は `.not. stabx` で +2。`cellstab` は tfshow が見る (偽なら `MKS.1/$$$ BMAGX/BMAGY` の 2 項目を表示に足す)。
- idp≠0 は `stabx = fam .or. qcell1 の stabx .and. (整数部の判定)`。idp=0 は qcell1 の stabx/staby そのもの (整数部の判定はない)。pysad は idp=0 を常に安定にしていた。
- `tffsa.f:1229-1247`: GO/CALC の最後に、全 idp の optstat で stabx が偽なら `Horizontal`、staby が偽なら `Vertical` (両方なら `Horizontal/Vertical`) + ` Unstable.` を出す。
- tfsetupcell (tffsa.f:1664): 上の #3。tdfun (qdcell.f:180) の `cell .and. ka>nfc0 .and. ka<=nfc0+2 .and. |df|<abmax=1e-8 → cycle`。
- tffsmatch: メッセージ回数 `nmes` (cellm の変数、<nmmax=4 のとき出す。反復表を 1 行出すたびと GO の終わりに 0 に戻る)。`lout` は qmdiag が `lfno=0` にするので和共鳴のメッセージは次の表示まで 1 回だけ。
  `rlist(inr)=r%r` (NetResidual)、`rlist(imr)=res2r(r)` (MatchingResidual)、`rlist(isl)=nstab` だが `isl=ktfsymbolz('UnstableLevel',14)-4` の `14` は名前の長さ 13 より 1 大きく、**UnstableLevel 自体は更新されない (常に 0)**。
  → pysad も `UnstableLevel=0` 固定、`NetResidual=r.r` (nstab を含まない)。
- resle (tffs.f:1644): `nstab` が小さい方が常に良く、同じなら `r1%r ≤ r2%r/tol`。pysad の `_Res.le` は同じ (読み直して一致を確認。matching02_r の原因は nstab の数え方・over の扱い・fmincube の上限だった)。

### FFS のフラグ (重要な落とし穴)
- FFS 内の `ON INTRES;` / `OFF INTRES;` は**何も変えず**、`Out[1]:= USE` を出すだけ (`FFSFlags[]` で確認: 値は変わらない)。FFS のフラグは `INTRES;` (ON) / `NOINTRES;` (OFF) と書く (`RFSW; NORADCOD;` と同じ)。
  → ベンチ 28b は `NOINTRES; NOHALFRES;` を使う。(`OFF INTRES` を使うと SAD の結果は 4 通りとも同じになる。最初に作った 28b はこれで無効なベンチだった。)

## 3. pysad の変更

- lattice.py: `_qcod_periodic` が収束 (`_LAST_QCOD_FOUND`)・`_solve_cod` (tsvd 退避)、`_qmdiag` (全分岐)、`_qcell_periodic_sad`、`Lattice.twiss(sad_stab=)` / `twiss_table(sad_stab=)`
  (`sad_stab={'intres':, 'halfres':, 'init':}` を渡すと不安定でも続行し、`Lattice._qstat = {stabx, staby, over, codfnd, tracex, tracey, msgs, stab}` を残す。渡さない従来の呼び出しは不変で、不安定なら RuntimeError)。
- match6.py: `_calc_optics` (メッセージと nmes)、`_optics_all` (idp=0 は qcell1 の判定、前の点が over なら後続は over)、`_tffscalc` (nstab/cellstab の数え方、TRX/TRY、idp=0 の奇数 mfitp の除外、cell 条件)、
  `_stale_rows`、`go()` (周期性条件)、`_log10`、`M_EULER`、`cellstab`/`last_stab`。
- ffscmd.py: `calc(quiet=)` (CELL は `sad_stab` で続行。安定のときだけ MARK へ書き戻す)、`_matcher_optics`、`go()` の後の `print_unstable` (`Horizontal/Vertical Unstable.`)、
  `print_fit_table` (BMAGX/BMAGY の追加、名前の 11 桁目の `f` を現在の FIT 位置の行に付ける)、`UnstableLevel=0`/`NetResidual`、`fn_calculate_optics` (CELL の stab・trace は qcell1 の値)。

表示で変わったこと (SAD に合わせた): 条件表の `f` の位置 (SAD は現在の FIT 位置の行に付ける。13b の `$$$ NX` は FIT MKC の後なので `f` が付かない)、13a-2 の `*****qcod---> Overflow ... 1/2` が出るようになった (SAD の bench_level13a.result の 78-79 行と同じ)。

## 4. ベンチ

`bench_level28a.sad` (ESRF pysad/tests/kek-ls_v03_12.lat 連結、NPARA=1。結果 `bench_level28a.result`):
- V: QDA* ×1.19 → 垂直だけ不安定 (BY=90、NY=30.2494)。CALC → `Vertical Unstable.` と `MKS.1/$$$ BMAGY 124207.53`。
  FIT NX NY / free QDA* / go → SAD の反復表 `5 0 4769 (DESCEND)`, `7 …1041 (NEWTON)` … (反復 1〜4 は nstab>0 で表に出ない。反復 2 で qcod overflow) まで全行一致、K1 は 1e-10 以下。
- H: QFA* ×0.9 → 水平だけ不安定 (BX=27043、AX=-9619)。GO は安定に戻れず `Unstable = 1 Residual = 2124`。反復表 `2 1 2157 (NEWTON)` まで一致。
- HV: QD ×0.95・QF ×0.9 → 両方不安定。7 変数 GO も安定に戻れず `Unstable = 2`。
`bench_level28b.sad`: INTRES/HALFRES の 4 通りで V/H/HV/V2 (QD×1.1: 安定) の CALC、`NOINTRES` では安定扱い (`Unstable.` なし、tracex=4/trace、分散は書かれる)、さらに `NOINTRES` で V の GO、`NOHALFRES` で H の GO。
`bench_level28c.sad` (V の emit のみ): pysad の emit6d は不安定な光学に未対応 (§6)。run_regression には入れない。

`run_regression.py 28a 28b` (`check_level28`): Print[] の行 (相対 1e-8、TWF_H だけ 1e-3、K1F_H 1e-6 — 不安定の縁の解は丸めが拡大される)、反復表、`Unstable.` / `Unstable = n` / qcod の報告の行の列が SAD と一致。

### 照合 (SAD との差)
- tune_sv_v30_24_1: 2 段目のあとの `USE RING; cell calc emit` が `Vertical Unstable.`、`BMAGY 1.069469`、BX 26.664801、BY .241658、NY 7.029326 で SAD と一致し、`ny not matched` に進む (スクリプトは最後まで走る)。
- matching02_r: 反復表 (6 → 4.280、9 → 4.164、12 → 4.072) と最終 `Residual = 4.072` が SAD と全行一致。デバッグ版 SAD ((WSL のホームの作業フォルダ)、`<work-scratch>/addS.sh addF.sh` で nstab と fmincube の入力を印字) で
  試行点の (r, nstab) が 1 反復目から一致 (r=3.7402426870 nstab=2、2.5941223246 nstab=1、…) することを確認した後、fmincube の `m_euler` の違いに行き着いた。
- 不安定の程度を 20 通り (QD/QF のスケール) 変えた CALC: β が 1e5 程度までは NX, NY, BX, BY が 1e-10 で一致、`Unstable.` の判定・メッセージも一致。
  成分が 1e10〜1e16 の強い不安定 (QF ×1.16) も `limitnan` で一致。1e18 を超えると SAD は COD の丸め雑音で NaN/`Closed orbit not found` になり再現できない。
- pf2 の他のスクリプト: matching_v30_12_r (最後の FITP の GO 含む)、matching_v31_23_r の反復表・Matched/Residual は SAD と一致のまま (v31_23 の `Matched. (1.3021E-12)` と SAD の `1.3022E-12` の最下位桁だけ既知)。
  (最初の版は v31_23 の最後の GO、DP=0.001 の FITP で、overflow した idp≠0 の点を条件なし (r=0) にして反復が SAD と違う経路を通り止まらなくなった。overflow した点は前回の rows で評価する (`_rows_idp`) ことで直した。)

## 5. 回帰

- 作業前 (`base_unst_pre/`) と作業後 (`base_unst_post2/`) の 13a/13b/14a/16a/22a/23a の `run_sad.py` 出力: 数値は不変 (違いは条件表の名前の 11 桁目 `f` の位置 (SAD に合わせた)、13a の `*****qcod---> Overflow ... 1/2` の 2 行 (SAD の bench_level13a.result の 78-79 行と同じ)、警告の行番号のみ)。
- pytest (`test_bend_feeddown` `test_coord` `test_sext_radiation` を `--ignore`): **153 pass** (146 + 新規 `test_unstable.py` 7 本)。
- `QUICK=1 python run_regression.py`: **37 ベンチ FAIL 0 (188 s)** (28a/28b を含む)。`run_regression.py 13a 16a`: FAIL 0 (13A-2 の TWB だけ既知 KNOWN、work_integration.md §2)。
- pf2: tune_sv_v30_24_1 (最後まで走り SAD と同じ表示)、matching02_r、matching_v30_12_r、matching_v31_23_r の反復表・Matched/Residual は SAD と一致。

## 6. 未解決・次の作業

1. **不安定な光学での EMIT** (`emit;`、`Emittance[]`): SAD の temit は不安定でも報告を出す (固有値が単位円上にない)。pysad の emit6d は別の値を返す
   (28c: SAD は Emittance Y = -1.177E-29、BY = 1.000000、`Imag.tune: -0 0.9885795 -0`。pysad は 2.886E-6 …)。tecalc の不安定時の扱いの移植が要る。
2. β が 1e30 を超える極端に不安定な光学 (QD ×0.9 など) は SAD の COD に 1e-20 の丸め雑音があり、r=Δcod²·β > conv で `Closed orbit not found` になる (pysad は COD が厳密に 0 で収束)。丸めの再現は不可能なので対象外。
   SAD が overflow を報告する強い不安定 (QF ×1.2) も、pysad は NaN にならずに値を返す。
3. `fam` (FITP の family 点) と idp≠0 の over の扱い (前の点が over のとき utwiss をコピー) は構造だけ移植 (SAD の utwiss(1,ix,·) のコピーは AX の列だけに見える点は未再現)。
4. `OFF COD` (orbitcal=偽) の qcod の分岐、`calopt` 偽の分岐は未対応。
5. `ON INTRES` (MAIN 流の書き方) は FFS 内では無効 (§2)。pysad の `ON`/`OFF` は従来どおりフラグを変える (SAD と違う)。FFS 内の `ON/OFF` の扱いは SAD に合わせて無効にするか未決。

## 7. 道具

- `<work-scratch>/make_28.py`: bench_level28a/b/c.sad の生成。`<work-scratch>/gen_scan.py`、`runsad.sh`: ESRF の K1 をスケールして CALC の結果を並べる探索。
- デバッグ版 SAD の追加の write 文: `<work-scratch>/addS.sh` (DBGS: r%nstab と rbest%nstab)、`addF.sh` (DBGM: fmincube の入力 f1 f2 g1 g2 rp0 dg)。
  `<work-scratch>/runsad.sh <sad> <名>` で実行。`PYSAD_MATCH_DEBUG=1` の pysad の印字と反復ごとに比べる (pysad の `iter=` は SAD の `DBGA iter` より 1 少ない)。
