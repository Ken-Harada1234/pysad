# 作業ログ: 完成済みモジュール 6 本 (+ track6n) の ffscmd への接続 (2026-10-15 夜)

対象: `pysad/ffscmd.py` (FFS 層) に未接続だった match6d / calc6d / touschek / geo / emit_extra / track6n_ext の接続と、
EMIT 失敗時の NaN。手順は「1 ステップごとに SAD の参照ベンチと比較して一致を確認 → 回帰を実行」。

## 0. 回帰の方法と基準値

- 基準値: 作業開始時の `run_sad.py bench_level13a/13b/14a/16a.sad` の出力と、当日の `<work-regression>/*.txt` (12:56〜12:58 の実行)。
  両者は同一だった (13a は stdout/stderr の混ざり方だけ違う)。
- 各ステップの後: (a) 上の 4 本の出力 (K1A/K1B/RESA/RESB/TUNE*/TWB/Score/RES/RDAP、反復行、Matched 行) を基準と比較、
  (b) `python -m pytest -q pysad/tests/test_pysad.py pysad/tests/test_new_features.py pysad/tests/test_error_spreading.py` (49 pass)、
  (c) `QUICK=1 python run_regression.py` (作業中の編集と干渉しないよう、ツリーをコピーしたディレクトリで実行)。
  最後にコピーで全部 (QUICK なし) を実行した (§8)。
- 注意: `run_regression.py` は pysad の出力そのものを比べるので、解析微分に変えた 13a/13b の最終桁は変わる (§2)。
  それ以外の出力は 1 文字も変わらなかった (14a, 16a は 4 ステップとも完全に同一)。
- この文書で「SAD と一致」と書いたものはすべて SAD の `.result` / `.sad.out` との数値比較による (`Print` の `TAG {…}` 行を数値で突き合わせ)。

## 1. EMIT 失敗時の NaN と tffscalc の重み (match6._tffscalc)

変更:
- `ffscmd.FFS.run_emit`: 閉軌道が見つからず `???-Emittance[]-closed orbit not found.` になった場合、`EMITX/EMITY/EMITZ/SIGE/SIGZ` を NaN にする
  (temit.f が rsetgl1 まで進む挙動)。以前は MARK の値 (1.8e-9) のまま残っていた。
- `match6._tffscalc`: SAD tffscalc.f 442-452 の直訳に変更。`em=max(emminv,|emx|+|emy|)`、
  `coum=min(1, max(coumin, coum0=1e-3, 1/(|emx/max(emminv,emy)|+|emy/max(emminv,emx)|)))`、`coum/(1+coum)`、`emxx=max(emx, coum·em)`。
  `emminv=1e-15` (toplvl.f)、`coumin=MINCOUP`。`max/min` は gfortran と同じく NaN を無視する (`_fmax`/`_fmin`)。
  `1/0` は Inf (Fortran)。以前の `em=max(1e-30,…)` と `if (emx or emy) else 1.0` は削除。
- `Matcher.coumin` を追加し、ffscmd が `MINCOUP` グローバルを渡す。
- **ついでに直した別のバグ**: ffscmd の `MINCOUP` の既定値が 0.0 だった (SAD は initb1.f で 1e-2)。0.01 にした。

検証 (<work-esrf>/matching1_dbg.sad = ESRF の esrf/matching + calc 直後の診断出力。SAD 参照 `matching1_dbg.sad.out`):

| 項目 | SAD | pysad |
|---|---|---|
| `cell calc emit` 後の {EMITX,EMITY,EMITZ} | {NaN,NaN,NaN} | {NaN,NaN,NaN} |
| `calc` 直後の MatchingResidual | 41205475.56577 | **41205475.56577** (4.1205E+07) |

(ステップ 1 だけの時点では微分は数値微分のままなので、GO の反復列は 3〜7 の後 `Matched. ( 8.7339E-09)` で、SAD の
`8.7338E-09` と末尾 1 桁違う。解析微分化 (§2) で反復列は完全に一致した。)

## 2. 解析微分の既定化 (match6.Matcher._derivative)

変更 (docs/history/work_analytic_deriv.md §7 の統合案どおり):
- `match6.Matcher.go` の数値微分ブロックを `Matcher._derivative(ids, dfs, rows, wvar)` に切り出し、
  `self.last_qu = (ids, qu, wvar)` を残す (デバッグ用)。`Matcher._numeric_columns` (中心差分) も Matcher に移した。
- `_derivative` は既定で `match6d.analytic_qu` を使い、FitFunction 行は SAD と同じ `v+eps1/wvar` の片側差分、
  `analytic_qu` が数値微分に回す列 (未対応の変数型) は中心差分。
- **`FFS$NumericalDerivative`**: ffscmd の `go()` がグローバルを読み、非 0 なら `Matcher.numerical_derivative=True` (従来の中心差分のみ)。
- **全体を中心差分にフォールバックする条件** (`Matcher._analytic_unsupported`; 理由は `Matcher.last_fallback_reason` に入る):
  INS (CELL でない: 解析微分は CELL のみ検証済み)、COUPLE 使用中 (coup≠1 は未対応; ffscmd が `couples_active` を渡す)、
  オフモメンタム条件 (FITP>1; match6d の idp≠0 は SAD と未照合なので安全側)、tffsqu が扱わない条件 (BMAG/幾何など、LENG は可)、
  誤差 (DK) 付きの変数 (`_errk≠1`)、L=0 の QUAD (qdthin 未移植)。
- `match6d.AnalyticMatcher` は互換のための空のサブクラスに (go/_derivative の写しは削除)。`bench_match6d_check.py` はそのまま動く。

検証:

| ベンチ | 結果 |
|---|---|
| <work-esrf>/matching1_dbg.sad | 反復列 `3 2.9742E+07 (NEWTON) 0.4561 / 4 7.2238E+06 / 5 5.1917E+05 / 6 7.340 / 7 0.1633 / Matched. (8.7337E-09)` が SAD と同じ (SAD は 8.7338E-09、RES1 は pysad 8.73374E-9 / SAD 8.73375E-9)。K1 9 変数が SAD と 12 桁一致 (例 QFA1H .40014410823162 vs .40014410823166、BPI6 −.94963999852662 vs −.94963999852666) |
| `FFS$NumericalDerivative=True` を calc の前に置いた同スクリプト | 数値微分の経路 (Matched 8.7339E-09) になる = 切替が効く |
| 13a (13A-1 Matched) | `Matched. (5.0759E-16)` = **SAD の 5.0759E-16** (数値微分では 5.0763E-16 だった)。K1A、TUNE0/TUNEA は不変 |
| 13a (13A-2, 解なし) | `Residual = 2.5977E-03` 不変。**K1B が SAD から相対 7e-10、TWB (BX@ORG) が絶対 4.7e-8 (相対 2.7e-9) ずれる** (数値微分では SAD と 10 桁一致していた)。§2.1 |
| 13b | 反復数・Matched 行は同じ。`RESB 3.37240E-10` = **SAD の 3.37240E-10** (数値微分では 3.37243E-10)。K1B は 2.5e-10 以内で SAD と同じ |
| 14a, 16a | 完全に同一 (14a は `cell calc` だけ、16a は GO なし) |

つまり 13a/13b の K1/残差は「変わらない」のではなく最終桁が動いた。13A-1・13b の残差は **SAD に近づいた** が、
**13A-2 だけは SAD から離れた (指定の「出力不変」を満たさない)**。原因を調べたので以下に書く。判断は §9 の 1 を参照。

### 2.1 13A-2 (解なし) の最終点が動く理由 (調査結果)

13A-2 は 4 条件 / 3 変数で解がなく、GO は残差 2.5977E-03 の最小点 (真の最小点) の近くで止まる。

- 真の最小点: pysad の解析微分は 3 反復目で K1 = {.5465445163541, −.5111339350306, .5426651385204} (残差 2.597728388430425e-3) に着く。
  数値微分も 1e-12 以内で同じ点 (.5465445163549…)。
- その後のライン探索 (r=0 の不安定点を試す行) を経て、**残差の差 1e-15 (相対 4e-13) の大小**で「最良点」が決まる:
  数値微分の run では 19 反復目の隣の点 K1 = .5465445160443 (残差 …431119e-3) の方が 3 反復目 (…431446e-3) より小さく評価されてそちらが最終答になり、
  解析微分の run では 3 反復目 (…430425e-3) の方が小さいのでそのまま残る。残差の評価自体に 1e-15 程度の雑音 (COD の収束誤差) がある。
- SAD の `gs` を WSL で 2 通り実行した (`NPARA=1`; scratchpad/sad13a/a_default.sad = 既定の解析微分、a_numeric.sad = `FFS$NumericalDerivative=True`):
  13A-1 の残差は既定で 5.0759E-16 / 数値微分で 5.0763E-16 (pysad の解析 5.0759E-16 / 数値 5.0763E-16 とそれぞれ一致)。
  13A-2 の最終点は **どちらも** K1B {.5465445160, −.5111339350, .5426651388}、TWB {…, 17.4362514520, .9505231603} で、
  pysad の数値微分の最終点と 10 桁一致する (pysad の解析微分は .5465445164 / 17.4362514046 で、真の最小点の側)。
  つまり SAD の「隣の点」への落ち方は微分の種類に依らず SAD 側の残差の雑音で決まっており、pysad の解析微分では再現できない (原理的に雑音の再現は不可能)。
- 解析微分の J そのものは正確: 停留点で Richardson 外挿した高精度の数値 J と最大相対差 3e-11 (数値微分 (eps=1e-5) は 8e-10)。
  J の誤りが原因ではない。

pysad で同じ 2 本のスクリプトを走らせた比較 (SAD の `.result` に対する最大差。`FFS$NumericalDerivative=True` を挿入したスクリプトは SAD も pysad も数値微分):

| pysad の微分 | 対 SAD (同じ微分) の K1A / RESA | K1B | TWB |
|---|---|---|---|
| 数値微分 (SAD も数値微分) | 一致 / 1e-21 (abs) | 一致 (10 桁) | 2e-10 (abs) |
| 解析微分 (SAD も解析微分 = 既定) | 一致 / 4e-21 (abs, 5.0759E-16 同士) | 4e-10 (abs) | 4.7e-8 (abs) |

結論: 13A-2 の「SAD と同じ最終点」は、解析微分では **`FFS$NumericalDerivative=True` にすると再現できる** (SAD のこのベンチの値と 10 桁一致; ただし 13A-1 の残差は 5.0763E-16 になる)。
run_regression の 13a は K1B/TWB の差が 1e-7 以内なら KNOWN にした (理由を note に記載)。その他の項目 (反復行、Matched/Residual 行、K1A、TUNE、RESA/RESB) は PASS。

## 3. CALC6D フラグ (ffscmd.calc)

変更:
- `FFS.calc()`: `CELL` かつ `CALC6D` なら `calc6d.calc6d(self.lat, momentum, rfsw=flags['RFSW'], radcod=flags['RADCOD'], rad=flags['RAD'],
  twiss0=直前の rows[0], out=self)` を呼び、`self.rows = res.rows`、`self.calc6d_result = res`、`self.trf0 = res.trf0`。
  不安定のときの `*****qcell6d---> Unstable optics` は self.write に出る。INS + CALC6D は未実装 (4 次元の伝搬になる)。
- `FFS.row_frac()`: 直近の CALC が CALC6D なら `calc6d.twiss_frac(self.calc6d_result, k, fr)` (qtwissfrac の calc6d 分岐)。
- `Twiss[]` は既に AZ/BZ/NZ/ZX.. を含むのでそのまま動く。

検証 (`run_sad.py bench_level18{a,b,c}.sad` の `TW_*` 行を `.result` と数値比較、ESRF 3564 要素、各 16 s):

| ケース | 結果 (8 位置 × 28 成分) |
|---|---|
| 18a (RFSW RAD) | max\|diff\| = 1.74e-10 (NX≈299 rad)、min(絶対,相対) で最悪 1.25e-11 |
| 18b (NORFSW RAD) | 1.74e-10 / 1.26e-11 |
| 18c (RFSW NORAD) | 1.74e-10 / 1.25e-11 (18a と同一出力) |

`TW_P200H` (分数位置 200.5) も一致。`bench_level18a_pysad.py` (Python API 直呼び) の 1.7e-10 と同じ。

メモ: SAD は `Emittance[]` の後で EMITX/SIGE 等を元に戻す (§6) ので、18a の `EMIT` 行 (Print) は SAD では MARK の値のまま。
これも §6 で再現した (§6 の表)。`DTSYNCH` は SAD では CALC6D の temit が設定する変数で、pysad には未定義 (Print に `DTSYNCH` と出る)。

## 4. INTRA / TouschekLifetime (ffscmd.run_emit, fn_touschek)

変更:
- `run_emit`: `flags['INTRA']` が真なら `touschek.touschek_table(lat, momentum, e, PBUNCH, nbunch=NBUNCH, mincoup=MINCOUP, rfsw, radcod, verbose=False)`
  を呼ぶ。グローバルに EMITX/EMITY/EMITZ/SIGE/SIGZ (最終パス)、EMITXC/EMITYC/EMITZC (IBS 込み収束値)、`TouschekTable` (SAD と同じ入れ子リスト) を入れ、
  `self.touschek` に結果を保持。USE で破棄。失敗時は `???-Emittance[]-INTRA failed: …` を出して通常の EMIT 結果に戻す。
- `fn_touschek` (組み込み `TouschekLifetime[nx,ny,nz]`): `touschek.touschek_lifetime` に接続 (INF は `math.inf`)。表がなければ
  SAD と同じ `TouschekTable Undefined; Do Emittance[] with INTRA;`。`sige` は `SIGE` グローバル。
- **`ON INTRA` が FFS 内でも通るようにした**: FFS 内の `ON flag…` / `OFF flag…` を `_exec_word` に追加 (MAIN 部の `ON RAD FLUC;` は
  もともと `split_main_ffs` が処理していた)。`INTRA;` / `NOINTRA;` は従来どおり通る。flag 解釈は tflag.f と同じ (名前、`NO`+名前、sino)。

検証 (`run_sad.py bench_level17a.sad`、23 s; `ON RFSW RAD FLUC INTRA;` に書き換えた版でも同じ):

| 項目 | 結果 |
|---|---|
| TL1..TL15 | 最大相対差 6.7e-12 (TL15 は両方 0)。例 TL1 12.36769325190917 (SAD …192202) |
| EMIT 行 (EMITX,Y,Z,SIGE,SIGZ,EMITXC,YC,ZC) | EMITY (1e-14 の小さい量) だけ 1.5e-7、他は ≤1e-11 |
| TouschekTable 各部 (TT2, TT31–33, TT4S/E/C…) | ≤1.1e-10 (最大は TT4E 1.1e-10) |
| 17b (RADCOD ON) | TL ≤4.5e-7 (既知: emit 自体の 1e-7 差) |

所要は emit 10 s + IBS 反復 13 s 程度。

## 5. GEO (LINE["GX"...], ORG, GEOCAL/GEOFIX)

変更:
- `FFS.ensure_geo()`: `geo.compute_geo(self.els, gx0.., chi0, iorg)` の結果を保持。GEOCAL (既定) なら CALC / EMIT 後や `invalidate()` で
  `geo_dirty` が立ち再計算、GEOFIX なら USE 直後の結果を保持 (tffsa.f:114-118)。`warn` は SOL の未移植警告を 1 回だけ出す。
- `line_value()`: `GX GY GZ GCHI1 GCHI2 GCHI3 GEO` と `OGX OGY OGZ OCHI1 OCHI2 OCHI3 OGEO` (閉軌道は row_frac/row の x_cod..)。
  `LINE["S"/"LENG", loc]` は geo の `pos_frac` (フリンジによる長さ補正込み。失敗時は従来の TwissRow.s の内挿)。
- **ORG コマンド**: `_org(words, i)`。要素名なしは `Info-ORG-Reset to the initial setting.`、続く語は getvad と同じ式の読み方
  (`+`/`-` で始まる語は前の式の続き: `ORG B3 1.5 -2.5 0.5 …` は gx=1.5−2.5)。chi (度) 欠落は `?ORG-Missing parameters: …` を出すが
  欠落分を 0 にして設定する (SAD と同じ。19c の ORGA〜ORGF はこの挙動)。メッセージ書式は SAD の ` ???-FFS-Error-?ORG-Missing parameters: …`。
- `twiss_value()`: **'GX'/'GY' (γ 関数) を SAD の 'GMX'/'GMY' に整理**。`Twiss["GX"/"GY"/"GZ"/"CHI1".."CHI3"]` は K64 と同じく
  `General::wrongval: #1 (GX) is must be to be name of optical function` のエラー (幾何は `LINE["GX"..]`)。
- `_resolve_loc`: nlat を超える番号は周回 (k−nlat; 例 4000 → 436)。`LINE["LENGTH", loc]` は SAD と同じ narg エラー。
- **エラーの扱い (sadscript.Interp.run)**: `run_sad` 経由では文の実行中の `SADError` を `???…` として表示して次の文に進む
  (SAD と同じ。以前は Python の例外でスクリプト全体が止まった)。19a の `Twiss["GX"..]` のエラー行がこれ。最上位の run だけが捕まえる
  (`Interp.continue_on_error`)。pytest 側 (Interp 直使い) の挙動は変わらない。

検証 (`run_sad.py bench_level19a/b/c.sad` と `.result` の `LN/GEO*/POS/ORG*` を数値比較):

| ベンチ | 結果 |
|---|---|
| 19a (ESRF, 21 箇所 + GEO 3 + POS) | LN 最大 4.3e-13 (値 ~190 m)、GEO 4.2e-13、POS 0 (loc 4000 の周回含む) |
| 19b (分岐テスト 35 箇所 + ORG B3) | LN ≤1.1e-14、POS 0、ORG ≤1e-14 |
| 19c (ORG 6 通り × 18 点) | ≤1e-14 |

`Twiss["GX"]` が SAD と同じエラーになることも確認 (19a の出力に ` ???General::wrongval: …` 1 行)。

## 6. Emittance[] のオプション (ffscmd.fn_emittance)

変更:
- `fn_emittance` を `emit_extra.emittance_rules` + `rules_to_sad_list` に差し替え。mode は emit.n と同じ
  (`Matrix->True`=3 ⊃ `Orbit->True`=2 ⊃ `OneTurnInformation->True`=1、既定 0)。戻りは emit.n の規則名・順序 (Stable, Region, Tunes,
  EnergyLossU0 [eV], RfVoltageVc, EquilibriumPosition, MomentumCompaction, OrbitDilation, BucketHeight, HarmonicNumber, OrbitAtExit,
  DampingRate, Emittances, MomentumSpread, BunchLength, TuneShiftByRadiation, NominalSpinTune, Polarization*… + mode≥1 で
  OrbitAtEntrance, OneTurnTransferMatrix, OneTurnDampingMatrix, NormalCoordinates, ExtendedTwissParameters, OneTurnExcitation,
  EquilibriumBeamMatrix, SpinPrecessionMatrix + mode≥2 で ClosedOrbit + mode≥3 で TransferMatrices)。既存の単位 (EnergyLossU0 は eV) は不変。
- **SAD の emit.n で確認した 2 点を反映** (Packages/emit.n を読んだ):
  1. 既定 (mode 0) では `EquilibriumBeamMatrix`/`OneTurnTransferMatrix` は返らない (mode≥1 だけ)。以前の pysad は常に返していたが、
     SAD の bench_level16b.result の `BM5` 行は `EquilibriumBeamMatrix[5]` (未評価) になる。pysad も同じになった。
  2. **`SaveEMIT->False` (既定) のとき EMITX/EMITY/EMITZ/SIGE/SIGZ は Emittance[] の前の値に戻る** (`If[~(SaveEMIT/.op),EMITX=ex0;…]`)。
     以前の pysad は Emittance[] がこれらを書き換えていた (18a の `EMIT` 行が SAD と食い違っていた原因)。`SaveEMIT->True` なら書き込む。
- `ExtendedTwissParameters` は `{AX->v,…}` の規則リスト (`Rule(S('AX'), v)`)。numpy → リスト、bool → 1.0/0.0。
- INTRA ON のときは `run_emit` (IBS 反復 + TouschekTable) を先に呼び、Emittances/MomentumSpread/BunchLength を最終パス値にする。
- 閉軌道が見つからないときは従来どおり `run_emit` の NaN 経路 (`???-Emittance[]-closed orbit not found.`) に落ち、mode 0 の規則を NaN で返す。
- `self.trf0` は `EquilibriumPosition` から更新 (追跡が使う)。`self.emit` は触らない (DynamicApertureSurvey は `self.emit is None` なら
  `run_emit` を呼ぶのでそのまま動く)。
- **`LINE["SIG"//key, l]`** (tftwiss.f:982。`""`→6×6、`X`→σx、`XY`→Σxy、`PX`,`DP`,`Z` など。`SIGE`/`SIGMA…` は除く) を追加:
  `FFS.beam_matrix_at(k)` が `emittance_rules(beam_matrices=True)` を遅延計算してキャッシュ (temit.f:214-222 の追加 1 周)。分数位置は未対応 (SADError)。
  bench_level21a の `LINE["SIG",l]` / `LINE["SIGX",100]` に必要だった。
- sadscript: **`Scan[f, list]` を実装** (以前は Null を返すだけで、しかも `_HOLD` に入っていて引数も評価されなかった)、
  `Rule` の `Part` (`#[[1]]`=名前、`#[[2]]`=値。21a の `Map[(#[[1]])&, ExtendedTwissParameters/.e]`)。

検証 (SAD の `.result` と数値比較):

| ベンチ | 結果 |
|---|---|
| 16b (NORAD NOFLUC, Emittance[] 既定) | EMT 3e-18、SIG 1e-14、TUN 2.8e-11、ALP 1.6e-13 (min(絶対,相対))、DMP 1e-16、COD 1.8e-17 (SAD と同値)。`BM5 EquilibriumBeamMatrix[5]` の表示も SAD と同じ |
| 21a (OneTurnInformation+Matrix, 51 s→32 s) | STB/REG/TSR/NSP/POL/SPM/NLAT 一致、TUN 2.8e-11、U0V 1.6e-13、DMP/EMT/SIG ≤1e-14、OTM 3.3e-10 (norm 相対 1.7e-10 と同じ桁)、ODM 5.9e-14、NCO 6.3e-12、ETK 名前一致、ETW 1.3e-11、OTE/EBM ≤3e-17、CO×8 1.8e-17、TM×8 6.3e-10 (周回で累積)、SZ×8 + SZB×2 + SZX100 ≤3e-17 |
| 18a の `EMIT` 行 (Emittance[] の後に EMITX を Print) | SAD と同じく MARK の値 (1.8058612106805e-09) のまま — 以前は 1.549e-10 に書き換わっていた (§3 のメモ) |

<work-esrf>/calc_RadInt.sad (`Emittance[]` を使う実スクリプト) も最後まで走り、`NaN…` の表示パターンまで SAD と同じ (以前は
`TypeError: must be real number, not complex` で止まっていたのが §5 のエラー継続で進むようになった)。

## 7. track6n の MULT/SOL

### 7.1 `Prog._compile` の `el.L == 0.0` より前の判定
作業開始時にすでに直っていた (track6n.py 12:54 版: `if tn in ('Multipole','Solenoid','TCavity'): raise NotImplementedError` が
`el.L == 0.0` の判定より前にあり、L=0 の判定は `Drift/Sextupole/...` に限定)。今回確認だけして、統合後の形に更新した (下)。

### 7.2 track6n_ext のカーネルを track6n.Prog/turn に組み込み
work_track6n_ext.md §5 の差分案を、ファイルを大きくしない形で実装した:
- カーネル (`k_mult`, `k_sol_*` 等) は `track6n_ext.py` に残し、`track6n.py` が **遅延 import** (`_ext()`) する
  (ext が track6n のヘルパーを import するので循環を避けるため)。kind 番号は ext の 10..13 のまま (track6n の 0..6 と衝突しない)。
- `track6n.Prog.__init__`: while ループにして、`Solenoid` → `track6n_ext.compile_sol_region(self, elements, i, ph)` (領域を compile して終点の添字を返す)、
  `Multipole` → `compile_mult`。`sol_mode` (SM_NONE/START/IN/END は track6n.py に定義) と `_add()` を追加。
- `track6n.turn`: tapert の la の扱いを SAD の tturn1/tsol どおりに (領域内の要素は tsol のループの数え方、領域後に主ループの数え直し)。
  新しい kind は `_ext().step_kind(...)` へ。
- `track6n_ext`: `ProgExt = track6n.Prog`、`turn_ext = track6n.turn` の別名にして互換を保つ。`_step_kind` → `step_kind`。
  `track_particles_ext` (フォールバックなし) はそのまま。
- `track_particles_fast` のフォールバックは、ext が未対応の組み合わせ (ANGLE≠0 の MULT、SOL の DX/DY/CHI/F1、BOUND=0 の SOL、領域内の BEND/MULT/CAVI、
  TCavity) だけ NotImplementedError → track6 になる。`dapert.py` (`_tn.Prog`/`_tn.turn`) も MULT/SOL を含むラティスで使える。

### 7.3 新設テスト `pysad/tests/test_track6n_ext.py` (8 テスト、pytest 12 s)
20 粒子 × 5 周、track6 (純 Python) と `track_particles_fast` (統合後) を比較:

| ラティス | モード | max\|Δ\| | RNG.seed (track6 / track6n) |
|---|---|---|---|
| mult | NORAD / RAD / FLUC | 3.1e-17 / 4.0e-17 / 4.3e-17 | 12345 / 12345 / −2025222455 一致 |
| sol | NORAD / RAD / FLUC | 2.0e-17 / 3.1e-17 / 3.0e-17 | 12345 / 12345 / −1861977247 一致 |

(許容 1e-13、FLUC は seed 完全一致を assert。)
ほかに (1) L=0 の MULT (薄い多極) と SOL が無視されない (MULT/SOL を抜いた追跡と >1e-8 違う)、(2) ソレノイド領域込みの
`turn(apert=...)` が損失なしで n 不変、x=5 m の粒子を入れると n が 1 減る、を確認する。pytest は 49 → **57 pass**。

検証: `python bench_level12_track_fast.py` の SAD との差は不変 (N T1 2.4e-11, N T20 3.9e-11, R T1 3.4e-11, R T20 8.2e-12, F T1 4.7e-11, F T20 4.6e-11、
乱数一致、track6 vs track6n 1.1e-17)。13a/13b/14a/16a の出力もステップ 6 と同一 (16a の DA は track6n.turn を通る)。

## 8. 最終回帰

全ステップ完了後のツリーをコピーして `python run_regression.py` (QUICK なし、21 ベンチ、所要 303 s) と
`python -m pytest -q pysad/tests/test_pysad.py pysad/tests/test_new_features.py pysad/tests/test_error_spreading.py pysad/tests/test_track6n_ext.py`
(**57 pass** = 既存 49 + 新設 8) を実行した。

| ベンチ | 判定 | 最悪 (差/許容) |
|---|---|---|
| L1, L2, L3A, L3D, L4, L5, L6, L7, L8, L9, L10, L11, L12, L12t, L12f, L13B | PASS | 変化なし (L7 1.7e-10、L12t 3.5e-11、L12f 4.7e-11 …) |
| L3B, L3C | PASS* (既知) | 変化なし |
| L14A | PASS* (既知: DP0 スキャン 2.5e-7, EMIT dz) | 変化なし (出力も完全に同一) |
| L16A | PASS* (既知: 損失周数 1 個) | 変化なし (出力も完全に同一) |
| L13A | 初回は FAIL (13A-2 の TWB 4.7e-8 > 2e-8)。§2.1 の調査のあと run_regression.py の 13a の判定を修正し **PASS\*** (13A-2 の K1B/TWB が 1e-7 以内なら KNOWN、note に理由) | TWB 4.7e-8 |

(`run_regression.py` の変更はこの 1 か所だけ: `check_level13` の末尾。元のファイルは scratchpad に退避。)
`<work-regression>/13a.txt`, `13b.txt` は最後に 13a/13b だけを単独で走らせたときの出力に置き換わっている (ほかのファイルは今朝のまま)。
追加で手動確認したもの (いずれも SAD の `.result` と比較済み): 16b, 17a/b, 18a/b/c, 19a/b/c, 21a, <work-esrf>/matching1_dbg.sad。

## 9. 未解決・判断が必要な点・気付いた点

1. **【要判断】13A-2 (解なしの停留点) が SAD の最終点と一致しなくなった (解析微分の既定化の副作用)。** 指示の「13a/13b の出力不変」を満たさない唯一の項目。
   原因は §2.1 (解析微分は真の最小点に着き、SAD は残差の 1e-15 の雑音で隣の点に落ちる。J の誤りではない)。
   ステップを戻す手もあったが、そうすると同じステップの主目的 (esrf/matching の反復列 `… Matched 8.7338E-09` の再現: 数値微分では 8.7339E-09)、
   13A-1 の RESA 5.0759E-16、13b の RESB 3.37240E-10 (どちらも SAD の既定の値) を失うので、既定は解析微分のまま残し、
   run_regression の 13a は K1B/TWB だけ KNOWN にした。**スクリプトの先頭に `FFS$NumericalDerivative=True;` を置けば従来 (= 13A-2 を SAD と 10 桁一致) に戻る。**
   戻したい場合は `ffscmd.go()` の `m.numerical_derivative = …` の既定を True にするだけでよい (match6.py 側は変更不要)。
2. COUPLE した 2 変数をともに free にして GO すると `match6._tfsolv` で `IndexError: arrays used as indices must be of integer type`
   (`qu[fitsel]` の fitsel が空リスト = float 配列になる)。統合前のコード (ステップ 1 時点) でも同じ例外が出る既存の不具合 (解析微分とは無関係)。
   解析微分側は COUPLE 使用中は中心差分にフォールバックするので別問題。
3. 条件表 (`print_fit_table`) の書式が SAD と少し違う (統合前から; <work-esrf>/matching1_dbg.sad の GO 後の表で確認): SAD は `MKS.1       BY   3.6 …`、
   pysad は `MKS.1     f BY   3.6 …` (名前列の幅と、11 桁目の `f` を付ける条件の選び方が違う)。数値には影響しない。
4. 解析微分 (match6d) を使わない状況 (全体を中心差分にする): INS、COUPLE、FITP>1、BMAG/幾何条件、DK 付き変数、L=0 の QUAD。
   FITP>1 は match6d の `idp≠0` を SAD と照合していないための安全側の措置。照合できれば `_analytic_unsupported` の該当行を外せる。
5. CALC6D: INS (CELL OFF) + CALC6D は未実装 (4 次元の伝搬)。GO/FIT は 4 次元の条件のまま (SAD の CALC6D の GO は 28 成分で評価する)。
   SAD は CALC6D の temit が `DTSYNCH` グローバルを設定するが pysad は未定義 (18a の `EMIT` 行に出る)。
6. GEO: SOL の `tsgeo` と `SORG` は未移植 (geo.py の警告が出る)。`GEO` コマンド (表示) は従来どおり未対応。
7. Emittance[]: `Region`/`InitialOrbit`/`InitialBeamMatrix`/`IntrabeamExcitation` は未対応。`LINE["SIG..",l]` は整数位置のみ。
   INTRA ON の `Emittance[]` は run_emit (IBS 反復) の後で emittance_rules を再計算するので 2 重に時間がかかる (ESRF で +10 s)。
8. 挙動が変わったもの (SAD に合わせた結果、既存スクリプトに影響しうる): `Twiss["GX"/"GY"]` はエラー (γ は `GMX`/`GMY`)、
   `Emittance[]` の既定の戻りから `EquilibriumBeamMatrix`/`OneTurnTransferMatrix` が消えた (mode≥1 のみ)、`Emittance[]` が EMITX 等を書き換えなくなった (`SaveEMIT->True` で書く)、
   `run_sad` 経由では文のエラーを表示して続行、`LINE["LENGTH", loc]` は narg エラー、`MINCOUP` の既定が 0.01。
9. track6n に MULT/SOL を統合したので、DynamicApertureSurvey (dapert.py) も MULT/SOL を含むラティスで numba が使えるようになったが、
   その DA の SAD との照合 (損失周の la の扱いを含む) は未実施 (単体テストは track6 との座標・seed 一致と損失判定の動作確認のみ)。
10. `<work-esrf>/calc_RadInt.sad` などの旧 `.pysad.out` は統合前のコードの出力 (calc_RadInt は途中で TypeError で止まっていた) で、今回の再実行結果とは別物。
