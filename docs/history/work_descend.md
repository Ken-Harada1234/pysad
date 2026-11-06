# DESCEND (最急降下) の分岐の調査 (2026-10-16)

担当: pysad/match6.py (match6d.py, tsvd.py は読んだだけで変更なし)。課題は docs/history/work_item4.md §E の
「tune_sv_v30_24_1 の 2 段目で pysad は SAD と 112 反復目まで一致するが、113 反復目の DESCEND (1.9297E-02 と 1.9683E-02) から分岐する」。

## 結論

**丸めの増幅ではなく、論理の違いだった。** pysad の `Matcher.go()` は tffsmatch.f の `wcal=chgmod` を実装しておらず、
`wcal` を最初の 1 回だけ真にしていた。SAD は **モード切替 (NEWTON↔DESCEND) の直後の評価ごとに** tffscalc の重み (twfit の `wfit`、`etamax`、`avebeta`) を
現在の光学から計算し直す。pysad は初期値の重みを使い続けたので、モード切替の後は残差・微分の行 (EX/EPX など光学に依存する重みの条件) が SAD と違った。
NEWTON だけで走る間は重みが初期のままで SAD と同じなので 112 反復目まで一致し、勾配 (tgrad) は重みの違いにそのまま効くので DESCEND の 1 歩目から 16% ずれた。
(docs/06_matching.md §4 の擬似コードにはすでに `wcal=chgmod` と書いてあったが、コードに入っていなかった。)

修正 (pysad/match6.py `go()`): do200 の先頭で `wcal = chgmod` とし、calc の後の `wcal = False` を削除。1 行の変更。

修正後は tune_sv_v30_24_1 の 2 段目の反復表 (NEWTON 2〜108、DESCEND 113 の 1.9297E-02 / 0.2500、NEWTON 116 の 1.8809E-02 / 0.3700、最終 Residual = 1.8569E-02) が
**SAD と全行一致** (反復の表示行はすべて同じ。差は行末の空白のみ)。

## 手順と根拠

1. 再現: `PYSAD_MATCH_DEBUG=1 python run_sad.py <work-pf2>/tune_sv_v30_24_1.sad_r_full.sad`。SAD 出力 (`tune_sv_v30_24_1.sad_r.sad.out`) と並べて 113 反復目 (SAD 1.9297E-02 / pysad 1.9683E-02) から違うことを確認。
2. SAD のデバッグ版: (WSL のホームの作業フォルダ) (WSL。(WSL のホームの作業フォルダ) を複製して rng 調査用に入れられていた 5 ファイル photons.f / random_plugin_sad.c / tbend.f / tquad.f / tturn.f を
   リポジトリ版に戻し、tffsmatch.f と tffscalc.f に write 文を追加。差分は `<work-scratch>/patch_tffsmatch.py` `patch2.py` `patch3.py`、ビルドは `make exe` で 1〜2 分。
   リポジトリの SAD は変更していない。デバッグ版を入れる前の出力が `.sad.out` と完全一致することを確認済み)。
   印字した量: 各反復の r, rbest, fact, crate, aimprv, smallf, badcnv, alate, amedcv, newton, chgmod (DBGA/DBGB)、dg, dval, wvar, bestval, wlimit (DBGC/D/W/V/L)、
   反復 105〜112 の df と qu 全成分 (DBGN/F/Q)、wcal のときの emx, emy, wfit, 適合点の βx, αx (DBGWP/F/B/A)。
   pysad 側は `PYSAD_MATCH_DEBUG=1` で同じ形式を印字 (`=2` で df と qu も。match6.py の `_DEBUG` ブロック)。
3. 調べた結果 (SAD と pysad を数値で突き合わせ):
   - NEWTON の 108 反復目まで: r の相対差 ≤ 5e-8 (丸め)。108 反復目の df と qu は相対 1e-7 / 1e-10 で一致。
   - 109 反復目で NEWTON→DESCEND に切り替わる (chgmod=真; 条件の `smallf/badcnv/amedcv` の値も SAD と pysad で同じ)。切替後の評価 (SAD の iter 110) の `df` と `qu` を比べると、
     条件 1, 2, 6, 7 (EX/EPX 系) の行が行ごとに一様に 7.6%, 0.23%, 9.1%, 8.2% ずれ、他の行は 1e-10。つまり重みだけが違う。
   - SAD の `wfit(1)` は反復 1 で 1.9187E+06 (βx=2.7436)、反復 110 で 2.0755E+06 (βx=2.3445)。pysad の値は 1.9187E+06 のまま (初期の重みの使い回し)。
   - `wcal` を修正すると 110 以降の DESCEND の試行点が SAD と一致: fact=1 で 8.9499710E-02 (SAD 8.9499712E-02)、fact=0.5 で 2.84654514E-02 (SAD 2.84654507E-02)、fact=0.25 で 1.92974074E-02 (SAD 1.92974074E-02)。
4. 論理の他の部分の直訳精度 (tffsmatch.f と match6.py を行ごとに読んだ): `fuzz`/`badcnv`/`alate`/`nretry` の論理、`fmincube` (tffindroot.f)、`tgrad`、`dg` と切替 (dg>0 → DESCEND)、
   `wlimit/nvara` の縮小、`chgmod` ブロックの `fact=min(1,2*fact)`、直線探索の `fact/2` と `a=fact/f1` の補間、`tsvd` の呼び出し — すべて SAD と同じで、他に違いは見つからなかった。
   SAD に 1e-15 の摂動を入れる実験は、原因が論理の違いと特定できたので行っていない。

## 回帰 (修正前の出力を保存して比較)

- bench_level13a / 13b / 14a の `run_sad.py` 出力: 修正前後で完全一致 (diff なし)。
- <work-esrf>/matching1_full.sad: 反復列・K1 とも完全一致。
- <work-pf2>/matching_v30_12_r_full.sad: NEWTON/Matched/Residual の表は同じ。最終の FIT (FREE S*、変数 2 個で SAD と別の局所解に着く既知の問題: EMITX 9.8323E-9 vs SAD 9.8358E-9) の
  SF/SD の K2 が 5e-8 (相対) 動いただけで、tune・残差は同じ。matching_v31_23_r_full.sad: 表・最終値とも同じ (警告 2 行が消えただけ)。
- pytest (`--ignore` 3 本): 146 pass。`QUICK=1 python run_regression.py`: 35 ベンチ FAIL 0 (196 s)。

## 残った点 (今回の担当外、記録のみ)

- tune_sv の続き (`USE RING; cell calc emit`): SAD は垂直方向が不安定 (Vertical Unstable) でもチューン (νy=7.0293) を返して続行し `ny not matched` を出すが、
  pysad は `periodic_optics` が `sin(μy)≈0` で ValueError、`qcod` の閉軌道が見つからず、スクリプトが止まる。qcell の不安定時の扱い (resume の「qcell の不安定時処理は未移植」) の別件。
  2 段目のフィットまでは SAD と一致するようになったので、以後の比較はこの別件待ち。
- matching02_r (NCELL、変数 1 個の解なし) の `fact` の列の違いは今回の修正では変わらない。SAD のデバッグ版で見ると、SAD は試行点の多くで `imprv=偽` (r=3.74 < rbest=4.44 でも不採用。
  `resle` は nstab を先に比べるので試行点が `nstab>0` (不安定) だったと推定、nstab は未印字)。pysad の不安定判定 (チューン境界付近の 1 周行列の扱い) が SAD と違う可能性が高く、上の qcell の件と同根かもしれない。

- **[10/16 追記] 上の 2 件は解消**: tune_sv の不安定ラインの続行と matching02_r の反復表は docs/history/work_unstable.md で SAD と全行一致にした。
  matching02_r の原因は不安定判定だけでなく、fmincube の上限 `m_euler` がオイラー定数 γ (0.5772) なのに pysad がネイピア数 e を入れていたこと (+ nstab の数え方・qcod overflow の扱い)。

## 道具 (再利用用)

- `<work-scratch>/patch_tffsmatch.py` `patch2.py` `patch3.py`: (WSL のホームの作業フォルダ) の tffsmatch.f / tffscalc.f へ write 文を入れる。`rebuild*.sh` でビルドと実行、`runsad.sh <sad> <名>` で任意の SAD スクリプトを実行 (gs はパスに空白があると動かないので WSL のホーム下 `(WSL のホームの作業フォルダ)` で実行)。
- `<work-scratch>/run/cmp.py` (反復ごとの r の比較)、`cmpq.py` (df と qu の比較)。
- pysad: `PYSAD_MATCH_DEBUG=1` で各評価の r/rbest/fact/変数と DBGC (dg, dval, wvar, bestval, wlimit)、`=2` で df と qu 全成分。
