# DESCEND (最急降下) の分岐の調査 (2026-10-16)

最終更新: 2026-11-16

担当: pysad/match6.py (match6d.py, tsvd.py は読んだだけで変更なし)。課題は docs/history/work_item4.md §E の
「FIT の 2 段目で pysad は SAD と 112 反復目まで一致するが、113 反復目の DESCEND から分岐する」。

## 結論

**丸めの増幅ではなく、論理の違いだった。** pysad の `Matcher.go()` は tffsmatch.f の `wcal=chgmod` を実装しておらず、
`wcal` を最初の 1 回だけ真にしていた。SAD は **モード切替 (NEWTON↔DESCEND) の直後の評価ごとに** tffscalc の重み (twfit の `wfit`、`etamax`、`avebeta`) を
現在の光学から計算し直す。pysad は初期値の重みを使い続けたので、モード切替の後は残差・微分の行 (EX/EPX など光学に依存する重みの条件) が SAD と違った。
NEWTON だけで走る間は重みが初期のままで SAD と同じなので 112 反復目まで一致し、勾配 (tgrad) は重みの違いにそのまま効くので DESCEND の 1 歩目から 16% ずれた。
(docs/06_matching.md §4 の擬似コードにはすでに `wcal=chgmod` と書いてあったが、コードに入っていなかった。)

修正 (pysad/match6.py `go()`): do200 の先頭で `wcal = chgmod` とし、calc の後の `wcal = False` を削除。1 行の変更。

## 手順と根拠

1. 再現: `PYSAD_MATCH_DEBUG=1` で pysad を実行し、SAD の出力と並べて 113 反復目から違うことを確認。
2. SAD のデバッグ版: WSL 内の SAD ソースの複製で、rng 調査用に入れられていた 5 ファイル photons.f / random_plugin_sad.c / tbend.f / tquad.f / tturn.f を
   リポジトリ版に戻し、tffsmatch.f と tffscalc.f に write 文を追加。ビルドは `make exe` で 1〜2 分。
   リポジトリの SAD は変更していない。デバッグ版を入れる前の出力が `.sad.out` と完全一致することを確認済み)。
   印字した量: 各反復の r, rbest, fact, crate, aimprv, smallf, badcnv, alate, amedcv, newton, chgmod (DBGA/DBGB)、dg, dval, wvar, bestval, wlimit (DBGC/D/W/V/L)、
   反復 105〜112 の df と qu 全成分 (DBGN/F/Q)、wcal のときの emx, emy, wfit, 適合点の βx, αx (DBGWP/F/B/A)。
   pysad 側は `PYSAD_MATCH_DEBUG=1` で同じ形式を印字 (`=2` で df と qu も。match6.py の `_DEBUG` ブロック)。
3. 論理の他の部分の直訳精度 (tffsmatch.f と match6.py を行ごとに読んだ): `fuzz`/`badcnv`/`alate`/`nretry` の論理、`fmincube` (tffindroot.f)、`tgrad`、`dg` と切替 (dg>0 → DESCEND)、
   `wlimit/nvara` の縮小、`chgmod` ブロックの `fact=min(1,2*fact)`、直線探索の `fact/2` と `a=fact/f1` の補間、`tsvd` の呼び出し — すべて SAD と同じで、他に違いは見つからなかった。
   SAD に 1e-15 の摂動を入れる実験は、原因が論理の違いと特定できたので行っていない。

## 回帰 (修正前の出力を保存して比較)

- bench_level13a / 13b / 14a の `run_sad.py` 出力: 修正前後で完全一致 (diff なし)。
- pytest (`--ignore` 3 本): 146 pass。`QUICK=1 python run_regression.py`: 35 ベンチ FAIL 0 (196 s)。

## 追記 (10/16)

- 不安定ラインの続行と反復表の違いは docs/history/work_unstable.md で直した。
  原因は不安定判定だけでなく、fmincube の上限 `m_euler` がオイラー定数 γ (0.5772) なのに pysad がネイピア数 e を入れていたこと (+ nstab の数え方・qcod overflow の扱い)。

## 道具 (再利用用)

- SAD のデバッグ版は tffsmatch.f / tffscalc.f へ write 文を入れて作る。gs はパスに空白があると動かないので WSL のホーム下で実行する。
- 反復ごとの r の比較、df と qu の比較は小さな比較スクリプトで行った (リポジトリには含まない)。
- pysad: `PYSAD_MATCH_DEBUG=1` で各評価の r/rbest/fact/変数と DBGC (dg, dval, wvar, bestval, wlimit)、`=2` で df と qu 全成分。

## 2026-11-16 の確認・追記

- 2026-11-16 確認: 変化なし。根拠: この話題の記述は work_unimpl.md と docs/22〜35 に出てこない (grep で確認)、関係するモジュールに 11/06 以後この話題の変更の記録がない。この文書は 11/06 時点の記録のまま。
