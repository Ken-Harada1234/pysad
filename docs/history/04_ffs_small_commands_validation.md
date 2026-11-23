# 04_ffs_small_commands の検証の記録

旧 35_ffs_small_commands.md の「gs と比べた確認 (2026-11-16 夜)」の節と、旧 31_delcor.md の「検証 (gs との比較)」の節の元の記述。

## 旧 35 gs と比べた確認

## gs と比べた確認 (2026-11-16 夜)
- `MATRIX` の名前だけの引数 (例 `MATRIX 要素名`): gs は `NAME.1` として扱い、ラティスに無い名前は読み飛ばす (to=$$$)。pysad も名前だけを位置にし、無い名前は読み飛ばす ([`test_ffsextra.py`](../pysad/tests/test_ffsextra.py))。
- `DUMP`: 誤差を入れたリングで gs と比べ、値が一致する。DK・DDK 列は Fortran の `1p g12.5` (0 は `0.0000`、小さい値は `4.22829E-04`) で、`ffscmd._gwd` が同じ書式にする (gs と一致。[`test_ffsextra.py`](../pysad/tests/test_ffsextra.py) の `test_dump_format_matches_gs`)。
- FIT のフラグ UNIPOL・RELW: gs は UNIPOL で変数の符号を変えず (tffsmatch.f:669)、RELW で重みが変わり、FIT の結果が変わる。pysad は `ffscmd.matcher_` で match6 のフラグとつなぎ、既定・UNIPOL・RELW の 3 通りが gs と全桁一致 ([`test_item4.py`](../pysad/tests/test_item4.py) に 3 件)。
- ソースの確認 (gs での動作は試していない): MON・MCOD・TCOD のコマンドは K64 の tffsa.f に無い (FFS では式として評価されて何も起きない。v2 モードは「K64 に無いコマンド」と警告して文全体を読み飛ばす。機能は実装しない)。
  フラグの効果: GAUSS は MAIN の TRACK・MEASURE の初期分布、CTIME は CPU 時間の表示、CONV は収束の判定・表示、EMIOUT は EMIT 報告の追加行 (pysad は spin の行だけ対応)、BIPOL は tffsmatch.f だけ。

- `DUMP pat`: gs の語の変数は `character*20` で、パターンは 20 文字で切れる。合う要素がなければ全要素。DDK 列は K0 のキーワードのある要素 (BEND・MULT・TCAVI) だけ (QUAD にはない)。
- `DISP` のモード R・O・P・GA・Z・A・B (BEAM、EMIT 後)・D・G・OG と REFERENCE・DREFERENCE (見出しの文字の置換は tfdisp.f の trefheader/tdrefheader のとおり)、`REFERENCE` コマンド、`Twiss["R…"/"D…"/"RALL"/"DALL"]` は gs と一致 ([`test_ffs_gs_match.py`](../pysad/tests/test_ffs_gs_match.py)。詳細は docs/03 §8)。
- `RENUMBER comp` と `LINE["MULT"]` (tfinimult の直訳。同名 1 台だけは 0)、`USE [NO]EXPAND`・`CALC|GO [NO]EXPAND` は gs と一致 (docs/03)。`NOKEEPEXP` (tffsa.f の `calexp = .not. tfvcomp()`) と RENUMBER の INDEX 属性 (2 以上なら同名の要素の番号をその値から振る。`LINE["INDEX",…]` も読める) は実装した (両モード。gs の実測と一致)。`tfvcomp` は、既定でないキーの成分の変数 (成分ごとの FREE) があれば真 (docs/03 §3.3)。
- `QF.3 K1 0.9;` のように `名前.番号` で始まる文は、gs では式として評価されて何も設定しない。pysad も何も設定しない。出力の `Out[n]:=` 行は再現しない。


## 旧 31 検証 (gs との比較)

## 検証 (gs との比較、`pysad/tests/test_delcor.py` 3 本。2 本は CELL と INS の組み合わせ、1 本は CAVI)

- CELL の 12 セル: 相関長 2 m (QF QD)、相関長 0 (B)、ADD とワイルドカード (Q*) の 3 回の DELCOR のあとのチューン 1e-14、閉軌道 1e-10。
- INS (4 セル、SEXT・BEND を含む、`*` で全要素): チューンと軌道が一致 (乱数の列が gs と同じ)。

