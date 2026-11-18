# FFS の小さなコマンド (`pysad/ffsextra.py`)

最終更新: 2026-11-18 (v2.1.0 の後) (現状のスナップショット。経緯は `history/` にある)

K64 の tfvars.f・tfattr.f・tfmat.f・tshow.f・tffsa.f の直訳。表の書式は autofg (`'12.9'` など) まで gs の出力と同じ (`tests/test_ffsextra.py`)。

| コマンド | 内容 |
|---|---|
| `VARIABLES` (`VARS`、`VAR`) | FREE した変数の表 (Now、Previous、Saved、Minimum、Maximum、結合先、係数)。`FFS["VAR"]` は同じ値のリスト |
| `ATTRIBUTE` (`ATTR`) | 要素の既定キーの値と限界値 |
| `MATRIX [SYMPLECTIC\|PHYSICAL] [from to]` | 転送行列 (CALC4D は 4×5、CALC6D は 6×6) |
| `SHOW` | 現在のマッチング条件 (`FFS["SHOW"]` はリスト) |
| `VARY keyword pattern` | 要素の既定キーを変える |
| `RECOVER` | FREE した要素の値を、直前の GO の前の値と交換 |
| `REVERSE` | 入口の AX AY EPX EPY R2 R3 DPX DPY の符号を反転 |
| `PRINT expr`、`EXECUTE "文字列"` | 式の値を `Out[n]:= 値` で出す、文字列を FFS コマンドとして実行 |
| `OUTPUT(OUT) file` / `APPEND(APP) file` / `CLOSE OUTPUT` | 出力 (FFS の表示と `Print`) をファイルへ切り替え・追記・戻す。ファイル名は引用符つき、またはファイル番号 |
| `INPUT(IN) file` / `READ file` | ファイルの文を続けて実行 (末尾で元の入力に戻る) |
| `STATUS` (`STAT`) | フラグ・特別変数の一覧 (twbuf と同じ 8 桁の桁揃え)。CTIME・DT は実行時間なので gs と一致しない |
| `CLOSE`・`TERMINATE`・`SUSPEND`・`RESUME` | 入力側は何もしない (入力はファイルを丸ごと実行するため) |

## gs と比べた確認 (2026-11-16 夜)
- `MATRIX` の名前だけの引数 (例 `MATRIX 要素名`): gs は `NAME.1` として扱い、ラティスに無い名前は読み飛ばす (to=$$$)。pysad も名前だけを位置にし、無い名前は読み飛ばす (`test_ffsextra.py`)。
- `DUMP`: 誤差を入れたリングで gs と比べ、値が一致する。DK・DDK 列は Fortran の `1p g12.5` (0 は `0.0000`、小さい値は `4.22829E-04`) で、`ffscmd._gwd` が同じ書式にする (gs と一致。`test_ffsextra.py` の `test_dump_format_matches_gs`)。
- FIT のフラグ UNIPOL・RELW: gs は UNIPOL で変数の符号を変えず (tffsmatch.f:669)、RELW で重みが変わり、FIT の結果が変わる。pysad は `ffscmd.matcher_` で match6 のフラグとつなぎ、既定・UNIPOL・RELW の 3 通りが gs と全桁一致 (`test_item4.py` に 3 件)。
- ソースの確認 (gs での動作は試していない): MON・MCOD・TCOD のコマンドは K64 の tffsa.f に無い (FFS では式として評価されて何も起きない。v2 モードは「K64 に無いコマンド」と警告して文全体を読み飛ばす。機能は実装しない)。
  フラグの効果: GAUSS は MAIN の TRACK・MEASURE の初期分布、CTIME は CPU 時間の表示、CONV は収束の判定・表示、EMIOUT は EMIT 報告の追加行 (pysad は spin の行だけ対応)、BIPOL は tffsmatch.f だけ。

- `DUMP pat`: gs の語の変数は `character*20` で、パターンは 20 文字で切れる。合う要素がなければ全要素。DDK 列は K0 のキーワードのある要素 (BEND・MULT・TCAVI) だけ (QUAD にはない)。
- `DISP` のモード R・O・P・GA・Z・A・B (BEAM、EMIT 後)・D・G・OG と REFERENCE・DREFERENCE (見出しの文字の置換は tfdisp.f の trefheader/tdrefheader のとおり)、`REFERENCE` コマンド、`Twiss["R…"/"D…"/"RALL"/"DALL"]` は gs と一致 (`test_ffs_gs_match.py`。詳細は docs/03 §8)。
- `RENUMBER comp` と `LINE["MULT"]` (tfinimult の直訳。同名 1 台だけは 0)、`USE [NO]EXPAND`・`CALC|GO [NO]EXPAND` は gs と一致 (docs/03)。`NOKEEPEXP` (tffsa.f の `calexp = .not. tfvcomp()`) と RENUMBER の INDEX 属性 (2 以上なら同名の要素の番号をその値から振る。`LINE["INDEX",…]` も読める) は実装した (両モード。gs の実測と一致)。`tfvcomp` は、既定でないキーの成分の変数 (成分ごとの FREE) があれば真 (v2.1.0 から。docs/03 §3.3)。
- `QF.3 K1 0.9;` のように `名前.番号` で始まる文は、gs では式として評価されて何も設定しない。pysad も何も設定しない。出力の `Out[n]:=` 行は再現しない。

## 限界・K64 との差
- `PRINT` の `Out[n]` は、USE が `Out[1]:= USE` を出したものとして 2 から数える (pysad は USE でこの行を出さない)。
- `MATRIX` の 5 列目 (分散の項) は η(l2) − M·η(l1)。要素が 1e-16 程度の雑音を持つ箇所は gs が `-4.798E-16` などと出すのに対し、pysad は 0 を出す。
- `MEASURE`・`TRCOD`・`LOGMEASURE`・`ORBITJITTER`・`PUTMEASURE`・`GETMEASURE`: **v2 モードだけ**実装した (`pysad/measure.py`。使い方は docs/03 の「v2 モードの MEASURE・TRCOD」)。k64 モードは未対応 (警告して無視)。`MEASURE` は gs が何も出力せず (B42、`TRCOD` は SIGSEGV)、**gs との数値比較は無い** (K64 のソースの読みと解析式・pysad の独立な経路で確かめた)。`SPLIT` は K64 でも obsolete。`SEED` も obsolete。
- v2 モードの台本の読み方 (MAIN の `OFF RFSW` が FFS に届く、`ORG` の負の値、MOMENTUM を設定しない計算の警告、`FFS;` の書き方、MAIN の代入 `NAME = 式;`) は docs/03 の同じ節に書いた。
- gs は `TERMINATE INPUT` の後で入力が尽きて止まるが、pysad は続ける。
