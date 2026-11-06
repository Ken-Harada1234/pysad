# 作業ログ: CodCorr (esrf COD 補正) と誤差の実現 (2026-10-16)

担当範囲: pysad/sadscript.py, ffscmd.py, lattice.py。`ユーザーのディレクトリ` では何も実行していない (<work-esrf> の写しと cc_test/ を使用)。
SAD 実行は `cc_test/runsad.sh` (WSL gs、NPARA=1、標準入力 /dev/null)。参照ベンチは bench_level31a / 31b (run_regression に追加)。

## 結論

| 項目 | 結果 |
|---|---|
| A. CodCorr の Dot 不一致 | 原因は項目 C。`tunex = Twiss["NX","***"]/2/Pi` が 1 要素リストになり、応答行列 `Responsex[j,i]` の各要素がリストになって Dot の長さが合わなかった。C を直して解消。続けて Eigensystem も直した (下記) |
| B. 誤差の実現 | 原因は terror.f の master グループ (同じ要素の分割片は 1 回の乱数を共有) と tfadjst.f (隣り合う磁石の鎖は先頭磁石の DX/DY を共有) が未実装だったこと。実装して bench_level31a (128 タグ) が 1e-10 (絶対) / 1e-7 (相対) 以内で一致 |
| C. Twiss/LINE のリスト形 | 一致数が 1 ならパターンの形に関わらずスカラー、0 個は {}、複数はリスト (tftwiss.f tflinestk)。bench_level31b (35 タグ) 一致 |
| D. FFS エラーで Do 中断 | 方針変更により打ち切り (試作は戻した。下記「未対応」) |

## A. 原因の特定

1. 最小再現: `CodCorr` 単体 (小さな行列) は pysad でも動く。実データでは `CalculateCODResponse` の `tunex` がリスト → `Cos[tunex*Pi-...]` がリスト → 応答行列がリスト入り → Dot で equalleng。
2. C を直すと CodCorr に入り、次に `Eigensystem` が numpy.linalg.eig の複素数を返して `LINE["K0",..] = ...` で float 変換エラー (対称行列の縮退で虚部の雑音)。
3. (最終版は下の「k23 比較」参照: n>8 の対称行列は eigh) SAD の Eigensystem は teigen.f (K. Oide) で、固有ベクトルは単位ベクトルではない (例 {{2,1},{1,3}} → {1,-.618}, {.447,.7236})。pysad には teigen.py (直訳) が既にあるので、実行列はそれを使うようにした (sadscript.py `_eigensystem`)。固有値の並び・規格化とも SAD と一致する (n=3〜5 のランダム対称行列で全桁一致)。複素固有値の対は 行 k = 実部+i虚部、行 k+1 = 実部-i虚部 (虚部 0 の成分は実数)。複素数の表示も SAD に合わせた (`(-I)`、`(1+I)`、`(1-2 I)` の符号は虚部の書式の中)。
   - 注意: teigen の 2x2 ブロック (実固有値の対) の並びは丸めに敏感で、三重対角などの特殊な対称行列 (整数行列など) では SAD (-Ofast) と並び・規格化が違うことがある。密な応答行列ではこの枝に入らない (40x50 の乱数行列で CodCorr の結果が SAD と 2e-10 で一致)。
4. SAD の qcod は一周行列の跡が NaN になると (`qtwiss1` の over) 後退せずに "Overflow & closed orbit not found" で終わり、twiss は前回の値のまま残る。pysad は NaN の後も Newton を続けて収束してしまい、SAD が失敗する誤差実現 (seeds[1] と seeds[4]) で値が違っていた。`_qcod_periodic` で NaN の跡を検出して打ち切り、calc() は直前の行 (`_last_rows`) を残すようにした。これで k1 テスト (esrf v05_61、seeds[1..4]) の 4 実現すべてが SAD と一致 (失敗する 1, 4 も含む)。
   - SAD で correct_cod の最初の実現 (k=1) の COD が 1e-17 だったのはこれ (誤差が効かないのではなく qcod が失敗して前回の twiss が残る)。

## B. terror.f の実現

- master (tmast.f): 磁石 (BEND/QUAD/SEXT/OCT/DECA/DODECA/MULT/CAVI/TCAV) の先頭は、同じ要素 (名前が同じ。間にゼロ長の別要素のみ許す) の最後の位置。2 番目以降は -先頭。DELX 等はグループの先頭ごとに乱数を 1 回だけ引き、グループ内の master≠0 の要素 (マーカーや L=0 の別要素は除く) に同じ値を入れる。pysad は全要素で個別に引いていたので乱数列がずれた。`_master_array`。
- tfadjst.f: 誤差コマンドの後 (UNIFORM 以外)、隙間なく並んだ磁石 (BEND/QUAD/SEXT/OCT/MULT/CAVI、間にゼロ長要素のみ) の鎖は先頭磁石の DX/DY のずれを共有する (例: `BPI3H HD BPI3H` の HD は BPI3H と同じ)。`_tfadjst`。
- LINE["K1"/"ANGLE",..] の読み出しは DK の誤差込みの実際の値 (SAD tflinek)。Element[] は設計値のまま。DK で効くのは主強さ (QUAD→K1、BEND→ANGLE) だけで BEND の K1 には掛けない。
- `LINE[key,"NAME.n"]=` で触ったキーは SAD の tftouch と同じく USE まで残り、CALC のたびに代表 (先頭の 1 台) の値に揃える (`_spread_pending` を使い切りにしない)。後から DELX 等で入れた非代表の誤差は次の CALC で代表の値に戻る。
- SAD の挙動メモ: 誤差を大きく入れると qcod が失敗する (ESRF で DELX 1e-4 Q* など)。bench_level31a は DELX 2e-5 などにして qcod が収束する範囲にした。

## 検証

- bench_level31a: DELX/DELY/DK/DTHETA/DDK、UNIFORM/COHERENT/ADD、グループ・鎖 (BPI3H/HD、SDA 4 連)、`LINE[..]=` の展開、`DELX 0` による戻し。Twiss の COD/ETA/NX/NY (9 位置)、LINE の DX/DY/K1/ROTATE/DROT/K2/ANGLE/K0 を 128 タグ比較。
- bench_level31b: `Twiss["NX","***"]`、`$$$`、`^^^`、名前、`NAME.*`、`Q*`、`RFC*` (一致 1 個)、一致 0 個、数値位置、`NAME+n`、リスト引数、`Twiss["*",..]`、LINE の同様の形。
- k1 テスト (cc_test/k1.sad: esrf v05_61 + correct_cod と同じ誤差セット + seeds[1..4]): COD/β/η/チューンが SAD と 1e-10 以内 (k=2,3 は実際の COD、k=1,4 は SAD の qcod 失敗 = 前回の twiss)。
- 回帰: run_regression QUICK=1 で 50 ベンチ FAIL 0 (31a/31b 追加)、pytest 164 passed (--ignore 3 本)。13a/13b/14a/16a/22a/23a/28a の出力は作業前後で diff (codcorr_base/before, after)。

## 4 本 (correct_cod, nomon, nomon_dq, print_emittance) の現状

- Dot で止まらなくなり、CodCorr に入って先まで走る。pysad は 1 実現 40 反復に 900 s 以上かかる (SAD は 110 実現で 270 s) ので、全 110 実現は完走しない (性能の項目。打ち切り 600 s)。実現数と反復数を絞った版 (correct_cod_nomon_k23.sad: k=2,3、反復 3) で SAD と比較した結果は下の「k23 比較」。

## 未対応・記録のみ

- D (FFS のエラーで外側の Do を中断): SAD は termes のエラー (すべて ' ???-FFS-Error-' で出力。出力ユニット 0 でも表示) で iffserr=1 → FFS[] の残りのコマンドを読まず、`FFS::error:  Error in FFS command 1 in FFS["..."]` で外側の文を中断する。FFS 文字列の中の REPEAT/UNTIL も SAD は動く (pysad は最上位のみ)。`DELX 1e-5 QXX*` (一致なし) や `QXX* K1 0.5` は SAD では無言 (pysad は '?Undefined' を出す)。試作は方針変更で戻した。
- SAD の teigen との並びの差 (三重対角など特殊な対称行列)。
- 性能 (4 本の全実現を走らせる場合)。

## k23 比較 (esrf v05_61、nomon、k=2,3、補正反復 3)

SAD は SAD 出力と同じ入力 (cc_test/kk.sad.out)、pysad は cc_test/kk.pysad。
- 補正前の COD/ETA の最大位置・最大値・RMS (k=2): 位置 (21, 83, 20, 163) と値 (.0033183117, .0023829386, .1429398595, .0768729489 ほか) が **SAD と一致** (誤差の実現は B で解決)。
- 補正の 1 回目 (CodCorr): SAD の `mat x/y rms` = 7.2e-6 / 1.64e-5。pysad は teigen 直訳で 3.1e-5 / 7.3e-5、numpy の厳密な対称固有値分解で 6.7e-6 / 3.9e-6。**一致しない**。応答行列 (160x160, 120x120 の tres.res) は縮退した固有値が多く、SAD の teigen (-Ofast) の固有ベクトルは丸めで決まる非直交の組で、CodCorr の結果 (非直交ベクトルを `V^T Λ^-1 V` に使う) が丸めに敏感になるため、pysad (どの実装でも) では再現できない。最終版は n>8 の実対称行列で厳密な固有値分解 (物理的に正しい擬似逆行列) を使う (x は SAD と 7%、y は 4 倍違う)。
- 補正後の CALC: SAD は "Closed orbit not found" (10 反復で未収束、Twiss は更新され「Unstable」)、pysad は一周行列の跡が NaN になり "Overflow" (前回の twiss を保持)。補正量が SAD と違うので比較対象にならない。
- 結論: A は「Dot で止まる」を解消し CodCorr 以降まで走るようになったが、補正の数値は SAD の teigen の丸めに依存するため SAD と一致しない (誤差の実現・補正前の COD は一致)。また 1 実現 40 反復に 900 s 超 (SAD は全体で 270 s)。
