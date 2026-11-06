# 作業ログ: SADScript / FFS / parser の実スクリプト互換 — 現状記録 (2026-10-15)

担当範囲: `pysad/sadscript.py`, `ffscmd.py`, `parser.py`, `run_sad.py`。前担当 2 名が中断した作業の「どこまで入っているか」を記録する。
コードは変更していない。方法は、`pysad_20261013.tar.gz` (最終更新 10/14 23:50 時点) を scratchpad の `old_pysad/` に展開して現在の `pysad/` と `diff -u` し、
動作確認を小さく行った。回帰 (bench_level13a/13b、pytest 49 件) は確認済みとのことなので再実行していない。

## 1. diff の規模

parser.py は 1640 → 1723 行、sadscript.py は 2404 → 2430 行、ffscmd.py は 1680 → 1756 行、run_sad.py は変更なし
(行末が CRLF/LF で入れ替わっているため、素の `diff -u` は 3000〜4800 行になる。`--strip-trailing-cr` を付けると実質差分は parser 172、sadscript 66、ffscmd 107 行)。
そのほか `lattice.py` (43 行、Twiss 書き戻し)、`track6n.py` (87 行、数値関数と la 計数) も変わっているが、前者は work_match_fixes.md、後者は他担当の範囲。

## 2. 依頼項目ごとの判定

| # | 依頼項目 | 判定 | 根拠 (diff / 動作確認) |
|---|---|---|---|
| 1 | 同型再定義は書いたキーワードだけ上書き (doelem.f) | 実装済み | `SADEnvironment.define_element` が同じ TYPE なら `old.update(params)`。`BEND B=(ANGLE=0.2)` を再定義しても L=1, E1=0.5, K0 が残るのを確認 |
| 2 | DECA/DODECA の無視 | 実装済み | キーワード表に DECA/DODECA ほか TCAV, UND, WIG, APERT, ST, COORD, MAP, BEAMBEAM, PROF, INS, SPCH, TEST を追加。`_build_element0` に Dodecapole と Solenoid を追加。`DECA D5=(…K4=1)`、`DODECA D6=(…K5=1)` が env に登録されるところまで確認 (build_lattice での生成は確認していない) |
| 3 | 単位 MV/MHZ/RAD/M/DEG の未登録 | 実装済み | `_SAD_UNITS` を initb1.f の表に置換 (RAD, M, T, GAUSS, JOULE, V/KV/MV/GV, HZ/KHZ/MHZ/GHZ を追加)。`VOLT=1MV` → 1.0e6、`FREQ=499.6 MHz` → 4.996e8 を確認 |
| 4 | -N*LINE の反転 | 実装済み | `_parse_beamline_refs_sad` を `[-] n * [-] name/(…)` に書き直し、`_expand_refs` で入れ子リストも反転。`LINE L=(Q -2*(B Q) B)` の展開が Q, Q, B, Q, B, B になる (期待どおり) |
| 5 | `_raw` にしか残らない属性 | 実装済み | `_apply_extra_params` を新設。QUAD/SEXT の ROTATE→ROTATION、DX/DY、BEND の K0/FRINGE/F1/F2/AE1/AE2、CAVI の PHI/DPHI/DVOLT、MULT の Kn/SKn、SOL の BZ を属性に反映。ROTATION=0.1, DX=1e-4, K0=0.01, FRINGE=1, PHI=0.2 を確認 |
| 5b | BEND E1/E2 の解釈 (推測ベースだった) | 実装済み (依頼外) | 常に ANGLE 比 (`E1*angle + AE1`)。tbende.f に合わせた。E1=0.5, ANGLE=0.2 → 0.1 を確認 |
| 6 | `ffscmd._TYPECODE` の MULT/SOL | 実装済み | Solenoid=20, Multipole=22 に訂正 (tfstk.f)。確認済み |
| 7 | RFSW の既定 (track.f:146-148) | 実装済み | `_FLAGS_DEFAULT['RFSW']` を True に。確認済み |
| 8 | 単一引用符文字列 `Twiss['*',1]` | 実装済み | `tokenize` と `tokenize_full` の両方で `'` を区切りとして扱う。`x=Twiss['*',1]` と `Twiss["*",2]` がどちらも str トークンになる |
| 9 | `Twiss["*"]`, `Twiss["*",k]` | 実装済み | `_TWISS_ALL` 28 個 (AX…ZPY) を返す。値の検証は未実施 |
| 10 | Twiss の端数位置を qfraccomp に合わせる | 実装済み (数値は未検証) | `frac_element` が qfraccomp の直訳 (BEND は F1/FRMD、E2=0、QUAD/CAVI/MULT の FB1/FB2 を処理)。`row_frac` は入口の閉軌道から転送行列で伝搬。逆向き要素の E1↔E2 入れ替えは未対応と注記されている |
| 11 | chrom の `?No variable.` (FREE SD*,SF*) | 現行では再現せず (原因は不明) | `_freefix` 自体は無変更。pf2 chrom を現行で走らせると変数が登録され収束する (work_match_fixes.md 参照) |
| 12 | esrf/pf2 実スクリプトを広く走らせて直す | 途中 | 10/15 08:50 の一括実行 (esrf 40 本) は完走 17 本、例外 12 本、打ち切り 11 本 (計 40 本)。例外のうち `unexpected character "'"` の 4 本は項目 8 で解消。USE の未定義ビームラインは SAD と同じ `Missing beamline in USE.` を出すようにした (disp_eta は TEST2 が esrf ラティスに無いので SAD でも同じ)。一括実行の再走は未 |

依頼外で入っていたもの: 複素数の SAD 表示 `(re+im I)`、`$FORM` 下の NaN 表示、`.5` のように先頭 0 を省く `autos`、Infinity の綴り、単項マイナスの優先度 (`-2^2=-4`)、
Sqrt 等の複素数対応 (cmath)、MAIN 部の `TYPE`/`PRINT` 読み飛ばし。

## 3. 動作確認

- **parser / tokenizer**: 上の表のとおり (項目 1〜8 は `python -c` 相当で確認)。
- **<work-esrf>/disp** (`cd <work-esrf>; python -u ../run_sad.py ../pysad/tests/kek-ls_v03_12.lat disp`): 約 15 s、例外なし。SAD 出力 `disp.sad.out` の Twiss 表 780 行を `#` 列で対応づけ、
  AX…DetR の数値を比べると表示精度の範囲で一致 (最大相対差 1e-12)。`$$$` は NX=47.5800, NY=17.6200。tune も SAD と一致 (-0.4199589, -0.3799494, -0.0040753)。
- **差が残る点**: (1) SAD は反転配置された要素名に `-` を付ける (`-BPI7.2` 等 380 要素、pysad は付けない)。(2) `Equilibrium position dz` は SAD 20.618669 mm、pysad 20.619466 mm、
  `Bunch Length` は 3.77333112 と 3.77333234 mm (2e-5 の差)。(3) SAD は `Beam tilt` の行を出す (pysad は出さない)。(4) Emittance Y の桁 (SAD -4.7e-38、pysad -4.4e-21) は 0 同然で意味がない。

## 4. 未着手・残り

`-` 付き要素名の表示、`Beam tilt` 行、DRAW/GEO/MATRIX 等の FFS コマンド (未対応のまま警告)、`ON INTRA` の FFS 内使用、`Repeat/Until`、`CalculateOptics`、`Twiss["SIGX"]`、
逆向き要素の qfraccomp、qfraccomp の数値検証 (端数位置の Twiss を SAD と並べる小ベンチ)。bench と pytest はこの記録の範囲では再実行していない。

> 追記 (10/16、項目 6 = docs/history/work_item6.md): §3 の「差が残る点」(1) `-` 付き要素名、(2) Equilibrium dz / Bunch Length、(3) Beam tilt の行は解決 (disp の 780 行すべて SAD と一致、dz 20.618669、Bunch Length 3.77333112)。(4) Emittance Y の桁は 0 同然で未対応のまま。
