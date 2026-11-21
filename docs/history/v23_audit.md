# v2.3.0 文書の総点検 (W2、2026-11-21)

依頼: 「解説、説明文書に抜け落ちや未更新、語句の不統一や個人名や個人ファイル名、設定がないかもう 1 回確認」。
対象: 追跡対象の文書・説明 (`*.md`、`example/**` の説明とヘッダ・コメント、`pysad/*.py` の docstring とコメント、`pysad/tests/**`)。`docs/history/` は個人情報・個人パスと事実の誤りだけ。`docs/html`・`local/` は追跡外で対象外 (リリースノートだけ用語を直した)。
挙動は変えていない。変えたコードは、`example/` のスクリプトのハードコードされたパス (環境変数に置き換え)、コメント 1 行 (`pysad/emit6d.py`)、試験の正規表現の書き方 (`test_ddba_lattices.py`) だけ。

## 1. 点検の観点

1. 個人名、メールアドレス、トークン
2. 個人パス・個人のフォルダ名・社内ホスト名・設定
3. 抜け落ち・未更新・矛盾 (版・件数・B 番号、モジュールの地図、API、フォルダの索引、リンク、古い記述)
4. 語句・表記の不統一
5. 例の `.md` の数値と出力ファイルの一致 (抜き取り)

## 2. 見つけたもの

### 2.1 個人名・メール・トークン

| 件数 | 内容 | 扱い |
|---|---|---|
| メールアドレス 0 件 | `@` を含む文字列は Python の行列積 `@` だけ | なし |
| トークン・パスワード・鍵 0 件 | | なし |
| 社内ホスト名・社内 URL 0 件 | 元のフォルダ名・社内のホスト名・社内の文書の名前は、DDBA の試験が「入っていないこと」を調べる正規表現の中だけ (文字列を分けて書いてある) | なし |
| 同僚の敬称つき 1 件 | `example/ddba/v07_14/modules/modules_errors.n` の「小林氏の逆行列法の式」 | 「逆行列法の式」に直した |
| 外部の人の敬称つき 2 件 | `docs/19_source_analysis.md` の「大見氏の DA」「大見氏の微分代数」 (SAD の `extensions/sad+`) | 敬称を外し「sad+ (DA (微分代数)、C++ 33 本)」「sad+ (微分代数の C++)」の形にする (下の 2.1 の補足) |
| 人名そのもの (ローマ字) 1 件 | `pysad/tests/test_ddba_lattices.py` の検査用の正規表現に、DDBA の作者名 (ローマ字と漢字) が入っていた | 文字列を分けて書き、ファイルに 1 つの語として載らないようにした |
| ログイン名 `coding` 7 ファイル | 試験の出力と台本 (`pysad/tests/scripts/*`) と `example/09_io_system/*.pysad.txt` の一時フォルダのパス `C:/Users/coding/AppData/Local/Temp/...` | `<temp>/...` に置き換えた |

補足 (許可の範囲の判断): SAD のソースの著者の表記 (`K. Oide`・`K. Ohmi`・`E. Forest`・KEK Report 92-14) は `pysad/__init__.py`・`teigen.py`・`und.py`・`phsrot.py`・`sadscript.py` のコメントと `docs/27_und.md` に残した (公開された報告書と K64 のソースの著者表記。ユーザーの判断が要る点は 4)。手法の名前に含まれる人名 (Piwinski、Bruck、Bassetti-Erskine、Oide-Yokoya、Yoshida、Haissinski、Dragt-Finn、Laskar-Robin、Bartolini-Schmidt) と、`example.oide.sad` (SAD 付属の例のファイル名) は残した。`K. Harada et al.` と KEK-LS の Raimondi 氏・Liuzzo 氏の記述は許可どおり。README 冒頭 4 行は触っていない (1 行目の「生出さん」以外に人名は無い。ほかに個人名が書かれていないことを確認した。2 行目以降に気になる点なし)。

`pyproject.toml` と `pysad/` の GitHub の URL (`Ken-Harada1234/pysad`) は公開の URL なので許可の範囲。

### 2.2 個人パス・設定

| 件数 | 内容 | 扱い |
|---|---|---|
| 文書中の `PYSAD_SAD_PACKAGES=D:/work/SAD-master/SAD-master/Packages` 15 ファイル 15 か所 (md) | docs/01、10、17、22、41、history/00_RESUME、v23_Q1・R5・S1・S2・T1・T2、work_example_fixes、example/v2_3/T1_README・T2_README | `<SAD の Packages のフォルダ>` に直した (環境変数名は残した) |
| 個人の作業フォルダの名前 5 か所 (md) | `D:\work\SAD\sad\pf2\sad2ele_*.sad` (v23_T1)、`D:\work\SAD-master\SAD-master\src` (v23_U1)、`D:\work\SAD\GUI` と `SAD-master\doc` (work_gui_tclbridge)、`D:\work\SAD-master\SAD-master` (work_wake)、`D:\work\Soft_FootPrint\Elegant_usage_ja.md` (T1_README) | `<SAD>`・「SAD 付属の GUI のフォルダ」などに直した。個人の文書の名前は削除し「elegant の公式文書」とした |
| コード中のハードコードされたパス 20 ファイル | `example/tools/*` (run_example・run_kekls・run_ddba_v07_14・canvas_render・gui_shot・ddba/gs_cmp)、`example/v2_3/*` (Q1・S1・S2・T1 の gs 比較、L2_common)、`example/v2_next/*` (P2・P5)、`example/ddba/v07_14/sextupole/sxt_common.py`、`example/tools/kek_ls_figs.py` のコメント | 新規 `example/tools/sadenv.py` (`PYSAD_SAD_PACKAGES` から Packages と gs の場所を決める。gs の場所だけ変えるときは `PYSAD_GS`) に置き換え。環境変数を設定した環境で同じ値になる (この環境で確認)。環境変数が無いときは、以前は個人のパスが既定だったが、いまはメッセージを出して止まる (setdefault だけの 9 本は、設定が無ければ GUI 系が skip になる従来の振る舞いに合わせて既定を外した) |
| 出力に出る Python のトレースバックのパス `D:\work\SAD-master\Claude outputs\...` 5 ファイル | `example/v2_next/*pysad_k64.txt`・`B1_sol_spin_bigbz.out.txt`、`example/v2_3/U2_verify.out.txt` | `<repo>\...` に置き換えた |
| `.gs.txt` に出る WSL のパス `/mnt/d/work/SAD-master/SAD-master/Packages/init.n` ほか | 約 200 ファイル、1450 行。gs が起動時に読み込むファイルの名前 | 許可の範囲 (既存の例と同じ)。個人のフォルダ名は `work` と配布フォルダ名 `SAD-master` だけなので、一般化はしていない (4 の判断) |
| WSL のホーム (`~/elegant/env.sh`・`~/sad_dbg`・`~/t1`・`~/t2`・`~/l2`・`~/u2`・`~/venvs/py314`) | 個人名を含まない慣例の場所。説明と WSL を呼ぶスクリプトの中 | 変更なし (報告のみ) |

### 2.3 抜け落ち・未更新・矛盾

| 観点 | 件数 | 代表例 | 扱い |
|---|---|---|---|
| (a) 件数・版の不一致 | 7 | `docs/01_overview.md` が 1691 件・1683 通過 (v2.3.0 の途中の値)、`docs/17_implementation_summary.md` が v2.1.0 の 1345 通過・`test_*.py` 118 ファイル・`test_v2_*.py` 54 本 (775 件)、`docs/history/00_RESUME.md` §3 の表が 578 件 (64 ファイル)。`pyproject.toml` の Changelog の URL が v2.1.0 のタグ | 現在の値 (1899 通過・8 skip・1907 件収集・158 ファイル・`test_v2_*.py` 73 本 1081 件) に直した。URL は `main` の CHANGELOG |
| (b) モジュールの地図 | 3 | `pysad/sxtopt.py` が docs/01・docs/12 に無い、`constants.py` が無い、`perfkern_bend.py` が無い。docs/12 に `diag_fit`・`diag_trace` などの旧い診断用 (非推奨・docs/01 に記載) と、`beambeam`・`coord`・`und`・`geo` ほかの要素別の docs への入口 (docs/01 の地図に書いてある) が無い | docs/01・docs/12 に追加。docs/12 の旧いモジュールは要素別の文書にあるので足さなかった (下の 3) |
| (c) API | 3 | docs/12 の `momentum_acceptance` の `check_every=0` (実際の既定は None = 1)、`resonance_terms` の引数が不足 (`rows`・`momentum_eV`・`order2`・`rdt4_options`)、`TrackSetup.from_source` の `fluc`。`compat.set_bend_k1_def`・`bend_dk_fix` が無い | 直した・追加した。`FrequencyMap[]` などの台本の関数、`adts(fourth=, optics=)`、`Formula->` は載っていた |
| (d) 索引 | 0 | `example/README.md` の一覧に全フォルダがある (`00_lattice` と `tools` には README が無い。一覧には載っている)。各フォルダの README の索引に、フォルダの全 `.sad` が載っている (DDBA の `*_full.sad` は表の「`*_full.sad`」の行で、v2_next の P5 は P5 の記録へのリンクの行で代表) | なし |
| (e) リンク切れ | 1 | `README.md` の `docs/html/` (追跡外の生成物。説明に「追跡しない」と書いてある)。`README_top.md` の 5 件は雛形で意図したもの | なし |
| (f) 古い記述 | 2 | docs/47 の v2.0.0 の節の N12 (`Twiss["DETR"]` が CALC の後で 0) が「未修正」のまま、CHANGELOG の v2.0.0 の「未解決」も同様 (v2.1.0 で修正済み。docs/05・docs/46 は修正済みと書く) | 「v2.1.0 で修正」と注記した |
| (f) 古い記述 | 1 | `README.md`・`CLAUDE.md`・docs/01 が `run_regression.py` を追跡対象のように書くが、ファイルはリポジトリに無い (作業者のローカル。`bench_level*.sad` も同じ) | README と CLAUDE.md に「リポジトリには入っていない」と注記 (docs/01 には書いてあった) |
| (g) B 番号 | 4 | `docs/46` の「B01〜B159」、`docs/17` の「B01〜B110」、`docs/11` の「B01〜B71」、`docs/history/00_RESUME.md` の「B1〜B62」 | B01〜B162 に直した (docs/21 の表は B1〜B162 の 162 行がそろっている)。docs/17 の「A1〜A58」は A1〜A60 |
| (h) 文書同士の矛盾 | 0 | `check_every` の既定 1、`adts()` の既定 (`fourth='r3'`・`optics='calc'`)、LMA と elegant の差 (0.5% 以内)、Touschek の式の差 (±5%) は、README・CHANGELOG・CLAUDE.md・docs/48・49・46・47 の間で一致 | docs/12 の `check_every=0` だけが食い違っていた ((c)) |

### 2.4 語句の不統一

| 項目 | 状況 | 扱い |
|---|---|---|
| 6 極 / 6 極、八極 / 8 極 | README・CHANGELOG・CLAUDE.md・docs では「6 極」が多数 (6 極は約 110、6 極は約 330)。pysad の docstring でも「6 極」が多数 (42 対 86)。example の md は「6 極」352・「6 極」242 | README・CHANGELOG・CLAUDE.md・docs (history と docs/13 を除く) の「6 極」「八極」を「6 極」「8 極」に統一 (24 ファイル)。リリースノート (ローカル) も同様。`pysad/*.py` の docstring と example の md は、example の md の多くが作成用のスクリプト (`example/tools/*`) の文字列から生成されること、出力と比べる試験があることから、変えていない (下の 3) |
| K64 モード / k64 モード | 628 対 1 | `pysad/emit6d.py` のコメント 1 か所を「k64 モード」に |
| 1/(m rad) / 1/(m·rad)、m rad / m·rad | 「m rad」が多数 (nm rad 89 対 nm·rad 0) | docs/49・50 と example/v2_3 の md 5 本の「m·rad」を「m rad」に (コード・出力・docs/history は変更なし) |
| タウシェック / Touschek | 2 件 (DDBA_26 の設計メモの要約) | 「Touschek」に (md 3、台本のヘッダ 1、`ddba_info.py` 1。DDBA の試験は通る) |
| KEK-LS / 大きなリング A、SP8 / 大きなリング B | 表記は「KEK-LS (大きなリング A)」「SP8 (大きなリング B)」で一貫 | なし |
| 周波数マップ / frequency map、LMA / 局所運動量アクセプタンス、RDT、ADTS、k64 モード、v2 モード | 日本語の文章では日本語を先に出し、略語は括弧で補う形で一貫 (「frequency map」は `freqmap.py` と `FMAP_summary.md` に括弧つきで 2 回) | なし |
| pysad の大文字小文字 | `PYSAD_` の環境変数とパッケージ名以外に「Pysad」「PySAD」なし | なし |

### 2.5 例の md の数値と出力 (抜き取り)

155 本の md の 2318 個の数値 (小数 4 桁以上) を、同じ名前の `.gs.txt`・`.pysad.txt`・`.json` の数値と、最後の桁の丸めの範囲で自動照合し、残り 263 個 (多くは解析式の列・単位換算・md の表の差の列) を手で調べた。不一致の訂正は 0 件。確かめたもの: `pf/tune_radiation/FFS_amp_dep_ts_PF.md` の Δνx 5σ・10σ・20σ・垂直 5σ (gs 出力 9.5983497550・9.6008159738・9.5982944172 と設計チューン 9.5975・5.2790 から、+0.000819・+0.003285・+0.000764 になる。md は gs の値による。pysad は 2e-5 違う)、`03_matching/FFS_FIT.md` の BX 7.73031・BY 1.36768 (gs 出力 7.730307・1.367681)、K1 0.321805・-0.643611、`04_emittance_radiation/FFS_RADINT.md` の I1 2.337645・I4 0.08023724 (gs 出力 8.0237245E-02)・I5、`FFS_EMITTANCE.md` の Jx 0.9651・U0 1.56327 MV (gs 出力 0.9651・1.5632657 MV)、`02_optics_cod/FFS_BEAMSIZE.md` の σx 1.28836 mm、`06_collective/flag_INTRA.md` の NOINTRA の 5.093949e-9、`03_matching/var_CONVERGENCE.md` の NX−6.7 (−0.007496)。

## 3. 直さなかったもの (理由)

- example の md の「6 極」「6 極」の混在 (約 600 か所): 多くは `example/tools/*` が文字列から生成し、md だけ直すと次の生成で戻る。作成用のスクリプトの文字列も直すと、生成物の差分が大きくなる。判断はユーザー (4)。
- `pysad/*.py` の docstring の「6 極」42 か所: コードの文字列 (出力・警告) との区別が要り、出力を比べる試験がある。コメントだけに限っても多い。
- docs/12 に載っていない旧いモジュール (`diag_fit`・`diag_trace` は旧い診断、`beambeam`・`coord`・`und`・`geo`・`multa`・`multae`・`radint`・`teigen`・`trptacc`・`tsvd`・`wakeopt`・`emit_extra`・`ffsextra`・`sadkw`): docs/01 の地図に載り、要素別の文書 (docs/22〜37) に API がある。重複を避けた。
- `.gs.txt` の `/mnt/d/work/SAD-master/SAD-master/Packages/...` (約 1450 行): 許可の範囲。
- `~/elegant`・`~/sad_dbg` ほか WSL のホームの慣例の場所: 個人名を含まない。
- docs/history の語句・数値の更新 (依頼どおり。個人パスと事実の誤りだけ直した)。
- `example/tools/pf_astrid/*` の「パスは整理前のまま」 (example/README.md の注記): 動くかどうかの確認はしていない。
- `local/` と `docs/html`: 追跡外。
- 例ごとの ○△× の件数: 合計 189・34・1・224 は表の小計と一致するが、フォルダの README の個別の判定の再集計はしていない。

## 4. ユーザーの判断が要るもの

1. SAD のソースの著者の表記 (`K. Oide`・`K. Ohmi`・`E. Forest`、KEK Report 92-14) をコメント・docstring に残すか。残した理由は、K64 のソースの出所の表記と公開の報告書の引用だから。消す場合は、`pysad/__init__.py`・`teigen.py`・`und.py`・`phsrot.py`・`sadscript.py` のコメントと `docs/27_und.md`。
2. 「大見氏」の敬称を外した 2 か所 (docs/19) で、人名のまま「大見氏の DA」と書かずに「sad+ (…)」にしたのでよいか。
3. example の md と `pysad/*.py` の「6 極」「6 極」を統一するか (統一する場合は、作成用のスクリプトの文字列と同時に)。
4. `.gs.txt` の WSL のパス (`/mnt/d/work/...`) を一般化するか。
5. `README.md` 冒頭 4 行は触っていない。
6. PyPI のビルドの作り直し (下の 5)。

## 5. リリース用のファイル

- **PyPI のビルド** (`D:\work\SAD-master\pysad-v2.3.0-build`、`release/pypi-v2.3.0`、先頭 `e70e13f`): main の README・CLAUDE.md・`pyproject.toml` (Changelog の URL)・`pysad/emit6d.py` (コメント 1 行)・docs を変えたので、**作り直しが要る** (README の PyPI 用の注記 1 コミットを、新しい main の上に載せ直し、同じ手順でビルドして `twine check`)。この点検では作り直していない (アップロードもしていない)。
- **リリースノート** (`local/wake/release_notes_v2.3.0.md`): 試験の件数・B 番号・例の数は変わっていない。「6 極」を「6 極」に直しただけ。
- push・タグ・GitHub のリリース・PyPI のアップロードはしていない。
