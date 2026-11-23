# tools: md を図入りの HTML にする

リポジトリ内の md を、目次つき・図入りの 1 枚の HTML にする。使うのは `markdown` パッケージだけで、外部通信・外部ライブラリの読み込みはない。HTML は md から作る生成物なので、md を直したら作り直す。

## 使い方

```
py -3.12 tools/build_html.py                  # 追跡対象の md をすべて変換
py -3.12 tools/build_html.py docs README.md   # 指定したファイル・フォルダの分だけ
py -3.12 tools/build_html.py --changed        # 変更された md だけ
py -3.12 tools/build_html.py --check          # 壊れたリンク・索引と入口のリンク先の検査だけ (HTML は書かない)
py -3.12 tools/build_html.py --require-math   # Node.js がない・数式の変換エラーがあるとエラー終了 (数式は「数式」の節)
py -3.12 tools/md_to_html.py README.md        # 1 ファイルだけ (-o 出力先、--fold-regex、--fold-level)
```

- 出力は、md と同じ名前の `.html` を同じフォルダに置く (名前に版・日付は付けない)。同名の `.html` が本ツールの生成物なら上書きする。そうでない既存の `.html` は、何もせず警告する。
- 引数なしの対象は `git ls-files '*.md'` のうち、作業ツリーに存在するもの。ただし `docs/history/` の下、`CHANGELOG.md`、`LICENSE` 類、`example/tools/kek_ls_src/` (生成用の断片) は除く。`--include-untracked` で、追跡されていない md も対象にする (コミット前の検査は `--check --include-untracked`)。
- `--changed` は、`git status` に出ている md と、HTML がない・HTML より新しい md を対象にする。ほかの引数と併用すると、その範囲の中で絞る。
- 変換の終わりに、全文検索の索引 `search_index.json`、検索の画面 `search.html`、入口ページ `index.html` をリポジトリの直下に作る (下の「全文検索と入口ページ」)。
- 終わりに、変換した件数・埋め込んだ図の枚数・警告・壊れたリンクの一覧を表示する。`--check` は壊れたリンク (md への相対リンク、画像、`#` のアンカー) があると終了コード 1。

## 変換の内容

- 左に固定の目次 (見出し 2〜3 階層)。表は見出し行を固定し、見出しのクリックで並べ替える。明暗の配色は OS の設定に従う。
- ページのタイトルは最初の h1。先頭に、元の md のファイル名だけの注記を付ける。
- 図: (a) `![説明](相対パス)` の画像 (png・jpg・gif・svg・webp) は base64 で埋め込む。1 枚が 3 MB を超えるときは埋め込まず相対リンクのままにして警告する。(b) 本文中にバッククォートで書かれた図のファイル名 (`fig_xxx.png` など) が、md と同じフォルダまたは `results/` にあれば、最初の言及を含む段落・箇条書き・表などのブロックの直後に図を置く。コードブロックの中の名前は無視する。画像記法で既に載せた図は重ねない。
- 文書間リンク `[..](xxx.md)`・`(xxx.md#anchor)` は、変換対象の md なら `.html` に書き換える。対象でなければそのまま。
- 見出しのアンカーは GitHub 風 (小文字化し、英数字・日本語・`_`・`-` 以外を除き、空白を `-` にする。同じ見出しは `-1`、`-2`)。日本語の見出しでも同じ。リンクの `#` が大文字小文字だけ違うときは、実在の ID に合わせる。
- コードブロックの色分け (下の「コードの色分け」)。各ページの先頭の注記の行に、検索と入口へのリンクを置く (`../search.html`・`../index.html` の相対パス)。
- 段落の直後に空行なしで始まる箇条書き・表、2 スペースの入れ子の箇条書きは、GitHub と同じように読めるよう前処理で整える (コードブロックの中は触らない)。

## 折りたたみ

既定では折りたたまない。`--fold-regex` (見出しの文字列に一致する正規表現) と `--fold-level` (既定 3) で、その階層の見出しの節 (同じか上の階層の見出しまで) を `<details>` にする。折りたたんだ節の中を指すリンクを開くと、その節が開く。

`build_html.py` では、`tools/html_config.json` に「リポジトリの根からの相対パス → オプション」を書く。

```json
{
  "example/ddba/v07_14/sextupole_v2/README.md": {"fold_regex": "^(ラウンド|揺らぎの確認|反省)", "fold_level": 3}
}
```

## 数式

md の中の TeX の数式は、HTML を作るときに KaTeX (Node.js) で HTML の部品にして埋め込む。できた HTML は、開くときに JavaScript が要らず、外部ファイルも外部通信もない単独のファイルで、オフラインの Chrome などでそのまま数式として見える。

書き方:

- インライン: `$\nu_x = 6.8$`。表示数式 (中央に独立した行、複数行も可): `$$ ... $$`。表の中にも書ける。
- インラインとみなすのは、同じ行の中で、開く `$` の直後と閉じる `$` の直前が空白でなく、中身が TeX の命令 (`\nu`)・`_`・`^`・`{}` を含むか、ASCII だけの短い式 (`$x$`・`$x=1$`・`$(x,y)$`・`$j+k+l+m=4$`・数だけの `$0.36$`) のもの。日本語を含む、閉じる `$` の直後が数字 (`$1〜$13`)、中身の前後が空白のもの (`$5 と $6`) は、金額として数式にしない。`$$$` (3 個以上連続。SAD の出力) も数式にしない。
- **コードの中は数式にならない**: ``` の囲み、字下げ 4 つ以上のコードブロック、バッククォートのインラインコード (`` `$x_1$` ``) の中の `$` は、そのまま。箇条書きの続きの字下げはコードとみなさない。
- 数式でない `$` を本文に書くときは、`\$` と書く (`\$5`)。シェルの `$HOME` のように、閉じる `$` がなければ、そのままでも数式にならない。
- 文書全体で数式を使わないときは、`tools/html_config.json` に `"math": false` を書く (SAD のヘルプの原文を貼った表など)。
- 変換できない数式は、赤字の TeX のまま残し、どの文書の何行目かを警告に出す。ビルドは止めない。

実装 (`mathtex.py`、`katex_render.js`):

- md を markdown に通す前に数式をプレースホルダーに置き換え (`_` や `*` が強調に誤変換されないように)、変換後に戻す。数式は文書ごとに 1 回の Node.js の呼び出し (標準入力に JSON、標準出力に JSON) でまとめて `katex.renderToString` にかける (`throwOnError` を使い、失敗した式だけ赤字にする)。出力は `htmlAndMathml` (見た目の HTML に、読み上げ・検索・コピー用の MathML を併記する。`html` だけだと読み上げられない。数式の部分は約 18% 増える)。`\lambdabar` は KaTeX になく、`katex_render.js` のマクロで `\bar{\lambda}` にしている。
- 数式のある HTML にだけ KaTeX の CSS を入れ、フォントは woff2 を base64 の data URI にして `@font-face` に直接入れる。入れるのは、数式の HTML の class から判定した使われているフォントの面だけ (Main の標準・Math の斜体・Size1〜3 など。フォントだけで文書あたり約 50〜120 KB)。数式のない HTML には CSS もフォントも入れない。長い表示数式は横にスクロールし、文字色は本文から継承するので暗い配色でも読める。
- 数式のための `<script>` は入らない (入っている `<script>` は表の並べ替えと折りたたみのためのもの)。

Node.js:

- 探す順は、環境変数 `PYSAD_NODE` (node 実行ファイルのパス) → PATH の `node` → `D:\work\tools\node-*\node.exe` (版が最新のもの)。公式の Windows 用 zip (`node-v*-win-x64.zip`) をどこかに展開して `PYSAD_NODE` にパスを入れるのがいちばん簡単。インストーラーは要らない。
- 見つからないときは、数式は TeX の文字のままにして、警告を 1 回出す (ビルドは続く)。`--require-math` (`build_html.py`・`md_to_html.py`) を付けると、Node.js が見つからないとき、または変換できない数式があるときに、エラー終了にする。
- KaTeX 0.19.1 (MIT ライセンス。https://katex.org/) の `katex.min.js`・`katex.min.css`・フォント (woff2 の 20 個) を `tools/katex/` に置いてある。ライセンスは `tools/katex/LICENSE`。ビルドで外部から取得することはない。

## コードの色分け

コードブロックは、ビルドのときに Pygments (`codehl.py`) で `<span class="…">` にする。開くときに JavaScript は要らず、CSS は色分けしたコードのあるページにだけ入る。明暗の配色は OS の設定に従い (`prefers-color-scheme`)、色相だけに頼らないよう、コメントは斜体、キーワードは太字にしてある。

- 対象: ` ```python `・` ```bash `・` ```sh `・` ```json `・` ```sad ` (SAD の台本)。
- 言語指定のないブロックは、SAD の台本と判断できるときだけ SAD として色分けする。判断は、空行を除く各行の 7 割以上が台本の形 (行末の `;`・FFS の命令・要素の定義・`Name[` の呼び出し・代入・コメント) で、3 割以上の行が `;` で終わり、数値の表・出力 (行頭の `$$$`・罫線・`|` の表)・コマンドの行 (`py`・`cd`・`set` など) がないこと。迷うものは色分けしない。出力の例・図・`text` は色分けしない。
- SAD の台本のレキサーは `sad_lexer.py` (簡易): コメント (`!` から行末、`(* *)`)、文字列、数、FFS の命令 (`USE`・`CALC`・`EMIT`・`DRAW` ほか)、要素 (`DRIFT`・`BEND`・`QUAD`・`SEXT`・`MARK`・`LINE`・`CAVI` ほか)、制御 (`Module`・`Do`・`If`・`With` ほか)、組み込み関数 (`Print`・`Fit`・`Show`・`Plot`・`TrackParticles` ほか)、記号。文法を解析するものではなく、書いたとおりの大文字小文字で一致したものだけに色が付く。
- 色分けしても文字列は変わらない (タグを除いて実体参照を戻すと元の文字列。合わなければ色分けしない)。コピーすると元のコードになる。
- 色分けの対象外: 箇条書き・番号つきリストの中に字下げして書いたコードブロック (Markdown の変換器がコードブロックとして扱わない)。
- Pygments がなければ色分けなしで HTML を作る (`pip install pygments`)。

## 全文検索と入口ページ

ビルドの最後に、次の 3 つをリポジトリの直下に作る (`site_index.py`)。対象の md は `docs/`・`example/`・`README.md`・`tools/README.md`・`pysad/tests/` の README・`physics_guide/` (`_deletable` と `docs/history/` は除く)。物理ガイドなどの md が増えれば、次のビルドで自動で入る。

- `search_index.json`: 索引。ページごとのパス・タイトル・要約・見出しと、節 (h1〜h3 の見出しごと) の本文 (数式は TeX の文字列、コードは文字列のまま。1 節 8000 字まで)。全角半角をそろえ (NFKC)、空白を 1 つに畳んである。
- `search.html`: 検索の画面。JavaScript は検索のためだけで、外部ライブラリはない。入力した語 (空白区切り) をすべて含む節を、大文字小文字・全角半角を区別せず、部分一致で探し (日本語・数式・コードにも効く)、見出し・タイトルに含まれるものを上位にして、周辺の文章つきで 60 件まで出す。範囲 (docs・example など) の絞り込みができる。`search.html?q=語` で検索語を渡せる。
- `index.html`: 入口ページ。目的 (docs/01_overview.md の最初の段落)、docs の目次 (docs/01_overview.md の目次の表の分類と、各 md の h1)、example の目次 (機能別フォルダの README と最初の段落)、試験・道具・物理ガイドへのリンク、検索へのリンク。md の見出しから自動で作るので、手で書かない。
- 一部の md だけを変換したときは、既存の索引のその分だけを入れ替える (入口と検索の画面は毎回作り直す)。`--check` は、索引の中のパスと、入口・検索のリンク先の存在を確かめる。
- 索引は `fetch` で読むので、`search.html` は Web サーバ経由で開く (GitHub Pages、または、リポジトリの直下で `py -3.12 -m http.server` を実行して `http://localhost:8000/search.html`)。ファイルを直接開く (`file://`) と、Chrome などは索引を読めない (画面にそのとおり表示する)。ほかの HTML (本文) はファイルを直接開いてよい。

## 図の作り方・SVG の方針

文書の図は、スクリプトで作る折れ線系の図は SVG、点の非常に多い図や画像は PNG にする。SVG は拡大してもぼけず、文字は読める。

- **保存の共通部分**: `tools/figures/figsave.py` の `save_figure`。`docs/figures/` の図 (`tools/figures/make_*.py`、`_common.py` の `save`) と、`example/` の図 (`example/tools/figlib.py` の `save`、`example/tools/kek_ls_figs.py`、`pf_history/make_figs.py`、`ddba/make_figs.py`) が使う。スクリプトには図の名前を `.png` で渡す。
- **SVG にする条件**: 図を SVG と PNG の両方に描き、SVG が PNG 以下、または 150 KB 以下なら SVG だけを置く (同名の PNG は消す)。そうでなければ PNG を置く。点の多い図 (周波数マップなどの散布図・画像・3D 図) は、`vector=False` を渡すか、SVG が大きくなるので PNG のままになる。
- **SVG の作り方**: `svg.fonttype = 'path'` で文字をパスにする (フォントのない環境でも崩れない。日本語フォントの設定はスクリプト側のまま)。ID は固定 (`svg.hashsalt`) で、同じ図なら同じ SVG になる。数値は小数 2 桁に丸め、作成情報・コメント・字下げを除く。文字をパスにするので、小さな図でも 40 KB 前後はある (同じ図の PNG の 1.5〜2.5 倍になることが多い)。PNG のほうが小さい図は、SVG が PNG 以下のものだけが SVG になる。
- **PNG の圧縮**: PNG は、画素を変えずに小さくして保存する (`optimize_png`: 不透明な RGBA は RGB に、256 色以下ならパレットに、フィルターは最良のものを選び、zlib 最大)。外部の道具は使わない。
- **再生成**: リポジトリの直下から `py -3.12 tools/figures/make_<名前>.py`、`py -3.12 example/<フォルダ>/figures_make_<名前>.py` など。`make_18_radlight.py` は `PYTHONPATH=pysad/tests/script_tests` が要る。環境変数 `PYSAD_FIG_SVG_ONLY=1` を付けると、PNG のままの図は書き直さない (SVG にできるかの確認用)。`make_09_emittance.py` は数分かかる。
- **SVG にした図**: `docs/figures/` の 06・07 (2 枚)・08・09・11・13・14・16・17・18、`example/` の 01 (3 枚)・02 (2 枚)・03・04・astrid (3 枚)・ddba (セルの光学)・kek_ls (誤差・ラティス・display の 3 枚・design の 2 枚)・pf (整合・履歴ラティスの 2 枚)。md の画像参照は `.svg`。
- **PNG のままの図**: 再生成のスクリプトがない図 (Tk の窓の画面、DRAW・OpticsPlot などの実行結果、周波数マップ・ADTS・RDT・LMA など長い計算の結果の図、DXF の図)、点の多い図、SVG が PNG より大きくなる図 (`KEKLS_ring_layout_optics`・`PF_history_layout`・`DDBA_overview`)。

## 限界

- 数式の検出は書き方の規則による (上記)。規則に合わない短い式 (たとえば `$x$` でなく `$ x $`) は数式にならない。KaTeX が解釈できない TeX (`\lambdabar` など、マクロに足していない記法) は赤字になる。
- 数式の見た目は Chrome で確認した。MathML は併記しているが、スクリーンリーダーでの読み上げは確かめていない。
- md の中の HTML タグはそのまま通す。
- 図の名前の言及は、同じフォルダと `results/` の直下だけを探す。
- 色分けと索引は、md の中の箇条書きに入れたコードブロックには効かない。SAD の台本の判定は、行の形による推定。
- 試験は `py -3.12 -m pytest pysad/tests/script_tests/test_tools_md_to_html.py -q`。

## GitHub Pages で公開する

リポジトリの HTML を、ダウンロードなしで、ブラウザからそのまま読めるようにする。サイトの URL は `https://ken-harada1234.github.io/pysad/` (入口ページ。全文検索は `…/pysad/search.html`)。リポジトリ側の準備 (下の 1〜2) はできている。GitHub の画面での設定 (3〜) は、リポジトリの所有者が行う。

1. **準備 (リポジトリ側)**: ルートの `.nojekyll` (空のファイル) が Jekyll の処理を止め、アンダースコアで始まるフォルダ・ファイルも、そのまま配る。入口 `index.html`・`search.html`・`search_index.json` と、各 md の `.html` を、`py -3.12 tools/build_html.py --include-untracked` で作り、コミットしておく。HTML 内のリンクは、すべて相対パスなので、`/pysad/` の下でもそのまま動く。
2. **大きさ**: 公開サイトの上限は約 1 GB (目安)。リポジトリの追跡ファイルの合計は約 110 MB (確認方法: `git ls-files -z | xargs -0 du -cb`)。余裕がある。
3. **リポジトリが public であることを確認する**: 無料のプランでは、public のリポジトリだけが公開できる。リポジトリのページの名前の横に `Public` と出ていれば public。または Settings → General → いちばん下の Danger Zone の Change repository visibility で確かめる。
4. **公開の設定**: リポジトリの Settings → Pages → Build and deployment → Source を **Deploy from a branch** にし、Branch を **main** (フォルダは **/ (root)**) にして **Save**。HTML が入っているブランチを選ぶ (HTML をコミットしたブランチを main に取り込んでおく)。
5. **公開まで**: 保存すると、Actions タブに `pages build and deployment` が走る。1〜数分 (長くても 10 分ほど) で、Settings → Pages の上部に「Your site is live at …」と URL が出る。初回は反映に少し時間がかかることがある。
6. **更新のしかた**: md を直したら、`py -3.12 tools/build_html.py --include-untracked` で HTML・索引を作り直し、md と HTML を一緒にコミットして main に push する。push のたびに自動で再公開される (数分)。
7. **非公開にする**: Settings → Pages → Build and deployment の **Unpublish site** (または Source を None に変更) で、サイトを止める。リポジトリを private にしても、無料のプランでは公開が止まる。止めたあとも、リポジトリそのもの (md・HTML のファイル) は残る。
8. **公開される範囲**: Pages はリポジトリの中の追跡ファイルをすべて配る (HTML だけではなく、台本・出力・ソースも URL を知っていれば読める)。公開してよい内容か、public にする前に確かめる。
