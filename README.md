- 本レポジトリ（作業メモ、解説文書を含む）は生出さんのソースをもとに、Claude code で全自動生成したものです。
- 基本的に、自分の書いたSADスクリプトに対する動作確認しかしておりません。
- うまくいかない場合は、そのままAIに渡して、作業メモやSADメモを読ませると、修正を続けることができます。（今後も随時修正すると共に、サンプルSADスクリプトを載せていきたいと思っています。）
- 下記、readmeもaiによる自動生成です。

# pysad — SAD (K64 版) の Python 直訳

加速器設計コード SAD (Strategic Accelerator Design) の K64 版 ([KatsOide/SAD](https://github.com/KatsOide/SAD)) の
光学計算・エミッタンス・マッチング・粒子追跡・ダイナミックアパーチャを Python に直訳したもの。
SAD の Fortran の式・分割数・反復順序・乱数の引き方を写し、SAD の数値 (印字精度 1e-10〜1e-11) を再現することを基準にしている。
SADScript のインタプリタと FFS のコマンド層を持ち、SAD のスクリプトをそのまま実行できる。

- 文書は `docs/` (`docs/01_overview.md` から。作業メモは `docs/history/00_RESUME.md` と `docs/history/work_*.md`)。
- SAD の仕様と癖の一覧は `docs/11_sad_quirks.md`、SADHelp との対応は `docs/14_pysad_vs_sadhelp.md`。
- SAD 本体はこのリポジトリに含まない。参照先は K64 版のソース (`src/*.f`, `Packages/*.n`) と CERN の SADHelp。

## 必要なもの

Python 3.12 と numpy, scipy, numba, pytest (`py -3.12 -m pip install numpy scipy numba pytest`)。

## 使い方

SAD の `gs` に渡すのと同じ順にファイルを並べる (複数ファイルは連結される。`Get[]` は最後のファイルのフォルダ基準)。

```
py -3.12 run_sad.py lattice.sad script.sad > script.pysad.out
py -3.12 run_sad.py lattice_and_script.sad              # ラティスとコマンドを 1 ファイルにしてあれば 1 つ
```

Python から:

```python
from pysad.ffscmd import run_sad_string
ffs, I = run_sad_string(open('script.sad', encoding='utf-8').read(), cwd='.')
print(I.get('EMITX'), ffs.rows[-1].betx)
ffs.run_text('fit NX 47.3; free QF* K1; go;')
```

下位の API (Lattice, Matcher, emittance6d, TrackParticles など) は `docs/13_pysad_api.md`。

## 検証

- `py -3.12 -m pytest pysad/tests -q` (`--ignore` する旧式の 4 本は `docs/01_overview.md` §5)。
  このリポジトリには pysad 本体と文書だけを入れている。SAD の参照値付きベンチ、回帰スクリプト、ユーザーの実スクリプトの照合は含めていない。
  誤差入りラティスを使う `test_elkern.py` の 2 件は環境変数 `PYSAD_ERR_LATTICE` (SAD スクリプトのパス) を指定したときだけ走る (無ければ skip)。

## 構成

| 場所 | 内容 |
|---|---|
| `pysad/` | 本体 (`docs/01_overview.md` §3 のモジュール地図) と `tests/` |
| `run_sad.py` | SAD スクリプトの実行 (`gs` と同じ引数) |
| `docs/` | 解説文書 (現状のスナップショット。`01_overview.md` から) |
| `docs/history/` | 履歴・作業メモ (`00_RESUME.md` は時系列の作業メモと残課題、`work_*.md` は項目ごとの経緯と照合の記録) |

## ライセンス

`LICENSE` は SAD の利用条件 (Copyright (c) 1986-2016, KEK Accelerator Laboratory。https://acc-physics.kek.jp/SAD/how-to-get-sources/licensing/ と K64 版の LICENSE に同じ文面) をそのまま置いたもので、文面は BSD 3 条項ライセンスと同じ。pysad は SAD の派生物として、この著作権表示と条件をそのまま適用する。
