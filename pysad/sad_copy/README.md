# pysad/sad_copy: SAD のソースからコピーしたファイル

**このフォルダのファイルは、pysad が書いたものではありません。KEK の SAD (Strategic Accelerator Design、K64 版のソース配布 `SAD-master`) から、そのままコピーしたものです。** pysad を単独の (SAD を別に入れなくても動く) パッケージにするために、pysad に同梱しています。

| フォルダ | 元の場所 (SAD のソース) | 中身 |
|---|---|---|
| `Packages/` | `SAD-master/Packages/` | SADScript で書かれた SAD のパッケージ (`*.n`、71 本。`init.n`、`correction.n` (CorrectOrbit・Steer・Monitor)、`bump.n`、`Tkinter.n`、`KBFrame.n` ほか) と `SADInspect.sad` |
| `KBFrame/` | `SAD-master/KBFrame/` | KBFrame (SAD の GUI の枠組み) の設定ファイル (`*.res`) と画像 (`*.gif`・`*.xbm`) |
| `LICENSE` | `SAD-master/LICENSE` | SAD のライセンス (KEK Accelerator Laboratory、BSD 3 条項)。**再配布の条件として、この著作権表示と免責条項を残す** |
| `MANIFEST.sha256` | (pysad が作成) | コピーした各ファイルの SHA-256 (コピー時点のもの。元の SAD と同じ内容かの確認に使う) |

## コピーのしかた・コピーしていないもの

- コピーは、ファイルの中身を変えていません (バイト単位で同じ)。pysad 側の変更は、`pysad/pkgload.py` と `pysad/tclbridge.py` が、環境変数 `PYSAD_SAD_PACKAGES` が無いときにこのフォルダを使うようにしたことだけです。
- 次のものは**コピーしていません**: `Packages/init.local.n` とその `.bak`・`.example` (その計算機だけの設定)、`OpticsPlot_230724a.n`・`OpticsPlot_230724b.n` (退避用の古い版)、`*.py` (SAD の Python 連携用で、pysad は使わない)、`CVS` フォルダ。SAD の Fortran のソース (`src/*.f`) と実行ファイル `gs` も、コピーしていません (pysad は SAD なしで動きます)。
- コピーした日: 2026-11-26。元の SAD は、K64 版のソース配布 (`SAD-master`。SADHelp 1.4.0.0k64 の系統) です。SAD の最新版との差は確認していません。

## 使われ方

- `CorrectOrbit`・`Steer`・`Monitor`・`MakeBump` などの SADScript のパッケージは、最初に使われたときに、このフォルダの `Packages/*.n` が読み込まれます (`pysad/pkgload.py`)。Tcl/Tk の GUI の接点 (`pysad/tclbridge.py`) も、このフォルダの `Packages/*.n` と `KBFrame/` を使います。
- 別の SAD の `Packages` を使いたいときは、環境変数 `PYSAD_SAD_PACKAGES` にそのフォルダを指定します (指定があれば、そちらが優先されます)。

## 著作権・ライセンス

Copyright (c) 1986-2016, KEK Accelerator Laboratory. All rights reserved. 条件は `LICENSE` のとおり (BSD 3 条項)。このフォルダのファイルの著作権は KEK にあり、pysad の作者のものではありません。SAD の公式の情報: <http://acc-physics.kek.jp/SAD/licensing.html>。

---

# pysad/sad_copy: files copied from the SAD source

**The files in this folder were not written by pysad.** They are verbatim copies from KEK's SAD (Strategic Accelerator Design, the K64 source distribution `SAD-master`): `Packages/` (SADScript packages, 71 `*.n` files and `SADInspect.sad`), `KBFrame/` (KBFrame resources), `LICENSE` (SAD's BSD 3-clause license, which must be kept on redistribution) and `MANIFEST.sha256` (SHA-256 of every copied file, made by pysad). Copied on 2026-11-26 so that pysad is self-contained. `init.local.n` (machine-specific settings), backup files, the Python helper scripts, the Fortran sources and the `gs` executable were not copied. The environment variable `PYSAD_SAD_PACKAGES` overrides this folder. Copyright (c) 1986-2016 KEK Accelerator Laboratory.
