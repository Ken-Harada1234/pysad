- 本リポジトリ（作業メモ、解説文書、実例を含む）は、生出さんのソースをもとに、Claude Code で全自動生成したものです。
- 動作確認は、基本的に自分で書いた SAD スクリプトに対してしかしていません。
- うまくいかない場合は、そのまま AI に渡して、作業メモや SAD メモを読ませると、修正を続けられます。（随時修正しています。サンプルの SAD スクリプト（それも自動生成）も載せました。）
- 以下の README も AI による自動生成です。

> **注意 / Notice**: PyPI のパッケージ名は `pysad-k64` です。PyPI にある異常検知ライブラリ `pysad` (Streaming anomaly detection) とは**無関係の別物**です。import 名はどちらも `pysad` なので、**同じ Python 環境に両方をインストールしないでください** (仮想環境を分けてください)。 / The PyPI package `pysad-k64` is **unrelated** to the streaming anomaly-detection library `pysad`. Both use the import name `pysad`; do not install both into the same environment.

# pysad — SAD (K64 版) の Python 直訳

> **v1 と v2 の違い (注釈)**: 以下の「直訳」は **v1 系 (v1.0.0) の方針**です。**v2 系 (v2.0.0 以降) は直訳ではありません。** v2 は K64 を出発点に、K64 のバグの修正、物理の誤りの修正、機能の追加を行った版で、**gs (K64) と結果が変わる箇所があります** ([`docs/47_v2_changes.md`](docs/47_v2_changes.md))。gs と同じ結果が要るときは、環境変数 `PYSAD_COMPAT=k64` (互換モード。v1 と同じ K64 の動作) を指定してください。なお、v2 でも、K64 の式・分割数・反復順序を写したコードの大部分は、そのままです。

加速器設計コード SAD (Strategic Accelerator Design) の K64 版 ([KatsOide/SAD](https://github.com/KatsOide/SAD)) の
光学計算・エミッタンス・マッチング・粒子追跡・ダイナミックアパーチャを Python に直訳したもの。
SAD の Fortran の式・分割数・反復順序・乱数の引き方を写し、SAD の数値 (印字精度 1e-10〜1e-11) を再現することを基準にしている。
SADScript のインタプリタと FFS のコマンド層を持ち、SAD のスクリプトをそのまま実行できる。

- 文書は [`docs/`](docs/) ([`docs/01_overview.md`](docs/01_overview.md) から。作業メモは [`docs/history/00_RESUME.md`](docs/history/00_RESUME.md) と `docs/history/work_*.md`)。
- SAD の仕様と癖の一覧は [`docs/11_sad_quirks.md`](docs/11_sad_quirks.md)、SADHelp との対応は [`docs/14_pysad_vs_sadhelp.md`](docs/14_pysad_vs_sadhelp.md)。
- 要素・機能別の文書は `docs/22_*.md`〜`docs/37_*.md` (GUI/Tcl、WAKE、BEAMBEAM、COORD、TCAVI、UND、MULT、SOL、MAP、DELCOR、GaussianCoulomb・BBBrem、特殊関数・Fit、SAD パッケージ、FFS の小コマンド、POL・CSR、描画)。
  大きな項目の解説文書は [`docs/38_wake_guide.md`](docs/38_wake_guide.md) (WAKE)、`39_beambeam_guide.md`、`40_csr_guide.md`、`41_tcltk_guide.md`、`42_graphics_guide.md`、`43_pol_guide.md`、`44_sol_guide.md`、`45_radlight_guide.md`。
- 現状 (2026-11-18、v2.1.0。PyPI は v2.0.0 まで公開): pytest 1345 通過・8 skip・失敗 0 (1353 件収集。`PYSAD_SAD_PACKAGES` つき。担当 E の 2026-11-18 の実行)、回帰 (`QUICK=1`) は 50 ベンチ FAIL 0。gs と合わない箇所・未対応は [`docs/history/work_unimpl.md`](docs/history/work_unimpl.md) と [`docs/history/00_RESUME.md`](docs/history/00_RESUME.md) §6 (例: SOL 領域内の四極 + RAD (gs の疑い)、最初の `Window[]` の親 (gs は Toplevel)、`SynchroBetaEmittance` の大きなリングの不一致、`Polarization` の過小評価)。未対応・未確認の最新の一覧は [`docs/46_unimplemented_and_unverified.md`](docs/46_unimplemented_and_unverified.md)。
- **互換モード (v2.0.0 から)**: 既定は `v2` (K64 のバグを直した動作)。gs (K64) と同じ結果が要るときは環境変数 `PYSAD_COMPAT=k64` (または `pysad.compat.set_mode('k64')`)。既存の試験と回帰は gs との比較なので k64 モードで走る。v1 との違いの一覧は [`docs/47_v2_changes.md`](docs/47_v2_changes.md)、K64 側のバグは [`docs/21_sad_bug_report.md`](docs/21_sad_bug_report.md) (B01〜B134)、変更履歴は [`CHANGELOG.md`](CHANGELOG.md)。
- SAD 本体はこのリポジトリに含まない。参照先は K64 版のソース (`src/*.f`, `Packages/*.n`) と CERN の SADHelp。

## インストール (PyPI)

```
pip install pysad-k64
```

PyPI: <https://pypi.org/project/pysad-k64/>。import 名は `pysad` です (パッケージ名と違います)。**v2.0.0** (K64 のバグの修正と機能の追加。互換モード `PYSAD_COMPAT=k64` で K64 と同じ動作) が最新です。ソースに忠実な翻訳の **v1.0.0** は `pip install pysad-k64==1.0.0` で入ります。

## 必要なもの

Python 3.12 と numpy, scipy, numba, pytest (`py -3.12 -m pip install numpy scipy numba pytest`)。

確認した Python の版 (v2.0.0、全試験): **3.12** と **3.10** は全通過 (937 通過・4 skip)。**3.14** は 935 通過・2 失敗 (3.14 に同梱の Tcl/Tk が 9.0 で、試験が 8.6 を期待していること、数値的に脆い 1 試験の閉軌道の収束の違い)。v2.1.0 では 2 件とも直し、3.14 で GUI 系と D2 の試験 (計 115 件) が通ることを確かめた (全試験は 3.14 では走らせていない。詳細は [`docs/46_unimplemented_and_unverified.md`](docs/46_unimplemented_and_unverified.md) §3)。3.13 は未確認。

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

下位の API (Lattice, Matcher, emittance6d, TrackParticles など) は [`docs/12_pysad_api.md`](docs/12_pysad_api.md)。

## 検証

- `py -3.12 -m pytest pysad/tests -q --ignore=pysad/tests/test_bend_feeddown.py --ignore=pysad/tests/test_coord.py --ignore=pysad/tests/test_sext_radiation.py` (`--ignore` する旧式の 3 本。手順は [`docs/history/00_RESUME.md`](docs/history/00_RESUME.md) §4)。
  GUI 系 (`test_plot.py`、`test_tclbridge.py`) は環境変数 `PYSAD_SAD_PACKAGES` (K64 の `Packages` の場所) が無いと skip。`PYSAD_SLOW=1` で遅い試験 (CSR MakeZL、約 7 分) も走る。
  このリポジトリには pysad 本体と文書だけを入れている。SAD の参照値付きベンチ、回帰スクリプトは含めていない。
  誤差入りラティスを使う `test_elkern.py` の 2 件は環境変数 `PYSAD_ERR_LATTICE` (SAD スクリプトのパス) を指定したときだけ走る (無ければ skip)。

## 構成

Git で追跡するのは [`pysad/`](pysad/)、[`docs/`](docs/)、`run_sad.py`、`README.md`、`LICENSE`、`CLAUDE.md`、`.gitignore` ([`docs/history/00_RESUME.md`](docs/history/00_RESUME.md) §9)。
ベンチ・回帰スクリプト・出力は入れない (`.gitignore`)。以下の一覧は 2026-11-18 の作業ツリー。

### トップ

| ファイル | 内容 |
|---|---|
| [`run_sad.py`](run_sad.py) | SAD スクリプトの実行 (`gs` と同じ引数: `py -3.12 run_sad.py lattice.sad script.sad`) |
| [`README.md`](README.md) | この文書 |
| [`CLAUDE.md`](CLAUDE.md) | AI アシスタント向けの読み順と要点 |
| [`CHANGELOG.md`](CHANGELOG.md) | 変更履歴 (リリースごとの内容、未対応・未確認の要点) |
| [`LICENSE`](LICENSE) | SAD の利用条件 (KEK Accelerator Laboratory) |
| [`.gitignore`](.gitignore) | 追跡しないもの (出力、ログ、キャッシュ) |
| [`pysad/`](pysad/) | 本体 (下の表) |
| [`pysad/tests/`](pysad/tests/) | 試験 (下の表) |
| [`example/`](example/) | 機能・関数ごとの短い例 (gs と pysad の画面表示つき。下の表) |
| [`docs/`](docs/) | 解説文書 (下の表。現状のスナップショット) |
| [`docs/history/`](docs/history/) | 作業履歴 (下の表) |

### `pysad/` (本体)

**SADScript と FFS**

| ファイル | 内容 |
|---|---|
| [`__init__.py`](pysad/__init__.py) | パッケージの入口 |
| [`parser.py`](pysad/parser.py) | SAD の MAIN レベル (ラティス定義、要素、LINE) と SADScript の読み取り |
| [`sadscript.py`](pysad/sadscript.py) | SADScript (SAD の Mathematica 風言語) のインタプリタ |
| [`sadlist.py`](pysad/sadlist.py) | SADScript のリスト・式の操作関数 |
| [`sadclass.py`](pysad/sadclass.py) | SAD の Class (`Class[...]`、`obj@member`、Constructor / Destructor) |
| [`sadelem.py`](pysad/sadelem.py) | 要素とビームラインに関する FFS 用の関数 (`LINE[...]`、`Element[...]` など) |
| [`sadkw.py`](pysad/sadkw.py) | 要素型ごとのキーワード一覧 |
| [`sadrandom.py`](pysad/sadrandom.py) | SAD の既定の乱数 (random_plugin_sad.c の直訳) |
| [`vecfn.py`](pysad/vecfn.py) | `Vector` / `ReleaseVector` |
| [`miscfn.py`](pysad/miscfn.py) | 小さな組み込み関数と `Packages/optimize.n` の翻訳 |
| [`specfn.py`](pysad/specfn.py) | 特殊関数・数学関数 |
| [`specfn_cx.py`](pysad/specfn_cx.py) | 複素の `Zeta`・`HurwitzZeta`・`PolyLog`・`HypergeometricU` (v2.1.0) |
| [`fitroot.py`](pysad/fitroot.py) | `FindRoot`・`Fit`・`FitGaussian`・`FitEmit`・`NIntegrate` |
| [`fileops.py`](pysad/fileops.py) | ファイル・ディレクトリ関数 (`Packages/FileSystem.n` の置き換え) |
| [`pkgload.py`](pysad/pkgload.py) | SAD の `Packages/*.n` を SADScript のまま読み込む仕組み |
| [`ffscmd.py`](pysad/ffscmd.py) | FFS のコマンド層 (`tffsa.f` の再現) と FFS 系の組み込み関数。`TrackParticles` もここ |
| [`ffsextra.py`](pysad/ffsextra.py) | FFS の小さなコマンド (VARIABLES・ATTRIBUTE・MATRIX・SHOW など) |
| [`maintrack.py`](pysad/maintrack.py) | MAIN レベルの `TRACK` コマンド |
| [`radint.py`](pysad/radint.py) | `RADINT` (放射積分) |
| [`nbcache.py`](pysad/nbcache.py) | numba のキャッシュの無効化 (`py -3.12 -m pysad.nbcache`。kernel を編集したら試験の前に走らせる。v2.1.0)。docs/12 §16.7 |
| [`compat.py`](pysad/compat.py) | 互換モード (`v2` が既定、`PYSAD_COMPAT=k64` か `compat.set_mode('k64')` で K64 (gs) の動作)。docs/47 |
| [`elements.py`](pysad/elements.py) | 要素の定義 |
| [`constants.py`](pysad/constants.py) | 物理定数 |

**光学・エミッタンス・マッチング**

| ファイル | 内容 |
|---|---|
| [`lattice.py`](pysad/lattice.py) | CELL/INS の光学計算 (qcell・qcod・qtwiss の直訳) |
| [`teigen.py`](pysad/teigen.py)、[`tsvd.py`](pysad/tsvd.py) | 固有値・固有ベクトル (teigen.f)、特異値分解ソルバ (tsvd) |
| [`turnkern.py`](pysad/turnkern.py)、[`perfkern.py`](pysad/perfkern.py)、[`perfkern_bend.py`](pysad/perfkern_bend.py)、[`elkern.py`](pysad/elkern.py) | 性能用の numba カーネル (結果は Python 版と同じ演算順) |
| [`emit6d.py`](pysad/emit6d.py) | `Emittance[]` (6 次元の平衡エミッタンス) |
| [`emit_extra.py`](pysad/emit_extra.py) | `Emittance[]` の戻り規則 (OneTurnInformation・Orbit・Matrix) |
| [`synchrobeta.py`](pysad/synchrobeta.py) | `SynchroBetaEmittance` |
| [`calc6d.py`](pysad/calc6d.py) | `CALC6D` (6 次元 CALC) |
| [`match6.py`](pysad/match6.py)、[`match6d.py`](pysad/match6d.py) | `FIT`/`FREE`/`GO` と解析微分 (tdgeo の GX.. を含む) |
| [`touschek.py`](pysad/touschek.py) | `INTRA`・`TouschekTable`・`TouschekLifetime` |
| [`trptacc.py`](pysad/trptacc.py) | TRPT (輸送線) の加速 |
| [`geo.py`](pysad/geo.py)、[`sgeo.py`](pysad/sgeo.py) | 要素の幾何 (GX・GY・GZ・CHI)、SOL 領域の幾何 |
| [`optics_fn.py`](pysad/optics_fn.py) | 光学・追跡の補助関数 (TransferMatrix など) |
| [`corr.py`](pysad/corr.py) | 軌道補正の組み込み (`ResponseMatrix`・`LinearSolveConditioned`) |

**追跡と物理の機能**

| ファイル | 内容 |
|---|---|
| [`track6.py`](pysad/track6.py) | `TrackParticles` の Python 版 (tturn の直訳)。POL のスピン追跡もここ |
| [`track6n.py`](pysad/track6n.py)、[`track6n_ext.py`](pysad/track6n_ext.py) | 同じ追跡の numba 版、MULT・SOL 領域の numba カーネル |
| [`dapert.py`](pysad/dapert.py) | `DynamicApertureSurvey` |
| [`wake.py`](pysad/wake.py)、[`wakeopt.py`](pysad/wakeopt.py) | ウェイク場 (追跡と光学) |
| [`beambeam.py`](pysad/beambeam.py)、[`bbbrem.py`](pysad/bbbrem.py)、[`gausscoul.py`](pysad/gausscoul.py) | BEAMBEAM 要素、BBBrem1、GaussianCoulomb |
| [`csr.py`](pysad/csr.py) | CSR (CSRTrack・CSRHaissin・CSRInit・CSRMatrix ほか) |
| [`spin.py`](pysad/spin.py) | POL (スピンと偏極。v2.1.0 でシンクロトロン側帯の項) |
| [`naff.py`](pysad/naff.py) | NAFF 風のチューン抽出 (`naff`・`naff_real`・`tune_xp`・`tunes_from_coords`。v2.1.0)。docs/12 §16.7 |
| [`radfield.py`](pysad/radfield.py) | `RadiationField`・`RadiationSpectrum` |
| [`measure.py`](pysad/measure.py) | `tlum`・`tluma`、MEASURE・TRCOD・LOGMEASURE・MAIN の `TRACK TURNS=-1` (v2 モード)、TRCOD FIXCOD など |
| [`coord.py`](pysad/coord.py) | COORD (座標変換) |
| [`tcav.py`](pysad/tcav.py) | TCAVI (横偏向空洞) |
| [`und.py`](pysad/und.py) | UND (アンジュレータ。v2 は磁場 `FBX`・`FBY` を含む) |
| [`phsrot.py`](pysad/phsrot.py) | PHSROT (位相の回転・減衰・拡散。v2 モード) |
| [`radtaper.py`](pysad/radtaper.py) | RADTAPER (閉軌道の局所の運動量に磁石の強さを合わせる) |
| [`spch.py`](pysad/spch.py) | 空間電荷 SPAC・PSPAC (SPCH)・WSPAC (v2 モード) |
| [`multa.py`](pysad/multa.py)、[`multae.py`](pysad/multae.py) | MULT (ANGLE≠0) の追跡、COD 写像と行列 |
| [`solmap.py`](pysad/solmap.py)、[`solemit.py`](pysad/solemit.py) | SOL 領域の Twiss と Emittance |

**GUI・描画**

| ファイル | 内容 |
|---|---|
| [`tclbridge.py`](pysad/tclbridge.py) | SAD の Tcl/Tk 接点を Python の tkinter で実装 |
| [`tkcanvas.py`](pysad/tkcanvas.py) | Tk の Canvas の組み込み関数 |

**旧 API (SAD と未照合の独自実装。文書の表の対象外)**: `ffs.py`、`matching.py`、`optics.py`、`tracking.py`。診断用: `diag_fit.py`、`diag_trace.py`。

### `pysad/tests/` (試験)

1345 通過・8 skip・失敗 0 (1353 件収集、2026-11-18。`PYSAD_SAD_PACKAGES` つき。v2.0.0 は 937 通過・4 skip)。参照は多くが SAD K64 (`gs`) の出力で、既存の試験は互換モード k64 で走る ([`conftest.py`](pysad/tests/conftest.py))。試験が走らせる SAD スクリプトと gs・pysad の画面表示は [`pysad/tests/scripts/`](pysad/tests/scripts/)、対応表は [`pysad/tests/OUTPUTS.md`](pysad/tests/OUTPUTS.md)。

| ファイル | 試験の対象 |
|---|---|
| [`conftest.py`](pysad/tests/conftest.py) | 試験の共通設定 (互換モードの既定を `k64` にする。フィクスチャ `v2_mode`・`k64_mode`) |
| `test_v2_A_*.py`、`test_v2_B_*.py`、`test_v2_C_*.py`、`test_v2_C2_*.py` | v2.0.0 の変更の試験 (29 ファイル 363 件: 放射・エミッタンス・スピン、特殊関数・スクリプト・FFS、追跡・集団効果・要素、残りの仕上げ)。v2 の動作は試験の中で `compat.set_mode('v2')` を明示する。一覧は [`docs/47_v2_changes.md`](docs/47_v2_changes.md) |
| `test_v2_A1_next.py`、`test_v2_A2_*.py`、`test_v2_B1_*.py`、`test_v2_B2_*.py`、`test_v2_B3_*.py`、`test_v2_C3_*.py`、`test_v2_D1_emit.py`、`test_v2_D2_*.py`、`test_v2_E_*.py` | v2.1.0 の変更の試験 (25 ファイル 412 件: 放射・F1・CAVI・MASS/CHARGE、スピン角・Polarization・RADINT、Emittance の F1・WSPAC・SOL、FREE・複素の特殊関数・Parallelize、側帯・光子・NAFF、MEASURE・RadiationField・NPARA、NORFSW・光子・POL + 加速、MAIN の式・Micado・nbcache、傾いたソレノイドの POL・ステアリングの F1)。補助の独立な計算は `_a2_indep.py`・`_c3_indep.py`・`_e_spinkick.py`。遅い試験は `PYSAD_SLOW=1` |
| [`test_pysad.py`](pysad/tests/test_pysad.py)、[`test_new_features.py`](pysad/tests/test_new_features.py) | 基本の数値検証、新機能 |
| [`test_item4.py`](pysad/tests/test_item4.py)、[`test_small_items.py`](pysad/tests/test_small_items.py)、[`test_error_spreading.py`](pysad/tests/test_error_spreading.py) | SADScript のリスト関数・tsvd・小項目 A〜E・誤差の spreading |
| [`test_evalfixes.py`](pysad/tests/test_evalfixes.py)、[`test_symbol_names.py`](pysad/tests/test_symbol_names.py)、[`test_symbolic.py`](pysad/tests/test_symbolic.py)、[`test_miscfn.py`](pysad/tests/test_miscfn.py) | gs に合わせた評価器の修正、シンボル名、記号微分 D、DownhillSimplex など |
| [`test_specfn.py`](pysad/tests/test_specfn.py)、[`test_fitroot.py`](pysad/tests/test_fitroot.py) | 特殊関数・数学関数、FindRoot・Fit・NIntegrate |
| [`test_fileops.py`](pysad/tests/test_fileops.py)、[`test_ffsextra.py`](pysad/tests/test_ffsextra.py)、[`test_sadrandom.py`](pysad/tests/test_sadrandom.py) | ファイル関数、FFS の小コマンド、乱数 |
| [`test_perf_kernels.py`](pysad/tests/test_perf_kernels.py)、[`test_elkern.py`](pysad/tests/test_elkern.py) | 性能用カーネルの回帰 (従来の経路と同じ結果) |
| [`test_emit_jac.py`](pysad/tests/test_emit_jac.py)、[`test_emit_report.py`](pysad/tests/test_emit_report.py)、[`test_emit_unstable.py`](pysad/tests/test_emit_unstable.py)、[`test_emit_extra.py`](pysad/tests/test_emit_extra.py) | Emittance・EMIT 報告・不安定な光学 |
| [`test_synchrobeta.py`](pysad/tests/test_synchrobeta.py) | SynchroBetaEmittance |
| [`test_geo_fit_sig.py`](pysad/tests/test_geo_fit_sig.py)、[`test_tdgeo.py`](pysad/tests/test_tdgeo.py) | 幾何の fit・ビームサイズ、tdgeo の解析微分 |
| [`test_unstable.py`](pysad/tests/test_unstable.py)、[`test_fam_points.py`](pysad/tests/test_fam_points.py)、[`test_trpt_acc.py`](pysad/tests/test_trpt_acc.py)、[`test_rotate.py`](pysad/tests/test_rotate.py)、[`test_cavrot.py`](pysad/tests/test_cavrot.py) | 不安定ライン、FITP の fam 点、TRPT の加速、ROTATE、空洞の DX/DY/ROTATE |
| [`test_bendk1_track.py`](pysad/tests/test_bendk1_track.py)、[`test_bend_feeddown.py`](pysad/tests/test_bend_feeddown.py)、[`test_sext_radiation.py`](pysad/tests/test_sext_radiation.py) | K1≠0 の BEND の追跡、フィードダウン、六極の放射 (後ろ 2 本は旧式で `--ignore`) |
| [`test_track6n_ext.py`](pysad/tests/test_track6n_ext.py)、[`test_main_track.py`](pysad/tests/test_main_track.py) | MULT・SOL の numba 追跡、MAIN の TRACK |
| [`test_coord.py`](pysad/tests/test_coord.py)、[`test_coord_elem.py`](pysad/tests/test_coord_elem.py)、[`test_tcav.py`](pysad/tests/test_tcav.py)、[`test_und.py`](pysad/tests/test_und.py)、[`test_mult_angle.py`](pysad/tests/test_mult_angle.py)、[`test_map.py`](pysad/tests/test_map.py) | COORD・TCAVI・UND・MULT ANGLE≠0・MAP (`test_coord.py` は旧式で `--ignore`) |
| [`test_sol_optics.py`](pysad/tests/test_sol_optics.py)、[`test_sol_emit.py`](pysad/tests/test_sol_emit.py)、[`test_sol_geo.py`](pysad/tests/test_sol_geo.py) | SOL の Twiss・Emittance・幾何 |
| [`test_wake.py`](pysad/tests/test_wake.py)、[`test_wake_optics.py`](pysad/tests/test_wake_optics.py) | WAKE の追跡、WAKEOPT |
| [`test_beambeam.py`](pysad/tests/test_beambeam.py)、[`test_bbbrem.py`](pysad/tests/test_bbbrem.py)、[`test_gausscoul.py`](pysad/tests/test_gausscoul.py)、[`test_measure.py`](pysad/tests/test_measure.py) | BEAMBEAM、BBBrem1、GaussianCoulomb、tlum |
| [`test_csr.py`](pysad/tests/test_csr.py) | CSR |
| [`test_spin.py`](pysad/tests/test_spin.py)、[`test_spin_track.py`](pysad/tests/test_spin_track.py) | POL (Emittance のスピン、TrackParticles のスピン追跡) |
| [`test_radfield.py`](pysad/tests/test_radfield.py)、[`test_radlight.py`](pysad/tests/test_radlight.py) | RadiationField、RADLIGHT |
| [`test_delcor.py`](pysad/tests/test_delcor.py)、[`test_corr.py`](pysad/tests/test_corr.py)、[`test_optics_fn.py`](pysad/tests/test_optics_fn.py) | DELCOR、軌道補正、光学の補助関数 |
| [`test_sadclass.py`](pysad/tests/test_sadclass.py)、[`test_gui_lang.py`](pysad/tests/test_gui_lang.py)、[`test_tclbridge.py`](pysad/tests/test_tclbridge.py)、[`test_plot.py`](pysad/tests/test_plot.py) | Class、GUI の言語機能、Tcl/Tk の接点、描画 (後 2 本は `PYSAD_SAD_PACKAGES` が必要) |

### `example/` (SAD の機能・関数ごとの短い例)

SAD の機能・関数を 1 つずつ短い例にしたものです。1 つの例は 4 ファイルで、`NAME.sad` (台本。サンプルのラティスを先頭に含み、そのまま走る)、`NAME.gs.txt` (gs で走らせた画面表示のリダイレクト)、`NAME.pysad.txt` (pysad で走らせた画面表示のリダイレクト)、`NAME.md` (説明。式と数値の照合、gs と pysad の比較、未確認・限界) です。ファイル名は SAD の機能・関数の名前そのままで、接頭辞は `FFS_` (コマンド)、`fn_` (関数)、`flag_` (フラグ)、`elem_` (要素)、`var_` (特殊変数)、`syntax_` (文法) です。図は同じフォルダの PNG を `.md` から埋め込んでいます。一覧と読み方は [`example/README.md`](example/README.md)。走らせ方: `py -3.12 example/tools/run_example.py example/02_optics_cod/FFS_CALCULATE.sad` (gs は WSL、pysad は `PYSAD_SAD_PACKAGES` つき)。例は既定の v2 モードで走らせた。gs と同じ結果が要るときは `PYSAD_COMPAT=k64`。

| フォルダ | 内容 |
|---|---|
| [`example/00_lattice/`](example/00_lattice/) | サンプルのラティス (`FODO_ring`: 24 セル、3 GeV、周長 105.6 m ほか) |
| [`example/01_lattice_elements/`](example/01_lattice_elements/) | ラティスと要素定義 (要素、LINE、USE、SPLIT、TYPE、ワイルドカード) |
| [`example/02_optics_cod/`](example/02_optics_cod/) | 線形光学と COD (CALCULATE、チューン、クロマティシティ、転送行列、ジオメトリ) |
| [`example/03_matching/`](example/03_matching/) | マッチング (FIT・FREE・FIX・GO、フラグ、OptimizeOptics) |
| [`example/04_emittance_radiation/`](example/04_emittance_radiation/) | エミッタンスと放射 (EMITTANCE、RADINT、RFSW・RAD・FLUC) |
| [`example/05_tracking/`](example/05_tracking/) | トラッキング (TrackParticles、DA、APERT、TRPT、JITTER) |
| [`example/06_collective/`](example/06_collective/) | 集団効果 (WAKE、BEAMBEAM、CSR、Touschek、IBS) |
| [`example/07_spin_solenoid_special/`](example/07_spin_solenoid_special/) | スピン・ソレノイド・特殊要素 (POL、SOL、MULT ANGLE、COORD、RADLIGHT) |
| [`example/08_scripting/`](example/08_scripting/) | SADScript (制御構文、反復子、`;`、リスト・文字列、関数定義、Module・Class、乱数) |
| [`example/09_io_system/`](example/09_io_system/) | 入出力・システム (Print・Write・Read、ファイル、`$FORM`、時間) |
| [`example/10_math_data/`](example/10_math_data/) | 数値・データ処理 (Fit、Spline、NIntegrate、特殊関数、行列、定数) |
| [`example/11_graphics/`](example/11_graphics/) | 描画 (Plot 系、OpticsPlot、BeamPlot、LaTeX) |
| [`example/12_gui_tkinter/`](example/12_gui_tkinter/) | Tkinter・GUI (`Tkinter.n`、`KBFrame.n`) |
| [`example/13_packages/`](example/13_packages/) | `Packages/*.n` (CorrectOrbit、MakeBump、optimize、D.n、TableForm) |
| [`example/pf/`](example/pf/) | PF (2.5 GeV) のラティスでの例集 (光学、マッチング、チューンシフト、放射積分、トラッキング、入射軌道)。索引は [`example/pf/README.md`](example/pf/README.md) |
| [`example/astrid/`](example/astrid/) | ASTRID2 (0.58 GeV) のラティスでの例集 (チューンサーベイ、MPW、位相空間、ラティスの作り方)。索引は [`example/astrid/README.md`](example/astrid/README.md) |
| [`example/v2_B/`](example/v2_B/README.md) | v2.0.0 の例 (特殊関数、文脈、共有メモリ・System・Fork、FFS の小項目)。既定の v2 モードで走る (pysad の出力のみ) |
| [`example/v2_C/`](example/v2_C/README.md) | v2.0.0 の例 (UND の磁場、PHSROT、SOL 領域内、ずれたソレノイド、スピン追跡、RADTAPER、LOSSMAP、RADLIGHT、MEASURE・TRCOD、空間電荷 ほか)。既定の v2 モードで走る (pysad の出力のみ) |
| [`example/v2_next/`](example/v2_next/) | v2.1.0 の例 (放射・F1・CAVI・陽子の追跡・PSPAC、スピン・MPW・RADINT、Emittance の F1・WSPAC・SOL、FREE・複素の特殊関数・Parallelize・文脈・Packages、偏極の側帯・光子・NAFF・幾何の fit、MEASURE・RadiationField・NPARA、NORFSW・光子・POL + 加速、MAIN の式・Micado・CSROYMatrix、ステアリングの F1・B108)。台本・gs の出力 (`.gs.txt`)・pysad の v2 (`.pysad.txt`)・k64 (`.pysad_k64.txt`)。A1・B1・B3・C3・D2 は README (`A1_README.md` ほか) |
| [`example/tools/`](example/tools/) | `.sad` の作成と、gs・pysad での実行 (`mk.py`、`run_example.py` ほか) |

pysad の試験が走らせる台本と、その gs・pysad の画面表示は [`pysad/tests/scripts/`](pysad/tests/scripts/) にあります (対応表は [`pysad/tests/OUTPUTS.md`](pysad/tests/OUTPUTS.md))。`example/` とは別のスクリプトなので、それぞれのフォルダに、それぞれの結果を置いています。

### `docs/` (解説文書。現状のスナップショットで、履歴は書かない)

| 文書 | 何の説明文書か |
|---|---|
| [`01_overview.md`](docs/01_overview.md) | pysad の全体像とモジュールの地図。最初に読む |
| [`02_lattice_main.md`](docs/02_lattice_main.md) | SAD の MAIN レベル (ラティス定義) と `parser.py` |
| [`03_ffs_commands.md`](docs/03_ffs_commands.md) | FFS のコマンドの解釈 (`tffsa.f`) と `ffscmd.py` |
| [`04_sadscript.md`](docs/04_sadscript.md) | SADScript 言語と `sadscript.py` |
| [`05_calc_optics.md`](docs/05_calc_optics.md) | CELL/INS の光学計算 (CALC) と `lattice.py` |
| [`06_matching.md`](docs/06_matching.md) | マッチング (FIT/FREE/GO)、解析微分、幾何の条件 |
| [`07_emittance.md`](docs/07_emittance.md) | `Emittance[]` と EMIT (`emit6d.py`) |
| [`08_tracking.md`](docs/08_tracking.md) | `TrackParticles` の粒子追跡 (`track6.py`・`track6n.py`) |
| [`09_dynamic_aperture.md`](docs/09_dynamic_aperture.md) | `DynamicApertureSurvey` (`dapert.py`) |
| [`10_benchmarks.md`](docs/10_benchmarks.md) | SAD (K64) との照合ベンチマーク一覧と回帰・試験の現状 |
| [`11_sad_quirks.md`](docs/11_sad_quirks.md) | SAD の癖・仕様・ビルド依存の挙動と pysad での扱い |
| [`12_pysad_api.md`](docs/12_pysad_api.md) | pysad の公開 API のリファレンス |
| [`13_sad_help_index.md`](docs/13_sad_help_index.md) | SADHelp (CERN 版) から機械的に抽出したコマンド・関数・変数の一覧 |
| [`14_pysad_vs_sadhelp.md`](docs/14_pysad_vs_sadhelp.md) | SADHelp の項目と pysad の対応表 |
| [`15_sadhelp_purpose_policy.md`](docs/15_sadhelp_purpose_policy.md) | SADHelp の項目の目的別分類と翻訳方針 |
| [`16_fortran_status.md`](docs/16_fortran_status.md) | SAD の Fortran ファイル別の翻訳状況と未翻訳の一覧 |
| [`17_implementation_summary.md`](docs/17_implementation_summary.md) | 実装のまとめ (モジュール表、準拠状況、未対応) |
| [`18_physics_guide.md`](docs/18_physics_guide.md) | SAD の物理計算ガイド (式と K64 の実装) |
| [`19_source_analysis.md`](docs/19_source_analysis.md) | SAD ソースの解析と pysad の対応表 |
| [`20_translation_candidates.md`](docs/20_translation_candidates.md) | pysad の翻訳候補と状況 (SAD の機能別、済み・見送りと理由) |
| [`21_sad_bug_report.md`](docs/21_sad_bug_report.md) | pysad の移植中に見つかった SAD K64 のバグ・癖の報告 (B 番号) |
| [`22_gui_tclbridge.md`](docs/22_gui_tclbridge.md) | Tcl/Tk の GUI (`Tkinter.n`・`KBFrame.n`) を pysad で動かす仕組み |
| [`23_wake.md`](docs/23_wake.md) | ウェイク場 (WakeFunction、LWAKE/TWAKE、WAKEOPT) の現状 |
| [`24_beambeam.md`](docs/24_beambeam.md) | BEAMBEAM (強弱ビームのビーム・ビーム相互作用) の現状 |
| [`25_coord.md`](docs/25_coord.md) | COORD (座標変換) の現状 |
| [`26_tcav.md`](docs/26_tcav.md) | TCAVI (横偏向空洞) の現状 |
| [`27_und.md`](docs/27_und.md) | UND (アンジュレータ) の現状 |
| [`28_mult_angle.md`](docs/28_mult_angle.md) | MULT (ANGLE≠0、曲げつき多極) の現状 |
| [`29_sol_optics.md`](docs/29_sol_optics.md) | SOL 領域の Twiss・Emittance・EMIT の現状 |
| [`30_map.md`](docs/30_map.md) | MAP (ExternalMap) の現状 |
| [`31_delcor.md`](docs/31_delcor.md) | DELCOR (相関のあるずれ) の現状 |
| [`32_gaussiancoulomb_bbbrem.md`](docs/32_gaussiancoulomb_bbbrem.md) | GaussianCoulomb と BBBrem1 の現状 |
| [`33_specfn_fit_misc.md`](docs/33_specfn_fit_misc.md) | 特殊関数、Fit・FindRoot・NIntegrate、DownhillSimplex などの現状 |
| [`34_sad_packages.md`](docs/34_sad_packages.md) | SAD の `Packages/*.n` を SADScript のまま使う (軌道補正・最適化・数学) |
| [`35_ffs_small_commands.md`](docs/35_ffs_small_commands.md) | FFS の小さなコマンド (`ffsextra.py`) |
| [`36_spin_csr.md`](docs/36_spin_csr.md) | POL (スピン) と CSR の現状 |
| [`37_plot_canvas.md`](docs/37_plot_canvas.md) | 描画 (ListPlot・Plot と Tk のキャンバス) の現状 |
| [`38_wake_guide.md`](docs/38_wake_guide.md) | 解説: WAKE (使い方、式、K64 の癖、限界) |
| [`39_beambeam_guide.md`](docs/39_beambeam_guide.md) | 解説: BEAM-BEAM (BEAMBEAM、BSTRL、光子、ルミノシティ、BBBrem) |
| [`40_csr_guide.md`](docs/40_csr_guide.md) | 解説: CSR (`CSR.n` と組み込み関数の使い方、物理、限界) |
| [`41_tcltk_guide.md`](docs/41_tcltk_guide.md) | 解説: Tcl/Tk の GUI (部品、使い方、ヘッドレス試験) |
| [`42_graphics_guide.md`](docs/42_graphics_guide.md) | 解説: グラフィック描画 (Plot 系、OpticsPlot、Canvas、gs との比較) |
| [`43_pol_guide.md`](docs/43_pol_guide.md) | 解説: POL (スピンと偏極、TrackParticles のスピン追跡) |
| [`44_sol_guide.md`](docs/44_sol_guide.md) | 解説: ソレノイド (SOL 領域の仕組みと対応範囲) |
| [`45_radlight_guide.md`](docs/45_radlight_guide.md) | 解説: RADLIGHT と放射場 (RadiationField・RadiationSpectrum) |
| [`46_unimplemented_and_unverified.md`](docs/46_unimplemented_and_unverified.md) | オリジナル SAD (K64) に対する pysad の未対応・未確認の一覧 (未翻訳、一部未対応、未確認、既知の差、gs 側の問題、優先度) |
| [`47_v2_changes.md`](docs/47_v2_changes.md) | pysad v2 の変更の一覧 (v2.0.0 と v2.1.0。v1 = K64 に忠実な翻訳との違い。互換モード、v2 で結果が変わるもの、新しい機能、pysad の不具合の修正、要判断、試験の状態) |
| [`sadhelp_summary.html`](docs/sadhelp_summary.html) | SADHelp 対応文書 (13・14・15) の要約 (リポジトリ外の `perf_work/sadhelp_html.py` で md から生成) |
| [`html/`](docs/html/) | 上の md を html にしたもの (`make_html.py` で生成。追跡しない)。手書きの `sad_wake_sigsegv_fixes_v001.html` (2026-11-14) と `_v002.html` (2026-11-16、B33 の差の原因と解消を追記) は K64 の wake まわりの SIGSEGV 3 件 (B31〜B33) の原因と修正の報告。差分は `docs/history/sad_patch/wake_sigsegv_fixes.patch` |

### `docs/history/` (作業履歴)

| 文書 | 内容 |
|---|---|
| [`00_RESUME.md`](docs/history/00_RESUME.md) | 引き継ぎメモ (方針・環境・現状・残課題・手順・落とし穴・履歴の索引)。最初に読む |
| [`work_unimpl.md`](docs/history/work_unimpl.md) | 11/14 以降の未実装要素の実装と、未解決の差 (最新) |
| [`work_session_20261116.md`](docs/history/work_session_20261116.md) | 2026-11-16 の作業 (POL のスピン追跡、RADLIGHT など) を話題別に |
| [`review_notes.md`](docs/history/review_notes.md) | docs のレビューメモ (2026-10-15) |
| [`work_wake.md`](docs/history/work_wake.md)、[`work_beambeam.md`](docs/history/work_beambeam.md)、[`work_gui_tclbridge.md`](docs/history/work_gui_tclbridge.md) | WAKE、BEAMBEAM、GUI の実装の経緯 (11/12〜11/14) |
| [`work_analytic_deriv.md`](docs/history/work_analytic_deriv.md)、[`work_descend.md`](docs/history/work_descend.md)、[`work_k0fit.md`](docs/history/work_k0fit.md)、[`work_match_fixes.md`](docs/history/work_match_fixes.md) | 解析微分、最急降下、K0 の fit 変数、match6 の修正 |
| [`work_calc6d.md`](docs/history/work_calc6d.md)、[`work_emit_extra.md`](docs/history/work_emit_extra.md)、[`work_emit_report.md`](docs/history/work_emit_report.md)、[`work_emit_unstable.md`](docs/history/work_emit_unstable.md) | CALC6D、Emittance の戻り規則、EMIT 報告、不安定な光学 |
| [`work_geo.md`](docs/history/work_geo.md)、[`work_geo_sig.md`](docs/history/work_geo_sig.md)、[`work_ins_offset.md`](docs/history/work_ins_offset.md)、[`work_bend_k0.md`](docs/history/work_bend_k0.md)、[`work_rotate.md`](docs/history/work_rotate.md) | 幾何、幾何の fit とビームサイズ、INS の初期軌道、BEND の K0、ROTATE |
| [`work_touschek.md`](docs/history/work_touschek.md)、[`work_track6n_ext.md`](docs/history/work_track6n_ext.md)、[`work_loss_rng.md`](docs/history/work_loss_rng.md)、[`work_loss_rng2.md`](docs/history/work_loss_rng2.md)、[`work_zdp.md`](docs/history/work_zdp.md)、[`work_trpt_acc.md`](docs/history/work_trpt_acc.md)、[`work_main_track.md`](docs/history/work_main_track.md) | Touschek、numba の追跡の拡張、粒子損失の乱数、z/dp、TRPT の加速、MAIN の TRACK |
| [`work_integration.md`](docs/history/work_integration.md)、[`work_item4.md`](docs/history/work_item4.md)、[`work_item6.md`](docs/history/work_item6.md)、[`work_small_items.md`](docs/history/work_small_items.md)、[`work_unstable.md`](docs/history/work_unstable.md) | モジュールの接続、SADScript と FFS の拡張、disp/emit の出力、小項目、不安定ライン |
| [`v2_changes_A.md`](docs/history/v2_changes_A.md)、[`B`](docs/history/v2_changes_B.md)、[`C`](docs/history/v2_changes_C.md)、[`C_measure`](docs/history/v2_changes_C_measure.md)、[`C_spch`](docs/history/v2_changes_C_spch.md)、[`C2`](docs/history/v2_changes_C2.md)、[`C2_measure`](docs/history/v2_changes_C2_measure.md)、[`C2_misc`](docs/history/v2_changes_C2_misc.md)、[`C2_spch`](docs/history/v2_changes_C2_spch.md) | v2.0.0 の作業の記録 (根拠の数値、試験、使える範囲、要判断)。まとめは [`docs/47_v2_changes.md`](docs/47_v2_changes.md) |
| [`v2next_A1.md`](docs/history/v2next_A1.md)、[`A2`](docs/history/v2next_A2.md)、[`B1`](docs/history/v2next_B1.md)、[`B2`](docs/history/v2next_B2.md)、[`B3`](docs/history/v2next_B3.md)、[`C`](docs/history/v2next_C.md)、[`D1`](docs/history/v2next_D1.md)、[`D2`](docs/history/v2next_D2.md)、[`E`](docs/history/v2next_E.md) | v2.1.0 の作業の記録 (根拠の数値、試験、使える範囲、要判断)。まとめは [`docs/47_v2_changes.md`](docs/47_v2_changes.md) |
| [`work_perf2.md`](docs/history/work_perf2.md)、[`work_perf3.md`](docs/history/work_perf3.md)、[`work_performance.md`](docs/history/work_performance.md) | 性能改善 (numba カーネル化) |

## ライセンス

`LICENSE` は SAD の利用条件 (Copyright (c) 1986-2016, KEK Accelerator Laboratory。https://acc-physics.kek.jp/SAD/how-to-get-sources/licensing/ と K64 版の LICENSE に同じ文面) をそのまま置いたもので、文面は BSD 3 条項ライセンスと同じ。pysad は SAD の派生物として、この著作権表示と条件をそのまま適用する。
