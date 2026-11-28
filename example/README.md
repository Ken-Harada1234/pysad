# example: SAD の機能・関数ごとの短い例 (gs と pysad の比較つき)

SAD の機能・関数を 1 つずつ短い例にしたものです。1 つの例は 4 つのファイルで 1 組です。

| ファイル | 内容 |
|---|---|
| `NAME.sad` | SAD スクリプト (サンプルのラティスを先頭に含み、そのまま走る) |
| `NAME.gs.txt` | SAD (gs) で走らせたときの画面表示のリダイレクト |
| `NAME.pysad.txt` | pysad で走らせたときの画面表示のリダイレクト |
| `NAME.md` | 説明 (何を見る例か、結果の読み方、物理・解析式と数値の照合、gs と pysad の比較、未確認・限界) |

ファイル名は SAD の機能・関数そのものの名前です。種別の接頭辞は `FFS_` (FFS コマンド)、`fn_` (関数)、`flag_` (フラグ)、`elem_` (要素)、`var_` (特殊変数)、`syntax_` (文法) です。図は同じフォルダの PNG を `.md` から埋め込んでいます。式は LaTeX (`$...$`) で書いてあり、GitHub で表示されます。

## フォルダ

| フォルダ | 分類 | 内容 |
|---|---|---|
| [`00_lattice`](00_lattice) | サンプルのラティス | `FODO_ring` (24 セル、3 GeV、周長 105.6 m) ほか |
| [`01_lattice_elements`](01_lattice_elements) | ラティスと要素定義 | 要素、LINE、USE、SPLIT、TYPE、ORG、ワイルドカードなど |
| [`02_optics_cod`](02_optics_cod) | 線形光学と COD | CALCULATE、チューン、クロマティシティ、転送行列、COD、ジオメトリ |
| [`03_matching`](03_matching) | マッチング | FIT・FREE・FIX・GO、フラグ、OptimizeOptics |
| [`04_emittance_radiation`](04_emittance_radiation) | エミッタンスと放射 | EMITTANCE、RADINT、RFSW・RAD・FLUC、放射積分の式 |
| [`05_tracking`](05_tracking) | トラッキング | TrackParticles、DA、APERT、TRPT、JITTER |
| [`06_collective`](06_collective) | 集団効果 | WAKE、BEAMBEAM、CSR、Touschek、IBS |
| [`07_spin_solenoid_special`](07_spin_solenoid_special) | スピン・ソレノイド・特殊 | POL、SOL、MULT ANGLE、COORD、RADLIGHT |
| [`08_scripting`](08_scripting) | SADScript | 制御構文、反復子、`;`、リスト・文字列、関数定義、Module・Class、乱数 |
| [`09_io_system`](09_io_system) | 入出力・システム | Print・Write・Read、ファイル、`$FORM`、時間 |
| [`10_math_data`](10_math_data) | 数値・データ処理 | Fit、Spline、NIntegrate、特殊関数、行列、定数 |
| [`11_graphics`](11_graphics) | 描画 | Plot 系、OpticsPlot、BeamPlot、LaTeX |
| [`12_gui_tkinter`](12_gui_tkinter) | Tkinter・GUI | Tkinter.n、KBFrame.n |
| [`13_packages`](13_packages) | Packages/*.n | CorrectOrbit、MakeBump、optimize、D.n、TableForm など |
| [`14_nonlinear_dynamics`](14_nonlinear_dynamics/README.md) | 非線形ダイナミクス | 周波数マップ (`frequency_map/`)、6 極の共鳴駆動項 (RDT)・振幅依存チューンシフト (ADTS)・クロマティシティ (`adts_rdt/`)、1 周の写像の正規形・K0≠0 の偏向磁石・共鳴の島 (`normal_form/`)、入門の図 (`primer/`)。解説は docs/24〜26。台本ごとのまとめは `adts_rdt/ADTS_summary.md`・`frequency_map/FMAP_summary.md` |
| [`15_lma_touschek`](15_lma_touschek/README.md) | 局所運動量アクセプタンス (LMA) と Touschek 寿命 | LMA の計算・信頼性・誤差と補正 (`lma/`)、elegant との LMA・Touschek 寿命の比較と Touschek の式 (`touschek_elegant/`)。解説は docs/26 |
| [`16_spin_polarization`](16_spin_polarization/README.md) | スピン追跡・偏極度・緩和率 | POL・Polarization・側帯・クロマティックな変調・領域内のスピン・加速中のスピン。解説は docs/16 |
| [`17_notebook`](17_notebook/README.md) | Jupyter の notebook から実行 | ラティスとコマンドをセルに貼って実行する `%%sad`・結果を Python で使う `sad.value()`・`pysad.run()`。出力つきの notebook `sad_in_notebook.ipynb` |
| [`pf`](pf) | PF (2.5 GeV) の例集 | 実機ラティスでの光学・マッチング・チューンシフト・放射積分・トラッキング・入射軌道。PF の履歴ラティス (1981・1986・1997・2005 年、[`pf/lattice/PF_history.md`](pf/lattice/PF_history.md)) の例を含む。索引は [`pf/README.md`](pf/README.md) |
| [`kek_ls`](kek_ls) | KEK-LS (3 GeV、20 セル、ESRF 型 HMBA) の例集 | 設計 (セルのチューン・6 極・入射セル・リングを閉じる・縮小・セル数・チューンサーベイ)、パラメータ表示 (εx・ξ・オプティクス・チューンシフト)、誤差入りの COD の補正。元の設計用スクリプトを整えた一式。索引は [`kek_ls/README.md`](kek_ls/README.md)、元との対応は [`kek_ls/ORIGINAL_MAP.md`](kek_ls/ORIGINAL_MAP.md) |
| [`ddba`](ddba) | 複雑なセルのラティスの設計手順の例 (DDBA、HLS ring の系列。5 GeV・2.5 GeV、周長 731〜757 m) | 設計メモの版 `DDBA_01`〜`DDBA_28` を、基本セルの構想・セルのマッチング・チューン点・6 極・等時性セル・リングへの展開・誤差と補正・最終の版、の段階に並べ直したもの。28 本のラティス (`cell calc emit` の周長・チューン・εx の表、gs との比較)、最終の版の例。最終の版 v07_14 のチューンのマッチング・チューンサーベイ・誤差の計算 (補正しながら誤差を大きくする手法)・パラメータ表示の一式は [`ddba/v07_14/`](ddba/v07_14/README.md)、6 極の検討 (RDT・ADTS の解析式による検証と最適化の試作。最適化は最終設定より DA の広い解を出せなかった) は [`ddba/v07_14/sextupole/`](ddba/v07_14/sextupole/README.md) (台本 `V0714_*.sad`、gs と pysad の比較、元のスクリプトとの対応 `ORIGINAL_MAP.md`)。索引は [`ddba/README.md`](ddba/README.md)、元との対応は [`ddba/ORIGINAL_MAP.md`](ddba/ORIGINAL_MAP.md) |
| [`astrid`](astrid) | ASTRID2 (0.58 GeV) の例集 | チューンサーベイ、MPW、位相空間、ラティスの作り方。索引は [`astrid/README.md`](astrid/README.md) |
| [`dxf_lattice`](dxf_lattice/README.md) | ラティスの CAD 図 (DXF) | `DrawMagnetsDXF[]` の平面図・磁石の立体。解説は docs/22 |
| [`_deletable`](_deletable/README.md) | 削除してよい作業ファイルの置き場 | 整理のあとに不要になった途中の出力。例の説明・試験は参照していない |
| `tools` | 道具 | `.sad` の作成と、gs・pysad での実行 (`mk.py`、`run_example.py` ほか。`pf_astrid/` は pf・astrid の作成用で、パスは整理前のまま) |

各フォルダの `README.md` に、そのフォルダの例の一覧と、gs と pysad の一致状況 (○ 一致、△ 差あり、× 未対応・失敗) があります。

## 走らせ方

```bash
py -3.12 example/tools/run_example.py example/02_optics_cod/FFS_CALCULATE.sad
```

gs は WSL で、pysad は `PYSAD_SAD_PACKAGES` (K64 版 SAD の `Packages` のフォルダ) つきで走ります。gs の場所は `PYSAD_SAD_PACKAGES` の隣の `bin/gs` で、別の場所のときは `PYSAD_GS` (WSL から見た絶対パス) で指定します ([`example/tools/sadenv.py`](tools/sadenv.py)。gs を使う `example/` のスクリプトに共通)。`NAME.gs.txt` と `NAME.pysad.txt` が作り直されます。

**v2 の注記**: 例は既定の v2 モードで走らせた (`NAME.pysad.txt` は v2 モードの出力)。gs と同じ結果が要るときは `PYSAD_COMPAT=k64` を付けて走らせる (K64 のバグも再現する)。v2 で K64 のバグを直したために gs と違う値になる例は、各 `.md` の「gs と pysad の比較」に「v2 の注記」として書いてある (v2 の変更の一覧は [`docs/28_v2_vs_k64.md`](../docs/28_v2_vs_k64.md))。主な理由は、`Emittance[]`・EMIT・追跡の Jx の端の項 (B75)、`MomentumCompaction` (B77)、`PolynomialFit` の x の平均が 0 のときの係数 (B18)、`InverseErf` (B47)、`For` の本体 (B71)、`RADLIGHT` の BEND の中の記録 (B65)、SPAC の空間電荷 (SPC1〜) です。各フォルダの `README.md` の ○△× は、v2 で違うものも含めて、`PYSAD_COMPAT=k64` では gs と一致するものを ○ としています (違いが残るものは行に書きました)。例外として `pf`・`astrid` のうち ASTRID2 の 5 例 (`MPW_v01_02_id`、`SetElement_BeamLine_Join`、`phase_x`、`amp_dep_ts`、`calc_RadInt`) は、v2 の既定との違いが表示の 7 桁目 (バケット高さ、αc) だけなので `.pysad.txt` を作り直していません。

## 読むときの注意

- 表示や数字の桁、テキストの微妙なずれは追っていません。基準は、完走することと、意味のある値が合うことです。
- gs と pysad の差が大きいものは、各 `.md` の「gs と pysad の比較」と「未確認・限界」に書いてあります。主な pysad の未対応と差は、`For`・`Module` のスコープ・`BeginPackage`・`Read` の区切り・APERT の追跡・RADTAPER・SPAC・LOSSMAP・Tkinter ウィジェットのメソッドなどです。
- gs が落ちる例 (RadiationField の SIGSEGV など) は、その旨を `.gs.txt` と `.md` に残しています。
- gs の窓の画像は撮れていません。描画と GUI の図は、Canvas の項目を書き出した再描画と、pysad の窓の切り取りです。
- SAD は `E1`・`E2` を ANGLE に対する比として読みます (端の角度 = E1×ANGLE)。`FODO_ring` の BEND はこの読みで解析式と合います。
- ラティスはサンプルで、特定の実機のものではありません。エネルギー 3 GeV、周長 105.6 m、RF 電圧は例ごとに `.md` に書いてあります。
- `pf` と `astrid` は、PF と ASTRID2 の実際のラティスを使った例集です。KEK の SAD チュートリアルの話題を、手元のラティスで作り直しました。
