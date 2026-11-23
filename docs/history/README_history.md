# README.md の履歴 (版ごとの注記・経過の記述)

README.md を現在の機能・使い方・構成・検証結果の説明に直したとき、書き換えた行の変更前の全文です。含まれていた節の見出しごとに、元の文脈のままの 1 行を入れています。


## 元の節: pysad — SAD (K64 版) の Python 直訳

> **v1 と v2 の違い (注釈)**: 以下の「直訳」は **v1 系 (v1.0.0) の方針**です。**v2 系 (v2.0.0 以降) は直訳ではありません。** v2 は K64 を出発点に、K64 のバグの修正、物理の誤りの修正、機能の追加を行った版で、**gs (K64) と結果が変わる箇所があります** ([`docs/47_v2_changes.md`](docs/47_v2_changes.md))。gs と同じ結果が要るときは、環境変数 `PYSAD_COMPAT=k64` (互換モード。v1 と同じ K64 の動作) を指定してください。なお、v2 でも、K64 の式・分割数・反復順序を写したコードの大部分は、そのままです。


## 元の節: インストール (PyPI)

PyPI: <https://pypi.org/project/pysad-k64/>。import 名は `pysad` です (パッケージ名と違います)。**v2.3.0** (K64 のバグの修正と機能の追加。互換モード `PYSAD_COMPAT=k64` で K64 と同じ動作。v2 は直訳ではありません。GitHub のリリース v2.3.0 と同じ) が最新です。ソースに忠実な翻訳の **v1.0.0** は `pip install pysad-k64==1.0.0` で入ります。


## 元の節: 必要なもの

確認した Python の版 (v2.0.0、全試験): **3.12** と **3.10** は全通過 (937 通過・4 skip)。**3.14** は 935 通過・2 失敗 (3.14 に同梱の Tcl/Tk が 9.0 で、試験が 8.6 を期待していること、数値的に脆い 1 試験の閉軌道の収束の違い)。v2.1.0 では 2 件とも直し、3.14 で GUI 系と D2 の試験 (計 115 件) が通ることを確かめた (全試験は 3.14 では走らせていない。詳細は [`docs/46_unimplemented_and_unverified.md`](docs/46_unimplemented_and_unverified.md) §3)。3.13 は未確認。


## 元の節: オリジナルの SAD にない機能 (pysad で追加)

pysad v2 系は直訳に加えて、**オリジナルの SAD (K64) にない解析機能**を持ちます。入門は [`docs/50_adts_rdt_primer.md`](docs/50_adts_rdt_primer.md)、台本ごとの数値のまとめは [`ADTS_summary.md`](example/v2_3/ADTS_summary.md)・[`FMAP_summary.md`](example/v2_3/FMAP_summary.md)、例の索引は [`example/v2_3/README.md`](example/v2_3/README.md) です。

| 機能 | 何を求めるか | 使い方 (SAD の台本 / Python) | 解説・例 | 版 |

|---|---|---|---|---|

| 周波数マップ | 初期条件の格子ごとのチューンの拡散 D の地図と共鳴線 | `FrequencyMap[]` / `pysad.freqmap` (NAFF ベース) | [docs/48](docs/48_frequency_map_lma.md)、[FMAP_summary](example/v2_3/FMAP_summary.md) | v2.3.0 |

| 局所運動量アクセプタンス (LMA) と Touschek 寿命 | 位置ごとの δ+(s)・δ−(s) (RF バケットとアパーチャを区別。`check_every` の既定は 1) と、それを使った Touschek 寿命 | `MomentumAcceptance[]`・`TouschekLifetimeLMA[]` / `pysad.lma`・`pysad.lmabucket` | [docs/48](docs/48_frequency_map_lma.md)、[L2](example/v2_3/L2_README.md)、[T1](example/v2_3/T1_README.md) | v2.3.0 |

| Touschek 寿命の独立な式 (Piwinski 型) | 結合・丸いビームまで正しい損失率。K64 の式との差は結合 1% で +0.1〜+6%、丸いビームで ±20% (B162)。v2 の `TouschekLifetime[INF, INF, nz]`・`TouschekLifetimeLMA[]` の既定 (`Formula->"K64"` で K64 の表) | `touschek.touschek_piwinski` | [docs/48 §4.6](docs/48_frequency_map_lma.md)、[T2](example/v2_3/T2_README.md) | v2.3.0 |

| 共鳴駆動項 (RDT) | 6 極 (1 次・2 次の 4 次の h_jklm)・8 極・偏向磁石の項。厚いレンズと薄レンズの両方 | `ResonanceDrivingTerms[]` / `pysad.resdrive`・`pysad.rdt2nd` | [docs/49](docs/49_resonance_terms_adts.md)、[ADTS_summary](example/v2_3/ADTS_summary.md) | v2.3.0 |

| 振幅依存チューンシフト (ADTS) | 6 極の 2 次・8 極・運動学・フリンジ・偏向磁石の項の ADTS (既定は物理的な項だけ。精密版は `fourth='r6'`・`optics='tracked'`) | `AmplitudeTuneShift[]` / `pysad.resdrive.adts` | [docs/49](docs/49_resonance_terms_adts.md)、[ADTS_summary](example/v2_3/ADTS_summary.md) | v2.3.0 |

| クロマティシティの項 | 自然クロマティシティ、6 極と偏向磁石の寄与 | `ChromaticityTerms[]` / `pysad.resdrive` | [docs/49](docs/49_resonance_terms_adts.md) | v2.3.0 |

| 数値の正規形 | 追跡の 1 周の写像から ADTS を取り出す | `pysad.nform` | [docs/49 §10](docs/49_resonance_terms_adts.md) | v2.3.0 |

| 振幅・δ に依る ν のべき級数 (J²・δ 依存・交差項) | ν(Jx, Jy, δ) を J・δ のべき級数 Σ c_abp Jx^a Jy^b δ^p として (J¹ = ADTS、J²、ξ の高次、振幅と δ の交差項)。模型の 1 周の写像を TPSA で作り Birkhoff の正規形で求める。線形 + 6 極・8 極に、オプションで偏向磁石の曲率の項・運動学的な項・端部のフリンジ (FODO・ASTRID2・大きなリング A で追跡と 1e-6 前後。入れていない効果は docs/49 §13.3) | `NormalFormTune[]` / `pysad.ampnf` | [docs/49 §13](docs/49_resonance_terms_adts.md)、[example/v2_4/](example/v2_4/README.md) | v2.4.0 |

| ν(J) の J の次数・5・6 次の駆動項・1 本の共鳴の島 | `ampnf` の J² (D ≥ 5)・J³ (D ≥ 7) を J の次数ごとに追跡と比較 (`nu(jorder=)`・`d2nu`)。純粋な 3〜6 次の共鳴の駆動項の表 (`drive_table`)。1 本の共鳴を残した正規形から固定点の作用と島の半幅 (`island_theory`。6 次で追跡と 0.04〜1.5%・島の端 0.3〜1.3%、5 次以下は島が広く外れる。対称性で禁止された共鳴にも有限の幅を返す)。K0≠0 の磁石の非線形の切り替え `k0_nonlinear` | `pysad.ampnf_ext`・`pysad.ampnf` | [docs/49 §15](docs/49_resonance_terms_adts.md)、[example/v2_4/](example/v2_4/README.md) | v2.4.1 |

| Touschek のバンチ内の δ の分布 (オプション) | 散乱される粒子の δ の分布 (対の重心の δ_c と相対運動量の縦成分 q_s) を入れた損失率。既定は入れない (K64・elegant・Piwinski と同じ δ=0 の基準)。ASTRID2 で受容が σδ の 9 倍のとき損失率が +4%。対のモンテカルロと 0.3% 以内 | `TouschekLifetimeLMA[BunchDelta->True, BunchQs->True]` / `pysad.touschek_bunch` | [docs/48 §4.4](docs/48_frequency_map_lma.md)、[example/v2_4/](example/v2_4/README.md) | v2.4.0 |

| `DynamicApertureSurvey` の `ExpandElementValues`・`ReferenceOrbit` | `ExpandElementValues->False` で成分ごとの値 (補正磁石の K0・位置ずれ) を揃えない。既定 (True) は gs と同じく揃えるので、補正した状態の DA は `False` を書く (DDBA v07_14 の誤差の例で gs と pysad が 29〜36σ で一致)。`ReferenceOrbit->{x0,px0,y0,py0,z0,dp0}` はその軌道のまわりで追跡 | `DynamicApertureSurvey[…, ExpandElementValues->False]` / `pysad.dapert` | [docs/09](docs/09_dynamic_aperture.md)、[`docs/history/v24_A.md`](docs/history/v24_A.md) | v2.4.0 |

| FFS の `DRAW` (matplotlib) | `DRAW [r1 [r2]] BX BY & EX LAT;` を Tk でなく matplotlib に描く (領域・`&` で枠を分ける・ラティス図と名前・クリックの InfoLabel)。保存は `PYSAD_DRAW_FILE` か `Draw$Option`。`OpticsPlot` の Graphics と点列が一致 (gs の Graphics と 1e-13 台)。matplotlib は任意 (`pip install pysad-k64[plot]`) | `DRAW …;` / `pysad.draw_mpl.draw` | [docs/42 §7b](docs/42_graphics_guide.md)、[example/11_graphics/fn_DRAW](example/11_graphics/fn_DRAW.md) | v2.4.1 (開発中) |

| `$DisplayFunction` の matplotlib 出口 (`MatplotlibDrawer`) | `Plot`・`ListPlot`・`OpticsPlot`・`GeometryPlot` などが作る `Graphics[…]` を Tk でなく matplotlib に描く (色・線の太さ・破線・点の記号・誤差棒・文字の位置揃え・`PlotRange`・`Scale`・`Frame`・`FrameTicks`・`GridLines`・`FrameLabel`・`PlotLabel`・凡例)。保存は `PYSAD_DISPLAY_FILE` (2 回目以降は番号つき)。未対応の要素は警告して飛ばす。Graphics の要素数・座標が matplotlib 側のデータと一致 (試験)。見た目は Tk の窓と比べていない。matplotlib は任意 | `PYSAD_DISPLAY=mpl` / `$DisplayFunction=MatplotlibDrawer;` / `pysad.mpl_display.render` | [docs/42 §7c](docs/42_graphics_guide.md)、[example/11_graphics/fn_OpticsPlot_mpl](example/11_graphics/fn_OpticsPlot_mpl.md) | v2.4.1 (開発中) |

| ラティスの CAD 図 (DXF) | ラティスの平面図 (基準軌道・偏向磁石・四極・6 極などの外形・MARK・要素名) を、テキストの DXF (R12 相当、外部ライブラリなし) に書き出す。レイヤーと色・`Region`・単位 (m・mm)・`Spec` (磁石の寸法)。gs でも pysad でも走る台本の関数と Python 関数の出力が一致 (gs と 8 本全行一致)。立体は未対応 | `DrawMagnetsDXF[]` / `pysad.dxf_lattice.write_dxf` | [docs/51](docs/51_lattice_dxf.md)、[example/dxf_lattice](example/dxf_lattice/README.md) | v2.4.1 (開発中) |

| 6 極の検討 (試作) | 6 極の K2L の最適化の試作と、解析式による検証の可否 (結論: 検証には条件つきで使える。最適化は最終設定より DA の広い解を出せなかった) | `pysad.sxtopt` | [example/ddba/v07_14/sextupole](example/ddba/v07_14/sextupole/README.md) | v2.3.0 (試作) |

| 6 極の最適化の項の指定 | 項ごとに objective/constraint/ignore・重み・目標値を指定して K2L を計算 (`SxtSpec`・`design_k2l`)、誤差感度の項、作動点の周りの共鳴の一覧。HLS ring v07_14 の例: K2L のランダム誤差 (3%・10%) の DA で評価 (ADTS だけの最適化は完成版の 0.6 倍、h00201 を 0 にした候補は誤差なしで 1.30 倍・3% で 1.17 倍、10% で同じ。gs で確認、LMA は未了) | `pysad.sxtopt` (`design_k2l`・`near_resonances`) | [example/ddba/v07_14/sextupole_v2](example/ddba/v07_14/sextupole_v2/README.md)、[docs/47](docs/47_v2_changes.md) | v2.4.1 (開発中) |

| 6 極の最適化の項の指定 | 項ごとに objective/constraint/ignore・重み・目標値を指定して K2L を計算 (`SxtSpec`・`design_k2l`)、誤差感度の項、作動点の周りの共鳴の一覧。HLS ring v07_14 の例: K2L のランダム誤差 (3%・10%) の DA で評価 (ADTS だけの最適化は完成版の 0.6 倍、h00201 を 0 にした候補は誤差なしで 1.30 倍・3% で 1.17 倍、10% で同じ。gs で確認、LMA は未了) | `pysad.sxtopt` (`design_k2l`・`near_resonances`) | [example/ddba/v07_14/sextupole_v2](example/ddba/v07_14/sextupole_v2/README.md)、[docs/47](docs/47_v2_changes.md) |

**v2 で K64 から改良・追加した主な点** (v2.0.0〜v2.4.0。一覧と条件は [`docs/47_v2_changes.md`](docs/47_v2_changes.md)):

- v2.4.0: パーサ・評価器・FIT の修正 (ユーザー定義の関数名 `At`、`FitFunction` を併用する `fit`、`ElementValues` とフィット変数の連動など)、放射・エミッタンス (tcod が収束しないときの Emittance[] (B163)、RADINT の ANGLE=0・K0≠0 の極、v2 の EMITDIV の既定を 4 に)、`DynamicApertureSurvey` の `ExpandElementValues`・`ReferenceOrbit`、`ampnf` (J²・δ 依存)。変更の一覧は [`CHANGELOG.md`](CHANGELOG.md) と [`docs/47_v2_changes.md`](docs/47_v2_changes.md)、朝に判断する項目の表は docs/47 の末尾。

- v2.3.0 の新しい例: PF の履歴ラティス ([`example/pf/lattice/PF_history.md`](example/pf/lattice/PF_history.md))、KEK-LS の例集 ([`example/kek_ls/`](example/kek_ls/README.md)。参考: KEK-LS Conceptual Design Report (CDR)、KEK、2017 年 5 月 22 日)、DDBA (HLS ring の系列) の複雑なセルの設計手順の例集 ([`example/ddba/`](example/ddba/README.md)・最終版 [`v07_14/`](example/ddba/v07_14/README.md)・6 極の検討 (試作) [`v07_14/sextupole/`](example/ddba/v07_14/sextupole/README.md)。参考: K. Harada et al., *Conceptual design of the Hybrid Ring with superconducting linac*, J. Synchrotron Rad. 29, 118–124 (2022)、DOI 10.1107/S1600577521012753)、ASTRID2 の例集 ([`example/astrid/`](example/astrid/README.md)。参考: *ASTRID2 Preliminary Design Report*、Aarhus University、ISA、2009 年 8 月)。

- v2.3.0 のパーサの修正 (pysad 自身の不具合): 逆向きの文脈の `-LINE` の二重反転、k64 モードの `1.7MV` などの数に隣接する単位 (gs は読む)。

- これらは **pysad のみ**で、gs (SAD 本体) では走りません。v2.3.0 の新機能 (上の表) は互換モード (`PYSAD_COMPAT=k64`) によらず同じ動作です。v2.0.0〜v2.2.0 で K64 の動作を変えたもの (互換モードで切り替わる) は docs/47 の (a) に区別して書いてあります。

- **PyPI は v2.3.0 まで公開**しています。v2.3.0 の新機能も `pip install pysad-k64` で入ります。SAD の台本を `py -3.12 run_sad.py lattice.sad script.sad` で動かすには、リポジトリを clone してください (PyPI のパッケージには `run_sad.py`・`example/`・`docs/` が入りません)。

- 検証の範囲と未確認は [`docs/46_unimplemented_and_unverified.md`](docs/46_unimplemented_and_unverified.md)・[`docs/47_v2_changes.md`](docs/47_v2_changes.md)・docs/48〜50・[`example/v2_3/`](example/v2_3/) の各まとめを参照してください。


## 元の節: 構成

ベンチ・回帰スクリプト・出力は入れない (`.gitignore`)。以下の一覧は 2026-11-18 の作業ツリー。


## 元の節: `pysad/` (本体)

| [`specfn_cx.py`](pysad/specfn_cx.py) | 複素の `Zeta`・`HurwitzZeta`・`PolyLog`・`HypergeometricU` (v2.1.0) |

| [`specfn_cx2.py`](pysad/specfn_cx2.py) | 複素の `ExpIntegralE`・`ExpIntegralEi`・`PolyGamma`・`WhittakerM/W`・`ProductLog`・`CosIntegral`・`FresnelS/C`・`GammaRegularized`・`Beta`・`Factorial`・`BesselJ` (複素の次数) など 15 関数 (v2.2.0) |

| [`nbcache.py`](pysad/nbcache.py) | numba のキャッシュの無効化 (`py -3.12 -m pysad.nbcache`。kernel を編集したら試験の前に走らせる。v2.1.0)。docs/12 §16.7 |

| [`synchrobeta.py`](pysad/synchrobeta.py) | `SynchroBetaEmittance` (v2.2.0 で始点の分散による過大を v2 で修正) |

| [`dapert.py`](pysad/dapert.py) | `DynamicApertureSurvey` (v2.2.0 で DAPWIDTH の窓を `trackd.f` の直訳に。既定 7) |

| [`spin.py`](pysad/spin.py) | POL (スピンと偏極。v2.1.0 でシンクロトロン側帯の項、v2.2.0 で厳密な緩和率・運動論的偏極・クロマティックな変調) |

| [`naff.py`](pysad/naff.py) | NAFF 風のチューン抽出 (`naff`・`naff_real`・`tune_xp`・`tunes_from_coords`。v2.1.0)。docs/12 §16.7 |

| [`trksetup.py`](pysad/trksetup.py) | 追跡の設定の共通部分 (v2.3.0。`TrackSetup`・`find_periodicity`・`select_locations`、閉軌道の求め直し、RF バケット) |

| [`freqmap.py`](pysad/freqmap.py) | 周波数マップ (v2.3.0。`frequency_map`、共鳴線の表。台本 `FrequencyMap[]`) |

| [`lma.py`](pysad/lma.py) | 局所運動量アクセプタンス (v2.3.0。`momentum_acceptance`、`touschek_lifetime_local`。台本 `MomentumAcceptance[]`・`TouschekLifetimeLMA[]`) |

| [`q1cmd.py`](pysad/q1cmd.py) | 周波数マップ・LMA の台本の関数の登録 (v2.3.0) |

| [`resdrive.py`](pysad/resdrive.py) | 6 極の共鳴駆動項 (R4 で偏向磁石の 3 次の項つき)・ADTS・クロマティシティ (v2.3.0。台本 `ResonanceDrivingTerms[]`・`AmplitudeTuneShift[]`・`ChromaticityTerms[]`) |

| [`adtsdirect.py`](pysad/adtsdirect.py) | 実際の強さの追跡 + NAFF + 多項式フィットで ADTS・チューンの δ 依存を直接求める (v2.3.0 R2。解析値 `resdrive.adts` の検証・基準値。`direct_adts`・`paraxial_tracker`) |

| [`bendnl.py`](pysad/bendnl.py) | 偏向磁石の 3 次の項 (曲率・硬い端・複合機能の本体・E1/E2) のリー生成関数を K64 の追跡の写像の 2 次のテイラー展開から取り出し、`resdrive.adts` の 2 次の摂動に渡す (v2.3.0 R3。K64 にない新機能) |

| [`rdttrack.py`](pysad/rdttrack.py) | 追跡から共鳴駆動項 (RDT) を取り出す (v2.3.0 R4。1 周の写像の 2 次のテイラー展開、閉軌道のヤコビ行列の δ 微分、3νx の近くの不動点。`resdrive.resonance_terms` の検証。K64 にない新機能) |

| [`fmapres.py`](pysad/fmapres.py) | 周波数マップの D の大きい領域を、RDT・ADTS から予想した共鳴の領域と比べる (v2.3.0 S1。台本ごとのまとめ [`example/v2_3/FMAP_summary.md`](example/v2_3/FMAP_summary.md)) |

| [`rdt2nd.py`](pysad/rdt2nd.py) | 4 次の共鳴駆動項 h_jklm (6 極の 2 次・8 極・運動学・端のフリンジ・偏向磁石の項。v2.3.0 R5。`resonance_terms` の 4 次の名前、台本 `Order2->`。K64 にない新機能) |

| [`nform.py`](pysad/nform.py) | 追跡の 1 周の写像の 3 次のテイラー展開からの数値の正規形 (J の 1 次の ADTS。v2.3.0 R5) |

| [`sxtopt.py`](pysad/sxtopt.py) | 6 極の K2L の最適化の**試作** (v2.3.0 R7。K2L の多項式モデル + 解析的な勾配 + SLSQP。結論: 検証には条件つきで使えるが、最終設定より DA の広い解は出せなかった。例 [`example/ddba/v07_14/sextupole/`](example/ddba/v07_14/sextupole/README.md)) |

| [`ampnf.py`](pysad/ampnf.py) | 振幅・δ に依る ν(Jx, Jy, δ) のべき級数 (v2.4.0。1 周の写像の TPSA + Birkhoff の正規形。J²・ξ の高次・振幅と δ の交差項。台本 `NormalFormTune[]`) |

| [`ampnf_ext.py`](pysad/ampnf_ext.py) | ampnf の拡張 (v2.4.1): 5・6 次の共鳴の駆動項の表 `drive_table`、1 本の共鳴を残した正規形 `resonance_normal_form` と島の理論 `island_theory`、追跡の島 `track_island` (検証用) |

| [`touschek_bunch.py`](pysad/touschek_bunch.py) | Touschek の損失率にバンチ内の δ の分布 (対の重心の δ_c と相対運動量の縦成分 q_s) を入れるオプション (v2.4.0。`touschek_piwinski_bunch`・`qs_ratio`) |

| [`spch_k64.py`](pysad/spch_k64.py) | PSPAC の K64 の形状関数 (k64 モード。v2.2.0。NPARA = 1 で gs と rms が 12 桁以上) |


## 元の節: `example/` (SAD の機能・関数ごとの短い例)

| [`example/v2_B/`](example/v2_B/README.md) | v2.0.0 の例 (特殊関数、文脈、共有メモリ・System・Fork、FFS の小項目)。既定の v2 モードで走る (pysad の出力のみ) |

| [`example/v2_C/`](example/v2_C/README.md) | v2.0.0 の例 (UND の磁場、PHSROT、SOL 領域内、ずれたソレノイド、スピン追跡、RADTAPER、LOSSMAP、RADLIGHT、MEASURE・TRCOD、空間電荷 ほか)。既定の v2 モードで走る (pysad の出力のみ) |

| [`example/v2_next/`](example/v2_next/) | v2.1.0・v2.2.0 の例 (P1〜P5: 偏極の緩和率と側帯、SynchroBeta・tradke・`face`、F1・AE の追跡、k64 の PSPAC・DAPWIDTH・WSPACSC・領域内のスピン・光子、MEASURE・DETR ≥ 1・複素の特殊関数・Tk 9.0。放射・F1・CAVI・陽子の追跡・PSPAC、スピン・MPW・RADINT、Emittance の F1・WSPAC・SOL、FREE・複素の特殊関数・Parallelize・文脈・Packages、偏極の側帯・光子・NAFF・幾何の fit、MEASURE・RadiationField・NPARA、NORFSW・光子・POL + 加速、MAIN の式・Micado・CSROYMatrix、ステアリングの F1・B108)。台本・gs の出力 (`.gs.txt`)・pysad の v2 (`.pysad.txt`)・k64 (`.pysad_k64.txt`)。A1・B1・B3・C3・D2 は README (`A1_README.md` ほか) |

| [`example/v2_3/`](example/v2_3/README.md) | v2.3.0 の例 (Q1: 周波数マップ・局所運動量アクセプタンス・gs との比較・台本の関数。Q2: 共鳴駆動項・ADTS・クロマティシティ・厚い ⇔ 薄い・追跡と gs の比較。R1〜R5: ADTS の追跡との検証・4 次の共鳴駆動項。S1: 周波数マップの検証。T1: elegant との比較。入門の図 `primer_*`。台本ごとのまとめ `ADTS_summary.md`・`FMAP_summary.md`。台本・出力・図) |

| [`example/11_graphics/fn_DRAW`](example/11_graphics/fn_DRAW.md) | v2.4.1 の `DRAW` (matplotlib): FODO と ASTRID2 の `DRAW BX BY & EX LAT;` の図 |

| [`example/11_graphics/fn_OpticsPlot_mpl`](example/11_graphics/fn_OpticsPlot_mpl.md) | v2.4.1 の `MatplotlibDrawer`: `OpticsPlot`・`Plot`・`ListPlot` の `Graphics[…]` を matplotlib に描いた図 (`$DisplayFunction=MatplotlibDrawer;`) |

| [`example/dxf_lattice/`](example/dxf_lattice/README.md) | v2.4.1 の `DrawMagnetsDXF[]` (ラティスの平面図の DXF): FODO・ASTRID2・KEK-LS の 1 セルと全周。台本・DXF・図・gs との一致 |

| [`example/v2_4/`](example/v2_4/README.md) | v2.4.0 の例 (A1: ampnf の ν(J, δ) と追跡の比較、A2: K0≠0 の偏向磁石の共鳴駆動項・ADTS・ξ、A3: Touschek のバンチ内の δ と対のモンテカルロ、A4: LMA の k64・結合・誤差、A5: 周波数マップと ν(J)、A7: 大きなリングの ampnf) |

| [`example/v2_4/I1_*`・`I2_*`](example/v2_4/I1_nu_of_J_orders.py) | v2.4.1 の ampnf の拡張: ν(J) の J の次数ごとの誤差 (FODO・ASTRID2・KEK-LS、gs の TrackParticles で基準を確認)、1 本の共鳴の島と追跡の島の比較 (非周期の FODO の 4〜6 次) |


## 元の節: `docs/` (解説文書。現状のスナップショットで、履歴は書かない)

| [`47_v2_changes.md`](docs/47_v2_changes.md) | pysad v2 の変更の一覧 (v2.0.0・v2.1.0・v2.2.0・v2.3.0。v1 = K64 に忠実な翻訳との違い。互換モード、v2 で結果が変わるもの、新しい機能、pysad の不具合の修正、要判断、試験の状態) |

| [`48_frequency_map_lma.md`](docs/48_frequency_map_lma.md) | 解説: 周波数マップと局所運動量アクセプタンス (LMA)。台本の関数 `FrequencyMap[]`・`MomentumAcceptance[]`・`TouschekLifetimeLMA[]` (v2.3.0) |

| [`49_resonance_terms_adts.md`](docs/49_resonance_terms_adts.md) | 解説: 6 極の共鳴駆動項・ADTS・クロマティシティ。台本の関数 `ResonanceDrivingTerms[]`・`AmplitudeTuneShift[]`・`ChromaticityTerms[]` (v2.3.0) |

| [`50_adts_rdt_primer.md`](docs/50_adts_rdt_primer.md) | **入門**: ADTS・RDT・周波数マップ・LMA を何をどうやって求めたか (v2.3.0。48・49 の前に読む。図と数値は `example/v2_3/primer_*`) |

| [`51_lattice_dxf.md`](docs/51_lattice_dxf.md) | 解説: ラティスの CAD 図 (DXF)。平面図の書き出し `DrawMagnetsDXF[]` / `pysad.dxf_lattice` (v2.4.1。使い方・描くもの・限界) |


## 元の節: オリジナルの SAD (gs) との一致度

| 全試験 (`pytest pysad/tests`、gs の出力との比較を含む、k64 が既定) | **2507 通過・13 skip・失敗 0 (2520 件収集。v2.4.1 の試験を含む)** (2026-11-23、v2.4.1 の最終の通し (1 時間 40 分)。その前の通しで落ちた 8 件は、四極の 4 次の合成 B226 と 6 極のフリンジ B291 による期待値の変化で、原因を切り分けて試験の期待を更新した。回帰 `QUICK=1` は 50 ベンチ FAIL 0。以前の記録: v2.4.0 は 2263 通過・9 skip・失敗 0 (2026-11-22)。`PYSAD_SAD_PACKAGES` つき、k64 が既定。skip の内訳: PYSAD_SLOW=1 が要る 4 件 ([`test_csr.py`](pysad/tests/test_csr.py)・[`test_v2_C3_measure.py`](pysad/tests/test_v2_C3_measure.py)・[`test_v2_D1_emit.py`](pysad/tests/test_v2_D1_emit.py)・[`test_v2_E_pol_kick.py`](pysad/tests/test_v2_E_pol_kick.py))、`PYSAD_ERR_LATTICE` が要る 2 件 ([`test_elkern.py`](pysad/tests/test_elkern.py))、POSIX だけ・fork が要る 2 件 ([`test_v2_B_script.py`](pysad/tests/test_v2_B_script.py)・[`test_v2_B2_script.py`](pysad/tests/test_v2_B2_script.py))) |

- **v2 で意図的に gs と違うもの**: v2 は K64 のバグ・物理の誤りを直すので結果が変わります ([`docs/47_v2_changes.md`](docs/47_v2_changes.md)。K64 側のバグは [`docs/21_sad_bug_report.md`](docs/21_sad_bug_report.md) の **B01〜B162 の 162 件**: v2.0.0 は B72〜B110、v2.1.0 は B111〜B134、v2.2.0 は B135〜B159、v2.3.0 は B160〜B162。B160 は誤りではなく K1 の暗黙の定義、B161 は CALC の行列の dk の不整合 (v2 で直した)、B162 は Touschek の損失率の式の適用範囲 (v2 は Piwinski 型が既定))。v2.3.0 の新機能 (上の節) は互換モードによらず同じ動作です。v2.3.0 で既存の v2 の結果を変えたのは `TouschekLifetime[INF, INF, nz]` の既定 (Piwinski 型。`Formula->"K64"` で K64 の表) と CALC の BEND の行列の dk (B161) だけです (`docs/47` の v2.3.0 の (a))。gs と同じ結果が要るときは `PYSAD_COMPAT=k64`。

詳細は [`docs/46_unimplemented_and_unverified.md`](docs/46_unimplemented_and_unverified.md)・[`docs/47_v2_changes.md`](docs/47_v2_changes.md)・[`example/README.md`](example/README.md)・[`pysad/tests/OUTPUTS.md`](pysad/tests/OUTPUTS.md)。

- 現状 (2026-11-21、v2.3.0。PyPI は v2.3.0 まで公開): 試験・回帰の件数は上の「オリジナルの SAD (gs) との一致度」の節。gs と合わない箇所・未対応は [`docs/history/work_unimpl.md`](docs/history/work_unimpl.md) と [`docs/history/00_RESUME.md`](docs/history/00_RESUME.md) §6 (例: SOL 領域内の四極 + RAD (gs の疑い)、最初の `Window[]` の親 (gs は Toplevel)、`Polarization` の追跡との比 (4.2 GeV・12 MV の 1 点)、偏向磁石の本体の B∥ (B136)、B108 の残る 2.7%、ユーザー定義の関数名 `At` が `@` と衝突して呼べない、`FitFunction` を `fit` と併用すると gs と別の解に落ちる)。v2.2.0 の要判断 (EMITDIV・DAPWIDTH・WSPACSC の既定など) は [`docs/47_v2_changes.md`](docs/47_v2_changes.md) の v2.2.0 の (d)、v2.3.0 の要判断・未解決 (共鳴線の予想の判別力、精密版の δ 微分、Touschek のバンチ内の δ 分布など) は同じく v2.3.0 の (d)。v2.3.0 の新機能 (周波数マップ・局所運動量アクセプタンス・Touschek の Piwinski 型の式・共鳴駆動項・ADTS・数値の正規形。K64 にない追加) の解説は [`docs/48_frequency_map_lma.md`](docs/48_frequency_map_lma.md)・[`docs/49_resonance_terms_adts.md`](docs/49_resonance_terms_adts.md) (入門は [`docs/50_adts_rdt_primer.md`](docs/50_adts_rdt_primer.md))。未対応・未確認の最新の一覧は [`docs/46_unimplemented_and_unverified.md`](docs/46_unimplemented_and_unverified.md)。

- **互換モード (v2.0.0 から)**: 既定は `v2` (K64 のバグを直した動作)。gs (K64) と同じ結果が要るときは環境変数 `PYSAD_COMPAT=k64` (または `pysad.compat.set_mode('k64')`)。既存の試験と回帰は gs との比較なので k64 モードで走る。v1 との違いの一覧は [`docs/47_v2_changes.md`](docs/47_v2_changes.md)、K64 側のバグは [`docs/21_sad_bug_report.md`](docs/21_sad_bug_report.md) (B01〜B162)、変更履歴は [`CHANGELOG.md`](CHANGELOG.md)。


## 元の節: [`pysad/tests/`](pysad/tests/) (試験)

1899 通過・8 skip・失敗 0 (1907 件収集。R7 の試験 [`test_ddba_sextupole.py`](pysad/tests/test_ddba_sextupole.py) 8 件を含む) (2026-11-21、v2.3.0 V3 の最終の通し。回帰 `QUICK=1` 50 ベンチ FAIL 0 (`QUICK=1`、133 s、k64)。skip の内訳: PYSAD_SLOW=1 が要る 4 件 ([`test_csr.py`](pysad/tests/test_csr.py)・[`test_v2_C3_measure.py`](pysad/tests/test_v2_C3_measure.py)・[`test_v2_D1_emit.py`](pysad/tests/test_v2_D1_emit.py)・[`test_v2_E_pol_kick.py`](pysad/tests/test_v2_E_pol_kick.py))、`PYSAD_ERR_LATTICE` が要る 2 件 ([`test_elkern.py`](pysad/tests/test_elkern.py))、POSIX だけ・fork が要る 2 件 ([`test_v2_B_script.py`](pysad/tests/test_v2_B_script.py)・[`test_v2_B2_script.py`](pysad/tests/test_v2_B2_script.py)))。以前: V1 (2026-11-20) は 1759 通過・8 skip、R5 の後は 1752 通過、R4 の後は 1718 通過、v2.2.0 は 1647 通過・9 skip、v2.1.0 は 1345 通過・8 skip、v2.0.0 は 937 通過・4 skip。参照は多くが SAD K64 (`gs`) の出力で、既存の試験は互換モード k64 で走る ([`conftest.py`](pysad/tests/conftest.py))。試験が走らせる SAD スクリプトと gs・pysad の画面表示は [`pysad/tests/scripts/`](pysad/tests/scripts/)、対応表は [`pysad/tests/OUTPUTS.md`](pysad/tests/OUTPUTS.md)。

| `test_v2_A_*.py`、`test_v2_B_*.py`、`test_v2_C_*.py`、`test_v2_C2_*.py` | v2.0.0 の変更の試験 (29 ファイル 363 件: 放射・エミッタンス・スピン、特殊関数・スクリプト・FFS、追跡・集団効果・要素、残りの仕上げ)。v2 の動作は試験の中で `compat.set_mode('v2')` を明示する。一覧は [`docs/47_v2_changes.md`](docs/47_v2_changes.md) |

| `test_v2_P*_*.py` (P1・P2・P3・P4・P5) | v2.2.0 の変更の試験 (19 ファイル: 偏極の厳密な緩和率・クロマティックな側帯・反転した BEND、SynchroBeta・tradke・RADCOD の微分・`face`、F1 と AE の追跡、k64 の PSPAC・DAPWIDTH の窓・FREE・領域内のスピンと光子・WSPACSC、k64 の MEASURE・DETR ≥ 1・複素の特殊関数・Micado・FITP・文法・Tk 9.0)。補助の独立な計算は `_p1_dk.py`・`_p1_models.py`・`_p2_indep.py`・`_p3_indep.py`・`_p3_gsdata.py`・`_p5_util.py`・`_p5_gscx.py` |

| [`test_v23_Q1_fmap_lma.py`](pysad/tests/test_v23_Q1_fmap_lma.py)、[`test_v23_Q2_resdrive.py`](pysad/tests/test_v23_Q2_resdrive.py)、[`test_v23_R1_rfsw_adts.py`](pysad/tests/test_v23_R1_rfsw_adts.py)、[`test_v23_R2_adts_direct.py`](pysad/tests/test_v23_R2_adts_direct.py)、[`test_v23_R3_bend_cubic.py`](pysad/tests/test_v23_R3_bend_cubic.py)、[`test_v23_R4_rdt.py`](pysad/tests/test_v23_R4_rdt.py)、[`test_v23_S1_fmapres.py`](pysad/tests/test_v23_S1_fmapres.py)、[`test_v23_R5_rdt4.py`](pysad/tests/test_v23_R5_rdt4.py)、[`test_v23_R6_bend_orbit.py`](pysad/tests/test_v23_R6_bend_orbit.py)、[`test_v23_S2_fmapcheck.py`](pysad/tests/test_v23_S2_fmapcheck.py)、[`test_v23_S3_fmappred.py`](pysad/tests/test_v23_S3_fmappred.py)、[`test_v23_L2_lma.py`](pysad/tests/test_v23_L2_lma.py)、[`test_v23_T1_elegant.py`](pysad/tests/test_v23_T1_elegant.py)、[`test_v23_T2_touschek.py`](pysad/tests/test_v23_T2_touschek.py)、[`test_v23_U1_bend_k1_field.py`](pysad/tests/test_v23_U1_bend_k1_field.py)、[`test_v23_U2_bend_k1_def.py`](pysad/tests/test_v23_U2_bend_k1_def.py)、PF 履歴・KEK-LS・DDBA の例の試験 ([`test_pf_history_lattices.py`](pysad/tests/test_pf_history_lattices.py)・[`test_kek_ls_example.py`](pysad/tests/test_kek_ls_example.py)・[`test_ddba_lattices.py`](pysad/tests/test_ddba_lattices.py)・[`test_ddba_v07_14.py`](pysad/tests/test_ddba_v07_14.py)・[`test_ddba_sextupole.py`](pysad/tests/test_ddba_sextupole.py) (R7、8 件)) | v2.3.0 の新機能と例の試験 (R6 の K0≠0・δ≠0・閉形式・非対角の正規形・`adts()` の既定と精密版、S2 の周波数マップの信頼性、S3 の共鳴線の予想、L2 の LMA の信頼性、T2 の Touschek の式、U2 の K1 の定義のオプションなどを含む。以下は R・S1 の内訳 (R5 の 4 次の共鳴駆動項・ADTS との整合・追跡の 3 次のテイラー展開・四極の積分器・数値の正規形・周波数マップの 14 件、S1 の周波数マップと RDT・ADTS の接続・誤差つき・δ・k64・MPW・gs の 12 件、R4 の共鳴駆動項の追跡との照合・偏向磁石の項 13 件、R3 の偏向磁石の 3 次の項 12 件を含む。以下は Q の計 35 件: 周波数マップ・共鳴線・LMA・周期性・Touschek への受け渡し・台本の関数・gs の出力との比較 (Q1 21 件)、共鳴駆動項・ADTS・クロマティシティ・厚いレンズ ⇔ 薄レンズ・追跡との照合 (Q2 14 件)) |

| [`test_v24_A1_ampnf.py`](pysad/tests/test_v24_A1_ampnf.py)、[`test_v24_A2_k0_bend.py`](pysad/tests/test_v24_A2_k0_bend.py)、[`test_v24_A3_touschek_bunch.py`](pysad/tests/test_v24_A3_touschek_bunch.py)、[`test_v24_A5_dapert_options.py`](pysad/tests/test_v24_A5_dapert_options.py)、[`test_v24_A7_misc.py`](pysad/tests/test_v24_A7_misc.py) | v2.4.0: `ampnf` (J¹ の収束・同じ模型の追跡との一致・ξ)、K0≠0 の偏向磁石の h・ADTS、Touschek のバンチ内の δ と q_s (対のモンテカルロとの一致)、`DynamicApertureSurvey` の `ExpandElementValues`・`ReferenceOrbit`、`errors_summary` の誤判定・`NormalFormTune[]` |

| [`test_v241_ampnf_ext.py`](pysad/tests/test_v241_ampnf_ext.py) | v2.4.1: `nu(jorder=)`・`d2nu`・J² の追跡との比較、6 次の島の理論と追跡の周期 6 の固定点・分離線、対称性で禁止された共鳴の駆動項、`k0_nonlinear` の独立性 (5 件) |

| [`test_v2_A1_next.py`](pysad/tests/test_v2_A1_next.py)、`test_v2_A2_*.py`、`test_v2_B1_*.py`、`test_v2_B2_*.py`、`test_v2_B3_*.py`、`test_v2_C3_*.py`、[`test_v2_D1_emit.py`](pysad/tests/test_v2_D1_emit.py)、`test_v2_D2_*.py`、`test_v2_E_*.py` | v2.1.0 の変更の試験 (25 ファイル 412 件: 放射・F1・CAVI・MASS/CHARGE、スピン角・Polarization・RADINT、Emittance の F1・WSPAC・SOL、FREE・複素の特殊関数・Parallelize、側帯・光子・NAFF、MEASURE・RadiationField・NPARA、NORFSW・光子・POL + 加速、MAIN の式・Micado・nbcache、傾いたソレノイドの POL・ステアリングの F1)。補助の独立な計算は `_a2_indep.py`・`_c3_indep.py`・`_e_spinkick.py`。遅い試験は `PYSAD_SLOW=1` |


## 元の節: [`docs/history/`](docs/history/) (作業履歴)

| [`v2_changes_A.md`](docs/history/v2_changes_A.md)、[`B`](docs/history/v2_changes_B.md)、[`C`](docs/history/v2_changes_C.md)、[`C_measure`](docs/history/v2_changes_C_measure.md)、[`C_spch`](docs/history/v2_changes_C_spch.md)、[`C2`](docs/history/v2_changes_C2.md)、[`C2_measure`](docs/history/v2_changes_C2_measure.md)、[`C2_misc`](docs/history/v2_changes_C2_misc.md)、[`C2_spch`](docs/history/v2_changes_C2_spch.md) | v2.0.0 の作業の記録 (根拠の数値、試験、使える範囲、要判断)。まとめは [`docs/47_v2_changes.md`](docs/47_v2_changes.md) |

| [`v23_R5.md`](docs/history/v23_R5.md)・[`v23_S1.md`](docs/history/v23_S1.md)・[`v23_Q1.md`](docs/history/v23_Q1.md)・[`Q2`](docs/history/v23_Q2.md)・[`T1`](docs/history/v23_T1.md) | v2.3.0 の作業の記録 (Q1: 周波数マップ・LMA、Q2: 共鳴駆動項・ADTS、T1: elegant との比較 (LMA・Touschek 寿命)、B160 (T1 が誤りの候補として見つけ、U1 が K1 の暗黙の定義に再分類。[`v23_U1.md`](docs/history/v23_U1.md)))。まとめは [`docs/47_v2_changes.md`](docs/47_v2_changes.md) |

| [`v24_A.md`](docs/history/v24_A.md)・[`v24_RAD.md`](docs/history/v24_RAD.md) | v2.4.0 の記録 (A: ADTS の J²・δ 依存 (ampnf)、K0≠0 の偏向磁石、Touschek のバンチ内の δ、LMA の対応範囲、DA の `ExpandElementValues`、周波数マップの予想・6 極検討の位置づけ、性能。RAD: 放射・エミッタンス (B163・RADINT の極・EMITDIV の既定)) |

| [`v2next_P1.md`](docs/history/v2next_P1.md)〜[`P5`](docs/history/v2next_P5.md) | v2.2.0 の作業の記録 (P1: 偏極、P2: 放射・エミッタンス、P3: F1・AE、P4: SOL・空間電荷・DAPWIDTH・FREE、P5: gs との検証)。まとめは [`docs/47_v2_changes.md`](docs/47_v2_changes.md) |

| [`v2next_A1.md`](docs/history/v2next_A1.md)、[`A2`](docs/history/v2next_A2.md)、[`B1`](docs/history/v2next_B1.md)、[`B2`](docs/history/v2next_B2.md)、[`B3`](docs/history/v2next_B3.md)、[`C`](docs/history/v2next_C.md)、[`D1`](docs/history/v2next_D1.md)、[`D2`](docs/history/v2next_D2.md)、[`E`](docs/history/v2next_E.md) | v2.1.0 の作業の記録 (根拠の数値、試験、使える範囲、要判断)。まとめは [`docs/47_v2_changes.md`](docs/47_v2_changes.md) |
