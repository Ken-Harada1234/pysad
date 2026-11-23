# script_tests/ (スクリプトテスト)

pysad の機能ごとの自動試験 (pytest)。198 本の試験ファイル `test_*.py`、収集 2534 件 (`PYSAD_SAD_PACKAGES` を設定し、旧式の 3 本を `--ignore` した `--collect-only` の値)。
各試験は、小さなラティスや SAD スクリプトを pysad で走らせ、(a) SAD (K64) で取って保存した出力、(b) 試験の中の解析式、(c) SAD・pysad のコードを使わない独立な数値計算 (`_*_indep.py` など) のどれかと比べる。SAD (K64) の実行は要らない。

用語: 「SAD 照合」= 保存した SAD (K64) の出力との比較。「独立」= 独立な数値計算や解析式との比較。試験は `conftest.py` により既定で互換モード k64 (`PYSAD_COMPAT=k64`) で走る。v2 の動作の試験は、試験の中で `v2_mode` フィクスチャか `compat.set_mode('v2')` を明示する。

SAD (K64) との照合ベンチは姉妹フォルダ [`../regression_tests/`](../regression_tests/README.md)。上位の説明は [`../README.md`](../README.md)。

## 目次

1. [実行方法](#実行方法)
2. [区分 1: 光学・CALC・FFS コマンド](#区分-1-光学calcffs-コマンド)
3. [区分 2: 要素ごとの試験 (BEND・MULT・SEXT・COORD・MAP・UND・TCAV・APERT・CAVI)](#区分-2-要素ごとの試験-bendmultsextcoordmapundtcavapertcavi)
4. [区分 3: エミッタンス・放射](#区分-3-エミッタンス放射)
5. [区分 4: 追跡・乱数・誤差](#区分-4-追跡乱数誤差)
6. [区分 5: DA・周波数マップ・LMA・Touschek・非線形力学](#区分-5-da周波数マップlmatouschek非線形力学)
7. [区分 6: スピン・偏光](#区分-6-スピン偏光)
8. [区分 7: SOL・ソレノイド](#区分-7-solソレノイド)
9. [区分 8: 集団効果 (wake・CSR・ビームビーム・空間電荷)](#区分-8-集団効果-wakecsrビームビーム空間電荷)
10. [区分 9: SADScript の評価器・パーサ・パッケージ](#区分-9-sadscript-の評価器パーサパッケージ)
11. [区分 10: 特殊関数・数学関数](#区分-10-特殊関数数学関数)
12. [区分 11: GUI・グラフィック (Tk・matplotlib・DXF)](#区分-11-guiグラフィック-tkmatplotlibdxf)
13. [区分 12: マッチング・FIT・軌道補正](#区分-12-マッチングfit軌道補正)
14. [区分 13: ラティス集・例の確認](#区分-13-ラティス集例の確認)
15. [区分 14: ツール・測定](#区分-14-ツール測定)
16. [部品・データ・scripts/](#部品データscripts)
17. [試験の足し方](#試験の足し方)

## 実行方法

リポジトリのルートから実行する (Python 3.12)。

```
set PYSAD_SAD_PACKAGES=<SAD の Packages フォルダ>        (PowerShell: $env:PYSAD_SAD_PACKAGES="...")
py -3.12 -m pysad.nbcache --force                          (numba のキャッシュを作り直す。初回・カーネルを直したあと)
py -3.12 -m pytest pysad/tests/script_tests -q -p no:cacheprovider `
    --ignore=pysad/tests/script_tests/test_bend_feeddown.py --ignore=pysad/tests/script_tests/test_coord.py --ignore=pysad/tests/script_tests/test_sext_radiation.py
py -3.12 -m pytest pysad/tests/script_tests/test_wake.py -q  (1 ファイルだけ)
```

- `--ignore` する 3 本は旧式のスクリプト形式で、import 時に `sys.exit` を呼ぶため pytest の収集を止める。
- `PYSAD_SAD_PACKAGES` が無いと GUI 系 (`test_plot.py`・`test_tclbridge.py` ほか) が skip する。
- 環境変数: `PYSAD_SLOW=1` (遅い試験も走らせる)、`PYSAD_ERR_LATTICE` (`test_elkern.py` の誤差入りラティス 2 件)、`PYSAD_PY314` (`test_v2_P5_tk90.py` の Tk 9.0 用 Python)、`V241_POL_CACHE` (偏光の比較の部品)。
- 全部で 1 時間 45 分ほど (FITP の five が約 10 分以上)。

## 区分 1: 光学・CALC・FFS コマンド

合計 17 ファイル 151 件。

| ファイル | 件数 | 何を確かめるか |
|---|---|---|
| [`test_pysad.py`](test_pysad.py) | 49 | ドリフト・四極・偏向の行列、Twiss、チューンの基本の数値検証 |
| [`test_optics_fn.py`](test_optics_fn.py) | 8 | TransferMatrix・BeamMatrix・NormalCoordinates・OrbitGeo などの光学関数 |
| [`test_ffs_gs_match.py`](test_ffs_gs_match.py) | 10 | FFS 内の MOMENTUM、Twiss["LENGTH"]、EXPAND/NOEXPAND の動作 (SAD 照合) |
| [`test_ffsextra.py`](test_ffsextra.py) | 8 | FFS の VARIABLES・MATRIX・SHOW・RECOVER・REVERSE・出力切替 |
| [`test_v2_A_ffs.py`](test_v2_A_ffs.py) | 9 | 4D の CALC が縦方向の列を戻すか、UnstableLevel、CalculateOptics の縦方向 |
| [`test_v2_B_ffs.py`](test_v2_B_ffs.py) | 14 | MON/MCOD/TCOD、LINE["LENGTH"]、ORG の負数、SeedRandom など FFS の小項目 |
| [`test_v2_C2_optics.py`](test_v2_C2_optics.py) | 7 | 加速空洞内の分散、CalculateOptics の TRPT 加速、RADTAPER の閉軌道 |
| [`test_trpt_acc.py`](test_trpt_acc.py) | 7 | TRPT の加速 (gammab)・MARK OFFSET・USE の大文字小文字・K=0 の四極 |
| [`test_unstable.py`](test_unstable.py) | 7 | 不安定なラインで CELL CALC が続行し nstab を判定する |
| [`test_geo_fit_sig.py`](test_geo_fit_sig.py) | 4 | 幾何の fit (GX/GY/GZ/CHI)・DISP G・ビームサイズ Twiss["SIG.."] |
| [`test_tdgeo.py`](test_tdgeo.py) | 2 | 幾何量の解析微分 (tdgeo) が有限差分と一致する |
| [`test_v2_B3_geofit.py`](test_v2_B3_geofit.py) | 4 | 四極の DX/DY の幾何微分と、SOL 領域内の四極のずれ |
| [`test_v2_P5_detr.py`](test_v2_P5_detr.py) | 11 | DETR≥1 (x-y モード入れ替え) の CALC・INS・Emittance[] (SAD 照合) |
| [`test_v2_A2_dx.py`](test_v2_A2_dx.py) | 10 | ずれた BEND の CALC と追跡の整合 (出口の幾何) |
| [`test_v24_H3_qcod_robust.py`](test_v24_H3_qcod_robust.py) | 1 | 大きなリングで周期閉軌道の探索 (qcod) が収束する |
| [`test_error_spreading.py`](test_error_spreading.py) | 0 | FFSState の誤差ばらまき関数の単体 (スクリプト形式) |
| [`test_new_features.py`](test_new_features.py) | 0 | MARK 先頭の初期軌道オフセット、save()/restore_saved() など (スクリプト形式) |

## 区分 2: 要素ごとの試験 (BEND・MULT・SEXT・COORD・MAP・UND・TCAV・APERT・CAVI)

合計 27 ファイル 442 件。

| ファイル | 件数 | 何を確かめるか |
|---|---|---|
| [`test_apert.py`](test_apert.py) | 13 | APERT 要素 (楕円・矩形・ずれ・回転・再入射) の追跡 (SAD 照合) |
| [`test_cavrot.py`](test_cavrot.py) | 2 | CAVI と VOLT 付き MULT のずれ・回転 (SAD 照合) |
| [`test_coord.py`](test_coord.py) | 0 | COORD 要素の各分岐の単体 (スクリプト形式。`sys.exit` で収集を止める) |
| [`test_coord_elem.py`](test_coord_elem.py) | 12 | COORD 要素の追跡・Twiss・エミッタンス (SAD 照合) |
| [`test_map.py`](test_map.py) | 6 | MAP 要素 (ExternalMap) の OPTICS・EMIT・TRACK (SAD 照合) |
| [`test_und.py`](test_und.py) | 4 | UND (アンジュレータ) の追跡・Twiss・エミッタンス (SAD 照合) |
| [`test_v2_C_und.py`](test_v2_C_und.py) | 9 | UND の磁場 FBX・FBY を解析式 (平面アンジュレータ) と比べる |
| [`test_tcav.py`](test_tcav.py) | 8 | TCAVI (横偏向空洞) の追跡とエミッタンス (SAD 照合) |
| [`test_mult_angle.py`](test_mult_angle.py) | 21 | ANGLE≠0 の MULT の追跡・エミッタンス・CALC (SAD 照合) |
| [`test_v24_A2_mult_angle_twiss.py`](test_v24_A2_mult_angle_twiss.py) | 1 | ANGLE≠0 の MULT の CALC (Twiss) が 10 桁で一致 |
| [`test_bendk1_track.py`](test_bendk1_track.py) | 2 | ANGLE=0 で K1≠0 の BEND の追跡と粒子損失 |
| [`test_bend_feeddown.py`](test_bend_feeddown.py) | 0 | Bend.constant_kick() の DX/DY フィードダウン (スクリプト形式) |
| [`test_sext_radiation.py`](test_sext_radiation.py) | 0 | 6 極の厚い放射キック: E0_GeV=0 で放射なし (スクリプト形式) |
| [`test_rotate.py`](test_rotate.py) | 16 | ROTATE (ROT/TILT) 付きの BEND・四極の行列・CALC |
| [`test_v241_A_bend_k1_k0.py`](test_v241_A_bend_k1_k0.py) | 36 | ANGLE=0 で K1≠0・K0≠0 の BEND が MULT と同じ結果 |
| [`test_v241_B_sext_fringe.py`](test_v241_B_sext_fringe.py) | 14 | 6 極・8 極などの端の非線形フリンジを CALC・Emittance[] へ入れる |
| [`test_v241_C_mult_edge_ykick.py`](test_v241_C_mult_edge_ykick.py) | 12 | MULT の端 (AE1・AE2) の縦方向の収束 |
| [`test_v241_quad_integrator.py`](test_v241_quad_integrator.py) | 31 | 四極本体の積分器の構造と、追跡の写像の正確さ |
| [`test_v244_quad_linear_fringe_f1.py`](test_v244_quad_linear_fringe_f1.py) | 68 | 四極・MULT の線形フリンジの係数 F1 の符号 (追跡と CALC・数値積分・gs の保存した出力、B323) |
| [`test_v250_chromaticity_bend.py`](test_v250_chromaticity_bend.py) | 19 | `Lattice.chromaticity()` の dν/dδ が偏向磁石 (なし・E1/E2・K1≠0・6 極つき) で追跡のヤコビアン・gs の `Twiss["NX",-1]` の差分と一致。`CHROMATICITY` コマンドの単位と BEND の寄与の欠落 (B325)、厚い 6 極の CALC と追跡のずれの切り分け・gs の TrackParticles (B326) |
| [`test_v2_P3_ae.py`](test_v2_P3_ae.py) | 24 | MULT の端の角 AE1・AE2 と K0 つき MULT のエミッタンス (SAD 照合) |
| [`test_v2_P3_f1.py`](test_v2_P3_f1.py) | 100 | ANGLE=0 の BEND と MULT の F1・FB1・FB2 の追跡 (SAD 照合) |
| [`test_v2_D2_edge.py`](test_v2_D2_edge.py) | 36 | 硬い端の極 (ANGLE=0 の BEND) の端の幾何を独立な基準と比べる |
| [`test_v2_B1_optics.py`](test_v2_B1_optics.py) | 42 | BEND の F1・FRINGE の Emittance[]・CALC・追跡 (SAD 照合と独立な厳密写像) |
| [`test_v2_E_steer.py`](test_v2_E_steer.py) | 19 | ANGLE=0 の BEND と MULT の F1 の Emittance[]・CALC (SAD 照合) |
| [`test_v2_A_tchge.py`](test_v2_A_tchge.py) | 8 | 据え付け誤差の座標変換 tchge の転送行列が写像のヤコビアン |
| [`test_v2_C_phsrot.py`](test_v2_C_phsrot.py) | 5 | PHSROT (位相空間の回転) がシンプレクティックでチューンを保つ |
| [`test_v23_U1_bend_k1_field.py`](test_v23_U1_bend_k1_field.py) | 7 | K1≠0 の BEND の K1 の定義を独立な数値積分と比べる |
| [`test_v23_U2_bend_k1_def.py`](test_v23_U2_bend_k1_def.py) | 14 | K1 の定義オプション ('log'/'lin') と CALC の dk |

## 区分 3: エミッタンス・放射

合計 24 ファイル 284 件。

| ファイル | 件数 | 何を確かめるか |
|---|---|---|
| [`test_emit_jac.py`](test_emit_jac.py) | 22 | emit6d の解析ヤコビアンが数値微分と一致し、tcod の直訳と合う |
| [`test_emit_report.py`](test_emit_report.py) | 37 | EMIT 報告の全行と、四極・薄い多極の解析行列 (SAD 保存値と比較) |
| [`test_emit_unstable.py`](test_emit_unstable.py) | 5 | 不安定な光学の EMIT / INS+CALC6D / tcod の回帰 |
| [`test_emit_extra.py`](test_emit_extra.py) | 0 | 空ファイル |
| [`test_synchrobeta.py`](test_synchrobeta.py) | 2 | SynchroBetaEmittance (SAD 照合) |
| [`test_small_items.py`](test_small_items.py) | 7 | 小さなリングの RADINT・DISP の位置範囲・RADINT の単精度の定数 |
| [`test_v2_A_radint.py`](test_v2_A_radint.py) | 13 | RADINT (放射積分) の曲げごとの半径・色収差の積分を独立な積分と比べる |
| [`test_v2_A2_radint.py`](test_v2_A2_radint.py) | 12 | RADINT・Emittance[]・追跡の整合を、大きなリングで自前の積分と比べる |
| [`test_v2_C3_radint.py`](test_v2_C3_radint.py) | 21 | 傾いた BEND・分散のある場所の I4 の独立な確認 |
| [`test_v2_A_damping.py`](test_v2_A_damping.py) | 8 | 放射減衰の分配 (Jx, Jy, Je) の端の角と勾配つき BEND |
| [`test_v24_rad_pole.py`](test_v24_rad_pole.py) | 4 | ANGLE=0・K0≠0 の極 (挿入光源) の RADINT と Emittance[] |
| [`test_v24_D_radbeta.py`](test_v24_D_radbeta.py) | 12 | 放射の β の因子と、BEND の T[z, δ] の厳密な行列 |
| [`test_v24_D3_rad_rf_sol.py`](test_v24_D3_rad_rf_sol.py) | 12 | RF 電圧・光子の β の因子と、SOL 内の四極の放射 |
| [`test_v2_C3_misc.py`](test_v2_C3_misc.py) | 16 | 結合リングの平衡エミッタンス、FITP>1 の最終解、FLUC の乱数列 |
| [`test_v2_P2_rad.py`](test_v2_P2_rad.py) | 22 | SynchroBetaEmittance の始点の分散 η≠0 など、放射・エミッタンスの残り |
| [`test_v2_D1_emit.py`](test_v2_D1_emit.py) | 22 | NORFSW+閉軌道のずれの Emittance[]、光子、WSPAC の最後のスライス |
| [`test_v2_P5_proton_beta.py`](test_v2_P5_proton_beta.py) | 6 | 陽子の放射の β の因子 (U0 が古典値/β³) |
| [`test_k64_tcod_notfound.py`](test_k64_tcod_notfound.py) | 3 | 閉軌道が収束しないときの Emittance[] の動作 |
| [`test_radfield.py`](test_radfield.py) | 2 | RadiationField・RadiationSpectrum (SAD 照合) |
| [`test_v2_C3_radfield.py`](test_v2_C3_radfield.py) | 8 | RadiationField を Liénard-Wiechert の式の数値計算と比べる |
| [`test_radlight.py`](test_radlight.py) | 1 | RADLIGHT (軌道の記録) (SAD 照合) |
| [`test_v2_C_radlight.py`](test_v2_C_radlight.py) | 3 | RADLIGHT の軌道が幾何と一致し、BEND 内を記録する |
| [`test_v2_B3_photons.py`](test_v2_B3_photons.py) | 6 | PHOTONS (光子の記録 PhotonList) と Stokes パラメータの範囲 |
| [`test_v2_A1_next.py`](test_v2_A1_next.py) | 40 | INS の Twiss["DETR"]、陽子の質量・電荷での追跡、ほか (解析式) |

## 区分 4: 追跡・乱数・誤差

合計 14 ファイル 160 件。

| ファイル | 件数 | 何を確かめるか |
|---|---|---|
| [`test_elkern.py`](test_elkern.py) | 37 | 要素ごとの numba カーネルが Python の経路と同じ結果を出す |
| [`test_perf_kernels.py`](test_perf_kernels.py) | 50 | 高速化カーネル (K1≠0 BEND・tchge・apply) が基準の Python 版と一致 |
| [`test_track6n_ext.py`](test_track6n_ext.py) | 8 | numba 追跡 (track6n) の MULT・SOL が純 Python と一致 (今回の収集では除外) |
| [`test_main_track.py`](test_main_track.py) | 8 | MAIN レベルの TRACK コマンド (規格化振幅・smear・チューン) |
| [`test_fourier_fit.py`](test_fourier_fit.py) | 4 | MAIN の TRACK の FOURIER 解析 (wjfit・wfres) (K64 の Fortran と照合) |
| [`test_v2_A_track.py`](test_v2_A_track.py) | 8 | 追跡の dv の項と、SOL 内の四極の放射が他粒子に依らない |
| [`test_v2_C_track.py`](test_v2_C_track.py) | 14 | tdusrn の NaN、最初の tapert、SOL 領域の未対応要素、LOSSMAP |
| [`test_v2_F2_photons_coord_tcavi.py`](test_v2_F2_photons_coord_tcavi.py) | 2 | PHOTONS と COORD・TCAVI を同じラティスに置いたときの光子の記録 |
| [`test_sadrandom.py`](test_sadrandom.py) | 3 | SAD の既定乱数と Random[]/GaussRandom[] の束縛 |
| [`test_delcor.py`](test_delcor.py) | 3 | DELCOR (相関のあるランダムなずれ) (SAD 照合) |
| [`test_v2_C_delcor.py`](test_v2_C_delcor.py) | 3 | DELCOR で MARK・CAVI へ書く先と、相関関数の形 |
| [`test_v2_C2_delcor.py`](test_v2_C2_delcor.py) | 7 | DELCOR の相関長≠0 の規格化 (各成分の rms が 1) |
| [`test_v2_B3_naff.py`](test_v2_B3_naff.py) | 12 | NAFF 風のチューン抽出が単一・複数周波数を正しく出す |
| [`test_v2_F3_measure_fifo.py`](test_v2_F3_measure_fifo.py) | 1 | PUTMEASURE/GETMEASURE を名前つきパイプ (FIFO) で受け渡す |

## 区分 5: DA・周波数マップ・LMA・Touschek・非線形力学

合計 29 ファイル 424 件。

| ファイル | 件数 | 何を確かめるか |
|---|---|---|
| [`test_v23_Q1_fmap_lma.py`](test_v23_Q1_fmap_lma.py) | 21 | 周波数マップ・局所運動量アクセプタンス (LMA)・共通の追跡設定 |
| [`test_v23_L2_lma.py`](test_v23_L2_lma.py) | 11 | LMA の check_every の既定、RF バケット、限界の種別 |
| [`test_v23_Q2_resdrive.py`](test_v23_Q2_resdrive.py) | 14 | 6 極の共鳴駆動項・ADTS・クロマティシティを追跡と比べる |
| [`test_v23_R1_rfsw_adts.py`](test_v23_R1_rfsw_adts.py) | 2 | RF つきラティスの 6 極 ×0 の ADTS (RFSW の ON/OFF) |
| [`test_v23_R2_adts_direct.py`](test_v23_R2_adts_direct.py) | 8 | ADTS の直接の検証 (追跡+NAFF+多項式フィット) |
| [`test_v23_R3_bend_cubic.py`](test_v23_R3_bend_cubic.py) | 12 | 偏向磁石の 3 次の項を入れた ADTS |
| [`test_v23_R4_rdt.py`](test_v23_R4_rdt.py) | 13 | 共鳴駆動項 h_jklmp の解析式を、追跡の 1 周写像から取り出した値と比べる |
| [`test_v23_R5_rdt4.py`](test_v23_R5_rdt4.py) | 14 | 4 次の共鳴駆動項 (6 極の 2 次を含む) |
| [`test_v23_R6_bend_orbit.py`](test_v23_R6_bend_orbit.py) | 15 | 偏向磁石の閉軌道まわりの分解、K0≠0 の磁石、δ≠0 の ADTS |
| [`test_v23_S1_fmapres.py`](test_v23_S1_fmapres.py) | 12 | 周波数マップと RDT・ADTS の接続、誤差つき・δ つきの周波数マップ |
| [`test_v23_S2_fmapcheck.py`](test_v23_S2_fmapcheck.py) | 19 | 周波数マップの警告・周回数の換算・誤差の与え方 |
| [`test_v23_S3_fmappred.py`](test_v23_S3_fmappred.py) | 12 | 共鳴線の予想 (保存した周波数マップの表を読むので速い) |
| [`test_v23_T1_elegant.py`](test_v23_T1_elegant.py) | 8 | elegant との比較 (変換、保存した elegant の LMA・寿命) |
| [`test_v23_T2_touschek.py`](test_v23_T2_touschek.py) | 8 | Touschek 寿命の式を独立な関数と比べる |
| [`test_v24_A1_ampnf.py`](test_v24_A1_ampnf.py) | 7 | 1 周写像の Birkhoff 正規形による ν(Jx, Jy, δ) |
| [`test_v24_A2_k0_bend.py`](test_v24_A2_k0_bend.py) | 8 | K0≠0 の BEND (挿入光源型) の解析式を追跡の写像と比べる |
| [`test_v24_A3_das_damping.py`](test_v24_A3_das_damping.py) | 3 | DynamicApertureSurvey の Damping・DampEnough (SAD 照合) |
| [`test_v24_A3_touschek_bunch.py`](test_v24_A3_touschek_bunch.py) | 4 | Touschek の損失率にバンチ内の δ 分布を入れるオプション |
| [`test_v24_A5_dapert_options.py`](test_v24_A5_dapert_options.py) | 2 | DynamicApertureSurvey の ExpandElementValues・ReferenceOrbit |
| [`test_v24_A7_misc.py`](test_v24_A7_misc.py) | 4 | Touschek の q_s の効果 (モンテカルロと一致)、誤差の集計 |
| [`test_v24_F1_rotated_bend.py`](test_v24_F1_rotated_bend.py) | 7 | ROTATE した磁石を RDT・ADTS・ν(J, δ) に入れる |
| [`test_v24_F2_defaults.py`](test_v24_F2_defaults.py) | 5 | v2 の既定値 (DAS・周波数マップ・Touschek) と k64 の違い |
| [`test_v24_G_k0_ampnf.py`](test_v24_G_k0_ampnf.py) | 154 | K0≠0 の磁石 (MPW・逆曲げ・傾き) を正規形 (TPSA) に入れる |
| [`test_v24_H1_coupled_h4.py`](test_v24_H1_coupled_h4.py) | 11 | 結合した光学 (ROTATE した磁石) の 4 次の h・δ≠0 |
| [`test_v241_ampnf_ext.py`](test_v241_ampnf_ext.py) | 5 | 正規形の J² の項・島の理論と追跡の比較 |
| [`test_v241_lma_long.py`](test_v241_lma_long.py) | 16 | 誤差つき LMA と 10^4 周級の放射追跡 (保存した gs の出力と比較) |
| [`test_v241_sxtopt_spec.py`](test_v241_sxtopt_spec.py) | 10 | 6 極の最適化 (sxtopt) の項の指定と作動点まわりの共鳴の警告 |
| [`test_v2_P4_dapwidth.py`](test_v2_P4_dapwidth.py) | 11 | DynamicApertureSurvey の DAPWIDTH の窓 (SAD 照合) |
| [`test_ddba_sextupole.py`](test_ddba_sextupole.py) | 8 | 6 極ファミリーの K2L に対する RDT・ADTS の多項式モデルと最適化 |

## 区分 6: スピン・偏光

合計 18 ファイル 186 件。

| ファイル | 件数 | 何を確かめるか |
|---|---|---|
| [`test_spin.py`](test_spin.py) | 8 | POL: PolarizationTime・SpinTune・SpinPrecessionMatrix (SAD 照合) |
| [`test_spin_track.py`](test_spin_track.py) | 5 | POL つきの TrackParticles (粒子ごとのスピン追跡) (SAD 照合) |
| [`test_v241_pol_fluc.py`](test_v241_pol_fluc.py) | 22 | FLUC つきスピン追跡を 4 つのリングで gs と比べる |
| [`test_v24_pol_emitdiv.py`](test_v24_pol_emitdiv.py) | 2 | POL の Emittance[] の EMITDIV への依存 |
| [`test_v24_D_spin_bpar.py`](test_v24_D_spin_bpar.py) | 4 | 磁石本体の運動量方向の磁場 B∥ の効果 (Maxwell を満たす基準) |
| [`test_v24_D3_spin_bpar_magnets.py`](test_v24_D3_spin_bpar_magnets.py) | 17 | 四極・6 極・MULT・ステアリング本体の B∥ の効果 |
| [`test_v2_A2_pol.py`](test_v2_A2_pol.py) | 7 | Polarization (Sokolov-Ternov) を独立な式と比べる |
| [`test_v2_A2_spin.py`](test_v2_A2_spin.py) | 9 | 500 周後のスピン角の違い (スライス角の誤差) の切り分け |
| [`test_v2_A_spin.py`](test_v2_A_spin.py) | 22 | スピンと SOL 内の放射: 運動量の共分散、端の縦磁場 |
| [`test_v2_B3_pol.py`](test_v2_B3_pol.py) | 7 | Polarization のシンクロトロン側帯 (スピン位相変調) の項 |
| [`test_v2_B3_radpol.py`](test_v2_B3_radpol.py) | 3 | RADPOL 追跡の緩和が Sokolov-Ternov の値と合う |
| [`test_v2_C2_spin.py`](test_v2_C2_spin.py) | 17 | SOL 領域・UND・PHSROT のスピン追跡 (Thomas-BMT 積分と比較) |
| [`test_v2_E_pol_kick.py`](test_v2_E_pol_kick.py) | 3 | δ の蹴りへのスピン軸の応答 G の検証 |
| [`test_v2_E_spinsol.py`](test_v2_E_spinsol.py) | 21 | SOL 領域 (傾いたソレノイド) の Emittance[] のスピン量 |
| [`test_v2_P1_pol_chrom.py`](test_v2_P1_pol_chrom.py) | 6 | 色収差による変調で、偏極の式が FLUC 追跡に近づく |
| [`test_v2_P1_rev.py`](test_v2_P1_rev.py) | 6 | 反転した BEND のあるリングの偏極度・時間が独立な式と合う |
| [`test_v2_P1_st.py`](test_v2_P1_st.py) | 9 | n が傾く・非平面リングの偏極の緩和率 |
| [`test_v2_P4_sol_spin.py`](test_v2_P4_sol_spin.py) | 18 | SOL 領域内の四極・MULT・CAVI のスピン (Thomas-BMT と比較) |

## 区分 7: SOL・ソレノイド

合計 9 ファイル 112 件。

| ファイル | 件数 | 何を確かめるか |
|---|---|---|
| [`test_sol_optics.py`](test_sol_optics.py) | 4 | SOL 領域の Twiss (CALC、INS、CELL) (SAD 照合) |
| [`test_sol_emit.py`](test_sol_emit.py) | 3 | SOL 領域の Emittance[] (SAD 照合) |
| [`test_sol_geo.py`](test_sol_geo.py) | 2 | SOL 領域の幾何 (ずれ・傾きがあるときの LINE["GX"] など) (SAD 照合) |
| [`test_v2_B1_sol.py`](test_v2_B1_sol.py) | 45 | SOL 領域内の要素を、ローレンツ力の積分と比べる (大きな Bz まで) |
| [`test_v2_C_sol.py`](test_v2_C_sol.py) | 32 | SOL 領域内の BEND・MULT・CAVI・APERT の追跡と放射 |
| [`test_v2_C2_sol.py`](test_v2_C2_sol.py) | 11 | 傾き・ずれた SOL (DX DY DZ CHI) と端の放射 F1、領域の幾何 |
| [`test_v2_C2_solrf.py`](test_v2_C2_solrf.py) | 9 | SOL 領域の CAVI・RF つき MULT の加速、wake、AUTOPHI |
| [`test_v2_F4_sol_bound0.py`](test_v2_F4_sol_bound0.py) | 2 | BOUND=0 の SOL が領域の外にあるときの追跡 (SAD 照合) |
| [`test_v2_P4_sol_photons.py`](test_v2_P4_sol_photons.py) | 4 | SOL 領域内の光子の記録 (PhotonList) |

## 区分 8: 集団効果 (wake・CSR・ビームビーム・空間電荷)

合計 15 ファイル 159 件。

| ファイル | 件数 | 何を確かめるか |
|---|---|---|
| [`test_wake.py`](test_wake.py) | 23 | WakeFunction の追跡 (SAD 照合) |
| [`test_wake_optics.py`](test_wake_optics.py) | 9 | WAKEOPT (光学の wake) と InitialOrbits (SAD 照合) |
| [`test_v2_C_wakecmd.py`](test_v2_C_wakecmd.py) | 3 | 旧式の WAKE コマンドが WakeFunction と同じ結果 |
| [`test_v2_C2_wakecmd.py`](test_v2_C2_wakecmd.py) | 6 | WAKE コマンドの表を複数行にまたがって書ける |
| [`test_csr.py`](test_csr.py) | 7 | CSR (CSRTrack・CSRHaissin・CSRSetupOY・CSRMatrix) (SAD 照合) |
| [`test_v2_C_csr.py`](test_v2_C_csr.py) | 2 | CSRSetupOY の位置 zj と運動量 dpj が別の値 |
| [`test_beambeam.py`](test_beambeam.py) | 24 | BEAMBEAM (強弱ビーム) の追跡・ルミノシティ (SAD 照合) |
| [`test_v24_bb_tilt_indep.py`](test_v24_bb_tilt_indep.py) | 5 | 傾いた強ビームの場を極座標の数値積分と比べる |
| [`test_v2_C_beambeam.py`](test_v2_C_beambeam.py) | 2 | BEAMBEAM の傾き回転が共変 |
| [`test_bbbrem.py`](test_bbbrem.py) | 1 | BBBrem1 (ビーム・ビーム制動放射の 1 事象) (SAD 照合) |
| [`test_gausscoul.py`](test_gausscoul.py) | 2 | GaussianCoulomb・GaussianCoulombU (SAD 照合) |
| [`test_v2_C_spch.py`](test_v2_C_spch.py) | 22 | 空間電荷 (SPAC・PSPAC・SPCH 要素) を解析式と比べる |
| [`test_v2_C2_spch.py`](test_v2_C2_spch.py) | 35 | MULT の RADIUS (SPAC) と WSPAC (ガウス分布) の空間電荷 |
| [`test_v2_P4_pspac_k64.py`](test_v2_P4_pspac_k64.py) | 14 | k64 モードの PSPAC (gs と全桁一致) |
| [`test_v2_P4_wspac_beam.py`](test_v2_P4_wspac_beam.py) | 4 | 追跡の WSPAC に自己無撞着なビーム行列を使うオプション |

## 区分 9: SADScript の評価器・パーサ・パッケージ

合計 18 ファイル 262 件。

| ファイル | 件数 | 何を確かめるか |
|---|---|---|
| [`test_evalfixes.py`](test_evalfixes.py) | 25 | SADScript の評価 (連鎖比較・位置指定・スレッド・空列) (SAD 照合) |
| [`test_example_fixes.py`](test_example_fixes.py) | 14 | Packages の遅延読み込み、記号式、有理数演算など (SAD 照合) |
| [`test_scriptgs.py`](test_scriptgs.py) | 17 | Module の字句スコープ、Return、Switch/Which、Out[n] (SAD 照合) |
| [`test_item4.py`](test_item4.py) | 9 | リスト関数、tsvd (SVD ソルバ)、制御構文、CalculateOptics の骨格 |
| [`test_sadclass.py`](test_sadclass.py) | 16 | SAD の Class (継承・メソッド・インスタンス変数) |
| [`test_symbol_names.py`](test_symbol_names.py) | 2 | 裸のシンボルを要素名に使う、要素の型名と同名の LINE |
| [`test_symbolic.py`](test_symbolic.py) | 3 | 記号微分 D とパターン照合 (Orderless・Flat) (SAD 照合) |
| [`test_fileops.py`](test_fileops.py) | 2 | ファイル・ディレクトリ関数、ARRAY・StringInsert など |
| [`test_environment_setenv.py`](test_environment_setenv.py) | 2 | Environment[名前, 値] が環境変数を設定する |
| [`test_miscfn.py`](test_miscfn.py) | 2 | DownhillSimplex・StringToStream など (SAD 照合) |
| [`test_gui_lang.py`](test_gui_lang.py) | 10 | Tkinter.n・KBFrame.n を動かすための言語機能 (リスト・文字列・Apply) |
| [`test_v2_B2_script.py`](test_v2_B2_script.py) | 16 | $Context、With、Visit、Bye など SADScript の小項目 (SAD 照合) |
| [`test_v2_B_script.py`](test_v2_B_script.py) | 25 | PolynomialFit・For・MessageName・複素 FindRoot など |
| [`test_v2_C2_script.py`](test_v2_C2_script.py) | 18 | FFS 文の読み方 (`FFS;NPARA=1;` など) |
| [`test_v2_D2_parser.py`](test_v2_D2_parser.py) | 58 | MAIN 部の要素定義の式 (`ANGLE=PI/8` など) の文法と値 (SAD 照合) |
| [`test_v2_P5_syntax.py`](test_v2_P5_syntax.py) | 16 | MAIN の要素定義の式の文法の細部 (SAD 照合) |
| [`test_v24_A_parser_fit.py`](test_v24_A_parser_fit.py) | 18 | パーサ・評価器・FIT まわり (関数名 At、NaN の FitFunction など) |
| [`test_v2_D2_packages.py`](test_v2_D2_packages.py) | 9 | CorrectOrbit の Micado と Packages の読み込み |

## 区分 10: 特殊関数・数学関数

合計 5 ファイル 88 件。

| ファイル | 件数 | 何を確かめるか |
|---|---|---|
| [`test_specfn.py`](test_specfn.py) | 2 | 特殊関数・数学関数を gs の出力 (specfn_ref.json) と比べる |
| [`test_v24_E_audit.py`](test_v24_E_audit.py) | 42 | 組込の数学・特殊・リスト関数を mpmath の参照値と比べる |
| [`test_v2_B_specfn.py`](test_v2_B_specfn.py) | 13 | Hypergeometric・Bessel・ExpIntegralEi などを独立な計算と比べる |
| [`test_v2_B2_specfn_cx.py`](test_v2_B2_specfn_cx.py) | 8 | 複素引数の Zeta・HurwitzZeta・PolyLog・HypergeometricU |
| [`test_v2_P5_specfn_cx.py`](test_v2_P5_specfn_cx.py) | 23 | 複素引数の特殊関数 (微分方程式の残差・漸化式・恒等式) |

## 区分 11: GUI・グラフィック (Tk・matplotlib・DXF)

合計 6 ファイル 111 件。

| ファイル | 件数 | 何を確かめるか |
|---|---|---|
| [`test_plot.py`](test_plot.py) | 8 | ListPlot・Plot と Tk キャンバスの描画 (SAD 照合。GUI 環境が要る) |
| [`test_tclbridge.py`](test_tclbridge.py) | 7 | Tcl/Tk の接点と Tkinter.n・KBMainFrame.n の動作 (要 `PYSAD_SAD_PACKAGES`) |
| [`test_v241_mpl_display.py`](test_v241_mpl_display.py) | 32 | `$DisplayFunction` = MatplotlibDrawer で Graphics を matplotlib に描く |
| [`test_v25_draw_mpl.py`](test_v25_draw_mpl.py) | 33 | FFS の DRAW (matplotlib) の構文解析と描画データ |
| [`test_dxf_lattice.py`](test_dxf_lattice.py) | 28 | ラティス平面図の DXF の読み戻し・偏向磁石の外形の幾何 |
| [`test_v2_P5_tk90.py`](test_v2_P5_tk90.py) | 3 | Tk 9.0 (Python 3.14 の venv) で GUI が動き窓の大きさが一致 |

## 区分 12: マッチング・FIT・軌道補正

合計 6 ファイル 35 件。

| ファイル | 件数 | 何を確かめるか |
|---|---|---|
| [`test_fitroot.py`](test_fitroot.py) | 2 | FindRoot・Fit・FitGaussian・NIntegrate・SelectCases (SAD 照合) |
| [`test_fam_points.py`](test_fam_points.py) | 6 | FITP の fam 点 (有限振幅ファミリー) の定義と自己検証 |
| [`test_v2_B2_free.py`](test_v2_B2_free.py) | 6 | 成分ごとの FREE (`FREE QF.1 K1;`) (SAD 照合) |
| [`test_v2_P4_free.py`](test_v2_P4_free.py) | 6 | 要素ごとの FREE の「何もしない変数」の数え方 (gs と全桁一致) |
| [`test_v2_P5_micado_fitp_oy.py`](test_v2_P5_micado_fitp_oy.py) | 12 | OYStability・Micado の変種・FITP>1 (SAD 照合) |
| [`test_corr.py`](test_corr.py) | 3 | 軌道補正 CorrectOrbit・Steer・最適化 (SAD 照合) |

## 区分 13: ラティス集・例の確認

合計 4 ファイル 50 件。

| ファイル | 件数 | 何を確かめるか |
|---|---|---|
| [`test_ddba_lattices.py`](test_ddba_lattices.py) | 24 | DDBA ラティス 28 本の `cell calc emit` の値が記録と一致 |
| [`test_ddba_v07_14.py`](test_ddba_v07_14.py) | 5 | DDBA の 1 式のラティス 2 本と誤差計算の縮小版が走る |
| [`test_kek_ls_example.py`](test_kek_ls_example.py) | 8 | KEK-LS の例のラティスの周長・チューン・εx と速い台本 |
| [`test_pf_history_lattices.py`](test_pf_history_lattices.py) | 13 | PF の履歴ラティス 4 本の `cell calc emit` の値 |

## 区分 14: ツール・測定

合計 6 ファイル 70 件。

| ファイル | 件数 | 何を確かめるか |
|---|---|---|
| [`test_tools_md_to_html.py`](test_tools_md_to_html.py) | 21 | md→HTML 変換ツール (tools/md_to_html.py・build_html.py。数式の KaTeX 変換を含む。描画の試験は Node.js がなければ skip) |
| [`test_packaging_cli.py`](test_packaging_cli.py) | 5 | 配布の内容: `pysad-run`・`pysad-share`、文書・例を `share/pysad/` に入れる `setup.py`、版番号の一貫性 |
| [`test_measure.py`](test_measure.py) | 1 | tlum が二つのガウスビームの重なり面積の解析値に近い |
| [`test_v2_C_measure.py`](test_v2_C_measure.py) | 19 | MEASURE・TRCOD・LOGMEASURE と MAIN の TRACK TURNS=-1 |
| [`test_v2_C2_measure.py`](test_v2_C2_measure.py) | 20 | FIXCOD・MVCOD・STNCOD・ALIGN・PUTMEASURE/GETMEASURE |
| [`test_v2_C3_measure.py`](test_v2_C3_measure.py) | 8 | MEASURE・TRCOD・ALIGN・FIXCOD の独立な検証 (大きなリング) |
| [`test_v2_P5_measure.py`](test_v2_P5_measure.py) | 10 | MEASURE・TRCOD の統計を gs と比べる |

件数 0 の 6 本: `test_bend_feeddown.py`・`test_coord.py`・`test_sext_radiation.py` (旧式。上記 `--ignore`)、`test_error_spreading.py`・`test_new_features.py` (スクリプト形式で、import 時に実行される)、`test_emit_extra.py` (空)。

## 部品・データ・scripts/

### 部品 (試験から import される)

| ファイル | 内容 | 使う試験 |
|---|---|---|
| `conftest.py` | 試験の共通設定: 互換モード k64 を既定にし、フィクスチャ v2_mode・k64_mode、Tk 窓・EMITDIV・環境変数の後始末 | 全部 |
| `__init__.py` | 空 (pysad.tests をパッケージにする。相対 import に必要) | - |
| `make_outputs_md.py` | scripts/ の台本と試験関数の対応表 OUTPUTS.md を作る | - |
| `_a2_indep.py` | 放射積分 I1〜I5・U0・Jx・減衰の独立な計算 (3×3 行列) | test_v2_A2_radint ほか |
| `_c3_indep.py` | Liénard-Wiechert 場、四極ずれの Jx、結合リングの平衡の独立な計算 | test_v2_C3_* ほか |
| `_d2_indep.py` | 硬い端の直線の極の厳密な写像 | test_v2_D2_edge |
| `rad_indep_d2.py` | 放射つき 1 周写像 (RK4)+薄い空洞+硬い端の極の独立な基準 | test_v2_D2_edge |
| `_e_spinkick.py` | δ の蹴りに対するスピン不変軸の応答 G を、追跡と線形モデルで測る | test_v2_E_pol_kick |
| `_p1_dk.py` | Derbenev-Kondratenko の式を閉軌道に沿った数値積分で直接評価 | test_v2_P1_st |
| `_p1_models.py` | 偏極の脱偏極率の線形の式と追跡の差を分解する補助 | (補助) |
| `_p2_indep.py` | ローレンツ力の数値積分による硬い端の BEND の厳密な写像 | test_v2_P2_rad |
| `_p3_gsdata.py` | gs の出力値 (3 粒子の x px y py z dp) の保存 | test_v2_P3_* |
| `_p3_indep.py` | 台本の組み立て (FODO リング・大きなリング A+補正磁石) と独立な基準 | test_v241_*, test_v2_P3_* ほか |
| `_p5_gscx.py` | 複素引数の特殊関数の gs の値 (一致点と gs の誤り点) | test_v2_P5_specfn_cx |
| `_p5_util.py` | example の台本と gs の出力を読み、pysad を互換モード指定で走らせる補助 | test_v2_P5_* |
| `_rad_beta_indep.py` | β<1 の放射の独立な計算 (遅延場の時間領域積分) | test_v24_D_radbeta |
| `_rad_indep.py` | ANGLE=0・K0≠0 の極の放射積分の独立な計算 | test_v24_rad_pole |
| `_spin_maxwell_indep.py` | Maxwell を満たす硬い端の扇形 BEND 内の Thomas-BMT 積分 | test_v24_D_spin_bpar |
| `_spin_maxwell_indep2.py` | 直線の磁石 (四極・6 極・MULT) の Maxwell 場内の Thomas-BMT 積分 | test_v24_D3_* |
| `_u1_indep.py` | 一様な磁場指数の扇形 BEND の独立な数値積分 (実験室系) | test_v23_U1/U2 |
| `_u2_indep.py` | ドリフト・四極・扇形 BEND のリングの独立な数値積分 | test_v23_U2 |
| `_v241_lma_tools.py` | 誤差つき LMA と長時間追跡の gs・pysad 比較の道具 (台本生成・gs 実行・解析) | test_v241_lma_long |
| `_v241_lma_errors_run.py` | 誤差つき LMA の gs・pysad の二分法を走らせて scripts/ に保存 | (生成用。直接実行) |
| `_v241_lma_errors_summary.py` | 誤差つき LMA の区間のまとめ (same/adj/diff) | test_v241_lma_long |
| `_v241_lma_k64_run.py` | 保存済みの台本を k64 モードで走らせ `.pysadk64.txt` に書く | (生成用) |
| `_v241_lma_long_run.py` | 10^4 周級の長時間追跡 (FLUC あり/なし) の gs・pysad 比較を走らせる | (生成用) |
| `_v241_lma_long_analyze.py` | 長時間追跡の出力から σδ・σz・ε を n ごとに解析 | test_v241_lma_long |
| `_v241_lma_long_compare.py` | 長時間追跡の gs・pysad (k64・v2) の比較表 | (生成用) |
| `_v241_pol_tools.py` | 偏光 FLUC 追跡の比較の道具 (台本生成・出力解析) | test_v241_pol_fluc |

`__init__.py` は空 (`pysad.tests.script_tests` をパッケージにして相対 import を通す)。

### データ

`kek-ls_v03_12.lat` (大きなリング A: 3 GeV、周長 603 m、3563 要素、20 セル) は、試験・ベンチ・example が共有するので、1 つ上の `pysad/tests/` に置いてある。

| ファイル | 内容 | 大きさ |
|---|---|---|
| `mkbump_ring.sad` | 24 セル・3 GeV の FODO に BPM・ステアラ・四極のずれを入れた誤差入りリング (軌道補正) | 7.6 KB |
| `co5_corr.sad` | 軌道補正の試験用の小リング (3 GeV) | 4 KB |
| `oo1_optimize.sad` | 最適化の試験用の小リング (3 GeV) | 0.6 KB |
| `data_bend_kernels.npz` | K1≠0 BEND の numba 化前の Python 版の出力 (基準) | 32 KB |
| `data_fourier_ref.json` | MAIN の TRACK の FOURIER 解析の K64 Fortran の参照値 | 17 KB |
| `data_tsvd_sad.npz` | SAD の tsvd の出力 (SVD ソルバの基準) | 1.7 KB |
| `specfn_ref.json` | 特殊関数の式と gs の値のリスト | 18 KB |
| `specfn_audit_ref.json` | 特殊関数の mpmath による高精度の参照値 (式・実部・虚部) | 46 KB |
| `symbolic_ref.json` | 記号微分の gs の出力 | 3 KB |

### scripts/ (1100 ファイル)

試験が走らせた SAD 台本と、その画面出力の写し。対応表は [`OUTPUTS.md`](OUTPUTS.md) (`make_outputs_md.py` で作る)。

- 名前: `<試験ファイル名 (test_ を除く)>__<試験関数名>[_2, _3, ...].{sad, gs.txt, pysad.txt}`。`.sad` = 台本、`.gs.txt` = 同じ台本を SAD (K64) で走らせた画面出力、`.pysad.txt` = pysad の画面出力。
- 内訳: txt 737、sad 358、json 5。`__` を持つもの 984 ファイル (48 種の試験ファイル)。持たない 116 ファイルは `v241_*` (偏光 FLUC・誤差つき LMA・10^4 周級の追跡の台本と出力、解析の json) と `v25_draw_gs_*` (FFS の DRAW)。
- 試験が `scripts/` を読むのは `test_v241_lma_long.py`・`test_v241_pol_fluc.py`・`test_v25_draw_mpl.py` の 3 本 (と部品)。ほかの試験の台本は試験の文字列のままで、`scripts/` は取り出した写し。

## 試験の足し方

1. ファイルは `test_<名前>.py`。ファイルの先頭の docstring に目的・条件 (エネルギー、周長など)・基準の出所を書く。
2. リポジトリのルートは 3 段上: `ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', '..'))`。共有ラティスは `os.path.join(HERE, '..', 'kek-ls_v03_12.lat')`。SAD との照合ベンチの基準 (`bench_*.result`) を読むときは `pysad/tests/regression_tests/bench/` (自動判定の対象外は `bench_extra/`) を指し、ファイルが無いときは `skipif` で飛ばす (飛ばした件数が増えていないか確認する)。
3. 他の試験・部品の import は `from . import _p3_indep as X`、`from pysad.tests.script_tests.test_xxx import ...` の形。
4. グローバルの状態: `conftest.py` が互換モード (既定 k64)、`emit6d.EMITDIV`、`PYSAD_*` の環境変数を、各ファイルの終わりに元へ戻し、Tk の窓も閉じる。試験の中でこれらを変えたら、ファイルの中で元に戻す。
5. v2 の動作は `v2_mode` フィクスチャ (k64 は `k64_mode`)。
6. 基準は、(a) 保存した SAD (K64) の出力、(b) 解析式、(c) SAD・pysad のコードを使わない独立な計算 (`_*_indep.py`) のどれか。結果を pysad の出力に合わせて基準を作らない。独立な計算は部品にまとめ、試験は `import` して使う。
7. `scripts/` に写しを置く場合は、上の命名規則で `.sad`・`.gs.txt`・`.pysad.txt` を並べ、`py -3.12 pysad/tests/script_tests/make_outputs_md.py` で対応表を作り直す。台本を `.sad` に取り出す専用のツールはこのリポジトリに無く、取り出し方は未確認。
8. 追加したら、この README の表 (区分・件数) を直す。
