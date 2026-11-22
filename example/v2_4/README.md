# example/v2_4: v2.4.0 の解析式・LMA・Touschek の検証 (A1〜A7)

記録は [`docs/history/v24_A.md`](../../docs/history/v24_A.md)、解説は [docs/49 §13](../../docs/49_resonance_terms_adts.md) (ampnf)・[docs/48 §4.4](../../docs/48_frequency_map_lma.md) (Touschek のバンチ内の δ)、入門は [docs/50 §7.7](../../docs/50_adts_rdt_primer.md)。
リポジトリ直下から `PYTHONUTF8=1 PYTHONIOENCODING=utf-8` で実行する。環境変数 `PYSAD_SAD_PACKAGES` (K64 版 SAD の Packages のフォルダ) は LMA (A4) に要る。

| 台本 | 内容 | 時間 | 出力 |
|---|---|---|---|
| [`A1_ampnf_vs_tracking.py`](A1_ampnf_vs_tracking.py) `fodo\|astrid` | [`pysad/ampnf.py`](../../pysad/ampnf.py) (1 周の写像のべき級数の正規形) の ν(Jx, Jy, δ) を、実際の追跡 (pysad k64、NAFF、δ = 0・±1e-3) と比べる。追跡の初期座標の J を正規形の J に直してから級数を評価 | FODO 1〜2 分、ASTRID2 1〜2 分 (D = 6、bend_cubic・kin・fringe) | `A1_fodo.out.txt`・`A1_astrid.out.txt` |
| [`A2_k0_bend_check.py`](A2_k0_bend_check.py) `[variant ...]` | K0≠0 の偏向磁石 (ANGLE≠0・E1・E2・K1・ROTATE) の共鳴駆動項・ADTS・ξ を、追跡の 1 周の写像のテイラー展開と比べる (FODO の 1 セルに 5 台のブロック)。変種 `k1_k00_leak`・`k1_kq`・`k1_k2` は、K0 と K1 の組み合わせのずれが閉軌道の漏れによることの切り分け | 7 秒 | 標準出力 ([`v24_A.md`](../../docs/history/v24_A.md) A2・A12 の表) |
| [`G_k0_ampnf.py`](G_k0_ampnf.py) `[mpw\|mpw_big\|mpw_ang\|mpw_edge\|tilt\|pi\|pi_neg]` | K0≠0 の偏向磁石 (3 極の挿入光源 K0 = a・(1, −2, 1) rad、ANGLE = 0・≠0・E1/E2・傾き・ROTATE = π) の ν(J, δ)・J¹・h を ampnf (TPSA、`pysad/ampnf_k0.py`) で求め、実際の追跡 + NAFF・追跡の 2 次のテイラー展開と比べる (D = 5) | 1 通り 2〜3 分 | [`G_k0_ampnf.out.txt`](G_k0_ampnf.out.txt) ([`v24_G.md`](../../docs/history/v24_G.md)) |
| [`G_independent_flow.py`](G_independent_flow.py) `[v2\|k64]` | K0≠0 の偏向磁石の K64 の追跡の式 (`tbendcore1`・`k_steer`) を、K64 の式を使わない独立な基準 (実験室の直交座標の荷電粒子の運動方程式の数値積分。Maxwell を満たす端の幅 w の場) と比べる。水平面 1e-12、Maxwell の項は w → 0 で一致 | 約 2 分 | [`G_independent_flow.out.txt`](G_independent_flow.out.txt) |
| [`G_real_scale.py`](G_real_scale.py) `[astrid\|astrid_tilt\|kekls] [D]` | 実機規模の K0≠0 (ASTRID2 + 2 T MPW 23 極、傾けた MPW、KEK-LS の 1 セルに MPW 23 極) の ν(J, δ)・J¹ を ampnf で求め、実際の追跡 + NAFF と比べる | ASTRID2 10 s (D = 5)、KEK-LS 144 s + 追跡 | [`G_real_scale.out.txt`](G_real_scale.out.txt) |
| [`A3_mc_bunch.py`](A3_mc_bunch.py) `astrid2\|kekls` | Touschek のバンチ内の δ ([`pysad/touschek_bunch.py`](../../pysad/touschek_bunch.py)) の独立な検証: 実験室系の 2 粒子の対のモンテカルロで、絶対値の判定 / kick の判定の損失率の比を、δ_c・q_s の内訳つきで数える | 数分〜15 分 (`gather` が Emittance[] を走らせる) | `A3_mc_bunch_astrid2.out.txt`・`A3_mc_bunch_kekls.out.txt` |
| [`A3_qs_3d.py`](A3_qs_3d.py) `astrid2\|kekls` | q_s の効果を入れた 3 次元の積分 (`touschek_bunch.qs_ratio`) と、上のモンテカルロの比較 | 数分 | `A3_qs_3d_astrid2.out.txt`・`A3_qs_3d_kekls.out.txt` |
| [`A3_qs_ring.py`](A3_qs_ring.py) | ASTRID2 全周の損失率の増え (δ_c だけ / δ_c + q_s を 6・12・24 スライスで補間) | 約 15 分 | `A3_qs_ring.out.txt` |
| [`A4_lma_support.py`](A4_lma_support.py) `k64\|coupl\|kekls [workers] [level] [nloc]` | LMA の k64 (v2 との比較)、強い結合 (誤差 level 5)、大きなリング A の誤差つき | 3〜6 分 (ASTRID2 8 位置、workers 2)、KEK-LS はさらに長い | `A4_k64.out.txt`・`A4_coupl.out.txt`・`A4_kekls.out.txt`・`A4_kekls_level1.out.txt` |
| [`A5_nu_of_J_grid.py`](A5_nu_of_J_grid.py) `astrid2\|kekls [D] [skip_den]` | 周波数マップの保存済みの表 ([`example/v2_3/Q1_fmap_*.table.txt`](../v2_3/README.md)) の各点のチューンを ampnf の ν(J) (J¹ のみ / 級数全体) と比べる。ロックした点の識別 (差の共鳴の島の位置予想) と、最小の共鳴の分母の項の影響 (`skip_den`) | ASTRID2 D = 6 で 1 分、D = 8 で 4〜5 分、KEK-LS D = 5 で 3 分 | `A5_astrid2.out.txt`・`A5_astrid2_D8.out.txt`・`A5_astrid2_D8_skip.out.txt`・`A5_kekls.out.txt` |
| [`A7_ampnf_kekls.py`](A7_ampnf_kekls.py) `[D]` | 大きなリング A (KEK-LS) の ampnf の全部入りの計測と、`adts().total`・CALC との比較 | D = 5 で約 2 分 | `A7_kekls.out.txt` |

DDBA v07_14 の 6 極検討への J²・δ 依存の追加: [`example/ddba/v07_14/sextupole/sextupole_j2_check.py`](../ddba/v07_14/sextupole/sextupole_j2_check.py) (`results/j2_correlation.out.txt`)。DDBA v07_14 の誤差の計算の DA (`ExpandElementValues->False`) の出力は [`example/ddba/v07_14/errors/`](../ddba/v07_14/errors/) の `*.expandfalse.*.txt`。

条件の要約: ampnf の模型 = 線形 + 薄い・厚い 6 極・8 極 + (オプション) 偏向磁石の曲率の 3 次の項と端の生成関数 + 運動学的な項 + 四極・K1≠0 の偏向磁石の端部のフリンジ (強さ K1/(1+δ))。入れていない: K0≠0・ROTATE の偏向磁石、結合、誤差、K1 の偏向磁石の 4 次以上、K64 の四極の積分器の台形則。
追跡は pysad の k64 (RF・放射なし、結合なし、1024 周、NAFF の Hann 窓)。LMA は v2 (RFSW・RAD ON、FLUC なし、LOSSAMPL ±10 mm の例の値、バケット判定あり)。

## H (結合した経路の RDT の拡張、KEK-LS の補正つき LMA。docs/49 §14.5・docs/48 §11、記録 `docs/history/v24_H.md`)

| スクリプト | 内容 | 所要時間 | 出力 |
|---|---|---|---|
| [`H1_coupled_h4.py`](H1_coupled_h4.py) `fodo\|astrid` | ロール・ずれ (DX・DY・DROTATE)・δ≠0・始点のある lat の 3 次・4 次の h を、模型の 1 周の写像 ([`pysad/ampnf_coupled.py`](../../pysad/ampnf_coupled.py)) と追跡の写像のテイラー展開から比べる。ADTS は h3・h4 の正規形・ampnf・追跡 + NAFF | FODO 約 3 分、ASTRID2 約 2 分 | `H1_fodo.out.txt`・`H1_astrid.out.txt` |
| [`H2_run.py`](H2_run.py) `種 ctrl niter ae fe te [workers [turns [lat]]]` | KEK-LS (v03_68・v03_12) に誤差を与え、COD を補正し、補正後の状態で `MomentumAcceptance[]` を台本から走らせる (台本 [`H2_kekls_lma_corrected.template.sad`](H2_kekls_lma_corrected.template.sad)・[`H2b_kekls_v03_12_lma_corrected.template.sad`](H2b_kekls_v03_12_lma_corrected.template.sad))。ctrl = 1 は補正後に `cell calc expand` (既定の CALEXP の対照)、niter = 0 は補正なし | 台本全体で 3〜5 分 (LMA 4 位置、workers 4) | `H2_kekls_*.out.txt` (v03_68)・`H2_kekls12_*.out.txt` (v03_12) |

条件: H1 は v2 モード (DX・DROTATE のある BEND の tchge は v2 の符号)、RF・放射なし。H2 は v2、RFSW・RAD ON (FLUC なし)、LOSSAMPL ±10 mm の例の値、周回数は自動 (シンクロトロン振動 2 周期)、バケット判定あり。

