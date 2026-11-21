# R2: ADTS の直接の検証 (追跡 + NAFF + 多項式フィット) と未確認の追及 (v2.3.0、2026-11-20)

記録は `docs/history/v23_R2.md`、解説は `docs/49_resonance_terms_adts.md` (§7)。ラティス: FODO (`example/00_lattice/FODO_ring.sad`、SF・SD の K2L = 4, −6、3 GeV、105.6 m)、ASTRID2 (`example/astrid/lattice/ASTRID2_v01_01.sad`、0.58 GeV、45.7 m)、ASTRID2 MPW (`ASTRID2_v01_02_id.sad` の `RINGMPW`)、大きなリング A (`pysad/tests/kek-ls_v03_12.lat`、3 GeV、603 m)。
条件 (すべて): RF 空洞は止める (pysad は rfsw=False、gs は `NORFSW;`)、放射なし、結合なし、δ=0 (δ の節を除く)、閉軌道 0、k64、1024 周、NAFF は Laskar の Hann 窓 (p=1) で線形のチューン ± 0.02 を探す。
J の格子: Jx, Jy ∈ {0, 1/5, …, 1}·Jmax の 35 点 × (φ, φ+π) の対 = 70 粒子 (Jmax = FODO 1e-7、ASTRID2 2e-8、大きなリング A 2e-10)。J はアクション (x = √(2Jβ) cos φ、単位 m rad、リング先頭)。フィットは全次数 ≤ 2 の多項式 (order=2)。

| ファイル | 内容 |
|---|---|
| `R2_adts_direct.py` | 直接の比較の本体。`pysad|gs|plot リング` (`--study` でフィットの次数・Jmax・NAFF の設定の走査)。gs は WSL、timeout 1800 s、1 本ずつ、終了後に sad1.exe の残りを確認 |
| `R2_adts_direct_{fodo,astrid,astrid_mpw,kekls}.pysad.txt`・`.json` | pysad (k64) の追跡のフィット、解析 (6 極 2 次・8 極 1 次・運動学・端部のフリンジ・合計) との表、パラキシャルな追跡 (解析と同じ模型) との照合、`--study` の走査 |
| `R2_adts_direct_{fodo,astrid,astrid_mpw,kekls}.gs.txt`・`.json` | gs の追跡のフィット (pysad と同じ粒子の差も) |
| `R2_fit_{fodo,astrid,astrid_mpw,kekls}.png` | νx−ν0 vs Jx (Jy=0)、νy−ν0 vs Jy (Jx=0)、νx−ν0 vs Jy (Jx=Jmax) と 3 次のフィット・解析の接線・残差 |
| `R2_unverified.py` | 未確認の節: `octupole`・`ablate`・`delta`・`resonance`・`coupling` |
| `R2_octupole_kekls.pysad.txt` | 8 極の強さ μ = 0〜8 倍の走査 (1 次の傾きが解析の 8 極 1 次と 2e-6〜2e-5、2 次の係数の μ 依存) |
| `R2_ablate_{fodo,astrid,kekls}.pysad.txt` | 端部のフリンジ (DISFRIN=1)・6 極・8 極を切った追跡と、解析 (切った分を除く) の比較 |
| `R2_delta_{fodo,astrid,kekls}.pysad.txt` | チューンの δ 依存 (追跡 + NAFF と CALC・6 極の 1 次の式)、オフモメンタムの ADTS (δ = ±δ0) |
| `R2_resonance.pysad.txt` | FODO の 3νx 共鳴線の近くの周期点の振幅: h30000 + ADTS の予測と追跡 (Newton 法) |
| `R2_coupling.pysad.txt` | ASTRID2 に四極の回転・DX・DY を入れた追跡の ADTS (結合・垂直分散・閉軌道を無視した解析との差) |

結果の要約 (1/(m rad)、(xx, xy, yy)。詳細は `docs/history/v23_R2.md` の表):

| リング | 解析の合計 | 追跡 pysad (k64) = gs | 相対差 |
|---|---|---|---|
| FODO | 397.94, −405.98, 569.07 | 408.99, −412.74, 519.95 | +2.8%, −1.7%, −8.6% |
| ASTRID2 | 522.48, 1504.32, 1388.64 | 470.96, 1304.81, 1485.80 | −9.9%, −13%, +7.0% |
| ASTRID2 MPW | 522.45, 1737.43, 1578.19 | 470.96, 1509.42, 1712.28 (gs 1712.26) | −9.9%, −13%, +8.5% |
| 大きなリング A | 50807.3, −40553.9, −6295.2 | 50709.5, −41042.5, −6527.8 | −0.19%, −1.2%, −3.7% |

解析の合計 = 6 極 2 次 + 8 極 1 次 + 運動学的な項 + 四極の端部のフリンジ (`resdrive.adts` の `.total`)。残りの差は 6 極 × (偏向磁石の曲率・硬い端) の交差項 (未実装)。gs と pysad の差は 1e-10 (ASTRID2 MPW は 1e-5)。
再現: `py -3.12 example/v2_3/R2_adts_direct.py pysad fodo --study`、`... gs astrid` (約 2 分)、`py -3.12 example/v2_3/R2_unverified.py ablate kekls` など。
