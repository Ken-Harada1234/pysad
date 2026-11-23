# 32 SAD ソース解析: gs との照合の桁 (旧 docs/19 §4 の「照合桁」の列)

旧 docs/19_source_analysis.md の冒頭と §4 から、gs (SAD K64 の WSL ビルド) との照合の桁をそのまま移したもの。本文 (docs/32_source_analysis.md §4) には、対応表と限界だけを置いている。
「照合桁」は SAD (K64) の印字精度 (10 桁、追跡は 5e-11) との最大差で、ベンチ名を併記する。ベンチの一覧は docs/27_benchmarks.md (旧 10_benchmarks)。

| コマンド / 関数 | 照合桁 (旧 19 §4 の原文) |
|---|---|
| `TrackParticles[]` | ≤5e-11 (Level 6〜11)、FLUC の乱数状態一致。スピン追跡は 4e-12 以下、RADLIGHT の軌道は差 0 (2026-11-16) |
| FFS CALC (4D) | 閉軌道 9-10 桁、Twiss ≤1e-10 (誤差・オフモメンタム付きは ≤3e-7、10_benchmarks §3.1) |
| FFS CALC6D | 1.7e-10 (Level 18a/b/c) |
| FFS EMIT / `Emittance[]` | 全量 ≤1.3e-11 (Level 16B)、OneTurn/Matrix 系 2.1e-10 (Level 21a)。POL の MULT・空洞・SOL は gs と 1e-11 (test_spin.py) |
| FFS RADINT | 7 桁 (Level 29a/29a2) |
| FFS CHRO / 色収差 | チューン ≤2.5e-7 (Level 3A, 14A) |
| FFS GO (FIT) | 反復数・残差・K1 一致 (Level 13)。解析版は SAD の反復列と K1 を 12 桁 |
| `CalculateOptics[]` | — |
| `DynamicApertureSurvey` | NORAD は地図・Score 完全一致 (Level 16A)。RAD+FLUC は損失周の乱数が未一致 |
| `TouschekLifetime` | ≤7e-12 (Level 17a)、RADCOD ON は 1.7e-8 (17b) |
| ジオメトリ (GX/GY/GZ/CHI) | 4.8e-13 (Level 19a/b/c) |
| `SynchroBetaEmittance` | gs と相対 1e-5 (2026-11-16) |
| `RadiationField` / `RadiationSpectrum` | RadiationField は gs が SIGSEGV (B49) で、等速運動の場と相対 6e-3。軌道は gs と差 0 |
| CSR 系 | 11 桁 (CSRTrack・CSRHaissin)、1e-9〜1e-13 (CSRSetupOY・CSRMatrix・CSRInit)。CSROYMatrix は gs が wrongtype (B60) で未照合 |
| `GaussianCoulomb(Fitted)` / BBBrem | docs/32 (旧。現在は docs/06・14) |
| DELCOR | チューン 1e-14、軌道 1e-10 (docs/31。現在は docs/04) |
| MAIN TRACK | SAD K64 の MAIN TRACK 自体は tinip で SIGSEGV のため、SADScript で同じ手順を書いた版 (bench_level27a/b) と照合 |
