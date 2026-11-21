# example/v2_3: 担当 Q2 の例 (v2.3.0、2026-11-20)

6 極の共鳴駆動項・振幅依存チューンシフト (ADTS)・クロマティシティ (`pysad/resdrive.py`、解説は `docs/49_resonance_terms_adts.md`、記録は `docs/history/v23_Q2.md`)。
(`example/v2_3/README.md` が Q1・Q2 の索引。ここは Q2 だけ。) ラティス: 小さな FODO (`example/00_lattice/FODO_ring.sad`、SF・SD に K2L = 4, −6)、ASTRID2 (`example/astrid/lattice/ASTRID2_v01_01.sad`、0.58 GeV、45.7 m)、
大きなリング A (`pysad/tests/kek-ls_v03_12.lat`、3 GeV、603 m、3563 要素)。誤差なし。条件: δ = 0、結合なし、RF・放射なし、6 極 2 次 + 8 極 1 次、解析式 (厚いレンズ) が既定。

| ファイル | 内容 |
|---|---|
| `Q2_resdrive_fodo.sad` / `.pysad.txt` | SAD の台本から `ResonanceDrivingTerms[]`・`AmplitudeTuneShift[]`・`ChromaticityTerms[]` を呼ぶ (`py -3.12 run_sad.py example/v2_3/Q2_resdrive_fodo.sad`)。gs には無い関数 |
| `Q2_rings.py` / `Q2_rings.pysad.txt` / `Q2_rings_elements.png` | 3 つのリングの共鳴駆動項 10 個 (要素ごとの寄与の上位)、ADTS (6 極・8 極、調和ごと、数値積分との照合、δ 依存)、クロマティシティ (数値微分と 6 極の 1 次の寄与)、作動点の確認。図: |h30000| と ADTS の累積寄与 |
| `Q2_thick_thin.py` / `Q2_thick_thin.pysad.txt` | **厚いレンズと薄レンズ列の収束の表**: n 分割 (1〜64) との差 (ADTS と全 10 項。n を倍にすると 1/4)、数値積分 (nslice = n) との一致、長さ → 0 の極限 |
| `Q2_adts_tracking.py` / `Q2_adts_tracking_{fodo,astrid,kekls}.pysad.txt` | pysad の追跡 (k64) の NAFF から ADTS を求めて解析値と比較 (6 極の強さ ×(+1, −1, 0) の C、8 極 ×(±1) の B) |
| `Q2_adts_gs.py` / `Q2_adts_gs_{fodo,astrid,kekls,kekls_oct0}.txt` | gs (SAD K64、WSL) の TrackParticles を 1 周ずつ回し (z は毎周 0 に戻す。R1 で台本に `NORFSW;` を足した。旧版の gs の出力 `Q2_adts_gs_*.txt` は RF ON のままなので ×0・×−1 が食い違う。RF なしの出力は `R1_adts_gs_*_norfsw.txt`)、同じ解析。`py -3.12 example/v2_3/Q2_adts_gs.py fodo\|astrid\|kekls [--turns 1024] [--cases 1:1,-1:1,0:1,0:-1] [--jb J]`。gs は timeout 付き、終了後に sad1.exe の残りを確認 |

結果の要約 (J はアクション、単位 1/(m rad)。数値の根拠・条件は上の各出力と docs/49 §6):

| リング | 解析 6 極 2 次 (xx, xy, yy) | pysad 追跡 C | gs 追跡 C |
|---|---|---|---|
| FODO | 323.41, −474.02, 496.07 | 323.42, −474.03, 495.96 | 323.42, −474.03, 495.96 |
| ASTRID2 | 440.30, 1294.57, 1127.55 | 440.97, 1296.00, 1125.57 (J = 2e-8) | 441.0, 1296.0, 1125.6 (NORFSW。RF ON の旧出力は 452.7, 1318.5, 1168.1) |
| 大きなリング A | 86841, −75099, −1610 | 86842, −75097, −1608 | 86842, −75097, −1608 (NORFSW。RF ON の旧出力は 121125, −46333, 22196) |

大きなリング A の 8 極 1 次: 解析 (−36905, 32480, −7154)、pysad 追跡 (−36905, 32480, −7154)、gs (−36890, 32480, −7154)。
**gs の差について (解決、R1)**: 旧版の gs の追跡は FFS の既定の RFSW ON のままで、RF つきのラティス (ASTRID2・大きなリング A) の空洞が 1 周内の z に応じて δ を与えていた (自然クロマティシティ × δ(J) の偽の ADTS。6 極 ×1 は ξ ≈ 0 で見えない)。`NORFSW;` を入れると gs は pysad (k64) と解析に一致する。`example/v2_3/R1_README.md`、`docs/history/v23_R1.md`。
