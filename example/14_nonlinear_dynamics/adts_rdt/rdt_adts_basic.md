# 6 極の共鳴駆動項・ADTS・クロマティシティの例

**この文書で分かること**: 6 極 (と 8 極) が作る共鳴駆動項 (RDT)、振幅依存チューンシフト (ADTS)、クロマティシティを、pysad の解析式 (`pysad/resdrive.py`) で 3 つのリングについて計算する例と、その解析値が追跡 + NAFF (チューンの高精度な周波数解析) から求めた ADTS と一致する範囲。
**前提**: 加速器物理 (ベータ関数、チューン、6 極の働き) を知っていること。pysad・SAD の使い方は知らなくてよい (コマンドは下にそのまま書いてある)。
**読む順**: 先頭の表 (何があるか) → 結果の要約の表 → 詳しい式・数値は `ADTS_summary.md`。解説は `docs/25_resonance_terms_adts.md`。この節の例の索引は `example/14_nonlinear_dynamics/README.md`。

用語: **共鳴駆動項 (RDT)** h_jklm は、6 極などの非線形磁石が共鳴 (例 3νx = 整数) をどれだけ駆動するかを表す複素数。**ADTS** は振幅 (アクション J) が大きいほどチューンが変わる量 ∂ν/∂J (単位 1/(m rad)、成分は xx = ∂νx/∂Jx、xy = ∂νx/∂Jy、yy = ∂νy/∂Jy)。J は x = √(2Jβ) cos φ で決まるアクション (単位 m rad)。**6 極 2 次**は、6 極の強さの 2 乗に比例して ADTS に効く寄与。

## 条件

| リング | ラティス | エネルギー・周長 |
|---|---|---|
| 小さな FODO | `example/00_lattice/FODO_ring.sad` (SF・SD に K2L = 4, −6) | 3 GeV、105.6 m |
| ASTRID2 | `example/astrid/lattice/ASTRID2_v01_01.sad` | 0.58 GeV、45.7 m |
| 大きなリング A | `pysad/tests/kek-ls_v03_12.lat` (3563 要素) | 3 GeV、603 m |

誤差なし。δ = 0、結合なし、RF・放射なし。ADTS は 6 極 2 次 + 8 極 1 次。解析式は厚いレンズ (数値積分でなく閉じた式) が既定。

## ファイル

| ファイル | 内容 |
|---|---|
| `fn_ResonanceDrivingTerms_fodo.sad` / `.pysad.txt` | SAD の台本から `ResonanceDrivingTerms[]`・`AmplitudeTuneShift[]`・`ChromaticityTerms[]` を呼ぶ (`py -3.12 run_sad.py example/14_nonlinear_dynamics/adts_rdt/fn_ResonanceDrivingTerms_fodo.sad`)。この 3 関数は pysad が加えたもので、元の SAD (gs) には無い |
| `rdt_rings.py` / `rdt_rings.pysad.txt` / `rdt_rings_elements.png` | 3 つのリングの共鳴駆動項 10 個 (要素ごとの寄与の上位)、ADTS (6 極・8 極、調和ごと、数値積分との照合、δ 依存)、クロマティシティ (数値微分と 6 極の 1 次の寄与)、作動点の確認。図は |h30000| と ADTS の累積寄与 |
| `rdt_thick_thin.py` / `rdt_thick_thin.pysad.txt` | **厚いレンズと薄レンズ列の収束の表**: 6 極を n 分割 (1〜64) した薄レンズ列との差 (ADTS と全 10 項。n を倍にすると差が 1/4)、数値積分 (nslice = n) との一致、長さ → 0 の極限 |
| `adts_tracking.py` / `adts_tracking_{fodo,astrid,kekls}.pysad.txt` | pysad の追跡 (互換モード k64) の NAFF から ADTS を求めて解析値と比較 (6 極の強さ ×(+1, −1, 0) から 6 極 2 次の係数 C、8 極 ×(±1) から 8 極 1 次の係数 B) |
| `adts_gs.py` / `adts_gs_{fodo,astrid,kekls,kekls_oct0}.txt` | SAD 本家 (gs、WSL) の TrackParticles で同じ解析を行う (検証用。`NORFSW;` で RF を止める。RF を止めない場合は `adts_gs_norfsw.md`) |

実行: `py -3.12 example/14_nonlinear_dynamics/adts_rdt/rdt_rings.py` (約 40 s)、`... adts_tracking.py fodo|astrid|kekls`、`... rdt_thick_thin.py`、gs は `... adts_gs.py fodo|astrid|kekls [--turns 1024] [--cases 1:1,-1:1,0:1,0:-1] [--jb J]` (gs は timeout 付き、終了後に WSL の sad1.exe が残っていないか確認する)。

## 結果の要約

6 極 2 次の ADTS 係数 C (J はアクション、単位 1/(m rad)、成分は (xx, xy, yy))。解析 = `resdrive` の厚い 6 極の閉じた式、追跡 = pysad の追跡 + NAFF (J = 2e-8 付近)。

| リング | 解析 6 極 2 次 (xx, xy, yy) | pysad 追跡 C |
|---|---|---|
| FODO | 323.41, −474.02, 496.07 | 323.42, −474.03, 495.96 |
| ASTRID2 | 440.30, 1294.57, 1127.55 | 440.97, 1296.00, 1125.57 (J = 2e-8) |
| 大きなリング A | 86841, −75099, −1610 | 86842, −75097, −1608 |

読み方: 6 極 2 次の ADTS は、解析式と追跡 + NAFF で FODO は 3〜4 桁、ASTRID2 で 0.2% 程度、大きなリング A で 1e-4 程度一致する。ASTRID2 の 1e-3 の差は、追跡の J = 2e-8 での J の 2 次以降の寄与による (J → 0 で近づく)。大きなリング A の 8 極 1 次は解析 (−36905, 32480, −7154) と追跡 (−36905, 32480, −7154) が一致する。

検証: gs (SAD の本家) の TrackParticles でも、RF を止めた条件 (`NORFSW;`) で同じ C が得られる (FODO・ASTRID2・大きなリング A で pysad と 1e-4 以内。詳細は `adts_gs_norfsw.md`)。gs の追跡を RF つきのまま (FFS の既定) にすると、RF つきのリングでは偽の ADTS が入る。
