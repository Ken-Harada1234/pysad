# K1≠0 の BEND の K1 の定義のオプションと CALC の dk の修正の例

試験は `pysad/tests/test_v23_U2_bend_k1_def.py`。すべて v2 モード。リポジトリ直下から `PYTHONUTF8=1 py -3.12 example/15_lma_touschek/touschek_elegant/bend_k1_def_xxx.py` (標準出力をファイルにリダイレクトするには `PYTHONUTF8=1` が要る)。

## 使い方

```python
from pysad import compat
compat.set_bend_k1_def('lin')     # B' が x に一様 (elegant CSBEND K2 = 0・MAD SBEND と同じ)。既定は 'log' (K64: (1 + hx)∂b/∂x 一様)
# 環境変数: PYSAD_BEND_K1_DEF=lin
```

v2 だけで有効 (k64 では無視)。追跡・CALC・`Emittance[]`・MULT (ANGLE≠0) の K1 に効く。'lin' = 'log' + 6 極 K2 = +hK1 (kick のポテンシャル K1 h (x³/3 − x y²/2))。δ = 0 の設計軌道の光学・行列は変わらない。elegant に変換するとき、'log' (K64 の定義) は K2 = −hK1 を足し、'lin' は K2 = 0 のまま。

## ファイル

| ファイル | 内容 |
|---|---|
| `bend_k1_def_verify.py` / `bend_k1_def_verify.out.txt` | (1) n = 0.4 の扇形磁石 16 台 (3 GeV 相当) の dν/dδ (CALC・追跡・独立な積分・解析式・elegant の保存値)、5 つの δ の追跡と積分のチューンの差、写像の差。(2) ASTRID2 の自然な色収差を独立な積分と比較。(3) CALC の dk の修正の影響 (n = 0.4 のリング、ASTRID2、大きなリング A) |
| `bend_k1_def_astrid2_chrom.py` / `bend_k1_def_astrid2_chrom.out.txt` | ASTRID2 の ξ・∂²ν/∂δ² を pysad ('log'・'lin'、フリンジあり・なし、6 極あり・なし) と elegant (K2 = −hK1・K2 = 0、EDGE_ORDER = 1・2) で。elegant は WSL (`~/elegant/env.sh`、serial、timeout 280 s、終了後 `ps` で残りなしを確認)。`sad2ele.py` の `astrid2_all.lte` の K2・EDGE_ORDER・6 極の K2 を書き換える |
| `bend_k1_def_lma_astrid2.py` / `bend_k1_def_lma_astrid2.out.txt` | ASTRID2 の 1 セル (25 位置) の LMA を 'log'・'lin' で (RF・放射つき、アパーチャ ±10 mm、4 workers、約 9 分) |
| `pysad/tests/_u1_indep.py`・`_u2_indep.py` | 独立な数値積分 (実験室系のローレンツ力、DOP853。pysad・SAD のコードを使わない)。`_u2_indep` はドリフト・厚い四極・複合機能の扇形磁石・端の縦方向の標準の集束 |

elegant との比較は保存値 (`bend_k1_elegant.out.txt`: n = 0.4 のリング) と、`bend_k1_def_astrid2_chrom.py` が走らせた elegant (ASTRID2) を使った。

## 結果の要約

n = 0.4 の扇形磁石 16 台 (3 GeV 相当、K1 = −0.061685 /m²、ρ0 = 2.5465 m):

| 量 | 'log' CALC | 'lin' CALC | 'lin' の積分 | 解析式・elegant K2 = 0 | 'log' の積分・gs・elegant K2 = −hK1 |
|---|---|---|---|---|---|
| dνx/dδ | −0.17213 | −0.60246 | −0.60251 | −0.60246 | −0.17213 |
| dνy/dδ | +0.21082 | +0.73786 | +0.73789 | +0.73786 | +0.21082 |

5 つの δ (±0.4%) のチューン (追跡・CALC) − 独立な積分: x 1.2e-8 以内、y 1e-11 以内。

ASTRID2 の自然な ξ (6 極なし): pysad 'log' (−6.2181、−11.8060) = elegant K2 = −hK1、pysad 'lin' (−6.2486、−11.0336) = elegant K2 = 0 (EDGE_ORDER = 2)、いずれも 4 桁以上一致。

∂²ν/∂δ² (ASTRID2、6 極つき、pysad 'lin' / elegant K2 = 0): x 76.92 / 71.53、y 79.70 / 27.50 (2.9 倍)。'lin' では y の差は解消しない。原因は pysad の K64 のフリンジと、elegant の `dnu/dp2` の精度 (K1 の偏向磁石のない FODO でも y が 23% 外れる)。独立な積分 (ASTRID2 自然、フリンジなし) は pysad と y 4.5% 以内。

LMA (ASTRID2、1 セル 25 位置): δ+ の差 4e-12、δ− の差は最大 2.6e-5 (相対 1.6e-3)。K1 の定義 (`log`・`lin`) に鈍感。

CALC の dk (`compat.set_bend_dk_fix`): v2 は独立な積分と一致する。n = 0.4 のリングで δ = 2% の ν_x の独立な積分との差は 7.5e-8 (k64 は 1.3e-4)。ASTRID2・大きなリング A では ξ に 4e-6、ν に 5e-7 以内。
