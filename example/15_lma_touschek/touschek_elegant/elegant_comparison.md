# elegant との比較の例 (局所運動量アクセプタンスと Touschek 寿命)

試験は `pysad/tests/script_tests/test_v23_T1_elegant.py`。この例は比較と、SAD → elegant の変換のスクリプトである。
実行はリポジトリ直下から `PYTHONUTF8=1 py -3.12 example/15_lma_touschek/touschek_elegant/xxx.py ...` (標準出力をファイルにリダイレクトするには `PYTHONUTF8=1` が要る)。すべて v2 モード。
elegant は WSL (Ubuntu、`~/elegant/env.sh`、elegant 2026.4.0)。使い方は elegant の公式文書。CPU は最大 4 コア (serial の elegant を `xargs -P 4`、または `mpirun -np 4 Pelegant`)。
gs (K64、2026-10-19 ビルド、WSL、NPARA = 1) は `gs_bend_k1.py` だけ (timeout 600 s、終了後に `ps` で sad1.exe の残りがないことを確かめた)。
使ったラティス: `example/astrid/lattice/ASTRID2_v01_01.sad` (ASTRID2、0.58 GeV、45.7 m、6 セル)、`pysad/tests/kek-ls_v03_12.lat` (3 GeV、603 m、20 セル)。

## 走らせるコマンド

```
# 1. pysad 側の基準値 (光学、IBS 込みのビーム、一様な受容の Touschek 寿命) → pysad_<名前>.ref.json
PYTHONUTF8=1 py -3.12 example/15_lma_touschek/touschek_elegant/pysad_ref.py astrid2
PYTHONUTF8=1 py -3.12 example/15_lma_touschek/touschek_elegant/pysad_ref.py kekls              # pysad の LMA の表 (3.) があれば pysad_local (局所の寿命) も入る
# 2. SAD → elegant の変換 (lte と要素の対応表 .map.json)。--u0 は RF の位相用 (省略なら pysad の値)
PYTHONUTF8=1 py -3.12 example/15_lma_touschek/touschek_elegant/sad2ele.py astrid2 example/15_lma_touschek/touschek_elegant
# 3. LMA: pysad (kek-ls は elegant と同じ位置で走らせる。ASTRID2 は保存済みの表 (`lma_astrid2.table.txt`) を使う) と elegant → 比較の表・図
PYTHONUTF8=1 py -3.12 example/15_lma_touschek/lma/lma_compare.py kekls --pysad --workers 2     # pysad の LMA を pysad_kekls.lma.txt に (52 位置、約 23 分)
PYTHONUTF8=1 py -3.12 example/15_lma_touschek/touschek_elegant/elegant_lma.py astrid2 --mpi --np 4            # elegant の LMA、全周 199 位置、Pelegant (約 9 分)
PYTHONUTF8=1 py -3.12 example/15_lma_touschek/touschek_elegant/elegant_lma.py kekls --np 4                    # 1500 周、14 位置、serial の elegant を位置ごとに 4 本同時
PYTHONUTF8=1 py -3.12 example/15_lma_touschek/touschek_elegant/elegant_lma.py astrid2 --mpi --np 4 --exact-k1 --cell   # 感度: K2 = −h K1 を足さない厳密な勾配 (1 セル)
PYTHONUTF8=1 py -3.12 example/15_lma_touschek/lma/lma_compare.py astrid2                         # 比較 (表 lma_compare_astrid2.table.txt、図)
PYTHONUTF8=1 py -3.12 example/15_lma_touschek/lma/lma_compare.py kekls
# 4. Touschek 寿命: elegant の touschekLifetime (一様な受容、elegant の LMA、pysad の LMA) → elegant_<名前>.touschek.txt
PYTHONUTF8=1 py -3.12 example/15_lma_touschek/touschek_elegant/elegant_touschek.py astrid2
PYTHONUTF8=1 py -3.12 example/15_lma_touschek/touschek_elegant/elegant_touschek.py kekls
PYTHONUTF8=1 py -3.12 example/15_lma_touschek/touschek_elegant/elegant_touschek_scan.py astrid2 1e6 0.01 0.1 1.0     # 式の差の切り分け (電流 1e6、結合の下限 MINCOUP)
PYTHONUTF8=1 py -3.12 example/15_lma_touschek/touschek_elegant/elegant_summary.py                                       # 寿命の表 (計算なし)
# 5. K1≠0 の BEND のオフモメンタムの光学 (gs・pysad・elegant・解析式)
PYTHONUTF8=1 py -3.12 example/15_lma_touschek/touschek_elegant/gs_bend_k1.py
PYTHONUTF8=1 py -3.12 example/15_lma_touschek/touschek_elegant/bend_k1_elegant.py
# 試験
PYSAD_SAD_PACKAGES=<SAD の Packages のフォルダ> py -3.12 -m pytest pysad/tests/script_tests/test_v23_T1_elegant.py -q      # elegant が無いと光学の 1 件は skip
```

## ファイル

| ファイル | 内容 |
|---|---|
| `sad2ele.py` / `astrid2.lte`・`kekls.lte`・`*_ma.lte`・`*_all.lte` (+ `.map.json`) | **SAD → elegant の変換**。pysad で展開した要素の並びから lte を書く。DRIFT → DRIF、QUAD → KQUAD、SEXT → KSEXT、OCT → KOCT、BEND → CSBEND (K1≠0 は `K2 = −K1·ANGLE/L`、`EDGE_ORDER = 2`。E1・E2 は SAD の値 × ANGLE)、CAVI → RFCA (`PHASE = 180° − asin(U0/V)`)、MARK・MONI・L = 0 で ANGLE = 0 の BEND → MARK。放射は SYNCH_RAD = 1・ISR = 0。`_ma`・`_all` は各要素の出口に MAXAMP (±10 mm、位置の名前 `AM` が momentum_aperture の対象) を足したもの。`.map.json` は要素の番号 ↔ s・名前 |
| `pysad_ref.py` / `pysad_astrid2.ref.json`・`pysad_kekls.ref.json`・`.ref.out.txt` | pysad の光学 (チューン、色収差、β・η の最大・最小、αc、U0)、自然・IBS 込みのビーム (EMITXC・EMITYC・EMITZC・SIGE・SIGZ)、一様な受容 a での Touschek 寿命、局所の LMA の寿命 (`pysad_local`) |
| `elegant_lma.py` / `elegant_astrid2.lma.txt`・`elegant_kekls.lma.txt`・`elegant_astrid2_exactk1.lma.txt`・`*_ma.ele` | elegant の `momentum_aperture` (RF・放射つき、ISR なし、MAXAMP ±10 mm、x = y = 0)。位置ごとの δ±・損失の周・位置。ASTRID2 は全周 199 位置 (Pelegant、1312 周)、kek-ls は 14 位置 (最初のセルの 11 位置 + 9 番目のセルの 3 位置。serial の elegant を位置ごとに 4 本同時、1500 周) |
| `lma_compare.py` / `lma_compare_astrid2.{table.txt,png,out.txt}`・`lma_compare_kekls.{table.txt,png,out.txt}`・`pysad_kekls.lma.txt` | LMA の比較 (同じ位置の pysad と elegant の δ±、相対差、限界の種類 = アパーチャ / RF バケット)。elegant の周期性の確認 |
| `elegant_touschek.py` / `elegant_{astrid2,kekls}.touschek.txt`・`*.aper_*.sdds` | elegant の `touschekLifetime`。入力は pysad の IBS 込みのビーム。受容: 一様な δ、elegant の LMA、pysad の LMA |
| `elegant_touschek_scan.py` / `elegant_touschek_scan_astrid2.out.txt` | 式の差の切り分け (電流 1e6 で IBS を小さくし、MINCOUP = 0.01・0.1・1 で一様な δ の寿命の elegant / pysad) |
| `elegant_summary.py` / `elegant_summary.out.txt` | 寿命の比較の表 (保存済みの結果から) |
| `gs_bend_k1.py` / `.out.txt`・`gs_bend_k1.sad` | 一様な磁場指数の扇形磁石 16 台のリング (K1≠0 の BEND) のオフモメンタムのチューン: gs・pysad の追跡 (1e-12 で一致)、解析式。B160 (K1 の暗黙の定義) の根拠 |
| `bend_k1_elegant.py` / `.out.txt` | 同じリングを elegant の CSBEND (K2 = 0、−h K1、−2 h K1) で |

## 結果の要約

- 光学: ASTRID2・kek-ls とも、チューン 1e-5 以内、αc 1.3e-4、U0 1e-5、σδ 5e-5、自然 εx は ASTRID2 で 2e-6 (kek-ls は pysad の分割 EMITDIV = 1 の既知の差 2.1%)。色収差は K1≠0 の BEND に K2 = −h K1 と EDGE_ORDER = 2 を足して合わせた (ξ の差 0.01〜0.03)。
- LMA (ASTRID2): 全周 199 位置の elegant と pysad の 1 セル 25 位置が相対差 0.5% 以内 (最大 7e-5)。アパーチャが限界の位置 (δ = 0.83%) も RF バケットで決まる位置 (δ = 1.60%) も同じ。
- LMA (kek-ls、周回数 491、アパーチャ ±10 mm): 最初のセルの 11 位置で、elegant と pysad の δ± の最大の差は 1.05e-4 (δ+)・8.95e-5 (δ−)、相対差は 0.35% 以内 (`lma_compare_kekls.out.txt`)。全位置が RF バケット (δ_rf = 2.994%) で決まる。別のセルの同じ位置との差 (elegant) は最大 7.8e-5。
- Touschek の式: 同じビーム・同じ一様な δ で elegant / pysad = 0.957〜1.025 (ASTRID2、δ = 0.4% → δ_rf)。pysad (SAD) の式は当てはめと等方の近似を含む。どちらが正確かは未確認。
- B160: K1≠0 の BEND のオフモメンタムの光学が gs = pysad と、解析式・elegant (B' 一様) で違う。「誤りではなく K1 の暗黙の定義 ((1 + hx)∂b/∂x 一様、B'' = −h B')」である (docs/29 A59)。elegant には K2 = −hK1 を足して合わせる。
