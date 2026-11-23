# Touschek 寿命の式の例 (pysad (K64) と elegant のどちらが正確か。独立な式での確定)

試験は `pysad/tests/script_tests/test_v23_T2_touschek.py` (8 件)。独立な式は `pysad/touschek.py` の `touschek_piwinski` (v2 の `TouschekLifetime[INF, INF, nz]`・`TouschekLifetimeLMA[]` の既定)。
実行はリポジトリ直下から `PYTHONUTF8=1 py -3.12 example/15_lma_touschek/touschek_elegant/touschek_formula_xxx.py ...`。v2 モード。elegant は WSL (Ubuntu、`~/elegant/env.sh`、elegant 2026.4.0。作業フォルダ `~/t2/<名前>`)。
ラティス: `example/astrid/lattice/ASTRID2_v01_01.sad` (ASTRID2、0.58 GeV)、`pysad/tests/kek-ls_v03_12.lat` (KEK-LS、大きなリング A。ESRF 型の HMBA、3 GeV)。`sad2ele.py` の lte (`<名前>_all.lte`、`elegant_comparison.md`) を使う。
CPU は pysad 1 本 + elegant 1 本 (上限 4 コア)。

## 結論

- **elegant の式 (Piwinski) が正確、K64 の式 (pysad の `TouschekTable`) が不正確**。同じビーム・同じ光学・同じ一様な δ で、独立な式 (2 粒子の Møller 散乱の積分、`touschek_piwinski`) との寿命の差:
  K64 は結合 1% で +0.1〜+6%、10% で +1〜+14%、100% (丸いビーム) で −14〜+21%。elegant は標本化を細かくすると −1.25〜+0.11% (36 点)。
- 独立な式は、3 次元の一般の積分 (1e-3 以内)・厳密な相対論の運動学の対のモンテカルロ (0.3% 以内)・平坦ビームの極限で Bruck・Piwinski の関数 D(ε) (記憶に基づく式、6 桁) と一致。
- K64 の差の原因: 当てはめの関数 (±1%) ではなく、x-y 面の運動量の広がりを平坦ビームの 1 次元の関数に ptrans² = tr(P) で入れる近似。elegant の誤差: 要素の出口だけの標本化 (粗いリングで寿命 +0.7〜+4.8%)。
- v2 の既定は Piwinski 型 (`Formula->"K64"` と k64 モードは K64 の式)。B162。

## 走らせるコマンド

```
# 1. 3 つの方法の寿命の表 (K64 / elegant / 独立な式)。結合 (MINCOUP)・電流 (PBUNCH)・受容を変える。elegant は光学を全要素 8 分割にして標本化を細かくする
PYTHONUTF8=1 py -3.12 example/15_lma_touschek/touschek_elegant/touschek_formula_scan.py astrid2 --sets "1e9:0.01,1e9:0.1,1e9:1,1e6:0.01,1e6:1,3e10:0.01"      # 約 5 分
PYTHONUTF8=1 py -3.12 example/15_lma_touschek/touschek_elegant/touschek_formula_scan.py kekls --sets "1e9:0.01,1e9:0.1,1e9:1"
# 2. 独立な式の検証: 3 次元の一般の積分と、実験室系の対のモンテカルロ (厳密な相対論) との比 (スライスごと)
PYTHONUTF8=1 py -3.12 example/15_lma_touschek/touschek_elegant/touschek_formula_validate.py astrid2        # 結合 1%・100% (1e9:0.01 1e9:1)
PYTHONUTF8=1 py -3.12 example/15_lma_touschek/touschek_elegant/touschek_formula_validate.py kekls
# 3. K64 の式の差の原因の切り分け (近似を 1 つずつ置き換える)
PYTHONUTF8=1 py -3.12 example/15_lma_touschek/touschek_elegant/touschek_formula_cause.py astrid2 1e9:0.01 1e9:0.1 1e9:1
PYTHONUTF8=1 py -3.12 example/15_lma_touschek/touschek_elegant/touschek_formula_cause.py kekls 1e9:0.01 1e9:1
# 4. 解析式との比較 (極限)、K64 の当てはめの誤差、バンチ内の δ の分布の見積もり
PYTHONUTF8=1 py -3.12 example/15_lma_touschek/touschek_elegant/touschek_formula_limits.py astrid2 kekls
# 5. 局所の LMA の寿命: K64 の表と Piwinski 型
PYTHONUTF8=1 py -3.12 example/15_lma_touschek/touschek_elegant/touschek_formula_lma_lifetime.py astrid2
PYTHONUTF8=1 py -3.12 example/15_lma_touschek/touschek_elegant/touschek_formula_lma_lifetime.py kekls
# 試験
PYSAD_SAD_PACKAGES=<SAD の Packages のフォルダ> py -3.12 -m pytest pysad/tests/script_tests/test_v23_T2_touschek.py -q
```

SAD の台本では (INTRA ON の `Emittance[]` のあと):
```
Print[TouschekLifetime[INF, INF, 0.004/SIGE]];                       ! v2: Piwinski 型。k64 モード: K64 の表
Print[TouschekLifetime[INF, INF, 0.004/SIGE, Formula->"K64"]];       ! 表
Print[TouschekLifetimeLMA[Formula->"Piwinski"]];                      ! MomentumAcceptance[] のあと
```
Python:
```
from pysad.touschek import touschek_piwinski
r = touschek_piwinski(ffs.touschek, a=0.006)               # 一様な受容
r = touschek_piwinski(ffs.touschek, ap=dpp, am=abs(dpm))    # 要素 l ごと (長さ nlat) の δ+・|δ−|
```

## ファイル

| ファイル | 内容 |
|---|---|
| `touschek_formula_scan.py` / `touschek_formula_scan_{astrid2,kekls}.out.txt`・`.json`・`touschek_formula_{astrid2,kekls}_div.ele` | 寿命の表 (pysad K64、elegant、独立な式)。elegant の入力は pysad の IBS 込みのビーム。`touschek_formula_*_div.ele` は elegant の Twiss (全要素 8 分割) |
| `touschek_formula_validate.py` / `touschek_formula_validate_{astrid2,kekls}.out.txt` | 1 次元の式と 3 次元の積分・厳密な運動学のモンテカルロの、スライスごとの比 (3D/1D = 0.9993〜1.0003、MC/1D = 0.9969〜1.0003) |
| `touschek_formula_check3d.py` | 3 次元の一般の積分 (q の方向は準モンテカルロ、動径は log の Gauss-Legendre、散乱角の半球の積分は θ の数値積分と φ' の解析)。q_s・非対角を入れる |
| `touschek_formula_mc.py` | 実験室系の 6 次元ガウス分布からの対のモンテカルロ。4 元運動量・重心系への Lorentz 変換・Møller の角分布・実験室系へ戻して δ_new = |P|/p0 − 1 を厳密に求める |
| `touschek_formula_cause.py` / `touschek_formula_cause_{astrid2,kekls}.out.txt` | K64 の差の原因: E0 (独立) / E1 (平坦の厳密な関数、λ1) / E2 (平坦の関数、tr P) / E3 (K64 の再現) / E4 (K64 の当てはめ、λ1) |
| `touschek_formula_limits.py` / `touschek_formula_limits.out.txt` | 平坦ビームの極限 = D(ε) (記憶に基づく式)、K64 の当てはめの誤差 (±1%)、a − δ_old の効き (概算) |
| `touschek_formula_lma_lifetime.py` / `touschek_formula_lma_lifetime.out.txt` | 局所の LMA の寿命 (K64 の表 / Piwinski 型): ASTRID2 9501.7 s → 9327.1 s (−1.8%)、KEK-LS 26875.7 s → 26862.5 s (−0.05%) |

## 結果の要約

ASTRID2 (IBS 込み、PBUNCH 1e9、結合 1%、一様な δ = 0.4%・0.6%・1.0%・1.6%) の寿命 [s]: pysad (K64) 453.5・1270・4970・19218、elegant 425.5・1195・4765・18753、独立な式 427.1・1200・4794・18912。
K64 / 独立 − 1 = +6.2・+5.8・+3.7・+1.6%、elegant / 独立 − 1 = −0.36・−0.42・−0.59・−0.85%。KEK-LS (一様な δ = 0.4%・0.6%・1.0%・2.0%) は K64 / 独立 − 1 = +3.8・+2.4・+1.0・+0.1%、elegant / 独立 − 1 = −0.37・−0.35・−0.29・−0.07%。
結合 100% (丸いビーム): K64 は ASTRID2 で −13.8 (0.4%) 〜 +15.8% (1.6%)、KEK-LS で −0.7 〜 +20.9% (2%)。elegant の粗い標本化 (要素の出口のみ) は ASTRID2 で +2.0〜+4.8%、KEK-LS で +0.7%。

使った条件・限界: kick の基準 (粒子のもとの δ = 0。K64・elegant・独立な式とも共通。a/σδ ≲ 10 で損失率が 5〜20% 増えうる概算が `touschek_formula_limits.py`)、非相対論的な Møller (0.3%)、q_s = 0 (3 次元で 1e-3)、誤差なし・結合の小さいリング。
