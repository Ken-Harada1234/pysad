# fn_Random_GaussRandom_SeedRandom
種別: 関数 (乱数) / フォルダ: 08_scripting / FFS 内の SADScript (ラティスなし)。gs は SAD 1.4.4.2k64、pysad は `run_sad_string`。

## 何を見る例か
`Random` (一様)、`GaussRandom` (正規、`GCUT` で打ち切り)、`ParabolaRandom` (放物線分布)、`SeedRandom[seed]` で列を再現、`SeedRandom[]` で状態を取り出して復元。

## スクリプト
```
SeedRandom[12345]; r1 = Random[5]; SeedRandom[12345]; r2 = Random[5];   (* r1 === r2 *)
SeedRandom[1]; u = Random[20000]; g = GaussRandom[20000]; p = ParabolaRandom[20000];
st = SeedRandom[]; a = Random[3]; SeedRandom[st]; b = Random[3];
```

## 結果の読み方
種を同じにすると同じ列 (`r1===r2`=1)。20000 個の一様乱数の平均 0.503181122013554、分散 0.083266 ($1/12=0.083333$)。`GaussRandom` の分散 0.99394 (標準偏差 1)、`GCUT` の既定は 1e+300 (打ち切りなし)、最大絶対値 4.03。`ParabolaRandom` の分散 0.19796 (理論値 $\int_{-1}^1x^2\cdot\frac34(1-x^2)dx=1/5$)、範囲 [-0.9957, 0.9860]。

## 物理・解析式
$[0,1)$ の一様分布の分散は $1/12=0.0833$ (表の 0.08327)、$\frac34(1-x^2)$ の放物線分布の分散は $1/5$ (表の 0.19796)、正規分布は 1 (0.9939)。標準誤差 $\sqrt{2/N}\approx0.01$ の範囲で一致。

## gs と pysad の比較
○ **乱数の列は gs と pysad で一致** (平均 0.503181122013554、`GaussRandom` の最大値 4.03258726603887 まで同じ)。2026-11-17 に 2 点を gs に合わせた: (1) `Random[2,3]` は 2×3 のテンソル ({2,3})。(2) `SeedRandom[状態]` は gs でも復元されない (gs の C 実装はリストの数を読まず、先頭が文字列なら乱数の種類の選択だけで状態は変わらない) ので、pysad も同じ (`a===b` は 0)。`SeedRandom["SAD"]` は状態を変えず現在の状態を返す。

## 未確認・限界
- 乱数の列が一致するのは、同じ種・同じ呼び出し順のときのみ。
- `SeedRandom["プラグイン名"]` は `fn_ParabolaRandom_ListRandom` を参照。
