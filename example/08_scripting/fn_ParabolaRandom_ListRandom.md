# fn_ParabolaRandom_ListRandom (ParabolaRandom, ListRandom, SeedRandom[plugin])
種別: 関数 (乱数) / フォルダ: 08_scripting / FFS 内の SADScript (ラティスなし)。gs は SAD 1.4.4.2k64、pysad は `run_sad_string`。

## 何を見る例か
`ListRandom[]` は使える乱数のプラグイン名のリスト (gs・pysad とも {"SAD"})、`SeedRandom["SAD"]` でプラグインを選ぶ。`ParabolaRandom` のヒストグラム。

## スクリプト
```
ListRandom[];   Scan[(SeedRandom[#]; SeedRandom[3]; Print[#, ": ", Random[3], " ", GaussRandom[2], " ", ParabolaRandom[2]]) &, ListRandom[]];
SeedRandom[3]; x = ParabolaRandom[40000]; h = Table[Count[Floor[(x + 1)*5], k], {k, 0, 9}];
```

## 結果の読み方
`SAD: {.534...,.0369...,.3325...} {.0999...,-1.3546} {.1569..., ...}`。ヒストグラム (10 区間、40000 個): {1194,3080,4517,5465,5963,5982,5336,4388,2999,1076}。放物線 $\frac34(1-x^2)$ を 10 区間に積分した比 (0.0280, 0.0768, 0.1128, ...) × 40000 = 1120, 3072, 4512, 5472, 5976, 5976, 5472, 4512, 3072, 1120 に近い。

## 物理・解析式
区間 $[x_a,x_b]$ の確率 $=\frac34\left[x-x^3/3\right]_{x_a}^{x_b}$。例: 中央 [-0.2, 0] は $\frac34(0.2-0.00267)=0.14800$ → 5920 (実測 5963/5982)。

## gs と pysad の比較
○ ヒストグラムは **完全に一致**。`SeedRandom["SAD"]` (プラグイン名の文字列) も gs と同じ (2026-11-17 確認。乱数の種類の選択だけで状態は変わらない)。`ListRandom[]` も一致。

## 未確認・限界
- 1 つしかないプラグイン ("SAD") しか見ていない。
