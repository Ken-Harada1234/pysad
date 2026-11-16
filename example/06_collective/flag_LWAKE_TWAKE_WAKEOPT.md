# flag_LWAKE_TWAKE_WAKEOPT
フラグ・06_collective・ラティス `00_lattice/wake_ring.sad` (RINGW)、3 GeV、`TRPT; INS; NORFSW NORAD NOFLUC`

## 何を見る例か
`LWAKE` (縦のウェイク)、`TWAKE` (横のウェイク) は追跡のウェイクを有効にするフラグ (`NO` を付けると OFF)。`WAKEOPT` は、さらに光学計算 (`CALC`) でもウェイクを入れる (`InitialOrbits` が必要)。4 通りの組み合わせの追跡と、`WAKEOPT` の光学を確かめる。

## スクリプト
```
FFS[{"NOLWAKE NOTWAKE", "LWAKE NOTWAKE", "NOLWAKE TWAKE", "LWAKE TWAKE"}[[k]]];   ! Do の中は FFS[...]
r = TrackParticles[{1, beam}, 3, 1][[2]];        ! Max[Abs[dp]], Max[Abs[px]]
LWAKE; TWAKE; WAKEOPT;
InitialOrbits = {{1e-3, 0, 0, 0, -0.01, 0}, {-1e-3, 0, 0, 0, 0.01, 0}};
FITP 2; CALC;
```
ウェイク関数の表は [fn_WakeFunction](fn_WakeFunction.md) と同じ ($W_0=10^9$、$W_t=10^{12}$、$N=10^{10}$)。`InitialOrbits` の各組は中心軌道からのずれ $\{x,p_x,y,p_y,z,\delta\}$。

## 結果の読み方
| case | LWAKE | TWAKE | max$|dp|$ | max$|p_x|$ |
|---|---|---|---|---|
| 1 | OFF | OFF | 0 | 0 |
| 2 | ON | OFF | $5.3405888\times10^{-10}$ | 0 |
| 3 | OFF | ON | 0 | $5.3422176\times10^{-11}$ |
| 4 | ON | ON | $5.3405888\times10^{-10}$ | $5.3422176\times10^{-11}$ |

WAKEOPT の光学 (2 軌道 + 中心軌道の 3 本): 要素 1 `ORIGIN` から 4 `SF.1` の `DX, DPX, DDP`:

| 要素 | DPX | DDP |
|---|---|---|
| ORIGIN (1) | 0 | 0 |
| WK (2、前半のキック後) | $8.9009813\times10^{-11}$ | $-4.450491\times10^{-11}$ |
| RF (3、WK の後) | $1.7801963\times10^{-10}$ | $-8.900981\times10^{-11}$ |
| SF.1 (4) | $1.7801963\times10^{-10}$ | $-8.900981\times10^{-11}$ |

## 物理・解析式
追跡のキックの最大値は [fn_WakeFunction](fn_WakeFunction.md) の式 ($\delta_{max}=-eNW_0/p_0c$、$p_{x,max}=eNW_t\langle x\rangle/p_0c$) どおり。
光学 (WAKEOPT) の中心軌道の横キックは、中心軌道 ($x=0$) が、ウェイクの源側になる軌道 1 ($x_1=10^{-3}$ m、$z=-0.01$) から受けるもの: 軌道の本数 $n_o=3$ (中心 + 2) で重みが $1/n_o$:

$$\Delta p_x=\frac{e\,N\,W_t\,x_1}{p_0c\,n_o}=\frac{1.602\times10^{-19}\times10^{10}\times10^{12}\times10^{-3}}{3\times10^9\times3}=1.7802\times10^{-10}.$$

これは gs・pysad の WK 通過後の DPX ($1.7801963\times10^{-10}$、前後で半分ずつ: 前半が $8.90\times10^{-11}$) と一致する (式 $1.78019626\times10^{-10}$)。DDP (縦) は式の期待値 $-1.78\times10^{-10}$ の半分になる (縦のキックは軌道の数が奇数のとき、軌道の対応が 1 つずれる K64 の癖 (docs/36 §6) が出ており、そのまま gs と同じ)。

## gs と pysad の比較
○ 一致 (追跡の 4 行、光学の 4 行とも印字桁まで同じ)。

## 未確認・限界
- フラグの切り替えは、`Do` の中では `FFS["LWAKE"]` のように書く必要がある (`LWAKE;` と書いても `Do` の中では効かない)。
- 光学側は 3 本の軌道、WK 1 か所のみ。加速のある TRPT では pysad は光学のウェイクを持たない (docs/36 §4)。この例のリングは `NORFSW` (RF なし) で行った。
- 縦のキック (DDP) は式と合わない (癖のため)。
