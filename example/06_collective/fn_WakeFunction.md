# fn_WakeFunction
種別: 関数 / フォルダ: 06_collective / ラティス: `00_lattice/wake_ring.sad` (新規。`FODO_ring.sad` の ORIGIN の直後に MARK `WK` を置いた `RINGW`)。3 GeV / フラグ: `TRPT; INS; LWAKE; TWAKE; NORFSW NORAD NOFLUC`、`PBUNCH=1e10 NBUNCH=1`

## 何を見る例か
`WakeFunction[Longitudinal, 要素名] = {{z1,w1},...}` (縦、z は m、w は V/C) と `WakeFunction[Transverse, 要素名]` (横、w は V/C/m) で要素ごとにウェイク関数の表を与える。ウェイクは `TRPT` の追跡 (`TrackParticles`) で効き、ここでは一定値の表を与えて、最大のキックを解析式と照合する。表の使い方 (`dz`、最後の行が使われない癖) も確かめる。

## スクリプト
```
USE RINGW; TRPT; INS;
PBUNCH = 1e10; NBUNCH = 1;
W0 = 1e9; Wt = 1e12;
LWAKE; TWAKE;
WakeFunction[Longitudinal, "WK"] = {{-1, W0}, {0, W0}, {1e-9, 0}, {1, 0}, {2, 0}};
WakeFunction[Transverse,   "WK"] = {{-1, Wt}, {0, Wt}, {1e-9, 0}, {1, 0}, {2, 0}};
r = TrackParticles[{1, beam}, 3, 1][[2]];      ! WK (位置 2) の直後まで
```
- 表は `dz = z(源の粒子) − z(受ける粒子)` の関数。ここでは $dz\le0$ で一定値、$dz>0$ で 0。**追跡では表の最後の 1 行が使われない** (K64 の癖、docs/13 §5.3) ので、0 の行を 2 本置いた。
- ビームは 100 粒子: $x_0=10^{-4}+10^{-5}\sin(0.37i)$ (平均 $1.000305\times10^{-4}$ m)、$z_0=0.01\sin(0.13i)$、他は 0。ウェイクは要素の前後で半分ずつ。

## 結果の読み方
| 行 | gs = pysad |
|---|---|
| `WakeFunction L:`、`T:` | 定義した表をそのまま返す |
| z 最大の粒子 | $dp=-5.340589\times10^{-10}$、$p_x=5.3422176\times10^{-11}$ |
| z 最小の粒子 | $dp=-2.864697\times10^{-11}$、$p_x=2.9540413\times10^{-12}$ |
| 平均 | $dp=-2.927789\times10^{-10}$、$p_x=2.9357286\times10^{-11}$ |

## 物理・解析式
全粒子から受けるウェイクの合計 (表が一定 $W_0$、全 $N=$ PBUNCH 粒子):

$$\delta_{max}=-\frac{e\,N\,W_0}{p_0c},\qquad p_{x,max}=\frac{e\,N\,W_t\,\langle x\rangle}{p_0c}.$$

(1 粒子あたりの電荷 $e$、バンチ電荷 $Q=eN$、ウェイク電位 $V=W_0Q$、$\Delta E=eV$ を $p_0c$ で割ったもの。符号は SAD の `CHARGE`(電子)。)

| 量 | 式の値 | gs | pysad |
|---|---|---|---|
| $\delta_{max}$ ($N=10^{10}$、$W_0=10^9$ V/C、$p_0c=3$ GeV) | $-5.340589\times10^{-10}$ | $-5.340589\times10^{-10}$ | 同左 |
| $p_{x,max}$ ($W_t=10^{12}$、$\langle x\rangle=1.000305\times10^{-4}$) | $5.342218\times10^{-11}$ | $5.3422176\times10^{-11}$ | 同左 |

z が最大の粒子はすべての粒子から受けて最大のキック、z が最小の粒子は自分自身 (同じスライス) の寄与だけを受ける ($dp=-2.86\times10^{-11}$ は最大の 5.4%)。平均は約 55% (粒子を z の順に並べたときに各粒子が受ける寄与の個数の平均と、自己項による)。

## gs と pysad の比較
○ 一致 (全 8 行、印字桁まで同じ)。

## 未確認・限界
- 一定値の表 (ステップ関数) と 100 粒子、WK 1 か所、$N=10^{10}$ の 1 バンチだけ。実際の因果的なウェイク (共鳴器、抵抗壁) は見ていない。
- どちらの端がバンチの先頭かという SAD の符号の取り方は、観測 (z が大きい粒子が全粒子分を受ける) から書いただけで、物理的な意味は確かめていない。
- 光学 (`CALC`) 側のウェイクは [flag_LWAKE_TWAKE_WAKEOPT](flag_LWAKE_TWAKE_WAKEOPT.md)。
