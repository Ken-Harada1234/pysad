# fn_AccelerateParticles
関数・05_tracking・ラティス `00_lattice/FODO_ring_V4MV.sad` (RF 4 MV、3 GeV)・`RFSW; RAD; NOFLUC; Emittance[]; NORAD`・120 周

## 何を見る例か
`AccelerateParticles[beam, mom, {n, nturn}]` はリングで運動量を周回ごとに変えながら追跡する (断熱減衰が自動で入る)。`mom` は周回番号 `n` の式で、ここでは `p0 + 1e6 n` (1 周あたり 1 MeV/c、120 周で 3.12 GeV/c)。

## スクリプト
```
em = Emittance[];          ! RF の同期位相 (RAD ON の平衡) を決める
NORAD;
beam = {1, {{0},{0},{y0},{0},{0},{0},{1}}};        ! {位置, 座標} の形で渡す
r = AccelerateParticles[beam, p0 + dE*n, {n, 120},
      Synchronize->((zs = Append[zs, #2[[2,5,1]]]; dps = Append[dps, #2[[2,6,1]]])&)];
```
- 内部 (`Packages/dapert.n`): 各周で `MOMENTUM=p1; u=p0/p1`、`px, py ← ×u`、`z ← z/v` ($v=u\sqrt{(p_1^2+m^2)/(p_0^2+m^2)}$)、$dp\leftarrow u\,dp-(1-u)$ として 1 周追跡する。
- `Synchronize->(...)&` は毎周呼ばれる (`#1`=周番号、`#2`=`{位置, ビーム}`)。式は `->((...)&)` と括弧で包む (括弧なしだと gs は Rule でないとエラー)。
- 戻りは `{1, ビーム}`。

## 結果の読み方
- `p0=3.0e9 p1=3.12e9 p0/p1=0.96154`: 最終の運動量と比。
- `inv0`、`inv1`: 垂直 ($y_0=1$ mm、分散なし) の Courant-Snyder 不変量の初期値と最終値 (ORIGIN の $\alpha_y,\beta_y$ で計算)。比 `inv1/inv0`。
- `z mean`、`zmax`、`zmin`、`dp max/min`: 毎周の $z$、$dp$ (Synchronize で記録) の統計。

## 物理・解析式
(1) 断熱減衰: $J_y\propto1/p$ なので $J_y(p_1)/J_y(p_0)=p_0/p_1$。(2) 加速ぶんの RF の平衡位相: 1 周の RF 電圧 $V\sin(\phi_s+kz)$ が放射損失 $U_0$ と加速 $\Delta E$ を補うので、$\sin(\phi_s+kz_{eq})-\sin\phi_s=\Delta E/V$。ここで $\sin\phi_s=U_0/V=1.5633/4=0.3908$ ($\phi_s=0.4015$ rad)、$k=2\pi f/c=10.47$ m$^{-1}$。

| 量 | 式 | gs | pysad |
|---|---|---|---|
| $J_{y,1}/J_{y,0}$ | $p_0/p_1=0.96154$ | 0.96076 | 0.96076 |
| $z_{eq}$ [m] | 0.02806 ($\sin\phi=0.6408$ から) | 中点 $(z_{max}+z_{min})/2=0.02932$、平均 0.03040 | 同じ |

粒子は $z=0$ から始めるので、平衡位置 $z_{eq}$ のまわりに振幅 $z_{eq}$ で振動し、$z_{max}\approx2z_{eq}=0.0561$ (観測 0.0586)。式との差 4〜8% は、$\sin$ の非線形 (振幅が大きい) による。

## gs と pysad の比較
○ 一致。全 4 行が同じ (印字桁まで)。

## 未確認・限界
- 3 GeV、1 MeV/周 × 120 周 (調整なし)。放射は OFF なので、実機のように $U_0$ を差し引く運動は見ていない。
- 水平方向は分散の項 (ηδ) が大きく (dp が ±2e-3)、不変量は垂直で評価した。
- `mom` を運動量の式 (eV/c) として渡す。エネルギーとの区別は MASS (電子) による。
