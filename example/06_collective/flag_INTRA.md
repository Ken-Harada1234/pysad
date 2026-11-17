# flag_INTRA
種別: フラグ / フォルダ: 06_collective / ラティス: `00_lattice/FODO_ring_V4MV_0p5GeV.sad` (新規。`FODO_ring_V4MV.sad` の運動量を 0.5 GeV に下げた版) / フラグ: `RFSW; RAD; FLUC; MINCOUP=0.01`

## 何を見る例か
`INTRA` を ON にすると、エミッタンス計算 (`Emittance[]`、`EMIT`) にバンチ内散乱 (IBS: 同じバンチの粒子どうしのクーロン散乱) が入る。バンチの粒子数 `PBUNCH` を変えて、平衡エミッタンス・エネルギー広がり・バンチ長がどう増えるかを見る。3 GeV では IBS が小さすぎて見えないので 0.5 GeV にした。

## スクリプト
```
NOINTRA;  em0 = Emittance[];            ! IBS なし
INTRA;
Do[ PBUNCH = {1e4, 1e8, 3e8, 1e9, 3e9}[[k]];
    em = Emittance[];
    Print[Emittances/.em, MomentumSpread/.em, BunchLength/.em], {k, 1, 5}];
```
- `MINCOUP`: IBS の計算で使う垂直エミッタンスの下限 = `MINCOUP*(EMITX+EMITY)` (理想リングは $\varepsilon_y=0$ で、IBS が非現実的に大きくなるのを避ける)。
- 運動量はラティス側 (MAIN レベル) で 0.5 GeV にした。`FFS` の中で `MOMENTUM=0.5e9` と書く方法は、gs は反映するが、当時の pysad は反映しなかったので避けた (現在の pysad は反映する。[flag_SPAC](../05_tracking/flag_SPAC.md))。

## 結果の読み方
| PBUNCH | $\varepsilon_x$ [m] | $\varepsilon_y$ [m] | $\varepsilon_z$ [m] | $\sigma_\delta$ | $\sigma_z$ [mm] |
|---|---|---|---|---|---|
| NOINTRA | 5.093949e-9 | 0 | 2.066340e-7 | 2.001562e-4 | 1.056 |
| $10^4$ | 5.100944e-9 | 1.0e-16 | 2.068383e-7 | 2.002551e-4 | 1.056 |
| $10^8$ | 1.698986e-8 | 1.64e-13 | 5.538140e-7 | 3.276815e-4 | 1.728 |
| $3\times10^8$ | 2.386539e-8 | 2.58e-13 | 7.544618e-7 | 3.824625e-4 | 2.017 |
| $10^9$ | 3.507525e-8 | 4.12e-13 | 1.081594e-6 | 4.579338e-4 | 2.415 |
| $3\times10^9$ | 5.019766e-8 | 6.19e-13 | 1.522903e-6 | 5.433839e-4 | 2.866 |

(gs の値。pysad も同じ。NOINTRA の $\varepsilon_y$ は gs で $5\times10^{-36}$、pysad で 0 (数値雑音の違い)。)

## 物理・解析式
式による照合は、次の極限と傾向だけ (IBS の Bjorken-Mtingwa 式の係数は検証していない)。

| 性質 | 式 | 結果 |
|---|---|---|
| $N\to0$ で IBS なしに戻る | $\varepsilon_x(N\to0)=\varepsilon_{x0}$ | $N=10^4$: $5.100944\times10^{-9}$ vs NOINTRA $5.093949\times10^{-9}$ (差 0.14%) |
| 増加は $N$ の増加より緩やか (放射減衰との釣り合い: 平衡では IBS 拡散 ∝ 減衰 × ε) | $\varepsilon_x\propto N^{p}$, $p<1$ | $N=10^8\to10^9$: $p=0.31$、$3\times10^8\to10^9$: $p=0.32$、$10^9\to3\times10^9$: $p=0.33$ |
| $\sigma_z=\sigma_\delta\cdot|\eta_c|c/\omega_s$ (縦の関係) | $\sigma_z/\sigma_\delta$ 一定 | $N=10^8$: 5.27 m、$3\times10^9$: 5.27 m (一定) |

縦の関係は、バンチ長がエネルギー広がりに比例して伸びる (RF が一定) ことの確認 ($\sigma_z/\sigma_\delta=5.27$ m は 6 行で変わらない)。

## gs と pysad の比較
○ 一致 (全 6 行が印字桁まで同じ)。

## 未確認・限界
- 0.5 GeV (電子)、RF 4 MV、理想リング (垂直エミッタンスが MINCOUP の下限)、粒子数 $10^4$〜$3\times10^9$。実機のパラメータではない。
- IBS の絶対値の式 (Piwinski、Bjorken-Mtingwa、Bane の高エネルギー近似) との照合は未実施。$N^{1/3}$ に近い傾きは計算結果から読んだもので、理論的な導出はしていない。
- pysad は `INTRA` を `SOL` のあるリングとは併用できない (docs/44)。
