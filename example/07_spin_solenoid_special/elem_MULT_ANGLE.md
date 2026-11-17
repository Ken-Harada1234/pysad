# elem_MULT_ANGLE
種別: 要素 (ANGLE 付きの MULT) / フォルダ: 07_spin_solenoid_special / ラティス: `00_lattice/multa_ring.sad` (新規。`FODO_ring.sad` の偏向磁石を ANGLE 付きの `MULT` に置き換えた `RINGMA0` と `RINGMA1`)。3 GeV / フラグ: `NORFSW NORAD NOFLUC`。参考: `docs/28_mult_angle.md`

## 何を見る例か
`MULT` は多極磁石だが、`ANGLE`≠0 を与えると曲げ磁石の上に多極成分 ($K_0$、$K_1$、$K_2$、…) をのせた要素になる (`tmulta`)。ANGLE 付きの MULT で、(1) 多極成分なしなら BEND と同じ結果になること、(2) $K_1\ne0$ の複合機能型磁石を、解析式 (行列の積) と照合する。

## スクリプト
```
MULT MA0 = (L=1.2 ANGLE=0.2617993877991494 E1=0.1308996938995747 E2=0.1308996938995747)
     MA1 = (L=1.2 ANGLE=0.2617993877991494 E1=0.1308996938995747 E2=0.1308996938995747 K1=0.05);
LINE RINGMA0 = (ORIGIN RF 24*CMA0);   ! BEND を MA0 に置き換え
LINE RINGMA1 = (ORIGIN RF 24*CMA1);   ! K1=0.05 付き
USE RINGMA1; CELL; CALC;  r = TrackParticles[...1 周...]
```
- `K1` は積分強度 (1/m)。長さ $L=1.2$ m なので勾配は $k_1=K1/L$。
- `E1 E2` は SAD では ANGLE に対する比として読まれる (端の角度 $e=E1\cdot$ANGLE。[fn_TrackParticles](../05_tracking/fn_TrackParticles.md))。

## 結果の読み方
| リング | NX | NY | $\beta_x$(ORIGIN) | $\eta_x$(ORIGIN) |
|---|---|---|---|---|
| BEND (RING) | 6.69632894 | 4.90605604 | 5.02518491 | 0.56522773 |
| MULT ANGLE K1=0 | 6.69632894 | 4.90605604 | 5.02518491 | 0.56522773 |
| MULT ANGLE K1=0.05 | 6.98599451 | 4.52286209 | 4.86338095 | 0.52740295 |

1 周の追跡 ($x_0=y_0=10^{-5}$): BEND と K1=0 の MULT は $x=1.713426\times10^{-5}$、$p_x=1.069066\times10^{-5}$、$y=1.58646\times10^{-6}$、$p_y=5.349766\times10^{-6}$ で完全に同じ。K1=0.05 では $x=1.186292\times10^{-5}$、$p_x=1.026738\times10^{-6}$、$y=-1.17827\times10^{-5}$、$p_y=1.381372\times10^{-6}$。

## 物理・解析式
複合機能型の扇形偏向 (曲率半径 $\rho=L/\theta=4.5837$ m、勾配 $k=K1/L$) の水平・垂直の伝達行列:

$$M_x=\begin{pmatrix}\cos\phi_x&\sin\phi_x/\sqrt{K_x}\\-\sqrt{K_x}\sin\phi_x&\cos\phi_x\end{pmatrix},\ K_x=\frac1{\rho^2}+k,\qquad K_y=-k$$

($\phi=\sqrt{K}L$。$K_y<0$ は双曲線関数)。両端に薄いエッジ $\begin{pmatrix}1&0\\\pm\tan e/\rho&1\end{pmatrix}$ (水平 +、垂直 −)、$e=E1\cdot\theta=0.03427$ rad。セル全体の行列の 24 乗からチューンと β を出す (Python)。

| 量 | 式 (K1=0.05) | gs | pysad |
|---|---|---|---|
| $\nu_x$ の小数部 | 0.9859945 | 0.98599451 | 同左 |
| $\nu_y$ の小数部 | 0.5228621 | 0.52286209 | 同左 |
| $\beta_x$(ORIGIN) | 4.8633810 | 4.86338095 | 同左 |
| $\nu_x$、$\nu_y$ (K1=0 の場合) | 0.6963289、0.9060560 | 0.69632894、0.90605604 | 同左 |

## gs と pysad の比較
○ 一致 (全 6 行が印字桁まで同じ。式とも 8 桁一致)。

## 未確認・限界
- $K_1$ だけ。$K_0$ (双極の過不足)、$K_2$ 以上、スキュー成分 (`SK`)、`BZ`、`VOLT` (最初の版では pysad 未対応と書いた。v2 の SOL 領域内の MULT と領域の外の RF つき MULT の記述が `docs/47_v2_changes.md` にあるが、この例では確かめていない)、`F1`、`TAPER`、`L=0` は見ていない。
- 追跡の 1 周の式による照合はしていない (チューンと β だけ)。
- 水平のエッジ項の符号 (+) と垂直 (−) は、$K_1=0$ のチューンが一致するように選んだ (SAD ソースでの確認ではない)。
