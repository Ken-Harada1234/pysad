# flag_POL_RADPOL
種別: フラグ / フォルダ: 07_spin_solenoid_special / ラティス: `00_lattice/FODO_ring_V4MV.sad` (3 GeV、RF 4 MV、平面リング) / フラグ: `RFSW; RAD; NOFLUC; POL`。参考: `docs/43_pol_guide.md`

## 何を見る例か
`POL` はスピン追跡と偏極の量 (偏極時間、スピンチューン、偏極ベクトル) の計算を有効にするフラグ (既定 OFF)。`RADPOL` (既定 ON) は、POL が ON のとき追跡に Sokolov-Ternov 効果 (放射によるスピンの偏極) を入れる。(1) `Emittance[]` の偏極の量を式と照合、(2) `TrackParticles` のスピンの `sy` が Sokolov-Ternov の式で立ち上がること、`NORADPOL` で止まることを確かめる。

## スクリプト
```
RFSW; RAD; NOFLUC; POL;
em = Emittance[];
Print[PolarizationTime/.em, SpinTune/.em, NominalSpinTune/.em, PolarizationVector/.em];
beam = {{0},{0},{0},{0},{zeq},{0},{0},{0},{1}};  ! 9 列: {x,px,y,py,z,dp, スピン角, sy/|s|, flag}
FFS["RADPOL"];  r = TrackParticles[{1, beam}, 1, 1, 500][[2]];
FFS["NORADPOL"]; r = TrackParticles[{1, beam}, 1, 1, 500][[2]];
```
- ビームは POL のとき 9 列 (10・11 列目は損失の予約)。7 列目はスピンの角度 $\arctan(s_z,s_x)$、8 列目は $s_y/|s|$。初期値は角度 0、$s_y=0$ (水平面内のスピン)。
- `zeq=0.038315657` m は RF の平衡位置 (`EMIT` の `Equilibrium position dz`)。ここから始めると dp が小さいままになる。
- `Do` の中のフラグの切り替えは `FFS["RADPOL"]` と書く。

## 結果の読み方
| 行 | gs | pysad |
|---|---|---|
| `PolarizationTime` [s] | 143.3652012269 | 143.3652012269 |
| `SpinTune` | −0.1918545590 | −0.1918545590 |
| `NominalSpinTune` | 6.8081482248 | 6.8081482248 |
| `PolarizationVector` | {1.5e-32, 1, 1.3e-31} | {0, 1, 0} |
| 500 周、RADPOL: $s_y/\lvert s\rvert$ | $1.1349035543\times10^{-6}$ | $1.1349035543\times10^{-6}$ |
| 500 周、NORADPOL: $s_y/\lvert s\rvert$ | $9.4\times10^{-21}$ (0) | $9.4\times10^{-21}$ |
| 500 周のスピン角 (両者とも) | 0.8793923838 | 0.8793923838 |

## 物理・解析式
スピンチューン (Thomas-BMT、平面リング): $\nu_s=a\gamma$、$a=0.00115965218$ (電子の異常磁気モーメント)。$\gamma=\sqrt{1+(p/m)^2}=5870.8536$ ($p=3$ GeV/c)。1 周のスピンの歳差角は $2\pi\nu_s$ で、`SpinTune` は小数部 $-(1-0.80815)=-0.19185$、`NominalSpinTune` は $a\gamma$ そのもの。
偏極時間 (Sokolov-Ternov): 偏向磁石の曲率半径 $\rho=L/\theta=1.2/0.2618=4.5837$ m、偏向の割合 $24\times1.2/105.6=0.2727$ なので $\langle\rho^{-3}\rangle=0.2727/\rho^3$:

$$\frac{1}{\tau_{ST}}=\frac{5\sqrt3}{8}\,r_e\,\lambdabar_C\,c\,\gamma^5\langle\rho^{-3}\rangle,\qquad P_{ST}=\frac{8}{5\sqrt3}=0.92376.$$

追跡の立ち上がりは $s_y(N)=P_{ST}\left(1-e^{-NT_0/\tau_{ST}}\right)$ ($T_0=C/c=3.5222\times10^{-7}$ s)。

| 量 | 式 | gs = pysad |
|---|---|---|
| $\nu_s=a\gamma$ | 6.80814822478 | 6.8081482248 (NominalSpinTune) |
| 歳差の小数部 | −0.19185178 | −0.1918545590 (差 $2.8\times10^{-6}$: 閉軌道・平衡の dp による) |
| $\tau_{ST}$ [s] | 143.365093 | 143.365201 (差 $7.6\times10^{-7}$) |
| $s_y$ (500 周、RADPOL) | $1.134825\times10^{-6}$ | $1.134904\times10^{-6}$ (差 $7\times10^{-5}$: 平衡位置と放射の離散化) |
| $s_y$ (NORADPOL) | 0 | $9.4\times10^{-21}$ |

偏極ベクトルは平面リングなので鉛直 {0, 1, 0}。`Polarization` (平衡偏極度) は K64 が NaN を返す (docs/43 B55) ので出していない。

## gs と pysad の比較
**v2 の注記 (2026-11-17、v2.0.0 の再実行、要確認)**: `.pysad.txt` は既定の v2 モードの出力に作り直した。スピンの追跡 500 周の `spin angle` が、gs と k64 モードの .8793923838 に対して v2 では 1.1515563959 になる (RADPOL・NORADPOL とも。`sy/|s|` は 1e-6 で同程度、`PolarizationTime`・`SpinTune` など `Emittance[]` の値は同じ)。原因は未特定: v2 は追跡の端の角の行路 (B75) や BEND の端の縦磁場の積分 (B58) を直したので、一周の軌道・スピンの歳差が変わった可能性があるが、この例では切り分けていない。`PYSAD_COMPAT=k64` では保存してあった旧い出力と同じで gs と一致。以下の「○ 一致」は k64 モードの結果。

○ 一致。上の表のすべての行が印字桁まで同じ。ただし `PolarizationVector` の成分は gs が $10^{-31}$ 程度の数値雑音、pysad が厳密な 0。

## 未確認・限界
- 平面リングだけ (偏極ベクトルは鉛直)。スピン共鳴、傾いた偏極ベクトル、ソレノイドのあるリング、`Polarization` (平衡偏極度) は見ていない。
- 500 周後のスピン角 (0.8794 rad) は、残る $dp$ の振動とスピンチューンの積算が入るので解析的に説明していない。
- 追跡は RAD NOFLUC (ゆらぎなし) の 1 粒子。FLUC のとき (乱数あり) の確認は docs/43 に gs 一致の記録がある。
- 偏極時間の式は、勾配 (四極) による補正の項を入れていない。この平面リングでは式と 7.6e-7 で一致した。
