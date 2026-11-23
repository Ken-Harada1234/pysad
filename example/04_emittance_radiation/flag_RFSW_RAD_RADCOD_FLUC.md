# flag RFSW, RAD, RADCOD, FLUC
種別: フラグ / フォルダ: 04_emittance_radiation / ラティス: `00_lattice/FODO_ring_V4MV.sad` (FODO_ring と同じ磁石列、RF 電圧だけ 4 MV。QF から始まる `RINGQ` を使う)。3 GeV、周長 105.6 m、24 セル、RF 500 MHz / フラグ: `RFSW`・`NORFSW`、`RAD`・`NORAD`、`RADCOD`・`NORADCOD`、`FLUC`・`NOFLUC`。

## 何を見る例か
4 つのフラグが `Emittance[]` (平衡の計算) と `TrackParticles` (追跡) のどちらにどう効くかを表にする。

## ラティスについて
- `FODO_ring.sad` の RF は 1.5 MV で、エネルギー損失 $U_0=1.563$ MV より小さく、バケットがない (同期位置 dz=150 mm、$\nu_s$ 虚数)。そこで RF 電圧を 4 MV にした `FODO_ring_V4MV.sad` を追加した。
- RADINT は先頭の要素が SEXT (L=0.2, K2=0) だと `DL=0. IN TCAL` で I2・I3・I5 が NaN になる (gs も pysad も同じ)。そのため先頭が QF の `RINGQ = (ORIGIN RF 24*CELLQ)` を使う (磁石の並びは同じ、位相だけ違う)。
- SAD の BEND の `E1`, `E2` は ANGLE に対する比なので、端の角は $\epsilon=E1\cdot\text{ANGLE}=0.0343$ rad (θ/2 ではない)。四極の `K1` は積分強度 (k=K1/L)。

## スクリプト
フラグの組ごとに `Emittance[]` を呼び、$U_0$、$V_c$、同期位置 dz、$\epsilon_x$、$\sigma_\delta$、$\sigma_z$ を `Print` する。最後に FLUC あり・なしで、原点から 200 粒子を 100 周追跡して $\delta$ の rms を比べる (`SeedRandom[7]`)。

## 結果の読み方 (gs。pysad も同じ)
| RFSW | RAD | RADCOD | FLUC | $U_0$ [eV] | dz [m] | $\epsilon_x$ | $\sigma_\delta$ | $\sigma_z$ [m] |
|---|---|---|---|---|---|---|---|---|
| 0 | 0 | 0 | 1 | 1563265.7 | 0.038316 | 1.7355E-7 | 0.001190 | 0.016315 |
| 0 | 1 | 0 | 1 | 1563265.7 | 0.038316 | 1.7355E-7 | 0.001190 | 0.016315 |
| 1 | 0 | 0 | 1 | 1563265.7 | 0.038316 | 1.7358E-7 | 0.001191 | 0.016060 |
| 1 | 1 | 0 | 1 | 1563265.7 | 0.038316 | 1.7358E-7 | 0.001191 | 0.016060 |
| 1 | 1 | 1 | 1 | 1563263.2 | 0.038316 | 1.7358E-7 | 0.001191 | 0.016060 |
| 0 | 1 | 1 | 1 | 1562441.8 | 3.0206E-4 | 1.7201E-7 | 0.001192 | 0.016334 |
| 1 | 1 | 0 | 0 | 1563265.7 | 0.038316 | 1.7358E-7 | 0.001191 | 0.016060 |

- `Emittance[]` は **RAD と FLUC に依存しない** (RAD を切っても、NOFLUC でも平衡値は同じ。これらは追跡 `TrackParticles` のフラグ)。効くのは RFSW と RADCOD。
- RFSW OFF でも $V_c$ と $\sigma_z$ は 0 でなく、空洞の積算量から作られる (docs/09 の tecalc)。RFSW ON との差は $\sigma_z$ で 1.6%、$\epsilon_x$ で 0.02%。
- RADCOD ON は閉軌道にエネルギー損失を入れる (RFSW ON で $U_0$ が 2.5 eV 小さい)。RFSW OFF と組み合わせると、同期位置が 0.0003 m になり、$U_0$ は 0.05%、$\epsilon_x$ は 0.9% 小さくなる (原因は未調査。RFSW OFF + RADCOD ON は docs/09 に別の閉軌道の扱いがある)。
- 追跡 (200 粒子、100 周、原点から): FLUC ON の $\delta$ rms は gs 4.2074E-4、pysad 4.5223E-4。NOFLUC は 2.6063E-4 (gs・pysad 同じ、原点が同期位置からずれているための同期振動)。量子励起の分は $\sqrt{4.2^2-2.6^2}\times10^{-4}=3.3\times10^{-4}$ (gs)、$3.7\times10^{-4}$ (pysad)。

## 物理・解析式
$\delta$ の励起の成長は $\sigma_\delta^2(n)=\sigma_{\delta0}^2\,[1-e^{-2nT_0/\tau_\epsilon}]$ ($\tau_\epsilon=0.6646$ ms、100 周 = 35.2 μs なので $1-e^{-0.1060}=0.1006$)。$\sigma_{\delta0}=1.1914\times10^{-3}$ なので $\sigma_\delta\simeq3.78\times10^{-4}$ で、gs・pysad の量子励起分 (3.3、3.7)E-4 と同じ桁 (200 粒子の統計誤差は 5% 程度、同期振動の位相の影響もある)。

## gs と pysad の比較
**v2 と k64 の違い**: `.pysad.txt` は既定の v2 モードの出力。v2 は `Emittance[]`・EMIT・追跡の Jx の端の項 (B75) と `MomentumCompaction` (B77) を K64 のバグに依らず正しい式で計算するので、この例の pysad の値のうち εx (約 2.4% 大きい)・εz・減衰率の X と Z・σδ・バンチ長・α・バケット高さ・放射によるチューンシフトが gs と違う。以下の表と本文の「pysad」の値と「一致」の判定は `PYSAD_COMPAT=k64` での結果 (2026-11-17 に確認)。gs と同じ値が要るときは `PYSAD_COMPAT=k64` を付ける。

- `Emittance[]` は、表の 7 行とも ○ 一致。
- 追跡の FLUC ON は △。同じ種 (`SeedRandom[7]`) でも、`TrackParticles` の量子励起の乱数列が gs と pysad で異なり、rms が 7% 違う。NOFLUC は一致。`Random[]`・`GaussRandom[]` の列そのものは一致する (08 の fn_Random_GaussRandom_SeedRandom) ので、差は放射の光子の引き方にある (原因は未調査)。

## 未確認・限界
- 追跡の粒子数・周回数が少なく (200 粒子、100 周)、統計の差かアルゴリズムの差かは判別していない。
- INTRA、POL、結合のあるラティスは見ていない。
