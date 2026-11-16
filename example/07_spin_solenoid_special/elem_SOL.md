# elem_SOL
要素・07_spin_solenoid_special・ラティス `00_lattice/sol_ring.sad` (新規: FODO_ring の 1 セル目の D1 (1 m) をソレノイド領域 `SA` ... `SB` にした RINGS。BZ=5 T)、3 GeV、`NORFSW NORAD NOFLUC`。参考: `docs/42_sol_guide.md`

## 何を見る例か
SOL は「ソレノイド領域の境界」を表す長さ 0 の要素で、`BOUND=1` の SOL を対にして置く。最初の SOL から次の BOUND≠0 の SOL までが領域で、領域の中は縦磁場 $B_z$ で運動する。1 つの領域 (`DRIFT` だけ) を通る 4×4 行列を追跡で求め、ソレノイドの伝達行列 (結合) の解析式と照合する。

## スクリプト
```
SOL SA = (L=0 BZ=5 BOUND=1 GEO=1) SB = (L=0 BZ=0 BOUND=1);
LINE CELLS = (SF DS QF D2 B D2 QD DS SD SA D1 SB);
LINE RINGS = (ORIGIN RF CELLS 23*CELL);
USE RINGS; CELL; CALC;                       ! 結合したリングのチューン
beam = (x, px, y, py をそれぞれ a=1e-6 にした 4 粒子);
r = TrackParticles[{pa, beam}, pb+1, 1, 1];   ! SA の位置から SB の次まで (領域を 1 回)
```
- `GEO=1` は対の片方 (`SA`) だけに付ける。両方に付けると K64 は `Duplicated GEO of Solenoid` を出して領域を扱わない。付けないと S 座標が途切れて周長を使う量が狂う (K64 の癖、docs/42 §5)。
- `L=0` にする (領域の長さは中の DRIFT・QUAD の長さの和)。領域内に置けるのは DRIFT と QUAD、領域内の SOL (BOUND=0) まで。
- この例は領域の中に四極を置かない。SOL 内の四極 (`tsolqu`) に `RAD` を組み合わせると、gs の結果は粒子の並びに依存する癖があり (docs/42 §5、pysad は再現しているが gs の現行ビルドと参照の保存値が食い違うことがある)、避けた。

## 結果の読み方
- `SA position=12 SB position=14 MOMENTUM=3e9`
- `with solenoid`: 領域のあるリングの `NX=6.69418751 NY=4.97467992`。ソレノイドなしは (6.696329, 4.906056)。結合のため整数部以外が変わる。
- `M[i,j]`: 領域の 4×4 行列の行 $i$ (列 $j$ は $x,p_x,y,p_y$ に $a$ を与えた粒子の応答)。

## 物理・解析式
一様な縦磁場 $B_z$ (長さ $L$) の中の運動 (力学的運動量 $p_x,p_y$ で、領域の外は磁場なし。端のフリンジ込み)。$K=B_z/(2B\rho)$、$B\rho=p/e=10.0069$ T·m (3 GeV/c)、$C=\cos KL$、$S=\sin KL$:

$$M=\begin{pmatrix}C^2&SC/K&SC&S^2/K\\-KSC&C^2&-KS^2&SC\\-SC&-S^2/K&C^2&SC/K\\KS^2&-SC&-KSC&C^2\end{pmatrix}$$

($x$ と $y$ の結合は $\pm SC$、$S^2$ の項)。$B_z=5$ T では $K=0.249827$ m$^{-1}$、$KL=0.249827$:

| 要素 | 式 | gs = pysad |
|---|---|---|
| $M_{11}=M_{33}=M_{22}=M_{44}$ ($C^2$) | 0.93887417 | 0.93887417 |
| $M_{12}=M_{34}$ ($SC/K$) | 0.95890728 | 0.95890728 |
| $M_{13}=-M_{31}$ ($SC$)、$M_{24}$、$-M_{42}$ | 0.23956098 | 0.23956098 |
| $M_{14}=-M_{32}$ ($S^2/K$) | 0.24467258 | 0.24467258 |
| $M_{21}=M_{43}$ ($-KSC$) | −0.05984881 | −0.05984881 |
| $M_{23}=-M_{41}$ ($-KS^2$) | −0.01527089 | −0.01527089 |

16 要素すべてが 8 桁で一致する (gs と pysad は印字桁まで同じ、式 (Python) とも 8 桁一致)。

## gs と pysad の比較
○ 一致。全 7 行が印字桁まで同じ (チューン、4×4 行列の 16 要素)。

## 未確認・限界
- 領域は DRIFT だけ (BZ=5 T、L=1 m、3 GeV)。BEND・MULT・CAVI が領域内にあるリング、BOUND=0 で BZ≠0 の SOL、DX・DY・CHI のずれは pysad 未対応 (docs/42 §1)。
- 結合したリングのチューン (6.69418751, 4.97467992) は、解析式で確認していない (gs と pysad の一致のみ)。
- RAD ON の領域内の四極、放射とスピン (POL) の組み合わせは見ていない。
- SOL を含む CALC6D と RADCOD つき CALC は pysad で拒否される (docs/42)。
