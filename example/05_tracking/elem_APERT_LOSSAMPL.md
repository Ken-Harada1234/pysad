# elem_APERT と LOSSAMPL
種別: 要素とグローバル変数 / フォルダ: 05_tracking / ラティス: `00_lattice/apert_ring.sad` (`FODO_ring.sad` に APERT を 2 個足した `RINGA`)、`00_lattice/FODO_ring.sad` (`RING`)。3 GeV / フラグ: `NORFSW NORAD NOFLUC`

## 何を見る例か
粒子を落とす 2 通りの方法: (1) ラインに置く `APERT` 要素 (楕円 `AX AY`、矩形 `DX1 DX2 DY1 DY2`)、(2) 全周のどこでも効く `LOSSAMPL` (|x|、|y| の上限、既定 1 m)。

## スクリプト
```
APERT AP0 = (AX=0.005 AY=0.005) AP1 = (DX1=-0.004 DX2=0.004 DY1=-0.004 DY2=0.004);
LINE RINGA = (ORIGIN AP0 RF 2*CELL AP1 22*CELL);     ! 00_lattice/apert_ring.sad
USE RINGA; CELL; CALC;
r = TrackParticles[{1, beam}, 1, 1, 1];              ! 7 粒子、1 周
USE RING; LOSSAMPL = 0.005; ...                      ! (B)
```
- (A) `AP0` は楕円の絞り (半径 5 mm): $(x/AX)^2+(y/AY)^2\le1$ を満たさない粒子が落ちる。`AP1` は矩形 ($DX1<x<DX2$ 等)。
- 粒子は $(x,y)$ [mm] = (0,0), (1,0), (3,0), (6,0), (−6,0), (0,4), (0,8)。

## 結果の読み方
| | (A) APERT in RINGA | (B) RING、LOSSAMPL=5 mm | (C) LOSSAMPL=1 m (既定) |
|---|---|---|---|
| gs の flag | {1,1,1,0,0,1,0} | {1,1,0,0,0,0,0} | すべて 1 |
| pysad の flag | {1,1,1,0,0,1,0} | {1,1,0,0,0,0,0} | すべて 1 |

## 物理・解析式
楕円の絞りの判定式 $(x/AX)^2+(y/AY)^2\le1$ で、(A) の結果が説明できる。AP0 での位置はほぼ初期値のままなので、7 粒子の判定は次のとおり。
- $(0,0)$、$(1,0)$: 通る。
- $(3,0)$: $0.36$、$(0,4)$: $0.64$ → 通る。
- $(\pm6,0)$: $(6/5)^2=1.44>1$、$(0,8)$: $(8/5)^2=2.56>1$ → 落ちる。

これは flag {1,1,1,0,0,1,0} と一致する。
(B) の LOSSAMPL は、判定点 (15 要素ごと) での位置で決まるので、解析式では示さない。

## gs と pysad の比較
(A)(B)(C) とも gs と一致 (○)。pysad の追跡 (numba 版 `track6n.py` の `k_apert_el`、Python 版 `track6.py` の `_apert_pass`) は、APERT 要素について `tapert.f` / `tapert1` の判定をそのまま行う。LOSSAMPL ((B)、(C)) は従来から一致している。
別に 16 粒子で、楕円・矩形・楕円+矩形・中心のずれ DX DY・ROTATE・全キーワード既定 0 (全粒子が落ちる)・APERT 2 個・3 周・JDPX/JDPY/DP の場合を gs と比べ、flag と落ちた時点の座標が一致した (`pysad/tests/test_apert.py`)。

## 未確認・限界
- ここでの比較は 3 GeV の FODO リング (RFSW/RAD なし、TrackParticles、1 周と 3 周) だけ。RAD・FLUC ありの追跡での落ち方、MAIN の TRACK の損失の報告 (`P. n lost in ... at l(APERT名)`)、TRPT は gs と比べていない。
- ソレノイド領域 (SOL BOUND≠0 の間) の中にある APERT は未対応 (pysad が拒否する。gs の tsol は tapert1 を呼ぶ)。
- JDPX/JDPY の再入射 (落とさず境界へ乱数で入れ直す) は、粒子が 11 個以下なら gs と同じ乱数列で一致するが、12 個以上では gs の値が一致しない (gs が別の乱数列を使うように見える。原因は未調査)。
- DP は JDPX/JDPY のどちらかが 0 でないときだけ効く (tapert.f のまま)。COUPLE は tapert が使わない。
