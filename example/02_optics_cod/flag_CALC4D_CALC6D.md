# flag_CALC4D_CALC6D (CALC4D と CALC6D)
種別: フラグ。フォルダ: 02_optics_cod。ラティス: `00_lattice/FODO_ring.sad` (RF 1.5 MV の `CAVI` 入り)。エネルギー 3 GeV。

## 何を見る例か
`CALC4D` (既定): `CALC` は 4 次元 + 分散 (4x5 行列) の光学。`CALC6D`: 縦方向も含めた 6 次元の完全な計算 (RADCOD を入れた軌道にもできる)。`MATRIX` の行列が 4x5 から 6x6 に変わり、`Twiss` の `NZ BZ AZ DZ DDP` が意味を持つ。
ここではどちらも `CALC` を走らせ、6 次元の 1 周行列から運動量圧縮率などを読む。

## スクリプト
```
USE RING; CELL; CALC;   em = Emittance[];   Print["CALC4D: Tunes ... ", Tunes/.em];
RFSW;   CALC6D;   CALC;
Print["CALC6D: NX NY NZ /2pi = ", Twiss["NX","$$$"]/2/Pi, " ", Twiss["NY","$$$"]/2/Pi, " ", Twiss["NZ","$$$"]/2/Pi];
MATRIX ORIGIN $$$;
```

## 結果の読み方
- CALC4D: `NZ` = 0、`BZ` = 1 (縦の情報なし)、チューンは (6.696329, 4.906056)。`Emittance[]` の `Tunes` の 3 つ目 (シンクロトロンチューン) は 0 (RF 1.5 MV は放射損失 1.56 MV 未満で平衡なし。elem_CAVI 参照)。
- CALC6D: NX NY は 4D と同じ。**`NZ` = -0.18567 は意味のある値ではない**: 1 周の 6x6 行列を見ると縦方向の ($z$, δ) ブロックは $\begin{pmatrix}1&-2.418\\0&1\end{pmatrix}$ で、RF による (6,5) 要素が 0 (RF の効果が入っていない、純粋なドリフト)。`CALC6D` は RF を入れない (RFSW が効くのは追跡と Emittance[] で、`CALC6D` の行列にはない) と思われる。
- 1 周行列の x ブロック: (1.713421, -4.742148; 1.069062, -2.375161)、(5,1)(5,2) と (1,6)(2,6) は分散の項。

## 物理・解析式
1 周行列は、チューンと運動量圧縮率で書ける:
$$\mathrm{tr}M_x=2\cos2\pi\nu_x,\quad \mathrm{tr}M_y=2\cos2\pi\nu_y,\quad M_{56}=-\alpha_c C\ (C=105.6\ \text{m},\ \alpha_c=0.0228977)$$
| 量 | 式の値 | gs | pysad |
|---|---|---|---|
| tr Mx | -0.66174 | -0.66174 | -0.66174 |
| tr My | 1.66158 | 1.66158 | 1.66158 |
| M56 = -α_c C | -2.418 | -2.418 | -2.418 |
| det Mx | 1 | 0.999999 | 0.999999 |

(α_c = 0.0228977 は `Emittance[]` の `MomentumCompaction` (elem_CAVI の gs の表示))

## gs と pysad の比較
全桁一致 (行列の雑音 1e-16 〜 1e-29 の要素だけ表示が違う。AZ の雑音: gs 6.0e-10、pysad -2e-16)。

## 未確認・限界
- RF の収束を含む 6D の周期解 (`CALC6D` のあとの `NZ`) は、この CAVI 設定では得られなかった。SAD でどの条件のときに RF が 6x6 に入るか (RADCOD、`Emittance[]` 経由など) は未確認。シンクロトロンチューンは `Emittance[]` (elem_CAVI) で確認済み。
- 6D の Twiss (BZ AZ DZ DDP の長さ方向の関数) の物理的な意味は、RF が入る条件を確かめてから見る必要がある。
