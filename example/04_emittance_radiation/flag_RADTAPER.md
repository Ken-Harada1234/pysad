# flag RADTAPER
種別: フラグ / フォルダ: 04_emittance_radiation / ラティス: `00_lattice/FODO_ring_V4MV.sad` (FODO_ring と同じ磁石列、RF 電圧だけ 4 MV。QF から始まる `RINGQ` を使う)。3 GeV、周長 105.6 m、24 セル、RF 500 MHz/ フラグ: RFSW, RAD, RADCOD, RADTAPER。

- **ラティスについて**: `FODO_ring.sad` の RF は 1.5 MV で、エネルギー損失 $U_0=1.563$ MV より小さく、バケットがない (同期位置 dz=150 mm、$\nu_s$ 虚数)。そこで **RF 電圧 4 MV の `FODO_ring_V4MV.sad` を新しく追加**した。
- RADINT は先頭の要素が SEXT (L=0.2, K2=0) だと `DL=0. IN TCAL` で I2・I3・I5 が NaN になる (gs も pysad も同じ)。そのため先頭が QF の `RINGQ = (ORIGIN RF 24*CELLQ)` を使う (磁石の並びは同じ、位相だけ違う)。
- SAD の BEND の `E1`, `E2` は ANGLE に対する比なので、端の角は $\epsilon=E1\cdot\text{ANGLE}=0.0343$ rad (θ/2 ではない)。四極の `K1` は積分強度 (k=K1/L)。

## 何を見る例か
RADTAPER を立てると、閉軌道の局所運動量 (DDP) に合わせて磁石が(ソレノイド以外)スケールされる。RADTAPER の前後で、閉軌道・$U_0$・チューンがどう変わるかを `Emittance[]` で見る。

## スクリプト
```
RFSW; RAD; RADCOD;
e0 = Emittance[];            (* RADTAPER なし *)
RADTAPER;
e1 = Emittance[];            (* RADTAPER あり *)
Print[Twiss["DDP", {1, 60, 121, 243}]];
(* 原点から 2 周 TrackParticles *)
```

## 結果の読み方 (gs)
| | RADTAPER なし | RADTAPER あり |
|---|---|---|
| OrbitAtExit x | -7.887E-5 | -2.198E-7 |
| OrbitAtExit δ | -2.624E-4 | 3.5126E-9 |
| $U_0$ [eV] | 1563263.2 | 1564904.0 |
| $\epsilon_x$ | 1.7358E-7 | 1.7372E-7 |
| Tunes | -.303702, -.093902, -.027625 | -.303633, -.094058, -.027619 |

DDP (閉軌道の運動量偏差、一周で $U_0$ ぶんだけ下がり RF で戻る) は {0, 3.9113E-4, 2.6068E-4, 3.5126E-9} (要素番号 1, 60, 121, 243)。RADTAPER なしの閉軌道は、損失エネルギーでずれた位置 (x=-79 μm、δ=-2.6E-4) を通り、RADTAPER ありでは磁石が運動量に追従するので軌道が原点近くに戻る (x=-0.2 μm)。

## 物理・解析式
エネルギー損失の比 $U_0/E_0=5.2\times10^{-4}$ に対し、一周の途中の運動量偏差は最大 $\delta\approx4\times10^{-4}$ (DDP の値)。磁石の強さを局所運動量に合わせてスケールするので、軌道のずれ (x=-79 μm) が 0.2 μm に減り、$U_0$ は 1641 eV (0.1%) 増え、チューンが $10^{-4}$ の桁で動く。$\Delta\nu\approx\xi\,\delta_{\rm local}$ (色収差 $\xi_x=-2.7$、$\delta\sim3\times10^{-4}$ で $10^{-3}$ の桁) より小さいのは、磁石が追従して補正されるため、と推定する (式の値としての照合は未実施)。

## gs と pysad の比較
**× pysad は未対応** (フラグを保持するだけで磁石をスケールしない)。pysad の値は RADTAPER なしと同じ ($U_0$ 1563263.2、OrbitAtExit x -7.887E-5、DDP 0)。`.pysad.txt` はそのまま出力 (エラー表示は出ず、値が変わらない)。docs/14 と docs/16 に「RADTAPER [未]」とある。

## 未確認・限界
- gs の RADTAPER ありの `Twiss["DDP"]` を要素番号で引いている (名前 `QF.1` では 0 が返った)。
- 追跡 (2 周) の結果は載せていない部分があるので .sad の `Print` を参照。gs と pysad の追跡は RADTAPER の有無で分かれるため比較できない。
- 運動量偏差が小さい (3 GeV、$U_0/E=5\times10^{-4}$) ので効果は小さい。大きな損失のリングでの確認はしていない。
