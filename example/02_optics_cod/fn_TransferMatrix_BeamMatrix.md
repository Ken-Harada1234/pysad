# fn_TransferMatrix_BeamMatrix (TransferMatrix と BeamMatrix)
種別: SADScript 関数。フォルダ: 02_optics_cod。ラティス: `00_lattice/FODO_ring.sad` (RF 6 MV に変更)。エネルギー 3 GeV。フラグ: CELL、CALC4D、EMIT 用に RFSW・RAD・NOFLUC。

## 何を見る例か
- `TransferMatrix[from, to]`: from から to までの転送行列 (CALC4D なら 4x4、CALC6D なら 6x6)。`to` が上流なら逆行列に相当 (リングでは 1 周回った行列)。オプション `Fold->`、`Symplectic->` (SADHelp)。
- `BeamMatrix[位置]`: その位置の 6x6 ビーム行列 $\Sigma$ (`EMIT` が必要)。

## スクリプト
```
m = TransferMatrix["SF.1","SF.2"];     Print["row1 = ", m[[1]]];     Print["row2 = ", m[[2]]];
mi = TransferMatrix["SF.2","SF.1"];    Print["... M*Minv = ", (m.mi)[[1,1]], " ", (m.mi)[[1,2]]];
RF VOLT 6.0E6; CALC;  RFSW; RAD; NOFLUC;  EMIT;
b = BeamMatrix[1];     Print["sqrt(b11) sqrt(b33) sqrt(b66) = ", Sqrt[b[[1,1]]], " ", Sqrt[b[[3,3]]], " ", Sqrt[b[[6,6]]]];
```

## 結果の読み方
`row1`、`row2` は 1 セル (SF.1 の入口から SF.2 の入口) の x ブロックの 1・2 行 (x-y 結合の 0 は 1e-16)。`M*Minv` の (1,1)=1、(1,2)=0 で逆行列を確認。
`BeamMatrix` は 6x6、`b11` = σx²。

## 物理・解析式
1 セルの周期行列 $M_{\rm cell}$ (`tools/ana_fodo_ring.py`):
| 要素 | 式の値 | gs | pysad |
|---|---|---|---|
| M11 | -2.3117 | -2.3117 | -2.3117 |
| M12 | 4.94192 | 4.94192 | 4.94192 |
| M21 | -1.1141 | -1.1141 | -1.1141 |
| M22 | 1.94912 | 1.94912 | 1.94912 |

ビーム行列は $\Sigma_{11}=\epsilon_x\beta_x+(\eta_x\sigma_\delta)^2$。**`BeamMatrix` の δ 方向の広がりは `EMIT` で求めた SIGE ではなく、特殊変数 `DP` (=0.01) を使う**ことが値から分かる: $\sqrt{b_{66}}$ = 0.01 = DP。
| 量 | 式の値 | gs | pysad |
|---|---|---|---|
| b11 (σδ=DP) | 3.282e-05 | 3.282e-05 | 3.282e-05 |
| b11 (σδ=SIGE) | 1.327e-06 | 3.282e-05 | 3.282e-05 |

(上段: EMITX·BX+(EX·DP)² を計算した値、下段: EMITX·BX+(EX·SIGE)²。BeamMatrix の b11 は上段と合う = DP を使っている。)

## gs と pysad の比較
**v2 の注記 (2026-11-17、v2.0.0 の再実行)**: `.pysad.txt` は既定の v2 モードの出力に作り直した。v2 は `Emittance[]`・EMIT・追跡の Jx の端の項 (B75) と `MomentumCompaction` (B77) を直したので、この例の pysad の値のうち εx (約 2.4% 大きい)・εz・減衰率の X と Z・σδ・バンチ長・α・バケット高さ・放射によるチューンシフトが gs と違う (K64 のバグを直したため。`docs/47_v2_changes.md`)。以下の表と本文の「pysad」の値と「一致」の判定は、v1 と同じ `PYSAD_COMPAT=k64` での結果 (旧い出力と同じことを 2026-11-17 に確認)。gs と同じ値が要るときは `PYSAD_COMPAT=k64` を付ける。

全桁一致 (行列、BeamMatrix の値)。EMIT の表の `alpha` の桁だけ違う (0.0222107 と 0.02221)。

## 未確認・限界
- CALC6D の 6x6 の `TransferMatrix`、`Fold->True`、`Symplectic->True` は未確認。RAD 付き 6D ビーム行列の平衡値 (`Emittance[...,Matrix->True]` の `EquilibriumBeamMatrix`) は 04_emittance_radiation の例で扱う。
