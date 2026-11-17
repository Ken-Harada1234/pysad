# fn_SingularValues (SingularValues)
種別: 関数 (SADHelp functions の Matrix-operations、`*` つき = SAD 独自。`Packages/Matrix.n`) / フォルダ: 10_math_data / ラティスなし

## 何を見る例か
3×2 の行列 $A=\begin{pmatrix}3&1\\1&3\\0&2\end{pmatrix}$ の特異値分解。`SingularValues[a]` は `{u, w, v}` を返す。**既定では `Inverse->True` で、`w` は特異値の逆数 (昇順)**。`Inverse->False` にすると特異値そのもの (降順)。

## スクリプト
```
s=SingularValues[a];                    (* {u, 1/w, v}: 1/w は昇順 *)
t=SingularValues[a,Inverse->False];     (* {u, w, v}: w は降順 *)
Transpose[t[[1]]].DiagonalMatrix[t[[2]]].t[[3]]     (* = A *)
Transpose[s[[1]]].DiagonalMatrix[s[[2]]].s[[3]]     (* 既定: 疑似逆行列の転置 (A^+)^T *)
```
オプションは `Tolerance` (既定 1e-8、これ以下の特異値を 0 扱い)、`Inverse`、`Sort`。

## 結果の読み方
`dims=` は {{2,3},{2},{2,2}}。`u` は 2×3、`w` は長さ 2、`v` は 2×2 で、行が特異ベクトル (行列の形は `A=Transpose[u].DiagonalMatrix[w].v` となる並び)。`recon=` は元の $A$ に戻る (丸め誤差程度のずれ)。`pinvT=` は $(A^+)^T$。

## 物理・解析式
$A=U\,\mathrm{diag}(w)\,V^T$、特異値は $A^TA=\begin{pmatrix}10&6\\6&14\end{pmatrix}$ の固有値の平方根。固有値は $12\pm2\sqrt{10}=18.325,\ 5.675$ なので $w=4.2807,\ 2.3823$。

| 項目 | 式 (numpy.linalg.svd) | gs | pysad | gs の相対差 | pysad の相対差 (参照値が 0 の行は絶対差) |
|---|---|---|---|---|---|
| 特異値 $w_1$ (Inverse->False) | 4.28071901908 | 4.280719019083 | 4.280719019083 | 5.9e-14 | 5.9e-14 |
| 特異値 $w_2$ | 2.3823191809 | 2.382319180896 | 2.382319180896 | 1.5e-13 | 1.5e-13 |
| 既定の出力 $1/w_1$ (小さい 1/w が先) | 0.233605615211 | .233605615211 | .233605615211 | 1.7e-12 | 1.7e-12 |
| 既定の出力 $1/w_2$ | 0.419759034818 | .419759034818 | .419759034818 | 1.2e-12 | 1.2e-12 |

## gs と pysad の比較
- 特異値・特異ベクトルは 12 桁で一致 (○)。特異ベクトルの符号は pysad が gs に合わせてある (各 u の先頭成分が負)。
- `recon` の (3,1) 要素は 0 の丸め (gs 2.2E-16、pysad 4.18E-17) で一致しない。

## 未確認・限界
特異値が `Tolerance` 以下のランク落ち行列、複素行列 (K64 は対応、pysad は未確認)、200 次元超の行列は試していない。この例でも、特異ベクトルの符号を gs に合わせた pysad の実装 (3 行列で確認) 以外の行列での符号の一致は未確認。
