# fn_OrbitGeo_GeoBase (OrbitGeo と GeoBase)
種別: SADScript 関数。フォルダ: 02_optics_cod。ラティス: `00_lattice/errors_ring.sad` (ずれのある四極を含むリング。flag_COD 参照)。エネルギー 3 GeV。フラグ: CELL。

## 何を見る例か
- `OrbitGeo[位置]`: **閉軌道**の幾何 {GX, GY, GZ} (設計軌道の幾何 `LINE["GX"..]` に、閉軌道のずれを座標変換して足したもの)。軌道がなければ設計軌道の幾何と同じ。
- `GeoBase[{chi1, chi2, chi3}]`: 座標軸の回転角 (幾何の CHI1 CHI2 CHI3) を 3x3 の変換行列 `{{xgx,xgy,xgz},{ygx,ygy,ygz},{zgx,zgy,zgz}}` (局所の x・y・z 軸が大域座標でどの向きか) にする。

## スクリプト
```
USE RING; CELL; CALC;
Print["design GX GZ at B.1 exit (position 8) = ", LINE["GX",8], " ", LINE["GZ",8]];
Print["OrbitGeo[8] (closed orbit) = ", OrbitGeo[8]];
Print["OrbitGeo[$$$] = ", OrbitGeo[LINE["LENGTH"]]];
Print["GeoBase[{0,0,Pi/2}] = ", GeoBase[{0, 0, Pi/2}]];
Print["GeoBase[{0,Pi/6,0}] = ", GeoBase[{0, Pi/6, 0}]];
Print["GeoBase[{Pi/6,0,0}] = ", GeoBase[{Pi/6, 0, 0}]];
```

## 結果の読み方
- 設計軌道: B.1 出口 (位置 8) の GX = 2.286339 m、GZ=0。`OrbitGeo[8]` = (2.286352, 0.156199, -0.00157): 閉軌道がずれた分、設計の (GX, GY, GZ) = (2.286339, 0.1561845, 0) から少し動く。**GZ = -0.001566 m は位置 8 での閉軌道の垂直変位 DY** (flag_COD の模型で 0.001566)。
- `OrbitGeo[$$$]` (リング終端): 閉じる軌道なので、始点の軌道のずれと一致 (gs: (-5.338e-05, 1.284e-04, -3.354e-04)。GY、GZ は始点の軌道のずれ (flag_COD の位置 1 の DX=-1.284e-4、DY=3.354e-4) の符号を変えた値 (GY=-DX、GZ=-DY))。
- `GeoBase`: 3 つの行が直交する単位ベクトル。例えば `{0,Pi/6,0}` の 3 行目 (0.866, 0, 0.5) は大きさ 1、2 行目と直交。

## 物理・解析式
設計軌道 (偏向 θ、弦 ρ sin θ): 位置 8 の入口の GX = 1.1 + ρ sin θ = 2.286339 ✓ (gs: 2.286339)。回転行列の単位ベクトル性・直交性は次の確認のとおり (gs/pysad の行列要素から計算):
| 行列 (gs) | norm(行1) | norm(行2) | norm(行3) | 行1・行2 | 行2・行3 |
|---|---|---|---|---|---|
| GeoBase[{0,Pi/6,0}] | 1 | 1 | 1 | 0 | 0 |


## gs と pysad の比較
- `LINE["GX"]`、`OrbitGeo[8]`、`GeoBase` は全桁一致。
- `OrbitGeo[$$$]` も一致 (2026-11-17 に対応: GX = -5.338e-05、GY = 1.2836e-4、GZ = -3.354e-4)。原因は `OrbitGeo` の式 (Optics.n) が一周の経路長のずれ `Twiss["DZ",$$$]` (= -5.338e-5) を進行方向に足すことで、pysad は最後の行 ($$$) の DZ を 0 にしていた。$$$ の DZ は周回の出口の値にした (CELL の 4D の CALC)。

## 未確認・限界
- 垂直の曲がり (CHI2 に出る) や SOL を含む幾何、`OGEO` キー (`LINE["OGEO",k]` は OrbitGeo と同じ量) は未確認。GeoBase の行・列の規約は SADHelp の記述だけで、式との照合は直交性の確認にとどまる。
