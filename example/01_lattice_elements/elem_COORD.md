# elem_COORD (座標変換 COORD)
種別: elem。フォルダ: 01_lattice_elements (pysad に実装済みの範囲: docs/25)。ラティスは START-COORD-DRIFT-FIN の輸送路 (INS)。エネルギー 3 GeV (結果は無関係)。

## 何を見る例か
`COORD C1 = (DX=1E-3)`: ずれ・回転の座標変換そのものを 1 つの要素にしたもの (軸外の要素の表現などに使う)。`DX DY` (ずれ)、`CHI1 CHI2 CHI3` (回転)、`DIR` (0 で順変換、≠0 で逆変換)。
ここでは (a) DX=1 mm の COORD で閉軌道 (INS では初期の軌道 `DX`) がどう変わるか、(b) 初期軌道 DX=1 mm の粒子に `CHI3=0.01` の回転 (x-y 面内) をかけた結果、(c) 幾何 (GX GY GZ GCHI3) への効果を見る。

## スクリプト
```
COORD C1 = (DX=1E-3);   COORD C2 = (CHI3=0.01);
MARK START = (BX=10 BY=10 AX=0 AY=0) FIN = ();   MARK START2 = (BX=10 BY=10 AX=0 AY=0 DX=1E-3);
LINE T1 = (START C1 D1 FIN);   LINE T2 = (START2 C2 D1 FIN);
FFS;
USE T1; INS; CALC;   Print["C1 (DX=1mm): orbit DX DPX at FIN = ", Twiss["DX","FIN"], " ", Twiss["DPX","FIN"]];
```
先頭 MARK の `DX` は初期の閉軌道 (INS で伝搬する軌道) の位置。Twiss の `DX DPX DY DPY` は**軌道**で、分散は `EX EPX EY EPY` (混同しないこと)。

## 結果の読み方
- `C1 ... orbit`: COORD で軌道が -1 mm にずれる (x → x − DX)。
- `C2 ... orbit`: 初期 DX=1 mm の軌道が回転で DY にも出る。
- 幾何の行: ずれ・回転が設計軌道の幾何に入る (GX=1 は D1 の長さ 1 m、`GCHI3`=0.01)。

## 物理・解析式
SADHelp の COORD の式 ($DIR=0$): $(x,y,z)_1=R_3R_2R_1\,(x-DX,\,y-DY,\,z)_0$。DX だけなら $x\to x-DX$、CHI3 だけなら $(x,y)\to(x\cos\chi_3-y\sin\chi_3,\ x\sin\chi_3+y\cos\chi_3)$ ($\chi_3=0.01$、初期 $x=1$ mm、$y=0$)。

| 量 | 式の値 | gs | pysad |
|---|---|---|---|
| C1: 軌道 DX = -DX_c | -0.001 | -0.001 | -0.001 |
| C1: GX(FIN) = 1 (D1 の長さ) | 1 | 1 | 1 |
| C2: DX = 1mm·cos χ3 | 1.000e-03 | 1.000e-03 | 1.000e-03 |
| C2: DY = 1mm·sin χ3 | 1.000e-05 | 1.000e-05 | 1.000e-05 |
| C2: GCHI3 = χ3 | 0.01 | 0.01 | 0.01 |


## gs と pysad の比較
全桁一致。

## 未確認・限界
- `DY DZ CHI1 CHI2`、`DIR≠0` (逆変換)、COORD を含むリングでの Emittance・追跡 (docs/25 に gs との照合が書かれているが、この例では確認していない)、COORD のパラメータを FIT 変数にする場合は未確認。
- 追跡の低速経路 (track6) は COORD に未対応 (docs/25)。
