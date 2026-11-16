# elem_DRIFT_MARK (DRIFT と MARK)
種別: elem (要素定義)。フォルダ: 01_lattice_elements。ラティスは 4 要素だけの輸送路 (リングではない)。エネルギー 3 GeV。フラグ: INS (初期値を MARK から与える)。

## 何を見る例か
`DRIFT` (長さ `L` だけの自由空間) と `MARK` (長さ 0 の印。**ラインの先頭の MARK は入射ビームの Twiss を持つ**) の例。
先頭の `MARK START = (BX=10 BY=10 AX=0 AY=0)` が初期条件、途中の `MID`、終端の `FIN` は光学関数を読む場所の目印。

## スクリプト
```
DRIFT D5 = (L=5);
MARK START = (BX=10 BY=10 AX=0 AY=0) MID = () FIN = ();
LINE T = (START D5 MID D5 FIN);
FFS;
USE T; INS; CALC;
DISP;
Print["BX START,MID,FIN = ", Twiss["BX",{"START","MID","FIN"}]];
```
要素の定義は `FFS;` の前 (MAIN 部) に置く。`USE T; INS; CALC;` は T を使い、周期解ではなく MARK の初期値から伝搬する (`INS`)。
`DISP;` は各要素入口の光学関数を表で出す (表の行は「その要素の入口」の値)。全文は `elem_DRIFT_MARK.sad`。

## 結果の読み方
- `Matched. ...` の次の 3 行は、ライン末尾 (`$$$`) の AX BX NX AY BY NY LENG。`LENG` = 10 m。
- `DISP;` の表: 行 `MID` の BX=12.5、`FIN` の BX=20 (どちらの平面も同じ)。`s(m)` の列で位置が分かる。
- `Print` の行は `Twiss[関数名, 位置の名前のリスト]` の戻り値。`LINE["S", ...]` は要素入口の s。

## 物理・解析式
ドリフトでは $\gamma_0=(1+\alpha_0^2)/\beta_0$ として
$$\beta(s)=\beta_0-2\alpha_0 s+\gamma_0 s^2,\quad \alpha(s)=\alpha_0-\gamma_0 s,\quad \psi(s)=\arctan\frac{s}{\beta_0-\alpha_0 s}$$
($\alpha_0=0$、$\beta_0=10$ m)。

| 量 | 式の値 | gs | pysad |
|---|---|---|---|
| BX(MID) (s=5 m) | 12.5 | 12.5 | 12.5 |
| BX(FIN) (s=10 m) | 20 | 20 | 20 |
| AX(MID) | -0.5 | -0.5 | -0.5 |
| AX(FIN) | -1 | -1 | -1 |
| NX(FIN) (rad) | 0.785398 | 0.785398 | 0.785398 |

(NX は位相進み [rad]。`FIN` で $\arctan(10/10)=\pi/4$)

## gs と pysad の比較
全桁一致 (出力の数値トークンはすべて同じ)。

## 未確認・限界
- 光学関数は 1 次 (線形) の話だけ。MARK の `OFFSET` や `GEO` キーワードは別の例 (FFS_SPLIT、FFS_ORG) で扱う。
- INS の初期値はエネルギーによらないので、3 GeV の仮定は結果に効かない。
