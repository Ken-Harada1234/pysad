# elem_SEXT_OCT_DECA_DODECA (6 極・八極・十極・十二極)
種別: elem。フォルダ: 01_lattice_elements。ラティス: 前半は 4 要素の輸送路 (INS)、後半は `FODO_ring.sad` と同じ光学のリング (定義はファイル内、RF なし)。エネルギー 3 GeV。

## 何を見る例か
多極磁石の要素の定義: `SEXT`(K2)、`OCT`(K3)、`DECA`(K4)、`DODECA`(K5)。**K2〜K5 は積分強度** (L を含む)。
(1) 輸送路で S3 O4 D5 D6 を並べると、線形の転送行列はドリフトと同じ (多極は 2 次以上の力だけ)。
(2) リングの SF/SD の K2 を変えても、設計軌道上 (COD=0) ではチューンも β も変わらず、**クロマティシティだけ**が変わる。

## スクリプト
```
SEXT  S3 = (L=0.2 K2=10);   OCT O4 = (L=0.2 K3=100);   DECA D5 = (L=0.2 K4=1000);   DODECA D6 = (L=0.2 K5=10000);
FFS;
USE T; INS; CALC;   TYPE S3 O4 D5 D6;   m = TransferMatrix["START","FIN"];
USE RING; CELL; CALC;   SF K2 5;   SD K2 -8;   CALC;   CHRO;
ps = LINE["POSITION","SF.1"] + {0,.25,.5,.75,1};   Print["SF.1 along s: BX = ", Twiss["BX",ps], ...];
```
`TYPE 名前` は要素の定義を表示。`SF K2 5;` は要素 SF の K2 を 5 にする (全 SF)。`CHRO` は要素ごとのクロマティシティ寄与 (次節)。全文は `.sad`。

## 結果の読み方
- 輸送路: `with S3 O4 D5 D6 : M11 M12 M21 M22 = 1 0.8 0 1` (L の合計 0.8 m のドリフト)。
- リング: K2=0 でも K2=(5,-8) でもチューンは 6.696329 / 4.906056 (同じ)。
- `CHRO;` の表 (1 要素 1 行、最後に `TOTAL`): x と y の列。SF.1 と SD.1 の行を下表で解析式と比べる。

## 物理・解析式
SAD の `CHRO` は四極と6 極の寄与を
$$\xi_{x}^{\rm (CHRO)}=-\frac{K_1}{L}\!\int\!\beta_x ds\ (\text{QUAD}),\quad +\frac{K_2}{L}\!\int\!\eta_x\beta_x ds\ (\text{SEXT}),\qquad
\xi_y^{\rm (CHRO)}=+\frac{K_1}{L}\!\int\!\beta_y ds,\ -\frac{K_2}{L}\!\int\!\eta_x\beta_y ds$$
で出す (これが $4\pi\xi$ の寄与: 実際の chromaticity は $\xi=\frac1{4\pi}\sum$)。6 極の符号は、SADHelp の式 ($-\frac{K_2}{L}\eta\beta$) とは逆で、コード・gs・pysad の値は上の符号 (物理的にも $\xi_x=-\frac1{4\pi}\oint(K_1-K_2\eta)\beta ds$ と合う)。
積分は FODO リングの周期解に沿って 2x2 行列で 200 分割して数値積分した (`tools/ana_fodo_ring.py`)。

| 行 | 式の値 (x, y) | gs (x, y) | pysad (x, y) |
|---|---|---|---|
| SF.1 (K2=5) | 15.9701, -6.78278 | 16.019, -6.788 | 16.019, -6.788 |
| SD.1 (K2=-8) | -4.8355, 18.8283 | -4.875, 18.804 | -4.875, 18.804 |

差は 0.3 % 以内 (SAD の積分の取り方と、こちらの台形則の差と思われる。原因までは調べていない)。

## gs と pysad の比較
- `CHRO` の表、チューン、β、行列はすべて一致。
- `TYPE D6` の `K5 =10000` も一致 (2026-11-17 に修正: pysad は DODECA の K5 を長さあたりの内部値 (50000) のまま返していた。`LINE["K5"]`・`Element["K5"]`・TYPE・ATTRIBUTE は積分強度 10000)。

## 未確認・限界
- OCT・DECA・DODECA の非線形の効果 (振幅依存のチューンシフトなど) は 05_tracking の例で扱い、ここでは確認していない。
- COD≠0 (ずれのある6 極) のときの線形光学への効き (フィードダウン) は扱っていない。
