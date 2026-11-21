# SetElement_BeamLine_Join_ASTRID2 (周期性のあるラティスの作り方: 6 極だけ全セル独立にする)
KEK の SAD チュートリアルの「SAD スクリプト上でラティスを定義する例」をもとにした例(`16*NCELL` のように作ったリングで、6 極だけ全セルで別の値にしたい。LINE を全部書き下さず、`SetElement`・`BeamLine`・`Join` で組む)。教材は別の 16 セルのリングだが、ここは ASTRID2 (0.58 GeV、周長 45.7044 m、6 セル、チューン 5.17999 / 2.14018) の要素で書き直した。ラティス: `lattice/ASTRID2_v01_01.sad` (`RING = (3*NCELL CAVI 3*NCELL)`)。フラグ: `RFSW RADCOD RAD COD` (`FLUC` は OFF)。

## 何を見る例か
`6*NCELL` のリングでは同じ要素 (SF・SD) を 6 回使うので、セルごとに K2 を変えられない。そこで6 極だけを `SetElement` でセルごとに別の名前の要素として作り (SDA01…SDB06、24 台)、6 極以外の部分を `BeamLine[...]` でくくり、`Join` でつないで 6 セルのリングを作る。元の `RING` と同じ光学になること (確認) と、1 つのセルの6 極だけを変えられること (応用) を見る。

## スクリプト
```
aa = {"01","02","03","04","05","06"};              ! セル番号の文字列 (名前の連結用)
Do[ id = aa[j];
    SetElement["SDA"//id, "SEXT", {"L"->lsd}, {"K2"->ksd}];   (* SDA01 ... : 要素を新しく作る *)
    SetElement["SFA"//id, "SEXT", {"L"->lsf}, {"K2"->ksf}];   (* SFA, SFB, SDB も同様 *)
,{j,1,6}];
bl1 = BeamLine[MKS, D14837, PM, D0100, QF3, D0180, PH, PV, D0180, BH, D0190];   ! 6 極 (SDA) の手前まで
bl2 = BeamLine[D0120, QF2, D0120];  bl3 = BeamLine[PM];  bl4 = BeamLine[D0120, QF2, D0120];
bl5 = BeamLine[D0190, BH, D0180, PH, PV, D0180, QF3, D0100, PM, D14837, MKE];
FFS["USE BeamLine[MKS];"];  bl = ExtractBeamLine[];            ! 出発点: 要素 1 つの BeamLine
Do[ id = aa[j];
    bl = Join[bl, bl1]; bl = Join[bl, BeamLine["SDA"//id]]; bl = Join[bl, bl2]; bl = Join[bl, BeamLine["SFA"//id]];
    bl = Join[bl, bl3]; bl = Join[bl, BeamLine["SFB"//id]]; bl = Join[bl, bl4]; bl = Join[bl, BeamLine["SDB"//id]];
    bl = Join[bl, bl5];  If[j == 3, bl = Join[bl, BeamLine[CAVI]]];                ! 3 セルのあとに RF
,{j,1,6}];
FFS["USE Join[bl];"];
```
- 1 セルの要素の並びは `NCELL` と同じ: `MKS D14837 PM D0100 QF3 D0180 HV D0180 BH D0190 SD D0120 QF2 D0120 SF PM SF D0120 QF2 D0120 SD D0190 BH D0180 HV D0180 QF3 D0100 PM D14837 MKE` (`HV = (PH PV)` はベンドの回転 0 と 90 度の 2 要素を `PH, PV` と並べた)。
- `Join` でつないだ `BeamLine` を `USE Join[bl]` でラインとして使う。
- **名前の連結に `ToString[j]` や `"0"//id` を使うと、名前が `SFA0` に切れて全セルが同じ名前になった** (gs・pysad とも。`SetElement` が前の要素を上書きし、`LINE["K2","SFA03"]` は空になって `Invalid set to LINE` で失敗した)。教材のように、文字列のリスト (`aa`) から取り出して連結すると、`SFA01`〜`SFA06` の別の名前になる。
- 元の `RING` が 200 要素なのに、新しいリングは 201 要素なのは、先頭に `USE BeamLine[MKS]` の出発点の `MKS` が 1 つ余計に入るため (教材の `BeamLine[MKNS]` から始める作り方と同じ)。

## 結果の読み方
`RING` と `NEWRING` で始まる行 (チューン、要素数、εx)。`xi (all same K2)` は全部同じ K2 のときのクロマティシティ (`DP0=±1e-4` の差分)。`formula:` と `SF of cell 3 -> 0 : xi` は、セル 3 の SF 2 台 (SFA03・SFB03) の K2 を 0 にしたときの ξ の変化 (式と計算)。

| 量 | 元の RING | 組み直した NEWRING |
|---|---|---|
| NX, NY | 5.179990, 2.140178 | 5.179990, 2.140178 (小数 8 桁まで同じ) |
| 要素数 | 200 | 201 |
| εx [m] | 1.268037e-8 | 1.268037e-8 (比 1.000000) |

## 物理・解析式
6 極 (薄肉) が与える色収差は
$$\xi_x^{\rm sext}=+\frac1{4\pi}\sum_i K_{2,i}\,\eta_i\,\beta_{x,i},\qquad \xi_y^{\rm sext}=-\frac1{4\pi}\sum_i K_{2,i}\,\eta_i\,\beta_{y,i}$$
($K_2$ は SAD の積分強度、$\eta$・$\beta$ は6 極の中点の値。符号は SAD の定義に合わせた。`mom_dep_ts_ASTRID2` の式と同じ。)
セル 3 の SF 2 台を外すと、ξ はその 2 台の寄与分だけ変わる:

| 量 | 式 | gs | pysad |
|---|---|---|---|
| 全部同じ K2 のときの $\xi_x$, $\xi_y$ (色収差補正済みで ~0) | | -2.10e-4, -2.79e-4 | 同じ |
| SF(セル 3) を 0 にしたときの $\Delta\xi_x$ | -2.0649 | -2.0664 (= -2.0666 - (-0.0002)) | 同じ |
| 同 $\Delta\xi_y$ | +2.9691 | +2.9692 (= 2.9689 - (-0.0003)) | 同じ |

式と計算が 0.1 % 以内で一致する。つまり、6 極をセルごとに別の要素にして、1 つのセルだけ変えられることが確かめられた (SFA03 だけ 0 に、他の 5 つのセルの `LINE["K2","SFA*"]` は 8.68 のまま)。

## gs と pysad の比較
○ 一致 (表のすべての値が全桁一致)。

## 未確認・限界
- 教材のリングは 16 セル (2 種類のセルを前半 10 / 後半 6 に使う特殊なリング) で、6 極は 5 種類 × 2 (反対向き) × 10 セル。ここは ASTRID2 の 6 セル (6 極 4 種類: SDA・SFA・SFB・SDB) に書き直した。教材の `-NCELLHA` (反転したライン) と `MULT` 要素への置換 (`SetElement[...,"MULT",{"L"->...},{"K2"->...}]`) は使っていない (ここは `"SEXT"`)。
- 作ったのは線形・色収差の確認まで。6 極を独立にしたあとの DA や多極成分の独立化 (`MULT`) は試していない。
- 要素名の長さ・文字の制約 (名前に `_` を含めるとき等) は調べていない。ここでは `SDA01` のような英数字 5 文字の名前を使い、`LINE["K2","SFA03"] = 0` で個別に値を変えられた。
