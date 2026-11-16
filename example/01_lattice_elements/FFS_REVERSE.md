# FFS_REVERSE (REVERSE)
種別: FFS コマンド。フォルダ: 01_lattice_elements。ラティスは START-ドリフト-QUAD-ドリフト-FIN の輸送路 (INS、初期 BX=BY=10 m、AX=1 AY=-1、EX=0.5 EPX=0.1)。エネルギー 3 GeV。

## 何を見る例か
SADHelp: 「`REVERSE` は入口の AX AY EPX EPY R2 R3 DPX DPY の符号を反転する」。進行方向を逆にして光学を伝搬させたいとき (出口の値から逆向きに追う) に使う想定のコマンド。
`REVERSE` の前後で、入口と出口の光学関数がどう変わるかを見る。

## スクリプト
```
USE T; INS; CALC;
Print["initial  AX AY EPX = ", Twiss["AX","START"], " ", Twiss["AY","START"], " ", Twiss["EPX","START"]];
Print["exit     AX BX AY BY = ", Twiss["AX","FIN"], ...];
REVERSE;   CALC;
Print["REVERSE: initial AX AY EPX = ", ...];   Print["REVERSE: exit    AX BX AY BY = ", ...];
```

## 結果の読み方
REVERSE の前後で出力は同じ (AX=1、AY=-1、EPX=0.1、出口 AX=-9.242220 BX=35.358212 AY=-105.9526 BY=676.90693)。

## 物理・解析式
なし。

## gs と pysad の比較
**一致するが、REVERSE は Twiss・MARK の値を変えない** (gs・pysad とも)。K64 のソース `tffsa.f` (`REV_ERSE` の分岐) を見ると、`axi ayi epxi epyi dpxi dpyi r2i r3i` という FFS 内の変数の符号を反転するだけで、MARK の保存値や Twiss の初期値を書き換えない。`FREE AXI` (入口値を変数にする) の変数にも反映されなかった (別に試した: 値は変化なし)。この K64 ビルドでは、REVERSE を単独で使っても光学の表示には現れない。SADHelp の記述どおりに働く例は作れなかった。

## 未確認・限界
- 上のソースの変数が FIT の `AXI` 等の条件式や `FFS$...` のどこで使われるかは追っていない。効く使い方 (マッチングの途中での使用など) は未確認。
