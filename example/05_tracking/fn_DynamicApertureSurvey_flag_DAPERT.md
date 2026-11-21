# fn_DynamicApertureSurvey と flag_DAPERT
種別: 関数とフラグ / フォルダ: 05_tracking / ラティス: `00_lattice/FODO_ring_V4MV.sad` (RF 4 MV、3 GeV) / フラグ: `RFSW; RAD; NOFLUC; EMIT; NORAD`。50 周

## 何を見る例か
`DynamicApertureSurvey[range, nturn, Output->6]` は、閉軌道・規格化座標・平衡エミッタンス (先に EMIT 等が必要。[multi-turn-tracking](multi-turn-tracking.md)) を使い、規格化振幅の格子で多周回追跡して、粒子が生き残る周数の地図 (1 周 = 1 文字) と Score を返す。`DAPERT` フラグはこの手順 (DAPERT 手続き) をマクロ粒子追跡に対して有効にするフラグで、`DynamicApertureSurvey` の内部が使う (利用者が ON/OFF する必要は通常ない)。

## スクリプト
```
RFSW; RAD; NOFLUC;
EMIT;                       ! RF 4 MV で U0=1.56 MV を補い、平衡エミッタンスを決める
NORAD;
DAPWIDTH = 51;              ! 追跡を始める粒子の窓の幅 (既定 7)
res = DynamicApertureSurvey[{{1,400},{1,1},{-0.5,0,0.5}}, 50, Output->6];
```
- `range = {xrange, yrange, zrange}`。既定の Axes `"ZX"` では zrange に値のリスト (ここでは $z$ 方向の規格化振幅 −0.5, 0, 0.5)、xrange は {最小, 最大} で 51 点に分けて走査する。値は初期振幅を $\sqrt{(\varepsilon_x+\varepsilon_y)}$ (ここでは $\sqrt{\mathrm{EMITX}+\mathrm{EMITY}}=\sqrt{1.7358\times10^{-7}}=4.17\times10^{-4}$) で割った規格化値。
- `Output->6` で標準出力に地図を出す。戻りは `{Score, 走査結果}`。

## 結果の読み方
```
 Turns = 50  Maximum number of particles = 306
      NZ     0----|----1----|----2----|----3----|----4----|----5
    -0.50 40 UUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUNIUaH1  1 .
     0.00 39 UUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUURMNja82  1 .
     0.50 39 UUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUsqNOA61  1 .
    Score:   118
```
(上は gs。文字は損失した周を表す符号で、`U` は 50 周生き残り、左から $x$ の規格化振幅が大きくなる。)
行頭の数 (40, 39, 39) は、$z$ の規格化振幅 (−0.5, 0, 0.5) ごとの score で、50 周生き残った $x$ 振幅の点の数 (51 点中)。戻り値の Score は各 $z$ の score の和 $40+39+39=118$ (SADHelp)。地図の左端から続く `U` が 50 周生き残る範囲で、右の文字は損失した周 (数字や英字) を表す。
$y$ の初期振幅は、`yrange={1,1}` により一定の規格化振幅 1 で与える ($x$ の振幅には比例しない)。EMITX=1.7358e-7 m、SIGE=1.1914e-3。

## 物理・解析式
動的口径は非線形で決まるので解析式はない。規格化振幅 $n$ の水平の初期振幅は、定義 $n=\sqrt{2J_x/(\mathrm{EMITX}+\mathrm{EMITY})}$ から $x_0=n\sqrt{(\mathrm{EMITX}+\mathrm{EMITY})\beta_x}$ ($\beta_x=5.025$ m、ORIGIN) となり、1σ 相当は $9.3\times10^{-4}$ m。51 点は $n=1$ から 400 を 51 等分 (間隔 7.98) なので、40 点目 ($n\approx312$) までが 50 周生き残り、$x_0\approx0.29$ m、角度 $\sim0.3$ rad に相当する。
ここまで大きいのは、このラティスの6 極が $K2=0$ で、非線形が厳密なドリフト・エッジの運動学だけだから。この限界は大振幅の運動学的非線形で決まり (解析式の裏付けは無い)。

## gs と pysad の比較
**v2 の注記 (2026-11-17、v2.0.0 の再実行)**: `.pysad.txt` は既定の v2 モードの出力に作り直した。v2 は `Emittance[]`・EMIT・追跡の Jx の端の項 (B75) と `MomentumCompaction` (B77) を直したので、この例の pysad の値のうち εx (約 2.4% 大きい)・εz・減衰率の X と Z・σδ・バンチ長・α・バケット高さ・放射によるチューンシフトが gs と違う (K64 のバグを直したため。`docs/47_v2_changes.md`)。以下の表と本文の「pysad」の値と「一致」の判定は、v1 と同じ `PYSAD_COMPAT=k64` での結果 (旧い出力と同じことを 2026-11-17 に確認)。gs と同じ値が要るときは `PYSAD_COMPAT=k64` を付ける。
動的口径の地図は EMIT のエミッタンスで決まる粒子の大きさを使うので、v2 の Score は 115 (gs・k64 は 118)、地図の端の数文字が変わる。

△ 概ね一致。Score (118) と、地図の 3 行は 1 文字 (または 2 文字) を除いて一致 (gs は `j`・`N`・`A`、pysad は `J`・`U` の箇所)。
損失の境界付近の粒子の周数が 1 周違うことによる差で、docs/08 の説明 (「損失境界の 1 粒子の周数が SAD と違う」、カオス的な境界) と同じ。走査結果の戻り値 (RES) も同じ理由で数カ所違う (数値は 1 周以内)。

## 未確認・限界
- 50 周、3 GeV、NORAD の短い追跡。放射減衰を入れる `Damping->True` は pysad 未対応。`DAPWIDTH` は 51 のときだけ gs と同じ窓にしてある (pysad は 51 以外の窓を再現しない)。
- 地図の `j`/`J` などの文字の意味 (どの周か) は SADHelp の DynamicApertureSurvey の記述 (損失した周数を 1 文字で示す) による。
