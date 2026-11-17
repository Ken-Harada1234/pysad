# fn_CorrectOrbit_Monitor_Steer (correction.n の CorrectOrbit・Monitor・Steer)

種別: パッケージ (Packages/correction.n)・フォルダ 13_packages・ラティス errors_ring_pkg (誤差入りの 24 セル、3 GeV、周長 105.6 m)

## 何を見る例か
`Monitor[{BPM の名前}]`、`Steer[{ステアラの名前}]` で補正系を定義し、`CorrectOrbit[kind, mon, str, opt]` で閉軌道のずれ (四極の据え付け誤差による) を最小二乗で補正する。ずれは BPM 48 台の rms で見る。

## スクリプト
```
mon=Monitor[Flatten[Table[{"BA"//ToString[k],"BB"//ToString[k]},{k,1,24}]]];
sh=Steer[Table["SH"//ToString[k],{k,1,24}]];   sv=Steer[Table["SV"//ToString[k],{k,1,24}]];
opt=CalculateOptics[1,LINE["LENGTH"],Twiss["*",1],FLAG["CELL"],2];
r=CorrectOrbit["Y",mon,sv,opt];     (* 補正量 r (24 個) を返し、LINE["K0"] にも足す *)
CALC;                               (* 軌道を計算し直す。CALC NOEXP では補正が効かない *)
```
`Twiss["DY",名前]` の BPM での値から rms を作る。

## 結果の読み方
`R:rmx0`・`R:rmy0` は補正前の水平・垂直の rms、`rmy1` は垂直補正後、`rmx1` は水平補正後、`rmy2` は水平補正後の垂直の rms。`nsteer=24`、`max kick` は最大のキック [rad]。gs は `Element values are not expanded` という注意を出すことがある。pysad も 2026-11-17 の再実行 (v2.0.0) で同じ注意と `NOCALEXP` の表示を出し、gs と同じになった (`NOKEEPEXP`・EXPAND の gs に合わせた対応による。最初の版の `.pysad.txt` は出していなかった。補正の数値は変わらない)。

## gs と pysad の比較
| 量 | gs | pysad |
|---|---|---|
| 水平 rms 補正前 → 後 [m] | 4.184880e-5 → 2.282492e-5 | 4.184880e-5 → 2.282492e-5 |
| 垂直 rms 補正前 → 後 [m] | 2.790881e-6 → 4.510214e-7 | 2.790881e-6 → 4.510214e-7 |
| 最大キック [rad] (垂直) | 2.212038e-5 | 2.212038e-5 |

- 全値が 9 桁以上一致 (○)。docs/34 の co5_corr の確認と同じ傾向。
- 垂直は rms が 6 分の 1 に下がる。**水平は 4.2e-5 → 2.3e-5 で、半分程度にしか下がらない**。BPM 48 に対しステアラ 24 の最小二乗なので完全には取れないはずだが、垂直ほど下がらない理由 (水平のステアラの位置と BPM の位相・分散の寄与) は調べていない。
- `CALC NOEXP` で補正が効かないことの確認は、docs/34 の記載 (gs・pysad とも効かない) による。この例では `CALC` を使った。

## 未確認・限界
`Solver->{"Micado",n}`・`CorrectLocal`・`CorStat` は試していない (試した Micado 6 本は垂直の rms が逆に悪化した: 補正後 2.5e-5。ステアラの選び方か設定の問題かは未確認で、例から外した)。2 平面を同時に補正する `{"X","Y"}` も試していない。
