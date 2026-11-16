# fn_MakeBump (bump.n の Bump・MakeBump)

種別: パッケージ (Packages/bump.n、correction.n から読まれる)・フォルダ 13_packages・ラティス errors_ring_pkg (誤差入りの 24 セル、3 GeV、周長 105.6 m)。pysad も動く (2026-11-17 に修正)

## 何を見る例か
閉バンプ (指定した位置の軌道だけを動かし、離れた所では軌道を元に戻す) を作る。`Bump[str, Condition->{{位置,"Y",値}}, Range->{始,終}, Plane->"Y"]` で条件を作り、`MakeBump[bump, str, opt]` でステアラの強さを解いて `LINE["K0"]` に足す。ここでは `BA5` の垂直軌道を 1 mm にする。

## スクリプト
```
y0=Twiss["DY","BA5"];
b=Bump[str,Condition->{"BA5","Y",1e-3-y0},Range->{"BA2","BA8"},Plane->"Y"];
x=MakeBump[b,str,opt];            (* {強さ, 使ったステアラの番号} *)
Print["R:ybump=",Twiss["DY","BA5"]];
```
(`Condition` の値は「今の軌道からの変化量」。`str` は水平と垂直の両方のステアラ `Steer[{SH1,SV1,...}]`。)

## 結果の読み方
`steerers used=` は使ったステアラの数、`R:ybump` は補正後の BA5 の $y$、`ybefore0`・`yafter0` と `ybefore`・`yafter` は範囲の外側の BPM (BA1、BA10) での $y$ の補正前と補正後。

## gs と pysad の比較
| 量 [m] | gs | pysad |
|---|---|---|
| BA5 の $y$ (補正前 → 後) | 3.0090e-5 → **1.0000141e-3** | 3.0090e-5 → **1.0000141e-3** |
| BA1 の $y$ (補正前 → 後) | -7.8194e-6 → -7.8075e-6 | -7.8194e-6 → -7.8075e-6 |
| BA10 の $y$ (補正前 → 後) | 9.6899e-6 → 9.7066e-6 | 9.6899e-6 → 9.7066e-6 |

- **○ 一致**。バンプの高さは目標の 1 mm に一致 (1.0000141e-3)。範囲の外側の BPM はバンプの約 1% しか動かず、ほぼ閉じている。使ったステアラは 12 本 (gs・pysad とも)。下位 1〜2 桁 (1e-17 程度) の雑音の差だけ。
- 以前の pysad は `???Set: list assignment length mismatch` で失敗した。`b=Bump[…]` は `bump.n` の `(a_=Bump[…])^:=…` (Set の UpValue) で、pysad は (1) `Bump` を自動読み込みの表に持たず、(2) 読み込んだ UpValue を組み込みの `Set` の定義として捨てていた。直した (docs/32、`work_example_fixes_20261117.md`)。
- gs の `Element values are not expanded` の注意は `CALC NOEXP` を内部で使うため (docs/32)。pysad も同じ文を出す。

## 未確認・限界
`Bump` の `Center`・`NumberOfSteers` の形 (3 バンプなど)、`CorrectOrbit` の `Bump->` オプション、水平のバンプは試していない。使ったステアラ 12 個の内訳 (水平・垂直) も見ていない。
