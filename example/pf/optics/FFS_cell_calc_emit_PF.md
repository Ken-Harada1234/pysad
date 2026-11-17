# FFS cell / calc / emit (PF): 実行の仕方と最初の三つのコマンド
種別: FFS コマンド (KEK の SAD チュートリアルの話題をもとにした)。加速器: PF 2.5 GeV リング (`lattice/v9_16_asym.sad`、周長 187.035 m、チューン 9.5975 / 5.2790)。フラグ: `ON LOG; ON RAD FLUC COD RFSW RADCOD;`

## 何を見る例か
ラティス台本の後ろに計算台本を付けて gs に渡す形 (`.sad` 一つ)、`USE`・`cell`・`calc`・`emit`・`SAVE ALL` の意味、周長・チューン・エミッタンスを読み出す。
gs への渡し方は `gs file.sad < /dev/null` (チュートリアルは `cat` でラティスと計算台本を 1 ファイルにして渡す)。トップレベル (要素・ライン定義) と FFS (`FFS;` 以降) の二段で、`stop` 1 回目で FFS が、2 回目で全体が終わる (ここでは `Exit;`)。

## スクリプト
```
ON LOG;                              ! 入力の復唱をしない
ON RAD FLUC COD RFSW RADCOD;         ! 放射・量子励起・COD・RF・放射による軌道
FFS;
$FORM = "15.10"; PageWidth = 1999;   ! 数値の桁と、行の折り返し無し
USE RING;  cell calc emit;           ! cell = 周期解、calc = 軌道と光学、emit = 平衡パラメータ
SAVE ALL;
Print["circumference [m] = ", LINE["S","$$$"]];  ...
```
- `cell` は周期解を求めるモードのフラグ (外すと `INS` = MARK の値から転送するモード)。`calc` は閉軌道のあと周期ツイス、`emit` は 6 次元 (RF 入り) の平衡量。
- コマンド (`cell` など) は大文字小文字を区別しない。関数 (`Print[]`、`Twiss[]`) は区別し、`[]` が付く。
- `$$$` はライン終点、`^^^` は始点。

## 結果の読み方
| 量 | 値 | 出所 |
|---|---|---|
| 周長 | 187.0354880 m | `LINE["S","$$$"]` |
| 要素数 | 756 | `LINE["LENGTH"]` |
| チューン | 9.5975304909 / 5.2789602116 | `Twiss["NX","$$$"]/2/Pi` |
| εx | 3.4538e-8 m rad | `EMITX` (`emit` の表: Emittance X) |
| エネルギー拡がり | 7.3366e-4 | `SIGE` |
| バンチ長 | 9.765 mm | emit 表 |
| U0 | 0.3991 MV | emit 表 |
| 運動量コンパクション率 | 6.563e-3 | emit 表 |
| 有効ハーモニック数 | 311.94 | emit 表 (RF の FREQ から求めた値で整数ではない。`v9_16_inj_noSF` は `HARM = 312`) |
| 減衰時間 | x 7.58 ms、y 7.82 ms、z 3.97 ms | emit 表 |
| 減衰分配 | Jx 1.0313、Jy 0.9999、Jz 1.9688 | emit 表 |

## 物理・解析式
放射積分からの εx・U0・α・Jx の照合は `pf/tune_radiation/calc_RadInt_PF` (式と数値の表)。

## gs と pysad の比較
**v2 の注記 (2026-11-17、v2.0.0 の再実行)**: `.pysad.txt` は既定の v2 モードの出力に作り直した。保存してあった旧い出力は `PYSAD_COMPAT=k64` の出力と同じ (確認済み)。v2 は K64 のバグを直したので (`docs/47_v2_changes.md`。主に `Emittance[]`・EMIT・追跡の Jx の端の項 B75 と `MomentumCompaction` B77)、この例の出力のうち `Emittance[]` の表示 (εx が数 % 大きい、εz・減衰率・α・バケット高さなど) と、それに依る量が gs と違う。以下の「一致」の判定は k64 モードの結果。gs と同じ値が要るときは `PYSAD_COMPAT=k64` を付ける。

○ 全項目が表示桁内で一致。差は εy の浮動小数点雑音だけ (gs 7e-38、pysad -3e-43 で、どちらも 0)。

## 未確認・限界
- 偏向磁石 28 台・誤差なしの設計ラティス。入射点 `INJ` や ID (`U*`) は MARK で、K 値は 0。
- 垂直エミッタンスは結合なしで 0。
