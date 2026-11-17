# FFS_BEAMSIZE (BEAMSIZE: ビームサイズ)
種別: FFS コマンド。フォルダ: 02_optics_cod。ラティス: `00_lattice/FODO_ring.sad` (RF 電圧だけ 6 MV に変更 — 1.5 MV では放射損失 1.56 MV に足りず平衡が無い)。エネルギー 3 GeV。フラグ: CELL、RFSW、RAD、NOFLUC (EMIT のため)。特殊変数: `MINCOUP`=0.01。

## 何を見る例か
`BEAMSIZE` (略 `BEAM`) は、現在の Twiss パラメータと **EMITX EMITY SIGE** (特殊変数) からビームサイズを計算する (5 次元。シンクロトロン振動が大きいときは不正確で、EMIT の CODPLOT を使う)。サイズは `LINE["SIGX", 位置]` などで読め、`DISP BEAM` (略 `DISP B`) で表になる。
`BEAMSIZE` はエミッタンスに $\max(\epsilon_k,\ {\rm MINCOUP}(\epsilon_x+\epsilon_y))$ を使う (結合の下限)。平面リングでは $\epsilon_y=0$ なので `MINCOUP` で垂直サイズを決める。
EMITX は先に `EMIT` で平衡値に更新する (これも `FFS_BEAMSIZE.gs.txt` の `Emittance X = ...` の表)。

## スクリプト
```
USE RING; CELL; CALC;
RF VOLT 6.0E6;  CALC;
RFSW; RAD; NOFLUC;
EMIT;                       ! 平衡エミッタンス・エネルギー広がりを EMITX EMITY SIGE に入れる
MINCOUP = 0.01;
BEAMSIZE;
DISP BEAM 1 8;
Print["sigx at QF.1, B.1 (mm) = ", LINE["SIGX",{"QF.1","B.1"}]*1E3];
```

## 結果の読み方
`EMIT` の表 (放射損失 U0=1.5633 MV、Vc=6 MV、Emittance X = 1.736E-7 m、Energy spread = 1.192E-3 …)。`EMITX EMITY SIGE SIGZ` の行 = (1.736e-07 m, 2.846e-38 m, 0.001192, 0.012812 m)。
`sigx at QF.1, B.1` の行が x のビームサイズ [mm]、`sigy` が y。

## 物理・解析式
$$\sigma_x^2=\epsilon_x\beta_x+(\eta_x\sigma_\delta)^2,\qquad \sigma_y^2=\max\bigl(\epsilon_y,\ {\rm MINCOUP}(\epsilon_x+\epsilon_y)\bigr)\,\beta_y$$
($\epsilon_x$=1.736e-07 m、$\sigma_\delta$=0.001192、MINCOUP=0.01)。β, η は gs の `Twiss`。

| 量 | 式の値 | gs | pysad |
|---|---|---|---|
| σx(QF.1) [mm] | 1.28826 | 1.28836 | 1.28836 |
| σx(B.1) [mm] | 1.06693 | 1.067 | 1.067 |
| σy(QF.1) [mm] | 0.0577287 | 0.057729 | 0.057729 |
| σy(B.1) [mm] | 0.0662391 | 0.066239 | 0.066239 |

(β、η は QF.1: BX=6.42693、EX=0.618669; B.1: BX=4.69262、EX=0.477274。)

## gs と pysad の比較
**v2 の注記 (2026-11-17、v2.0.0 の再実行)**: `.pysad.txt` は既定の v2 モードの出力に作り直した。v2 は `Emittance[]`・EMIT・追跡の Jx の端の項 (B75) と `MomentumCompaction` (B77) を直したので、この例の pysad の値のうち εx (約 2.4% 大きい)・εz・減衰率の X と Z・σδ・バンチ長・α・バケット高さ・放射によるチューンシフトが gs と違う (K64 のバグを直したため。`docs/47_v2_changes.md`)。以下の表と本文の「pysad」の値と「一致」の判定は、v1 と同じ `PYSAD_COMPAT=k64` での結果 (旧い出力と同じことを 2026-11-17 に確認)。gs と同じ値が要るときは `PYSAD_COMPAT=k64` を付ける。
なお、保存してあった旧い `.pysad.txt` には `[pysad] DISP のモード B は未対応` と書いてあったが、これは「`DISP BEAM` は 2026-11-17 に対応」の前の出力で古かった。今回の再実行で `DISP BEAM` の表が出ている。

- サイズ・Twiss・`EMIT` の表の主要な値は一致 (EMITX の相対差 1e-6 以下)。表示の差: 垂直エミッタンスの雑音 (gs 2.8e-38、pysad 0 など)、EMIT の表の `alpha` の桁 (0.0222107 と 0.02221)。
- `DISP BEAM` (`DISP B`、EMIT 後) は 2026-11-17 に対応 (gs と同じ見出しと桁。射影 Twiss の列 AXp BXp EMITXp… を含む。雑音 1e-16 の列は pysad が 0)。`LINE["SIGX",...]` は pysad でも gs と同じ値を返す (ただし pysad は EMIT を先に実行していないと `BEAMSIZE` が `EMIT を先に実行してください` と止まる)。

## 未確認・限界
- 6 次元のビームサイズ (EMIT + CODPLOT の `SIZE...`)、結合のあるリングの回転角 (`Rot(deg)`)、結合のあるリングの `DISP BEAM` は未確認。
- RF 電圧を 6 MV に変えたことでエミッタンス自体は変わらない (放射積分で決まる。04_emittance_radiation の例で扱う)。
