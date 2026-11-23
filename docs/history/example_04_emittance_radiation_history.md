# example/04_emittance_radiation の md から取り除いた履歴の記述

例の md を「現在の動作と結果の説明」にするため、版・経過・以前の状態の記述を取り除いた。取り除いた元の行を、ファイルごと・元の行番号つきで残す (行番号は整理前のもの)。この節は、整理前の md (リポジトリの先頭の版) と現在の md の差分から再構成した。

## FFS_EMITTANCE.md

元の 66 行目:
> **v2 の注記 (2026-11-17、v2.0.0 の再実行)**: `.pysad.txt` は既定の v2 モードの出力に作り直した。v2 は `Emittance[]`・EMIT・追跡の Jx の端の項 (B75) と `MomentumCompaction` (B77) を直したので、この例の pysad の値のうち εx (約 2.4% 大きい)・εz・減衰率の X と Z・σδ・バンチ長・α・バケット高さ・放射によるチューンシフトが gs と違う (K64 のバグを直したため。`docs/47_v2_changes.md`)。以下の表と本文の「pysad」の値と「一致」の判定は、v1 と同じ `PYSAD_COMPAT=k64` での結果 (旧い出力と同じことを 2026-11-17 に確認)。gs と同じ値が要るときは `PYSAD_COMPAT=k64` を付ける。

元の 72 行目:
> - $J_x$ の RADINT 対 EMIT の 2.4% の差の原因は、端の項の取り扱いと推定した段階 (ソース未確認、他の BEND 形状での確認なし)。どちらが物理的に正しいかは、追跡で測っていない。**2026-11-17 に v2.0.0 で決着** (`docs/47_v2_changes.md` の B75): 解析式 Jx = 1 − I4/I2 と厳密な幾何 + Robinson の式が一致し、RADINT が正しく K64 の EMIT が誤り (端の角の蹴りが行路を変えない)。v2 の EMIT の Jx は 0.9429 (この例の再実行) で、RADINT の 0.9415 に近づいた。k64 モードは gs と同じ 0.9651。この例のリング以外 (矩形磁石でない BEND など) の確認は記録の条件に限る。

## FFS_RADINT.md

元の 51 行目:
> **v2 の注記 (2026-11-17、v2.0.0 の再実行)**: `.pysad.txt` は既定の v2 モードの出力に作り直した。v2 は `Emittance[]`・EMIT・追跡の Jx の端の項 (B75) と `MomentumCompaction` (B77) を直したので、この例の pysad の値のうち εx (約 2.4% 大きい)・εz・減衰率の X と Z・σδ・バンチ長・α・バケット高さ・放射によるチューンシフトが gs と違う (K64 のバグを直したため。`docs/47_v2_changes.md`)。以下の表と本文の「pysad」の値と「一致」の判定は、v1 と同じ `PYSAD_COMPAT=k64` での結果 (旧い出力と同じことを 2026-11-17 に確認)。gs と同じ値が要るときは `PYSAD_COMPAT=k64` を付ける。

## README.md

元の 10 行目:
> | flag_RADTAPER | RADTAPER の前後 | ○ (v2 で対応。k64 モードで gs と一致。v2 の既定では εx が 2.4% 違う: Jx の修正 B75) |

元の 15 行目:
> - RADINT の $J_x$ (0.9415) と EMIT の $J_x$ (0.9651) は 2.4% 違う。$I_4$ の端の項を 2 倍にすると EMIT に近づく (FFS_EMITTANCE.md)。**v2.0.0 で原因が分かり直した**: K64 の `Emittance[]`・EMIT が端の角の蹴りの行路を入れず $I_4$ の端の項が 2 倍になっていた (B75。RADINT が正しい)。v2 の既定では EMIT の $J_x$ は 0.9429 で RADINT に近い (gs と違う値になる)。`PYSAD_COMPAT=k64` では gs と同じ 0.9651。この節の他の例の pysad の出力 (`.pysad.txt`) は v2 のもの。

元の 16 行目:
> - 例 (`FFS_EMITTANCE`、`FFS_RADINT`、`fn_Emittance`、`flag_RFSW_RAD_RADCOD_FLUC`、`fn_SynchroBetaEmittance`) の判定は k64 モードでの gs との一致。v2 の既定では εx・減衰率・α などが B75・B77 の修正の分だけ gs と違う。

## flag_RADTAPER.md

元の 37 行目:
> **○ 一致 (k64 モード)。v2 の既定では εx・εz が gs と 2.4%・1% 違う (意図した差)。** 2026-11-17 の再実行 (pysad v2.0.0)。RADTAPER は v2 で実装した (最初の版では「フラグを保持するだけで磁石をスケールしない」未対応だった)。

元の 52 行目:
> - **εx・εz の差は RADTAPER ではなく v2 の修正による**: RADTAPER なしの εx も 1.7358E-7 (gs) と 1.7767E-7 (v2) で違う。K64 の `Emittance[]` の Jx の端の項 (B75) を v2 で直したため (`docs/47_v2_changes.md`)。k64 モードでは gs と 5 桁以上一致する。

## flag_RFSW_RAD_RADCOD_FLUC.md

元の 35 行目:
> **v2 の注記 (2026-11-17、v2.0.0 の再実行)**: `.pysad.txt` は既定の v2 モードの出力に作り直した。v2 は `Emittance[]`・EMIT・追跡の Jx の端の項 (B75) と `MomentumCompaction` (B77) を直したので、この例の pysad の値のうち εx (約 2.4% 大きい)・εz・減衰率の X と Z・σδ・バンチ長・α・バケット高さ・放射によるチューンシフトが gs と違う (K64 のバグを直したため。`docs/47_v2_changes.md`)。以下の表と本文の「pysad」の値と「一致」の判定は、v1 と同じ `PYSAD_COMPAT=k64` での結果 (旧い出力と同じことを 2026-11-17 に確認)。gs と同じ値が要るときは `PYSAD_COMPAT=k64` を付ける。

## fn_Emittance.md

元の 42 行目:
> **v2 の注記 (2026-11-17、v2.0.0 の再実行)**: `.pysad.txt` は既定の v2 モードの出力に作り直した。v2 は `Emittance[]`・EMIT・追跡の Jx の端の項 (B75) と `MomentumCompaction` (B77) を直したので、この例の pysad の値のうち εx (約 2.4% 大きい)・εz・減衰率の X と Z・σδ・バンチ長・α・バケット高さ・放射によるチューンシフトが gs と違う (K64 のバグを直したため。`docs/47_v2_changes.md`)。以下の表と本文の「pysad」の値と「一致」の判定は、v1 と同じ `PYSAD_COMPAT=k64` での結果 (旧い出力と同じことを 2026-11-17 に確認)。gs と同じ値が要るときは `PYSAD_COMPAT=k64` を付ける。

## fn_SynchroBetaEmittance.md

元の 31 行目:
> **v2 の注記 (2026-11-17、v2.0.0 の再実行)**: `.pysad.txt` は既定の v2 モードの出力に作り直した。v2 は `Emittance[]`・EMIT・追跡の Jx の端の項 (B75) と `MomentumCompaction` (B77) を直したので、この例の pysad の値のうち εx (約 2.4% 大きい)・εz・減衰率の X と Z・σδ・バンチ長・α・バケット高さ・放射によるチューンシフトが gs と違う (K64 のバグを直したため。`docs/47_v2_changes.md`)。以下の表と本文の「pysad」の値と「一致」の判定は、v1 と同じ `PYSAD_COMPAT=k64` での結果 (旧い出力と同じことを 2026-11-17 に確認)。gs と同じ値が要るときは `PYSAD_COMPAT=k64` を付ける。


# 台本 (`.sad`・`.py`) の冒頭コメントから、機能別フォルダへの移動時に取り除いた履歴の記述

v2_next・v2_C・v2_B の作業フォルダから機能別のフォルダへ移したとき、台本 (`.sad`・`.py`) の冒頭コメントにあった作業符号 (「v2next A1 項目 4」など) と、旧版の動作を語る文 (「旧 pysad は…」など) を取り除いた。取り除いた行を、元の行と直したあとの行の対で残す。

- `var_MASS_CHARGE_proton_ring.sad` 1 行目:
  - 元: `! v2next A1 項目 2: TrackParticles・Emittance[]・CALC6D の MASS・CHARGE (陽子)。旧 pysad は MASS を追跡に渡さず、電子として追った。`
  - 後: `! TrackParticles・Emittance[]・CALC6D の MASS・CHARGE (陽子)。`
- `flag_RADCOD_closure.sad` 1 行目:
  - 元: `! v2next A1 項目 7: RADTAPER・RADCOD の閉軌道が、追跡 (RAD) の 1 周で閉じるか。Emittance[] の閉軌道から始めた 1 周の動き (z は 0 から)。`
  - 後: `! RADTAPER・RADCOD の閉軌道が、追跡 (RAD) の 1 周で閉じるか。Emittance[] の閉軌道から始めた 1 周の動き (z は 0 から)。`
- `FFS_RADINT_I4V_coupled.sad` 1 行目:
  - 元: `! v2next A2 項目 3: 結合のあるリング (45 度傾いた四極 4 台、積分強度 0.05) の RADINT の I4V。`
  - 後: `! 結合のあるリング (45 度傾いた四極 4 台、積分強度 0.05) の RADINT の I4V。`
- `FFS_RADINT_pf.sad` 452 行目:
  - 元: `! ---- v2next A2 項目 3: PF リング (2.5 GeV、187 m) の RADINT・Emittance[]・EMITDIV と Polarization ----`
  - 後: `! ---- PF リング (2.5 GeV、187 m) の RADINT・Emittance[]・EMITDIV と Polarization ----`
- `elem_BEND_F1_emittance.sad` 1 行目:
  - 元: `! v2next B1 項目 1: Emittance[]・EMIT・CALC6D の BEND の F1・FRINGE。K64 の tbende の F1 (tblfre、tbendef1 の 2 断片、alcorr の長さの補正) と`
  - 後: `! Emittance[]・EMIT・CALC6D の BEND の F1・FRINGE。K64 の tbende の F1 (tblfre、tbendef1 の 2 断片、alcorr の長さの補正) と`
- `elem_BEND_F1_emittance.sad` 2 行目:
  - 元: `! 周長 (tfgeo の pos: BEND の F1 の長さ補正) を入れた。旧 pysad は F1 を黙って無視した (U0 が 5%、チューンが 0.006 違った)。`
  - 後: `! 周長 (tfgeo の pos: BEND の F1 の長さ補正) を入れた。`
- `FFS_RADINT_k1_bend_vdisp.sad` 1 行目:
  - 元: `! v2next C3 項目 5 (B107): K1≠0 の偏向磁石 + 結合による Y 方向の分散。RADINT の I4U + I4V (= 物理的な I4) と減衰分配 (DU, DV)。`
  - 後: `! (B107) K1≠0 の偏向磁石 + 結合による Y 方向の分散。RADINT の I4U + I4V (= 物理的な I4) と減衰分配 (DU, DV)。`
- `FFS_RADINT_quad_offset.sad` 1 行目:
  - 元: `! v2next C3 項目 5 (B107b): 据え付けのずれた四極の RADINT。軌道は四極の中心からずれ、四極が放射する (h = k x_o、I4 の項 2 k² x_o η)。`
  - 後: `! (B107b) 据え付けのずれた四極の RADINT。軌道は四極の中心からずれ、四極が放射する (h = k x_o、I4 の項 2 k² x_o η)。`
- `fn_SynchroBetaEmittance_coupled.sad` 1 行目:
  - 元: `! v2next C3 項目 3: SynchroBetaEmittance を結合のあるリングで走らせる。24 セルの FODO (3 GeV、RF 4 MV、6 極なし) に skew 四極 (K1 = 0.02、45 度) を 1 台。`
  - 後: `! SynchroBetaEmittance を結合のあるリングで走らせる。24 セルの FODO (3 GeV、RF 4 MV、6 極なし) に skew 四極 (K1 = 0.02、45 度) を 1 台。`
- `flag_NORFSW_orbit.sad` 1 行目:
  - 元: `! v2next D1 項目 1: NORFSW; RAD の Emittance[] と閉軌道のずれ (12 セル FODO、3 GeV、周長 45.6 m、QF に DX = 1e-4)。K64 の tcod は RFSW OFF で x px y py だけ閉じ δ・z は動かさない。`
  - 後: `! NORFSW; RAD の Emittance[] と閉軌道のずれ (12 セル FODO、3 GeV、周長 45.6 m、QF に DX = 1e-4)。K64 の tcod は RFSW OFF で x px y py だけ閉じ δ・z は動かさない。`
- `flag_EMITDIV_rect.py` 2 行目:
  - 元: `"""v2next P2 項目 2: v2 の Emittance[] と RADINT の差 (矩形磁石、SP8 (大きなリング B) で 1.4e-3) の原因 = EMITDIV の既定 (1) の離散化の誤差 (k64・v2・扇形・矩形とも同じ)。`
  - 後: `"""v2 の Emittance[] と RADINT の差 (矩形磁石、SP8 (大きなリング B) で 1.4e-3) の原因 = EMITDIV の既定 (1) の離散化の誤差 (k64・v2・扇形・矩形とも同じ)。`
- `flag_RADCOD_derivative.py` 2 行目:
  - 元: `"""v2next P2 項目 4: RADCOD ON の 1 周 (または要素) の解析の転送行列 T + D と、写像の数値微分 (4 次の中心差分、刻み 2e-5) の差。`
  - 後: `"""RADCOD ON の 1 周 (または要素) の解析の転送行列 T + D と、写像の数値微分 (4 次の中心差分、刻み 2e-5) の差。`
- `fn_SynchroBetaEmittance_dispersion_start.sad` 1 行目:
  - 元: `! v2next P2 項目 1: SynchroBetaEmittance (B128) はリングの始点の分散 η ≠ 0 のとき Emittance[] の数倍〜40 倍の値を返す (K64 の式の誤り)。`
  - 後: `! SynchroBetaEmittance (B128) はリングの始点の分散 η ≠ 0 のとき Emittance[] の数倍〜40 倍の値を返す (K64 の式の誤り)。`
- `fn_SynchroBetaEmittance_rings.py` 2 行目:
  - 元: `"""v2next P2 項目 1: SynchroBetaEmittance と Emittance[] の比 (AzimuthalModes = 3)。実行: py -3.12 example/04_emittance_radiation/fn_SynchroBetaEmittance_rings.py (出力は fn_SynchroBetaEmittance_rings.out.txt)`
  - 後: `"""SynchroBetaEmittance と Emittance[] の比 (AzimuthalModes = 3)。実行: py -3.12 example/04_emittance_radiation/fn_SynchroBetaEmittance_rings.py (出力は fn_SynchroBetaEmittance_rings.out.txt)`
- `flag_RAD_tradke_exact.py` 2 行目:
  - 元: `"""v2next P2 項目 3: 放射反作用の角度の項 ddpx の符号 (B111)。K64 の tradke (Emittance[]・閉軌道) は ddpx = +0.5 dpx dg、tradkn (追跡) は -0.5 dpx dg (dg < 0)。`
  - 後: `"""放射反作用の角度の項 ddpx の符号 (B111)。K64 の tradke (Emittance[]・閉軌道) は ddpx = +0.5 dpx dg、tradkn (追跡) は -0.5 dpx dg (dg < 0)。`
- `var_MASS_proton_beta_p1p5.sad` 1 行目:
  - 元: `! v2next P5 項目 8: 陽子 (MASS = 938.272088 MeV) の放射の β の因子 (B112)。Emittance[] の EnergyLossU0 を、古典式 U0 = (2/3) r_p m c² β³ γ⁴ Σ θ²/L と比べる。`
  - 後: `! 陽子 (MASS = 938.272088 MeV) の放射の β の因子 (B112)。Emittance[] の EnergyLossU0 を、古典式 U0 = (2/3) r_p m c² β³ γ⁴ Σ θ²/L と比べる。`
- `var_MASS_proton_beta_p3.sad` 1 行目:
  - 元: `! v2next P5 項目 8: 陽子 (MASS = 938.272088 MeV) の放射の β の因子 (B112)。Emittance[] の EnergyLossU0 を、古典式 U0 = (2/3) r_p m c² β³ γ⁴ Σ θ²/L と比べる。`
  - 後: `! 陽子 (MASS = 938.272088 MeV) の放射の β の因子 (B112)。Emittance[] の EnergyLossU0 を、古典式 U0 = (2/3) r_p m c² β³ γ⁴ Σ θ²/L と比べる。`
- `var_MASS_proton_beta_p30.sad` 1 行目:
  - 元: `! v2next P5 項目 8: 陽子 (MASS = 938.272088 MeV) の放射の β の因子 (B112)。Emittance[] の EnergyLossU0 を、古典式 U0 = (2/3) r_p m c² β³ γ⁴ Σ θ²/L と比べる。`
  - 後: `! 陽子 (MASS = 938.272088 MeV) の放射の β の因子 (B112)。Emittance[] の EnergyLossU0 を、古典式 U0 = (2/3) r_p m c² β³ γ⁴ Σ θ²/L と比べる。`
- `flag_RADTAPER_optics.sad` 先頭付近:
  - 元: `! C2d: RADTAPER`
  - 後: `! RADTAPER`

## 手で書き直した記述
- 「以前の記録の…」「以前の pysad は…」「記録は docs/history/v2next_*.md §…」「(.pysad_k64.txt)」の言及など、複数行にまたがる記述は個別に書き直した (`.sad` 側に残るのは現在の動作の説明だけ)。

## 削除した k64 モードの出力
- この移動で、旧 `v2_next` の `*.pysad_k64.txt` のうち、このフォルダに移る例の分を削除した (k64 モードの出力は `PYSAD_COMPAT=k64` で再現できる。`flag_PSPAC_k64_shape` だけは k64 オプションの例として残した)。
