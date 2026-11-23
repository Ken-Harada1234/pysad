# 01_overview の検証の記述 (gs との比較。本文から移したもの)

本文 (docs/01_overview.md) から、gs (K64) との比較の細部を移した。元の文脈つき。

## 1. GUI・描画の gs との確認 (旧 01 §1 の「対象外」の項)

11/16 夜の確認: `ListPlot`・`Plot` のオプション 10 例で戻りの Graphics 式は gs と文字列で完全一致 (見た目ではなく式の比較。`Options[g]` は gs が全オプションを返し pysad は `{}`)。**差あり**: 最初の `Window[]` の親が、gs では Toplevel、pysad ではルート窓。`OpticsPlot` の text 項目数が gs と違う (FODO 24 セルで gs 69、pysad 74。Tk のフォントによる `RemoveOverlap` の差)。`w=Window[]`・`e=Entry[w]` は `Widget[…]` 値になり、`e[Insert]`・`e[Get]`・`e[Delete]` は gs と同じく働く (docs/22・41)。

## 2. モジュール表の行にあった gs との一致の記述 (旧 01 §3)

- 行 `synchrobeta.py`: | synchrobeta.py | SynchroBetaEmittance1 (シンクロ・ベータ結合のエミッタンス。運動量のずれごとの閉軌道と tturne を取り方位角モード展開を解く。gs と emitx 相対 4e-7〜1e-5) | tfemit.f tfemits, temits.f ほか |
- 行 `coord.py / tcav.py`: | coord.py / tcav.py / und.py | COORD (追跡・6 次元写像・幾何、docs/25)、TCAVI (追跡 RAD FLUC まで、Emittance[] は L=0 のみ gs と一致、docs/26)、UND (K64 は磁場を設定できず磁場なしのドリフト。v2 モードは `FBX`・`FBY` の磁場つきで追跡・Twiss・Emittance[] に入る、docs/27) | tcoord.f, tcoorde.f, qcoord.f, ttcav.f, ttcave.f, undulator.f |
- 行 `spch_k64.py`: | spch_k64.py | PSPAC の K64 の形状関数 (k64 モードだけ。CIC・「ほぼ最近接」・5 点差分 + 双 3 次補間・3 次スプライン。NPARA = 1 で gs と rms が 12 桁以上) | spch.c・psn.c・tspch_.c |
- 行 `radfield.py`: | radfield.py | RadiationField・RadiationSpectrum。RADLIGHT の軌道記録 (`TrackParticles` が `{結果, 軌道{t,x,y,z}}` を返す) は `ffscmd._RadLight` (k64 モードは gs と差 0、v2 モードは 3 次元の軌道で BEND の中も記録。docs/45) | tlinit.f, tlstore, tlresult |
- 行 `csr.py`: | csr.py | CSRTrack・CSRHaissin・CSRSetupOY・CSROYMatrix・CSRMatrix (mode 1000 = csrhelmat を含む)・CSRConvert (= tfcsrvec)・CSRInit と実FFT・spline1 (docs/36。`Packages/CSR.n` の MakeZL などは gs と一致、ScanTrack・OYStability は未照合) | csrtrack.f, csroy.f, csrtest.f, tcftr.f, spline.f |

## 3. 旧 01 §6「到達点と未対応」の「SAD と一致しているもの」(未対応の表は docs/31 へ)

## 6. 到達点と未対応

SAD と一致しているもの (詳細は 10_benchmarks.md):
- 線形光学 (CELL/INS、CALC6D、ROTATE/ずれ付きの全要素、ANGLE=0 BEND の K0/K1、TRPT の加速と MARK OFFSET、幾何 GX..、RADINT) ≤1e-10 (誤差・オフモメンタム付きは ≤1e-7)。
- Emittance[]/EMIT (RFSW × RADCOD の 4 通り、INTRA、不安定な光学、報告の全行) ≤1e-11〜1e-9、Touschek ≤7e-12。
- マッチング (解析微分、DESCEND、fuzz、tsvd、2 段目フィット、INS、フック) は反復列・残差・K1 を SAD と同じ桁で再現。
- 追跡 (NORAD/RAD/FLUC、粒子損失) ≤5e-11、DA の地図と Score、誤差の実現 (DELX.. の乱数列)。
- SADHelp の対応 (14_pysad_vs_sadhelp.md の 11/16 時点): FFS コマンド 56 中 45 実装 (未対応 2、無視 9)、関数 94 中 46 (この数は表が数える SADHelp 記載の関数の数で、特殊関数・Fit・Packages の関数 (docs/33, 32) は SADHelp の外に増えている)、特別変数 48 中 44 を参照、フラグ 60 中 59 を登録。
- 要素・機能 (各文書の「限界」を必ず見る): COORD・TCAVI・UND・MULT ANGLE≠0・SOL の Twiss/Emittance[]/追跡・MAP・DELCOR・GaussianCoulomb・BBBrem・BEAMBEAM・WAKE・CSR・POL のスピン (Emittance[]・EMIT・TrackParticles)・RADLIGHT・SynchroBetaEmittance・特殊関数・Fit/FindRoot・SAD パッケージ (CorrectOrbit・MakeBump ほか)・FFS の小さなコマンド・描画 (docs/22〜35)。gs との一致の桁は各文書の表にある。
