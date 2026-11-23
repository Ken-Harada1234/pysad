# 14_pysad_vs_sadhelp.md から取り除いた履歴の記述

14_pysad_vs_sadhelp.md (現在の仕様の解説) から、版・経過に結びついた記述を取り除いた。下は取り除く前の元の行を、節の見出しつきでそのまま並べたもの (行番号は整理前のもの)。

## # pysad と SADHelp (1.4.4.2k64 Updated: 10/19/2026) の対応表 — SAD 側の定義場所つき

整理前の 3 行目:

````text
最終更新: 2026-11-17 (現状のスナップショット。SADHelp は 1.4.4.2k64 Updated: 10/19/2026 のまま。pysad 側の列は 2026-11-16 の Interp.builtins と FFS の表から再照合した)
````

## ## 2. SADScript / FFS の関数 (94、pysad 実装 86)

整理前の 78 行目:

````text
数え方: 94 は `functions` 節の分類の中の関数の数 (CERN の SADHelp 1.4.0.0k64 の目次で 2026-11-17 に数え直した)。分類に入らない `DownhillSimplex` (pysad 実装) を加えると 95。FFS-dedicated-functions は関数 28 と下位項目 3 (見出しの 31)。
````

## ## 3. 特別変数 (48、pysad 参照あり 44)

整理前の 246 行目:

````text
| `ElementValues` | Fortran tffsa.f | Couple.n (SetCoupledElements) を CALC の前に実行 (gs と一致)。VISIT の退避・BYE の復元は v2.1.0 で実装 (docs/04) |
````

整理前の 265 行目:

````text
| `MOMENTUM` | Fortran ActTra.f,bb.f,bbstrhl.f,…; .n dapert.n | 参照あり (ffs.py, ffscmd.py, optics_fn.py)。初期値は K64 の initb1.f で 30 GeV (以前の文書の 1 GeV は誤り。docs/21) |
````

## ## 4. フラグ (60、pysad 登録 59)

整理前の 347 行目:

````text
| `WSPAC` | Fortran tffs.f (fname/sino) | 実装 (v2 モード: ビーム行列のガウス分布の力が `TrackParticles` に入る。Emittance[]・EMIT・CALC6D・CALC (RFSW ON) へも v2.1.0 で対応 (Emittance[]・EMIT は k64 も gs と全桁)。k64 モードの追跡は値を保持するだけ) |
````

## ## 6. 手で直した行 (機械照合を作り直したら戻す)

整理前の 371 行目:

````text
- §1: MATRIX・RECOVER・REVERSE・SHOW・STATUS・VARY (ffsextra.py)、ELSEIF (sadscript.py)、MEASURE・RENUMBER・SEED・READ の注 (2026-11-17 に MEASURE の行を v2 の実装に直した)。
````

整理前の 373 行目:

````text
  Graphics の行 (Plots.n などを動かす)、ListRandom・MemoryCheck・TracePrint (限定の代用)、Fork・OpenShared・Shared・SharedSize・System (2026-11-17 に v2 の実装に直した)、LOSSMAP・PSPAC・RADTAPER・SELFCOD・SPAC・WSPAC (同)、
````

## ## 7. 未確認・限界

整理前の 379 行目:

````text
- 「実装」は名前が登録されているという意味で、挙動が K64 と一致するかは要素ごとの文書による。Interp.builtins の登録名は 2026-11-16 に実行時に集めた (462 個。11/06 は 281 個)。
````

整理前の 383 行目:

````text
- 要素のキーワード表 (§5) は 2026-11-16 の機械照合のまま。v2 で足したキーワード (`UND` の `FBX`・`FBY`、`PHSROT` の `PX PY PZ`・`Dij`・`Bij`、`SOL` の `DZ`・`CHI1..3`・`AUTOPHI`、`CAVI`・`MULT` の `AUTOPHI`、`MARK` の `JDX JDPX JDY JDPY` の読み方など) は反映していない (要確認)。
````

