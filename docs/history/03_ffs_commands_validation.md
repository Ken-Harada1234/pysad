# 03_ffs_commands の検証の記録 (K64 との比較の詳細)

docs/03_ffs_commands.md から移した、gs との比較の数値の元の記述。見出しは移す前の節。

## §2 BIPOL の行

````text
QF K1・QD K1 の 2 変数で NX・NY を 0.2 に fit した例で既定 K=−.0068645845/−.0057261776、UNIPOL K=.3289522375/−.3239879365 と gs と全桁一致 ([`test_item4.py`](../pysad/tests/test_item4.py))。BIPOL は tffsmatch.f だけで効く
````

## §2 ABSW の行

````text
`ffscmd.matcher_` で `match6` の `absweit` につないである。同じ fit 例で RELW の K=−.0024047842/−.0107636499 が gs と全桁一致 ([`test_item4.py`](../pysad/tests/test_item4.py))
````

## §3.3 成分ごとの FREE の検証

````text
小さな FODO (3 GeV、周長 18.4 m) で B.3 だけの K1 で NX を 1.30 にすると反復 3 回 (残差 2.1729E-03・4.8421E-08・3.2593E-17)、NX 8.168140905042474 が gs と全桁。大きなリング A (3 GeV、周長 603 m、3563 要素。`FREE BPI6.2 K1; FREE BPI7.5 K1` と `FREE QDA2.1 K1; FREE QDA3.4 K1; FREE QDA4`) でも反復の行・チューン (47.6000342385、17.6000441754)・K1 が gs と一致 (2.8 秒)。未確認:
````

## §5 EMIT のあとの USE の検証

````text
gs と全桁一致 (小リングで EMIT 前 δ = (−3.0771610e-4、6.82902209e-4、−1.023444e-2)、EMIT 後・USE 後 (−4.1421124e-5、9.50018237e-4、−9.97733e-3))。大きなリング A (CELL CALC・RAD・EMIT の 5 周追跡 → USE RING → 5 周追跡で 1.3e-11)・PF リングでも一致。
````

## §7 SynchroBetaEmittance の検証

````text
gs と emitx 相対 4e-7〜1e-5 (12 セルのリング + CAVI で確認)。2026-11-16 に大きなリング (24 セルの平面リング + CAVI VOLT=2e6、AzimuthalModes->7) の 4 点で emitx が gs と相対 1.4e-7〜1.4e-6 で一致 (emity は雑音: gs −2.5e-19 など、pysad 1.3e-20 など)。
````

## §7 TrackParticles の l_track の検証

````text
pysad は 11/16 に `TrackRun.l_track0` と `k_apert` の lcur==0 で再現し、gs (1 回目・2 回目) と全桁一致 (15 値、[`test_small_items.py`](../pysad/tests/test_small_items.py))。(元の K64 は B33 でこの状況に SIGSEGV するため、基準は gs を直した版。)
````

