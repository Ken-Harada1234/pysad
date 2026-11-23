# 33 Fortran 別の翻訳状況: gs との照合の桁の目安 (旧 docs/16 §0 から移した)

docs/33_fortran_status.md の本文には、各行の「照合」欄の個別の桁を残している。ここには、旧 docs/16_fortran_status.md §0 にあった全体の目安の段落を元のまま移す。ベンチ (Level) の一覧は docs/27_benchmarks.md (旧 10_benchmarks)。

照合桁の目安 (10_benchmarks): CALC 光学 ≤1e-10 (誤差・オフモメンタム・結合経由は 1e-7〜3e-7)、Emittance[] ≤1.3e-11 (Level 3C の結合は 7.8e-7)、追跡 ≤5e-11 (Level 7/9/11 の z・trf0 は ≤1.7e-10)、マッチング K1 ≤1e-10、DA は地図・Score 一致。

旧 docs/20 の冒頭にあった試験の件数の行 (「試験は `py -3.12 -m pytest pysad/tests -q` (旧式の 3 本を除く) で 937 通過・4 skip (941 件収集、2026-11-17。GUI 系は `PYSAD_SAD_PACKAGES` が無いと skip。k64 モードが既定)、回帰 (`QUICK=1`) は 50 ベンチ FAIL 0 (2026-11-17)」) は、旧 16 §0 の「試験の現状」と同じ内容なので、docs/33 §0 の 1 つにまとめた。
