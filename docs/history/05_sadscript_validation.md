# 05_sadscript の検証の記録 (gs との突き合わせ。旧 04_sadscript.md から移したもの)

本文 (docs/05_sadscript.md) では 1〜3 行にした gs (K64) との比較の詳細を、元の文脈つきで残す。

## 旧 §5.1 (数値・行列・乱数)
- Spline: 「境界条件は確認済み: 自然スプライン `Spline[{{x,y},..}]` と `Spline[..,0,0]` の形・`Sort[list, Greater]`・`Scan` が gs と全桁一致。`Spline[data,0,0]` は gs も pysad も未評価のまま。」
- 複素行列の演算: 「`LDSolve` の終点 (cos 1, −sin 1) と 3 通りの呼び出しが gs と一致 ([`test_evalfixes.py`](../../pysad/tests/test_evalfixes.py) 2 件)。`Packages/Optics.n` の `TransferMatrix`・`BeamMatrix`・`NormalCoordinates`・`OrbitGeo` は gs と表示桁で一致 (gs は 1e-16 の雑音の要素を出し pysad は 0)。`PeriodicR` は呼び出し方が不明で比較できなかった (gs も引数エラー)。」

## 旧 §6 (特殊な機能)
- 文脈: 「K64 (init.n の定義と gs の実測 15 行の調査) の規則を写した。」「`$Context`・`$ContextPath` への代入、`WithContext`、入れ子の文脈 (両モード。gs で実測し、台本 `example/v2_next/B2_context.sad` の 19 行が gs と全行一致)」「この 51 個の名前 (gs の `Context[名前]` の実測)」「gs の `Names` が列挙しないので 51 個しか調べていない」「SAD の Packages の読み込み … は変えていない (`PYSAD_SAD_PACKAGES` つきで GUI・描画・パッケージ系の試験が通る)。」
- MessageName: 「例 `B2_symbol_message.sad` (12 行が gs と全行一致)。」「`General::mytag` の継承は gs の `h::mytag` の実測と同じ」
- Parallelize: 「(両モード。… `Packages/shared.n` の `Fork`・共有メモリ・`Wait4` の動作を gs で確認)」「WSL (Python 3.14) で、A (5 点)・B (リストの結果)・C (1/0 → Infinity)・E (x は 0 のまま)・F (子の Print) が gs と同じ。FODO の CALC を 8 通りの K1 で逐次と並列の結果がビット単位で同じ。4 つの `Sleep[1]` が 1.02 秒。」(後半 2 文は本文にも残した)
- VISIT/BYE: 「台本 `B2_visit_bye.sad` の 8 行が gs と全行一致。」`BYE ALL` は gs が SIGSEGV で比べられない。
- Packages の関数: 「小さな台本を gs と行ごとに突き合わせた (docs/34)。」

## 旧 §7 (テスト方法と確認)
SAD のベンチでの確認 (docs/history/00_RESUME.md、<work-olddocs>/sad_level3_4_benchmark_results.md Level 14 節):

- `bench_level13a.sad` / `bench_level13b.sad` (FIT/FREE/GO + Print/LINE/Twiss/MatchingResidual) をそのまま実行し、反復行・`Matched` 行・条件表・Print の出力が SAD (`*.result`) と同じ (K1 ≤1e-10、残差 5.0763e-16 vs 5.0759e-16)。
- `bench_level14a.sad` (= 20 セルのラティス + 色収差スキャンの台本: Get, OpenWrite, Table, Do, DP0 スキャン + `FFS["cell calc;"]`, PolynomialFit) で 21 点のチューンが ≤2e-8、PolynomialFit の係数は SAD の式どおり (相対差 ≤1e-5)。
- Element / BeamLine 系は bench_level24f (164 タグ一致)。
- 試験の件数の行: 「pytest 941 件を収集、2026-11-17 に 937 通過・4 skip。k64 モードが既定」。
