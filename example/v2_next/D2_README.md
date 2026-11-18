# v2next D2 の例 (パーサ・GUI・パッケージ・整理)

記録: `docs/history/v2next_D2.md`。gs は SAD K64 (1.4.4.2k64、2026-10-19 ビルド、WSL)。`.gs.txt` が gs の出力 (台本の echo つき)、`.pysad.txt` が pysad v2、`.pysad_k64.txt` が k64 モード。

| 例 | 内容 | gs と pysad |
|---|---|---|
| `D2_main_expr` | MAIN の要素定義の式 (`ANGLE=PI/8`、`SIN(0.1)`、`2*PI/16`、`45 DEG`、前に代入した変数と関数、`1/0` は Infinity) | 13 個の ANGLE ほか全値が gs と一致 (v2・k64 とも同じ) |
| `D2_periodicr` | Optics.n の `PeriodicR` (4x4 の一周行列から R・A・B) | 6 行列 x 12 成分が 3e-16 以内 (NaN になる 2 行列も同じ) |
| `D2_correctlocal` | `CorrectLocal` (局所バンプ。領域 {BA3, BA8} の外の軌道は変わらない) | 24 BPM の DY が 1e-14 以内、補正量 12 個も同じ |
| `D2_micado` | `CorrectOrbit` の Micado (`Solver->{"Micado",n}`) の補正後 rms | k64 は gs と 3e-8 以内 (どちらも n に対して単調に下がらない)。v2 は n が増えるほど単調に下がる (K64 の pmicad の誤り) |
| `D2_csroymatrix` | `CSROYMatrix` を手で作った実数リストで | 36 要素が 4e-20 (相対) で一致。`CSRSetupOY` の出力の 4 番目の Head が gs では数 (B60 の原因) |

リングは 12 セルの FODO (3 GeV、周長 52.8 m)。四極の据え付け誤差 DX・DY は式 `1E-4*SIN(k*1.7)` で与える。
