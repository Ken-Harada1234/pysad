# 15_sadhelp_purpose_policy.md から取り除いた履歴の記述

15_sadhelp_purpose_policy.md (現在の仕様の解説) から、版・経過に結びついた記述を取り除いた。下は取り除く前の元の行を、節の見出しつきでそのまま並べたもの (行番号は整理前のもの)。

## # SADHelp (1.4.4.2k64 Updated: 10/19/2026) の項目の目的別分類と翻訳方針

整理前の 3 行目:

````text
最終更新: 2026-11-17 (現状のスナップショット。pysad の列は docs/14_pysad_vs_sadhelp.md の 2026-11-16 の照合による。RENUMBER・ElementValues の 2 行は 2026-11-17 に手で直した。MEASURE・Fork・OpenShared・Shared・SharedSize・System の行と「現状の補足」は 2026-11-17 に v2 の実装に直した)
````

## ## フック・特殊 (変数) (6) — 方針: 直接翻訳

整理前の 347 行目:

````text
| 変数 | `ElementValues` | Fortran tffsa.f | Couple.n (SetCoupledElements) を CALC の前に実行 (gs と一致)。VISIT の退避・BYE の復元は v2.1.0 で実装 (docs/04) | 直接翻訳 |
````

## ## 現状 (2026-11-17) の補足

整理前の 361 行目:

````text
- v2 モード (v2.0.0 から) で実装した: `Fork` (POSIX で本当に分ける)・`OpenShared`・`Shared`・`SharedSize` (fork した子と共有)、`System` (既定では実行せず、許可があれば実行)、`MessageName` の関数形、MEASURE・TRCOD ほか (K64 のソースの読み。gs が動かない B42)、`Polarization` (K64 が NaN)、PHSROT (K64 で動作しない)。k64 モードでは従来どおり (限定・未対応)。docs/20、47。
````

