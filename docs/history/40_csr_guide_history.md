# 40_csr_guide の履歴 (整理で本文から取り除いた記述)

本文 (`docs/40_csr_guide.md`) を版のないスナップショットにするとき、取り除いた版・経過の記述を、元の文脈つきで残したもの。元の行をそのまま引用する (引用した行のうち、履歴にあたるのは版・経過を述べた部分)。

## 元 3 行 (節: 40. CSR (コヒーレント放射) の解説)

> 最終更新: 2026-11-18 (v2.1.0 の後)

## 元 158 行 (節: 4. K64 の癖)

> ([`docs/21_sad_bug_report.md`](21_sad_bug_report.md) の B59・B60。B60 は v2.1.0 で原因を特定した。pysad の `CSRSetupOY` は両モードで正しい形の出力を返し、gs の wrongtype は再現しない。B59 は k64 モードが K64 どおりに再現し、v2 モードは直した。)

## 元 161 行 (節: 4. K64 の癖)

> - **B60**: `CSROYMatrix` が、`CSRSetupOY` の出力を渡すと `General::wrongtype` を返す (この K64 ビルドで再現)。**手で作った正しい形の実数リストなら gs でも動く** (以前の「手で作った実数リストでも」は誤り)。原因 (D2 の調査、v2.1.0): `tfcsroysetup` (csroy.f) が 4 番目 (zj) の実数リストへ `c_f_pointer(c_loc(rlist(kazjn)), zjn, [mphi,nj])` で書くとき、先頭の語 `rlist(kazjn)` (リストの頭 = head) から書き始める (要素は `rlist(kazjn+1)` から書くべきで、ほかの 4 つの配列は `rlist(k+1:…)` で正しい)。head が zjn(1,1) で上書きされるので、`tfcsroymat` の `tfreallistq(las%dbody(4), lazj)` が偽になる (B59 の `dpjn` と `zjn` が同じ記憶域も同じ箇所)。gs で `CSRSetupOY` の出力の各要素の `Head` を調べ、4 番目が実数 (1.2952814481588496e-5) で List でないことを確認した。登録名は `CSRSetupOY`・`CSROYMatrix` だが、ソース内の関数名は `tfcsroysetup`・`tfcsroymat` で綴りが違う。`OYStability` と `OYScan` は gs で完走しない (K64 の誤り)。pysad は両モードで正しい形の出力を返して動く (k64 で wrongtype を再現するかは要判断)。手で作った入力 (nr=64、nj=3、nphi=2) の 6x6 の行列 (36 要素) は gs と相対 4e-20 で一致 ([`example/v2_next/D2_csroymatrix.sad`](../example/v2_next/D2_csroymatrix.sad)、`test_csroymatrix_hand_made_input_matches_gs`)。

