# 29_sol_optics から取り除いた履歴の記述

`docs/29_sol_optics.md` を現在のスナップショットだけの解説にするため取り除いた文章。元の文脈 (節の見出し) つきで、文章はそのまま。

## # SOL 領域の Twiss (CALC / INS) と Emittance[] / EMIT

取り除いた部分: 「最終更新: 2026-11-18 (v2.1.0 の後) (現状のスナップショット。」

元の行:

> 最終更新: 2026-11-18 (v2.1.0 の後) (現状のスナップショット。経緯は `history/` にある)

## ## 限界・未確認

取り除いた部分: 「は v2.1.0 で対応した (両モード。k64 は gs と 1e-9。docs/44 の CALC6D の節)**」

元の行:

> - **CALC6D・RADCOD つきの CALC、Touschek (INTRA) と SOL の組み合わせは v2.1.0 で対応した (両モード。k64 は gs と 1e-9。docs/44 の CALC6D の節)**。TouschekLifetime の値は SOL のあるリングで未確認 (gs が NaN)。BOUND=0 の SOL (BZ≠0) の CALC・RADTAPER + SOL・WSPAC + SOL の値は未確認。

## ## v2 モードで追加された範囲 (要約。詳細・根拠・限界は docs/44_sol_guide.md §7a)

取り除いた部分: 「| v2.1.0 で対応 (k64 も gs と 1e-9) |」

元の行:

> | CALC6D・RADCOD つきの CALC・Touschek (INTRA) | v2.1.0 で対応 (k64 も gs と 1e-9) | 対応 (両モード同じ) |

## ## v2 モードで追加された範囲 (要約。詳細・根拠・限界は docs/44_sol_guide.md §7a)

取り除いた部分: 「領域内の MULT・CAVI の POL (Emittance[]) は v2.1.0 で確認 (docs/44 の v2.1.0 の節。bz = 20 まで)。」

元の行:

> v2 の領域内の BEND・MULT・CAVI の条件と確認の範囲: 電子、0.5〜3 GeV、bz = 0.3〜3 1/m、ローレンツ力の数値積分と 1e-14〜1e-8。gs との数値比較は無い。bz ≥ 3 1/m を超える場、領域内の MULT・CAVI の POL (Emittance[]) は v2.1.0 で確認 (docs/44 の v2.1.0 の節。bz = 20 まで)。加速する CAVI の POL (Emittance[]) は未確認。領域内の四極の DX・DY (追跡・Emittance[]・幾何) は v2.1.0 で対応 (gs と全桁)。

## ## v2 モードで追加された範囲 (要約。詳細・根拠・限界は docs/44_sol_guide.md §7a)

取り除いた部分: 「(追跡・Emittance[]・幾何) は v2.1.0 で対応 (gs と全桁)。」

元の行:

> v2 の領域内の BEND・MULT・CAVI の条件と確認の範囲: 電子、0.5〜3 GeV、bz = 0.3〜3 1/m、ローレンツ力の数値積分と 1e-14〜1e-8。gs との数値比較は無い。bz ≥ 3 1/m を超える場、領域内の MULT・CAVI の POL (Emittance[]) は確認済み (docs/44。bz = 20 まで)。加速する CAVI の POL (Emittance[]) は未確認。領域内の四極の DX・DY (追跡・Emittance[]・幾何) は v2.1.0 で対応 (gs と全桁)。

## # SOL 領域の Twiss (CALC / INS) と Emittance[] / EMIT

取り除いた部分: 「(下の「v2 モードで追加された範囲」に要約)」

元の行:

> **モードによる違い**: この文書の本文 (作り方・検証・限界) は **k64 モード** (領域内は DRIFT・QUAD・L=0 の要素・SOL。gs と比較した範囲) の説明。**v2 モード** (既定) は、領域内の BEND・MULT・CAVI・K1=0 の四極、領域内の空洞の TRPT の加速・wake・AUTOPHI、ずれて傾いたソレノイド (`DX DY DZ CHI1..3`・`F1`) も Twiss・Emittance[]・EMIT・追跡に入る。これらは gs ではなくローレンツ力の数値積分などで確かめたもので、使える範囲・使えない範囲・根拠は docs/44_sol_guide.md §7a に書いた (下の「v2 モードで追加された範囲」に要約)。

## ## v2 モードで追加された範囲 (要約。詳細・根拠・限界は docs/44_sol_guide.md §7a)

取り除いた部分: 「## v2 モードで追加された範囲 (要約。」

元の行:

> ## v2 モードで追加された範囲 (要約。詳細・根拠・限界は docs/44_sol_guide.md §7a)
