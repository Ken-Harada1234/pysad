# v2.5.0 記録: K1 ≠ 0 の MULT の FB1・FB2 の端が途中にも入る (B340) と、MULT・BEND の FRINGE の対応

日付: 2026-11-24。試験 `pysad/tests/script_tests/test_v250_mult_k1_fringe.py` (9 件)。gs の保存した出力 `pysad/tests/script_tests/scripts/v250_mult_k1_fringe.{sad,gs.txt}` (1 回の実行に BEND・MULT 8 ケース)。

## 1. 前の作業 (B331) で残った 2 つの食い違いの切り分け

B331 の作業で、同じ磁場・同じ長さの MULT と BEND (ANGLE ≠ 0) に (a) FRINGE = 1・2 で x の 3e-5 級、(b) K1 ≠ 0 で y の 1e-5 級の食い違いが見つかっていた。gs で、3 GeV/c、L = 1 m、ANGLE = 0.3、NORAD・NORFSW、7 粒子 (x・px・y・py・δ の単独と複数の同時励起、δ = 1e-2 を含む) を追跡して調べた (台本は 1 回の gs の実行に複数の LINE をまとめた)。

### 1.1 (a) FRINGE = 1・2: 仕様の違い (誤りではない)

`tffs.f:1263-1277` (BEND): FRINGE ≤ 0 のとき fb1 = 0 (FRINGE = -1 以外)・fb2 = 0 (FRINGE = -2 以外)。FRINGE > 0 (1・2・3) は fb1・fb2 が両方残る (F1 を両端に使う。SADHelp の BEND の FRINGE の説明「両端」と同じ)。`tffs.f:1346-1364` (MULT): FRINGE = 1 は fb1 だけ、2 は fb2 だけ、3 は両方。したがって BEND の FRINGE = 1・2 と MULT の FRINGE = 1・2 は別の磁場。対にするのは BEND の -1・-2・3 と MULT の 1・2・3 (BEND の -1・-2 はヘルプにない)。

gs (F1 = 0.5、E1 = 0.5・0.3 の有無、DISFRIN = 0・1):

| 比較 | 差 (最大) |
|---|---|
| F1 = 0 の BEND と MULT (FRINGE = 0〜3、E1・E2 あり・なし) | 0 |
| BEND の FRINGE = 1 (F1 = 0.5) と MULT の FRINGE = 1 | 6.0e-5 (E1 = E2 = 0)、5.7e-5 (E1・E2 あり) |
| BEND の FRINGE = 2 と MULT の FRINGE = 2 | 7.5e-5、7.4e-5 |
| BEND の -1 と MULT の 1 | 1.26e-6 (E1 = E2 = 0)、2.5e-7 |
| BEND の -2 と MULT の 2 | 1.27e-6、2.6e-7 |
| BEND の 1・2・3 と MULT の 3 | 2.54e-6、5.1e-7 (E1・E2 あり) |

対にした差は B331 の二重補正の分だけで、v2 は差 0 (1e-15 以内)。pysad の k64 は gs の BEND・MULT と 1e-15 以内。

## 2. (b) K1 ≠ 0: K64 の誤り (B340)

### 2.1 症状

gs で ANGLE = 0.3、K1 = 0.3、E1 = 0.5、E2 = 0.3 の BEND と MULT:

| F1 (FB1 = FB2) | BEND と MULT の差 (最大) |
|---|---|
| 0 (FRINGE = 0 または 3) | 2.2e-11 (成分は z、他は 1e-13) |
| 0.5、両端 (BEND 3・MULT 3) | 1.41e-4 (y 7e-5、py 1.1e-4。E1 = E2 = 0 では 1.47e-4) |
| 0.5、入口のみ (BEND -1・MULT 1) | 7.0e-5 |
| 0.5、出口のみ (BEND -2・MULT 2) | 7.0e-5 |

DISFRIN の 0・1、E1・E2 の有無によらず同じ大きさ。入口のみ・出口のみの合計が両端 (1.4e-4) で、「余計な端が 1 つ」の大きさが 7e-5。

### 2.2 原因

`tmulta.f:148,236` は、MULT を「半スライス (入口) → 中間 (ndiv − 1 枚) → 半スライス (出口)」の `tbend` で作り、前半に (fb1, fb2, mfr = -1)、後半に (fb1, fb2, mfr = -2) を渡す。K1 = 0 の `tbend` (`tbend.f:400,415`) は `(mfring > 0 or -1) and fb1` のように mfr で端を選ぶので、前半は入口の fb1 だけ、後半は出口の fb2 だけが効き正しい。K1 ≠ 0 は `tbend.f:595` で `tbendi` に行き、`tbendi.f:73` (fb1 ≠ 0 なら入口) と `tbendi.f:159` (fb2 ≠ 0 なら出口) は mfr を見ない。前半の `tbendi` が fb2 を、後半の `tbendi` が fb1 を、それぞれ半スライスの反対側の端に入れてしまい、MULT の途中 (前半の終わり、後半の始め) に余計な端のキック (x の δ 依存・py の y と y³・z) が入る。

### 2.3 修正 (v2 のみ)

`multa.k_multa` が、前半の `_tbend` に fb2 の代わりに 0、後半に fb1 の代わりに 0 を渡す (`b2` = v2)。K1 = 0 (k_bend0) は mfr で選ぶので結果は変わらない。v2 の MULT は BEND と 1e-12 以内で一致 (E1・E2・DISFRIN・FRINGE 1・2・3 / BEND -1・-2・3)。k64 は gs と 1e-15 以内。確認は「修正後の v2 MULT が gs の BEND と同じ追跡になる」こと (独立した照合の基準は gs の BEND。v2 の BEND と gs の BEND の差 3e-9〜5.6e-7 は別の v2 の修正 (E1・E2 の端の幾何など))。

## 3. 限界・未確認

- 確認は ANGLE = 0.3・K1 = 0.3・NORAD・3 GeV/c。K2 以上を持つ MULT、RAD ON (`tbendi` の first・last による放射の端の断片 `f1r`・`f2r`) は未確認。RAD ON では pysad は `_tbend` が first = last = True で呼ぶので、v2 の修正で前半の f2r・後半の f1r の断片もなくなる (本来の MULT の中間に端の断片はない)。
- CALC の `tmultae` → `tbende` が同じ渡し方かは未確認 (pysad は FB1・FB2 の CALC が未対応)。
- BEND の FRINGE = -1・-2 はヘルプにない (ソースと gs で確認)。
