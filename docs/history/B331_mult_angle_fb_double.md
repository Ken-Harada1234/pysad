# v2.5.0 記録: ANGLE ≠ 0 の MULT の FB1・FB2 による長さ L' の補正の二重適用 (B331)

日付: 2026-11-23。試験 `pysad/tests/script_tests/test_v250_mult_angle_fb_double.py` (12 件)。gs の保存した出力 `pysad/tests/script_tests/scripts/v250_mult_angle_fb_double__{bend,mult,multc}.{sad,gs.txt}`。

## 1. 症状

K64 は ANGLE ≠ 0 の MULT の FB1・FB2 (フリンジの長さ) による長さの補正 L' = L − ((φ FB1)² + (φ FB2)²)/(48 L)·sin(w/2)/sin(φ/2) (w = φ − ψ1 − ψ2) を、`tffs.f:1368` (tpara) で `p_L_MULT` に入れ、`tturn.f` がその値を `tmulta` の al に渡す。`tmulta.f:82-90` はスライス長 aln = al/ndiv から同じ式 (al はすでに L') をもう一度引く。BEND は `tffs.f:1282` の 1 回だけ (`tbend.f` に補正は無い)。

## 2. gs での確認

3 GeV/c、NORAD・NORFSW、L = 1 m、ANGLE = 0.3、E1 = 0.5、E2 = 0.3、FRINGE = 3、TrackParticles の 5 粒子 (基準、px = 1e-3、py = 1e-3、δ = 1e-2 を 1 つずつ)。

| 比較 | 差 (最大) |
|---|---|
| gs の BEND (F1 = 0.5) と gs の MULT (FB1 = FB2 = 0.5) | x 2.78e-7、px 3.4e-9、y 1.4e-9、py 7.9e-9、z 2.46e-8 |
| gs の BEND と gs の MULT (L = 1.00018814134845: 二重の補正の分だけ先に伸ばした) | 全成分で 0 (ビット単位) |
| FB = 0 の MULT と BEND | 0 |

L = 1.00018814134845 は、c(L) = ((φ FB1)² + (φ FB2)²)/(48 L)·sin(w/2)/sin(φ/2) として L − c(L) − c(L − c(L)) = 1 − c(1) (BEND の L') を解いた値 (c(1) = 1.882e-4 m)。二重でなければこの伸ばしで差は出ないので、差は二重の補正の分だけ。別の小さい例 (ANGLE = 0.1、F1 = 0.2、E1 = E2 = 0、δ = 1e-3) では x の差が 8.32e-10 で、二重の分の短縮 1.67e-5 m × (1 − cos φ)/φ·δ = 8.3e-10 の予想と一致。pysad の k64 は gs と 1e-15 以内 (3 つの台本)。

## 3. 修正 (v2 のみ)

`multa.k_multa` (numba) に引数 `b2` を足し、v2 では `aln` から 2 回目の補正を引かない (al は p_L として L' になっている)。`track6n` の呼び出しが `_cav_b2()` (= `compat.v2()`) を渡す。v2 の MULT (FB1 = FB2 = F1) と BEND は、E1・E2・FRINGE = 3・DISFRIN・DX・ROTATE で 1e-15 以内に一致。FB が効かない場合 (FB = 0・FRINGE = 0) と K2 のキックは変わらない。

## 4. 限界・未確認

- 確認は FRINGE = 3・K1 = 0・K2 以上なし・RAD OFF・3 GeV/c。
- FRINGE = 1・2 (片側だけ) の MULT と BEND の x の 3e-5 級の食い違いは、仕様の違い (BEND は FRINGE ≠ 0 で F1 を両端に使う。対にすれば差はこの項の分だけ) で、誤りではない。K1 ≠ 0 の y の食い違いは別の誤り (B340、`B340_mult_k1_fringe.md`)。
- CALC `tmultae` (`tmulta.f:327-333`) も同じ式を持つが、pysad は FB1/FB2 の CALC を未対応で、gs の CALC では確かめていない。
- L' の式そのもの (ヘルプの BEND の F1 の定義) の導出は未確認。「1 回が正しい」の根拠は、同じ磁場・同じフリンジの BEND (ヘルプの定義) との一致。
