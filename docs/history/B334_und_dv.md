# v2.5.0 記録: UND の z の更新の −dv·L (B334)

日付: 2026-11-23。試験 `pysad/tests/script_tests/test_v250_und_dv.py` (9 件)。gs の保存した出力 `pysad/tests/script_tests/scripts/v250_und_dv__{und,drift}_{3gev,10mev}.{sad,gs.txt}`。

## 1. 症状

K64 `undulator.f:109` の z の更新は `z(i) = z(i) − 0.5*(fgx²+fgy²)*dpds²/ds` (横の運動による経路の長さ) だけで、速度の差による −dv·ds の項が無い。DRIFT (`tdrift`) は `z += dpz·al1 − dv·al` で持つ。K64 の UND は磁場を設定できない (B37) ので、磁場なしの UND は DRIFT と同じはずが、δ ≠ 0 の粒子の z が dv·L (dv = −δ(2+δ)/[h(h + p h0)]) だけ違う。R56 = L/γ² の欠落。

## 2. gs での確認

NORFSW・NORAD、UND (L = 1 m、KX = 0.1、QY = 0.2、POLE = 10、SLICE = 4、PHI = 30) と DRIFT (L = 1 m)、5 粒子 (x = 1e-4、x = 2e-4、px = 1e-4、px = 2e-4、すべて 0 の順で δ = 1e-3・1e-3・1e-2・−1e-2・0)。

| p | z(UND) − z(DRIFT) (δ = 1e-3, 1e-3, 1e-2, −1e-2, 0) |
|---|---|
| 3 GeV/c | −2.897e-11、−2.897e-11、−2.858e-10、+2.945e-10、0 |
| 10 MeV/c | −2.6005e-6、−2.6005e-6、−2.5659e-5、+2.6439e-5、0 |

自由飛行の式 L(h0/p0 − h/p) (基準粒子より早い δ > 0 で z > 0) の符号を変えたものと相対 2e-3 以内で一致 (横の運動の項は UND と DRIFT で同じ大きさで、差には残らない)。k64 の追跡は gs と全桁一致。

## 3. 修正 (v2 のみ)

`und.k_und`・`und.k_und_field` の最後 (dv を更新したあと) に `S[i,4] -= S[i,6]*(nslice*ds)` (引数 `b2`、呼び出しは `compat.v2()`)。放射つきの python 版 (`track6` の刻みごとの経路) は刻みごとに `c[4] -= p.dv*ds`。v2 の磁場なしの UND の z が DRIFT と 1e-14 以内、FBX = 1e-9 T の弱い磁場でも DRIFT の z と 1e-12 以内。

## 4. 限界・未確認

- gs で確認できたのは磁場なしの UND (K64 は磁場を設定できない)、3 GeV/c と 10 MeV/c、L = 1 m。v2 の磁場つき UND の z の絶対値は、磁場を 0 に近づけた極限しか確かめていない。
- 放射つき: numba 版は最後に全長 L に最終の δ の dv を掛ける (放射で δ が変わる分は 2 次の差)。python 版は刻みごと。
- UND の Twiss・Emittance[] (v2 の `und.emit`・`solmap.step` は追跡の写像の数値微分) は追跡と同じ写像なので R56 が入る (別に確かめていない)。K64 は UND を CALC で恒等として扱う (k64 は変えない)。
