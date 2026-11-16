"""measure.py の tlum / tluma: ガウス分布の二つのビームの重なり面積の推定が解析値 4π·σx·σy·… に近いこと (参照の gs は MEASURE が動かないため取れない)"""
import numpy as np


def test_tlum_gaussian_overlap_is_close_to_analytic():
    from pysad.measure import tlum
    rng = np.random.default_rng(3)
    n = 40000
    sx, sy = 1.0e-3, 2.0e-4
    x = rng.normal(0.0, sx, (2, n))
    y = rng.normal(0.0, sy, (2, n))
    val, nmesh = tlum(x, y, 0.0, 0.0, sx, sy)
    assert nmesh == 200
    # 同じガウス分布 2 つ: Σ bin1·bin2 = n²Δ²/(4π σxσy)、Δ² = σxσy(6/nmesh)² なので tluma = σx·σy (実効面積)
    assert abs(val / (sx * sy) - 1.0) < 0.1
