# 16_spin_pol の検証の詳細 (gs との比較)

docs/16_spin_pol.md の「検証」の詳細。旧 docs/43_pol_guide.md §7 と旧 docs/36_spin_csr.md §1.3 を元の文脈つきで移した。

## 旧 43 §7 gs との一致の度合い


試した格子: エネルギー 3.0〜3.08 GeV (aγ = 6.8〜6.99)、縦ずれ DY 0〜2 mm、ROTATE 0〜0.03 rad。

| 量 | 平面リング | 縦ずれ・傾き・6 極のあるリング |
|----|----|----|
| PolarizationTime | 全桁 (123.8154010596169 s) | 13 桁 (130.35269523782 s) |
| SpinTune | 13 桁 | 12 桁 |
| PolarizationVector | 12 桁 | 10 桁以上 |
| SpinPrecessionMatrix・脱偏極ベクトルの表 | 一致 (書式まで) | 一致 |
| Emittances・DampingRate (POL 付き) | 11〜12 桁 | 11〜12 桁 |
| TrackParticles の戻り (7・8 列) | 4e-12 以下 | 4e-12 以下 (NOFLUC・FLUC) |

- MULT のリング、SOL (ドリフトと四極のある領域) のリングは偏極時間・スピンチューン 1e-11、偏極ベクトル 1e-10 (`atol`) で一致。追跡は平面リング、ずれ・回転のあるリング、MULT (K1 K2 K3)、K1 のある BEND、SEXT、薄い MULT、SOL 領域内のドリフトで 4e-12 以下 (試験の許容は `tol=1e-8` 相対)。
- すべて MOMENTUM 3 GeV 前後、周長は十数 m〜数十 m の小さな試験リング (周長に SOL の長さは入らない)。実機の格子での照合はしていない。
- 試験は 13 件が通る (`pytest pysad/tests/test_spin.py pysad/tests/test_spin_track.py`、2026-11-16 に確認、約 12 秒)。


## 旧 36 §1.3 gs との一致


| 量 | 平面リング | 縦ずれ・傾き・6 極のあるリング |
|----|----|----|
| PolarizationTime | 全桁 (123.8154010596169 s) | 13 桁 (130.35269523782 s) |
| SpinTune | 13 桁 | 12 桁 |
| PolarizationVector | 12 桁 | 10 桁以上 |
| SpinPrecessionMatrix・脱偏極ベクトルの表 | 一致 (書式まで: Fortran の 1p g15.7) | 一致 |
| Emittances・DampingRate (POL 付き) | 11〜12 桁 | 11〜12 桁 |

エネルギー 3.0〜3.08 GeV (aγ=6.8〜6.99)、縦ずれ DY 0〜2 mm、ROTATE 0〜0.03 rad の格子を試した。
