# 18_radlight の検証の詳細 (gs との比較)

docs/18_radlight.md の「検証」の詳細。旧 docs/45_radlight_guide.md §5 の gs との一致の表と、旧 docs/36_spin_csr.md §4 の 2026-11-16 の確認を元の文脈つきで移した。

## 旧 45 §5 gs との一致の度合い


| 項目 | 試験 | 条件 | 一致 |
| --- | --- | --- | --- |
| 軌道の記録 (RADLIGHT) | [`pysad/tests/test_radlight.py`](../pysad/tests/test_radlight.py) (参照は gs の出力) | 1 GeV、BEND を含む 3 セル、2 粒子、複数周、始点 ls≠1 | 22 点 × 4 本 (t, x, y, z) × 2 粒子で相対 1e-12 以下。実際は差 0 |
| `RadiationSpectrum` | `test_radfield.py::test_radiation_spectrum_matches_gs` (参照は gs の出力) | 60 点のガウス × cos、λ 1〜5 mm | 相対 2e-6 以内の試験 (gs と全桁一致) |
| `RadiationField` | `test_radfield.py::test_radiation_field_of_uniform_motion` | 等速直線運動 (β=0.8)、観測点 {0.3, 0.1, 2.0} | gs 参照なし。Liénard-Wiechert の速度場と内部の点で相対 6e-3 (スプライン微分の離散化誤差) |


## 旧 36 §4 (RADLIGHT) 2026-11-16 夜の確認


2026-11-16 夜の確認:
- 軌道の「進行方向」は局所 z 軸で、`geo[0]=[[0,0,1,0],[-1,0,0,0],[0,-1,0,0]]` から局所 z 軸の大域成分は (1,0,0) なので、進行方向の距離は大域 x (軌道の 2 行目) に出る。
- RADLIGHT → RadiationField → FindPeak → RadiationSpectrum の流れ (B を ANGLE=0 K0=0.1 の STEER にして 1 周): pysad は最後まで通る (軌道 8 点 {4,8}、f {10,8}、ピーク位置 3、スペクトル {3,1001})。gs は RadiationField で SIGSEGV (B49) のため通しの比較はできない。軌道 (TrackParticles の戻り) の部分だけ gs と一致。
- 限界: 複数区間 (複数周) の軌道は区間の始めの点が重なる (t が同じ点が 2 つ) ため、pysad の `RadiationField` は `x must be strictly increasing sequence` で失敗する。gs も RadiationField が SIGSEGV で落ちるので比較できない。
未確認: 大きなリング・多数粒子での規模。解説は docs/18_radlight.md。
