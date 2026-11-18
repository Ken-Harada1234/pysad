# example/v2_next: 担当 B3 の例 (v2next、2026-11-18)

記録は `docs/history/v2next_B3.md`。実行は `py -3.12 run_sad.py B3_xxx.sad` (既定は v2。`PYSAD_COMPAT=k64` で K64 の動作)。gs の出力は `example/tools/run_example.py` で取った (`.gs.txt`)。

| 台本 | 項目 | pysad v2 (`.pysad.txt`)・k64 (`.pysad_k64.txt`) と gs (`.gs.txt`) の関係 |
|---|---|---|
| `B3_pol_sideband.sad` | 項目 1: 4.2 GeV の 12 セル FODO の Polarization | **gs・k64 は P = 1.27e-4 (spdepol が脱偏極率を 1000 倍大きく出す)、v2 は 0.119 (側帯の項込み。線形の式だけなら 0.205)**。偏極時間 gs 23.06 s、v2 23.08 s (緩和率の定義の違い) |
| `B3_spin_mc.py` | 項目 1: 線形の軌道 + 厳密なスピン回転の MC | 側帯の式の検証。`py -3.12 example/v2_next/B3_spin_mc.py 4.2 12 3000 30000 11 out.json` (3000 粒子 × 30000 周 約 8 分) |
| `B3_radpol_roundoff.py` (`.pysad.txt`) | 項目 2: ASTRID2 の RADPOL の緩和が 1.58 倍速い原因 | 出力: P0 = 1、0.5、0.25 から始めた dP/周と式の比 (1.578、0.934、0.979) = 倍精度の丸め |
| `B3_photons_bend.sad` | 項目 3: BEND・QUAD の光子 (PHOTONS) | k64 は gs と光子 358 個が全列一致 (E は 3e-8)。v2 は偏光 xi1 が ±1 以内 (K64 は −3735 などの誤り) |
| `B3_geofit_quad_dx.sad` | 項目 4: SOL 領域内の四極の DX・DY を変数にした幾何の fit | gs と v2 が 2 回の Newton で収束し、QF DX = 2.788484326e-3、QD DY = 2.667980563e-3 (1e-12 で一致)。k64 は SOL 自体のずれが未対応で動かない |
| `B3_naff_tunes.py` (`.out.txt`) | 項目 5: NAFF でチューンを取り出す | 4 つのリング (12 セル FODO、大きなリング A、PF、ASTRID2) の追跡のチューンが CALC と 1e-9 以内。wtunem (FFT) は 1e-4 |

`B3_photons_bend.gs.txt` は gs の出力のうち光子の行 (`PA`・`PB`) を試験 `test_v2_B3_photons.py` が読む。
