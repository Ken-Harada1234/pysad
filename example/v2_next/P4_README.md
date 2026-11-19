# example/v2_next: 担当 P4 の例 (SOL・空間電荷・DAPWIDTH・FREE の残り、2026-11-19)

記録: `docs/history/v2next_P4.md`。試験: `pysad/tests/test_v2_P4_*.py`。gs の値は K64 (2026-10-19 ビルド) の実測。**gs の追跡の例は `NPARA = 1`** (gs の既定は粒子を子プロセスに分け、PSPAC は各プロセスが自分の粒子だけで電位を作る)。
実行: `py -3.12 run_sad.py <台本>.sad` (v2 が既定。`PYSAD_COMPAT=k64` で K64 の動作)。gs は `gs <台本>.sad` (WSL)。

| 例 | 内容 | gs との一致 | ファイル |
|---|---|---|---|
| `P4_free_noop_variable` | 成分ごとの FREE の何もしない変数 (`FREE QF.3 K1`)。VAR の表・エラー・フィットの変数の数 | 全桁 (両モード) | `.sad`・`.gs.txt`・`.pysad.txt` (v2)・`.pysad_k64.txt` |
| `P4_dapwidth_window` | `DynamicApertureSurvey` の DAPWIDTH の窓 (7・3・12・51)。RAD + FLUC の最終 seed (DAPWIDTH = 2) | NORAD の地図・Score が一致 (境界の 1 点 1 周差)。**seed は k64 で 4/5 が一致** (v2 は放射の式が違うので合わない) | 同上 |
| `P4_dapwidth_ringA` | 大きなリング A (`kek-ls_v03_12.lat`、3 GeV、603 m) の DAPWIDTH の窓 (7 と 51): 地図 (gs と pysad) | 地図・Score 63 が一致 (境界の 1 点 1〜2 周差) | `.out.txt` |
| `P4_pspac_k64_shape` | k64 の PSPAC (K64 の形状関数)。5 MeV・1e11・縦 4 セルの粗いメッシュ | **k64 が全桁一致**。v2 は rms x が 4.5% 違う | 同上 |
| `P4_pspac_astrid2` | 大きなリングの確認: ASTRID2 (0.58 GeV、45.7 m) に SPCH 3 か所、PBUNCH = 1e14、10 周 | **k64 が 11 桁以上一致** (v2 は x が 8% 違う) | 同上 |
| `P4_wspac_tracking_beam` | 追跡の WSPAC のビーム行列: 既定 (ゼロ電流) と v2 のオプション `WSPACSC = 1` (空間電荷つき・MINCOUP の下限) の試験粒子のチューンシフト | gs (−2.938e-5、−1.210e-3) に既定が x 1.5%・y 0.7%。WSPACSC = 1 (下限つき) は y が 3 倍違う (K64 の追跡には下限が入らない)。k64 は追跡に WSPAC を入れない | 同上 |
| `P4_sol_photons` | v2: SOL 領域内の DRIFT・MULT・QUAD の光子。要素ごとの光子数・エネルギーの合計 (MULT は古典値 285 eV) | gs は MULT・DRIFT を記録し QUAD は記録しない。k64 は領域内の MULT が未対応 | 同上 |
| `P4_sol_spin_quad_mult.py` | 領域内の QUAD・MULT のスピンの Thomas-BMT の積分との差 (v2: 1e-6 以下) | gs は使わない | `.out.txt` (v2)・`.out_k64.txt` |

## 未確認・限界
- 大きなリング A (603 m) は DAPWIDTH だけ。PSPAC・SOL 領域・WSPAC の追跡は ASTRID2 (0.58 GeV) までの確認。PF リング・8 GeV のリングは使っていない。
- DAPWIDTH: NPARA > 1・`Damping->True` は未対応。損失の境界の粒子は gs と 1〜数周違うことがある。
- PSPAC (k64): gs の NPARA > 1 の再現、1 セッションの 2 回目の `TrackParticles` (gs は PSPAC が働かない) は再現しない。PSPACNZ < 3 はエラー。
- WSPAC の追跡: 強い空間電荷 (0.2 GeV で PBUNCH ≳ 1e9) は反復が収束しない。
- SOL 領域内の光子: BEND・SEXT・OCT・傾いた SOL の領域は未対応。領域内の加速 (POL + TRPT) は実装しない (記録 §2.6)。
