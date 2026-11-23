# v2.5.0 記録: SOL 領域内の QUAD の ΔH の運動量の回転の順序 (B333)

日付: 2026-11-23。試験 `pysad/tests/script_tests/test_v250_sol_quad_rotation.py` (7 件)。gs の保存した出力 `pysad/tests/script_tests/scripts/v250_sol_quad_rotation__{p05,p10}.{sad,gs.txt}`。独立なレビューの指摘 (ΔH の回転の更新順序) を gs で確認した。

## 1. 症状

`tsolqu.f:122-123` (`tsolqum` の `:302-303` も同じ):

    pxi  = a22*pxi + a24*py(i)
    pyi  = -a24*pxi + a22*py(i)

2 行目の pxi は 1 行目で更新された値。正しい回転は py' = −a24·px_old + a22·py。K64 の写像は (px, py) → (a22 px + a24 py, −a24(a22 px + a24 py) + a22 py) で、行列式 a22² = cos²(wℓ) (w = bz(−dpz)/(1+dpz)、ℓ = スライスの長さ)。位置 x・y の更新 (更新前の pxi) は正しい。1 スライスの角の 2 乗 (wℓ)² の分、運動量の面積が縮み、シンプレクティックでない。

## 2. gs での確認

3 GeV/c、NORAD・NOFLUC、`SOL SA=(L=0 BZ=13 BOUND=1) SB=(L=0 BZ=0 BOUND=1)` (bz = 1.3 /m)、`QUAD EL=(L=0.4 K1=0.144)` (K1 = 0.36 /m²)、LINE (ST SA EL SB)。点 (x, px, y, py) と各成分 ±1e-6 の 13 粒子の TrackParticles 1 要素から 6×6 ヤコビアンを作る (δ = z = 0)。

| 点 (x, px, y, py) | 量 | gs | k64 (pysad) | v2 |
|---|---|---|---|---|
| (1e-3, 0.05, −5e-4, 0.02) | 4×4 ブロックの \|JᵀSJ − S\| | 1.54e-7 | 1.54e-7 | 4.4e-10 |
| 同 | det J − 1 | −1.504e-7 | −1.504e-7 | −6.4e-11 |
| (1e-3, 0.1, −5e-4, 0.05) | 4×4 ブロック | 3.39e-6 | 3.39e-6 | 5.3e-10 |
| 同 | det J − 1 | −3.234e-6 | −3.234e-6 | −5.2e-10 |

k64 の 13 点の追跡は gs と全桁一致 (最大差 0)。v2 の残りは中心差分の丸め (1e-16/1e-6 の桁)。レビューの見積り (px = 0.1 で 4 回の回転の θ² の和が 3.0e-6) と det の大きさが合う。

## 3. 修正 (v2 のみ)

- `track6._tsolqu` (python): v2 は `pyi = -a24*pxi + a22*py` を `pxi` の更新前に計算。
- `track6n_ext.tz_step1` (numba): 引数 `v2f` (既定 False) を足し、`tsolqu_bz_all` が `rs[1]` (v2 の印) から渡す。`solemit` (SOL 領域の Emittance[] の写像) の呼び出しは `compat.v2()`。
- bz = 0 の分岐と、スライスの最後の半ステップ (`tsolqu_bz_final1`、元から更新前の pxi) は変わらない。

## 4. 限界・未確認

- 1 構成 (bz = 1.3 /m、QUAD 1 台、3 GeV/c) の 2 点。bz = 10 /m・px = 0.1 でレビューは 8.8e-5 を報告 (この記録では確かめていない)。
- v2 でも (y, δ)・(x, δ) のシンプレクティック性の欠陥が残る (px = 0.1 で 7.7e-6・1.7e-6、px = 0.05 で 8.4e-7・1.0e-8)。回転の更新順序とは別の原因 (z の更新の近似などと思われる。未調査)。
- `tsolqum` (K0 つき) と SOL 内の MULT の K1 は、同じ式を使う pysad の経路 (`tz_step1`) で直るが、gs では個別に確かめていない。
