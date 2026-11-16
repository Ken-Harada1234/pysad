# GaussianCoulomb と BBBrem1

最終更新: 2026-11-16 (現状のスナップショット。経緯は `history/` にある)

## GaussianCoulomb[x, y, sigx, sigy] / GaussianCoulombU[x, y, sigx, sigy]

K64 の `tspac.f` (`bbkick1`、`cwerf`)、`tintrb.f` (`twspu`)、`bint.f` (`rombint`、`polint`)、`gamma.f` (`gamma0log`)。`pysad/gausscoul.py`。

- `GaussianCoulomb`: 大きさ (sigx, sigy) のガウス分布が (x, y) に作る力 `f_x + i f_y` (複素数。虚部が 0 なら実数で返す: gs と同じ)。Yokoya の `CWERF` (Gautschi) をそのまま写した。
- `GaussianCoulombU`: ポテンシャル。a+b<30 はロンバーグ積分 (k=5 の多項式外挿)、それ以上は漸近展開、r=1 は `gamma0log`。許容は (1e-10, 1e-13)。
- 検証 (`pysad/tests/test_gausscoul.py`): 15 点 (軸上・x 軸上・y 軸上・負の座標・σx<σy・σx≈σy (1e-5)・大きい引数 (漸近展開)・小さい引数) の実部・虚部が gs と 12 桁で一致 (表示の 12 桁が全部同じ)。
- **`GaussianCoulombFitted[x, y, sigx, sigy]`** = `{fx, fy, fu, fxx, fyy, fxy}` (`twspfu`) は実装済み。21×61×61 の表 (`twspu` と `bbkick1` から) を `spline1` で作り、3 次のスプラインで補間する (最初の呼び出しで表を作る、数秒)。gs と相対 1e-9 以内で一致 (試験 `test_gaussian_coulomb_fitted_matches_gs`)。K64 の `twspfuinit` は `spline1` の境界値に `ktaloc` の記憶を読む (uxx の 2 回目の呼び出し) が、0 とみなして gs と一致した。`sigx < sigy` のとき K64 は fxx と fyy を入れ替えない (そのまま写した)。

## BBBrem1[sqrt(s) (eV), ds/s]

K64 の `tbbbrem.f` (`bcube1`、`tbcube1`)。`pysad/bbbrem.py`。ビーム・ビーム制動放射の 1 事象 (Kleiss-Burkhardt)。
戻り値は `{(-1,-1,-1,1)*p2 (eV の 4 元運動量), 近似断面積 (mbarn×1e-31)}`。乱数は SAD の `tran()` (nran=3)。

- 検証 (`pysad/tests/test_bbbrem.py`): `SeedRandom[3]` のあと 6 事象、sqrt(s)=10,20,30,40 GeV の 4 事象で、乱数の列も含めて相対 1e-15 (最後の桁の丸め) で一致。
- パッケージ `Packages/bbbrem.n` (`BBBrem` クラス) は `Parallelize`・`LBoost`・`Rot3P` などで `BBBrem1` を使う。pysad の SADScript インタプリタで `Parallelize` が動くかは確かめていない (未対応とみなす)。
- 事象を棄却する (t < tmin) 場合に K64 が出す診断の print は写していない。

解説: docs/37_beambeam_guide.md (BEAM-BEAM。GaussianCoulomb・BBBrem はビーム・ビームの計算と同じ系)
