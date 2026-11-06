# 作業ログ: touschek.py (TouschekTable / TouschekLifetime の直訳) の解説と照合

日付: 2026-10-15 / 対象: `pysad/touschek.py` (861 行)。コードは未変更。照合は `bench_level17_pysad.py a|b` を実際に再実行した結果。

## 1. 何をするモジュールか

SAD で `INTRA` フラグを立てて `emit` (または `Emittance[]`) を実行すると、temit が 1 周追跡 (tturne) の各要素で `tintrb` を呼び、
(a) バンチ内散乱 (IBS) の拡散を beam(1:21) に足し、(b) 最後の収束パスでは Touschek の損失率表 (touckl / touckm / toucke) を積算する。
収束後に `TouschekTable` (SAD のリスト) が組み立てられ、`TouschekLifetime[nx, ny, nz]` (Packages/emit.n) がその表をスプライン補間して寿命を返す。
touschek.py はこの経路全体を直訳し、emit6d.py の `_State` / 要素写像 / `tecalc` を再利用している。

`touschek_table(lattice, momentum_eV, em, pbunch, nbunch, mincoup, rfsw, radcod)` が `emit` (INTRA ON) 相当で、
戻り値の `'table'` が `TouschekTable` と同じ構造
`[id, [dp(120), 1/τ(120)], [[nx(34), ny(34), nz(34)], mx(34×25), my(34×25)], toucke(nlat×120)]`。
`touschek_lifetime(table, nx, ny, nz, sige)` が `TouschekLifetime`。INF は `None` か `math.inf`。

## 2. SAD ルーチン ↔ Python 関数の対応表 (touschek.py のコメントによる)

| SAD (ファイル/名) | Python | 内容 |
|---|---|---|
| temit.f `touschek_table` モジュール (ntouckx=34, ntouckz=25, ntouckl=120) | 定数 `NTOUCKX/NTOUCKZ/NTOUCKL` | 表の大きさ。dp = 0.004〜0.242 (0.002 刻み) |
| temit.f `initialize_tampl` | `initialize_tampl()` | 振幅表 tampl(34,3) (x,y: 1〜88 と ×10、z: 4.6〜) |
| tintrb.f `touckf` / emit.n `TouschekF` | `touschek_f(x)` | C(x)/x の近似 (定数 EEULER, TA, TB)。x=0 で 1e200 |
| tintrb.f `fintrbu` | `fintrbu(u, a1, a2)` | IBS 積分の被積分関数 |
| bint.f | `bint(f, a, b, epslon, epsabs)` | 台形則細分 + Simpson 外挿。収束判定まで同じ反復をたどる |
| teigens33 (thess3+tqr3) | `_eigs33` | 3×3 固有値分解。LAPACK (eigh) は絶対誤差で σy² が 1% 狂うため、相対精度の出る 2 側 Jacobi 法で代用 (直訳ではない唯一の箇所) |
| sols33 | `_sols33` | 3×3 連立方程式 (対角で釣り合わせてから解く) |
| tmulbs (calint=T で beam(22:42) も変換) | `_IntraState.apply` | B2 ← J B2 Jᵀ |
| photons.f:1518 tradke (dbeam を beam(1:21) と beam(22:42) に加える) | `_IntraState.tradke` | B2_new = B_new + J (B2 − B) Jᵀ |
| temw / touschek_table / tmacro の変数 | `_IntraContext` | r, ri, eemx/y/z, tampl, touckl, touckm, toucke, id |
| tturne.f 先頭 (calint で表を 0 に) | `_IntraContext.reset_tables` | |
| tintrb.f (wspac なし、optics=.false.) | `_tintrb(st, al, al1, l, ctx)` | スライス中心のビーム行列から IBS 拡散 (bmi) と Touschek 表の積算 |
| tturne.f calint 分岐 (要素ごと tintrb、ドリフトの 2 分割) | `_tturne_intra` | emit6d.tturne の要素分岐 + tintrb 呼び出し |
| temit.f `eintrb` | `_eintrb` | IBS 反復の加速 |
| temit.f converge ループ + `tintraconv` | `touschek_table` | 収束判定 (ITMAX=100, RESIB=3e-6, DCMIN=1e-6, MINCOUP の下限)、収束後に caltouck=T でもう 1 周して表を作る |
| spline.f `spline` (mode=0) | `_spline_ddy` | 自然スプライン (ddy は y''/6) |
| spline.f `tffindindex` / D.n SplineData | `_find_index`, `SadSpline` | 区間探索と 3 次補間 (範囲外は端の区間で外挿) |
| emit.n `Touschek` クラス (Setup) | `TouschekLifetimeCalculator.setup` | spz, spxz, spyz スプラインの準備、範囲 (dpmin/max, nxmin/max, …) |
| emit.n `TouschekLifetime[nx,ny,nz]` | `TouschekLifetimeCalculator.lifetime`, `touschek_lifetime` | nx, ny の INF の組合せで 3 通りに分岐 (両 INF: dp 受容、片方 INF: 振幅-nz 2 次元補間)。両方有限は NaN |

なお docstring が参照する `work_touschek.md` は今回この文書として初めて作成した (固有値ルーチンを Jacobi にした経緯はコメントの記述のみ)。

## 3. 照合: bench_level17a / 17b

SAD 側: `bench_level17a.sad` は `RFSW; RAD; FLUC; INTRA; MINCOUP=0.019;` に `PBUNCH = 200/2.998/1.602*LINE["S","$$$"]/992*1e8`
(= 2.5318e9 粒子/バンチ)、`emit`、TouschekTable の各部と `TouschekLifetime` 15 ケース (TL1..TL15) を Print。`17b` は RADCOD を加えただけ。
ラティスは pysad/tests/kek-ls_v03_12.lat (nlat = 3564)。pysad 側は `bench_level17_pysad.py` をそのまま実行 (IBS 反復 it=5 で収束)。

所要時間: 17a は emittance6d 9.7 s + touschek_table 19.4 s、17b は 18.4 s + 14.3 s。

### 17a (RADCOD OFF): 全項目 1e-9 台以下 (ただし EMITY を除く)

| 項目 | 最大相対差 |
|---|---|
| EMITX / EMITZ / SIGE / SIGZ / EMITXC / EMITYC / EMITZC | 5.0e-12 / 2.2e-12 / 7.7e-13 / 1.9e-12 / 5.0e-12 / 1.7e-12 / 2.2e-12 |
| EMITY (≈1.2e-14、結合 0 近傍の小さな値) | 1.5e-7 |
| 表 [[2,1]] dp | 2.0e-16 |
| 表 [[2,2]] 1/τ | 3.3e-11 |
| 表 [[3,1]] 振幅 | 0 |
| 表 [[3,2]] mx / [[3,3]] my | 8.6e-12 / 6.7e-12 |
| 表 [[4]] 要素和 (120) | 3.3e-11 |
| 表 [[4]] dp 和 (非零 2901 要素) | 1.1e-10 (零要素は pysad も 0) |
| [[4,1]] / [[4,2]] / [[4,100]] / [[4,-1]] | 0 / 0 / 3.2e-11 / 7.8e-10 |
| TouschekLifetime TL1..TL15 | 全 15 ケース ≤ 6.7e-12 (TL15 は両者 0) |

例: TL1 [INF,INF,3] pysad 12.36769325191 / SAD 12.36769325192、TL9 [1000,INF,200] 4.614551171303e5 / 4.614551171326e5。
id、TouschekTable の長さと形 (4 要素、{3564, 120}) も一致。EMITY の 1.5e-7 は絶対値が 1e-14 の最終パス値で、
y 面が結合のみで決まる小さい量の丸め由来と見られる (EMITYC = 4.6e-12 のほうは 1.7e-12 で一致している)。

### 17b (RADCOD ON): 1e-7 台で一致

| 項目 | 最大相対差 |
|---|---|
| EMITX/Y/Z, SIGE, SIGZ, EMITXC/YC/ZC | 7.9e-8 〜 2.9e-7 |
| 表 [[2,2]] 1/τ, [[3,2]] mx, [[3,3]] my | 3.3e-7, 6.1e-7, 3.5e-7 |
| 表 [[4]] dp 和 (2901 要素) / [[4,-1]] | 1.1e-5 / 7.7e-6 |
| TL1..TL15 | 全 15 ケース ≤ 4.5e-7 (TL5 4.4e-7, TL6 4.5e-7 が最大) |

原因 (推定、ただし docs/10_benchmarks.md の Level 4 の既知事項と整合): RADCOD ON の emittance6d 自体が SAD に対し ≤1.1e-7
(SAD の tcod の収束誤差内、z と trf0 の分け方の違い) であり、IBS 反復は元の平衡エミッタンスから出発するため、その差がそのまま残る。
touschek.py 側の誤差ではない可能性が高い (17a で 1e-11 台なので表の計算自体は正しい)。1e-9 まで詰めるには tcod の収束基準 (codw0, 1e-4) を SAD と同じにする必要があるが、
SAD 側の収束が途中打ち切りなので「SAD の数値の再現」は原理的に 1e-7 程度が限界。

## 4. ffscmd.py への接続案 (今回は実装しない)

現状: `ffscmd.py` の `fn_touschek` は `SADError('TouschekLifetime is not implemented yet')` を投げ (877 行)、`B['TouschekLifetime']` は登録済み (222 行)。
`run_emit` (751 行) は INTRA フラグを無視して `emittance6d` のみ呼ぶ。`flags['INTRA']` は既に存在 (62 行)。

1. `run_emit`: `emittance6d` の後に `if self.flags['INTRA']:`
   `from .touschek import touschek_table` → `out = touschek_table(self.lat, self.momentum(), e, PBUNCH, mincoup=MINCOUP, rfsw=..., radcod=...)`。
   PBUNCH / MINCOUP / NBUNCH は `self.I.globals` から読む (SAD の既定: MINCOUP=1e-2、PBUNCH=0 なら IBS も Touschek も 0)。
2. 結果の反映: SAD は INTRA のとき EMITX/Y/Z, SIGE, SIGZ (最終パス) と EMITXC/YC/ZC (IBS 込み収束値) をグローバルにセットするので、
   `out['emittance_x']` 等と `out['EMITXC']` 等を `g[...]` に書き、`g['TouschekTable'] = out['table']` を保存。
   `self.touschek = out` にも持たせる (`Emittance[]` の戻り規則に影響しない)。
3. `fn_touschek(args)`: `nx, ny, nz` を評価 (`INF` は `sadscript.INF = float('inf')` なのでそのまま `touschek_lifetime` に渡せる)。
   表が未作成なら SAD と同じ `TouschekTable Undefined; Do Emittance[] with INTRA;` を `SADError` で投げる。
   `sige` は `self.I.get('SIGE')`。`touschek_lifetime(self.touschek['table'], nx, ny, nz, sige)`。
4. print_emit: INTRA 時の表示 (IBS 込みエミッタンス行の有無) は SAD の出力を確認してから決める。
5. 注意: 追加コストは 15〜20 s (IBS 反復 5 回)。`touschek_table` の `verbose=True` は反復ごとにログを出すので、FFS 経由では `verbose=False` か `self.write` 経由にする。
   結合 (EMITY) の最終パス値は 1e-7 相対しか合わないので、EMITY を使う後段の照合には EMITYC を使う。

## 5. 判定

touschek.py は RADCOD OFF で SAD と 1e-9〜1e-11 で一致 (TouschekLifetime 15 ケースは 7e-12 以下)、RADCOD ON では SAD の tcod 収束誤差 (1e-7) の範囲。
次の手は上記 ffscmd.py 接続と、PBUNCH=0 / MINCOUP 既定値での動作確認 (本ベンチは MINCOUP=0.019 のみ)。
