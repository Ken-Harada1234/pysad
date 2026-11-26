# 12. pysad 公開 API リファレンス

**この文書で分かること**: pysad を Python から使うときの公開 API (関数・クラス・引数・戻り値) の一覧。§0 は単位・座標・位置の約束、§1〜15 はモジュール別の API、§16 は追加のモジュール。
**前提**: SAD のラティスと FFS の用語 (docs/01・02・03)。例はコピーして実行できる。
**読む順**: §0 → 使うモジュールの節 (§1 parser、§2 elements …)。互換モードは §16.6。

(§1〜15 に無い API は §16。互換モード `pysad.compat` と v2 のモジュールは §16.6、`naff`・`specfn_cx`・`nbcache` は §16.7、周波数マップ・LMA・共鳴駆動項は §16.8)

> 場所: `pysad/`。全体像と文書の案内は `01_overview.md`、各計算の直訳メモは 07〜11、SAD の癖は `29_sad_quirks.md`。
> 本書の例はすべて `PYTHONIOENCODING=utf-8 python -u` で実行して動作を確認してある (ラティスは `<試験用ラティス>.lat` と本文中の小さな FODO リング)。
> [`pysad/__init__.py`](../pysad/__init__.py) の `from pysad import ...` は旧 API (§15) を export しているので、**各モジュールを直接 import する** (`from pysad.parser import load_sad_file` など)。
> **互換モード**: 既定は `v2` (K64 のバグを直した動作)。gs (K64) と同じ値が要るときは `PYSAD_COMPAT=k64` か `pysad.compat.set_mode('k64')` (§16.6)。本書の「K64 の動作」の記述は k64 モードのもの。モードで変わる箇所は docs/28。

## 0. 共通の約束 (単位・座標・位置)

| 項目 | pysad | SAD | 備考 |
|---|---|---|---|
| 運動量 `momentum_eV` | eV/c | `MOMENTUM` (`3GEV` = 3e9) | 各関数に明示的に渡す。`FFS` 層は `MOMENTUM` グローバルを使う |
| QUAD/BEND `K1` | 長さあたり [m⁻²] = SAD の K1 / L | 積分値 [m⁻¹] | parser が割る。L=0 なら SAD の値のまま。元の値は `el._raw['K1']` |
| SEXT `K2`, OCT `K3`, DECA `K4` | 長さあたり (SAD の値 / L) | 積分値 | 同上 |
| BEND `ANGLE` | rad | rad | |
| BEND `E1`, `E2` | 端部角 [rad] = SAD の E1 × ANGLE (+ AE1) | ANGLE に対する比 | SAD と定義が違うので注意 |
| 回転 | 属性名 `ROTATION` [rad] | `ROTATE` / `DROTATE` | BEND の `ROTATE` は設計ロール角 (`Bend.ROTATE`)。ANGLE=0 でも ROTATE≠0 の BEND (VD/HD 等) は `Drift` にせず `Bend(ANGLE=0)` のまま保持する |
| CAVI `VOLT`, `FREQ`, `PHI` | eV, Hz, rad | 同じ | `HARM` は `build_lattice` で `FREQ = HARM·c/C` に変換 |
| 粒子座標 | `(x, px/p0, y, py/p0, z, δ)` の行ベクトル、配列形 `(np, 6)` | TrackParticles の列 | z = −cΔt [m]、δ = Δp/p0 |
| Twiss 表の位置 | `rows[0]` = 入口 `*`、`rows[i]` = i 番目要素の出口 (= i+1 番目の入口)、`rows[-1]` = `$$$` | `twiss(l)` (l = 1..nlat) | pysad の添字 i = SAD の l − 1。長さは nelem+1 |
| 位相 `mux`, `muy` | rad (累積) | `NX`, `NY` (rad) | チューン = `mux/2π` |
| 要素位置の指定 (`Matcher.locate`, FFS) | `'QF'` (最初の 1 台), `'QF.3'` (3 台目), `'$$$'`/`'***'` (末尾), `'^^^'` (先頭) | 同じ | |

SAD の要素ファミリー vs 誤差 (インスタンス) の区別は `Lattice.set_element` (全インスタンス) と `Lattice.set_error` (1 台) で行う。

## 1. parser.py — MAIN 部 (要素定義・LINE) の読み込み

対応: `toplvl.f` / `doelem` / `doline` / `expln` (要素定義、`LINE`、`N*`、`-LINE` の逆順、再定義は後勝ち)。

| 関数 / クラス | 戻り | 説明 |
|---|---|---|
| `load_sad_file(path, verbose=False)` | `SADScript` | .lat/.n ファイルを読む |
| `load_sad_string(source, verbose=False)` | `SADScript` | 文字列から読む (ファイル不要) |
| `load_sad_files(*paths, verbose=False)` | `SADScript` | 複数ファイルを順に同じ環境に読む |
| `SADScript.build_lattice(name)` | `Lattice` | LINE 名 (大文字) を展開して `Lattice` を作る。周長が決まった後に CAVI の HARM → FREQ |
| `SADScript.element_names()` / `beamline_names()` | `list[str]` | 定義済みの要素名 / LINE 名 |
| `SADScript.get_variable(name)` / `set_variable(name, v)` | | MAIN 部の変数 (例 `MOMENTUM`) |
| `SADScript.define_element(name, etype, **params)` / `define_beamline(name, refs)` | | Python から要素・LINE を追加 |
| `SADScript.env` | `SADEnvironment` | 生の定義辞書 `env._elements[name]` (`'TYPE'` と SAD のキーワード値)、`env._beamlines` |

```python
from pysad.parser import load_sad_file
sc = load_sad_file(r'<試験用ラティス>.lat')
lat = sc.build_lattice('RING')          # Lattice('LATTICE', nelem=3563, C=603.120 m)
qf = [e for e in lat.elements if e.name == 'QFA1H'][0]
qf.K1, qf._raw                          # 長さあたりの K1、{'L': .125, 'K1': ...} (SAD の積分値)
```

## 2. elements.py — 要素クラス

`@dataclass` の `Element` を基底に `Drift, Quadrupole, Bend, Sextupole, Octupole, Decapole, Dodecapole, Multipole, Cavity, TCavity, Solenoid, Monitor, Marker`。共通属性 `name, L, K0, DX, DY, ROTATION` (`Element.ROTATE` は垂直キック用で SAD の ROTATE ではない)。`make_element(etype, name, **kw)` は SAD 型名 (`'QUAD'` 等) から生成する。`transfer_matrix(delta)` は 6×6 線形行列 (旧 API 用。CALC/追跡は lattice.py / track6.py の非線形写像を使い、この行列は使わない)。

| クラス | 固有属性 (単位) | 追跡・CALC で使う SAD キーワード (属性として setattr 可) |
|---|---|---|
| `Quadrupole` | `K1` [m⁻²] | `DX, DY, ROTATION, F1, F2, FRINGE, DISFRIN, K0` |
| `Bend` | `ANGLE` [rad], `K1` [m⁻²], `E1, E2` [rad], `HGAP, FINT` (未使用) | `K0, DX, DY, ROTATION (DROTATE), F1, FRINGE, FB1, FB2` |
| `Sextupole` / `Octupole` / `Decapole` | `K2` / `K3` / `K4` (長さあたり) | `DX, DY, ROTATION, K0` |
| `Cavity` | `VOLT` [eV], `FREQ` [Hz], `PHI` [rad] | `L` (≠0 も可), `_harm` |
| `Multipole` | `KnL` 辞書 | `K0, K1, K2, K3, DX, DY, ROTATION` (ANGLE=0 のみ) |
| `Solenoid` | `BZ` [T] | 追跡のみ |
| `Marker` | — | 先頭 MARK の `AX, BX, …, EMITX, EMITY, DP, SIGZ, SIGE` は `env._elements` から読む (INS の入口 Twiss、`Matcher` の既定) |

## 3. lattice.py — CALC (閉軌道 + Twiss)

対応: `qcell.f qcell1` → `qtwiss.f` (`qcod` 4D 閉軌道 Newton+SVD、`qtwiss`/`qmat2twiss` 結合込み伝搬)、`INTGRL.f`。

| メソッド | 戻り | 説明 |
|---|---|---|
| `Lattice(elements, name='LATTICE')` | | `lat.elements` (list), `lat.circumference`, `lat.nelem`, `lat[i]`, `len(lat)` |
| `twiss_table(delta=0.0, cod=None, full=True)` | `list[TwissRow]` (nelem+1) | **CELL; CALC** 相当。δ = DP0。閉軌道の初期値は直前の結果 (SAD の twmov)。`cod` で初期値を与える |
| `twiss(delta=0.0, full=True)` | `list[dict]` | 同じ計算の辞書版 (キー `ax, bx, nx, ay, by, ny, ex, epx, ey, epy, r1..r4, detr`) |
| `tune(delta=0.0)` | `(νx, νy)` | `Twiss["NX","$$$"]/2/Pi` |
| `find_cod(delta=0.0)` | `ndarray(6,)` | 入口の閉軌道 (`qcod`)。`lat._cod_prev` に残る |
| `ins_twiss(elements, twiss0, cod0=None, delta=0.0)` (モジュール関数) | `list[dict]` (各要素に `'cod'`) | **INS** 相当。`twiss0` は `ax, bx, ay, by, ex, epx, (ey, epy, r1..r4)`、`cod0` は `(x, px, y, py)` |
| `set_element(family, key, value)` | 台数 | ファミリー全台に `setattr` (値は pysad 単位: K1 は長さあたり) |
| `set_error(instance, key, value)` | `Element` | `'QD.3'` 形式で 1 台だけ。SAD の `QD.3 DX 1E-4` |
| `transfer_matrix(delta)`, `partial_transfer_matrix(i_from, i_to, delta)` | 6×6 | 線形行列の積 (旧 API 系) |
| `chromaticity()`, `summary()`, `radiation_integrals()`, `emittance(energy_eV)` | | `chromaticity()` は dν/dδ (`twiss_table(±1e-4)` の位相進みの差分。BEND・6 極を含む。追跡のヤコビアン・gs の `Twiss["NX",-1]` の DP0 差分と一致。3 GeV の 1 リング構成で確認、docs/30 B325)。他は旧 API 系で SAD 未照合。FFS の `CHROMATICITY` は QUAD・SEXT・MULT だけ・単位 4π·ξ |

`TwissRow` の属性 (SAD の `twiss(l, 0, 1:28)` と同じ順): `s, name, alfx, betx, mux, alfy, bety, muy, etx, etpx, ety, etpy, r1, r2, r3, r4, detr, x_cod, px_cod, y_cod, py_cod, z_cod, delta_cod, alfz, betz, muz, zx, zpx, zy, zpy`。`to_dict()` で `ax, bx, …` の別名も得られる。4D CALC では `alfz..zpy` は 0 (これらは `calc6d` が埋める)。

```python
tab = lat.twiss_table()                       # 3563 要素で 0.1 s
end = tab[-1]
print(end.mux / 2 / math.pi, end.muy / 2 / math.pi)       # 47.580034, 17.620044
mkc = [r for r in tab if r.name == 'MKC']     # 同名要素の全出口 (40 行)
print(mkc[0].s, mkc[0].betx, mkc[0].etx)
tab2 = lat.twiss_table(delta=0.01)            # DP0=0.01 のオフモメンタム
```

## 4. emit6d.py — Emittance[] (temit)

対応: `temit.f` (temit/tecalc/synchm), `tturne.f` (tquade/tbende/tthine/tcave/tmulte), `photons.f tradke`。放射込み 6D 一周写像 (T: シンプレクティック部, D: 減衰部, B: 量子励起) → 固有値分解 → 平衡ビーム行列。

`emittance6d(lattice, momentum_eV, charge=1.0, amass=ELMASS, close_z=True, rfsw=False, radcod=False, fshift=0.0) -> dict` (`fshift` = SAD の FSHIFT (RFSW ON のとき dv に足す)。`emittance_rules`・`calc6d` にも同じ `fshift=`)

| キー | 内容 | SAD の Emittance[] 規則名 / EMIT 表示 |
|---|---|---|
| `emittance_x`, `emittance_y`, `emittance_z` | m | `Emittances`, `EMITX/EMITY/EMITZ` |
| `momentum_spread`, `bunch_length` | σδ, σz [m] | `MomentumSpread`, `BunchLength` (`SIGE`, `SIGZ`) |
| `momentum_compaction` | αp | `MomentumCompaction` |
| `energy_loss_eV` | U0 [eV] | `EnergyLossU0` |
| `rf_voltage_eff`, `trf0`, `bucket_height` | Vc [eV], 同期位置 [m], ΔV/P0 | `RfVoltageVc`, `EquilibriumPosition`, `BucketHeight` |
| `closed_orbit` | `ndarray(6,)` 6D 閉軌道 (RFSW+RADCOD は tcod の入口の値、他は一周後の出口) | `OrbitAtExit` |
| `tunes`, `damping_rate` | 長さ 3 (x, y, z)。tunes は固有値の偏角/2π (分数部、符号付き) | `Tunes` (Real tune), `DampingRate` |
| `beam_matrix`, `one_turn_matrix` | 6×6 | `EquilibriumBeamMatrix`, `OneTurnTransferMatrix` |
| `normal_r`, `normal_ri` | 6×6 正規モード → 物理 (と逆) | temw の r/ri (DA の初期粒子配置に使う) |
| `cd`, `omegaz`, `heff`, `omega0`, `h0`, `ext_twiss` | tecalc の複素 cd (1:3 = 減衰 + i チューンシフト、4:6 = log λ)、シンクロトロン角周波数、実効ハーモニック数、ω0、h0、`tfetwiss(ri, cod)` の 28 個 | EMIT 報告 (`print_emit`) の Extended Twiss, fs, h, Imag.tune, Damping time, Tune shift, 分配数 |

`rfsw=True` が SAD の `RFSW; emit`。`radcod=True` は RADCOD (放射込みの閉軌道)。`trf0` は追跡 (`track_particles` の `trf0`) に渡す。下位関数 `tturne(elements, cod0, ph, trf0, rfsw, radcod)` と `tecalc(state, circ, rfsw)` も公開。

```python
from pysad.emit6d import emittance6d
em = emittance6d(lat, 3.0e9, rfsw=True)      # 3563 要素で約 1.6 s
em['emittance_x'], em['momentum_spread'], em['bunch_length']   # 1.549e-10, 6.947e-4, 3.773e-3
```

下位 API: `closed_orbit6d(els, ph, rfsw, radcod, close_z=True, norf_trf0=False) -> (cod, trf0, found, st_full)` (temit の tcod 部分。`emittance6d` と `calc6d._closed_orbit` が共用)、
`_tcod_radcod(els, ph, cod, rfsw)` (tcod.f の直訳、RADCOD ON)、`beam_tilt_sizes(beam_matrix) -> (tilt, sigx, sigy)` (EMIT の Beam tilt / Beam size)、`tturne(..., want_db=True, want_jac=True)`
(`want_jac=False` は cod だけ追う一周、`want_db=False` は D・B なし = SAD の irad=6)、`emit6d._USE_ANALYTIC_JAC` (False で解析ヤコビアンを使わず数値微分)。

## 5. match6.py — FIT / FREE / FIX / GO (tffsmatch)

対応: `tffsmatch.f` (Newton/DESCEND、fuzz、fmincube、tweigh、tfsolv の SVD)、`tffscalc.f` (運動量点ループ、twfit の重み、nstab)、`qdcell.f tdfun`、`tgetfv.f`。微分は既定で解析微分 (match6d.py、§6)。`numerical_derivative=True` で数値微分。

`Matcher(lattice, cell=True, dp0=0.0, dpmax=None, emx=None, emy=None, itmax=40, convergence=1e-9, exponent=2.0, offmw=1.0, out=print)` — `dpmax`/`emx`/`emy` の既定は先頭 MARK の `DP`/`EMITX`/`EMITY` (SAD と同じ、DP の既定 0.01)。

| メソッド | SAD | 説明 |
|---|---|---|
| `free(pattern, key=None)` | `FREE QF* K1` | パターンは SAD のワイルドカード、カンマ区切りで複数。key 省略時は型の既定 (QUAD K1, BEND ANGLE, SEXT K2, DRIFT L, CAVI VOLT…) |
| `fix(pattern='*', key=None)` | `FIX` | 変数から外す |
| `fit_at(loc=None, loc1=None)` | `FIT loc [loc1]` | 以後の条件の位置 (省略で `$$$`) |
| `fit(fun, value, npoints=None, loc=None, loc1=None)` | `FIT [loc] NX 47.6 [n]` | `fun` は `'AX','BX','NX',…,'LENG','CHI1'` (`_FUNS`)、末尾 `M` で最大フィット。`value` は数値または `(min, max)`。NX/NY/NZ は 2π 倍、CHI は deg→rad |
| `fit_clear()`, `fitp(n)` | `FIT;`, `FITP n` | |
| `go(calc_only=False)` | `GO` / `CALC` | 反復して `_Res` (`.r` 残差, `.nstab` 不安定数) を返す。途中経過は `out` へ SAD と同じ書式 |
| `values()` | | `{'QF K1': 値}` — **SAD 単位 (積分値)** |
| `iterations`, `matching_residual`, `unstable_level` | `MatchingResidual`, `UnstableLevel` | |
| `fit_value`, `fit_weight`, `fit_function` | `FitValue`, `FitWeight`, `FitFunction` | コールバック (関数を代入) |

```python
from pysad.match6 import Matcher
m = Matcher(lat, out=print)
m.free('QF,QD')                  # 既定キー K1
m.fit('NX', 2.18); m.fit('NY', 1.84)
res = m.go()                     # ' Matched. ( 2.5279E-16) DP = 0.01000 ...' を印字
res.r, res.nstab, m.iterations, m.values()   # {'QF K1': 0.25215, 'QD K1': -0.20305}
lat.twiss_table()[-1].mux / 2 / math.pi      # 2.18000000
```

## 6. match6d.py — 解析微分のマッチング (tffsqu / qdcell)

対応: `qdcell.f` (tffsqu, qdcell, dgammag), `qdtwis.f`, `qdquad.f`, `qdbend.f`, `qdmdia.f`, `qgettr.f`, `tftmat.f`, `qtent.f`, `qchg.f`。SAD の既定経路 (`FFS$NumericalDerivative=False`) の直訳。

| API | 説明 |
|---|---|
| `AnalyticMatcher(lattice, ...)` | `Matcher` のサブクラス (引数・メソッドは同じ)。解析微分は `Matcher` の既定 (`Matcher._derivative`) なので、互換のための別名。`numerical_derivative=True` で数値微分 |
| `compare_derivatives(m, verbose=True) -> (qu_analytic, qu_numeric, ids)` | 現在の変数・条件で解析微分と数値微分を並べて比較 (検証用) |
| `analytic_qu(m, ids, rows, wvar) -> (qu, numeric_cols)` | 微分行列 qu[条件, 変数] と、解析微分できず数値に落ちた列 |

制限: CELL のみ検証済み。変数は QUAD K1/ROTATE, BEND ANGLE/K1 (他は数値微分列にフォールバック)。L=0 の QUAD、MARK の変数は微分 0。

## 7. track6.py / track6n.py — TrackParticles (tturn)

対応: `tftrack.f`, `tturn.f`, `tquad.f`/`tsolqu.f`, `tthin`, `tbend.f`/`tbendi.f`, `tcav.f`, `tmulti.f`, `tsol`, `photons.f` (tradk / tradkfn)。要素外側・粒子内側のループ、FLUC の乱数を引く順序も SAD と同じ。

`track_particles(elements, coords, nturn, momentum_eV, trf0=0.0, rfsw=True, rad=False, fluc=False, seed=None, charge=1.0, amass=None, progress=0) -> ndarray(np, 6)` (track6: 純 Python)
`track_particles_fast(...)` (track6n: numba、同じ引数・戻り。MULT/SOL を含むと track6 へフォールバック)

| 引数 | 意味 |
|---|---|
| `elements` | `lat.elements` (始点 = 終点 = ラティスの先頭 = `TrackParticles[{1, beam}, 1, 1, nturn]`) |
| `coords` | `(np, 6)` の `(x, px/p0, y, py/p0, z, δ)` |
| `trf0` | `emittance6d` の `trf0` (SAD の global trf0。emit を先に実行するのと同じ) |
| `rfsw, rad, fluc` | FFS フラグ RFSW / RAD / FLUC (FLUC は RAD と併用) |
| `seed` | `SeedRandom[seed]`。None なら共有乱数 `track6.RNG` の現在状態から続ける |
| `progress` | n 周ごとに経過を印字 (0 で無し) |

戻りは nturn 周後の座標 `(np, 6)`。発散すると `FloatingPointError` (何周目かを含む)。損失判定 (tapert) は TrackParticles でも行う (落ちた粒子は alive=0、座標は落ちた時点の値)。

```python
import numpy as np
from pysad.track6n import track_particles_fast
x0 = np.array([[1e-3, 0, 1e-4, 0, 0, 0]])
x = track_particles_fast(lat.elements, x0, 20, 3.0e9, trf0=em['trf0'], rfsw=True, rad=True, fluc=False)
# 20 周 0.5 s (3563 要素、初回は numba のコンパイル込み)。track6.track_particles と 1e-15 で一致
x = track_particles_fast(lat.elements, x0, 10, 3.0e9, trf0=em['trf0'], rad=True, fluc=True, seed=12345)
```

## 8. dapert.py — DynamicApertureSurvey (trackd)

対応: `Packages/dapert.n`, `tfdapert.f`, `trackd.f` (スロット順、tapert の損失判定、Score、62 進の地図)。

`dynamic_aperture(elements, range3, r1, nturn, phi, ivar1, ivar2, momentum_eV, emit, emx, emz, trf0, rfsw, rad, fluc=True, alost=1.0, zlost=100.0, out=None, progress=0) -> (score, result)`

| 引数 | 意味 (`Axes->"ZX"` 既定なら ivar1=3 (Z), ivar2=1 (X), ivar3=2 (Y)) |
|---|---|
| `range3` | `[[a2min, a2max], [a3min, a3max]]` (ivar2, ivar3 の規格化振幅、51 点を同じ j で動かす) |
| `r1` | ivar1 の値のリスト (Z なら δ の値そのもの) → 行になる |
| `phi` | `[PhaseX, PhaseY, PhaseZ]` (既定 `[0, 0, -π/2]`) |
| `emit` | `emittance6d` の戻り (`normal_r`, `closed_orbit` を使う) |
| `emx`, `emz` | `sqrt(\|EMITX\|+\|EMITY\|)`、RFSW なら `sqrt(\|EMITZ\|)`、さもなければ `\|SIGE\|` |
| `alost`, `zlost` | `LOSSAMPL`, `LOSSDZ` |
| `out` | 地図の出力先 (`sys.stdout.write` など)。None で表示しない |

戻り: `score` (= `Σ nsc`) と `result = [[[a2min,a2max],[a3min,a3max],[a1min,a1max]], [[a1, nsc, [ntloss×51]], ...]]` (SAD の `ResultOfDAPERT`)。`ntloss` は損失までの完了周数 (生存は nturn)。未対応: `Damping->True`、`DAPWIDTH` の窓。

```python
from pysad.dapert import dynamic_aperture
emx = math.sqrt(abs(em['emittance_x']) + abs(em['emittance_y'])); emz = math.sqrt(abs(em['emittance_z']))
score, result = dynamic_aperture(lat.elements, [[0, 40], [0, 20]], [-0.02, 0.0, 0.02], 50,
                                 [0.0, 0.0, -math.pi/2], 3, 1, 2.5e9, em, emx, emz, em['trf0'],
                                 rfsw=True, rad=False, out=sys.stdout.write)
#  Range    Xmin:   0.000 Xmax:  40.000 ... / -0.02 48 UUUU...  /    Score:   144
```

### 8b. maintrack.py — MAIN の TRACK コマンド (ActTra/track/tracka)

対応: `ActTra.f`, `track.f`, `tracka.f`, `tinip` (trackd.f), `tsmear.f`, `wjfit.f`, `wtunem.f`, `tcftr.f`, `tapert.f` の損失報告 (docs/history/work_main_track.md)。

```python
from pysad.parser import parse_main_commands          # MAIN 部の ON/OFF, EMIT, `NAME=値`, TRACK, STOP を順に読む
from pysad.maintrack import run_track                  # run_track(ffs, spec, write) → dict
# 普通は run_sad.py で MAIN だけのスクリプトを実行する (MAIN に TRACK 文があると命令を書かれた順に実行)
```

戻りの dict: `alive` (生存数)、`lost_turn[ID]`、`final` (ID 順の最終座標、落ちた粒子は落ちた時点)、`smear` ((周, ID, [2Jx,2Jy,2Jz,σx,σy,σz]) のリスト)、
`tunes` (`TunesByTracking` と同じ: 粒子ごと `[生存, [面 x,y,z ごとに 8 本の [周波数, 振幅]]]`)、`losses` ((ID, 周, 要素番号))、`emx`, `emz`, `codin`。
`FFS.track_results` に TRACK ごとの dict が溜まる。`PYSAD_PROGRESS=100` で 100 周ごとに進捗を表示。

## 9. sadrandom.py — SAD の乱数

対応: `random_plugin_sad.c` (合同法 seed ← 48828125·seed mod 2³², Box–Muller)、`drndsr.f`。

| API | 説明 |
|---|---|
| `SadRandom(seed=17)` / `.restore(seed)` | `SeedRandom[seed]` (奇数化)。FFS 開始時は 17 |
| `.uniform0()` (= `.tran()`) / `.uniform()` | [0,1) / (0,1) の一様乱数 (`Random[]`) |
| `.gauss(n=1, gcut=1e300)` | `GaussRandom[n]` (GCUT 棄却)。n=1 なら float、それ以外 list |
| `.dump()` | `(seed, filled, next)` = `SeedRandom[]` の戻りの数値部 |
| `drndsr(rng)`, `tdusrn(rng, anp)` | 放射光子のエネルギー分布 / ポアソン光子数 (FLUC で使う) |

追跡・DA・FFS は 1 本の共有乱数 `pysad.track6.RNG` (`SadRandom`) を使う (SAD も 1 系列)。`track_particles(seed=…)` は `RNG.restore(seed)` と同じ。

## 10. sadscript.py — SADScript インタプリタ

対応: `tfetok.f` (字句), `tfeval.f` / `tfstk.f iprior` (構文・優先順位), `tfefun*` (組み込み), `autofg.f` (`$FORM`), `Packages/Fit.n` (PolynomialFit)。

| API | 説明 |
|---|---|
| `Interp(out=None)` | インタプリタ。`.globals` (dict), `.builtins` (name → f(args)), `.cwd` (Get/OpenRead の基準) |
| `Interp.run(src, mode='ffs')` | 文字列を実行して最後の値を返す。`'ffs'` では先頭語が FFS コマンドの文を `ffs_handler` に渡す、`'script'` は SADScript のみ |
| `Interp.get(name, default)` / `.assign(name, v)` | グローバル変数の読み書き (実数は float、True/False も 1.0/0.0、リストは list、文字列 str) |
| `fmt_real(x, form='')` | `autofg` の直訳: `''` は最短表現、`'15.10'`, `'S15.10'`, `'F10.5'`, `'E12.5'` |
| `fmt_g(x, w=11, d=4)` | Fortran `1PGw.d` (GO の残差表示) |
| `parse(src)`, `tokenize(src)`, `strip_comments(src)` | 構文木 (`Expr`, `Sym`, `Rule`, `Func`) |
| `SADError` | 実行時エラー |

単独で使うより `ffscmd.FFS.install(interp)` で `Twiss[]` 等を登録した状態 (= `run_sad_string`) で使う。

```python
from pysad.sadscript import Interp, fmt_real
I = Interp(); I.run('a = Table[i^2, {i, 3}]; Print[a];', 'script'); I.get('a')   # [1.0, 4.0, 9.0]
fmt_real(1.2345678901234e-3, '15.10')   # '    .0012345679'
```

## 11. ffscmd.py — FFS コマンド層 (tffsa) と run_sad

対応: `tffsa.f` (1 語の解釈順: tfprint → tflag → tgetfv → terror → 要素名 → コマンド → tfgetv), `tftwiss.f` (Twiss[]/LINE[]/Element[]), `tflag.f`, `terror.f`, `tftype.f`。

| API | 説明 |
|---|---|
| `run_sad(paths, out=None, cwd=None) -> (FFS, Interp)` | 複数ファイルを連結して実行 (`gs lattice script` と同じ)。`cwd` 既定は最後のファイルのディレクトリ |
| `run_sad_string(src, out=None, cwd=None) -> (FFS, Interp)` | 文字列版。`FFS;` 行で MAIN 部と FFS 部に分け、MAIN 部の `ON/OFF flag;` をフラグに反映 |
| `split_main_ffs(src) -> (main, ffs_src, flags)` | 上の分割だけ |
| `FFS(script_env, interp, out=None)` | `parser` の `env` と `Interp` を結び、組み込み `Twiss, LINE, Element, FFS, Emittance, TrackParticles, SeedRandom, GaussRandom, Random, DynamicApertureSurvey, FFSCalc, BeamLineName` のほか、`CalculateOptics, TouschekLifetime, Master, RenumberElement, ResponseMatrix, ParabolaRandom, SynchroBetaEmittance1` などを登録 (`install`)。全体の一覧は §16.1 |
| `FFS.execute(stmt)` / `FFS.run_text(text)` | FFS コマンド 1 文 / コマンド列 (`FFS["..."]`) |
| `FFS.use(name)`, `FFS.calc()`, `FFS.run_emit(quiet)`, `FFS.go(fit)` | `USE`, `CALC` (→ `FFS.rows`: `TwissRow`、`rows[l-1]` = SAD の `twiss(l)`), `EMIT`, `GO`/`CALC` |
| `FFS.twiss_value(key, k, fr=0.0)`, `FFS.line_value(key, k, fr)` | `Twiss[key, k+fr]`, `LINE[key, k+fr]` (k は 1 始まり、fr は要素内の割合) |
| `FFS.set_family(fam, key, v)`, `FFS.set_value(el, key, v)` | `QF K1 v` (SAD 単位) / 1 台だけ |
| `FFS.flags` (dict), `FFS.lat`, `FFS.els`, `FFS.nlat`, `FFS.emit`, `FFS.trf0`, `FFS.matcher`, `FFS.rng` | 状態 (SAD のフラグ・ラティス・Emittance 結果・乱数) |

FFS 層に接続済みの機能: `CalculateOptics`、`CALC6D` フラグ (CELL のみ)、`INTRA` + `emit` + `TouschekLifetime[]`、`LINE["GX"..]`/`ORG`、`Emittance[OneTurnInformation/Orbit/Matrix->True]` と `LINE["SIG..",l]`、解析微分 (既定)。詳細は docs/history/work_integration.md。

```python
import io
from pysad.ffscmd import run_sad_string
script = open(LAT, encoding='utf-8').read() + '''
FFS;
USE RING; CELL; CALC;
Print["NX = ", Twiss["NX","$$$"]/2/Pi, "  BX at MKC: ", Twiss["BX","MKC"]];
Print["K1 of QFA1H: ", LINE["K1","QFA1H"], "  Element: ", Element["K1","QFA1"]];
FFS["fit NX 47.6; fit NY 17.6; free QFA4 QDA5; go;"];
RFSW; emit;
Print["EMITX = ", EMITX, "  SIGE = ", SIGE];
'''
buf = io.StringIO()
ffs, I = run_sad_string(script, out=buf)     # 約 11 s (CALC + GO + emit)
print(buf.getvalue())                        # SAD と同じ書式の GO 表・EMIT 報告・Print の出力
I.get('EMITX'), ffs.flags['RFSW'], ffs.twiss_value('BX', ffs.nlat), ffs.rows[-1].betx
```

出力の抜粋 (SAD と同じ書式):

```
 Matched. (  0.000    ) DP = 0.01000  DP0 = 0.00000  ExponentOfResidual = 2.0  OffMomentumWeight =   1.000  CALEXP
$$$         AX    #######  #  -2.08E-14 $$$         BX    #######  #   4.643203 $$$         NX    #######  #  47.580034
NX = 47.580034238511544  BX at MKC: 1.0923652872679401
Energy loss per turn U0 =  .2960749 MV  Effective voltage    Vc = 1.4000000 MV
Emittance X            = 1.5494E-10 m   Emittance Y            = 1.2747E-21 m
```

### run_sad.py (<repo> 直下) — コマンドライン

```
cd "<repo>"
PYTHONIOENCODING=utf-8 python -u run_sad.py <試験用ラティス>.lat <台本ファイル>
```

引数のファイルを順に連結して `run_sad_string` に渡す。SAD の `gs` に複数ファイルを渡すのと同じ。進捗表示は環境変数 `PYSAD_PROGRESS=n` または SADScript 変数 `PYSAD$PROGRESS` (DA の n 周ごと)。

pip で入れたときのコマンド (`pysad/cli.py`): `pysad-run lattice_file [script_file ...]` (上と同じ)、`pysad-share` (文書・例のフォルダの表示)、`pysad-warmup [--check | --clear] [--quiet]` (numba のカーネルを先に機械語へ変換してキャッシュを作る。最初の実行の約 1 分の待ちがなくなる。`--check` はキャッシュが有効なら終了コード 0、`--clear` は pysad のキャッシュを消す。`pysad/warmup.py`)。

## 12. calc6d.py — CALC6D (qcell6d)

対応: `qcell.f qcell6d` → `temit` (tcod, tturne, tecalc) → `tsetetwiss` / `tfetwiss` (`tffs.f`) / `tfadjustn`、分数位置は `qtwiss.f qtwissfrac`。

| API | 説明 |
|---|---|
| `calc6d(lattice, momentum_eV, rfsw=False, radcod=False, rad=True, charge=1.0, amass=ELMASS, twiss0=None, radiation=None, out=None, with_emit=True) -> Calc6DResult` | `CELL; CALC6D; CALC;`。SAD と同じく RAD フラグに依らず放射減衰込み |
| `Calc6DResult` | `rows` (nelem+1 の `TwissRow`、28 成分すべて: `alfz, betz, muz, zx..zpy` も入る), `stable`, `closed_orbit`, `trf0`, `r`, `ri`, `tunes` (νx, νy, νz の分数部), `emit` (tecalc の dict = Emittance[] 相当), `messages` |
| `twiss_frac(res, l, fr) -> TwissRow` | 要素 l (1 始まり) の入口から割合 fr の位置 (`Twiss[key, 200.5]`) |
| `tfetwiss(r, cod, normi)`, `etwiss2ri(tw)` | 拡張 Twiss 28 成分 ↔ 正規化行列 |

```python
from pysad.calc6d import calc6d, twiss_frac
r6 = calc6d(lat, 2.5e9, rfsw=True)
r6.stable, r6.tunes, r6.rows[-1].betz, r6.rows[-1].muz / 2 / math.pi
twiss_frac(r6, 3, 0.5).betx        # 3 番目の要素の中央
```

## 13. touschek.py — INTRA (IBS + Touschek 寿命表)

対応: `temit.f` (INTRA の converge ループ, tintraconv), `tintrb.f` (IBS 拡散と touckl/touckm/toucke の積算), `Packages/emit.n` (Touschek クラス, Spline)。

| API | 説明 |
|---|---|
| `touschek_table(lattice, momentum_eV, em, pbunch, nbunch=1.0, mincoup=1e-2, rfsw=True, radcod=False, charge=1.0, amass=ELMASS, verbose=True) -> dict` | `INTRA; emit` 相当。`em` は `emittance6d` の戻り (`closed_orbit`, `trf0` を使う)、`pbunch` は `PBUNCH` |
| 戻りのキー | `'table'` (SAD の `TouschekTable` と同じ入れ子リスト `[id, [dp(120), 1/τ(120)], [[nx, ny, nz], mx, my], toucke]`), `'EMITXC','EMITYC','EMITZC'` (IBS 込み), `'emittance_x', …` (最終 tecalc), `'iterations'`, `'converged'`, 生の表 `'touckl','touckm','toucke','tampl'` |
| `touschek_lifetime(table, nx, ny, nz, sige) -> float` | `TouschekLifetime[nx, ny, nz]` [s]。`nx`, `ny` は `None`/`math.inf` で INF。両方 INF なら `nz·sige` の運動量受容に対する寿命 |
| `TouschekLifetimeCalculator(table)` / `.lifetime(nx, ny, nz, sige)` | emit.n の Touschek クラス (スプラインを保持) |
| `touschek_f(x)`, `SadSpline(x, y)` | TouschekF、`Spline[]` の直訳 |

```python
from pysad.touschek import touschek_table, touschek_lifetime
tt = touschek_table(lat, 2.5e9, em, 1e10, mincoup=1e-2, rfsw=True, verbose=False)
tt['iterations'], tt['converged'], tt['EMITXC']
tau = touschek_lifetime(tt['table'], None, None, 10.0, em['momentum_spread'])   # 受容 10 σδ の寿命 [s]
```

## 14. 通しの使用例 (Python から直接)

下の例は `load_sad_string` の FODO リング (170 要素) で、各段を順に実行したもの (全体で数秒)。実ラティスでも同じ呼び方で動く (§3〜§7 の時間を参照)。

```python
import sys, math, numpy as np
sys.path.insert(0, r'<repo>')
from pysad.parser import load_sad_string
from pysad.emit6d import emittance6d
from pysad.match6 import Matcher
from pysad.track6n import track_particles_fast
from pysad.dapert import dynamic_aperture

SRC = '''
MOMENTUM = 2.5 GEV;
DRIFT  D1 =(L=0.3)  D2=(L=0.8);
BEND   B  =(L=1.5 ANGLE=0.22439947 E1=0.5 E2=0.5);
QUAD   QF =(L=0.3 K1=0.26)  QD=(L=0.3 K1=-0.24);
SEXT   SF =(L=0.2 K2=1.0)   SD=(L=0.2 K2=-1.5);
CAVI   RF =(VOLT=1.0E6 HARM=100 L=0);
MARK   MK =();
LINE CELL=(QF D1 SF D1 B D2 QD D2 SD D1 B D1)
     RING=(MK 14*CELL RF);
'''
lat = load_sad_string(SRC).build_lattice('RING')          # 1. ラティスを読む
tab = lat.twiss_table()                                    # 2. CELL; CALC → Twiss 表
print(tab[-1].mux / 2 / math.pi, tab[-1].muy / 2 / math.pi)   # 2.0779, 2.2566
em = emittance6d(lat, 2.5e9, rfsw=True)                    # 3. RFSW; emit
print(em['emittance_x'], em['momentum_spread'], em['trf0'])   # 1.770e-6, 9.93e-4, 0.0823
m = Matcher(lat, out=print)                                # 4. FREE QF QD; FIT NX 2.18 NY 1.84; GO
m.free('QF,QD'); m.fit('NX', 2.18); m.fit('NY', 1.84)
res = m.go(); print(res.r, m.values())                     # 2.5e-16, {'QF K1': 0.25215, 'QD K1': -0.20305}
x = track_particles_fast(lat.elements, np.array([[1e-3, 0, 0.5e-3, 0, 0, 0]]),   # 5. TrackParticles 100 周
                         100, 2.5e9, trf0=em['trf0'], rfsw=True, rad=False)
emx = math.sqrt(abs(em['emittance_x']) + abs(em['emittance_y'])); emz = math.sqrt(abs(em['emittance_z']))
score, result = dynamic_aperture(lat.elements, [[0, 40], [0, 20]], [-0.02, 0.0, 0.02], 50,   # 6. DA
                                 [0.0, 0.0, -math.pi / 2], 3, 1, 2.5e9, em, emx, emz, em['trf0'],
                                 rfsw=True, rad=False, out=sys.stdout.write)
print(score)                                               # 144
```

同じことを SAD スクリプトとして `run_sad_string` / `run_sad.py` に渡す場合は §11 の形 (`FFS;` 以降に `USE RING; CELL; CALC; FIT ...; FREE ...; GO; RFSW; emit; DynamicApertureSurvey[{{-.02,0,.02},{0,40},{0,20}}, 50, Output->6];`) で書く。

## 15. 旧 API (非推奨)

`ffs.py`, `matching.py`, `tracking.py`, `optics.py` (`FFSState`, `Match`, `track_particles`(旧), `dynamic_aperture`(旧), `periodic_optics` 等) と `parser.SADScript.ffs` / `run_ffs` は独自実装で SAD と未照合 (tracking.py は放射係数に疑い)。[`pysad/__init__.py`](../pysad/__init__.py) がこれらを export しているが、新規コードでは使わず、本書のモジュールを直接 import すること。`lattice.py` の `chromaticity()` / `summary()` / `emittance()` も同じ系統。

## 16. その他の API (§1〜15 に無いもの)

組み込み関数は `Interp.builtins` (名前 → f(引数のリスト)) に登録される。ここでは前節までに無い物を書く。
各項目の gs との一致の程度と限界は、対応する要素別の文書 (docs/13〜22) と docs/27 にある。ここの関数名と引数の形は `pysad/*.py` のソース (登録箇所とドキュメント文字列) で確認したもので、実行して値を確かめたものではない。

### 16.1 組み込み関数の登録場所 (一覧)

| モジュール | 登録される名前 (`B['...']`) |
|---|---|
| `sadscript.py` (216 個) | 数学・文字列・リスト・入出力・制御の基本関数、`Read`・`Skip`・`ReadString`・`ReadList`、`GaussianCoulomb`・`GaussianCoulombU`・`GaussianCoulombFitted`・`BBBrem1` ほか |
| `sadlist.py` (88 個) | `Append, Apply, Cases, Complement, Count, Delete, DeleteCases, Dimensions, Dot, Drop, Extract` などのリスト関数 |
| `ffscmd.py` (22 個) | `Twiss, LINE, Element, FFS, Emittance, TrackParticles, SeedRandom, GaussRandom, Random, DynamicApertureSurvey, TouschekLifetime, CalculateOptics, Master, RenumberElement, ResponseMatrix, ParabolaRandom, FFSCalc, BeamLineName, SynchroBetaEmittance1, Set` ほか |
| `sadelem.py` (9) | `Element, Element$, SetElement, FLAG, BeamLine, WriteBeamLine, ExtractBeamLine, PrintBeamLine, GetMAIN` |
| `specfn.py` (22) | `HypergeometricPFQ(Regularized), LegendreP, Fourier, InverseFourier, Multinomial, BitwiseAnd/Or/Xor, Complex, Union1` ほか (特殊関数は 458 式で gs と照合: docs/06) |
| `fitroot.py` (5) | `FindRoot, Fit, FitEmit, FitGaussian, NIntegrate` |
| `vecfn.py` (2) | `Vector, ReleaseVector` |
| `miscfn.py` (26) | `DownhillSimplex, StringToStream, TemporaryName, ListRandom, Fork, Short, Inequality, ToInputString, MessageName, SingularValues, MatrixBalance` ほか |
| `fileops.py` (35) | `OpenShared, ReadShared, WriteShared, Parallelize, FileNames, FileQ, CopyFile, RenameFile, DirectoryName, Environment` ほか |
| `optics_fn.py` (9) | `BeamMatrix, GeoBase, NormalCoordinates, OrbitGeo, SymplecticJ, SymplecticInverse, TransferMatrix, AccelerateParticles, SurvivedParticles` |
| `csr.py` (7) | `CSRInit, CSRMatrix, CSRConvert, CSRTrack, CSRHaissin, CSRSetupOY, CSROYMatrix` |
| `radfield.py` (2) | `RadiationField, RadiationSpectrum` |
| `corr.py` (1) | `LinearSolveConditioned` (CorrectOrbit・MakeBump などは SAD の `Packages/*.n` を `pkgload.py` 経由でそのまま実行) |
| `sadclass.py` (2) | `Class, At` (SAD のクラスの機構) |
| `tclbridge.py`・`tkcanvas.py` | `Tcl*`・`Tk*`・`X*` (GUI: docs/20)、`TkCanvas*` (描画: docs/21) |

(表の数は `B['名前'] = ...` の形で代入している箇所を数えたもの。`B.update(...)` や関数経由で登録される名前は数えていないので全数ではない。)

### 16.2 新しい関数

| 関数 | 形式 | 説明 |
|---|---|---|
| `SynchroBetaEmittance1` | `SynchroBetaEmittance1[nus0, nus1, dnus, mphi]` → `{{nus, emitx, emity, emitxp, emityp, conv}, …}` | K64 の `tfemits`/`temits` (シンクロ・ベータ結合の平衡エミッタンス)。`Packages/emit.n` の `SynchroBetaEmittance` が呼ぶ。`USE` 後に使う。実装は `synchrobeta.py` (`synchro_beta_emittance`)。gs との差・未確認は docs/27。24 セルの平面リング + CAVI で emitx が gs と 1.4e-7〜1.4e-6、結合のあるラティスは gs が SIGSEGV で比較不能 (pysad は NaN の ValueError で失敗) |
| `CSRInit`, `CSRMatrix`, `CSRConvert` | `CSRMatrix[pipe, {dx,dy,f0..f3,tilt}, omega, rho, center, mode]`、`CSRConvert[v または m, pipe]` | CSR の行列。mode 22・33・23・32・110・220・330・1000 (らせん断面 `csrhelmat`)。`CSRConvert` は K64 の `tfcsrvec` と同じ。docs/15 |
| `CSRTrack`, `CSRHaissin`, `CSRSetupOY`, `CSROYMatrix` | `CSRTrack[z, dp, zl, {dz, ndiv, inten, r56, r65, frf, phirf, damp, sigmae}, nturn]` ほか | 縦方向のマクロ粒子追跡 (z・dp を**その場で更新**する)、Haissinski 方程式、OY の方法。`Packages/CSR.n` の `Track` などが呼ぶ |
| `RadiationField`, `RadiationSpectrum` | `RadiationField[traj, obs]`、`RadiationSpectrum[{τ, f}, {λ1, λ2, dλ}]` | `traj = {{t..},{x..},{y..},{z..}}` (RADLIGHT の軌道) と観測点 `{x,y,z}` から放射場 (`radfield.py`、K64 の tlfield/tlcalc/tlspect)。`RadiationSpectrum` は NPARA の fork 並列なし (逐次) |
| `GaussianCoulomb`, `GaussianCoulombU`, `GaussianCoulombFitted` | `GaussianCoulombFitted[x, y, sx, sy]` (実数 4 つ、リストを返す) | 横方向ガウス分布のクーロン場。`Fitted` は K64 の `twspfu` (フィット式、`gausscoul.py`)。docs/06 |
| `BBBrem1` | `BBBrem1[e, ds/s]` (ds/s は 0〜1 の外でエラー) | 制動放射のモンテカルロ (`bbbrem.py`、docs/14) |
| `Skip` | `Skip[unit, type, n]` | 型 (Character・Real・Word など) の項目を n 個読み飛ばし、最後に読んだ項目を返す。EOF なら `EndOfFile`。n≤0 や Record などは gs と同じく未評価で返す |
| `Read` | `Read[unit, type]` | 型の読み方を gs に合わせた。`Real` で数でない語は式として評価した値を返す。`Word` は語の文字列 |
| `FindRoot`, `Fit`, `NIntegrate`, `FitGaussian`, `FitEmit` | | 数値微分で gs と 5e-10〜1e-9。ベクトルの独立変数は実装。複素 `FindRoot` は v2 モードで実装 (複素の Newton 法。k64 モードは TypeError。gs は壊れた値: B61)。docs/05・06 |
| `DownhillSimplex` ほか | | `optimize.n` の Amoeba を移植 (gs と 1e-13) |

注: `MessageName[f, "usage"]` の関数形は、gs は f を無視して `Undefined MessageName` を出す (B70)。v2 モードは f の本文を返す (`f::usage` と同じ。無ければ `General::tag`)。`f::usage` の形は一致。

### 16.3 既存の関数の拡張

| 関数 | 変更 |
|---|---|
| `TrackParticles` | (1) POL ON かつ RAD ON のとき、粒子ごとにスピンを追跡する (`track6.py` の `sprot`・`bsi`)。入力は 9〜11 列 `{x,px,y,py,z,dp, 角度, sy/\|s\|, flag}`、戻りの 7・8 列は `atan(sz,sx)` と `sy/\|s\|`。この場合は高速追跡 (track6n) でなく Python の追跡 (track6) を使う。k64 モードは損失の処理が無く、v2 モードは落とす。(2) RADLIGHT が ON のとき `{結果, 軌道}` を返す: 軌道は粒子ごとに `{t, x, y, z}` の 4 本の表 (`ffscmd._RadLight`)。k64 モードは K64 と同じで BEND の中を記録せず横方向の座標も入れない。v2 モードは 3 次元の軌道 (BEND の中も `RADLIGHTDIV` 点)。(3) TRPT かつ LWAKE/TWAKE のときウェイクを適用 (`wake.py`、docs/13)。(4) `LOSSMAP` が ON のとき戻りに落ちた要素・落ちた周の 2 行を足す (両モード。POL のときは 10・11 行目)。(5) v2 モードで `SPAC`・`PSPAC`・`WSPAC` が ON のとき空間電荷を入れる (`spch.make`。docs/10 §9.7)。(6) `RADTAPER` が ON のとき直前の Emittance[] の係数でテーパーを入れる |
| `Emittance[]` | POL のスピン量 (偏極時間・スピンチューン・偏極ベクトル。MULT・空洞・SOL (ドリフトのある領域) まで。`Polarization` は v2 モードで有限値)、SOL 領域を含む Twiss・Emittance[] (`solmap.py`、`solemit.py`。v2 モードは領域内の BEND・MULT・CAVI・ずれ・F1 も)、COORD・TCAVI (L=0)・UND (v2 は磁場つき)・PHSROT (v2)・MULT ANGLE≠0・ExternalMap の要素。`RADTAPER`・`EMITDIV` が効く |
| FIT の幾何の条件 (`GX GY GZ CHI1..3`) | DRIFT の L・BEND の ANGLE・L は K64 の `tdgeo` の解析微分 (`match6d.tdgeo_term`)、他の変数は中心差分 |
| `obj@name[...]` | クラスのメンバーでない名前は式のまま (K64 と同じ) |
| `ElectronMass` (スクリプト) | 510998.95 (K64 の値) |

### 16.4 `run_sad` と `run_sad_string` の現状

§11 の記述 (`run_sad(paths, out=None, cwd=None)`、`run_sad_string(src, out=None, cwd=None)`、どちらも `(FFS, Interp)` を返す) は 2026-11-16 のソース (`ffscmd.py`) で確認した。補足:

- `run_sad_string` は `FFS;` の行で MAIN 部と FFS 部に分ける。MAIN 部に `TRACK` 文があるときは、MAIN の命令 (ON/OFF, EMIT, 代入, TRACK, STOP) を書かれた順に実行する (`ffs.main_exec`)。
- 実行中のエラーは `I.continue_on_error = True` で、SAD と同じく表示して次の文へ進む。`SystemExit` は捕まえて戻る。
- 初期の `MOMENTUM` は 30 GeV (K64 の `initb1.f`)。
- SAD の `Packages/*.n` (CorrectOrbit, MakeBump, OptimizeOptics, `CSR.n` の MakeZL・SkinDepth など) は `pkgload.py` 経由で SADScript のまま読み込まれる (docs/06)。

### 16.5 限界・未確認

- 引数の検査 (引数の数・型のエラー文) は関数ごとに違うので、細部は各モジュールのソースを見ること。実行例は作っていない。
- `B['名前']` の数え上げは全数ではない (§16.1 の注)。
- 未翻訳: `wjfit/wfres/wfbin/wiord/wjdraw` (MAIN の TRACK の Fourier 解析。gs も SIGSEGV で参照が無い。チューンの抽出は `pysad.naff`: §16.7)。k64 モードでは MEASURE/TRCOD、SOL 領域内の MULT・BEND・CAVI、PHSROT、空間電荷、`Polarization` の有限値は動かない (v2 モードで実装。§16.6)。`measure.py` の `tlum/tluma`・MEASURE 系は gs と数値比較ができない (gs で MEASURE が動かない)。

### 16.6 互換モード `pysad.compat` と v2 のモジュール・API

**互換モード** ([`pysad/compat.py`](../pysad/compat.py))

| 関数・変数 | 内容 |
|---|---|
| `compat.mode()` | 現在のモード (`'v2'` か `'k64'`)。環境変数 `PYSAD_COMPAT` を毎回読む (`set_mode` で上書きしていなければ)。環境変数が `k64`・`v1`・`compat` なら k64、それ以外 (既定) は v2 |
| `compat.set_mode(m)` | `'v2'` または `'k64'` を設定する。`None` で環境変数に戻す。それ以外は ValueError |
| `compat.v2()`・`compat.k64()` | 真偽。コードは `if compat.v2(): <修正した動作> else: <K64 の動作>` と書く |
| SADScript の `Compat$V2[]` | 互換モードが v2 なら 1 (prelude の `PolynomialFit` が使う) |

- 既定は v2。`py -3.12 run_sad.py lattice.sad script.sad` も v2 で走るので、gs と数値を比べるときは `PYSAD_COMPAT=k64` を付ける。
- 既存の試験と回帰は k64 で走る ([`pysad/tests/script_tests/conftest.py`](../pysad/tests/script_tests/conftest.py)、`pysad/tests/regression_tests/run_regression.py`)。v2 の試験 (`test_v2_*.py`) はフィクスチャ `v2_mode` で明示する。
- 追跡のカーネルは `rs[1]` (v2 の印) を読むので、乱数状態 `rs` は必ず `track6n.make_rs()` で作る (`dapert.py`・`solmap.py`・`maintrack.py` も同じ)。
- 互換モードで切り替えない変更 (pysad 自身の不具合の修正、K64 に無い機能・gs と同じ結果が出る機能) の一覧は docs/28。

**v2 のモジュール** (`pysad/*.py`)

| モジュール | 内容 | 文書 |
|---|---|---|
| `compat.py` | 互換モード (上) | docs/28 |
| `phsrot.py` | PHSROT (v2 だけ)。`matrices(params, rad, fluc)` (M0・減衰・拡散)、`k_phsrot` (numba の追跡カーネル)、`annotate`、`step`、`emit` (Emittance[] の T・D・B) | docs/10 §9.7 |
| `radtaper.py` | RADTAPER。`solve(els, momentum_eV, dp0, fshift, charge, amass, tol, maxit)` が閉軌道の δ の分布と係数を外側の反復で求める。`scaled`・`factors`・`tapered` が磁石の強さに係数を掛ける | docs/10 §9.7 |
| `spch.py` | 空間電荷 (v2 だけ)。`make(ffs, tr, npart)` (`TrackParticles` から。SPAC・PSPAC・WSPAC)、`SpaceCharge` (追跡に挿入するフック)、`Solver3D`・`kick3d` (SPAC の 3 次元)、`Solver2D`・`kick2d` (PSPAC)、`WSpace`・`gauss_field`・`gauss_potential`・`make_wspace` (WSPAC)、`SpaceChargeError` | docs/10 §9.7 |
| `measure.py` | MEASURE・TRCOD・LOGMEASURE・ORBITJITTER・PUTMEASURE・GETMEASURE (v2 だけ)。`ttinit(mark, np0, rng, ...)` (初期分布)、`ttstat`・`ttstat_lines` (統計の表示)、`Tracker`、`cmd_measure`・`cmd_trcod`・`cmd_align`・`cmd_orbitjitter`・`cmd_putmeasure`・`cmd_getmeasure`、`trcod_fix`・`trcod_move`・`trcod_restore`・`trcod_reference` (基準の軌道)、`track_trpt` (MAIN の `TRACK … TURNS=-1`)、`close_measure_files`。`tlum`・`tluma` は従来のまま | docs/03 |

**v2 の主な API**

- `emit6d.emittance6d(lat, p, rfsw=True, radcod=True, pol=True)['spin']` に、v2 の平衡偏極度の量が入る: `equpol` (偏極度)、`depol_rate` (1 周あたりの脱偏極率)、`taudep` (脱偏極時間 [s])、`spin_spread` (平衡分布での n の広がり)、`spin_orbit` (入口の ∂n/∂u、3x6)、`equpol_k64` (K64 の式の参考値)。docs/16 §5.3。
- `miscfn.allow_system(True)`: `System[cmd]` の実行を許可する (環境変数 `PYSAD_ALLOW_SYSTEM=1` でも)。許可がなければ v2 は実行しない。docs/05。
- `track6.SolRegionError`: SOL 領域内に置けない要素 (v2)。
- 環境変数: `PYSAD_COMPAT`、`PYSAD_ALLOW_SYSTEM`、`PYSAD_SYSTEM_SHELL` (許可があるときのシェル。sh・bash の実行ファイルを書くと `<shell> -c cmd`)、`PYSAD_SAD_PACKAGES` (別の SAD の `Packages` のフォルダ。あれば最優先。無ければ pysad に同梱した SAD のコピー `pysad/sad_copy/Packages` を使うので、通常は要らない。GUI の接点はこの変数があるときも有効になる)、`PYSAD_GUI` (`1` で、環境変数なしでも同梱のコピーから Tcl/Tk の GUI の接点 (Tkinter.n・KBFrame.n・Canvas.n。Plot などを画面に描く) を有効にする。無いときは画面を開かず、`$DisplayFunction` は Identity)。
- FFS の変数・フラグ (v2): `SPACSTEP`・`SPACNX/NY/NZ`・`PSPACNX/NY/DX/DY/NZ/DZ/NTURNCALC`・`WSPACBEAM`・`WSPACRADIUS`・`PBUNCH`、`RADLIGHTDIV`、`EMITDIV`・`EMITDIVB`・`EMITDIVQ` (両モード)。
- `ffs.track_results` (MAIN の `TRACK … TURNS=-1` の戻り): entrance / exit の (重心, 分散行列)、final、alive (v2)。

**使えない範囲・未確認**: k64 モードでは上の v2 だけの機能は動かないか無視される (それぞれの文書の「k64 モード」の欄)。`spch.py` の電子以外の粒子の式は 100 MeV の陽子で確かめた (`TrackParticles` は `MASS`・`CHARGE` を追跡に渡す。SAD スクリプトからの空間電荷の使用は未確認)。上のモジュールの公開関数名はソースの `def`・`class` から書いたもので、引数の細部は各ファイルの docstring を見ること。

**同梱のコピーの未確認・限界**: 同梱の SAD のコピー (`pysad/sad_copy`) を元の SAD と比べたのは、3 つの例 (`fn_CorrectOrbit_Monitor_Steer`・`fn_CorrectLocal`・`fn_CSROYMatrix`) の出力の一致だけ。SAD の最新版との差は確認していない。GUI の接点を同梱のコピーだけ (環境変数なし、`PYSAD_GUI=1`) で動かした確認は、1 つの試験 (Tk のバージョンとフォント系の取得) だけ。

### 16.7 互換モードと切り替えないモジュール (`pysad/*.py`)

互換モード (`pysad.compat`) と切り替えない追加 (どちらのモードでも同じ関数)。使い方・範囲は docs/28。

| モジュール | 内容 | 文書 |
|---|---|---|
| `naff.py` | NAFF 風のチューン抽出 (Laskar の数値的な周波数解析)。周波数の単位は 1/周回、符号は exp(+2πiνn) の向き。`naff(z, nfreq=1, window=1, fmin=-0.5, fmax=0.5, pad=8, tol=1e-14)` (複素信号 z_n = Σ a_k exp(2πi ν_k n) の周波数 ν_k と複素振幅 a_k を大きい順に nfreq 個。Hann 窓つき、FFT のピークをゼロ詰めで探し、窓つきの内積の最大化 (Brent 法) で精密化し、見つけた成分を Gram-Schmidt で引いて次へ。振幅は窓つき最小二乗)、`naff_real(x, nfreq=1, window=1, remove_mean=True)` (実数の信号。解析信号にして `naff`。周波数は 0〜0.5)、`twiss_normalize(x, px, alpha, beta)`、`tune_xp(x, px, alpha=0.0, beta=1.0, window=1, nfreq=1)` ((x, px) の列と Twiss から z = x̂ − i p̂ を作り符号つきのチューン ν (−0.5〜0.5) を返す。SAD の追跡の座標は px = dx/ds なので x̂ = x/√β、p̂ = (α x + β px)/√β に直す)、`tunes_from_coords(coords, twiss=None, window=1)` (`track_particles` の座標列 (ターン × 6) から (νx, νy, νz))。確認: 既知の信号で周波数 1e-11、追跡の座標列 (1024 周) が CALC と 12 セル FODO 3e-9・大きなリング A 6e-10・PF 1.4e-10・ASTRID2 5e-11。MAIN の TRACK の `wtunem` (FFT の放物線補間) は 4e-5〜2e-4。限界: 周波数の間隔が約 2/N 未満の成分は分離できない、放射減衰のある列 (放射を切って使う)、共鳴線の近傍 | docs/28、[`example/05_tracking/fn_NAFF_tunes.py`](../example/05_tracking/fn_NAFF_tunes.py) |
| `warmup.py` | numba のカーネルの事前コンパイル (コマンド `pysad-warmup`、`py -3.12 -m pysad.warmup` でもよい)。小さな FODO を土台にした埋め込みの台本 10 本 (光学・エミッタンス・追跡・放射、BEND の K1・フリンジ、ステアリング、MULT・OCT・APERT、COORD・TCAVI、UND・PHSROT、空洞の横方向電圧、BEAMBEAM・WAKE、SOL 領域、TRPT の加速) を 1 回ずつ走らせ、カーネル 195 本のうち 191〜192 本の機械語をキャッシュに作る (残りの 3〜4 本は、この入力では呼ばれず、使うときに変換される。SOL の領域内の粒子が早く落ちるかどうかで 1 本変わる)。`--check` はキャッシュが有効か (CPU・numba の版・ソースが一致するか) を調べて終了コード 0/1、`--clear` は pysad のキャッシュだけを消す。キャッシュの場所は numba の規則 (`NUMBA_CACHE_DIR` があればその下、なければ `pysad/__pycache__`、書き込めなければユーザーのキャッシュ)。`pysad.nbcache` の依存先の記録も更新する。変換の時間は 1 台の PC (Python 3.12・numba 0.68) で 80〜90 秒 (キャッシュは 27.5 MB。作り直しの 2 回目は 4 秒)。実行時は `pysad-warmup` と同じ `NUMBA_CACHE_DIR` にすること。別の CPU・numba の版では使われず再変換される | docs/01 |
| `specfn_cx.py` | 複素の特殊関数 (mpmath が無いので自前)。`zeta_c(s)`・`hurwitz_c(s, a)` (Euler–Maclaurin)、`polylog_c(s, z)` (冪級数・μ = ln z の展開・倍角の公式の再帰)、`hyperu_c(a, b, z)` (積分表示・多項式・漸化式と Kummer の変換)。SADScript の `Zeta`・`HurwitzZeta`・`PolyLog`・`HypergeometricU` が複素の引数のときに呼ぶ。範囲は docs/06 | docs/06 |
| `specfn_cx2.py` | 複素の特殊関数の追加。scipy の複素対応 (`exp1`・`expi`・`digamma`・`sici`・`fresnel`・`lambertw`・`hyp1f1`・`hyp2f1`)、整数次の ExpIntegralE の漸化式・級数、`polygamma` (Hurwitz ゼータ経由)、複素の次数の BesselJ の級数 (\|x\| ≤ 25)、Whittaker W (`specfn_cx.hyperu_c`)。`specfn.install` が複素の引数で切り替える (実数の引数では切り替えない)。両モード |
| `spch_k64.py` | `K64Pspac` (k64 モードの PSPAC。`SpaceCharge._kick2d` から呼ぶ)。K64 の形状関数: 横は CIC、縦は「ほぼ最近接」、電位は 5 点差分 (境界 φ = 0) の厳密解と双 3 次補間、線密度は 3 次スプライン。`PSPACNZ < 3` はエラー |
| `nbcache.py` | numba のキャッシュの無効化。`py -3.12 -m pysad.nbcache [--force]` (`purge(dry_run=False, force=False)`、`stale_modules()`、`cached_modules()`)。numba の `cache=True` の索引 (`.nbi`) は関数を定義したファイル 1 つの内容のハッシュでしか無効にならないので、呼ばれる側 (別のファイルの njit 関数や、`emit6d.ELMASS` などの大域の定数) を編集しても呼ぶ側のキャッシュが古いまま残る。`__pycache__/nbdeps.json` にモジュールごとの依存先の内容のハッシュ (改行の違いは無視) を覚え、変わっていたら呼ぶ側の `.nbi`・`.nbc` を捨てる。最初の 1 回は記録だけ (最初から古いものがあるかもしれないときは `--force`)。**kernel (`elkern.py`・`perfkern*.py`・`turnkern.py`・`track6n*.py` など) を編集したら、試験の前に 1 回走らせる** (自動では呼ばない) | docs/28、`test_nbcache_purges_caller_when_callee_file_changes` |

**その他の主な API (1)**

- `spin.sideband_depolarization(...)`: 偏極の脱偏極率にスピンの歳差の位相の変調 (シンクロトロン側帯) を厳密に取った項を入れる。`emit6d.emittance6d(..., pol=True)['spin']` に `depol_rate` (側帯込み)・`depol_rate_linear` (線形)・`spin_spread_sb`・`sideband_error`。`spin.SIDEBAND = False` で線形の式だけに戻る。`spin.SOL_BPAR = False` で K64 の縦磁場の足し方 (B130) に戻る。docs/16。
- `track6n.bend_links` (tfgeo の inext・iprev の再現。B113)、`track6n.tent_all`・`texit_all`・`k_cav_q` (CAVI の DX・DY と V1 V20 V11 V02)、`track6n.edge_exact1` (K1≠0 の BEND の追跡の厳密な端。B117)。
- `emit6d._tbende_core`・`_tblfre`・`_tbendal` (BEND の F1・FRINGE の Emittance[])、`emit6d._tsteee`・`_tmulte` (ANGLE=0 の BEND・MULT の F1)。
- `match6d.sol_dcod` (領域内の四極の DX・DY の解析微分)、`sgeo._quad_step`・`solemit._tquade` (領域内の四極のずれ)。
- `ffscmd._mass_charge()` (電子と同じ値のときは何も渡さない。`TrackParticles`・Emittance[]・CALC6D に `MASS`・`CHARGE` を渡す)、`ffscmd.FFS._npara_chunks`・`_fn_track_npara` (`NPARA=n;` の fork の分け方)。
- `radfield._drop_repeated_points` (複数区間の `RadiationField`)。`fileops` の `Parallelize` (POSIX の fork。`PYSAD_PARALLELIZE=serial` で逐次)。
- 環境変数: `PYSAD_PARALLELIZE` (`serial` で `Parallelize` を逐次にする)。

**その他の主な API (2)**

- `compat.set_bend_k1_edge('face'|'radial'|None)`・`compat.bend_k1_edge()`・環境変数 `PYSAD_BEND_K1_EDGE`: K1≠0 の BEND の端の模型 (v2 だけ。k64 は常に `radial`)。実装は `lattice._face_kt`・`_face_kick_cod`・`_face_kick_jac`・`_bend_k1_cod_matrix`、`emit6d._tbende_core`、`track6n.k_bendi`・`track6._tbendi`。
- `compat.set_bend_k1_def('log'|'lin'|None)`・`compat.bend_k1_def()`・`compat.bend_k1_lin()`・環境変数 `PYSAD_BEND_K1_DEF` (B160): K1≠0 の BEND の K1 の定義。'log' (既定) は K64 (中央面の磁場が B0 + (K1/h) ln(r/ρ0)、(1 + hx)∂b/∂x が一様)、'lin' は B' 一様 (elegant・MAD と同じ)。v2 だけ。k64 は常に 'log' (docs/29 A59、docs/28)。
- `compat.bend_dk_fix()`・`compat.set_bend_dk_fix(True|False|None)` (B161): CALC・`Emittance[]` の K1≠0 の BEND の解析行列の dk を kick の導関数に直す修正。v2 は既定で True、k64 は常に False。`set_bend_dk_fix` は測定・試験用の強制 (None で既定に戻す)。
- `spin.ST_EXACT`・`KINETIC`・`STEP_RATE`・`CHROMATIC` (v2 の Polarization)、`CHROM_KMAX`・`CHROM_NEL`。`spin.sideband_depolarization(..., C6=...)`、`spin._chromatic_matrix`。`spin_results` の `sdamp_rel` (厳密な緩和率)・`sdamp_rel_b2` (A2 の Σ r (b·n)²)・`kinetic`、不安定なリングの P = NaN (`depol_error`)。`track6._spin_rot_cs` (BEND の ROTATE のスピン)。
- `synchrobeta.tesumb_v2`・`tefsetup(jump=)` (v2 の SynchroBeta)。`perfkern.tradke_db` の `exactd` (v2 の tradke の微分)。`radtaper` (v2 の BEND の係数 (δ_l + δ_l+1)/2)。
- `track6._tsteer_ens`・`_tblfri`・`_tmulk`・`_tmulti_ens` (F1・FB・FRINGE・AE の追跡)、`lattice._mult_params`・`_bend_as_mult` (MULT の AE)。F1 つきの MULT を含むリングは `track6n_ext.compile_mult` が NotImplementedError にして track6 に落とす。
- `dapert.dynamic_aperture` (trackd.f のスロット管理。DAPWIDTH の既定 7)。FFS の大域変数 `WSPACSC` (追跡の WSPAC のビーム行列)。`measure._k64_measure`・`_k64_trcod` (k64 の MEASURE・TRCOD)。`lattice._exch_phase` (入れ替わる INS の位相の整数部。v2)。`parser.Token.pos/end`。
- FFS: `NPARA=n;` (FLUC の追跡の子プロセスの分け方)。

### 16.8 周波数マップ・LMA・共鳴駆動項・ADTS のモジュール

K64 にない追加。互換モードによらず同じ関数 (docs/28)。解説は docs/26・docs/25。SAD の台本の関数は pysad だけ (gs には無い)。

| モジュール | 内容 | 文書 |
|---|---|---|
| `trksetup.py` | `TrackSetup.from_source(src, 'RING', rad=True, rfsw=True, fluc=False, alost=None)` (追跡の設定。`.track`、`.tunes_full()`、閉軌道の求め直し、RF バケット `bucket_zlimit`)、`find_periodicity(els, rows)` (セルの数・各セルの最初の要素・セルの長さ)、`select_locations(st, i0, i1, ...)`、`lower_priority` | docs/26 §2 |
| `freqmap.py` | `frequency_map(st, loc, x=, y=, dp=, nturn=1024, workers=1, norm='data')` → `FMapResult` (`D`・`nux1`・`nuy1`・`nux2`・`nuy2`・`lost_turn`・`x`・`y`・`dp`、`save`・`load_table`・`plot`)、`resonance_lines`・`resonance_check` (整数・半整数を含む)・`line_distance` | docs/26 §3 |
| `lma.py` | `momentum_acceptance(st, locs, nturn=None, nsynch=2, check_every=None, workers=2, ...) (`check_every` の既定 None は 1 = 要素ごとの損失の判定、0 は SAD の tapert のまま。`TrackSetup` が 0 でない値を持つときはその値)` → `LMAResult` (`s`・`dp_plus`・`dp_minus`・`touschek_input()`・`save`・`plot`)、`touschek_lifetime_local(touschek_table, result, setup=, period=, s0=)`、`acceptance_at`、`linear_momentum_limit` | docs/26 §4 |
| `q1cmd.py` | 台本の関数 `FrequencyMap[位置, nturn, X->, Y->, DP->, Output->, Plot->]` (既定は RFSW・RAD・FLUC OFF)、`MomentumAcceptance[位置, Cell->, Distance->, CheckEvery->, Verify->, TurnCheck->, Workers->, Output->, Plot->]` (`CheckEvery` の既定は 1)、`TouschekLifetimeLMA[Formula->]` (既定は Piwinski 型) | docs/26 §3.3・§4.3・§4.5 |
| `lmabucket.py`・`fmapcheck.py`・`fmapres.py`・`fmappred.py` | LMA の RF バケットとアパーチャの区別 (`limit_types`、`bucket_1d`)、周波数マップの信頼性 (`check_setup`・`check_result`・`D_norm`・`feeddown_estimate`)、RDT・ADTS から予想した共鳴領域との比較 (`ResonancePredictor`、`apply_alignment_errors`、`coverage_stats`)、予想の作り直し (`ResonancePredictorS3`、`evaluate_prediction`・`auc_score`) | docs/26 §3.5・§3.6・§4.5・§10、docs/25 §12 |
| `rdt2nd.py`・`nform.py` | 4 次の h_jklm (`rdt4`、`resonance_terms` に 4 次の名前を渡すと回る)、正規形の補正 (`nf_correction_h4`・`rdt4_normal_form`)、`adts_from_rdt4`、`fourth_order_adts`、数値の正規形 (`nform.adts_from_map`・`adts_at_delta`・`adts_delta_derivative(method='analytic'\|'nform')`) | docs/25 §10・§11 |
| `resdrive.py` | `resonance_terms(lat, terms=, method='analytic', nslice=8, delta=0, ref=None, rows=None, bend_cubic=True, bend_nslice=16, momentum_eV=3e9, order2=None, rdt4_options=None)` (引数の `None` は既定の動作。h_jklm、要素ごとの寄与。`bend_cubic=True` が既定で偏向磁石の 3 次の項の分 `.total_bend` を足す、`.total_sext` は 6 極のみ。δ≠0 では入れない)、`adts(lat, method=, include_octupole=True, delta=0, fourth=None, optics=None)` (6 極 2 次 + 8 極 1 次、`by_harmonic`・`per_element`。既定は物理的な項だけ = `fourth='r3'`・`optics='calc'`。`fourth='r6'` は四極の積分器の台形則と偏向磁石の 4 次の項 (`.bend4`) を足す精密版、`optics='tracked'` は δ≠0 の光学を追跡の閉軌道から作る。docs/25 §11)、`chromaticity_numeric(lat, order=3, dmax=2e-3)`・`chromaticity_sextupole(lat)`・`chromaticity_natural(lat)` (自然クロマティシティの式 −(1/4π)∮(K1+h²)β + 端)・`chromaticity_bend(lat)` (偏向磁石の 3 次の項の分散による分)、`slice_thin`・`adts_tracking` (検証用)、`kinematic_adts`・`quad_fringe_adts` (`adts().kin`・`.fringe`・`.total`)。`adts(lat, bend_cubic=True, bend_nslice=16)` (偏向磁石の 3 次の項: `.bend`・`.cross`・`.total`、`bend_cubic=False` で偏向磁石の 3 次の項を入れない値)、`bend_cubic_adts`。`bendnl.py`: `bend_decomposition(el, nslice)`・`bend_atoms`・`extract_map`・`Taylor2` (追跡の写像から偏向磁石の 3 次のリー生成関数を取り出す)。`rdttrack.py` (追跡から RDT を取り出す): `ring_map(elements, momentum, row0, amp=1e-5)`・`rdt_from_map(rm)` (1 周の写像の 2 次のテイラー展開から p=0 の h_jklm)・`delta_quadratic_total`・`delta2_orbit_h10002` (δ の項)・`tracked_rdt_table`・`periodic_points`・`predict_periodic_J` (3νx 共鳴線の近くの周期点)。`bendnl.py`: `harmonic_amplitudes_shift`・`bend_atoms_rdt` (h への偏向磁石の寄与、δ の項は分散による軌道のずれの展開)。`adtsdirect.py`: `direct_adts(lat, momentum_eV, jmax, delta=0, order=2, nturn=1024)` (追跡 + NAFF + 多項式フィットの ADTS)・`tune_vs_delta`・`fit_tune_surface`・`fit_delta_curve`・`paraxial_tracker`。台本 `ResonanceDrivingTerms[]` (入力 `BendCubic->`・`BendSlices->`)・`AmplitudeTuneShift[]` (`Kinematic->`・`QuadFringe->`・`Total->`。入力 `BendCubic->`・`BendSlices->`、出力 `BendCubic->`・`SextBendCross->`)・`ChromaticityTerms[]` (出力 `NaturalPart->`・`BendPart->`) | docs/25 |
| `sxtopt.py` | [試作] 6 極ファミリーの K2L の多項式モデル (`SxtModel`: RDT・ADTS の解析式を k の線形・二次形式で厳密に組み、`check(k)` で解析式の直接の値と比べる) と SLSQP の最適化。台本の関数なし。検証用で、LMA の予測には使えない ([`example/ddba/v07_14/sextupole/`](../example/ddba/v07_14/sextupole/)) | [`example/ddba/v07_14/sextupole/README.md`](../example/ddba/v07_14/sextupole/README.md) |

限界・未確認 (詳細は docs/26 §7・docs/25 §5・§11.8・docs/31): 周波数マップ・LMA の誤差つき・結合・FLUC は調べた条件だけ (docs/31)。k64 モードの LMA、PF リング、SP8 (大きなリング B) は未確認。共鳴駆動項・ADTS は結合・垂直分散・端部・運動学的な項・6 極×8 極・10 極以上を入れない。
