# 08. TrackParticles — SAD の tturn 系粒子追跡と pysad の track6.py / track6n.py

最終更新: 2026-11-06 (現状のスナップショット)。SAD ソースの行番号は 2026-10 時点。

> 照合対象: `<SAD>\src` (tftrack.f, tturn.f, tapert.f, tdrift.f, tquad.f, tsolqu.f, tbend.f, tbendi.f, tcav.f, tmulti.f, tsol.f, photons.f, drndsr.f, random_plugin_sad.c, tgauss_.c, tfRandom_.c) と `pysad/track6.py`, `pysad/track6n.py`, `pysad/sadrandom.py`。
> 追跡 (TrackParticles) は CALC/Emittance[] の経路 (tturne: tquade, tbende, tthine, tcave, tradke) とは**別のルーチン**を使う。pysad もそれに合わせ、光学計算 (lattice.py / emit6d.py) とは独立に track6.py を書いている。

---

## 1. TrackParticles[] の仕様 (tftrack.f)

```
TrackParticles[{loc, beam}, dest, nbegin, nend]
```

| 引数 | 内容 | 根拠 |
|---|---|---|
| `{loc, beam}` | `loc` は開始要素 (番号または名前)、`beam` は粒子リスト | tftrack.f:74-86 |
| `dest` | 終了要素 (省略時 nlat = ラティス末尾) | tftrack.f:64-73 |
| `nbegin` / `nend` | 開始・終了ターン番号。`mt = nend − nbegin + 1` 周 (既定は 1 周) | tftrack.f:40-61 |

- `beam` = `{{x…},{px…},{y…},{py…},{z…},{dp…},{flag…}[,{lost-pos…},{lost-turn…}]}`。POL OFF なら 7〜9 列、POL ON なら 9〜11 列 (phis, sz が加わる) でなければエラー (tftrack.f:106-116, tfsetparticles のコメント)。
- **px, py は px/p0** (基準運動量で割った値) で渡し、同じ形で返る。内部では別の座標を使う (§2)。
- `flag` が 0 の粒子は最初から損失扱い。`tfsetparticles` が `kptbl` (粒子 ID ↔ 配列位置の対応表、損失位置、損失ターン) を初期化し、**死んだ粒子を配列の末尾へ寄せる** (tftrack.f:480-)。生きている粒子数 `npa` だけが追跡に使われる。
- 戻り値は `{位置, 粒子リスト}`。位置は区間が順方向なら `dest`、逆向きなら `loc`。flag 列は生存 1 / 損失 0。LOSSMAP ON のときは損失要素番号と損失ターンの列が加わる。
- `NP` (FFS 変数) は粒子数の上限ではなく、`DynamicApertureSurvey` 等が使う既定の粒子数。TrackParticles は渡した粒子をそのまま追跡する。

### 1.1 周回ループ (tftrack.f:243-281)

1. `mt > 1` なら最初に `tturn(npa, ls, nlat, …)` (開始要素から 1 周の終わりまで)。
2. 残り `mt − 1` 周は `tturn(npa, 1, nlat, …)`。`npa ≤ 0` で打ち切り。
3. 最後に `tturn(npa, ls, ld, …)` で終了位置まで。`ld ≤ ls` なら `normal = .false.` にして省略。

`TrackParticles[{1, beam}, 1, 1, N]` は N 周 (終点 = 始点なので最後の部分周回なし)。RFSW OFF のときは各周の終わりで `z(1:np) = 0` (tturn.f:142-145)。

### 1.2 NPARA > 1 での fork (tftrack.f:159-196)

`NPARA` (Packages/init.n の既定は 4) が 1 より大きいと `tftrack` は粒子を分割して `itffork()` で子プロセスに配る。

```
npara = min(NPARA, npz/npparamin + 1, ne*npz/npnlatmin + 1)     npparamin = 9, npnlatmin = 3000
```

- **粒子が 9 個以上**で分岐が起こる (8 個以下なら `npz/9 + 1 = 1`)。
- 子プロセスごとに `kseed = kseed + 2` として `tfaddseed(kseed)` を呼び、`SeedRandom[現在の種 + kseed]` で乱数の種をずらす (tftrack.f:629-)。したがって FLUC の乱数列を SAD と pysad で一致させるには **`NPARA = 1` をスクリプトに入れる**必要がある (親子で別の乱数列になるため)。
- wake 有効時は並列化しない (tftrack.f:133)。

---

## 2. 内部座標と要素ループ (tturn.f)

### 2.1 tconvm — 入出力の変換 (tconvm.f)

| 向き | 変換 |
|---|---|
| 入力 (idir < 0) | `px ← px/(1+δ)`, `py ← py/(1+δ)`, `dv = −δ(2+δ)/h1/(h1 + (1+δ)h0) + dvfs` (h1 = √(1+p1²), p1 = p0(1+δ)) |
| 出力 (idir ≥ 0) | `px ← px(1+δ)`, `py ← py(1+δ)` |

内部の運動量は **px/p (全運動量で割った値)** で、`dv` は粒子ごとの状態変数 (速度差 1 − v/v0 に相当)。`dvfs` は `tsetdvfs` (tturne.f:859) が `merge(FSHIFT, 0, rfsw)` で決める定数 (RFSW ON のときだけ FSHIFT)。放射でエネルギーが変わると `dv` も更新される (§4)。

**dvfs が dv に足される場所 (追跡)** (詳細は docs/history/work_small_items.md E): tconvm/tconv/tinip の入力変換、tradkn (photons.f) と trad/tsdrad/tserad/tqrad の放射後、tcav の h1 (`p1r h0/(1−dv+dvfs)`) と dv、tcavfrin、tmulti の空洞部。tcav の `dvn` (dzn) は含まない。pysad は `_Phys.dvfs` (`track_particles(…, fshift=)`、FFS の TrackParticles / DA / MAIN の TRACK は FSHIFT を渡す)。numba 版は状態 S を (n, 11) にして最後の列に dvfs を持つ。bench_level29e (NORAD/RAD/FLUC、FSHIFT=1e-6) が座標 ≤1.9e-11。

### 2.2 tturn / tturn1 — 要素外側・粒子内側

`tturn(np, lb, le, x, px, y, py, z, g, dv, sx, sy, sz, kptbl, n, normal)` → `tturn1` が `do l = lbegin, lend` で**要素を外側**に回し、各要素ルーチンが粒子配列 `x(np)…` をまとめて処理する。pysad の `_tturn1(parts, …)` も同じ構造 (粒子ごとに要素を回すのではない)。これが乱数の消費順 (§4) とソレノイド内四極の癖 (§3.7) に効く。

要素の分岐 (tturn.f:286-552):

| 要素型 | ルーチン | 備考 |
|---|---|---|
| DRIFT | `tdrift_free` | DISKIN≠0 なら近軸式 |
| BEND | `tbend` | K0 に RANK×gauss を加える |
| QUAD | `tquad` | |
| SEXT/OCTU/DECA/DODECA | `tthin` (nord = 6/8/10/12) | |
| MULT | `tmulti1` → `tmulti` (segment あれば `tmultiseg`) | |
| SOL | `tsol` (ソレノイド領域全体) | |
| CAVI | `tcav` | VOLT+DVOLT, DPHI, APHI≠0 なら COD の z で位相補正 |
| MARK | 何もしない | |
| APERT/MAP/BEAM/INS/COORD/UND/WIG | 専用ルーチン (pysad 未対応) | |

### 2.3 tapert — 15 要素ごとの損失判定

- `tturn` で `la = 0` にして始めるので最初の要素の前に必ず `tapert`。以後、`select` を抜けた要素ごとに `la = la − 1` し、`la ≤ 0` で `tapert` を呼んで `la = 15` に戻す (tturn.f:231-240, 553)。要素ループの後にも 1 回 (tturn.f:576)。`exit` で抜ける MARK 等では減らない (tsol 内のループは毎要素減る)。pysad は `Prog.la_count` で同じ数え方 (DRIFT/MARK/MONI は la を減らさない) をする。
- `tapert` はまず `limittrack` (tturn.f:7-18): `limitnan(a, x)` = `merge(x, merge(a, sign(x,a), |a| ≤ x), isnan(a))` — NaN は +x に、|a| > x は ±x に丸める。x, y は `xlimit = 3·LOSSAMPL`、px, py は 0.9999、z は 1e10。
- 判定 (ax = ay = 0 の場合): `−alost < x < alost` かつ `−alost < y < alost` かつ `|z| ≤ zlost` を満たさない粒子を落とす。`alost = LOSSAMPL` (既定 1 m)、`zlost = LOSSDZ` (既定 100 m)。コメントにあるとおり NaN を落とすために「条件を否定した」形で書かれている。
- 落ちた粒子は `kptbl(i,4) = l_track`, `kptbl(i,5) = kturn` を記録し、**末尾の生きている粒子と入れ替えて `np` を減らす** (tapert.f: Sweaping)。粒子の並びが変わるので、以後の乱数の消費順や §3.7 の癖に影響する。

---

## 3. 要素ごとのルーチンと式の要点

### 3.1 ドリフト `tdrift_free` (tdrift.f:10-32)

`dpz = pxy2dpz(px,py) = √(1−px²−py²) − 1`, `al1 = L/(1+dpz)`:
`x += px·al1`, `y += py·al1`, `z += dpz·al1 − dv·L`。近軸近似なしの厳密ドリフト。

### 3.2 四極 `tquad` (tquad.f:1-116)

L = 0 → `tthin(nord=4)`、K1 = 0 → ドリフト (bz≠0 なら `tdrift_solenoid`)。それ以外:

1. `tsolrot` (DX, DY のずらし → 回転 theta2)。theta2 = ROT + `akang` で K1 < 0 なら π/2 を足す (焦点面の入れ替え)、`ak = sign(K1L, L)`。回転は chi3≠0 のときだけ行う (pysad の `rot1` も th2=0 のときは回さない。x=inf, y=inf のとき 0·inf=NaN になるのを避けるためで、落ちる粒子の乱数消費 (§7) に効く)。
2. 入口ハードエッジフリンジ `ttfrin(nord=4, +ak)` (DISFRIN=0 かつ mfring ≠ 2)。akk = ak/L/4 の x³, y³ 次の写像 (tquad.f:325-)。
3. 線形フリンジ (FRINGE=1 or 3): `a = f1in/p, ea = e^a, b = f2in/p` で `x ← ea·x + b·px`, `y ← y/ea − b·py`, `px ← px/ea`, `py ← py·ea`, z も補正。`f1 = −|akk·F1²|/24`, `f2 = |akk|·F2` (tsetfringep, tffs.f:952)。
4. 本体: 放射なし `tsolqu`、放射あり `tsolqur` (入口で f1in ≠ 0 なら `tradk(f1in)`)。
5. 出口の線形フリンジ (FRINGE=2 or 3)、`ttfrin(−ak)` (mfring ≠ 1)、放射ありで f1out ≠ 0 なら `tradk(f1out)`、`tsolrot` の戻し。

`tsolqu` (tsolqu.f): `eps = 0.2·EPS` (EPS=0 なら 0.2)、`ndiv = 1 + int(|L·hypot(ak, bz)|/eps)`、各スライスは「半ドリフト → 厚い四極の閉じた式 (c1, s1, ch2, sh2) → 半ドリフト」。
`tsolqur` (放射あり): `ndiv = max(ndiv, ndivrad(…))`、**スライスを外側、粒子を内側**に回し、各スライスの後に全粒子まとめて `tradk(np, …, aln, 0)`。

### 3.3 薄い／厚い多極 `tthin` (tquad.f:118-322)

TENT (DX/DY/ROT の入口変換) → `ttfrin(nord)` (DISFRIN=0; 六極 nord=6 は専用式、それ以外は一般次数) → 本体 → `ttfrin(−ak)` → 放射ありなら `tradk(L)` 1 回 → TEXIT。
本体 (L ≠ 0): drift(α₁L) → ½キック → drift(αL) → 補正キック (β, γ 係数) → drift(αL) → ½キック → drift(α₁L)、α = 1/√12, β = 1/6 − α/2, γ = 1/40 − 1/(24√3), α₁ = 1/2 − α の 4 次積分。キックは `kord = nord/2 − 1`, `aki = K/kord!/(1+δ)`, `cx = (x − iy)^kord`, `px −= aki·Re cx`, `py −= aki·Im cx`。K は積分値 (pysad の Sextupole.K2 は長さあたりなので K2·L を渡す)。

### 3.4 偏向 `tbend` (tbend.f:531-673) と `tbendi`

- ANGLE = 0: K1 = 0 → `tsteer`、phib = ANGLE + K0 = 0 → `tquad`、それ以外 → `tmulti`。
  - ROT の渡し方 (tbend.f:565-): tsteer / tquad / tmulti には **theta+dtheta (ROT + DROT)** を渡す (CALC の tsteee は theta だけ)。tsteer は L=0 → tthin(nord=2, ak=K0): `px −= K0/(1+δ)` (TENT/TEXIT)、K0=0 → ドリフト、L>0 → 入口エッジ tan(AE1) + 本体 (ndiv 分割、放射は各スライスの後に全粒子 tradk)。
    垂直ステアリング (`ROTATE=-90 DEG`) は TENT で 90° 回した座標系で水平キックになる。BEND の ROT は `Bend.ROTATE` (+ROTATION) で、track6 `_tsteer_ens`/`_tbend_ens`、track6n の K_STEER/K_STEER0/K_BEND が読む。
  - ANGLE=0, K1≠0 の追跡 (tbend.f:566-590): K0=0 → tquad (theta+dtheta)、K0≠0 → tmulti。pysad は `track6._bend_k1_proxy` が QUAD/MULT の代理を作る (track6n は K0=0 を K_QUAD)。SAD の tquad 呼び出しは theta2 に akang を含めないので K1<0 の BEND は |K1| の集束四極として追跡される (11_sad_quirks A38、pysad も再現)。K1≠0 かつ K0≠0 は SAD が CALC で落ちるので追跡の照合はない。
- K1 ≠ 0 → `tbendi` (tbendi.f): `tbshift`/`tbrot` → `ttfrin(4)` → FB1 のフリンジ → **粒子外側**で `ndiv = 1 + int(|ANGLE|/1e-3)` (放射ありは `ndivrad` と max) のスライスを `tbendiinit/tbendicorr/tbendibody` で回す。**`rhoe = rhob·(1+δ)` は入口で 1 回だけ決め、放射で δ が変わっても更新しない**。放射は分割ごとに `tradk1`、最後に全粒子 `tradk`。
- 純粋な偏向: `tbshift(DX, DY, 回転)` → DROT/CHI2 があれば `tbrot(alg = L/2, phig = ANGLE/2)` → `tbendcore` (入口エッジ psi1 → 扇形 → 出口エッジ psi2 を 1 ステップで写す閉じた式。フリンジ `dxfr = fb²/ρ/24`, `dyfr = fb/ρ²/6`) → 戻し。
- 放射あり `tbendr`: `ndiv = min(1024, ndivrad(phir, 0, 0, EPS))`、入口・出口のフリンジ領域 fb/2 は `iprev/inext = 0` のときだけ別スライス (n1 = −1, n2 = ndiv + 2)。各スライスで `tbendcore` → 全粒子 `tradk`。

### 3.5 空洞 `tcav` (tcav.f)

- `phis = 0 (TRPT) / phic (APHI≠0) / w·trf0`。`trf0` は Emittance[]/CALC が決める同期位相 (pysad は `track_particles(trf0=…)` で与える)。
- `offset` の切り替え: `rad .or. trpt .or. autophi` なら `offset = sin(dphis), offset1 = 0`、それ以外は `offset = sin(dphis) − sin(phis), offset1 = sin(phis)`。RAD OFF では平均のエネルギー損失ぶんのキックを消す (放射損失がないので同期位相のぶんだけ打ち消す)。
- `tlim = 1e4 (RFSW) / 0 (NORFSW)`。`t = clamp(−z·h1/p1, ±tlim)`、`dh = max(oneev − h1, V·(−2 sin(ph) cos(ph − dphis) + offset)·wsn)`、`z = −t·p2r·p0/h2 − dzn`。h1 は `dv ≤ 0.1` のとき `p1r·h0/(1 − dv + dvfs)` から取る。
- 長さのある空洞: `vnominal = 0` (リング) なら `ndiv = 1 + int(min(|w·L|, …))`、`tcavfrin` (DISFRIN=0, mfring ≥ 0) が入口・出口で x, y に比例した横キックとエネルギー変化。

### 3.6 多極 `tmulti` (tmulti.f:37-)

ANGLE ≠ 0 → `tmulta` (pysad 未対応)。それ以外: `tsolrot` → `ttfrins(2n+2, akr[n])` (dofr(n): n ≤ 2 は非零で常に真) → [`tsolqu` / `tsolqum` (放射なし) / `tsolqur` (放射あり) の半幅 → n ≥ 2 のキック]×ndiv → 半幅 → `ttfrins(−)` → `tradk` → 戻し。`akr[n] = K_n·cr1^(n+1)` (cr1 は theta2 の回転)。`ndiv0` は tpara が `ampmaxm = 0.05` から決める。**`tsolqum` の z の dv 項は `−0.25(ap − dv)·aln0` で `tsolqu` の `−dv·aln0` と違う** (SAD の不整合、pysad は `_tsolqu(qum=True)` で再現)。

### 3.7 ソレノイド領域 `tsol` (tsol.f:510-)

`tsol` は BOUND ≠ 0 の SOL から次の BOUND ≠ 0 の SOL までを内部ループで追跡する。`bzs = tfbzs(k)` = charge·(BZ+DBZ)/(p0/c)。
入口: `trots` (CHI/DX/DY/DZ の座標変換) → F1 ≠ 0 なら `tradks` → DISFRIN=0 なら `tsfrin(bzs)`。出口: `tsfrin(−bzs)` → `trots`。
領域内: DRIFT は `tdrift_solenoid` (`px += bzp·y/2` で回転してから戻す; 放射は粒子外側の `tradkf1/tradk1`)、BEND は双極キック付き `tdrift`、QUAD は `tquad(bz = bzs)` → `ttfrin` の bz 分岐と `tsolqu/tsolqur` の bz ≠ 0 分岐 (tzsetparams)、CAVI は `tmultiacc`。

**癖 (tsolqu.f:465, 519-531)**: `tsolqur` の bz ≠ 0 分岐では `bzp = bz/(1+δ_i)` が各スライスの粒子ループ内で粒子ごとに設定されるが、スライス後の最後の半ステップ (`do i = 1, np` で `xsincos(r·bzp)`) は `bzp` を計算し直さないので、**最後のスライスの最後の粒子の bzp が全粒子に使われる**。RAD ON のソレノイド内四極の結果は粒子の並びに依存する (差は δ_i − δ_last に比例)。pysad は要素外側・粒子内側の構造 (`_tsolqu_bz_ens` → `_tsolqu_bz_final(p, aln0, bzp_last)`) でこれを再現している (Level 9 で確認)。

---

## 4. 放射 `tradk` (photons.f, module kradlib)

```
tradk(np, …, al, phi0):  al ≠ 0 なら  rfluct ? tradkfn : tradkn
                          その後  pxr0 = px, pyr0 = py, zr0 = z   (基準点の更新)
```

| 項目 | 式 (tradkn: FLUC なし) | tradkfn (FLUC あり) |
|---|---|---|
| 曲がり角 | 基準点 (pxr0, pyr0) を phi0 だけ回した px0 と現在の px, py の外積 `(ppx, ppy, ppz)` から `theta = asinz(|pp|)` | 同じ (ただし `dpz0` を回転後に再計算) |
| 光子数 | `anp = anrad·h1·theta`, `anrad = 5/(2√3)·rcratio` | 同じ `anp` を `tdusrn` に渡す |
| 臨界エネルギー | `uc = cuc·h1³/p0·theta/al1`, `al1 = al − z + zr0`, `cuc = 1.5·rclassic/rcratio` | `rho = al1/theta`, `uc = cuc·h1³/p0/rho` |
| エネルギー損失 | `de = cave·anp·uc`, `cave = 8/(15√3)` | `de = dph·uc` (dph = Σ 光子エネルギー/uc) |
| δ の更新 | `dg = −de/(1 + 2de)`, `g = max(gmin, g + dg)`, `gmin = −0.999999` | 同じ |
| 角度・位置 | `ddpx = −½·dpx·dg`, `x += ddpx·al1/3` | `ddpx = −r1·dpx·dg`, `x += r2·ddpx·al1` |
| dv, z | `dv = −g(1+pr)/h2/(h2 + p2) + dvfs` (tconvm と違う形だが SAD のまま), `z ← z·p2/h2·h1/p` | 同じ |

- **基準点 (setr0)**: 各要素ルーチンが入口で `pxr0 = px, pyr0 = py, zr0 = z` を置き、`tradk` の最後に更新する。放射のキックは「基準点からの角度変化」から計算するので、スライスごとに `tradk` を呼ぶ構造と一体。pysad は `_Particle.setr0()` / track6n の S[:,7:10]。
- `ndivrad(phi, ak, bz, eps)` = `int(1 + h0/100·aka)` (TRPT なら ×10)、`aka = (|phi| + hypot(ak,bz)·0.02)/eps` (tffs.f:641)。
- `tdusrn` (drndsr.f): `p = tran()` を 1 回引き、`p ≤ pn` になるまで光子を 1 個ずつ加える (ポアソン分布の累積)。光子ごとに `drndsr()` (横谷の逆関数近似、tran 1 回) と `r = tran()` を引く。光子 0 個でも p は消費する。`pn > 0.9999` で打ち切り。
- **NaN の粒子**: `anp` が NaN だと `p ≤ pn` も `pn > cutoff` も偽なので `npmax = 10000` 回ループし、1 粒子で 20001 個の乱数を消費する。落ちる粒子の乱数消費 (§7) に関係する。
- `max(gmin, g + dg)` は SAD のコンパイル結果では g+dg が NaN のとき gmin を返す。pysad は `gfloor` でこれに合わせる (g=NaN のままだと次の要素の tdusrn が 20001 個の乱数を引いて SAD とずれる)。

---

## 5. 乱数 (random_plugin_sad.c, tgauss_.c, tfRandom_.c)

| 項目 | 内容 |
|---|---|
| 生成法 | 32 bit 符号付き合同法 `seed ← 48828125·seed` (int32 のオーバーフローで mod 2³²)。乗数 A = 5¹¹ |
| `tran()` / `Random[]` | RANDOM_UNIFORM0: `(seed + 2³¹)/2³²` ∈ [0, 1) |
| `tgauss()` / `GaussRandom[]` | Box–Muller: `φ = seed·π/2³¹`, `ρ = √(−2 log(seed/2³² + 0.5))`。`ρ cos φ` を返し `ρ sin φ` を `next` に保持 (`filled`)。`|値| > GCUT` は棄却 |
| `SeedRandom[n]` | `seed = |n/2|·2 + 1` (奇数化、負なら符号付き)、`filled = 0`, `next = 0`。起動時の初期値 **17** |
| `SeedRandom[]` の戻り | `{"SAD", seed, filled, next}` (内部状態を uint32 4 語で dump: seed, filled, next の下位/上位 → 数値 4 個 + 名前で 5 要素。tfRandom_.c:177-195)。pysad の `SadRandom.dump()` が対応し、FFS 層は `{"SAD", seed, filled, next, 0}` で表示 (seed は符号なし表示) |

`tfaddseed(kseed)` (NPARA > 1 の子プロセス) は `SeedRandom[]` の 2 番目の要素 (seed) に kseed を足して `SeedRandom[seed + kseed]` を呼ぶ。

---

## 6. pysad の実装

### 6.1 track6.py (純 Python、粒子オブジェクト)

`_Particle`: `c = (x, px/(1+δ), y, py/(1+δ), z, δ)`, `dv`, 基準点 `pxr0, pyr0, zr0`。`out()` で px/p0 に戻す。

| SAD | track6.py | 備考 |
|---|---|---|
| tconvm | `_Particle.__init__` / `out` | dv は `ph.dv` (dvfs を含む) |
| tturn1 | `_tturn1(parts, elements, trf0, rfsw, rad)` | 要素外側・粒子内側。SOL は `_tsol_region` に委譲 |
| tdrift_free | `_tdrift_free` | |
| tquad / ttfrin / tqlfr / tsolqu / tsolqur | `_tquad_ens`, `_ttfrin4`, `_tqlfr`, `_tsolqu(qum)`, `_tsolqur_ens` | `_tsetfringep` が F1/F2 の係数 |
| tthin / ttfrin(6, 一般) | `_tthin_sext`, `_tthin_gen`, `_ttfrin6`, `_ttfrin_gen`, `_ttfrin_n` | SEXT/OCT/DECA/DODECA |
| tbend / tbendcore / tbendr / tbshift / tbrot | `_tbend_ens`, `_tbend0_ens`, `_tbendcore`, `_tbshift`, `_tbrot`, `_tbendal` | |
| tbendi (bendib) | `_tbendi`, `_BendI`, `_tbendiinit/_tbendicorr/_tbendibody` | 粒子外側 |
| tcav / tcavfrin | `_tcav`, `_tcav_drift`, `_tcavfrin` | |
| tmulti / ttfrins | `_tmulti_ens`, `_ttfrins` | ANGLE=0 のみ |
| tsol / tsfrin / tdrift_solenoid / tsolqu(bz) | `_tsol_region`, `_tsfrin`, `_tdrift_solenoid`, `_ttfrin4_bz`, `_TZ`, `_tsolqu_bz_ens`, `_tquad_sol` | BOUND≠0 のみ |
| tradk / tradkn / tradkfn / tradk1 | `_tradk`, `_tradk_core` (FLUC 分岐) | 乱数は `RNG` (`SadRandom`) |
| tran / drndsr / tdusrn / SeedRandom | sadrandom.py `SadRandom.uniform0`, `gauss`, `drndsr`, `tdusrn` | |

`track_particles(elements, coords, nturn, momentum_eV, trf0, rfsw, rad, fluc, seed, …)` が `TrackParticles[{1, coords}, 1, 1, nturn]` 相当。進捗表示 (`progress=`) と非有限値での即時停止を持つ。放射あり要素は `*_ens` でスライス外側・粒子内側にして SAD と同じ順で乱数を引く。

### 6.2 track6n.py (numba)

- 粒子配列 `S(n, 11) = (x, px/(1+δ), y, py/(1+δ), z, δ, dv, pxr0, pyr0, zr0, dvfs)` に対する `@njit` カーネル (`k_drift`, `k_quad`, `k_thin`, `k_bend0`, `k_bendi`, `k_cav`, `k_steer0`, `k_apert`) を track6 と同じ演算順序で転記。`NJ = dict(cache=True, error_model='numpy', nogil=True)`: `error_model='numpy'` でゼロ除算や `sqrt(負)` が例外ではなく Fortran 風に NaN/inf になる。
- `Prog(elements, ph)`: 要素を `(種類コード, パラメータ組)` の表に前処理 (フリンジ係数、回転の cos/sin、BEND の w/sqwh なども事前計算)。`turn(prog, S, n, ph, trf0, rfsw, rad, fluc, rs, apert)` が 1 周。`apert = (kid, lost_turn, kturn, alost, zlost)` を渡すと最初の要素・15 要素ごと・周末に `k_apert` (limitnan + 判定 + 末尾の生きている粒子との入れ替え) を行い新しい n を返す。`k_apert` は損失の記録も持ち、MAIN の TRACK の損失報告 (§8) が使う。
- 乱数は `rs = np.array([seed], int64)` を持ち回り、`lcg` は `(48828125·seed) % 2³²` を符号付きに折り返す (int32 のオーバーフローの再現)。`tran`, `drndsr`, `tdusrn` も `@njit`。終了時に `_t6.RNG.seed` へ書き戻すので track6 と乱数の状態を共有する。
- `fmax/fmin` と `NAN_MODE`: gfortran の `MAX/MIN` (x86-64 の maxsd/minsd, `a > b ? a : b`) は NaN が混じると後ろの引数を返すが、Python の `max/min` は前の引数を返す。発散した粒子の挙動を合わせるために `NAN_MODE = 1` で Fortran 流にしている (`tcav` の `max/min`、`sqrt1` の `max(1e-200,·)` など)。ただし `tradk` の `max(gmin, g + dg)` だけは SAD のコンパイル結果が逆 (NaN で gmin) なので `gfloor` を使う (サイトごとに違うので、新しい `max/min` を足すときは SAD のデバッグ版の再同期テストで確かめる)。`track6n.TRACE` は 1 要素ごとのデバッグ用フック。
- 数学関数は mathfun の直訳: `sqrt1`/`pxy2dpz` (NaN の出方を含む)、`asinz` (`asinz(NaN) = −π/2`)、`norm2_3` (±inf の扱い)、`xsincos`/`xsincosh` (級数付き)。発散した粒子の値 (inf/NaN) の出方を SAD と揃えるためで、落ちる粒子の乱数消費 (§7) に効く。
- `track_particles_fast(...)`: `track_particles` と同じ仕様。`Prog` が `NotImplementedError` を出す要素 (MULT, SOL, L=0 四極, ANGLE=0 K1≠0 BEND など) があれば **track6 にフォールバック**する。
- 性能: track6 の約 17 µs/粒子・要素に対し track6n は **176 ns/粒子・要素** (約 100 倍)。500 粒子 × 10 周 (RAD FLUC, ESRF 3563 要素) 3.1 秒。DynamicApertureSurvey (153 粒子 × 100 周) は track6 で 6.5 分、track6n で数秒。

---

## 7. 照合結果

| ベンチ | 内容 | 結果 |
|---|---|---|
| Level 6 | 小リング 5 粒子 × 1/10/100 周、NORAD / RAD (NOFLUC) | ≤5e-11 (SAD の印字精度 1e-10) |
| Level 7 | K1≠0 BEND (tbendi)、誤差 DELX/DELY/DTHETA (tsolrot, TENT/TEXIT, tbshift, tbrot) | ≤5e-11 (z のみ 1.7e-10、trf0 の微差) |
| Level 8 | 四極 F1/FRINGE=3、長さのある空洞 L=0.5、MULT K2 | ≤5e-11 / ≤5.9e-10 / ≤1.3e-10 |
| Level 9 | ソレノイド領域 (BZ=2 T, D1+QD)、tsolqur の bzp の癖 | NORAD ≤5e-11、RAD ≤1.5e-10 |
| Level 10 | RAD+FLUC, SeedRandom[17] | ≤5e-11、追跡後の seed も一致 |
| Level 11 | MULT K0/K1/K2 + DX/DY/ROTATE (tsolqum の dv 項) | ≤5.9e-11 |
| Level 25 (a,b,c,d,e,f,h) | ROTATE: 垂直ステアリング (tsteer L>0 と L=0)、skew 四極 (L=0 の SQ は tthin(nord=4))、ROTATE 付き BEND (K1 あり)、六極〜十極・MULT、組み合わせ | 20 周 NORAD ≤7e-16、RAD NOFLUC ≤4.4e-13 |
| Level 12 | 実ラティス ESRF v03_12 (3563 要素、NLAT=3564) 5 粒子 × 1/20 周、NORAD/RAD/FLUC | ≤4.7e-11、乱数状態一致 (track6n も同一、track6 との差 4e-18) |
| Level 29e | FSHIFT=1e-6 の追跡 (ESRF、NORAD/RAD/FLUC、5 粒子) | 座標 ≤1.9e-11 |
| Level 16j | 51 粒子 × 20 周 NORAD (DA の初期分布) | 5e-16 |
| Level 16g/16k (NPARA=1) | 51/32 粒子 RAD+FLUC、損失なし | 最終 seed まで完全一致 |
| Level 16l_j22/25/27 | 単粒子 (損失あり) を TrackParticles で 1,2,4,10 周 | 12 件すべて最終 seed 一致 |
| Level 16a | DynamicApertureSurvey NORAD 153 粒子 × 100 周 | SAD の地図 153 文字が全て一致 |
| Level 20 small | DynamicApertureSurvey RAD+FLUC 153 粒子 × 100 周 (esrf リング) | 地図・Score 63 一致。最終 seed は不一致 (カオス、下記)。再同期テストは 30 周で seed 全一致 |

**落ちる粒子の乱数消費** (詳細は `docs/history/work_loss_rng.md`, `docs/history/work_loss_rng2.md`): 粒子が落ちる周は、発散した粒子が tapert で除かれるまで光子の乱数を引き続ける (θ=π/2, anp≈94 なら約 100 個、NaN なら 1 個; §4) ので、その間の消費数を SAD と揃える必要がある。pysad が SAD に合わせている点は、`sqrt1`/`pxy2dpz`/`asinz`/`norm2`/`xsincos` の inf/NaN の出方 (§6.2)、`tapert` の la の数え方 (§2.3)、`rot1` が th2=0 のとき回さないこと (§3.2)、tradkfn の `max(gmin, g+dg)` が NaN で gmin を返すこと (`gfloor`, §4)。
判定は**再同期テスト** (`<work-scratch>/r2_resync.py`) で行う: SAD のデバッグ版 (各要素の全粒子の状態と seed を印字) の状態から pysad の次の 1 要素を走らせて seed・座標を比べる。30 周 (小規模 DAsv の 153 粒子、93 030 記録) のすべてで seed が SAD と一致する。
**自由に走らせた最終 seed は一致しない**: 落ちる粒子は発散中に光子を引き続け、その有無が 1e-11 の入力差 (丸め誤差の増幅) で逆になる。
SAD 自身も初期座標の相対 1e-16 の摂動で 1 周後の seed が変わる (12 回中 2 回、1e-10 では 12 回中 12 回) ので、この基準での一致は期待できない。
判定は再同期テストの seed 一致で行う。

**未対応**: PHOTONS (tdusrnpl / tphrec の光子リスト)、POL (スピン sprot)、NPARA > 1 の fork、tmulta (ANGLE≠0 MULT)、BOUND=0 の SOL、SOL の DX/DY/DZ/CHI と F1、wake・空間電荷、APERT/MAP/BEAM/INS/COORD/UND/WIG 要素、LOSSMAP の出力列。

---

## 8. MAIN レベルの TRACK コマンド (詳細は docs/history/work_main_track.md)

SAD の MAIN の `TRACK USE=RING NPART= TURNS= NX=(..) NY=(..) NZ=(..);` (`ActTra.f` → `track.f` → `tracka.f`) は、emit の正規モード行列 r と閉軌道から σ 単位の規格化振幅で初期分布を作り、
`tturn` (TrackParticles と同じ追跡) を TURNS 周回して、周ごとの規格化座標から smear (2J の周平均とゆらぎ) と tune (FFT のピーク) を表示する。
`pysad/maintrack.py` に直訳 (`run_track`)。追跡は `track6n.turn` をそのまま使い、損失の報告 (`P. id lost in n-th turn at l(name), amplitudes: x y z`) は `k_apert` の記録から出す。

| 項目 | 内容 |
|---|---|
| 初期分布 | `emx=√\|EMITX+EMITY\|`、`emz=√\|EMITZ\|` (RFSW ON) か `\|SIGE\|`。物理座標 = r·(NX,0,NY,0,0,NZ)×(emx,emx,emx,emx,emz,emz) + codin (+TDXI/TDYI/TDZI) |
| smear | nsmear=`max(100,TURNS/50)/100*100`。n が nsmear の倍数か n=TURNS+1 のとき表 (nsm=`mod(n-1,nsmear)+1` > 1 のとき)。和は常に nsmear で割る |
| tune | zbuf (x+i px, y+i py, z+i dp、周の頭、長さ nts=TURNS+1 以上の 2 のべき) を `tcftr` (exp(+i…)) で FFT、パワーを (.5,1,.5) で平滑して極大を放物線補間、8 本 |
| K64 | MAIN の EMIT/TRACK は動かない (EMIT は何もせず TRACK は `tinip`→`tmap` で SIGSEGV)。参照は SADScript で同じ手順を書いた `bench_level27a.sad` / `27b` |
| 照合 | 27a: RFSW ON/RAD OFF、RAD ON/FLUC OFF、RAD ON/FLUC ON の 6 粒子 250 周。27b: 粒子が落ちるケース。smear・最終座標・tune・損失周が SAD と一致 (座標・tune ≤5e-11、smear・振幅 相対 ≤7e-9) |
| 未対応 | FOURIE ON の多項式フィット、`ResultsOfTracking` フック、TURNS<0 (trpt)、NPARA>1、PHSPLOTS の描画 (無視) |
