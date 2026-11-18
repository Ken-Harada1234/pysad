# 08. TrackParticles — SAD の tturn 系粒子追跡と pysad の track6.py / track6n.py

最終更新: 2026-11-18 (v2.1.0 の後) (現状のスナップショット)。SAD ソースの行番号は 2026-10 時点。

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
- 戻り値は `{位置, 粒子リスト}`。位置は区間が順方向なら `dest`、逆向きなら `loc`。flag 列は生存 1 / 損失 0。gs は LOSSMAP ON のとき損失要素番号と損失ターンの 2 列 (8・9 列目) を戻りに足す。pysad は**両モードで足す** (POL のときは 10・11 行目。生存は 0 と nend+1。9.7 の LOSSMAP。`example/05_tracking/flag_LOSSMAP`)。
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
- 子プロセスごとに `kseed = kseed + 2` として `tfaddseed(kseed)` を呼び、`SeedRandom[現在の種 + kseed]` で乱数の種をずらす (tftrack.f:629-)。したがって FLUC の乱数列を SAD と pysad で一致させるには **`NPARA = 1` をスクリプトに入れる**か、**`NPARA=n;` で同じ分け方を再現する** (v2.1.0。pysad の既定は NPARA=1)。子 k (k = 1, 2, …、条件 `ipn + np1 < np` かつ `npr < npara`) は乱数の種 = 親の種 + 2k で粒子 `[ipn, ipn + np1)` を追い、親は残りを元の種で追う (有効な並列数 = min(NPARA, np/9 + 1, ne·np/3000 + 1))。pysad の分け方は RAD + FLUC の追跡だけ (POL・RADLIGHT・PHOTONS・WAKE は分けない)。FODO リング (np = 200、1・10・100 周) の rms dp と最終の乱数の種が NPARA=1 と 4 の両方で gs と全桁一致、np = 16 (有効 npara = 2) の粒子ごとの dp も 12 桁一致、PF リング (760 要素、六極つき、np = 40) も一致 (k64 モード。v2 は放射の式の修正で rms が 6e-6〜1e-3 違うが、乱数の引き方・最終の種は同じ。例 `example/v2_next/C3_fluc_nparallel.sad`)。
- wake 有効時は並列化しない (tftrack.f:133)。wake は実装済み ([23_wake.md](23_wake.md))。

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
| COORD / TCAVI / UND / MAP / MULT (ANGLE≠0) / BEAMBEAM | 専用ルーチン | pysad は `track6n` の `K_COORD`・`K_TCAV`・`K_UND`・`K_MAP`・`K_MULTA` と BEAMBEAM で対応 (§9)。track6 (低速) は COORD・TCAVI・UND を扱わない |
| APERT | `tapert1` → `tapert` (§2.4) | pysad は `K_APERT` (numba) と `_apert_pass` (track6) で対応 |
| INS/BEAM/WIG | 専用ルーチン | INS は docs/history/work_unimpl.md (gs と一致を確認)。BEAM・WIG は未確認 |

### 2.3 tapert — 15 要素ごとの損失判定

- `tturn` で `la = 0` にして始めるので最初の要素の前に必ず `tapert`。以後、`select` を抜けた要素ごとに `la = la − 1` し、`la ≤ 0` で `tapert` を呼んで `la = 15` に戻す (tturn.f:231-240, 553)。要素ループの後にも 1 回 (tturn.f:576)。`exit` で抜ける MARK 等では減らない (tsol 内のループは毎要素減る)。pysad は `Prog.la_count` で同じ数え方 (DRIFT/MARK/MONI は la を減らさない) をする。
- `tapert` はまず `limittrack` (tturn.f:7-18): `limitnan(a, x)` = `merge(x, merge(a, sign(x,a), |a| ≤ x), isnan(a))` — NaN は +x に、|a| > x は ±x に丸める。x, y は `xlimit = 3·LOSSAMPL`、px, py は 0.9999、z は 1e10。
- 判定 (ax = ay = 0 の場合): `−alost < x < alost` かつ `−alost < y < alost` かつ `|z| ≤ zlost` を満たさない粒子を落とす。`alost = LOSSAMPL` (既定 1 m)、`zlost = LOSSDZ` (既定 100 m)。コメントにあるとおり NaN を落とすために「条件を否定した」形で書かれている。
- 落ちた粒子は `kptbl(i,4) = l_track`, `kptbl(i,5) = kturn` を記録し、**末尾の生きている粒子と入れ替えて `np` を減らす** (tapert.f: Sweaping)。粒子の並びが変わるので、以後の乱数の消費順や §3.7 の癖に影響する。**`l_track` はセッション最初の TrackParticles では 0** なので、その最初の tapert では最初から LOSSAMPL の外にいる粒子の印が 0 (生きている) になり落ちない (詳細は §9.5)。


### 2.4 APERT 要素 (tapert1)

- キーワードの既定は全部 0: `DX1 DX2 DY1 DY2` (矩形)、`AX AY` (楕円)、`DX DY` (中心のずれ)、`ROTATE` (回転)、`JDPX JDPY DP`、`COUPLE` (tapert は使わない)。光学では L=0 のドリフト、追跡だけで効く (SADHelp の APERT)。
- `tapert1` (tapert.f) は `xl=min(DX1,DX2)`, `xh=max(DX1,DX2)`, `yl, yh` も同様にして `tapert` を呼び、`tturn` の icAPRT の case は呼んだあと `la = 15` に戻して `exit` (la は減らさない。tturn.f:417-423)。
- 判定 (落ちない条件、tapert.f:78-174): `eli = AX≠0 かつ AY≠0`。`eli` で `xh==xl` または `yh==yl` なら矩形は ±1e100 (楕円だけ)。`ROTATE≠0` なら `xa = cosθ(x−DX) − sinθ(y−DY)`, `ya = sinθ(x−DX) + cosθ(y−DY)`、0 なら `xa=x−DX, ya=y−DY`。落ちない条件は `(xa/AX)²+(ya/AY)² ≤ 1` (eli のとき) かつ `xl < xa < xh` かつ `yl < ya < yh` かつ `|z| ≤ LOSSDZ`。**全キーワード 0 の APERT は矩形が空なので全粒子が落ちる**。LOSSAMPL の判定は同じ式で ±LOSSAMPL の矩形・楕円なし・ずれなし。
- `JDPX` か `JDPY` が 0 でないとき (tapert.f:32-60) は落とさず、`(x/xa)²+(y/ya)² > 1` (xa=min(|xh|,|xl|)、ya も同様。0 なら制約なし) の粒子を `tran_array` の 2 個の乱数で境界へ入れ直す: `x = xa cos φ1, y = ya sin φ1, px = JDPX cos φ2, py = JDPY sin φ2` (φ = 2π·乱数)、`g += DP`。`DP` はこの再入射のときだけ効く。dv は更新しない。
- 落ちた粒子: `kptbl(i,4)=l_track` (= APERT 自身の要素番号)、`kptbl(i,5)=kturn` を記録して §2.3 と同じ入れ替え。TrackParticles の戻りの flag は 0、座標は落ちた時点の値 (`limittrack` の頭打ちを含む)。MAIN の TRACK の報告 `P. id lost in n-th turn at l(APERT名), amplitudes: x y z` は、TRPT でなければ出し、TRPT のときは l_track が APERT のときだけ出す (tapert.f:182-187)。
- pysad: parser が APERT を `Drift(L=0)` に `_apert=True` と各キーワードの属性を付けて作る。`track6n.k_apert_gen` (numba、tapert.f 全体の直訳。LOSSAMPL の判定 `k_apert` も同じ関数)、`k_apert_el` (APERT 要素)、`Prog._compile` の `K_APERT`、`turn`。Python の追跡 (track6) は `_apert_pass` で同じ判定をし、落ちた粒子を末尾へ寄せる (`track_particles(alive=…, zlost=…)`)。track6 には LOSSAMPL の判定 (最初の要素・15 要素ごと) はない。
- gs との一致: 楕円・矩形・楕円+矩形・DX DY・ROTATE・全部 0・APERT 2 個・3 周・JDPX JDPY DP の flag と落ちた時点の座標 (`tests/test_apert.py`、3 GeV の FODO リング、RFSW/RAD なし、粒子 16 個)。
- 未確認・限界: ソレノイド領域内の APERT (tsol.f:789 は tapert1 を呼ぶ) は、k64 モードでは未対応で NotImplementedError、v2 モードでは領域内の APERT で粒子が落ちる (9.7。スピン追跡の経路だけは未対応)。RAD・FLUC つき、TRPT、MAIN の TRACK の報告は gs と比べていない。JDPX/JDPY の乱数は粒子が 11 個以下なら gs と一致、12 個以上では gs が別の乱数列を使うように見え一致しない (原因は未調査)。`DAPERT` (動的口径) の追跡は APERT 要素を通すが、gs と比べていない。

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
- K1 ≠ 0 → `tbendi` (tbendi.f): `tbshift`/`tbrot` → `ttfrin(4)` → FB1 のフリンジ → **粒子外側**で `ndiv = 1 + int(|ANGLE|/1e-3)` (放射ありは `ndivrad` と max) のスライスを `tbendiinit/tbendicorr/tbendibody` で回す。**`rhoe = rhob·(1+δ)` は入口で 1 回だけ決め、放射で δ が変わっても更新しない** (k64 モード。v2 モードはスライスごとの放射のあとに更新する: 9.7)。放射は分割ごとに `tradk1`、最後に全粒子 `tradk`。
- 純粋な偏向: `tbshift(DX, DY, 回転)` → DROT/CHI2 があれば `tbrot(alg = L/2, phig = ANGLE/2)` → `tbendcore` (入口エッジ psi1 → 扇形 → 出口エッジ psi2 を 1 ステップで写す閉じた式。フリンジ `dxfr = fb²/ρ/24`, `dyfr = fb/ρ²/6`) → 戻し。
- 放射あり `tbendr`: `ndiv = min(1024, ndivrad(phir, 0, 0, EPS))`、入口・出口のフリンジ領域 fb/2 は `iprev/inext = 0` のときだけ別スライス (n1 = −1, n2 = ndiv + 2)。各スライスで `tbendcore` → 全粒子 `tradk`。

### 3.5 空洞 `tcav` (tcav.f)

- `phis = 0 (TRPT) / phic (APHI≠0) / w·trf0`。`trf0` は Emittance[]/CALC が決める同期位相 (pysad は `track_particles(trf0=…)` で与える)。
- `offset` の切り替え: `rad .or. trpt .or. autophi` なら `offset = sin(dphis), offset1 = 0`、それ以外は `offset = sin(dphis) − sin(phis), offset1 = sin(phis)`。RAD OFF では平均のエネルギー損失ぶんのキックを消す (放射損失がないので同期位相のぶんだけ打ち消す)。
- `tlim = 1e4 (RFSW) / 0 (NORFSW)`。`t = clamp(−z·h1/p1, ±tlim)`、`dh = max(oneev − h1, V·(−2 sin(ph) cos(ph − dphis) + offset)·wsn)`、`z = −t·p2r·p0/h2 − dzn`。h1 は `dv ≤ 0.1` のとき `p1r·h0/(1 − dv + dvfs)` から取る。
- TRPT で VOLT≠0 の空洞は `vnominal ≠ 0` で基準運動量が上がる。`track6n.cav_full` が tcav.f の ndiv 分割と p0、h0 の更新を写している ([23_wake.md](23_wake.md) §2a)。
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
| 曲がり角 | 基準点 (pxr0, pyr0) を phi0 だけ回した px0 と現在の px, py の外積 `(ppx, ppy, ppz)` から `theta = asinz(\|pp\|)` | 同じ (ただし `dpz0` を回転後に再計算) |
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
| `tgauss()` / `GaussRandom[]` | Box–Muller: `φ = seed·π/2³¹`, `ρ = √(−2 log(seed/2³² + 0.5))`。`ρ cos φ` を返し `ρ sin φ` を `next` に保持 (`filled`)。`\|値\| > GCUT` は棄却 |
| `SeedRandom[n]` | `seed = \|n/2\|·2 + 1` (奇数化、負なら符号付き)、`filled = 0`, `next = 0`。起動時の初期値 **17** |
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
- `Prog(elements, ph)`: 要素を `(種類コード, パラメータ組)` の表に前処理 (フリンジ係数、回転の cos/sin、BEND の w/sqwh なども事前計算)。`turn(prog, S, n, ph, trf0, rfsw, rad, fluc, rs, apert)` が 1 周。`apert = (kid, lost_turn, kturn, alost, zlost)` を渡すと最初の要素・15 要素ごと・周末に `k_apert` (limitnan + 判定 + 末尾の生きている粒子との入れ替え)、APERT 要素の位置で `k_apert_el` を行い新しい n を返す。`k_apert` は損失の記録も持ち、MAIN の TRACK の損失報告 (§8) が使う。
- 乱数は `rs = np.array([seed], int64)` を持ち回り、`lcg` は `(48828125·seed) % 2³²` を符号付きに折り返す (int32 のオーバーフローの再現)。`tran`, `drndsr`, `tdusrn` も `@njit`。終了時に `_t6.RNG.seed` へ書き戻すので track6 と乱数の状態を共有する。
- `fmax/fmin` と `NAN_MODE`: gfortran の `MAX/MIN` (x86-64 の maxsd/minsd, `a > b ? a : b`) は NaN が混じると後ろの引数を返すが、Python の `max/min` は前の引数を返す。発散した粒子の挙動を合わせるために `NAN_MODE = 1` で Fortran 流にしている (`tcav` の `max/min`、`sqrt1` の `max(1e-200,·)` など)。ただし `tradk` の `max(gmin, g + dg)` だけは SAD のコンパイル結果が逆 (NaN で gmin) なので `gfloor` を使う (サイトごとに違うので、新しい `max/min` を足すときは SAD のデバッグ版の再同期テストで確かめる)。`track6n.TRACE` は 1 要素ごとのデバッグ用フック。
- 数学関数は mathfun の直訳: `sqrt1`/`pxy2dpz` (NaN の出方を含む)、`asinz` (`asinz(NaN) = −π/2`)、`norm2_3` (±inf の扱い)、`xsincos`/`xsincosh` (級数付き)。発散した粒子の値 (inf/NaN) の出方を SAD と揃えるためで、落ちる粒子の乱数消費 (§7) に効く。
- `track_particles_fast(...)`: `track_particles` と同じ仕様。MULT と SOL 領域は `track6n_ext` のカーネル (`K_MULT`、`K_SOL_FRIN`・`K_SOL_DRIFT`・`K_SOL_QUAD`)、COORD・TCAVI・UND・MAP・MULT (ANGLE≠0) は `K_COORD`・`K_TCAV`・`K_UND`・`K_MAP`・`K_MULTA` で処理する。
`Prog` が `NotImplementedError` を出す要素 (BOUND=0 の SOL、ANGLE=0 で K1≠0 かつ K0≠0 の BEND、L=0 の BEND、ANGLE=0 の BEND (tsteer) の F1/FB1/FB2、CAVI の AUTOPHI、PHSROT (K64 で動かない B38) など) があれば **track6 にフォールバック**する (track6 が対応しない COORD・TCAVI・UND・MAP は NotImplementedError のまま)。
- 性能: track6 の約 17 µs/粒子・要素に対し track6n は **176 ns/粒子・要素** (約 100 倍)。500 粒子 × 10 周 (RAD FLUC, 3563 要素) 3.1 秒。DynamicApertureSurvey (153 粒子 × 100 周) は track6 で 6.5 分、track6n で数秒。

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
| Level 12 | 実ラティス 試験用のラティス (3563 要素、NLAT=3564) 5 粒子 × 1/20 周、NORAD/RAD/FLUC | ≤4.7e-11、乱数状態一致 (track6n も同一、track6 との差 4e-18) |
| Level 29e | FSHIFT=1e-6 の追跡 (実ラティス、NORAD/RAD/FLUC、5 粒子) | 座標 ≤1.9e-11 |
| Level 16j | 51 粒子 × 20 周 NORAD (DA の初期分布) | 5e-16 |
| Level 16g/16k (NPARA=1) | 51/32 粒子 RAD+FLUC、損失なし | 最終 seed まで完全一致 |
| Level 16l_j22/25/27 | 単粒子 (損失あり) を TrackParticles で 1,2,4,10 周 | 12 件すべて最終 seed 一致 |
| Level 16a | DynamicApertureSurvey NORAD 153 粒子 × 100 周 | SAD の地図 153 文字が全て一致 |
| Level 20 small | DynamicApertureSurvey RAD+FLUC 153 粒子 × 100 周 (実ラティス) | 地図・Score 63 一致。最終 seed は不一致 (カオス、下記)。再同期テストは 30 周で seed 全一致 |

**落ちる粒子の乱数消費** (詳細は `docs/history/work_loss_rng.md`, `docs/history/work_loss_rng2.md`): 粒子が落ちる周は、発散した粒子が tapert で除かれるまで光子の乱数を引き続ける (θ=π/2, anp≈94 なら約 100 個、NaN なら 1 個; §4) ので、その間の消費数を SAD と揃える必要がある。pysad が SAD に合わせている点は、`sqrt1`/`pxy2dpz`/`asinz`/`norm2`/`xsincos` の inf/NaN の出方 (§6.2)、`tapert` の la の数え方 (§2.3)、`rot1` が th2=0 のとき回さないこと (§3.2)、tradkfn の `max(gmin, g+dg)` が NaN で gmin を返すこと (`gfloor`, §4)。
判定は**再同期テスト**で行う: SAD のデバッグ版 (各要素の全粒子の状態と seed を印字) の状態から pysad の次の 1 要素を走らせて seed・座標を比べる。30 周 (小規模 DynamicApertureSurvey の 153 粒子、93 030 記録) のすべてで seed が SAD と一致する。
**自由に走らせた最終 seed は一致しない**: 落ちる粒子は発散中に光子を引き続け、その有無が 1e-11 の入力差 (丸め誤差の増幅) で逆になる。
SAD 自身も初期座標の相対 1e-16 の摂動で 1 周後の seed が変わる (12 回中 2 回、1e-10 では 12 回中 12 回) ので、この基準での一致は期待できない。

**未対応** (両モード): SOL 領域内・COORD・TCAVI・MAP のあるラティスの光子 (PHOTONS の PhotonList は BEND・QUAD・MULT・SEXT・OCT・磁場つきドリフト・BEAMBEAM + BEND で対応: 9.7 (b))、BOUND=0 の SOL、wake のうち MULT の F1・FRINGE (23_wake.md)、BEAM/WIG 要素。k64 モードではさらに、SOL の DX/DY/DZ/CHI と F1、SOL 領域内の BEND・MULT・CAVI (§9)、空間電荷 (v2 モードは SPAC・PSPAC・WSPAC を入れる。9.7)。
(POL のスピン追跡、tmulta、COORD・TCAVI・UND・MAP、WAKE、BEAMBEAM は §9 で対応済み。LOSSMAP は両モードで対応。)

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
| 未対応 | FOURIE ON の多項式フィット、`ResultsOfTracking` フック、NPARA>1、PHSPLOTS の描画 (無視)。TURNS<0 (trpt) は k64 モードでは未対応 (リングの追跡として扱う)、v2 モードでは `TURNS=-1` が輸送路の 1 回通過 (`Statistics at Entrance` / `Exit`。docs/03 の「v2 モードの MEASURE・TRCOD」)。空間電荷・MEASURE 系は MAIN の TRACK (リング) では使えない |

- TRPT の放射の分割: gs の ndivrad (tffs.f) は TRPT で `int(1+h0/100*aka*10)` (そうでなければ 10 倍なし)。track6n は TRPT のとき ndivrad に渡す h0 を 10 倍にして合わせる。RAD FLUC のあとの乱数の状態が gs と一致する (pysad/tests/test_beambeam.py)。

---

## 9. 追加した要素・機能の追跡の範囲と既知の差

### 9.1 要素ごとの対応 (gs と照合した範囲は各文書)

| 要素・機能 | 追跡の実装 | 一致の範囲 | 文書 |
|---|---|---|---|
| WAKE | `track6n.turn`・`TrackRun.set_wake` (TRPT かつ LWAKE/TWAKE のときだけ有効)、CAVI の分割 `cav_full`、`mult_full` | docs/23 の範囲 (RAD ON は wake と別に 1e-6 程度ずれる記述あり。docs/23 §参照) | docs/23 |
| BEAMBEAM | `pysad/beambeam.py` (弱ビームが強ビームから z スライスごとにキック)。光学には効かない | docs/24 | docs/24 |
| COORD | `K_COORD`。`tturn` の icCOORD は `exit` なので tapert までの要素数 la を減らさない。`dv` は更新しない | 200 粒子 3 周で相対 0 | docs/25 |
| TCAVI | `K_TCAV` (`tcav.apply_tcav`)。放射は RAD と FLUC の両方が ON のときだけ。RANKICK・RANPHASE は NotImplementedError | 200 粒子 3 周、RAD FLUC まで相対 0 | docs/26 |
| UND | `K_UND`。K64 は磁場を設定できず磁場なしのドリフト (B37)。v2 モードは `FBX`・`FBY` つきの磁場のカーネル (9.7) | k64: 20 粒子 3 周で相対 0。v2: 解析式 (gs は比較できない) | docs/27 |
| MAP | `K_MAP` → `ffscmd._map_track` (`ExternalMap["TRACK", l, nt, x]`)。落ちた粒子もフラグ≠0 なら復活する (K64 の `doinject`) | 10 桁。MAIN の TRACK・DA・Touschek は MAP を恒等で通す | docs/30 |
| MULT (ANGLE≠0) | `K_MULTA` (`multa.k_multa`)。放射・FLUC・TRPT の ndivrad ×10 も BEND と同じ | 11 通り 20 粒子 3 周で相対 0 | docs/28 |
| SOL | `K_SOL_*`。BOUND≠0 の領域。K64 は SOL (BOUND≠0) の長さを周長 `circ` に入れない (`emit6d._emit_circ`、追跡も同じ) | Level 9 (§7) | docs/29 |

### 9.2 POL のスピン追跡 (`track6.py`、docs/36)

`POL` かつ `RAD` が ON のときだけ、`TrackParticles` が粒子ごとのスピン (`_Particle.s`、`bsi`) を `photons.f` の `sprot` (tradk の中) の直訳で追う。RAD OFF ではスピンは変わらない (K64 も同じ)。
POL ON・RAD ON の `fn_track` は高速追跡 (`TrackRun`) を使わず Python の追跡 (track6) に切り替える。**k64 モードでは LOSSAMPL の損失の処理 (§2.3 の tapert・flag) がなく粒子は落ちない** (APERT 要素では落ちる。§2.4)。v2 モードでは track6 の経路でも落とす (9.7)。
入力は 9〜11 列 `{x,px,y,py,z,dp, 角度, sy/|s|, flag}`、戻りの 7・8 列は `atan(sz, sx)` と `sy/|s|` (スピン = (st·cos 角度, sy, st·sin 角度)、st=√(1−sy²))。
`bsi` (磁場の縦成分の積分) は四極 (tsolqur の各スライスの両端)、BEND (tbendcore の入口・出口、K1 のある BEND は tbendi の各スライス)、MULT (tmulti)、tdrift の ak0 付き、SOL 内のドリフトと四極 (bz·al) に入る。四極の ROTATE は tsolrot の chi3 のスピン回転。
gs と 4e-12 以下 (NOFLUC・FLUC): 平面リング、ずれ・回転のあるリング、MULT (K1 K2 K3)、K1 の BEND、SEXT、薄い MULT、SOL 領域内のドリフト (`test_spin_track.py`)。
未確認: SOL (bz≠0) 内の四極のスピン (下の差 ③ のため)、空洞・MAP のスピン。v2 モードでは SOL 領域内 (BEND・MULT・CAVI・ずれ・F1)、UND、PHSROT のスピン追跡が通る (9.7)。

### 9.3 RADLIGHT の軌道 (`ffscmd._RadLight`、docs/36)

`TrackParticles` は RADLIGHT のとき RAD と TRPT を ON にし、`{{ls または ld, 座標}, 軌道}` を返す。軌道は粒子ごとの {t, x, y, z} の 4 本の表で、`RadiationField` の入力になる。
K64 の `tlinit`・`tlstore`・`tlresult` の直訳で、gs と差 0 (2 粒子・複数周・始点 ls≠1、`test_radlight.py`)。K64 の癖: 記録するのは各区間 (tturn) の始めと、L≠0 の BEND 以外 (型 < CAVI) の要素の後だけで、**BEND と曲げのある MULT の中は記録しない**。
(k64 モード) `tlinit` が geo0 の軸 1・2 を設定しないので**横方向の座標は軌道に入らない** (進行方向に沿った点だけが並ぶ)。位置は straight (dphi=0) の積み上げで、BEND で向きが変わらない。
そのため、B を ANGLE=0 の STEER にする流れでないと、意味のある場にならない。**v2 モードは軌道が実際の 3 次元の軌道になり、BEND の中も記録する** (9.7、docs/45)。

### 9.4 SOL 領域の追跡の限界

**k64 モード**: 領域内に置ける要素は DRIFT、QUAD (K1≠0)、L=0 の要素、SOL。BEND・MULT・CAVI は未対応 (NotImplementedError)。領域内の直線 BEND は K64 の `tsol.f` どおりに `tdrift(al, bzs, phiy, phix, …)` を呼ぶ直訳を試したが、
gs は偏向をほぼ無視し (BZ=0.5 T で ANGLE=0.1 でも位置・角度の変化が 1e-5)、pysad の直訳 (ANGLE·距離だけ曲がる) と合わないので取り下げた。原因は未確認 (K64 側の癖かもしれない)。
**v2 モード**: 領域内の BEND・MULT・CAVI・APERT と、ずれて傾いたソレノイド (`DX DY DZ CHI1..3`・`F1`) に対応した。解析式 (ローレンツ力の数値積分) で確かめたもので、gs との比較ではない (9.7、docs/44)。

### 9.5 損失の処理

- 損失の判定・入れ替えは `tapert` (§2.3) のとおり。高速追跡 (`track6n.turn` の `apert`) と MAIN の TRACK (§8) にある。track6 (Python) 側は APERT 要素の損失だけ扱う (§2.4。LOSSAMPL の判定はなく、POL のスピン追跡でも LOSSAMPL では落ちない)。
- **POL のスピン追跡の損失**: k64 モードは処理なし (9.2)。v2 モードは track6 の経路でも落とす (9.7)。TRPT で最初から LOSSAMPL の外の粒子の扱いは次項。
- **セッション最初の TrackParticles の `l_track=0`** (gs と照合済み。**k64 モードの動作。v2 モードは最初の tapert で落とす (B33)**、9.7): K64 の `tapert` は落とす印に `kptbl(i,4)=l_track` を書く。`l_track` は tmacro のセッション変数で、セッション最初の TrackParticles では 0 (CALC をしても変わらない)。このため**最初の tapert (tturn の最初の要素) では印が 0 = 「生きている」になり、最初から LOSSAMPL の外にいる粒子は落ちず**、座標は 3·LOSSAMPL に頭打ちされて追跡を続け、後の tapert で落ちる。2 回目の TrackParticles からは `l_track≠0` で最初の要素で落ちる。`TrackRun.l_track0` と `k_apert` の lcur==0 で再現し、gs (1 回目・2 回目) と、元の K64 (B33 で SIGSEGV) を直した SAD (LOSSAMPL=1e-4、3 粒子 x0={1e-3,1e-6,2e-3}、NORAD、TRPT、2 周、15 値) と全桁一致 (`test_small_items.py`)。元の K64 は落ちるので「正しい基準」は直した版。**B29**: QUAD に `RAD=1` を書くと COD・エミッタンスが 0 になるのはバグではなく仕様 (K64 は `krad = (RAD_QUAD == 0) .and. al /= 0`: tturn.f:360、tsol.f:677、tturne.f:572。要素の RAD が 0 以外ならその要素の放射を止める)。
- COORD は `tturn` の `exit` なので la を減らさない (9.1)。MAP は落とす粒子をフラグ 0 で返せ、フラグ≠0 なら復活する。

### 9.7 v2 モードで変わった追跡 (k64 モードは K64 のまま)

試験は `test_v2_C_*.py`・`test_v2_C2_*.py`、記録は `docs/history/v2_changes_C*.md`、例は `example/v2_C/`。根拠は解析式・独立な数値積分 (scipy の solve_ivp によるローレンツ力と Thomas-BMT 方程式の積分など)・保存則・極限。gs は K64 が正しく動く部分の参照にだけ使った。条件: 電子、3 GeV・周長 31〜106 m の FODO リング、1〜3 GeV の短い輸送路、ソレノイド 1〜3 T (bz = 0.3〜3 1/m)。**実機の大きなリングでは確かめていない。**

**(a) v2 モードで結果が変わるもの (K64 のバグ・物理の修正。k64 モードは K64 のまま)**

| 項目 | k64 モード | v2 モード | 根拠 |
|---|---|---|---|
| B33 セッション最初の TrackParticles | 最初の tapert は l_track = 0 のため、LOSSAMPL の外の粒子に落とす印が付かず追跡が続く (§9.5) | 最初の tapert で落ちる (座標は 3·LOSSAMPL に頭打ち、落ちた周 1、落ちた要素は始点の 1 つ前) | 「絞りの外から始まる粒子は落ちる」。2 回目以降の K64 の動作と同じ |
| B10 NaN の期待光子数 | anp が NaN だと 10000 回まわり乱数を 20001 個引く | anp が NaN・±Inf なら光子を出さず、乱数は 1 個 | 有限の anp では K64 と同じ乱数列 |
| B34 BEAMBEAM の傾き | cost = sin(傾き角)、sint = 傾き角 (回転にならない) | (cos t, sin t) の回転 | 回転に対する共変性 (t = 0.3 で 1e-9)。強ビームに x-y 結合 (R1..R4 など) があるときだけ結果が変わる |
| B37 UND の磁場 | 磁場なしのドリフト | `FBX`・`FBY` の磁場 (docs/27) | 平面アンジュレータの解析式 |
| B38 PHSROT | 行列がすべて 0 で粒子の座標が 0 | 意図された写像 M = N·Rot(2πPX, 2πPY, 2πPZ)·(1−D)·N⁻¹ + 拡散 (docs/36 の PHSROT) | M0 のシンプレクティック性 1e-13、固有値 = exp(±2πi P) |
| B35 SOL 領域内の未対応の要素 | SIGSEGV (pysad は NotImplementedError) | `SolRegionError` (K64 の tsol.f が出すはずの文 `Only DRIFT, BEND, … are supported within SOL`。長さ 0 の SEXT なども黙って通さない) | tsol.f の select case の case default |
| B69 SOL 領域内の BEND・MULT・CAVI・APERT | gs が領域内の BEND の偏向をほぼ無視。pysad は未対応 | docs/44 | ローレンツ力の数値積分と 1e-14〜1e-10 |
| SOL のずれ・回転・F1 | 未対応 | docs/44 | 実験室系の合成と 2e-12〜1e-8 |
| B81〜B83 スピン追跡 (sprot の縦磁場 bsi、tsolqur、刻みの角) | K64 の式 | SOL 領域・UND の中は縦磁場を運動量の向きに置く (N6)。領域内の四極の bsi と放射の基準 (N7)。歳差の角を厳密な θ に (N8。**すべての要素の追跡のスピンに効く**) | Thomas-BMT の積分と 1e-7〜1e-5 (K64 の式は 6e-4〜4e-3)、平面リングの 1 周の回転角 = 2π aγ |
| B84 wake の重み | 空洞が 2 つ以上に分割されると wake の蹴りの合計が (1.5n − 0.5)/n 倍 | 台形の重み (合計 1) | docs/23 |
| B85 AUTOPHI | TRPT でないとき位相の差が消える | 「位相の基準 = 設計軌道の到着時刻」 (docs/44) | gs (TRPT) と全桁。リングの追跡の AUTOPHI は K64 と違う |
| B14・B08・B39 (SOL 内の MULT・四極) | docs/07 | docs/07 | docs/07 |
| B111 tradke の反作用の角度の符号・RADTAPER の追跡の重み | tradke (Emittance[]・閉軌道) と tradkn (追跡) で符号が逆。Emittance[] の閉軌道から追跡すると 1 周で閉じない (x 4e-7 m) | 符号を反転。RADTAPER の追跡の BEND の重みを (3δ_l + δ_l+1)/4 に | 閉軌道の追跡の 1 周が 1e-9 m。v2 の Emittance[] (RADCOD) は gs と 1e-4〜1e-3 違う |
| B112 放射の古典半径 | `re0 = elradi` (電子の値。質量によらない) | `re0 = elradi·elmass/amass` | 30 GeV/c の陽子の U0 = 6.88e-4 eV (古典式) と 1e-6。電子は両モードで同じ |
| B117 K1≠0 の BEND の端 (tbendi) | 薄い写像 | `DISFRIN=0` かつ端の角 ≠ 0 で CALC と同じ厳密な tbedge | 厳密な写像 (ローレンツ力の積分) と 1e-10 |
| B123・B132 光子の偏光 ξ1・エネルギー列 | ξ1 が pp で 2 回割られ ±1 を超える。E 列に βγ を掛ける | pp で割らない。γ を掛ける | \|ξ1\| ≤ 1、k64 は gs と全列一致 (9.7 (b)) |
| B131 SOL 領域内の CAVI のスピン | 歳差しない (`krad = .false.`) | 放射・スピンの刻みに通す | Thomas-BMT と 1e-9 台 |
| SOL 領域内の DRIFT の刻み | 1 刻み (bz = 9 で放射が 1.7 倍外れ、スピンが 1e-2) | ndiv = 1 + int(\|bz·al\|/0.02) | Thomas-BMT と 2.3e-7〜8.9e-7、らせんの放射損失と 1e-4 |
| B83 (N8) スピンの刻みの角 (補足) | example/07 の 500 周のスピン角 0.8794 rad は刻みの誤差の積算 | 厳密な θ (1.1516 rad) | EPS=1E-5 で gs と v2 が 1e-8 (A2) |
| B76 K1≠0 の BEND の放射 | rhoe を入口の値のまま | rhoe = rhob·(1+δ) を更新 | 減衰の和 Jx+Jy+Je = 4 (docs/07) |
| B58 K1≠0 の BEND のスピン (tbendi) | 端の縦磁場に y/ρ がなく k·x·y を全スライスで足す | 入口で +(y/ρ + k x y)、出口で −(…) | docs/07 |
| 空間電荷 (SPAC・PSPAC・WSPAC) | 追跡に入らない | `TrackParticles` に入る (9.7 (c)) | 解析式 |
| RADLIGHT の軌道 (B65) | 横方向の座標が入らず、BEND の中を記録しない | 3 次元の軌道、BEND の中も記録 | docs/45 |

**(b) モードによらない追加・修正 (gs と一致する、または pysad 自身の不具合の修正)**

- **LOSSMAP**: `LOSSMAP;` で `TrackParticles` の戻りの flag の後に 2 行 (落ちた要素の番号、落ちた周) が付く (POL のときは 10・11 行目)。落ちた要素は判定 (tapert) のときの l_track なので、15 要素ごとの判定点か APERT の番号。gs の例 (`example/05_tracking/flag_LOSSMAP`) で要素 {0,0,242,242,242}、周 {31,31,4,1,1} と一致 (`test_lossmap_rows_match_gs`)。k64 モードでも足す。v2 は l_track = 0 で落ちた粒子にも落ちた周を返す。使えない: MAP 要素 (ExternalMap) が落とした粒子の要素番号 (0 になる)、空間電荷 (spch) が落とした粒子の要素番号。POL なしで高速追跡に乗らない要素の並びは `LOSSMAP: … 記録できない` の SADError。
- **スピン追跡の損失**: POL + RAD の TrackParticles も tapert で落ちる (最初の要素・15 要素ごと・区間の終わり。LOSSMAP の 2 行も付く)。gs: 落ちた粒子の要素 74・周 1・座標 .022037636775、生存粒子のスピンが一致 (k64 で 1e-12)。k64 モードの最初の呼び出しの B33 の癖 (復活) は近似 (gs は要素 60、pysad は 15)。
- **RADTAPER**: `RFSW; RAD; RADCOD; RADTAPER;` がすべて ON のとき、Emittance[]・EMIT が磁石の強さ (BEND の K0+ANGLE・K1、QUAD、SEXT〜DODECA、MULT) を閉軌道の局所の運動量に合わせる。続く TrackParticles も同じ係数 (直前の Emittance[] の結果) を使う。`Twiss["DDP"]` に閉軌道の δ。実装は K64 の tcod の反復の直訳ではなく、同じ不動点 (δ(1) = DP0、一様な係数 dptaper を未知数にして 6 次元で閉じる) を外側の反復で求める (`radtaper.py`)。gs と比べるときは gs 側で `CODCONV = 1e-4;` を入れる (既定の判定では軌道の残りが 1% 違う)。gs (CODCONV=1e-4): x −2.173197e-7 (pysad −2.17319317e-7、6 桁)、U0 10 桁、DDP 9 桁、チューン全桁、追跡 2 周 (FLUC) x 4.808978e-5 (pysad 4.808984e-5)。閉軌道の z は 0 になり位相の基準は `EquilibriumPosition` (trf0) に入る。`CELL; CALC6D; RFSW; RAD; RADCOD; RADTAPER; CALC;` でも効く。使えない: 4 次元の CALC のテーパー (K64 の qtwiss に式が無く、RFSW・RAD・RADCOD・RADTAPER ON の CELL CALC (4D) は RADTAPER なしの CALC と全桁同じ: gs で確認。実装しない。テーパーの光学は CALC6D で出る。k64 も CALC6D にテーパーを入れる: gs と 8e-8)、SOL 領域内の要素、RFSW OFF。見つからないときは `???-RADTAPER-closed orbit with taper not found` を出してテーパーなしで続ける。追跡 (RAD、NOFLUC) を Emittance[] の閉軌道から始めたとき、**k64 モードは 1 周で x が 4e-7 m、z が 2e-6 m 動く (K64 の tradke と tradkn で反作用の角度の符号が逆: B111。RADTAPER ありで x 2e-6 m、z 1e-5 m)。v2 は符号を直し、RADTAPER の追跡の BEND の重みも tturne と同じ (3δ_l + δ_l+1)/4 にして、x・px が 1e-9 m、z が 1e-8 m で閉じる** (3 GeV、周長 100〜600 m の電子リング、U0/E ~ 1e-4〜5e-4。小リング x 1.4e-9 m、大きなリング A (RADCOD) −4.5e-11 m、PF 2.7e-11 m、ASTRID2 (MPW 入り) −5.4e-9 m。20 周の追跡も閉軌道の周りに止まる)。v2 の Emittance[] (RADCOD) は gs と相対 1e-4〜1e-3 ずれる (k64 は全桁一致)。未確認: U0/E が大きいリング (8 GeV、周長 1.4 km)、SOL 領域内の tradke (bzh ≠ 0) と追跡 (`tsolqur` など) の整合、RADTAPER の重みのどちらが物理的か (「閉じる」を根拠にした選択)。試験 `test_radtaper_matches_gs`・`test_v2_C2_optics.py`。例 `example/v2_C/C4_loss_and_taper.sad`・`C2d_taper_optics.sad`。
- **MULT のスキュー成分** (SK0 SK1 SK2 …、ANGLE = 0): pysad は黙って無視していた (`lattice._mult_params` が `KsL` を読んでいなかった)。`akr(n) = cmplx(K_n, SK_n)·cr1^(n+1)` で CALC・Emittance[]・追跡のすべてに入る。gs: 追跡 4 成分が 2e-13、Twiss (BX BY R1 R2 DX DY) が 1e-9 (`test_mult_skew_components_match_gs`)。
- **ResponseMatrix の TRPT**: 応答に sqrt(gammab(stt+1)/gammab(mon)) を掛ける。gs と全桁 (`test_response_matrix_trpt_matches_gs`)。`CalculateOptics` の TRPT の加速、加速区間の分散の gr の向きの修正は docs/05。
- **`MASS`・`CHARGE`** (v2.1.0。pysad の不具合の修正で、両モード): `TrackParticles`・Emittance[]・CALC6D に渡る (以前は電子の質量・電荷 1 のまま)。FFS の中で与える (`CHARGE` は FFS に入るたびに 1 に戻る: gs と同じ)。運動学 (ドリフトの z = L (β1/β0 − 1)、1e-12)、ソレノイドのラーモア回転 (電荷 1・−1・2、5e-12 m)、Qs (30 GeV/c の陽子、1e-4) が解析式と合い、k64 は gs と全桁一致 (陽子の小リング、大きなリング A の 3 GeV/c の陽子)。放射の古典半径は v2 だけ質量に比例 (B112)。使えない・未確認: 放射の β の因子 (β < 0.99 の粒子の放射は U0 の 1/β² 倍)、TouschekLifetime・Intrabeam・RADINT の古典半径 (電子の値のまま)、陽子での確認は 3 GeV/c・30 GeV/c のリングだけ。例 `example/v2_next/A1_proton_ring.sad`。
- **光子 (PhotonList。RAD・FLUC・PHOTONS ON)** (v2.1.0): `tdusrnpl`・`radangle1`・`tphrec`・`tphotonlist` の直訳。行 = [E (eV)、gx gy gz、ux uy uz、xi1 xi2 xi3、粒子番号、要素番号]。対応: BEND (K1=0: `tbend` の刻みごと、K1≠0: `tbendi` の刻みごと + 最後。ROTATE・DX・DY・DROTATE)、QUAD、MULT (K0・K1・K2・DX・DY・ROTATE)、SEXT・OCT (DX・ROTATE)、磁場つきドリフト、BEAMBEAM + BEND (TRPT つきでも)。k64 は gs と全列一致 (光子 358 個: BEND と QUAD、325 個: MULT・QUAD・BEND、71 個: BEAMBEAM + BEND + TRPT と追跡後の乱数の状態。E 列は 1.4e-8 違う: K64 が γ でなく βγ を掛けるため。B132)。**v2 は ξ1 の誤りを直す** (B123: |ξ1| ≤ 1) し、E 列に γ を掛ける。解析式 (v2): 光子数 (5/(2√3)) α γ θ が統計の 4σ 以内 (24 セル FODO・3 GeV・100 粒子で −0.10%)、⟨u⟩/u_c = 8/(15√3) = 0.3079 (0.3059 ± 0.0023)、⟨u²⟩/u_c² = 11/27、スペクトルの KS 検定、位置が円弧の上 (2 mm 以内)、SEXT・OCT・MULT の光子数。BEND・MULT・SEXT・OCT・QUAD・BEAMBEAM があるラティスは Python の追跡 (`track6`) になり遅い (24 セルで 1 粒子 1 周 1 ms 程度)。使えない・未確認: **SOL 領域内の要素・ドリフト (磁場なし) の光子は記録しない** (乱数だけ K64 と同じに引く)、COORD・TCAVI・MAP があるラティス (従来の高速追跡で BEND があっても光子なし)、POL の追跡との併用、QUAD の光子の位置は gs と未照合 (BEND は全桁)、強い六極・八極で不安定になる条件。track6 の `_ndivrad` が TRPT の 10 倍を見ていなかったのを直した (POL・RADPOL・光子の追跡の TRPT の結果が gs に近づく)。例 `example/v2_next/B3_photons_bend.sad`・`D1_photons_mult.sad`・`D1_photons_bb_bend.sad`。
- **BEND の F1・FRINGE の RAD の追跡** (v2.1.0。pysad の不具合の修正。B113): K64 の tfgeo の inext・iprev を再現した (`track6n.bend_links`)。先頭が長さ 0 の MARK の輸送路の最後の BEND は自分自身に繋がり、F1 の放射の断片が無い。以前の pysad は断片を常に入れていた。gs と 1e-6 (8 桁表示): 1 個の BEND (INS) の 1 周 (δ = −2.17142054679e-5)、隣り合う BEND の 12 セルのリング、K1≠0 の BEND の対 (`tbendi`)、POL (Python の追跡)。大きなリング A・PF・ASTRID2 は F1 なしで一致。
- **CAVI の DX・DY・ROTATE と `V1 V20 V11 V02`** (v2.1.0。pysad の不具合の修正): 追跡の CAVI に K64 の TENT・TEXIT を入れ、パーサが V1・V20・V11・V02 を捨てていたのを直した。L≠0 の CAVI + DX・DY の差 1e-3 が解消 (gs と全桁。大きなリング A で空洞に DX = 2e-4・DY = −3e-4、QDA2 を 0.03 rad 回して RAD の 5 周追跡が差 0)。使えない: WakeFunction のある CAVI・SOL 領域内の CAVI の V1 V20 V11 V02 (NotImplementedError)。
- **EMIT のあとの `USE`**: `trf0` (RF の位相の基準) を `USE` で持ち越す (K64 の仕様。EMIT・Emittance[]・CALC6D のどれのあとでも gs と全桁)。
- **POL + TRPT + 加速** (v2.1.0。v2 だけ): SOL 領域の外の VOLT≠0 の CAVI (DX・DY・ROTATE も) を track6 で通す。断熱加速の解析式と 5e-7 (1 GeV、扇形の偏向磁石の間に 0・0.2・0.5 GeV)。k64 は従来どおり明示のエラー。使えない: SOL 領域内の加速、RF つき MULT・TCAVI の加速、リングの何周もの TRPT (K64 の TRPT は輸送路用)、スピン追跡の損失の処理。
- **PSPAC のメッシュの警告** (v2.1.0。v2 だけ): 横の rms が `PSPACDX`・`PSPACDY` の 1.5 セル未満 (キックが 10% 以上小さい)、縦の rms が `PSPACDZ` の 1.5 セル未満 (`PSPACNZ ≥ 4` のとき)、粒子の 1% 以上がメッシュの外のとき `???-PSPAC-Warning-…` を出す (結果は変えない)。メッシュが十分細かいとき (a ≥ 4 DX、一様なバンチの lb ≥ 10 DZ または `PSPACNZ<4`) は解析式と 1〜2%、gs と 1%。粗いメッシュは gs と違う (K64 は電荷を横 CIC・縦「ほぼ最近接」で置き、v2 は 2 次の B スプライン)。gs の x と y の 16% の非対称は再現しなかった。gs は PSPAC の状態を C の static 変数に持つので、1 セッションで TrackParticles を繰り返すと 2 回目から PSPAC が働かない (gs と比べるときは 1 セッション 1 回)。例 `example/v2_next/A1_pspac_resolution.sad`。
- **RADTAPER と `DP0`**: CELL CALC6D のあと `DP0` は先頭の `DDP` (閉軌道の δ) の別名になり、RADCOD の CALC6D のあとの 4D の CALC はそのエネルギーで光学を計算する (k64。gs と一致)。k64 の CALC6D も RADTAPER を入れる (gs と 8e-8)。

**(c) 空間電荷 (SPAC・PSPAC・WSPAC。v2 モードだけ。`pysad/spch.py`)**

k64 モードは今までどおり空間電荷を追跡に入れない (`test_spac_k64_mode_unchanged`、`test_mult_radius_k64_mode_unchanged`)。K64 の SPAC は B89 (力が 1/γ 倍に弱い)・B90 (SPAC のキックの後 px を (1+δ) でもう一度割る)・B92 (メッシュ) の誤りを持ち、それをビット単位で再現する翻訳は行っていない。PSPAC は式が正しく、v2 の値が gs と 3〜4 桁合うがビット単位ではない。gs での SPAC の確認 (1 m のドリフト、`RADIUS=0.02`、半径 1 mm の一様な円・長さ 0.1 m の一様なバンチ、2000 粒子): 5 MeV・1e11 個で gs の rms x 0.6226 mm に対し解析式 (包絡方程式 a'' = K/a、K = 2 r_e λ/(β²γ³)) 1.639 mm、v2 の pysad 1.623 mm。gs は K/γ の式 (0.645 mm) に近い (K64 の誤り B89)。

| フラグ | 入る場所 | モデル |
|---|---|---|
| `SPAC` | `DRIFT` (`RADIUS>0`)、SOL 領域内のドリフト、MULT (`RADIUS≠0`、ANGLE=0、VOLT=0、SOL 領域の外) | バンチの静止系の 3 次元の電位のキック。ドリフトを `SPACSTEP` (既定 0.05 m) ごとに割る。円筒の壁 (半径 RADIUS) は像の電荷。`RADIUS<0` は半径 \|RADIUS\| の絞りだけ。MULT は [本体の半分 → (多極のキック → 絞り → 空間電荷のキック) × ndiv → 本体の半分] で、真空ダクトの軸は磁石の軸 |
| `PSPAC` | `SPCH` 要素 | 横 2 次元のポアソン方程式 (PSPACNX×PSPACNY の矩形メッシュ、境界で電位 0) の電位と縦の線密度 λ(z) の積。前の SPCH からの距離の分のキックを SPCH の位置で入れる。`PSPACNZ<4` はバンチに合わせた 64 ノードを自動 (K64 の既定では gs が落ちる: B93。SPAC は粒子数が少ないと gs が落ちる: B98) |
| `WSPAC` | 各要素の入口 | ビーム行列のガウス分布の力 (強-弱)。ビーム行列は Emittance の平衡値、または `WSPACBEAM = {εx, εy, σz[, σδ]}`。`WSPACRADIUS` で縦の力の電位の基準。`SELFCOD` で力の中心を粒子の重心に |

設定変数: `PBUNCH` (既定 1e10)、`SPACSTEP`、`SPACNX/NY/NZ`、`PSPACNX/NY/DX/DY/NZ/DZ/NTURNCALC`、`WSPACBEAM`、`WSPACRADIUS` (詳細は `docs/history/v2_changes_C_spch.md`・`v2_changes_C2_spch.md`)。フラグは FFS の最上位に書く (`Do[...]` の中では切り替わらない)。粒子には z の広がりが要る。`RFSW` を ON にする (`NORFSW` は区間ごとに z を 0 にするので v2 は空間電荷をエラーにする。`RADIUS<0` の絞りだけなら `NORFSW` でも走る)。例 `example/05_tracking/flag_SPAC` (作り直した。以前の「gs でも SPAC の効果が出ない」は `NOSPAC;`・`SPAC;` を `Do[...]` の中に書いてフラグが切り替わっていなかったため)、`example/v2_C/C2s_mult_radius.sad`・`C2s_wspac.sad`。

- 確かめた範囲: 電子 5 MeV〜50 MeV、1e8〜1e12 個/バンチ、ビーム半径 0.1〜7 mm、真空ダクト 3 mm〜3 cm。SPAC・PSPAC・WSPAC のチューンシフトが Laslett の式と 1〜5% で一致 (SPAC・PSPAC はマクロ粒子の雑音)。gs との比較は SPAC の包絡 (上)、PSPAC (リング 3 周後の rms x/y が 1.776/1.042 mm と 1.777/1.042 mm、4 桁)、WSPAC (Δν が 1% 以内) の 4 例。
- 使えない・未確認: **追跡の WSPAC のビーム行列は空間電荷なし・MINCOUP なしの平衡値** (Emittance[]・EMIT・CALC6D・CALC (RFSW ON) には空間電荷が入る: docs/07。ビームサイズ `LINE["SIG…"]` は従来どおり)、ANGLE≠0・VOLT≠0 の MULT と SOL 領域内の MULT の SPAC (メッセージを出して空間電荷なしで追う)、QUAD・BEND・SEXT の中の SPAC (偏向磁石のあるリングでは SPAC は全周には入らない。全周なら PSPAC + SPCH か WSPAC)、複数バンチ (`NBUNCH`)、壁の抵抗、円でないダクト、MAIN の TRACK・`DynamicApertureSurvey`・スピン (`POL`) の追跡、ウェーク・放射との組み合わせ (走るが独立な根拠なし)、バンチの静止系で相対論的な運動 (γ·x' ≳ 1)。電子以外の粒子: 式は 100 MeV の陽子で 0.5% 確認した。`TrackParticles` は `MASS`・`CHARGE` を追跡に渡す (v2.1.0。SAD スクリプトからの空間電荷の使用は未確認)。マクロ粒子は 1000 個以上を勧める (チューンシフトを数 % で見るなら 4000 個以上)。

**(d) SOL 領域・UND・PHSROT の追跡のつなぎ**

- v2 の SOL 領域: 領域内の BEND (直線の座標のままの一様な双極磁場)・MULT・CAVI・APERT、ずれて傾いたソレノイド。使える範囲・使えない範囲は docs/44。
- UND・PHSROT は高速追跡 (track6n) と、POL + RAD のスピン追跡 (track6。UND は磁場つき、PHSROT は位相空間の写像だけでスピンは動かさない) の両方で通る。docs/27、36。
- TrackParticles 以外: MAIN の `TRACK`・DA は UND・PHSROT・領域内の要素を扱うが、確かめたのは `TrackParticles` と Emittance[]・CALC。

### 9.8 既知の差 (gs と合わない。いずれも `docs/history/work_unimpl.md`、スピンなしでも出る・未調査)

1. **SOL (BZ≠0) 領域内の四極 + RAD の追跡**: 保存済みの参照 `bench_level9.sad` (gs の出力、2026-10-14) は pysad と一致するが、現在の gs (2026-10-19 ビルド) は同じ入力で別の値を返す (`S R T1 P1` の x が −.0027867898 から −.0014539780)。NORAD は両方で一致。
   **保存済みの参照と現在の gs のどちらが正しいかは断定できないが、2026-11-16 夜の切り分け (単一粒子の 1 周追跡、2.5 GeV、RAD NOFLUC、emit のあと) で、pysad ではなく現在の gs が誤っている疑いが強くなった**。SOL の BZ を 0 → 0.01 T に変えると:
   (1) 領域内がドリフトだけ (`SB D1 SE`): gs = pysad (BZ=0・0.01 とも 9 桁以上で一致)。
   (2) 領域内が四極 (`SB QD SE`): BZ=0 は gs = pysad。BZ=0.01 で gs は x=6.33e-4、y=−5.6e-16 (横の結合なし)、pysad は x=4.04e-4、y=1.1e-8 (BZ に連続)。
   (3) `SB D1 QD SE`、BZ=0.01: gs の結果が四極なしの `SB D1 SE` と 8 桁で一致する (x=1.194038e-4、px=−5.363477e-4)。**gs は RAD ON かつ BZ≠0 の領域内で四極の収束・発散力を入れていないように振る舞う**。BZ=0.01 と 0.2 でも gs の x はほぼ同じ (2.95e-5、BZ=0 の −2.23e-4 とは大きく違う) ので BZ に連続でない。
   (4) 現在の gs の 2 つのビルド (10-19 の Linux6 と、11-14 に -O0 で作り直した `~/sad_dbg`) は同じ値 (x=2.954816161E-5、10 桁一致) で、最適化の違いではなく現在の K64 ソースから作った gs の振る舞い。保存済みの参照 (10-14) を出した 10-15 の古い実行ファイル (`obj/AMD64-Linux5/sad1.exe`) は symbol lookup error で今は動かず、同じ版かは確認できない。
   結論: pysad が壊れているのではなく、現在の gs が SOL 領域内の四極 + RAD を物理的に誤って扱う疑いが強い (pysad は物理と保存済みの参照に合う)。**K64 ソースのどこが原因かは未特定**。§7 の Level 9 の RAD の行 (≤1.5e-10) は保存済みの参照との比較である。
