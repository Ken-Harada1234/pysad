# 11. SAD の癖・仕様・ビルド依存の挙動 (pysad での扱い)

> 対象: K64 版 SAD (github KatsOide/SAD、WSL Ubuntu で `bin/gs` 実行)。ソースは `<SAD>\src` と `Packages/`。
> 行番号は 2026-10 の照合時点。「ベンチ」は `<repo>\bench_level*.sad/.result` と `*_pysad.py/.result`。
> pysad は「再導出せず直訳」が方針なので、下記の癖は (バグに見えるものも) そのまま再現している。

## A. 数値計算の癖 (結果の値に効くもの)

| # | 現象 | 根拠 (ソース / ベンチ) | pysad での扱い |
|---|---|---|---|
| A1 | **tsolqur (RAD ON, bz≠0) は最後の粒子の bzp を全粒子に使う。** 放射あり・ソレノイド内の四極は最後の半ステップで `tz%tz0%bzp` (最後に tzsetparams を呼んだ粒子の値) を全粒子に使うため、結果が粒子の並び順に依存する (差 ∝ δ_i − δ_last)。放射なしの tsolqu は粒子ループ内で完結するので起きない | tsolqu.f tsolqur (tzsetparams :268、bzp :273-)。Level 9 で同じ座標の粒子セットを変えると SAD 自身が 1.15e-6 変わるのを確認 | track6 を要素外側・粒子内側 (`_tturn1(parts,…)`) に組み直し、`_tsolqu_bz_final` に最後の粒子の bzp を渡して再現。SAD との比較は同じ粒子セットで行う |
| A2 | **tsolqum (放射なし MULT) の z の dv 項が tsolqu と違う。** tsolqu は `−0.25·(… + ap·aln0 + …) − dv·aln0`、tsolqum は `−0.25·(… + (ap − dv)·aln0 + …)` で dv に 0.25 が掛かる。δ=1e-3 で 1 周 5e-10 の z のずれ | tsolqu.f:86-88 vs :250。Level 11 (観測 6.7e-10) | `_tsolqu(qum=True)` で SAD のまま |
| A3 | **Cavity PHI の既定値は 0** (π/2 ではない)。同期位相は trf0 として自動で決まる (RADCOD OFF: `trf0=asin(U0/Vc)/w`) | initb1.f ky_PHI_CAVI=4、tturne.f:47-57、tcod.f。Level 4・12 | `elements.Cavity.PHI=0.0` に修正 (旧値 π/2)。trf0 は emit6d で決める |
| A4 | **4 次元 CALC では空洞は RFSW に依らずドリフト。** qcav が `rfsw = rfsw .and. trpt` としてから tcave を呼ぶので、リング (CELL) の CALC では RFSW ON でも空洞は長さ L のドリフト。追跡 (tturn) では CAVI は RFSW に関係なく tcav が呼ばれ、RFSW OFF は周回末に z=0 にするだけ (TCAV だけがドリフト) | qcav.f:19-20、tturn.f:142-145, 443-486 (bench_level14a_dbg) | `lattice.RFSW` 既定 False、CALC の空洞はドリフト。追跡 (track6) は tcav |
| A5 | **tbendi (K1≠0 BEND) の rhoe は入口で 1 回だけ決め、放射で δ が変わっても更新しない** (module bendib の変数) | tbendi.f / tbend.f tbendiinit。Level 7 (track6n 転記時に発覚) | 同じく入口で固定 |
| A6 | **薄い多極 (tthine) の中のドリフトにも dv 項** (`z −= L·dv(δ)`, tdrife0 dvon=.true.)。抜けていると TransferMatrices の (5,6) が六極 1 台ごとに 1.3e-9 ずれ、νs/α_c が 3e-6 ずれる | tthine.f、bench_level16b/16e | lattice の `_tthine_cod` に dv 項を追加 |
| A7 | **E1/E2 は常に ANGLE に対する比** (`psi1=E1*ANGLE+AE1`)。E1=0.5 で ANGLE/2、E1=1 で ANGLE。逆向き要素では E1↔E2 を入れ替える | tffs.f:1261-1262、tbende.f:448、INTGRL.f:91 | lattice/emit6d/track6 の Bend.E1 は rad (= E1×ANGLE+AE1)。parser.py:349-354 が SAD と同じ `E1 = E1*angle + AE1`, `E2 = E2*angle + AE2` に変換する。逆向き通過は parser.py:385-386 で E1↔E2 を入れ替える |
| A8 | **K1/K2/K3 は積分値** (QUAD K1 = K1·L [1/m]、SEXT K2 も積分値)。L=0 なら薄レンズ強度 | tffs.f:1316、tquad.f:22-24、tbend.f:566-573 | parser が長さで割って per-length に保存。**TYPE/LINE[] で戻すときは割り算の丸めを避けるため `el._raw` の元の値を返す** (ffscmd.py:357-359) |
| A9 | **BEND の K0 は ANGLE に加算される双極成分** (phib=ANGLE+K0、rho0=L/ANGLE)。ANGLE=0 の BEND は K1=0 なら tsteer、K0=0 なら tquad、他は tmulti に分岐。INTGRL (RADINT) は K0 を見ない | tbend.f:565-611、tbende.f:390-423 | lattice/emit6d に phib/drhob、track6 に tsteer (L=0 は薄いキック) |
| A10 | **CELL の分散は一周行列の δ 列から** (`EX=0.5(a13+a12·a23−a22·a13)/dcosμ`)、x_cod/δ の差分ではない。結合時は s 行列で対角化したモード分散 (EX は PEX と違う) | qcell.f:267-279、tgfun.f:239 | lattice.twiss() を qcell1 の結合版に直訳 (EY が 8e-10 で一致) |
| A11 | **6D 行列の z 行はシンプレクティック条件で組む** (`trans1(5,j) = −t1j t26 + t2j t16 − t3j t46 + t4j t36`)。数値微分のままだと一周行列が 6.6e-5 非シンプレクティックになり、軌道がずれたときの εy が 2% 狂う | tqfrie.f、tbende.f、tsolque.f。Level 3C | emit6d `_State.piece` で同式 |
| A12 | **Emittance[] の閉軌道は z を閉じる δ を選ぶ** (6D tcod)。RF なしでも誤差で周長が変わると δ≠0 (3C: δ=−4.85e-8) | tcod.f | emit6d に close_z (割線法、tcod の直訳ではない) |
| A13 | **qcod は収束後に行列を評価し直さない** (最後に 1 周させた出発点の行列)。収束判定 `conv=min(1e-10, 1e-18(1+(δ/1e-3)²+nlat))` が緩く、off-momentum の βx は初期値の取り方だけで相対 1e-6 動く。初期値は前回の安定な CALC の閉軌道 (先頭 MARK に保存) | qtwiss.f:810-970、tffscalc.f:93,118-120。Level 2 | `_qcod_periodic` 直訳、`Lattice._cod_prev` (twmov 相当)。SAD と同じ順序で twiss() を呼ぶと一致 |
| A14 | **不安定時も 1/cosμ に置き換えて続行** (INTRES/HALFRES が立っているときだけ不安定扱い) | qcell.f:209-248 | `lattice._qmdiag` / `_qcell_periodic_sad` が qcell.f の直訳 (history/work_unstable.md) |
| A15 | **四極のハードエッジフリンジ tqfrie は FRINGE=0 でも既定で入る** (切るのは DISFRIN=1) | tquade.f、Level 2/5 | `_tqfrie_cod`、DISFRIN 対応 |
| A16 | **追跡は CALC/emit と別ルーチン**: tsolqu は eps=0.2 分割 (tsolque は 0.1)、tthin は α/β/γ の 4 次積分 + ttfrin 六極フリンジ (tthine は L/6-2L/3-L/6)、tbendcore は閉じた式、tradkn は `dg=−de/(1+2de)` (tradke と別式)。MULT は追跡 tmulti にフリンジ ttfrin が入り、emit tmulte には入らない | tquad.f、tsolqu.f、tbend.f、photons.f、tmulti.f/tmulte.f。Level 6/8 | track6/track6n は追跡側、lattice/emit6d は CALC 側をそれぞれ直訳 |
| A17 | **BEND のロール誤差 DROTATE は中点まわり** (LGEO=L/2, ANGLGEO=ANGLE/2)。回転の向きは `x' = cosθ·x − sinθ·y` (tsolrote/tchge)。DTHETA は四極では ROTATE、BEND では DROTATE に入る | tfgeo.f:212-213、tsol.f、tchge.f、terror.f。Level 3C/7 | `_tbrote_cod`、`_tchge_cod`/`_bend_element_cod_matrix` (旧 `_bend_err_wrap`)、回転の符号を統一 |
| A18 | **四極ロールの feed-down は座標変換で入る** (tsolrote でずらす→回す→tsolque; tsolque の dx0=ak0x/ak は四極では 0) | tquade.f:30-80 | `_quad_offset_map` 内で (ずらす→回す) |
| A19 | **数値微分の刻み**: pysad 側の注意だが、eps=1e-6 は K1=0 BEND の式の分岐境界に当たり (1,6) 成分に 5.6e-7 の相対誤差 | bench_level16b | `_num_jacobian` eps=1e-7 |
| A20 | **RADCOD ON の tcod は z を一周ごとに trf0 へ移しながら反復し、収束判定が緩い** (`r=Σ(dcod1/codw)²<1e-4`、z は 1e-7 まで)。z と trf0 (Equilibrium dz) の最終値は反復の履歴で決まり、cod の出口と入口も 1e-9 (dp) ずれる (物理量は z+trf0 の和と、その cod での一周行列だけ) | tcod.f、bench_item6_emit.sad | `emit6d._tcod_radcod` が直訳。dz は SAD と 12 桁一致 |
| A21 | **RFSW OFF でも Vc・h・Bucket height・εz・σz・dz が 0 にならない**: tturneg/tecalc が空洞の積算量 (vc0, dvcacc …) から作り、`ωz=√|αp 2π heff Vc/pgev cos φrf|ω0/2π`、trf0 は tturne の update (T(5,6) の符号と U0・Vc) | temit.f tecalc、tturne.f、bench_item6_tiny.sad | `emit6d.tecalc` / `emittance6d` / `closed_orbit6d(norf_trf0)` |
| A22 | **emit の結果は軌道の連続関数ではない**: 量子励起・放射の分割数 (tquade の `nq`, tbende の `ndiv` の `int()`) が cod の x に依存するので、軌道が 1e-14 動くだけで εx 等が 2e-8、σxx が 1e-5 動く。SAD と 1e-11 で一致させるには数値微分の雑音 (刻み 1e-7 で BEND 本体に 1e-8〜1e-7) を無くす必要があった | bench_item6_emit.sad | 解析ヤコビアン (`emit6d._jac_*`) |
| A23 | **EMIT の表示は `autofg(x,'10.7')(1:10)` で最後の桁が切り捨て**: dz=20.61866962 → `20.618669`、f0=497069.3289 → `497069.32`。f0 は `c·β/C` | temit.f temitprint | `fmt_real(x,'10.7')` |
| A24 | **DISP は反転配置 (`-`) の要素名に `-` を付ける** (名前欄は `' -'//name`、`max(10,lname+2)` 桁)。`LINE["NAME"]`・`Twiss` の返値には付かない | tfdisp.f:438-448 | `ffscmd.disp_cmd` |
| A25 | **tchge の入口で th≠0 かつ回した座標系のずれ≠0 (または dz, dtheta≠0) のとき、tbrote/tbgrote/tdrife1 が回転行列 `trans1` を「上書き」する**: COD は回るが一周行列に回転が入らない。出口は tmultr5 で掛けるので回転が入る → 一周行列が R^T·M になり、結合 (R1..R4)・チューンが物理とずれる。適用先: BEND (ROTATE と DX/DY/DROTATE)、薄い多極 (SEXT/OCT/DECA と L=0 の四極)。phig≠0 では tdrife1 (BEND 中点のドリフト) が行列に入る。四極・MULT の tsolrote にはない | tchge.f (tchge の `transa=>trans1` → tbrote、tbgrote の trans 代入)、bench_level25f / 25g | `lattice._tchge_jac` (入口の癖)、`_thin_cod_matrix`、emit6d `_tbende`/`_tthine_gen`。この癖のある要素を含むリングの Emittance[] (δ, チューン) は SAD の tcod の収束で決まるため pysad と一致しない |
| A26 | **BEND の theta の渡し方は経路ごとに違う**: CALC/emit の tsteee・tquade は theta だけ (DROTATE は渡らない)、追跡の tsteer・tquad は theta+dtheta、tmulte/tmulti は theta+dtheta (+akang)、ANGLE≠0 は theta (DROTATE は tbrote で BEND 中点まわり)。tsteee/tthine は K0=0 (ak=0) だと tchge を呼ばない (ずれ・回転とも無効) | tbende.f、tbend.f:565-、tsteee.f、tthine.f | `lattice._bend_chg_args` |
| A27 | **ANGLE=0 かつ K1≠0 かつ K0≠0 の BEND を定義すると SAD K64 (WSL) は FFS 初期化で SIGSEGV (exit 139)** (ROTATE の有無・L に依らない。K1 のみ、K0 のみは動く)。参照値が取れない | bench 25e の初版で確認 (tmp25/bis3.py) | 照合不可。CALC/emit は tmulte (`_bend_as_mult`)、追跡は NotImplementedError |
| A28 | **QUAD の L=0 は tthine(nord=4, ak=K1 の積分値, dx, dy, theta)**: K1<0 の −90° 回転 (akang) は入らず、フリンジなし、放射なし。pf2 の SQ (`L=0 K1 ROTATE=45 deg`) | tquade.f:25-29、tquad.f:25-30、bench 25h | CALC/emit/追跡とも実装 |
| A29 | **DTHETA は DROT のある型 (BEND, MULT) では DROT、なければ ROT に入る** (設計値 + 誤差)。LINE["ROTATE",BEND] は ROT (設計ロール) で DROT とは別 | terror.f:268-285 | `ffscmd._apply_error`, `_AttrMap` |
| A30 | **EMIT 報告 (temitprint) の Extended Twiss は EPX の欄に EY の値を出す** (format 9001 への引数の並びの誤り: `vout9(mfitey)` が 2 回)。`ZX/ZPX/ZY/ZPY/AZ/BZ/PSIZ` はそのまま。`Emittance[]` の ExtendedTwissParameters は正しい (EPX は EPX) | temit.f:1015-1017 (引数), 1021 (format), bench_level26a | `ffscmd._print_emit_exttwiss` が同じ誤りを再現 (報告の行を SAD と一致させるため) |
| A31 | **RFSW OFF の報告でも RF の量 (Vc, h, fs, バケット高さ) は 0 にならない** (空洞の積算量から作る)。`Imag.tune` は log|λ|/2π (0 同然の雑音で `-0.0000000` と `0.0000000` が混じる)。z 方向 (RFSW OFF の Jordan ブロック) の `Real tune` の 0 の符号 (`-0.0000000`) は teigenc の丸めで決まり、pysad は再現できない (値は 0) | bench_level26a | 差は符号だけ |
| A32 | **RFSW OFF でも RADCOD ON なら tcod の rfsw=.false. の枝 (im=5) が放射込みの閉軌道を求める**: x,px,y,py を閉じ、z のずれは trf0 (Equilibrium dz) に移る。出口の δ は −U0/P0 (RF がないので閉じない)。U0 (ESRF で 2.96054e5 vs 2.96075e5 eV)、α_c、チューンが RADCOD OFF と 1e-5 違う。旧 pysad は RFSW OFF のとき RADCOD を無視していた | tcod.f, bench_level26a R3 | `emit6d.closed_orbit6d` の新しい枝 (`_tcod_radcod(rfsw=False)`)、`emittance6d` の最後の一周も radcod=True |
| A33 | **FSHIFT (Orbit dilation)**: RFSW ON のとき `dvfs=FSHIFT` が全要素の `dv` に足される (tsetdvfs)。報告の `Orbit dilation dl = -dleng·1e3 mm`、`dleng = FSHIFT·周長` (tfgeo.f:134)。FSHIFT=1e-6 (周長 603.12 m で dl=−0.60312 mm) で ESRF の閉軌道の δ が −4.3e-3、U0 が 0.2943 MV に変わる | bench_level26c | `_Phys.dvfs` (emittance6d/emit_extra/calc6d の `fshift=`)。追跡 (track6/track6n) は未対応 |
| A34 | **RADINT (INTGRL.f) の癖**: ①tbendebody の bendib の rhob/rho0 を INTGRL は設定しないので、直前に tbend 系が設定した値 = **最後に処理した曲がりのある BEND の L/ANGLE が全 BEND に使われる** (ESRF は反転した ARCAH2 の BID5=20.757 m。ほかの値だと I9 が 1e-3 ずれるので一致で裏づけ)。②六極のブロックは DL を設定せず直前の BEND/QUAD の最後のスライス長を使う (I9, I10=γ·DL が六極でも足される)。DKL2 は代入されず六極の CKL0X/CKL1X は 0。L=0 の六極は TCAL の後で DL=0.2 (単精度)。③`1e27`、`1.E-6`、`0.2`、`2.3026` は単精度 (ELOSS が相対 1.15e-8 小さい)。④CAVI は VOLT/HARM/FREQ を最後のもので上書き。⑤CKL1X=DKL1·DX/DR は軌道の向きだけ。⑥MASS の既定は elmass=0.51099895e6 (pysad の既定 0.51099906e6 は誤りだった) | bench_level29a, 29a2 | `pysad/radint.py` (stale_rho は最後の BEND の L/ANGLE)。history/work_small_items.md A |
| A35 | **FFS 内の `ON`/`OFF` の語は無効**: `Out[1]:= USE` を出すだけで、後ろの語は普通に次の文として読まれる。`ON INTRA;` ≡ `INTRA;` (立てる)、`OFF INTRA;` も INTRA を**立てる**、`ON NORAD NOFLUC;` は RAD, FLUC を切る (`ON;` も `Out[1]:= USE`)。MAIN 部の ON/OFF は有効 | bench_level29d | `ffscmd._exec_word`。history/work_small_items.md D |
| A36 | **`CALC6D;` / `CALC4D;` は単にフラグを変えるだけで再計算しない** (Twiss[] は直前の CALC のまま)。**4D の CALC は縦方向の列 (AZ BZ NZ ZX ZPX ZY ZPY) を作り直さず**、直前の 6D の CALC の**先頭行**の値を全行に残す (tffscalc.f:87-90 は 4D では 15 列だけ入れ直し、qtwiss1 が縦の列を直前の行から写す) | bench_level29c | `FFS._z_stale`。history/work_small_items.md C |
| A37 | **DISP の位置・範囲** (tfdisp.f/ielm.f): 位置を省くと前回の範囲 (id1, id2。USE で 1..nlat に戻る)。`loc1==loc2` はその要素から `$$$` まで (2 行)、`loc2<loc1` は周回。数式の位置は `int(v+0.499)`。`NAME` は代表、`NAME.n` は複数あるときの n 番目 (1 台だけの `RFC.1` は未定義エラー)。`*` を含む語は数式でなくパターン (`disp 2*47 100` は何も出ない) | bench_level29b | `ffscmd.disp_cmd`。history/work_small_items.md B |
| A38 | **ANGLE=0 で K1≠0 (K0=0) の BEND の追跡は K1 の符号を失う** (tbend.f:573 の tquad 呼び出しが theta2=theta+dtheta を渡し、akang (K1<0 の π/2) を含めない。tquad は ak=sign(ak0,al) で符号を theta2 に任せる) → K1<0 の BEND は |K1| の集束四極として追跡される。CALC (tquade) は正しく K1<0。SAD で QUAD K1=+0.5 と BEND ANGLE=0 K1=−0.5 の TrackParticles が全桁同じ、QUAD K1=−0.5 とは違う | <work-perf>/bendk1_*.sad | track6 `_bend_k1_proxy` が |K1| で tquad を呼ぶ (再現)。history/work_trpt_acc.md 10 節 |
| A39 | **EMIT コマンドは GAUSS フラグを ON にする** (`UNIFORM; emit; ?GAUSS` → 1)。`Emittance[]` は変えない。以後の BeamMatrix[] の dp² は GAUSS の値 (UNIFORM なら dp²/3) | <work-perf>/hiprio9_gauss.sad, hiprio11.sad | `ffscmd` の EMIT が `flags['GAUSS']=True`。history/work_trpt_acc.md 11 節 |
| A40 | **INTRA ON の EMIT 報告は 2 段** (temit.f): 1 周目の temitprint が零電流の報告、収束後に Touschek の表 (運動量アクセプタンスごとの寿命と、Nx/Nz・Ny/Nz ごとの寿命/100s の表: 9104/9105/9131) を出し、最後に `Parameters with intrabeam scattering:` + `Particles/bunch/1e10 = … Minimum coupling = … %` + IBS 込みの EMITXC/EMITYC/EMITZC・SIGE・SIGZ の 9102 ブロック。Emittance[] は表を出さない。EMITX などのグローバルは IBS 込みの値 | <work-userscripts>/pf2 touscheck_r, disp_bunchlength_r | `ffscmd.run_emit` → `print_emit(e0)` → `_print_touschek` → `_print_emit_intra` |
| A41 | **PBUNCH の既定値は 1e10** (initb1.f RsetGL)。INTRA の EMIT で PBUNCH を設定しないスクリプトも Touschek の表が出る (NBUNCH=1, MINCOUP=1e-2, GCUT=1e300, LOSSAMPL=1, LOSSDZ=100, DAPWIDTH=7 も同じ表) | <work-userscripts>/pf2 disp_bunchlength_r | `ffscmd` のグローバル既定値 |

## B. 乱数と並列 (再現には順序が重要)

| # | 現象 | 根拠 | pysad での扱い |
|---|---|---|---|
| B1 | **乱数は 32 bit 符号付き合同法** `seed ← 48828125·seed` (int32 で折り返し)。`Random[]`=(seed+2³¹)/2³²、`GaussRandom[]` は Box–Muller (φ=seed·π/2³¹, ρ=√(−2 log(seed/2³²+0.5)))、ρcosφ を返し ρsinφ を次回に回す、\|値\|>GCUT は棄却 | src/random_plugin_sad.c:17,107-128 | `pysad/sadrandom.py` (15 桁一致) |
| B2 | **SeedRandom[n] の種は `|n/2|·2+1` (奇数化、負なら負)**。FFS 開始時の種は 17。seed 17 の GaussRandom 列 0.617158…, 1.657394…, −0.476296… は誤差コマンド DELX/DELY/DTHETA の値そのもの | random_plugin_sad.c:54-68。Level 7/10 | `SadRandom.restore` |
| B3 | **SeedRandom[] の戻りは `{"SAD", seed (符号なし = seed mod 2³²), filled, next, 0}` の 5 要素** | Level 16 | `dump()`、ffscmd の SeedRandom[] (符号なし表示) |
| B4 | **FLUC の乱数を引く順**: 要素外側・粒子内側。放射ありの BEND (tbendr)・四極 (tsolqur)・MULT はスライスごとに全粒子の本体→全粒子の tradk。tbendi とソレノイド内ドリフトは粒子外側・スライス内側。粒子ごとに tdusrn で p=tran() 1 回、光子ごとに drndsr 1 回 + r 1 回 (光子 0 個でも p は引く) | photons.f tradkfn、drndsr.f。Level 10/12 | `*_ens` ルーチンで同順。乱数の最終状態まで一致 |
| B5 | **NPARA の既定は Packages/init.n で 4** (initb1.f の 1 ではない)。TrackParticles は粒子 9 個以上 (npparamin=9) で fork し乱数列が分かれる。NPARA>1 だと子プロセスが出力を重複して書く | Packages/init.n:219、tftrack.f。Level 16 | 参照値は `NPARA=1;` を入れて取る。pysad は `NPARA` 1 固定 |
| B6 | **損失粒子の NaN/inf と乱数**: 粒子が落ちる周は発散中の値が NaN/inf/巨大値のまま tradkfn が光子の乱数を引く。θ=π/2 (有限の巨大値) なら anp≈94 で約 100 個、ppa が NaN (inf の 2 成分の norm2) なら θ=−π/2 → anp<0 で 1 個、anp=NaN (g=NaN など) なら tdusrn が 10000 回ループして 20001 個。gfortran -Ofast の max/min の NaN の返り値はサイトごとに違う (tradkfn の `max(gmin, g+dg)` は NaN で gmin、asinz の `max(-1,x)` は −1、sqrt1 の `max(1e-200,·)` は後者)。tsolrot は chi3=0 なら回さない (0·inf=NaN を作らない)。発散中の分岐は丸めの 1 ulp で逆になるので、自由走行の最終 seed は SAD と一致しない (SAD 自身も 1e-16 の摂動で変わる) | drndsr.f, photons.f (tradkfn), tsol.f (tsolrot), tbend.f (tbendcore の asinz)。bench_level16f/16l/20 small | track6n: `sqrt1`/`pxy2dpz`/`asinz`/`norm2_3`/`gfloor`/`rot1`、`Prog.la_count`。SAD デバッグ版との再同期テストで 30 周 seed 一致 (history/work_loss_rng.md, history/work_loss_rng2.md) |
| B7 | **DynamicApertureSurvey の粒子スロット順** は trackd と同じ (i=1..n1p, j=51→1)、損失時は末尾の生きている粒子と入れ替える (tapert)。損失判定は最初の要素、以後 15 要素ごと、周末。n1 の上限は 200 (help の 256 ではない)、ntloss は 51 点 (help の 50 ではない)、DAPWIDTH はグローバル変数 | trackd.f、tapert.f、tfdapert.f:56。Level 16 | dapert.py + track6n.k_apert。NORAD の DA は地図 153 文字一致 |

## C. FFS / SADScript の仕様 (スクリプトの書き方に効くもの)

| # | 現象 | 根拠 | pysad での扱い |
|---|---|---|---|
| C1 | **LINE[key,"QD.1"]=v は CALC で全 QD に展開** (QD.1 は代表 klp)。`LINE["DX","QD.5"]=v` は CALC で代表の値に上書きされて効かない。1 台だけの誤差は FFS コマンド `DELX σ QD.1` (QD* で全台に別々の乱数) | tfsetv.f tffsadjust:95-148、terror.f。Level 1/3B | `set_element(fam,key,v)` (ファミリー) / `set_error('QD.1',key,v)` (1 台)。ffscmd の DELX… は SAD の乱数で同じ値 |
| C2 | **CELL では入口 COD が先頭 MARK に書き戻される** (DUMP の ORG.1)。先頭 MARK の AX/BX/…/DDP が前回の解の保存場所 | tffscalc.f:118-120, 540- | `_cod_prev` |
| C3 | **ドリフトなしで続く要素は上流の要素と一緒に動く** (SF だけに誤差を入れられない) | SAD 実行 (Level 3B で DQS を挟んだ) | 未実装 (ベンチ側で回避) |
| C4 | **FIT; は条件を消さない。** 引数なしの FIT は位置を `$$$` (終端) に戻すだけ。消すのは `REJECT TOTAL` | tffsa.f:509-531、tfrej.f:15-19 | ffscmd.py:1072 同じ |
| C5 | **FFS["cmd"] は出力ユニット省略で表示なし** (outfl=0)。`FFS["go",6]` の 6 は出力ユニット (−1 で現在の出力先)。戻り値は計算結果のリスト | tffsa.f:1807-1812 | ffscmd.py:700- `fn_ffs` |
| C6 | **PolynomialFit は平均 0 (xa=0) のとき xl=1 になる** (`If[xa, Table[...], 1]`) ので係数の次数補正が壊れる → ユーザーは MyPolyFit を使用 | Packages/Fit.n:37-56 | sadscript.py:877- が Fit.n をそのまま評価 (癖も再現) |
| C7 | **tffsmatch の fuzz は NaN でモードを変えない。** 残差が悪化すると crate<0 → log10 が NaN → gfortran の min/max は NaN を無視して badcnv=0 → DESCEND に切り替えず fmincube の直線探索。再現しないと早く "Residual =" で止まる | tffsmatch.f:255-289。esrf/matching | match6.py:113-121 (`_log10` が NaN を返し fuzz が 1) |
| C8 | **MAXITERATION は反復上限そのものではない** (alate のランプと DESCEND の 10·itmax)。収束は `res2r(r) ≤ CONVERGENCE·max(nfcol,1)`。OffMomentumWeight の分母は help の √(点数) ではなく √\|mfitp\| | tffsmatch.f:98,162-177、tffscalc.f:336-345 | match6 直訳 |
| C9 | **マッチングの微分は既定で解析的** (tffsqu/qdcell)。数値微分 (FFS$NumericalDerivative=True) だと解が複数ある問題 (esrf/matching) では別の解に着く | tffsmatch.f:320-393 | `match6d.py` (解析微分) が Matcher の既定。esrf/matching の反復列と K1 を 12 桁再現 |
| C10 | **LINE RING=(10*(A B)) の入れ子は構文エラー** → `LINE CELLP=(A B); LINE RING=(10*CELLP);` | lread.f。Level 13 | parser も同じ制約で書く |
| C11 | **FFS で使う LINE の先頭は MARK 必須** ("The first element must be a MARK element.") | tffs.f:2988-2991 | 同じ前提 |
| C12 | **Mainlevel では SADScript (Print[] 等) は使えない**。FFS 進入時に tfinitn で初期化。FFS に入ると RFSW=ON、CALC6D=OFF、TRPT=OFF に強制 | toplvl.f、track.f:104,146-153 | run_sad.py: MAIN 部は parser、`FFS;` 以降は sadscript/ffscmd |
| C13 | **MON/MCOD/TCOD はこのビルドでは何もしない** (語が式として評価されるだけ)。**後ろの語は次の文として読まれる**ので `MON 0 0 0 0 PM*` は `0` が式、`PM*` が値の無い要素値コマンドで `?Missing value for PM*` (FFS[] の中なら FFS::error で Do を抜ける) | SAD 実行で確認、<work-perf>/mon_test.sad | `ffscmd` は語を 1 つ飛ばすだけ。数値の語は式として無視 |
| C14 | **Emittance[] の EnergyLossU0 / RfVoltageVc は eV**。TransferMatrices は先頭から各要素入口までの**累積**行列。NormalCoordinates は ri (物理→正規) | Packages/emit.n:47,63,76。bench_level16b/16e | ffscmd.py:729 |
| C15 | **`type` の出力書式**: 値はタブ位置 26, 31, 36, … (次の位置まで空白、最低 1 つ)、数値は最短表現で先頭の 0 を書かない (.413) | tftype.f。Level 16 で書式一致 | ffscmd.py:1563- `type_cmd`, :83 `_shortest` |
| C16 | **Twiss[...,"QD.1"] は要素の入口**。`Twiss[...,6]` は 6 番目の要素の入口 (table[5]) | Level 1/5 | 同じ |
| C17 | **`FFS["QD.1 DX 1E-4"]` はフィット条件** (閉軌道 DX を 1e-4 にする) として解釈され、誤差設定ではない | Level 2 | 誤差は DELX… |
| C18 | **USE の引数に変数は使えない、Do[] の反復子に文字列リストは使えない (多分)、InputForm は関数ではない** (ToString の書式)。桁数は `$FORM="15.10"`。Do[] の中の FFS コマンドは FFS["…"] で包む | Level 3A/5 (ユーザー指摘) | sadscript は同じ制約 |
| C19 | **DUMP は 5 桁しか出ない** → 10 桁は `LINE["DX","QD.1"]` で取る | Level 7 | — |
| C20 | **色収差は CHRO / RADINT の値を使わない** (BEND 成分なし)。DP0 スキャン + 多項式当てはめで求める。RADINT と Emittance[] は手法が違い値も違う (L2: 529 nm vs 460 nm) | INTGRL.f、temit.f。Level 2/3A | emit6d は Emittance[]、`radiation_integrals` は RADINT |
| C21 | **COD フラグはエミッタンス計算時のフラグ** (temit の calcod)、CALC の閉軌道には関係しない。RADCOD は追跡では taper 条件にしか現れず、emit/COD 経路で放射損失を閉軌道に反映する | tturn.f:317-321、photons.f:1420-1430 | emit6d の radcod 分岐 |
| C22 | **QUAD 定義に RAD=1 を書くと COD・エミッタンスが 0 になる** (原因不明) → 切り分けには DampingRate を使う | Level 3C (sadrun/c3ye_nq) | — |
| C23 | **ExpandElementValues->True (既定)** で Emittance[]/DynamicApertureSurvey は Element["EXPAND"] を呼ぶ | emit.n:25,34、dapert.n:17,29-30 | ffscmd 内で同様 |
| C24 | **SADScript の `1GEV` は `1*GEV`** (数値の直後の英字は掛け算のシンボル)。GEV/MM/DEG は未定義なので `MOMENTUM = 1GEV` は MOMENTUM=1 になる (ラティス定義の `L=3MM` や FFS コマンド行は単位として読む)。MAIN/FFS の SADScript では `1E9` と書く | <work-perf>/hiprio.sad, unit_test.sad | pysad は字句解析で単位を掛けて 1e9 にする (意図的な拡張)。SAD 互換の結果が要るなら `1E9` |
| C25 | **FFS["…"] の中のコマンドがエラー (`?No variable.` 等) を出すと `???FFS::error:  Error in FFS command n in FFS["…"]` を出して、それを含む最上位の文 (Do など) を打ち切り、次の文へ進む**。エラー行は出力ユニット 0 (表示なし) でも表示される。`LINE[key,"NAME"]=v` で要素が無いときも `???General::invset:  Invalid set to LINE in LINE[…]` で同様に打ち切る。pysad は以前は続行していた (dx_sv/dy_sv/chrom_b/disp_b/matching21/amp_dep_ts_x_octa の出力差の原因) | <work-perf>/err_test.sad | `ffscmd.fn_ffs` (`_RecOut` でエラー行を検出して SADError)、`set_line`、`sadscript.run` が `... in 文` を出す |
| C26 | **qcod が溢れて閉軌道が見つからないとき SAD は止まらず壊れた値で続行する** (`*****qcod---> Overflow & closed orbit not found`、Twiss は直前の値や 1e200、チューンは 1/cosμ 置換の値)。pysad は直前の twiss があればそれを、無ければ NaN の表で続行する (値は SAD と違う) | <work-userscripts>/esrf amp_dep_ts_dp (dp=−0.04) | `ffscmd.calc`、`match6._stale_rows`、`lattice.nan_rows` |
| C27 | **`LINE[key,"FAM"]=v` (番号なし) の直後の読み出しは代表 (FAM.1) だけが変わる** (他の台は次の CALC で揃う)。pysad は直ちに全台に入れる (読み出しが CALC 前だけ違う) | <work-perf>/use_test.sad (`LINE["K2","SF"]=5` → SF.1=5, SF.3=2) | `set_family` (差は記録のみ) |
| C28 | **USE / VISIT はラティス定義から作り直す**ので、FFS で変えた要素値 (LINE[…]=、`QF K1 v`) は USE/VISIT/BYE で定義の値に戻る (SAVE しない限り) | <work-perf>/use_test.sad | `FFS.use` が `env.build_lattice` で作り直す (同じ) |

## D. ビルド依存・環境の挙動

| # | 現象 | 根拠 | pysad での扱い |
|---|---|---|---|
| D1 | **WSL ビルドは FITP>1 の GO で再起動する** ("*** Run time Environment" からスクリプト再実行、NPARA=1 でも)。参照値が取れない | bench_level15a、esrf/chrom | match6/ffscmd の FITP>1 (idp ループ、重み、INTRES 判定、fam 点) は直訳のまま SAD と未照合 (自己検証のみ) |
| D2 | **gfortran の min/max は NaN を無視する** (C7, B6 の根拠)。Python の max/min は前の引数を返す | tffsmatch.f fuzz、tapert/photons | match6: `_log10`→NaN で fuzz=1。track6n: NAN_MODE |
| D3 | **esrf スクリプトは datas.txt に OpenAppend で書く** → ユーザーのディレクトリでは実行しない (<work-esrf> にコピー) | esrf/*.sad | run_sad.py は <work-esrf> で |
| D4 | **Bash ツールの heredoc はバックスラッシュが 1 段外れる** → Python パッチは Edit ツールで | 作業メモ | — |
| D5 | 非コンパイルのファイル (pqcell.f, qquad.f, qthin.f, tdfun.f, tquads.f 等 122 本) を参照しない。qquad/qthin の実体は qtwiss.f 内、tdfun は qdcell.f 内 | mk/sad.obj.mk OBJOBS、obj/AMD64-Linux5 | 直訳元はコンパイル対象のみ |

## E. 未解決・要確認

- B6: 損失周の乱数消費は論理が再同期テストで一致するが、最終 seed はカオスで一致しない (history/work_loss_rng2.md)。
- C3: ドリフトなしで続く要素の連動 (未実装、ベンチ側で回避)。
- D1: FITP>1 の照合 (SAD 側が落ちるため別ビルドか別環境が必要)。
