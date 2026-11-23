# 44_sol_guide の履歴 (整理で本文から取り除いた記述)

本文 (`docs/44_sol_guide.md`) を版のないスナップショットにするとき、取り除いた版・経過の記述を、元の文脈つきで残したもの。元の行をそのまま引用する (引用した行のうち、履歴にあたるのは版・経過を述べた部分)。

## 元 3 行 (節: 44. ソレノイド (SOL) の解説)

> 最終更新: 2026-11-18 (v2.1.0 の後)

## 元 16 行 (節: 1. 何ができるか)

> | POL のスピン (Emittance[]・EMIT・追跡) | 対応 (k64 モード: 領域内は DRIFT・QUAD。v2 モード: 追跡のスピンは領域内の BEND・MULT・CAVI・ずれ・F1 も。Emittance[] の POL の量は k64 では領域内が K64 の式、v2 は v2.1.0 で修正: B130) | [`docs/36_spin_csr.md`](36_spin_csr.md)、docs/43 §8a |

## 元 17 行 (節: 1. 何ができるか)

> | CALC6D、RADCOD つきの CALC | **v2.1.0 で対応 (両モード)**。k64 は gs と 1e-9 | `ffscmd._check_sol_optics`、`calc6d._turn` |

## 元 18 行 (節: 1. 何ができるか)

> | Touschek (INTRA) との組み合わせ | **v2.1.0 で対応 (両モード)**。INTRA + SOL の Emittances は gs と εx・εy 2e-8、εz 9e-8 (相対)。TouschekLifetime は gs・pysad とも NaN (この条件で確認できず) | `touschek._tturne_intra` |

## 元 100 行 (節: Twiss (CALC / INS) (`solmap.py`、K64 の `tsole1`))

> `CALC` のたびに `solmap.annotate` が領域内の要素に `_solmap` を付ける。領域内の Twiss の座標は領域内の正準運動量での値 (K64 の `tsole1` と同じ)。幾何 (GX など) は従来どおり。

## 元 109 行 (節: CALC6D)

> **v2.1.0 で対応した (両モード。pysad の不足を gs に合わせた: gate なし)**。SOL (BZ≠0 または BOUND≠0) を含む CALC6D と RADCOD つきの CALC は、`_check_sol_optics` が領域を annotate して通し、`calc6d._turn`・`touschek._tturne_intra` が領域内の要素 (`_solmap`) を `solemit.emit` で通す (`emit_extra` と同じ)。K64 は `qtwiss` の 4D が領域を `qsol` (`tsole1(.. enarad=.false.)`)、CALC6D が `qtwissfrac` の `tsole1(…, .true., .true.)`、INTRA の `tturne` が領域内の要素を `tsole1` で通す。確認 (gs、SOL 領域のあるリング、3 GeV、BZ = 2 T、回転四極、QFE の DX・DY): CALC 4D 2.73475252978800・1.87531224259842 (pysad k64 2.73475252978794)、CALC6D 2.73495437251036・1.87554377922481・NZ −.01173531080122 (2.73495437251039 など。1e-9、NZ は絶対 1e-9)、RADCOD つき CALC6D 2.73497629185364 (1e-10)、RADCOD の CALC6D のあとの 4D の CALC 2.73465818826401 (pysad 2.73465818826743。CELL CALC6D のあと `DP0` を先頭の `DDP` にする)、β (1・30・70 番) も 1e-9。INTRA + SOL の Emittances は {2.706114556595e-7、2.245979910482e-8、2.778701736555e-5}。未確認: TouschekLifetime の値 (gs が NaN)、BOUND=0 の SOL (BZ≠0) の CALC (NotImplementedError)、RADTAPER + SOL、WSPAC・INTRA + SOL の併用の値。例 [`example/v2_next/B1_sol_calc6d_intra.sad`](../example/v2_next/B1_sol_calc6d_intra.sad)。

## 元 114 行 (節: 幾何 (`sgeo.py`、K64 の `tsgeo.f`))

> 軌道が軸上 (DX DY DPX DPY が 0) のときは直進と同じで従来の経路。`LINE["S"|"GX"|"GZ"|"GCHI2"|"GCHI3", k]` が gs と一致する (§6)。

## 元 164 行 (節: 7. 限界・未確認)

> 9. k64 モードでは領域内の MULT・BEND・CAVI の `tmulte1` と、`tsgeo` の領域幾何、CALC6D、Touschek は未翻訳。`tdgeo` の領域内の変数は v2 モードで中心差分に対応 (§7a)。CALC6D・Touschek (INTRA) は v2.1.0 で両モードとも対応 (上の CALC6D の節)。

## 元 205 行 (節: ずれて傾いたソレノイド (`DX DY DZ CHI1..3`・`F1`))

> - **使えない**: k64 モード (CALC6D・RADCOD つきの CALC は v2.1.0 で k64 も対応)。追跡を領域の途中から始める・途中で終える (`TrackParticles` の区間の端が領域の中: `ソレノイドの終端がない` のエラー。従来どおり)。領域より上流の幾何に軸まわりの回転 (GCHI3≠0) があるラティスで、追跡を先頭以外から始めた場合 (区間ごとに幾何を作るので CHI3 が入らない)。領域内の APERT のスピン追跡 (NotImplementedError)。領域内の四極の DX・DY のある POL の追跡 (track6) は NotImplementedError。領域内の光子は記録しない (乱数だけ K64 と同じに引く)。

## 元 211 行 (節: スピン追跡と Emittance[] の POL)

> - **Emittance[]・EMIT の POL の量 (v2.1.0 で v2 は修正。B130)**: k64 は領域の中が K64 の式のまま: 軸上のソレノイドでは追跡と 2e-5 で合うが、ずれて傾いたソレノイドではスピンチューンで 1.5e-4〜1.1e-3、n0 の向きで 2.5e-3〜1.5e-2 rad 違う (縦磁場の積分を z 軸の向きに足す)。**v2 は B⊥ と B∥ を運動量の向きに分けて足す**: SpinTune の追跡との差 1e-13〜5e-7、n0 の差 0〜1e-5 rad (傾き 3・10 mrad、ずれ 2・4 mm、BZ = 2 T・0.8 m の領域。小さな FODO・ASTRID2・MPW 入り・PF・大きなリング A)。回転した縦磁場の歳差の閉じた式と追跡が 5e-7 (傾き 1〜50 mrad・1〜3 T・0.3〜1 GeV)。領域内の DRIFT・QUAD の放射・スピンの刻みは bz·L/n ≤ 0.02 (DRIFT の刻みは B1 で追跡にも入れた: bz = 4.5 で 1.5e-3、9 で 8.7e-3 ずれていたスピンが 2.3e-7・8.9e-7、らせんの放射損失が 1e-4)。使えない・未確認: bz ≥ 3 1/m の QUAD・MULT の D は 1e-4 台 (放射の分割が EPS 由来)、領域内の QUAD の端の ψ は追跡と同じ模型 (hard edge の Bs なし。軌道が数 mm の範囲で確認、cm のずれは未確認)、領域内の CAVI の POL (Emittance[]) は VOLT=0 だけ確認、領域の外の B∥ (垂直な軌道の傾きの中の偏向磁石の B_y y')。docs/43 の「傾いてずれたソレノイドの POL の量」。

## 元 213 行 (節: v2.1.0 の追加 (領域内の四極のずれ・bz ≥ 3・領域内の CAVI のスピン))

> ### v2.1.0 の追加 (領域内の四極のずれ・bz ≥ 3・領域内の CAVI のスピン)

## 元 217 行 (節: v2.1.0 の追加 (領域内の四極のずれ・bz ≥ 3・領域内の CAVI のスピン))

> - **領域内の四極の `DX`・`DY`** (B3。両モード。gs と同じ): K64 の `tquad` は `tsolrot` で座標を (dx, dy) ずらす。ソレノイドの一様磁場のベクトルポテンシャルのゲージが変わるので正準運動量もずれる: x → x − dx、y → y − dy で px → px + (bz/2)dy、py → py − (bz/2)dx (出口は逆)。追跡のカーネル (`k_sol_quad`)・`solemit`・`sgeo` に入れた。gs と `LINE["GX"…]` が全桁一致 (領域内の DX・DY・ROTATE つきの 2 台)。**幾何の fit** (QF の DX、QD の DY で M2 の GY = 0.0008、GZ = −0.0035): gs と v2 (解析微分 `match6d.sol_dcod`: 要素の写像 c_out = F(c − t) + t、t = (DX, −(bz/2)DY, DY, (bz/2)DX) より d c_out/d DX = (I − J) ∂t/∂DX) は 2 回の Newton で収束し、QF DX = 2.78848432670e-3 (gs) / 2.78848432669e-3 (v2)、QD DY = 2.66798056318e-3 (両方)。解析微分は幾何の有限差分と 2e-3。K64 の CHI3 の微分の式は領域の終端で GCHI3 が 0 に戻される (tsgeo の `tfderotgeo`) ので幾何の変化と合わない (9e-4 を返すが幾何は動かない): CHI3 を使う fit は K64 でも収束しない (B124)。使える範囲: 領域内の四極 (K1≠0・F1=FRINGE=0)。SEXT・OCT・MULT は `_QDX_TYPES` に入れたが領域内のずれは未確認。四極の F1・FRINGE と領域内の BEND・MULT のずれは未対応。k64 モードは SOL 自体のずれ (trots) が未対応でずれた SOL の fit は使えない。例 [`example/v2_next/B3_geofit_quad_dx.sad`](../example/v2_next/B3_geofit_quad_dx.sad)。

## 元 218 行 (節: v2.1.0 の追加 (領域内の四極のずれ・bz ≥ 3・領域内の CAVI のスピン))

> - **bz ≥ 3 1/m** (B1): 領域内の DRIFT・BEND・MULT (K0・K1・K2) の写像 (力学的運動量) がローレンツ力の積分と bz = 6・10・20 1/m で 1.2e-8 以内 (多極の分割 EPS の分)。放射損失 −(2/3) r_e γ³ (K0/L)² L が bz = 6・20 で 4e-5。**DRIFT の刻みを足した** (v2): 領域内の BEND・MULT は bz·al/ndiv ≤ 0.02 に刻んでいたが DRIFT は 1 刻みだったため、軸を外れた粒子 (δ = 1e-3・−5e-3、0.4 m) のスピンが Thomas-BMT の積分と bz = 4.5 (0.1 GeV) で 1.5e-3、9 で 8.7e-3 (0.3 GeV の bz = 3 で 1.3e-4) ずれ、放射も bz = 9 で 1.7 倍外れていた (`track6n_ext.k_sol_drift`・`track6._sol_drift_v2`。ndiv = 1 + int(\|bz·al\|/0.02)。結果 2.3e-7・8.9e-7・4.5e-8、放射 1e-4)。k64 は 1 刻み (K64 のまま)。QUAD (K1 = 0.3、0.4 m) のスピンは bz = 4.5・9 で 2.0e-5・9.0e-5、MULT (K0 + K1) は 1.3e-4〜1.6e-4、BEND は 6e-7・9e-7 (放射の分割が EPS 由来で bz·al/ndiv ≤ 0.2 のため。v2 の既存の試験値を変えないため足していない)。

## 元 222 行 (節: v2.1.0 の追加 (領域内の四極のずれ・bz ≥ 3・領域内の CAVI のスピン))

> - **領域内の四極 (K1<0 の −90°、ROTATE) のスピンの座標系** (E。pysad の不具合。両モード): 回していなかったのを直した (track6 の `_tquad_sol`・`solemit._tquade`。K64 の `tsolrote` は `srot = rr·srot`)。k64 の Emittance[] の srot が gs と全桁一致 (BZ = 1 T、K1<0・ROTATE つきの四極)。

## 元 224 行 (節: v2.2.0 の追加 (記録 P4。試験 `test_v2_P4_sol_spin.py`・`test_v2_P4_sol_photons.py`))

> ### v2.2.0 の追加 (記録 P4。試験 `test_v2_P4_sol_spin.py`・`test_v2_P4_sol_photons.py`)

## 元 226 行 (節: v2.2.0 の追加 (記録 P4。試験 `test_v2_P4_sol_spin.py`・`test_v2_P4_sol_photons.py`))

> - **領域内の QUAD・MULT のスピン (B146。v2)**: POL のとき領域内の QUAD を bz·al/ndiv ≤ 0.02 に刻む (Emittance[] の `solemit._tquade` と同じ。以前は刻み bz·al/ndiv ≈ 0.2 の弦の誤差で bz = 4.5・9 で 2e-5・9e-5)。端の縦磁場の積分 ψ (QUAD の k x y、MULT の k0 y + k1 x y + …) は端のフリンジを入れるとき (DISFRIN=0) だけ入れる (以前は DISFRIN に関係なく足し、MULT の K0 つきで 1.6e-4 ずれた)。Thomas-BMT の積分 (0.4 m、DISFRIN=1、軸上と軸を外れた粒子・3 つのスピンの向き) と bz = 0.3・3・4.5・9 1/m で **9.3e-7 以下**。D は追跡 1e-8・Emittance[] 2e-6、ASTRID2 (RF の前に BZ = 0.5 T・0.9 m の領域に QUAD・MULT) の Emittance[] の SpinTune が追跡の 1 周のスピン行列と 1.2e-7。領域内の CAVI (VOLT = 0・0.2 MV) の SpinTune は追跡と 1.6e-7 だが、位相 0 なので SpinTune がエネルギー変化にほぼ鈍い弱い試験。k64 は変えない。**未確認**: DISFRIN=0 の ψ を硬い端の B_s と軌道のフリンジ写像を合わせた厳密な基準で確かめていない、ずれ・ROTATE を持つ MULT の ψ、CAVI の加速・位相が 0 でないとき、大きなリング A・PF・8 GeV (0.1〜1 GeV・ASTRID2 までのリングで確認)。

## 元 228 行 (節: v2.2.0 の追加 (記録 P4。試験 `test_v2_P4_sol_spin.py`・`test_v2_P4_sol_photons.py`))

> - **領域内の加速 (POL + TRPT)** は実装しない (v2.2.0 の判断): 刻みごとの基準運動量 p0(s)・bz p0/p(s)・歳差の aγ(s) の作り直しが要り、RF の電場を含む運動方程式の積分の独立な基準が無い。実装するなら `sol_mult_full` の刻みで `_RadHook` が基準運動量を更新し `_tradk_sol` に p(s) を渡す。現状は明示のエラー (`test_pol_with_trpt_acceleration_is_an_explicit_error`)。

