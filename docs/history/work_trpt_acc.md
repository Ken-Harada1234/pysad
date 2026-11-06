# 作業ログ: USE の大文字小文字、TRPT の加速 (gammab / tcave / tmulte)、MARK OFFSET、K=0 の QUAD/MULT (2026-11-06)

残課題「USE LinacPFBT (disp_bt_r)」。調べると USE 自体は 1 行の不具合で、その先に線形加速器 (13.2 MeV から) の加速と MARK の OFFSET が
未実装という本体があった。disp_bt_r (linac-pfbt.sad、3546 行の表) は最終的に **全列が印字精度で SAD と一致**。

## 1. USE の大文字小文字 (直接の原因)

linac-pfbt.sad の `LINE LinacPFBT=` は大文字小文字混じり。parser は書かれたとおりに登録するが、ffscmd の USE / use() /
sadelem の ExtractBeamLine は名前を大文字にして引いていた → 「Missing beamline in USE.」。
`Env.resolve_beamline(name)` (書かれたとおり → 大文字 → 大文字小文字を無視) を足し、USE / use() / build_lattice / ExtractBeamLine が
それで引く。`ffs.name` は登録名のまま (BeamLineName[] も)。MAIN の `TRACK USE=` / `EMIT USE=` の比較は大文字小文字を無視。

## 2. TRPT の加速 (pysad/trptacc.py、新規)

SAD の仕組み (読んだ所):
- tfgeo.f tfgeo1 (acconly): `trpt` のとき CAVI/MULT ごとに **tfacc** で基準運動量 p (単位 amass) を進め `gammab(i+1)` に入れる。
  tfacc: v=|charge|·VOLT/amass、dh=max(oneev−h1, −v·sin(PHI·sign(charge)))、h2=h1+dh、p2=p1+dh(h2+h1)/(p2'+p1)。
- qtwiss.f: CAVI は qcav → **tcave** (`rfsw = rfsw .and. trpt`: RING では RFSW に依らずドリフト = 11_sad_quirks A4)、
  MULT (VOLT≠0) は qmult1 → **tmulte** (`acc = vc≠0 .and. rfsw`)。
- tcave.f / tmulte.f の加速: ndiv = 1+int(min(|w L|, sqrt((v(1/h0+1/h1))²/3/eps))) (tcave eps=1e-2、tmulte eps=eps00m·EPS=0.005·EPS、≤300)、
  入口 tcavfrie (L≠0, DISFRIN=0, FRINGE≠2) → [ドリフト (最初 L/2ndiv、以後 L/ndiv) + 基準粒子の速度補正 tgetdvh(dgb)·L → キック]×ndiv →
  ドリフト L/2ndiv + tgetdvh → 出口 tcavfrie(−v)。キックは dh=max(oneev−h1, −va(sinφ+offset1))、φ = w t + PHI·sign(charge) − trf0·w、
  t = −z/v1、va = vn + v20a x² + v02a y² (v20a = v02a = vn (w(½/p0+½/p1))²/4: RF の横集束)。trpt では offset1=s0=0、trf0=0 (qtwiss.f:87)。
  解析ヤコビアン trans1 は Fortran のとおり。
- 要素の末尾 (dhg≠0): **rg2 = gammab(l)/gammab(l+1)** を trans の px, py, δ の行と cod(2),(4) に掛け、cod(6)=(cod(6)+1)rg2−1
  (運動量の規格化を新しい基準に移す)。pgev=gammab(l+1)·amass → tphyzp (以後の p0, h0)。
- **qcopymatg**: 4×5 (x px y py; δ 列) は gammab(k)≠gammab(k+1) なら sqrt(gammab(k+1)/gammab(k)) 倍 (断熱減衰を外してユニモジュラに)。
- **qmat2twiss** (非結合の枝): trpt なら入口の分散 ex epx ey epy に gr = gammab(l)/gammab(l+1) を掛けてから伝搬 (結合の枝では掛けない)。
- 各要素は軌道の cod(6) を運動量のずれに使う (加速の後は cod(6) が DP0 から 1e-8 程度ずれる)。

pysad:
- `trptacc.gammab_profile` (tfacc)、`trptacc.acc_step(el, cod, gb0, gb1, …)` → (6×6 ヤコビアン, 出口 cod): emit6d の `_State`/`_Phys`/`_tdrife0`
  を使い、tcavfrie と加速キックは解析ヤコビアン (`_cavfrie_piece`, `_acc_kick_piece`)、末尾の rg2 のブロック。CAVI と MULT(VOLT, K=0) を
  同じ関数で扱う (ndiv の eps と上限が違うだけ)。MULT の VOLT と K_n/ANGLE/BZ の併用、加速要素の DX/DY/ROTATE は NotImplementedError。
- `lattice.ins_twiss(..., acc=)` → `_ins_twiss_acc`: 要素を順に進め (加速要素は acc_step、それ以外は `_el_cod_step(el, cod, cod[5])`)、
  加速要素の行列は sqrt(gb1/gb0) 倍、`_transport_twiss(tw, M, gr=gb0/gb1)` (非結合の枝で分散に掛ける)。
- `FFS._trpt_acc()`: TRPT かつ RFSW のとき `trptacc.make_acc` (p0 = MOMENTUM/MASS、CHARGE、omega0 は全長から)。`_ins_rows` が渡す。
  RING (TRPT OFF) では None → 従来どおり空洞はドリフト。
- `Cavity`/`Multipole` の VOLT/PHI/FREQ/HARM は parser が属性に付けている (MULT は `_build_element` の ('ANGLE','E1','E2','VOLT','FREQ','PHI','HARM','BZ'))。

## 3. MARK の OFFSET (qcell.f:318-345, tffs.f tffsmarkoffset / tffselmoffset)

OFFSET≠0 の MARK は位置 xp = OFFSET + l (逆向きなら (OFFSET−.5)(−1)+.5) の Twiss を自分の行にする: lx=int(xp)、fr=xp−lx、
fr=0 なら twiss(l)=twiss(lx)、そうでなければ qtwissfrac1(lx, fr)。xp が MARK を指すならその MARK の OFFSET を辿る (最大 256 回、?Recursive OFFSET)。
linac-pfbt.sad では 5 個 (OFFSET=−.5 が 4 個: 直前の QUAD の中央 = 対称点、−0.1618 が 1 個)。
pysad: `FFS._mark_offsets(rows)` (calc の 4D 経路の末尾、CELL/INS 両方)。`row_frac(k, fr, rows=)` に rows 引数を足した。
s は `_rows_pos` が tfgeo の pos (geo.py: OFFSET MARK の位置と BEND/MULT のフリンジの長さ補正 `L − ((φ·fb1)²+(φ·fb2)²)/L/48·sin(…)/sin(φ/2)`)
に置き換える (OFFSET MARK か FRINGE≠0 かつ F1/FB1/FB2≠0 の BEND/MULT があるときだけ。linac-pfbt の BMR01 (F1=.06 FRINGE=1) で 2.1e-5 m)。
6D (CALC6D) の経路にはまだ入れていない。

## 4. ついでに直した pysad 自身の不具合

- **K1=0 の QUAD と K が全部 0 の MULT が軌道のずれた所でドリフトになっていなかった**: `_tquade_cod` の K1=0 の枝と `_mult_sequence` の
  nmmax<0 の枝が `_tqente_cod` (tsolque の中の非線形部分だけ: x += r·px の r は px² 程度) を使っていた。SAD は tdrife0 (厳密ドリフト)。
  設計軌道上 (cod=0, δ=0) は transfer_matrix の近道で正しかったので回帰では出なかった。`_exact_drift_cod` に置き換え、
  `_el_cod_step` は K 無しの MULT (`_mult_is_drift`) を DRIFT と同じ枝 (解析行列) で扱う (数値微分をやめた)。
- DISP の Value 欄: SAD (tfdisp.f:293-302) は既定キーのある要素なら値が 0 でも tdtrimz (` .0`)、MARK や既定キーの無い要素だけ ` 0`。
  pysad は値が 0 なら ` 0` にしていた (L=0 の DRIFT、ANGLE=0 の BEND、K1=0 の QUAD で違っていた)。

## 5. 検証

| 項目 | 結果 |
|---|---|
| pf2 disp_bt_r (linac-pfbt.sad、TRPT、CAVI 2 本 + 加速管 MULT 130 本、OFFSET MARK 5 個、3546 行) | 表の **13 列すべてが印字精度で SAD と一致** (<work-perf> の cmp。EX/EPX/EY/EPY の 1e-12 以下と DetR の 1e-33 の雑音を除く)。条件表の行も同じ。0.8 s |
| <work-perf>/linac_tiny.sad (20 MeV、CAVI 8 MV、MULT 12 MV、K=0 の MULT/QUAD、F1 付き BEND、OFFSET=−.5/1.25 の MARK、DX=1e-4 の軌道) | 要素ごとの DX..DDP, AX..BY が B1 (F1 付き BEND、δ=1.5e-8 の軌道) の前まで 1e-10 で一致、その後 1e-9〜3e-8 (BEND のフリンジの経路、追わない)。OFFSET MARK の s (4.8, 6.825) と Twiss が一致。DZ は 1e-11 の差 (tdrife0 の dv 項: 追わない) |

テスト: `pysad/tests/test_trpt_acc.py` (6 件: USE の大文字小文字、加速の SAD 値、OFFSET MARK、RING では加速なし、K=0 の QUAD/MULT のドリフト、DISP の Value 欄)。
pytest: 256 passed (18 s)。回帰 `QUICK=1`: FAIL 0 (50 ベンチ、113 s)。

## 6. 変更したファイル

- `pysad/trptacc.py` (新規)。
- `pysad/parser.py`: `Env.resolve_beamline`、build_lattice。`pysad/sadelem.py`: ExtractBeamLine。`pysad/maintrack.py`: TRACK USE= の比較。
- `pysad/lattice.py`: `_transport_twiss(gr=)`、`ins_twiss(acc=)`/`_ins_twiss_init`/`_ins_twiss_acc`、`_mult_is_drift`、`_el_cod_step` の DRIFT 枝、
  `_tquade_cod` K1=0、`_mult_sequence` nmmax<0、`_dv_drift`。
- `pysad/ffscmd.py`: USE/use()、`_trpt_acc`、`_mark_offsets`/`_mark_offset`/`_elm_offset`/`_rows_pos`、`_mark_offset_value`、calc の末尾、
  `row_frac(rows=)`、`_disp_show` の Value 欄、EMIT USE= の比較。
- `<work-perf>/`: patch8_usecase.py、patch9_trpt.py (差し替えの記録)、linac_tiny.sad / linac_tiny2.sad (+ .sad.out = SAD、.pysad.out)、regression_quick5/6.log、pytest_full5.log。
  <work-pf2>/disp_bt_r.pysad.out は一致した出力に更新。

## 7. 裸の要素名シンボル (disp_q_K1_r、同日)

`qlist = {Q1A,Q2A,…}; Print[LINE["K1", qlist[i]]]` の Q1A は未定義のシンボル。SAD の tflinestk は tfgetstrns で文字列でもシンボルでも名前を取る
(tftwiss.f: LINE / Twiss / Element に共通)。pysad は `_resolve_loc` (LINE/Twiss の位置引数) と `fn_element` が `Sym` を名前にする。
ついでに: v31_24 の `TEST = (-N2 2*NCELL)` は TEST が要素の型名 (icTEST) と同じなので `_parse_line_block` が LINE ブロックを終えてしまい、
`[pysad parser] Warning line 222: bad operand type for unary -` が出て TEST が登録されなかった → 次のトークンが `=` なら LINE の続きとして読む。
結果: K1 14 個と K2×2 2 個が SAD と同じ (<work-pf2>/disp_q_K1_r.pysad.out)。テスト `pysad/tests/test_symbol_names.py`。

## 8. FITP の fam 点 = MatchingAmplitude / NFAMP の有限振幅ファミリー (同日、ユーザー指示「自己検証で十分」)

SAD (tffamsetup.f, tffsa.f:1183, tffscalc.f:147-215, qdcell.f tdfun, qcell.f qcell1 の fam):
- `MatchingAmplitude = {{dp, xamp, yamp}, …}` と `NFAMP` (≤12) があると、GO/CALC のたび (tffsa) に tffamsetup が FITP の各 dp 点 k
  (dpm(k) ≤ dp < dpm(k+1)、dpm は隣り合う dp の中点) のまわりに、先頭 MARK の Twiss (twmov(1,·)) の正規モード振幅
  x0 = sqrt(βx em) (em = |EMITX|+|EMITY|)、位相 θ = 2π(i−1)/NFAMP の軌道 dx(i) = (cos θ x0, sin θ px0 − α/β x, …) を
  tfnormaltophysical (R1-R4) で物理座標にして xamp/yamp 倍した軌道のずれ dfam を NFAMP 本ずつ (x と y) 作る。
  添字は nfr+1, −(nfr+1), nfr+2, … と交互 (nfa=merge(−nfa, −nfa+1, nfa≥0))、jfam = 親の k、kfam = 1..NFAMP (x) / −1..−NFAMP (y)、dp は親と同じ。
- tffscalc: fam 点は親 i2=jfam の入口の軌道 + dfam を初期軌道に qcell1(fam=.true.) → qcod は codfnd=.true. で呼ばれるので閉軌道を探さず、
  与えた軌道のまわりの一周行列で周期解を作る (INTRES/HALFRES でも不安定にせず 4/trace に置き換え)。親が over なら utwiss をコピー。
  fam 点の stabx/staby は常に真 (`optstat(ix)%stabx = fam .or. …`)。
- tdfun: 条件は通常の点 −mp..mp と「親 jfam がその範囲にある fam 点」で評価。軌道の条件 (DX..DPY) は fam 点 (|idp|>nfr) では評価しない。
  重みは idp≠0 の式 (OffMomentumWeight)。FitValue の {id, dp} の id は kfam。
- 1 点の条件 (nfr=0) では dpm の範囲が空なので fam 点はできない (FITP>1 が前提)。

pysad: `Matcher._fam_define` (tffamsetup)、`_fam_points` (tffscalc の fam ループ)、`_eval_idps` / `_fam_dps` (tdfun)、`_normal_to_physical`、
`lattice._qcod_periodic(fixed=)` / `Lattice.twiss(fixed_cod=)` / `twiss_table(fixed_cod=)` (閉軌道の初期値・先頭 MARK への書き戻し・memo に触れない)、
`FFS._fam_setup` (MatchingAmplitude, NFAMP, EMITX/EMITY, 先頭 MARK の Twiss)、`_matcher_optics_fixed` (CELL は周期解、INS は `_ins_rows(cod0=)`)、
FitValue の id。ファミリーの定義は GO/CALC ごとに 1 回 (`_fam_def_done`)。

自己検証 (`pysad/tests/test_fam_points.py`、6 件): 本数と添字・kfam・jfam・dfam の式、NFAMP≤0 と nfr=0 と形の誤りの扱い、fam 点の入口の軌道 =
親の閉軌道 + dfam、六極ありで振幅依存のチューン、残差に fam 点が入る、fam 点の計算が閉軌道の初期値を変えない、周期解のチューンが
その軌道の一周行列の正規モードの跡と一致 (y のずれは BEND で x-y 結合を生む)、INS は軌道だけずれる、FitValue の id に kfam。
SAD(WSL) は FITP>1 の GO で再起動するので SAD との照合はない。

## 9. CAVI/SOL の ROTATE (同日)

残課題「CAVI/TCAV/MONI/SOL の ROTATE」(10/16 の ROTATE 全面対応で残したもの)。調べると:
- ユーザーのラティスの CAVI は VOLT/HARM だけ (DX/DY/ROTATE なし)、TCAV・SOL は無く、MONI に ROTATE もない。
- SAD: tcave.f (emit / CALC6D) は入口・出口で tchge(dx,dy,0,theta,…) (回してからずらす)、tmulte.f (VOLT 付き MULT) は tsolrote (ずらしてから回す)。
  追跡の tcav.f は dx,dy,theta を受け取るが**使わない** (無視)。pysad の track6 もそのまま。
- **実装**: emit6d `_tcave` と trptacc `acc_step` に tchge (CAVI) / tsolrote (MULT、行列は数値微分) を入れた。ROTATE だけなら軸対称な空洞では
  結果は変わらず、DX/DY と組むと tcavfrie の集束 (x·dpt) が回った座標でかかって x-y 結合を生む。
- **SAD との一致**: <work-perf>/cavrot.sad (CAVI L=0.5 VOLT=2 MV HARM=24 DX=2e-3 DY=−1e-3 ROTATE=30°、CELL CALC6D + emit): Twiss["DX".."R2"] と
  EMITX/EMITY/EMITZ/SIGE/SIGZ が 1e-11 で一致 (EMITY=1.74e-7、Beam tilt .2927 rad が結合の証拠。外すと EMITY=0)。
  <work-perf>/linac_rot.sad (TRPT の CAVI と加速管 MULT に DX/DY/ROTATE): 要素ごとの軌道と Twiss が B1 の前まで 1e-9 (DZ は 1.5e-10)。
- 残り (使用なし): SOL の ROTATE (SOL 自体の幾何 tsgeo が未移植)、TCAV (要素は読むが CALC/追跡の写像が未実装)、MONI の ROTATE (読み値の回転)。

テスト: `pysad/tests/test_cavrot.py` (2 件)。

## 10. ANGLE=0 K1≠0 BEND の追跡、TrackParticles の粒子損失、LINE["GAMMABETA"]、plot_inj_r (同日、残課題 A)

- **ANGLE=0, K1≠0 の BEND の追跡** (tbend.f:566-590): K0=0 なら tquad(al, K1, bz=0, dx, dy, theta+dtheta, theta+dtheta, krad, chro=T, fringe,
  f1in=f2in=f1out=f2out=0, mfring=0, eps)、K0≠0 なら tmulti (akm(0)=K0, akm(1)=K1、theta+dtheta)。pysad: `track6._bend_k1_proxy` が QUAD / MULT の
  代理を作って `_tquad_ens` / `_tmulti_ens` に渡す (track6n は K0=0 を K_QUAD、K0≠0 は track6 へ)。
- **SAD K64 の癖 (11_sad_quirks A38)**: tbend の tquad 呼び出しは theta2 に akang (K1<0 の π/2) を含めないので、K1<0 の BEND は |K1| の
  集束四極として追跡される (QUAD K1=+0.5 と全桁同じ)。CALC は正しい。pysad は |K1| で再現 (`test_bend_k1_sign_is_lost_only_in_tracking`)。
  K0≠0 かつ K1≠0 は SAD が CALC で落ちる (A27) ので追跡の照合はない。
- **TrackParticles の粒子損失** (tftrack.f / tturn.f tapert): 以前は 1 周で非有限になると FloatingPointError で止まっていた。`track6n.TrackRun` が
  tturn を区間ごとに呼び、tapert (LOSSAMPL=1, LOSSDZ=100、最初の要素・15 要素ごと・周末) で落ちた粒子を末尾へ寄せて追跡から外し、結果は粒子の
  番号順に alive (落ちたら 0) と落ちた時点の座標 (limittrack: |x|,|y|≤3·LOSSAMPL、|px|,|py|≤0.9999、|z|≤1e10) で返す。初期フラグ 0 は追跡しない。
  numba の表に無い要素 (SOL 等) は従来の track6 (損失なし)。
- **LINE["GAMMABETA"] / ["GAMMA"]** (tftwiss.f:837): gammab(k) (p/mc、TRPT の加速分布込み、分数位置は補間)。GAMMABETA は未実装で 0 を返していた。
- **plot_inj_r**: 原因は `emix = gemx/LINE["GAMMABETA","^^^"]` が INF → 初期ビーム INF → 追跡が発散。両方を直して完走 (0.5 s)。
  1000 粒子の PAR/INI (GaussRandom の既定シードの乱数まで) と生存数 1000 が SAD と同じ、XMEAN −.0118282 と −.0118276 は FLUC の乱数の差 (da_xy_r と同じ既知)。
  決定論的な 5 粒子 (<work-perf>/plot_inj_trk.sad: x=30 mm 注入 + KK3/KK4 のバンプ) は 1 周 (全粒子生存)、3 周 (注入粒子 3 個が落ちる、座標と alive)、
  RAD FLUC + SeedRandom[17] の 1 周とも SAD と全桁同じ。ListPlot/Tk 以降は対象外 (SAD は core dump)。

| 照合 | 結果 |
|---|---|
| <work-perf>/bendk1_{plain,rot,dxy,drot}.sad (BEND ANGLE=0 K1=−0.5 に ROTATE/DX,DY/DROTATE、NORAD 3 周と RAD 1 周、5 粒子) | 全変種・全粒子が SAD と同一 (生存粒子は 12 桁、落ちた粒子は抑えた座標 3/.9999/1e10 と −.999999 まで同じ) |
| <work-perf>/bendk1_quad.sad (等価な回した・ずらした QUAD) | 同一 |
| pf2 plot_inj_trk / plot_inj_cnt | 上記 |

テスト: `pysad/tests/test_bendk1_track.py` (2 件)、`test_trpt_acc.py` に GAMMABETA (1 件)。

## 11. SADHelp 照合で未対応だった物理関数・変数の優先度 高・中 (同日、docs/15_sadhelp_purpose_policy.md の続き)

ユーザー指示「物理関数で未対応のもので対応した方がよいものをリストアップして」→「優先度高、中を順に進めて」。SAD 側の定義を読んで pysad に入れ、
<work-perf>/hiprio*.sad (RING: 7 セル + RFC、QD に ROTATE=3 DEG、MOMENTUM=1GEV) で SAD K64 と照合した。新規モジュール `pysad/optics_fn.py`。

| # | 項目 | SAD 側 | pysad | SAD との照合 |
|---|---|---|---|---|
| 高1 | **TransferMatrix[p1, p2, Fold->, Symplectic->]** | Optics.n (.n → Python で書き直し): 要素行列を p1..p2 で積む。p1>p2 は Fold (周回) か逆行列、CALC6D なら 6×6 | `optics_fn.transfer_matrix` (lattice の el_mats を積む、名前・番号・Fold・Symplectic) | 4D 4 例と 6D 1 例が 1e-9 |
| 高2 | **BeamMatrix[p]** | Optics.n: NormalCoordinates の r と Diag(emx, emx, emy, emy, σz², dp²) で r.D.rᵀ。GAUSS off (UNIFORM) なら dp² → dp²/3 | `optics_fn.beam_matrix` (4D は EMITX/EMITY/DP、6D は emit の 6×6) | 4D (UNIFORM, GAUSS) と 6D が 1e-9 |
| 高2' | **NormalCoordinates[p]** | Optics.n: 4D は R 行列と Twiss から r (4×4 + 分散)、6D は etwiss2ri (NormalCoordinate6) | `optics_fn.normal_coordinates4/6` (後者は calc6d.etwiss2ri) | 4D/6D が 1e-9 |
| 高3 | **SurvivedParticles[beam]** | Fortran tfinitn.f / tftrack.f: フラグ 1 の粒子だけ残した {x,px,y,py,z,dp,flag} | `optics_fn.survived_particles` | 5 粒子 → 3 粒子、座標 1e-9 |
| 高4 | **EMIT の結果変数** DTSYNCH / PHICAV / EFFVCRATIO / EFFVC / EFFRFFREQ / OMEGA0 | tturne.f:60-64 (emit の後に trf0, vceff, wrfeff を書き出す)、OMEGA0 = 2πc/周長 (USE 時と emit 時) | `run_emit` が emit6d の戻り値 (`wrfeff` を追加) から書き出す、`_set_omega0` | 6 値が 1e-8 |
| 高5 | **EMITXE / EMITYE / EMITZE** | temit.f:456-467: sqr2 の反復の後に emitn を max(emit, EMIT?E) で下から押さえる | `emit6d.tecalc(emit_floor=)` | EMITYE=1e-7 で EMITY=1.742e-7 (SAD と 1e-8) |
| 高6 | **OpticsProlog / OpticsEpilog** | tffscalc.f:64 (変数を要素に入れた直後) / :398 (戻る直前) で定義されていれば評価する | `Matcher.optics_prolog/epilog` フック (`I.run(name, 'script')`) | PROLOG/EPILOG の表示と Twiss の値が SAD と同じ |
| 中7 | **SynchroBetaEmittance** | emit.n の定義は SAD K64 にない (tfinitn id 1038 の Fortran 実体も見つからず)、SADHelp にも本文なし | **保留** (実体不明) | — |
| 中8 | **AccelerateParticles[beam, p(n), {n, nmax}]** | dapert.n (HoldAll): 各周で MOMENTUM を p(n) にして TrackParticles、MOMENTUM は変えたまま返す | `optics_fn.OpticsFunctions.fn_accelerate_particles` (`sadscript._HOLD` に追加) | 生存粒子は 6e-10、落ちた粒子の座標は SAD と異なる (下) |
| 中9 | **OrbitGeo[p] / GeoBase[{χ1,χ2,χ3}] / SymplecticJ[n]** | Optics.n の tforbitgeo (DISP OG の値)、回転角→3×3 基底、emit.n の J 行列 | `optics_fn.orbit_geo / geo_base / symplectic_j` (SymplecticInverse も) | 1e-9 |
| 中10 | **FFS$SHOW[]** (= FFS["SHOW"]) | Fortran tffs.f: 現在の fit 条件を {{n1, n2, fun, goal, \|mfitp\|-1, scale}, …} の平らなリストで返す | `ffs_show_list` (Matcher.conditions から) | 4 条件が同じ |
| 中11 | **InitialOrbits / ElementValues** | tffsa.f | **保留** (ユーザースクリプトに使用なし) | — |
| 中12 | **ParabolaRandom[n]** | C tfRandom_.c: 一様乱数 2 個の和から放物線分布 | `sadrandom.parabola` | SeedRandom[5] の 6 個が 1e-11 |

ついでに直したもの:
- **`?FLAG`** (`g = ?GAUSS;`) をパーサに追加 (`FFS$FlagQ`)。これまでは「unexpected token '?'」。
- **SADScript の単位付き数値** `1GEV`, `30DEG`, `3MM`: 字句解析で単位を掛ける (1GEV → 1e9)。**SAD 自身は `1GEV` を `1*GEV` と読み、GEV は未定義シンボルなので
  `MOMENTUM = 1GEV` は MOMENTUM=1 になる** (hiprio.sad で確認、11_sad_quirks C24)。pysad は意図どおり 1e9 にする (意図的な拡張、SAD 互換が要るなら `1E9` と書く)。
- **EMIT コマンドは GAUSS フラグを ON にする** (SAD で `emit; ?GAUSS` が 1。Emittance[] は変えない)。pysad の EMIT も立てるようにした (11_sad_quirks A39)。
  これがないと EMIT 後の BeamMatrix の dp² が 1/3 になって SAD と合わなかった。
- Print の数値書式は SAD (`$FORM` 依存の固定小数) と pysad (repr) で違うので、テストは値で照合する。

AccelerateParticles の落ちた粒子: SAD は周ごとに TrackParticles を呼び直すので、前の周で落ちた粒子 (flag 0) も次の周の入力に入り、tapert が抑えた座標
(±3·LOSSAMPL 等) から再び追跡される。pysad も同じ手順だが、落ちた粒子の座標の内部表現が違う (5e-2 の差) ので、生存粒子だけ照合した。

| 照合 | 結果 |
|---|---|
| <work-perf>/hiprio.sad〜hiprio11.sad (TM1-4, BM1, NC, GB, OG1, SJ, OM, EV, EYE, SHOW, SV, PR, TM6, BM6, NC6, AP) | 上の表。AP は生存粒子のみ |

テスト: `pysad/tests/test_optics_fn.py` (8 件)。
