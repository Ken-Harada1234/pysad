# v2 の変更 (担当 C2): 空間電荷の続き — MULT の RADIUS (SPAC) と WSPAC

2026-11-17。前任 C の `v2_changes_C_spch.md` (SPAC の DRIFT、PSPAC + SPCH) が「残した項目」を仕上げたもの。
対象: `pysad/spch.py` (追加)、`pysad/track6n.py` の `turn` (3 行)、`pysad/ffscmd.py` の `TrackParticles` (条件 1 行)、`pysad/parser.py` (1 語)、
試験 `pysad/tests/test_v2_C2_spch.py` (35 本、約 60 秒)、例 `example/v2_C/C2s_mult_radius.sad`・`C2s_wspac.sad`、作り直した例 `example/05_tracking/flag_SPAC.*`。
どれも v2 モードだけで動く。**k64 モードは今までどおり、空間電荷を追跡に入れない** (親の決定。§5)。

値の根拠は gs (K64) ではなく、解析式と独立な数値計算 (包絡方程式の数値積分、場の積分表示の数値積分、自前の 2×2 行列のリング)。gs の値は「K64 は X」として並べる。

## 1. MULT の RADIUS の空間電荷 (SPAC)

### K64 の動作 (tmulti.f)
- `SPAC` が ON で MULT の `RADIUS≠0`・`L≠0` のとき `spac1` が立ち、分割を `ndiv = max(ndiv, nint(|L|/(alstep·EPS)), nint(|bz L|/1.5/EPS))` に増やす (`alstep` = 5 cm、`EPS` = 0 は 1)。
- 分割のループは「本体 (四極 + 双極) の半分 → [多極 (n≥2) のキック、`spkick(al*wm)`] → 本体 → … → 本体の半分」。最後に `tapert(radius, radius)`。
- **`spkick` の呼び出しが壊れている** (SPC3、tmulti.f:222・541): `call spkick(np,x,px,y,py,z,g,dv,sx,sy,sz,al*wm,radius,alx,kturn,kptbl)` は `sx,sy,sz` を余分に渡しており、`spkick` の引数 (spkick.f) と合わない。
- もう 1 つ: `spkick` に渡す長さ `al*wm` は、最後の回だけ `wm = wi/2` なので、合計が L の (ndiv − 1/2)/ndiv 倍になる (多極のキックは `akrm = akr*wi` で ndiv 回とも全量。空間電荷だけ最後が半分)。

### v2 の動作 (`spch.SpaceCharge._mult`)
- 対象: `SPAC` ON、`MULT` の `RADIUS≠0`、`ANGLE = 0`、`VOLT + DVOLT = 0`、SOL 領域の外。強さ 0 の MULT (ドリフトとして追う) は DRIFT と同じ経路 (ビット単位で同じ結果)。
- 分割: `ndiv = max(MULT の分割, nint(|L|/(SPACSTEP·EPS)), 放射の分割)` (上限 1000。`SPACSTEP` は既定 0.05 m)。
- 順序: 入口 (DX・DY・ROTATE、フリンジ) → 本体の半分 → [多極 (n≥2) のキック、半径 |RADIUS| の絞り、**空間電荷のキック (長さ L/ndiv)**] を ndiv 回 (間に本体 1 つ分) → 本体の半分 → 出口のフリンジ → 絞り → 出口。空間電荷のキックは ndiv 回とも L/ndiv (合計 L)。
- 場は DRIFT と同じ `Solver3D`・`kick3d` (バンチの静止系の 3 次元のポアソン方程式、円筒の壁)。座標は MULT の中の系 (DX・DY・ROTATE の後) なので、**真空ダクトの軸は磁石の軸**。
- 本体・多極のキック・フリンジ・放射は、今までの MULT のカーネル (`track6n_ext` の `mult_qu_all`・`mult_kick_all`・`mult_fringe_all`、`track6n.tradk_all`) をそのまま呼ぶ (カーネルは変更していない)。
- `RADIUS<0` は絞りだけ (DRIFT と同じ)。`NORFSW` で `RADIUS>0` はエラー (DRIFT と同じ)。
- 入れられない MULT (`ANGLE≠0`、`VOLT≠0`、SOL 領域の中) は、メッセージ ` SPAC: 次の MULT は RADIUS≠0 だが空間電荷を入れない (…): 名前 (理由)` を出して、空間電荷なしで追う。
- `Element["RADIUS", "名前"]`・`LINE["RADIUS", k]` が MULT でも定義の値を返す (`parser.py`。今までは 0。モードによらない)。

### 根拠と試験

| 試験 | 比べる相手 | 条件 | 一致 |
|---|---|---|---|
| `test_mult_k0_equals_drift` | 同じ長さ・同じ RADIUS の DRIFT (C が検証済み) | 5 MeV、1e11 個、1 m | ビット単位 |
| `test_mult_zero_charge_equals_plain_mult` | PBUNCH = 0 で、RADIUS の無い同じ MULT。六極・スキュー四極・DX・ROTATE・放射 (RAD NOFLUC) つきも | SPACSTEP を大きくして分割を同じにする / 既定 | 1e-15 / 四極だけ 1e-8、六極つき 1e-6 (分割の違い) |
| `test_mult_quad_kv_envelope` | 四極の収束力つきの KV の包絡方程式 X'' = −kX + 2K/(X+Y) + ε²/X³、Y'' = +kY + … の数値積分 (X = 2.1214 mm、Y = 6.8133 mm。K = 0 では 1.1510、4.2997 mm) | k = 3 /m²、L = 1 m、K = 5.98e-6、ε = 2e-6 m、4096 粒子 | 2·rms 0.07%・0.2% (0.5% で判定)、傾き 1% |
| `test_mult_kick_step_convergence` | SPACSTEP = 0.25、0.125、0.0625、0.03125 m の結果の差 | 同上 | 差が 1/4 ずつ (2 次)。最後の差 < 0.2% |
| `test_mult_negative_radius_and_unsupported` | RADIUS<0 は絞りだけ。ANGLE≠0 はメッセージ | | 厳密 |
| `test_ring_model_matches_laslett_formula` | 線形レンズのモデル (下) と 1 次の摂動の式 Δν = −(1/4π)∮βk ds | 4 セルの FODO、ν = (0.7856, 0.6468) | 1%。ドリフトだけなら全周の 80% |
| `test_ring_tune_shift_mult_full_circumference` | 自前の 2×2 行列のリングに、全周 1.25 cm ごとに線形レンズ kp/(σ(σx+σy)) を入れたチューン | 5 MeV、周長 8 m、εx = 1e-6、εy = 0.5e-6 m、σz = 5 cm、2.5e8 個、4096 粒子 | MULT: 追跡/モデル = 0.967 (x)、0.942 (y)。QUAD (ドリフトだけ) はドリフトだけのモデルに合い、MULT の 0.80 倍 (7% で判定) |
| `test_three_models_agree_on_tune_shift` | 同じビームで SPAC (MULT、全周)・PSPAC・WSPAC | 同上 | 下の表 |
| `test_mult_radius_k64_mode_unchanged` | k64 モードで SPAC・WSPAC が ON でも NOSPAC とビット単位で同じ | | 厳密 |
| `test_mult_radius_readable` | `Element["RADIUS", ..]`、`LINE["RADIUS", k]` | | 厳密 |

**リングで「DRIFT だけでは全周の 80% しか入らない」問題**: 四極を `MULT` (RADIUS = 3 cm) で書くと全周に入る。4 セルの FODO (周長 8 m、ドリフト 6.4 m) のチューンシフト (追跡は 1 周の行列の対角和から):

| 方法 | Δνx | Δνy | 線形レンズのモデル (同じ場所にレンズ) |
|---|---|---|---|
| SPAC、四極が MULT (全周) | −0.004135 | −0.006314 | −0.004275、−0.006706 (全周) |
| SPAC、四極が QUAD (ドリフトだけ) | −0.003308 | −0.005062 | −0.003421、−0.005354 (ドリフトだけ = 全周の 80%) |
| PSPAC + SPCH (全周の分を SPCH にまとめる) | −0.004146 | −0.006789 | −0.004263、−0.006725 |
| WSPAC (要素の入口) | −0.004263 | −0.006726 | −0.004263、−0.006726 |
| Laslett の式 (1 次の摂動、全周) | −0.004264 | −0.006672 | |

SPAC の追跡/モデルが 0.94〜0.97 なのは、マクロ粒子の雑音とメッシュ (C の記録と同じ傾向): 16384 粒子では 0.982、0.976、メッシュを 48×48 にすると 0.96〜0.98。中心の勾配 (チューンシフト) は、ビームの中心付近の少数の粒子で決まるので雑音に敏感で、分布によっては大きめにも出る (準乱数 4000 個で +5〜6%、1000 個で +10〜15%: `flag_SPAC` の (B))。

### 使える範囲と条件
- `TrackParticles`、電子 5〜20 MeV、1e8〜1e12 個/バンチ、ビームの大きさ 0.5〜7 mm、ダクト 2〜3 cm で確かめた (式は DRIFT と同じなので、C の §6 の範囲がそのまま当てはまる)。
- MULT は K0 (双極のキック)・K1・SK1・多極 (K2 以上)・DX・DY・ROTATE・フリンジ・放射つきで走る (PBUNCH = 0 で空間電荷なしの MULT と一致することを確かめた)。
- キックの間隔は `SPACSTEP` (既定 5 cm)。短い四極 (5 cm 以下) では 1 キック。2 次で収束する。
- 条件は DRIFT と同じ: `RFSW` ON、粒子に z の広がり、ビームの外縁が 0.9·RADIUS 以内、1000 粒子以上 (チューンシフトを数 % で見るなら 4000 粒子以上)。

### 使えない範囲
- `ANGLE≠0` の MULT (曲がった軌道の中の空間電荷。カーネルが別 (`multa`) で、分割の間にキックを挟む形にしていない。直線の 3 次元の場をそのまま使ってよいかも未確認)、`VOLT≠0` の MULT (`tmultiacc`: 加速で β γ が途中で変わる)、SOL 領域の中の MULT (`K_SOL_MULT`)。いずれもメッセージを出して空間電荷なしで追う。
- `QUAD`・`BEND`・`SEXT` などの MULT 以外の磁石の中 (K64 にも無い)。偏向磁石のあるリングでは、SPAC は全周には入らない (全周に入れるなら PSPAC + SPCH か WSPAC)。
- 放射 (RAD) と空間電荷を同時に入れた MULT は、走るが値の確認はしていない (PBUNCH = 0 の一致だけ)。FLUC つきは乱数を引く順が MULT の分割で変わる。
- K64 との一致 (K64 は呼び出しが壊れているので比べる相手が無い)。

## 2. WSPAC (ビーム行列のガウス分布の力。強-弱のモデル)

### K64 の動作 (ソースを読んだ結果)
- SADHelp: 「WSPAC: 強-弱のモデルの空間電荷。ビームの大きさは EMIT (CODPLOT つき) で求めておく。全方向ガウス分布、粒子数は PBUNCH。光学・エミッタンス・追跡に効く」。「SELFCOD: WSPAC の追跡で、力の中心を、その時の粒子の重心にする (既定は EMIT の閉軌道)」。
- 追跡 (tturn.f:270-278、584-591): 各要素 l を追う前に、`sspac = (pos(l)+pos(l+1))/2` (要素の中点) として、長さ `sspac − sspac0` の分を `twspac(…, cod = twiss(l, DX..DDP), beamsize(:,l), l)` で蹴り、`sspac0 = sspac`。区間の終わりに `pos(lend) − sspac0` の分を `beamsize(:,lend)` で蹴る。つまり**キックは要素の入口で、入口のビーム行列を使い、長さは「前の要素の中点からこの要素の中点まで」**。
  - `sspac0` の初期値は `pos(lbegin−1)` (区間の始めの 1 つ前の要素の入口。`lbegin = 1` では `pos(0)`)。
  - SOL 領域の中 (tsol) では呼ばれない (領域の長さの分は、領域の次の要素の入口でまとめて入る)。
  - `wspaccheck`: ビーム行列が無い (`ifsize = 0` か 6 次元でない) と `WSPAC without beam matrix!  You need EMIT with CODPLOT.` を出して、追跡せずに戻る。
- `twspac` (tintrb.f:378): ビーム行列 Σ から、z で条件を付けた横の分散 `xx1 = Σxx − Σxz²/Σzz`、`yy1`、`xy1` → 主軸の向き θ と σ1、σ2。粒子ごとに、中心 (`SELFCOD` なら粒子の重心、そうでなければ閉軌道 `cod`) からのずれを主軸に回し、
  `f = PBUNCH·r_c·L/(p1²·h1·√(2π σz²))` (`p1 = (1+δ_cod)·βγ`、`h1 = γ`: **係数は 1/(β²γ³) で正しい**)、`az = f·exp(−(z−zc)²/2σz²)`、
  `Δpx = −az·fx` (`fx` は Bassetti-Erskine の `bbkick1` と同じ符号: 同符号の電荷で負 → 外向きのキック)、`Δδ = −(z−zc)/σz²·az·fu` (`fu` は中心で 0 の電位 `twspu`)。
  `fx, fy, fu` は 60×60×20 の表 (`twspfu`) のスプライン補間で、表の外は `bbkick1`・`twspu` を直接呼ぶ。
- 光学と Emittance (tintrb.f:78-93 `twspace`、temit.f): 線形のレンズ `akx = 2 N r_c L/(p1² h1 (σx+σy) √(2π)σz)/σx` を転送行列に入れ、ビーム行列が収束するまで反復する。

### v2 の動作 (`spch.WSpace`、`spch.make_wspace`、`SpaceCharge.wkick`)
- `WSPAC` が ON の `TrackParticles` で、各要素の入口に K64 と同じ区切りのキックを入れる。区間の始めは区間の始点 (`pos(lbegin)`) から数える。SOL 領域の中の要素でも入口ごとに入れる。
- ビーム行列の取り方:
  1. 既定 (K64 と同じ): **Emittance の平衡のビーム行列** (各要素の入口。`ffs.beam_matrix_at(l)` = `emit_extra` の `BeamMatrices` = `LINE["SIG", l]` の CALC6D の値)。EMIT を先に実行していなくても、追跡のときに計算する (`RFSW`・`RADCOD` の今の値で)。
  2. `WSPACBEAM = {εx, εy, σz}` か `{εx, εy, σz, σδ}` (v2 で追加): Twiss (β、α、分散、R1〜R4) とエミッタンスから作る (`tfbeam` と同じ形。MINCOUP の下限は掛けない。x-z の相関は 0)。放射平衡の無い系 (低エネルギーのリング、輸送路、垂直のエミッタンスが 0 の平らなリング) で使う。
- 閉軌道 (力の中心): Twiss の `DX DPX DY DPY DZ DDP` (各要素の入口)。`SELFCOD` なら粒子の重心 (x、y、z)。
- 場: Bassetti-Erskine の式を複素誤差関数 (`scipy.special.wofz`) で直接計算する (`gauss_field`。K64 の表の補間はしない)。円 (|σx−σy| ≤ 1e-6 σx) は 2(1−e^(−r²/2σ²))/r²。電位は 1 次元の積分表示を 64 点の Gauss-Legendre で積分する (`gauss_potential`)。
- キック (L は上の区切りの長さ、r_c = q² r_e m_e/m、(x1, y1) は主軸の座標、F = ∇U は単位の線密度のガウス分布の場で遠方 2/r):

      A(z) = PBUNCH·r_c·L·exp(−(z_b − zc)²/2σz²) / (β²γ³ √(2π) σz)
      Δ(px/p0) = A·F_x (主軸から戻す)、Δ(py/p0) も同様 (外向き)
      Δδ = (z_b − zc)/σz² · A · (U_wall − U(x1, y1))

  - `z_b = β1 z/β_i` (バンチの速度での位置 = −β1 c × 到着時刻。K64 は粒子の z をそのまま使う。差は δ/γ² の程度)。こうすると、キックは (x, px, y, py, −ct, E) で厳密に電位 A(z_b)·(U_wall − U) の勾配になり、シンプレクティック。到着時刻は変えない (C の SPAC と同じ扱い)。
  - `U_wall`: 既定 0 (K64 と同じ: 軸上で電位 0。縦の力は「軸上の粒子との差」だけ)。`WSPACRADIUS = R` (v2 で追加) を与えると U_wall = U(R) (= [U(R,0)+U(0,R)]/2) で、半径 R の円筒の真空ダクトの中の軸上の電位 (g 因子 ≈ 2 ln(R/σ) − ln 2 + γ_E) が入る。
- エラー (分かりやすい文で止める): ビーム行列が使えない (`WSPAC: 要素 1 の入口のビーム行列が使えない (σz² = 0)。… WSPACBEAM = {εx, εy, σz, σδ} を与えること`)、垂直の大きさが 0 (`… ビームの大きさが 0 以下 … 垂直のエミッタンスが 0 のリングでは、結合を入れるか WSPACBEAM … を与えること`)、`WSPACBEAM` の形が違う、追跡できない要素の並び。`SPAC`・`PSPAC` も ON のときは「力が二重に入る」とメッセージを出す (止めない)。

### 線形のチューンシフトの式 (導出)
1. 2 次元のガウス分布 (線密度 λ) の中心付近の場: E_x = λ·2x/(σx(σx+σy)) (ガウス単位。Bassetti-Erskine の式の中心での展開。円なら λx/σ²)。
2. 実験室系で、電場の力と磁場の力の和は e E_x (1 − β²) = e E_x/γ²。運動方程式 γ m β² c² x'' = e E_x/γ² から x'' = k_x x、
   **k_x = 2 r_c λ/(β²γ³) · 1/(σx(σx+σy))**、k_y = 2 r_c λ/(β²γ³) · 1/(σy(σx+σy))。バンチの中心では λ = N/(√(2π)σz)。`twspace` の `akx` と同じ。
3. 1 周の行列 M (tr M = 2 cos μ) に、β の場所で薄いレンズ [[1,0],[k L,1]] を掛けると tr = 2 cos μ + β k L sin μ、つまり μ' = μ − β k L/2。
   足し合わせて **Δν_u = −(1/4π) ∮ β_u(s) k_u(s) ds** (Laslett のチューンシフト。1 次の摂動)。円いビーム (σ² = βε) では Δν = −r_c λ C/(4π β²γ³ ε) (C は周長)。

### 根拠と試験

| 試験 | 比べる相手 | 一致 |
|---|---|---|
| `test_gauss_field_and_potential_vs_integral[6 通り]` | 場と電位の 1 次元の積分表示 F_x = 2x∫_0^∞ exp(−x²/(2σx²+t) − y²/(2σy²+t))/((2σx²+t)^{3/2}(2σy²+t)^{1/2}) dt、U = ∫_0^∞ [1 − exp(…)]/√(…) dt の数値積分 (quad)。σy/σx = 0.3、3.3 (縦長)、0.999、1、0.01、0.95。点は 0.0001σ〜20σ | 1e-8 (実際は場 1e-10、電位 1e-10) |
| `test_gauss_field_limits` | 中心の勾配 2/(σx(σx+σy))、遠方 2/r、円の式、∇U = F (差分。σy/σx = 1e-3 でも) | 1e-9〜1e-5 |
| `test_wspac_kick_formula_tilted_beam` | x-y と x-z の相関のあるビーム行列で、キックが上の式 (主軸・条件つき分散・A(z)·F (quad)・縦) どおり。線形の極限が `twspace` の式 | 1e-7、縦 1e-5 |
| `test_wspac_kick_is_symplectic` | Jᵀ S J = S (差分のヤコビ行列、(x, px, y, py, −ct, E)) | 1e-6 |
| `test_wspac_longitudinal_wall_potential` | `WSPACRADIUS`: 軸上の電位 U(R) = Ein(R²/2σ²) を、SPAC の 3 次元のソルバー (Solver3D、独立) の軸上の電位/線密度と比べる。頭 (z>0) は加速、尾は減速。WSPACRADIUS = 0 では軸上の Δδ = 0 | 1% |
| `test_wspac_ring_tune_shift` | 自前の 2×2 行列のリングで、各要素の入口に同じ長さの線形レンズを入れたチューン。Laslett の式。PBUNCH に比例 | 0.1% (モデル)、3% (式)、2% |
| `test_three_models_agree_on_tune_shift` | 同じガウスビームで SPAC (−0.004135, −0.006314)、PSPAC (−0.004146, −0.006789)、WSPAC (−0.004263, −0.006726) | 互いに 8% 以内 (SPAC・PSPAC はマクロ粒子の雑音) |
| `test_wspac_emittance_beam_matrix_ring` | Emittance のビーム行列の経路。0.5 GeV の電子リング (周長 35.3 m、ν = (2.126, 1.437)、弱いスキュー四極、2e10 個)。4×4 の 1 周の行列の固有値のチューンの変化と、1 次の式 −(1/4π)Σβ k (Twiss の β、WSpace の σ・傾き) | Δνx = −3.008e-5 (式 −3.0005e-5)、Δνy = −1.1965e-3 (式 −1.1962e-3)。1% で判定 |
| `test_wspac_flat_ring_needs_vertical_size`、`test_wspac_errors_and_messages` | エラーの文 | |
| `test_wspac_selfcod_and_offset` | 中心からずれた集団: 既定は外向きの正味のキック、SELFCOD では重心のキックが 0 | 1e-3 |
| `test_example_wspac_script` | 例のスクリプトが走り、追跡と式が合う。gs の値とも 1.5% | 1% |

**gs (K64) との比較** (根拠ではなく照合): 例 `C2s_wspac.sad` の前半 (0.5 GeV のリング、2e10 個) を gs で、`CODPLOT; WSPAC; EMIT;` の後に NORAD で追跡すると Δνx = −2.938e-5、Δνy = −1.210e-3。v2 は −2.963e-5、−1.196e-3 (差 0.8%、1.2%)。K64 は EMIT がビーム行列を空間電荷込みで解き直す分だけ違う。gs で EMIT を `NOWSPAC` で行ってから `WSPAC` で追跡すると NaN になる (K64 は「WSPAC つきの EMIT」が前提)。

### 使える範囲と条件
- `TrackParticles` の中の、ガウス分布のビームが作る空間電荷の力 (強-弱: ビームの大きさは固定で、追跡する粒子は力を受けるだけ)。粒子数によらない (試験粒子 1 個でも使える)。速い (1 キック 1 ms/1000 粒子、10 ms/4096 粒子)。
- 全要素 (DRIFT・四極・偏向・空洞…) の入口で入るので、リングの全周に入る。RADIUS も SPCH 要素も要らない。
- 確かめた範囲: 電子 5 MeV (WSPACBEAM)、0.5 GeV (Emittance のビーム行列)、Δν = −3e-5〜−0.025、σy/σx = 0.02〜3、結合のあるビーム (傾き)。
- `NORFSW` でも走る (粒子の z が毎周 0 になるので、全粒子がバンチの中心の線密度を感じる)。

### 使えない範囲・限界
- **光学 (CALC) と Emittance への空間電荷** (K64 の `twspace`・`qwspac`、Emittance の反復): 未対応 (`emit6d`・`emit_extra`・`calc6d` は変更しない方針)。したがって v2 の WSPAC のビーム行列は**空間電荷なし (ゼロ電流) の平衡値**で、空間電荷によるビームの大きさの変化は入らない (上の gs との 1% の差)。Twiss のチューンにも空間電荷のシフトは出ない (追跡でだけ見える)。
- ビームの大きさは追跡の間じゅう固定 (強-弱)。エミッタンスの増大やコヒーレントな振動は扱えない (それには SPAC・PSPAC)。
- 垂直のエミッタンスが 0 のリング (結合の無い平らなリング) は Emittance の経路では使えない (エラー)。`WSPACBEAM` を与える。
- `WSPACBEAM` は σz が一定、x-z の相関なし、MINCOUP なし。加速 (TRPT) で運動量が変わる輸送路ではエミッタンスの断熱減衰を入れていない (β²γ³ は空洞の後の値を使う)。
- 縦の力は、既定では K64 と同じ「軸上で 0」。`WSPACRADIUS` の U_wall は円いビーム・円筒のダクトの近似 (扁平なビームでは軸の平均)。
- SOL 領域の中の要素でのキック、`SELFCOD` と閉軌道のあるリング、放射つきの追跡との組み合わせは、走るが値の確認はしていない。
- MAIN の `TRACK`、`DynamicApertureSurvey`、スピン追跡は対象外 (C の SPAC と同じ)。

## 3. 設定 (FFS の変数。C の表への追加)

| 変数 | 既定 | 意味 |
|---|---|---|
| `SPACSTEP` | 0.05 | SPAC のキックの間隔 [m]。DRIFT と MULT で共通 |
| `WSPACBEAM` | なし | `{εx, εy, σz}` か `{εx, εy, σz, σδ}`。与えると WSPAC のビーム行列を Twiss から作る (無ければ Emittance の平衡値) |
| `WSPACRADIUS` | 0 | WSPAC の縦の力の電位の基準 (真空ダクトの半径 [m])。0 は K64 と同じ (軸上で 0) |

フラグ: `SPAC`・`PSPAC`・`WSPAC`・`SELFCOD` (FFS の最上位に書く。`Do[...]` などの式の中では切り替わらない)。

## 4. 例

- `example/v2_C/C2s_mult_radius.sad` (19 秒): (1) 四極 (MULT) の中の KV ビームの包絡 — 追跡 2·rms = 2.109、6.779 mm、包絡方程式 2.121、6.813 mm (全粒子。バンチの端の分 0.5% 小さい)。(2) リングのチューンシフト — SPAC (MULT、全周) −0.004463、−0.007066、PSPAC −0.004343、−0.006568、WSPAC −0.004263、−0.006726、Laslett (全周) −0.004264、−0.006671、(ドリフトだけ) −0.003414、−0.005331。SPAC の +5% は準乱数 4000 個のビームの雑音。
- `example/v2_C/C2s_wspac.sad` (3 秒): 0.5 GeV の電子リングで Emittance のビーム行列を使う WSPAC。追跡 Δν = (−2.963e-5, −1.196e-3)、スクリプトの中で計算した Laslett の式 (−2.989e-5, −1.193e-3)。SELFCOD (重心のキック 0.200 nrad → 0)、平らなリングのエラーと `WSPACBEAM`。
- どちらの `*.pysad.txt` も **v2 モードの出力** (Windows、2026-11-17)。k64 モードでも走る (空間電荷なしの値になる)。
- `example/05_tracking/flag_SPAC.*` を作り直した (§6)。

`example/v2_C/README.md` に足す行の案 (親がまとめる):

    | `C2s_mult_radius.sad` | MULT の RADIUS の空間電荷 (SPAC)。四極を MULT で書くと全周に入る | 四極の中の KV ビームの包絡が包絡方程式と 0.5% で一致、リングのチューンシフトが SPAC・PSPAC・WSPAC・Laslett の式で一致 (SPAC・PSPAC は粒子の雑音で 5%) |
    | `C2s_wspac.sad` | WSPAC (ビーム行列のガウス分布の力) と SELFCOD、`WSPACBEAM` | 追跡のチューンシフトが Laslett の式と 1% で一致、gs (K64) とも 1%。平らなリングはエラーで知らせる |

`example/05_tracking/README.md` の flag_SPAC の行の案: `| [flag_SPAC](../../example/05_tracking/flag_SPAC.md) | 空間電荷 (ドリフトの包絡とリングのチューンシフトを解析式と比べる。RADIUS<0 は絞りだけ) | ○ pysad v2 は解析式と一致 (k64 モードは未対応)。gs は力が 1/γ 倍 (K64 の誤り) |`

## 5. k64 モード
親の決定: 空間電荷の k64 モードは K64 の動作のままにする (= 今までどおり、SPAC・PSPAC・WSPAC は追跡に入らない)。変えていない (`test_mult_radius_k64_mode_unchanged`、C の `test_spac_k64_mode_unchanged`)。
モードによらない変更は 1 つ: MULT の `RADIUS` が要素の属性になった (`parser.py`。`Element["RADIUS", ..]` が gs と同じく定義の値を返す。RADIUS を読む処理は spch.py だけなので、追跡・光学の値は変わらない)。

## 6. 例 `flag_SPAC` の作り直し
- 誤り: `NOSPAC;`・`SPAC;` を `Do[...]` の中に書いていて、フラグが切り替わっていなかった (式の中ではフラグの命令にならない)。そのため「gs でも効果が出ない」と書いていた。
- 直した内容: フラグを最上位に書く。`RFSW` にする (v2 は NORFSW の空間電荷をエラーにする)。ラティス `00_lattice/spac_ring.sad` は空洞の電圧を 0 にし (RFSW で 20 MeV のビームを乱さない)、1 m のドリフトの輸送路 `TL` を足した。`PBUNCH = 1e13` (Δν ≫ 1) をやめ、(A) ドリフトの包絡 (1e12) と (B) リングのチューンシフト (1.5e9、Δν ≈ −0.012) を解析式と比べる。元の RADIUS<0 の絞りの例は (C) として残した。
- 結果 (v2 の pysad / gs / 解析式): (A) rms 0.716 mm / 0.506 mm / 0.720 mm (K/γ なら 0.506 mm)。(B) Δν = (−0.01180, −0.01323) / (−1.85e-4, −2.08e-4) / (−0.01109, −0.01254) (1/γ 倍なら −2.83e-4、−3.20e-4)。(C) 生存 84/200 で一致。
- `.gs.txt`・`.pysad.txt` を `run_example.py` で作り直した (pysad は v2 モード。約 60 秒、gs は 11 秒)。`.md` に、どちらのモードの出力か、SPC1、解析式との比較、未確認・限界を書いた。試験 `test_example_flag_spac_script` が (A) と (C) を走らせる。

## 7. C が「未確認」としたもの

| 項目 | 結果 |
|---|---|
| TRPT の加速との組み合わせ | **確かめた** (`test_spac_with_trpt_acceleration`)。0.5 m のドリフト (5 MeV) → 空洞 (+5 MV、長さ 0、10 MHz) → 1 m のドリフト (10.01 MeV/c)。包絡方程式を 2 区間で積分 (空洞で傾きとエミッタンスが p1/p2 倍、2 つ目のドリフトは空洞の後の β²γ³): 2·rms = 3.848 mm (空間電荷なし 2.235 mm)、追跡 3.845 mm (0.07%)。空洞の前の運動量の K を使うと 4.829 mm なので、空洞の後の p0 を使っていることが分かる。長さのある空洞の中・加速しながらの空間電荷は入らない (空洞は DRIFT でない) |
| 電子以外の粒子 | **式は確かめた、SAD スクリプトからは使えない**。`TrackParticles` (ffscmd.fn_track) は `TrackRun` に `MASS`・`CHARGE` を渡しておらず、追跡は常に電子の質量 (空間電荷と関係のない pysad の既存の制限)。`TrackRun(amass=陽子)` と `SpaceCharge(amass=陽子)` を直接作る試験 (`test_spac_proton_direct`) で、100 MeV の陽子 (β = 0.428、γ = 1.107) の KV ビームの包絡が K = 2 r_p λ/(β²γ³) の包絡方程式と 0.5% で合う。足りないのは、fn_track が MASS・CHARGE を TrackRun に渡すこと (親の判断。他の追跡にも影響する) |
| 「小さい RADIUS>0 で gs が落ちる」 | **再現し、原因を特定した** (SPC10)。RADIUS の大きさではなく、**SPAC の中の粒子数が少ないこと**が原因。`spkick.f` の `spmesh` は縦を `nza` = 10 枚 (np < 4000) のスライスに `npz = ceil(np/nza)` 個ずつ分け、`dz2 = zz(itab(np)) − zz(itab((nza−1)·npz+1))` を計算する。np < 9·ceil(np/10)+1 では添字が np を超えて配列の外を読み、等しいと最後のスライスが 1 粒子で dz2 = 0 になる。gs (1 m のドリフト、RADIUS = 3 cm) で確かめた: np = 2、15、33 は SIGSEGV、5 は Bus error、9・10・11・19・25・28・81・82 は NaN で全粒子が失われ、1・20・30・40・50・60・90・100・120・150・200・300・2000 は正常 (規則どおり)。リング (RADIUS = +0.4 mm、200 粒子) では、絞りで粒子が減る途中でこれに当たって落ちる (NORFSW・PBUNCH = 1e9、RFSW・1e9、RFSW・1e13 の 3 通りで SIGSEGV・Bus error)。v2 は粒子数によらず走る (`test_spac_small_particle_number`: 2、5、15、33、81 個) |
| ウェーク・放射との組み合わせ | 未確認のまま。走るが、独立な根拠 (空間電荷とウェークを同時に含む解析解) を用意できていない。足りないもの: それぞれ単独で検証済みのキックが、同じ追跡で干渉しないことを示す試験 (例: 順序を入れ替えても 2 次の差に収まること) |
| SOL 領域の中の四極・MULT | 未対応のまま (空間電荷はドリフトだけ。MULT はメッセージを出す) |
| バンチの静止系で相対論的な運動 | 未確認のまま (静電場の近似の限界。γ·x' ≳ 1 で崩れる。確かめるには、遅延ポテンシャルを含む場の計算が要る) |

## 8. 変更したファイルと箇所

| ファイル | 箇所 |
|---|---|
| `pysad/spch.py` (LF) | 追加: `_ein`、`gauss_field`、`gauss_potential`、`WSpace`、`mult_unsupported`、`twiss_beam_matrix`、`make_wspace`、`SpaceCharge.wkick`・`_sckick`・`_mult`。変更: `SpaceCharge.__init__` (MULT の RADIUS、`wspace`)、`begin`、`handles`、`element`、`end`、`make` (WSPAC、MULT のメッセージ)、冒頭の説明 |
| `pysad/track6n.py` (CRLF) | `turn`: `spch.ctx = (rad, fluc, rs, tm)` (1 行)、要素の前の `spch.wkick` (2 行) |
| `pysad/ffscmd.py` (CRLF) | `fn_track`: `spch.make` を呼ぶ条件に `WSPAC` を足した (1 行と注釈) |
| `pysad/parser.py` (CRLF) | `_apply_extra_params`: Multipole の属性の一覧に `'RADIUS'` (1 語) |
| 新規 | `pysad/tests/test_v2_C2_spch.py`、`docs/history/v2_changes_C2_spch.md`、`example/v2_C/C2s_mult_radius.sad`・`.pysad.txt`、`example/v2_C/C2s_wspac.sad`・`.pysad.txt` |
| 作り直し | `example/05_tracking/flag_SPAC.sad`・`.md`・`.gs.txt`・`.pysad.txt`、`example/00_lattice/spac_ring.sad` |

`track6n_ext.py`・`lattice.py`・`emit6d.py`・`emit_extra.py`・`calc6d.py` などは変更していない (読むだけ、または既存の関数を呼ぶだけ)。

## 9. 要判断 (親へ)
1. `TrackParticles` が `MASS`・`CHARGE` を追跡に渡していない (電子以外は追えない)。空間電荷の式は粒子の種類に対応しているが、入口が無い。
2. 光学・Emittance への WSPAC (K64 の `twspace` の反復) を入れるか。入れると、WSPAC の追跡のビーム行列が自己無撞着になり、Twiss にも空間電荷のチューンシフトが出る (`emit6d`・`calc6d` の変更が要る)。
3. `docs/46` (「SPCH・SPAC は読むだけ」「gs が落ちる: SPAC の小さい RADIUS>0」)、`docs/21` (K64 のバグ: SPC1〜SPC10 を B 番号にするか)、`example/05_tracking/README.md` の flag_SPAC の行、`example/v2_C/README.md` の更新 (§4 に案)。
4. `parser.py` の 1 語 (MULT の RADIUS を属性にする) はモードによらない変更。k64 の追跡・光学の値は変わらないが、`Element["RADIUS", MULT]` の値が 0 から定義の値に変わる (gs と同じになる)。
