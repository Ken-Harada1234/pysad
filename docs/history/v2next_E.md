# v2next の記録: 担当 E (D1 の未着手 3 項目: 傾いたソレノイドの POL、ANGLE=0 の BEND・MULT の F1、脱偏極率の追跡との差、2026-11-18)

対象: `pysad/spin.py` (`tradke_spin`・`SOL_BPAR`・`SPIN_ANGLE_EXACT`)、`pysad/solemit.py`、`pysad/emit6d.py` (`_tsteee`・`_tmulte`・`_tquade`)、`pysad/lattice.py` (`_mult_params`・`_mult_sequence`・`_tbfrme_cod`・CALC の経路)、
`pysad/elkern.py` (`_tsteee` の edgefix)、`pysad/track6.py` (`_tquad_sol`)、`pysad/track6n_ext.py` (領域内の CAVI)。numba のカーネル (`elkern.py`) を編集したので `py -3.12 -m pysad.nbcache` でキャッシュを掃除した。
試験: `test_v2_E_spinsol.py` (20 件)、`test_v2_E_steer.py` (19 件)、`test_v2_E_pol_kick.py` (2 件 + `PYSAD_SLOW=1` の 1 件)。既存の試験の変更は 1 件 (§6)。
例 (`example/v2_next/`): `E_steer_f1.sad`・`E_mult_f1.sad` (gs の出力 `.gs.txt`、pysad の v2 `.pysad.txt`・k64 `.pysad_k64.txt`)、`E_b108_poles.py`、`E_sol_tilt_spin.py`、`E_pol_kick_response.py` (各 `.out.txt`)。
使ったリング: 12 セル FODO (3 GeV・0.5 GeV・4.2 GeV)、大きなリング A (`pysad/tests/kek-ls_v03_12.lat`、3 GeV、603 m、3563 要素)、PF リング (`example/pf/lattice/v9_16_asym.sad`、2.5 GeV、187 m)、
ASTRID2 (`ASTRID2_v01_01.sad`)・ASTRID2 + MPW (`ASTRID2_v01_02_id.sad`、0.58 GeV)。SP8 (大きなリング B、8 GeV) は使っていない。

## 1. 一覧

| 項目 | 結論 | モード | 試験 |
|---|---|---|---|
| (1) 傾いてずれたソレノイドの Emittance[] の POL の量 | **実装した**。SpinTune の追跡との差 1.5e-4〜1.1e-3 → 1e-13〜5e-7、n0 の向き 4e-3〜1.5e-2 rad → 0〜1e-5。大きなリング A・PF・ASTRID2・ASTRID2 + MPW でも同じ。閉じた式 (回転した縦磁場の歳差) とも一致 | v2 だけ (k64 は K64 の式のまま) | `test_v2_E_spinsol.py` |
| (1') 見つけた pysad 自身の不具合 3 件 | 軌道ゼロの四極の Emittance[] (POL) のスピン軌道結合 D が 0 / 領域内の四極 (K1<0・ROTATE) でスピンの座標系を回さない (Emittance[] と track6) / 領域内の CAVI のスピンが追跡で歳差しない。前 2 件は gs (k64) と全桁一致に直った | 前 2 件は両モード、CAVI は v2 | 同上 |
| (2) ANGLE=0 の BEND (補正磁石) と MULT の F1 の Emittance[]・CALC | **実装した**。k64 は gs と全桁一致 (FODO 14 通り、大きなリング A・PF・ASTRID2・ASTRID2 + MPW の挿入) | 両モード (追加) | `test_v2_E_steer.py` |
| (2') B108 | **K64 の誤りと確定** (端の項が約 2 倍)。v2 だけ直した。独立な厳密な幾何 (円弧と直線のロビンソン) と 3〜4% で一致 | v2 だけ | 同上 |
| (3) 脱偏極率の追跡との差 (2.3 倍) | **確定的な蹴りへの応答で再現した**: 全座標 (x・px・py・δ) の蹴りへのスピン軸の応答は、追跡が線形の軌道モデルの 2.07 ± 0.24 倍 (δ の応答 1.9 倍、水平 1.3 倍、垂直 1.0 倍)。振幅 (2σ 以上、シンクロトロン + ベータトロン) とともに増え、低エネルギー (ASTRID2) では 1.0。**式には入れられず**、限界を表にした (§4) | — | `test_v2_E_pol_kick.py` |

## 2. 項目 1: 傾いてずれたソレノイドの POL の量

| 欄 | 内容 |
|---|---|
| 状況 | D1・B1: Emittance[] の POL の量は SOL 領域の中だけ K64 の式。ずれて傾いたソレノイド (GEO=1 の SOL に `DX DY DPX DPY`) で SpinTune が追跡と 2.3e-4 違う |
| 原因 | K64 の `tradke` は、運動量の向きの磁場の積分 (縦磁場 bz·al、端の ψ) を **z 軸の向きの磁場** として足す。運動量が z 軸から傾くと (傾いた軌道・軸を外れた粒子)、(1) 向きの変化から作る B⊥ (すでに運動量に垂直) に垂直成分が二重に入り、(2) δ≠0 で (1+a)B∥ の歳差が 1/(1+δ) 倍にならない。追跡 (`track6._sprot` の `bpar`、C2) は直してあった |
| v2 の動作 | `spin.tradke_spin`: 領域内では B⊥ = 向きの変化 (角そのものにする: `asin`) と、B∥ = [bz·al + (bx px + by py)·al/pz + ψ]/(1+δ) を運動量の向きに分けて足す (勾配も)。`solemit._trk` が各刻みの bz・横磁場 (bx, by)・その勾配・端の ψ (`add_bsi`) を渡す: DRIFT・BEND (一様な双極)・QUAD (端の ψ = k x y)・MULT/CAVI/RF つき MULT (多極の横磁場、`_mult_fld`)・SOL の端 (F1 の放射、B∥ = 0)。あわせて v2 だけ: 領域内の DRIFT・QUAD の放射・スピンの刻みを bz·L/n ≤ 0.02 にする (B1 で追跡に入れたものと同じ。弦の誤差)、向きの変化の大きさを sin θ でなく θ にする (領域の外にも効く: N8、平面リングで SpinTune が 2.6e-5 変わる) |
| 根拠 (独立) | (a) **閉じた式**: 領域は一様な縦磁場 bz ẑ_s なので、運動量は B のまわりを回り (角速度 −B、路程 L/p̂z)、スピンの歳差は Ω = (1+aγ)B − (aγ−a)(B·p̂)p̂。p̂ と一緒に回る系では Ω_rot = aγB − (aγ−a)(B·p̂0)p̂0 が定常で、M = R(−B l)R(−Ω_rot l)。出口の座標系の行列は r_出口·M·r_入口 (trots の回転)。傾き 1〜50 mrad・1〜3 T・0.3〜1 GeV の 5 通りで追跡と 5e-7 以内 (放射の刻み)、Emittance[] の srot は追跡と 1e-8 以内。(b) Thomas-BMT の数値積分 (厳密): DRIFT・BEND・QUAD (K1<0・ROTATE つき)・MULT (K0+K1、ROTATE つき)・CAVI・RF つき MULT を BZ = 2 T・0.5 GeV (bz = 1.2 1/m、歳差 2.4 rad) で、軸上のスピン回転が 3e-5 以内、軌道座標についての微分 (D) が追跡 1e-6 + 1%、Emittance[] 3e-6 + 1% 以内 |
| 数値 (SpinTune = Emittance[] − 追跡の 1 周の行列の回転角、n0 の向きの差) | 追跡は放射の損失を 1e-9 倍にして歳差だけ見る。リング先頭に BZ = 2 T・0.8 m の領域。**傾き 3 mrad・ずれ 2 mm / 傾き 10 mrad・ずれ 4 mm**: 小さな FODO (0.5 GeV): 新 −8e-16 / −4e-16、K64 式 −3.8e-4 / −7.3e-4 (n0 4.9e-3 / 1.4e-2 rad)。ASTRID2: 新 1e-14 以下、K64 式 +4.1e-4 / +9.1e-4 (4.9e-3 / 1.5e-2 rad)。ASTRID2 + MPW は同じ。PF: 新 5e-13、K64 式 −3.9e-4 / −1.1e-3 (2.5e-3 / 8.1e-3 rad)。大きなリング A: 新 4e-13、K64 式 −3.8e-4 / −1.07e-3 (3.6e-3 / 1.2e-2 rad)。軸上は新・K64 式とも一致 (`E_sol_tilt_spin.out.txt`)。新と追跡の一致が 1e-13 になるのは、両方が同じ刻みの式を使うため (連続の値との差は上の閉じた式で 5e-7) |
| 試験 | `test_tilted_solenoid_spin_matrix_equals_closed_form_in_tracking_and_emittance_v2[5]`、`test_ring_spin_tune_and_n0_with_tilted_solenoid_v2[3]`、`test_ring_with_tilted_solenoid_k64_formula_is_off_in_spin_tune_and_n0` (`spin.SOL_BPAR = False` で旧式を再現)、`test_solenoid_region_spin_orbit_coupling_emittance_equals_thomas_bmt_v2[9]`、`test_solenoid_region_quad_spin_frame_rotation_k1neg_and_roll_v2`、`test_one_turn_depolarization_vectors_match_gs_k64[2]` |
| k64 モードの扱い | `spin.tradke_spin` の新しい枝・`solemit._trk` 等は `compat.v2()` と `spin.SOL_BPAR` で分けた。k64 は K64 のまま |
| 使える範囲 | v2 の Emittance[]・EMIT の POL (SpinTune・n0・Polarization・脱偏極ベクトル) の、SOL 領域 (ずれ・傾き・F1 つき。領域内は DRIFT・QUAD・BEND・MULT・CAVI・RF つき MULT・SOL) があるリング。1〜3 T・0.3〜3 GeV |
| 使えない・未確認 | bz ≥ 3 1/m の QUAD・MULT の D は B1 のまま (放射の刻み EPS 由来で 1e-4 台)。領域内の QUAD の端の ψ は追跡と同じ模型 (hard edge の Bs なし。軌道が数 mm の範囲で確認、cm のずれは未確認)。領域の外の B∥ (垂直な軌道の傾きの中の偏向磁石の B_y y'、K64 も追跡も入れていない: 12 セル FODO で SLIM の δ 応答が追跡と 13% 違う原因の候補 — §4.5) は未対応。SP8 (大きなリング B)・8 GeV は使っていない |

### 2.1 見つけて直した pysad 自身の不具合 (gate なし)

| 不具合 | 見つけ方 | 直し | gs との一致 |
|---|---|---|---|
| (a) `emit6d._tquade`: 軌道が正確に 0 の四極は、POL ON でも tradke を通らずスピン軌道結合 D (srot(:,4:9)) が 0 | 平面リングの要素ごとの D を追跡の微分と比べた (四極が `dt` 非 0、`dm` 厳密に 0) | `calpol` のときは通す (K64 の `tquade` は常に `tsolque` → `tradke`) | 斜め四極つきの完全なリング (EMIOUT): D の 18 成分が gs と全桁 (7 桁) 一致 (`_GS_PERFECT`) |
| (b) 領域内の四極 (`K1<0` は −90°、`ROTATE`) の座標の回転でスピンの座標系を回さない (`track6._tquad_sol`、`solemit._tquade`)。K64 の `tsolrote` は `srot = rr·srot` | track6 のスピンを Thomas-BMT と比べて K1<0 で 2.4e-3 ずれた。領域の外の `_tquad_ens` は回していた | 両方に入れた (`_spin_rot_z`、`solemit._spin_rot`) | SOL 領域 (BZ = 1 T) に K1<0 と ROTATE つきの四極: k64 の Emittance[] の srot が gs と全桁一致 (`_GS_SOL`)。回さないと x 列が 1.7 違う |
| (c) 領域内の CAVI (VOLT=0 でも) が追跡 (v2、track6) でスピンを歳差させない (`radok = False`: K64 の tsol は CAVI を krad = .false. で呼ぶ) | リング全体で追跡と Emittance[] の R が 0.12 違った | v2 は CAVI も放射・スピンの刻みに通す (`track6n_ext._compile_sol_v2`)。Thomas-BMT と 1e-9 台 | (k64 は K64 のまま) |

## 3. 項目 2: ANGLE=0 の BEND・MULT の F1、B108

### 3.1 F1 (`tsteee`・`tmulte`)

| 欄 | 内容 |
|---|---|
| 状況 | D1: `emit6d._tsteee` は F1・FB1・FB2 で NotImplementedError、MULT は F1・FRINGE で NotImplementedError (FB1・FB2 は黙って無視)、CALC (`lattice`) は F1 を無視。K64 の `tsteee.f`・`tmulte.f`・`tbfrie.f` を読んで直訳 |
| BEND (tsteee) | fb1 = F1 + FB1、fb2 = F1 + FB2、`mfring = FRINGE`。入口: `mfring > 0 か −1` で `tblfre(dxfr1 = −fb1²/ρb/24, dyfr1 = fb1/ρb²/6, dyfra1)`、f1r = fb1/2、断片 n = −1, 0。出口: `mfring > 0 か −2` で f2r = fb2/2、断片 ndiv+1, ndiv+2、本体のあと `tblfre(dxfr2 = +fb2²/ρb/24)`。本体の長さ alc = al·rbc、断片の長さ・放射は `tbendal` (B1 の `_tbendal`)。`tbedge` は `mfring ≠ −2` (入口)・`≠ −1` (出口)。放射の刻みは最後の断片の長さ alr。CALC は同じ `_tsteee` (放射なし) を使う (`lattice._bend_core_cod_matrix`)。FRINGE = −1・−2 は F1 なしでも片側のエッジを省く (高速な numba の `k_bend` は両側だったので、その場合 Python の経路に回す) |
| MULT (tmulte) | `tsetfringepe`: F1・F2 は L≠0 かつ K1≠0 のときだけ (f1in = F1 + F1K1F、f2in = F2 + F2K1F、f1out = F1 + F1K1B、f2out = F2 + F2K1B)。入口: `tqfrie` (mfring ≠ 2) → `tbfrme(ak0/al1, fb1)` (mfring = 1, 3、それ以外 (≠ 2) は fb なし) → `tqlfre(al1, ak1, f1in, f2in)` (mfring = 1, 3) → `tradke(f1in)` か `tsetr0`。出口は逆順 (`tqlfre(−f1out, f2out)`・`tbfrme(−ak0, fb2)`・`tqfrie`・`tradke(f1out)`)。`tbfrme` の線形フリンジ (`tbfrmle`) を足した。放射の基準点の `tsetr0` は入口のフリンジの全部の後 (旧実装は最初の写像の後で、tqfrie と tbfrme の両方があると K64 と順序が違った。`'radstart'` の印で直した) |
| 根拠 (gs、k64) | 12 セル FODO (3 GeV): ステアリング (L = 0.2 m、K0 = ±8e-4 を交互に) の F1 のみ・F1 + FB1 + FB2 + AE1 + AE2・FRINGE = 1/2・−1/−2・DISFRIN = 1・K0 = ±3e-3・F1 なしの FRINGE = −1/−2、MULT (K1 = 0.7・−0.6) の F1・F1 + F2・FRINGE = 1/2・K0 つき・K0 + FB1 + FB2・K2 つき・F1K1F/F1K1B・DISFRIN = 1・ROTATE + DX・DY: CALC のチューン、Emittance[] のチューン (RF なし・RF あり)・U0・運動量圧縮・同期位置・減衰率・εx・εz が gs と全桁 (表示 10 桁)。F1 の効果は U0 で −4.4 eV (ステアリング)、運動量圧縮で 2e-3 (MULT)。**大きなリング** (補正磁石 2 台 (K0 = ±5e-4、F1 + FB1 + AE1 + AE2) または MULT 2 台 (K1 = ±0.05、K0、F1 + F2 + FB1 + FB2) をリングの先頭に挿入): 大きなリング A・PF・ASTRID2・ASTRID2 + MPW (ステアリングのみ) で gs と全桁 (εx は 4e-9 以内) |
| 試験 | `test_steering_bend_f1_matches_gs_k64[6]`、`test_mult_f1_matches_gs_k64[8]`、`test_steering_and_mult_f1_change_the_result_v2`、`test_steering_calc_tunes_equal_emittance_tunes_without_rf_and_rad` (gs を使わない: 6D の写像のチューン = CALC のチューン) |
| 使えない・未確認 | **TrackParticles の F1 (ANGLE=0 の BEND の `tsteer`、MULT の `tmulti`)** は未対応のまま (track6n・track6 が NotImplementedError)。MULT の AE1・AE2 (端の角: `tmulke` の apsi) は pysad が MULT では読まない (黙って無視) ので、`tbende` から tmulte に回る K1≠0・K0≠0 の ANGLE=0 BEND の AE も無視。`tbfrmle` の行列は数値微分 (K64 の `trans1(2,1)` の式は gs と未照合)。逆向きの要素 (`direlc < 0`) の入れ替えは未対応。**ステアリングに ROTATE と DX/DY を同時に付けると gs と CALC のチューンが 1e-4 違う** (F1 に関係なく既存: K64 はこのとき xy を回す行列を返す癖があり、pysad も似た値を返すが一致しない。未調査) |

### 3.2 B108 (ANGLE=0 の BEND の端の角)

| 欄 | 内容 |
|---|---|
| 状況 | C: `tsteee.f` は端の角 apsi1・apsi2 (AE1・AE2) を `tbedge` で蹴りに入れるが `tradke` に行路の補正 (B75 の tin・tout) を渡さない。独立な値が無かった |
| 独立な値 | 硬い端の一様磁場の極 (円弧と直線の厳密な幾何。面は参照軌道に垂直な面から AE だけ傾く、K0<0 は鏡映) と、ロビンソンの定理 Je = 2 + dln(Σ B_i² l_i(δ))/dδ、Jx = 3 − Je (U ∝ E²∮B²dl。l_i は運動量 δ の閉軌道が磁場の中を通る路程)。24 セルの FODO (3 GeV) に、極 3 個の組 (K0・−2K0・K0。軌道が極の外で 0 に戻る) を各セルに入れた。pysad の光学・放射のコードを使わない |
| 結果 (AE なしとの差 ΔJx) | K0 = 0.005 (軌道の傾き 5 mrad) で AE = (0.3, 0.3): 厳密 +0.01440、v2 +0.01487 (1.032)、**k64 +0.02992 (2.08 倍)**。(−0.2, 0.4): 厳密 +0.00574、v2 +0.00596 (1.037)、k64 +0.01238 (2.16)。K0 = 0.0025 でも同じ比 (v2 1.03〜1.04、k64 2.06〜2.10)。くさび (0.3, −0.3) は 1 次が打ち消し、v2 −0.00068 対 厳密 −0.00061。**K64 は端の項が約 2 倍 = B75 と同じ誤り** (`E_b108_poles.out.txt`) |
| v2 の動作 | `_tsteee` (Python と numba の `elkern._tsteee`) の最初の `tradke` に `tin = tan(AE1)`、最後に `tout = tan(AE2)` (ANGLE≠0 の `tbende` と同じ)。k64 は K64 のまま |
| 試験 | `test_steering_edge_angle_term_v2_matches_exact_geometry_and_k64_doubles_it[2]`、`test_steering_edge_angle_wedge_cancels_to_first_order` |
| 使える範囲 | ANGLE=0 の BEND (K1 = 0) の AE1・AE2 つきの Emittance[]・EMIT・追跡の放射以外。K0 が小さく軌道の傾きが 5 mrad 程度の極 (残る差 3〜4% は EMITDIV = 20 の有限の分割と 2 次の項) |
| 使えない・未確認 | K1≠0・K0≠0 の ANGLE=0 BEND と MULT の端の角 (上: AE を読まない)。軌道の傾きが 30 mrad の極 (K0 = 0.03) では独立な値と v2 が 0.02〜0.03 違う (2 次の項。傾き 5 mrad の表は 3〜4%)。追跡 (`tsteer`) の端の放射は未確認 |

## 4. 項目 3: 脱偏極率と FLUC 追跡の差 (4.2 GeV・12 MV で追跡が式の 2.3 倍) の切り分け

### 4.1 方法 (FLUC なし・放射の損失なしの追跡)
B3 の差は FLUC 追跡の統計 (480 粒子) でしか測れなかった。E は**確定的な蹴りへの応答**で分けた: 放射の光子 1 個 = 座標への小さな蹴り h。軌道 u の粒子 A と、u に蹴り ±h を足した B± を同じスピン (n0) から K = 1500 周追い (放射の損失は `emit6d.CAVE` を 1e-9 倍にして消す。RF の同期粒子は損失が無くても加速しないことを確かめた)、
ζ = s·(m + i l) (m, l ⟂ n0) の、**スピンチューンの線の複素振幅** (NAFF、k=0 に外挿) の差 (a_B+ − a_B−)/(2h) = G = ∂n/∂(座標) を得る (B± の円錐の大きさ = Δn。軌道の線 νx・νy・νs の組み合わせとは分離する)。脱偏極率は Σ_i ½ Σ_jk B_i,jk Re(G_j G_k*) (B_i は要素の拡散行列)。
同じ初期条件の組 (ベータトロン振幅だけ / シンクロトロン振幅だけ / 両方。振幅は平衡のビーム行列から、モードごとの不変量が指数分布になるよう抽出し f 倍) を、追跡 (M2: 非線形な軌道 + 厳密なスピン) と、線形の軌道 + 要素ごとの線形のスピン軌道結合の厳密な回転 (M1 = B3 の MC と同じモデル) の両方に通し、比 M2/M1 を見る (`pysad/tests/_e_spinkick.py`、`example/v2_next/E_pol_kick_response.py`)。
方法の検算: M1 の閉軌道上の δ の蹴りへの応答は SLIM の N0 の δ 列 (n に垂直な成分) と 0.02% で一致 (2.7283 対 2.7278)。

### 4.2 結果

| 条件 (12 セル FODO、四極の縦ずれ ±3 mm) | 蹴り | 比 <\|G2\|²>/<\|G1\|²> (M2 追跡 / M1 線形の軌道) |
|---|---|---|
| 4.2 GeV・12 MV、δ だけ。ベータトロンのみ f = 0.5 / 1 / 2 (各 10 試料) | δ | 0.75 ± 0.04 / 0.92 ± 0.06 / 1.27 ± 0.12 |
| 同、シンクロトロンのみ f = 0.5 / 1 / 2 | δ | 0.95 ± 0.09 / 0.96 ± 0.08 / 1.70 ± 0.25 |
| 同、両方 f = 0.5 / 1 / 2 | δ | 0.95 ± 0.08 / 1.08 ± 0.13 / **2.10 ± 0.26** |
| 4.2 GeV・6 MV (f = 1、12 試料): ベータトロンのみ / シンクロトロンのみ / 両方 | δ | 0.87 / 1.17 / 1.08 |
| 4.2 GeV・24 MV (f = 1、12 試料): ベータトロンのみ / シンクロトロンのみ / 両方 | δ | 1.24 / 2.06 ± 1.0 / 1.29 |
| ASTRID2 (0.58 GeV、60 台でなく 12 台の四極の縦ずれ rms 1 mm)、ベータトロン f = 1 (8 試料) / 両方 (3 試料) | δ | 1.02 / 1.01 |
| **4.2 GeV・12 MV、両方 f = 1、偏向磁石の直後で 4 座標 (x, px, py, δ) の蹴り、同じ 12 試料の組** | x, px, py, δ + 拡散行列 | **2.07 ± 0.24** (δ の応答 1.9 倍、x・px 1.3 倍、py 1.0 倍) |

- **ゼロ振幅 (閉軌道) の応答も追跡は線形モデル (Emittance[] の N) より 13% 小さい** (G = 2.36 対 2.73。追跡の刻みを 12 倍にしても同じ = 追跡の刻みの粗さではない)。一周の D の x・px 列の n に垂直な成分が追跡と 10〜20% 違う (`kr_test4`)。原因の候補は、領域の外の偏向磁石の中の縦磁場成分 B∥ = B_y y' (垂直な軌道の傾きがあるとき。K64 も追跡も入れていない)。**要判断**。
- 線形のモデルでは偏向磁石の放射の脱偏極への寄与は δ だけでなく **px の蹴り (拡散行列の px 成分は δ の 2%、N の px 列は δ 列の 13 倍) が支配的** (δ だけで見ると線形の寄与を 3 倍過小に見る。`rate_prof`)。δ の蹴りだけの比 (上の上 3 段) は、全座標の比 (最下段) の指標として不十分。
- 振幅に対する応答の増え方が急で (M1 でも <G²> が f = 1 → 2 で 4.5 倍)、脱偏極率は 2σ 以上の裾の寄与が大きい。追跡 / 線形の比は振幅とともに増える (δ のみ両方: 0.95 → 1.08 → 2.1)。
- 側帯の和 (ν_spin + kν_s、k = −4..4) の応答は中心の線の 12% 増しにとどまる (M1)。

### 4.3 切り分けの結論
1. **B3 の 2 倍の差は、追跡のスピン軸の蹴りへの応答 (全座標) が線形の軌道モデルの 2.07 ± 0.24 倍であることで再現できた** (FLUC 追跡の 2.1 倍と一致)。確定的な追跡で見えるので、乱数の統計や FLUC の実装の問題ではない。
2. 原因の座標は **δ の応答 (1.9 倍) と水平 (x・px) の応答 (1.3 倍)**。垂直 (py) は同じ。ベータトロン振幅だけでは δ の応答の比は 1σ で 0.9、2σ で 1.3、シンクロトロン振幅だけ・両方では 2σ で 1.7〜2.1 に増える。**振幅 (特に 2σ 以上) の増加とともに軌道の非線形が効く** = 振幅依存のチューン・高次の共鳴・ベータトロン側帯・分散つきのスピン軌道結合のうち、**シンクロトロン振幅 (νs の変調) とベータトロン振幅の組み合わせ**が支配的で、ベータトロンだけでは効かない (1σ)。どの非線形かまでは分けていない (小さなリングに6 極が無いので、色収差・運動学的な項・RF の sin の非線形が候補)。
3. 低エネルギー (ASTRID2、0.58 GeV) では比は 1.0 (非線形の効果なし)。B3 の 3 GeV・2.9 GeV の追跡 / 式の比 (0.54〜0.87) とも合う。**高エネルギー (aγ σδ/νs ≳ 0.5)・高 RF 電圧のリングで効く**。
4. **式には入れられなかった**: 応答の比は振幅の関数で、平衡分布で平均すると 1σ の少数試料の値 (δ のみの回帰: 平衡分布の平均 0.8〜1.0、全座標の比: 2.1 ± 0.24) が試料数・蹴りの座標の取り方で変わる (誤差が大きい)。解析式に入れるには 2σ 以上の裾の応答を多数試料 (数百) で測る必要がある。

### 4.4 表 (限界の記録): v2 の Polarization の脱偏極率 / FLUC 追跡
| リング・エネルギー・RF | v2 の式 | 追跡 | 比 (式/追跡) | 出典 |
|---|---|---|---|---|
| 12 セル FODO 4.2 GeV・12 MV | 側帯込み 4.46e-8 | (10.4 ± 1.6)e-8 | 0.43 | B3 §2.3。E の応答の比 2.07 ± 0.24 と整合 |
| 同 4.2 GeV・6 MV / 24 MV | 3.50e-8 / 4.93e-8 | 未取得 (線形軌道 MC 3.83e-8 / 5.50e-8) | (応答の比は 1.08 / 1.29。全座標は未測定) | B3・E |
| 12 セル FODO 3 GeV・1.5 MV / 2.903 GeV・1.5 MV | 6.26e-11 / 1.20e-9 | 3.4e-11 / 1.05e-9 (±22〜24%) | 1.8 / 1.14 | B3 §2.4 |
| ASTRID2 (0.58 GeV) | 1.43e-12 | (1.7 ± 0.3)e-12 | 0.84 | B3 §2.2。E の応答の比は 1.0 |
| 大きなリング A (3 GeV)・PF (2.5 GeV) | 8.8e-14・3.5e-9 | 未取得 | — | 追跡 1 粒子 1 周が 0.4 s (A) と遅く、応答の測定 (3〜9 粒子 × 数百周 × 試料) は PF で 1 試料 1 時間以上。E は ASTRID2 だけ取った |

### 4.5 要判断・残した項目 (項目 3)
- 全座標の応答の比 2.07 ± 0.24 は 12 試料 (位置 1 か所、4.2 GeV・12 MV のみ)。6 MV・24 MV・PF の全座標の応答、位置依存、試料数 100 以上は未実施。
- ゼロ振幅の線形応答の 13% の差 (§4.2): 領域の外の偏向磁石の B∥ = B_y y' を含めると変わるか未確認 (K64 も追跡も入れていない)。
- 追跡の δ の応答が 2σ 以上で急に増える原因 (共鳴か非線形のどれか) は未特定。

## 5. 全試験・回帰 (2026-11-18)
全試験 (`--ignore` の 3 本を除く、`PYSAD_SAD_PACKAGES` つき、4 組に分けて並行): **1342 通過・5 skip・失敗 0** (E の追加 41 件 + `PYSAD_SLOW` の 1 件 skip を含む)。回帰 (`QUICK=1`): **50 ベンチ FAIL 0** (159 s)。
既存の試験の変更 1 件: `test_v2_B1_optics.py::test_emittance_f1_unsupported_is_still_explicit_for_steering_bend` → `test_emittance_f1_of_steering_bend_is_implemented_by_e` (ステアリングの F1 を実装したので)。

## 6. 親が更新する文書
- docs/21: B108 を「確定」(K64 が端の項を 2 倍にする。AE1・AE2 つきの極、v2 で行路を x tan(e) だけ短くする)。新しい項目 (番号は親): (i) `tquade` の軌道ゼロで POL のスピン軌道結合が抜ける (pysad 自身、両モード修正)、(ii) 領域内の四極の srot の回転が抜ける (同)、(iii) K64 の `tradke` は縦磁場の積分を z 軸の向きに足す (傾いた SOL で SpinTune が 1.5e-4〜1.1e-3、n0 が 2.5e-3〜1.5e-2 rad ずれる。v2 で修正)、(iv) 領域内の CAVI の spin が歳差しない (K64 の tsol は krad = .false.。v2 は物理的に歳差させる)。
- docs/44 §7a・docs/46・docs/47: 「Emittance[] の POL の量は領域の中だけ K64 の式」を「v2 は修正済み」に。「ANGLE=0 の BEND・MULT の F1 は未対応」を「Emittance[]・CALC は実装 (gs と一致)、追跡 (TrackParticles) は未対応」に。B108 の「要判断」を解消。docs/07 の B108 の記述も同様。
- 数値: 傾いたソレノイドの SpinTune の差 (§2)、F1 の gs との一致 (§3)、B108 の比 (厳密 / v2 / k64 = 1 / 1.03 / 2.08)。
