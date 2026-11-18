# v2next 担当 A1 の記録 (2026-11-17): すぐ直せるものと gs との未解決の差 8 項目

試験: `pysad/tests/test_v2_A1_next.py` (40 件)。例: `example/v2_next/A1_*.sad` (pysad の v2 の出力 `.pysad.txt`、k64 の出力 `.pysad_k64.txt`、gs の出力 `.gs.txt`)。
gs の確認は `local/wake/gsrun.sh`。互換モードの規則は `local/wake/v2_brief.md`。ここに書いた gs の値は 2026-11-17 の gs (10-19 ビルド) の出力。

確認に使ったリング (小): 3 GeV の FODO リング (周長 41.8〜105.6 m、RF つき)。(大): 「大きなリング A」(`pysad/tests/kek-ls_v03_12.lat`、3 GeV、周長 603 m、3563 要素、RF つき)、PF リング (`example/pf/lattice/v9_16_asym.sad`、2.5 GeV、187 m)、ASTRID2 (`example/astrid/lattice/ASTRID2_v01_02_id.sad` の `RINGMPW`、0.58 GeV、MPW 入り)。大きなリング B (8 GeV) は使っていない (下の「残した項目」)。

## 一覧

| 項目 | 結果 | モード |
|---|---|---|
| 1 N12: INS の `Twiss["DETR"]` | 直した (pysad の不具合)。gs と全桁一致 | 両モード (gate なし) |
| 2 TrackParticles の MASS・CHARGE | 直した (pysad が渡していなかった)。gs と全桁一致。**K64 の放射の古典半径 (電子の値) の誤りを v2 で直した** | MASS・CHARGE の受け渡しは両モード。古典半径は v2 だけ |
| 3 `emit6d._tquade` の解析ヤコビアン | **記述の症状は再現しなかった** (Python 版・numba 版・数値微分が 1e-12 で一致)。関数単体の z の行 (列 x px y py) の欠落を v2 で補った | v2 だけ (結果は変わらない) |
| 4 BEND の F1 + FRINGE + RAD の追跡 | 直した (K64 の tfgeo の inext・iprev の再現)。gs と一致 | 両モード (gate なし。K64 の動作) |
| 5 L≠0 の CAVI + DX・DY | 直した (CAVI の TENT・TEXIT が無かった。V1 V20 V11 V02 も入れた)。gs と全桁一致 | 両モード (gate なし。K64 の動作) |
| 6 EMIT のあとの USE | 直した (trf0 が USE で 0 に戻っていた)。gs と全桁一致 | 両モード (gate なし。K64 の仕様) |
| 7 RADTAPER の追跡の閉軌道 | 原因を特定して v2 で閉じた (K64 の tradke の符号)。3 GeV の小リングで 4e-7 → 8e-10 m | v2 だけ。k64 は gs と同じ |
| 8 PSPAC の x と y が gs と 16% 違う | 16% の再現はできず、原因を **メッシュの分解能** と特定。警告を足した | v2 だけ (PSPAC 自体が v2 だけ) |

## 項目 1: N12 `Twiss["DETR"]` が CALC の後で 0

| 欄 | 内容 |
|---|---|
| 状況 | `ffscmd._rows_from_tws` が `TwissRow` を作るときに `detr` を渡していなかった。INS の CALC (`_ins_rows`) と InitialOrbits・wake の行が通る。CELL の CALC (リング) は別の経路 (`calc`、`calc6d`) で `detr` が入る |
| K64 の動作 | `twiss(i, mfitdetr)` = 結合行列の行列式 (DETR < 1 は R1R4 − R2R3、DETR ≥ 1 は 1 + XYTH − R1R4 + R2R3 の形) |
| v2 の動作 | k64 と同じ (gate なし)。`_rows_from_tws` に `detr=tw.get('detr', 0.0)` を足した |
| 根拠 | gs: 45 度回した四極 (K1=0.8) を入れた結合リングで、CELL CALC の後の INS CALC の M2 での DETR = 0.0021332033027、R1 = −0.0800626854892 (pysad 0.0021332033028)。DETR = R1R4 − R2R3 の恒等式 (1e-12) |
| 大きなリング | 大きなリング A (3 GeV、603 m、QDA2 を 0.03 rad 回し、空洞に DX=2e-4・DY=−3e-4): CELL CALC の後と INS CALC の後の `Twiss["DETR", {1, 500, 1000, 2000}]`、R1 とも gs と全桁 (10 桁) 一致 (DETR(1) = −0.0293428076、他 0.0123505545)。PF リング (Q2 を回す): 全桁一致 (DETR 0.0256283385、0.1315129060) |
| MEASURE への効果 | `measure.mark_values` が読むのは先頭 MARK の生の値 (CELL CALC の `_twmov_to_mark` が書く)。INS の CALC の後の `Twiss["DETR"]` ではないので、MEASURE の値は元から正しかった。結合リング (K1=0.8) の M2 の MEASURE の σx = 4.3381e-5 m は `Twiss["SIGX"]` = 4.3328e-5 m と 0.12% (NP = 20000 の統計) で一致。`Twiss["SIGX"]` など INS の行から作る量は DETR を使うので、DETR ≥ 1 (モードが入れ替わる) の INS では旧 pysad の値が違った可能性がある |
| 試験 | `test_n12_ins_calc_keeps_detr[k64、v2]` |
| 使える範囲 | 結合したリングの INS の CALC の DETR と R の関係、gs との一致 (DETR < 1) |
| 使えない・未確認 | **DETR ≥ 1 (モードが入れ替わる) の INS と MEASURE は未確認**。DETR < 1 のリング (R1 = 0.3 まで、K1 = 0.5〜5 の QS を調べた) は 0.24〜0.43 で入れ替わらなかった。入口を DETR > 1 の生の値で与える輸送路は gs と pysad で入口の規格化 (R が 1e-17、DETR 1.375 と、R がそのまま DETR = 1) が違い、比較できなかった |

## 項目 2: TrackParticles の MASS・CHARGE

| 欄 | 内容 |
|---|---|
| 状況 | `fn_track` が `TrackRun` に `amass`・`charge` を渡さず (電子の質量・電荷 1 のまま)、`Emittance[]` (`emittance6d`、`emittance_rules`)、CALC6D (`calc6d`、`calc6d_ins`) の呼び出しにも渡していなかった。`ffscmd._mass_charge()` (電子と同じ値のときは何も渡さない) を足し、`TrackRun`・track6 の経路・`emittance6d`・`emittance_rules`・`calc6d` に渡した |
| K64 の動作 | `tfsetparam.f`: `charge = CHARGE`、`amass = MASS`、`p0 = pgev/amass`、`brhoz = pgev/c`、RF は `v = vc/amass |charge|`、`phic = (phi+dphi)·charge`、ソレノイドは `bzs = charge·BZ/(p0 c/amass…)`。`crad = sign(rclassic (c/amass)²/p0/1.5, charge)`、`rclassic = charge²·re0`、**`re0 = elradi` (電子の古典半径。質量によらない)**。BEND の ANGLE・四極の K1 は幾何なので質量・電荷に依らない (磁気剛性は p/(charge c)) |
| v2 の動作 | MASS・CHARGE の受け渡しは両モード同じ (K64 の動作)。**放射**: v2 は `re0 = elradi·elmass/amass` (粒子の古典半径)。k64 は K64 のまま (電子の値) |
| 根拠 (運動学・磁気剛性) | ドリフト 10 m、陽子 3 GeV/c: z = L (β1/β0 − 1) (δ = +1e-3: 8.89799109e-4 m、δ = −1e-2: −9.03335291e-3 m) で 1e-12 以内。z は δ と正準に共役な β c (t0 − t) 型の座標 (K64 の dv の式の積分)。px だけの粒子は z = −L px²/2。電子の質量のままだと 1e-3 の δ で 3e-10 m しか出ない。ソレノイド (BOUND=1 の SOL 領域、L=2 m、BZ=1 T、3 GeV/c): 厚いソレノイドの行列 x_out = C² x0、y_out = −S C x0 (k = charge BZ c/(2p)、C = cos kL) と 5e-12 m 以内、電荷 −1 で回転の向きが逆、電荷 2 で k が 2 倍 |
| 根拠 (Qs) | 30 GeV/c の陽子 (γ = 32、周長 45.6 m、h = 24、V = 2 MV): 1 周の z-δ 写像 (TrackParticles の 1 周の中心差分): δ の傾き M[δ,z] = 2π h V/(C β² E) と 1e-3 以内、η は RF なしの 1 周の z の変化 (η = 1.44289482/C)、Qs = sqrt(h η V/(2π β² E)) と 1e-4 で一致、行列式 1 (1e-9)。Emittance[] の Qs (−0.00283934) は gs と全桁一致 |
| 根拠 (放射) | 陽子の U0: 古典式 U0 = (2/3) r_p mc² β³γ⁴ Σ θ²/L = 6.8800e-4 eV (r_p = r_e m_e/m_p)。v2 の追跡の 1 周の δ の変化 (RF なし): −U0/(β³ p c) と 1e-6 で一致 (質量 1×・2× 陽子で確認)。K64 の tradkn は光子数 (5/(2√3)) α γ θ と臨界エネルギーから δ を Δδ = −(エネルギー)/(p c) だけ下げる超相対論的な式で β の因子が無く、U0_K64 = U0/β³ (γ = 32 で +0.15%)。Emittance[] の EnergyLossU0 も同じ (6.88996e-4 = U0/β³)。k64 は 1.2651 eV (mp/me = 1836 倍。gs と全桁一致) |
| gs との比較 | 小リング (30 GeV/c の陽子): Emittance[] (Tunes、U0、EquilibriumPosition、BucketHeight、HarmonicNumber、OrbitAtExit など) と、1 周・10 周・100 周の追跡 (NORAD・RAD) が k64 モードで gs と全桁一致 (表示 6〜9 桁)。大きなリング A (陽子、3 GeV/c): Emittance[] (Tunes、U0 = 5.7641486711e-5 eV、EquilibriumPosition、MomentumCompaction、BucketHeight) が 1e-9 以内、50 周の追跡 (NORAD・RAD) が 1e-10 以内で一致 |
| 試験 | `test_track_proton_drift_z_exact[k64、v2]`、`test_track_solenoid_larmor_charge[電荷 1・−1・2 × 両モード]`、`test_proton_rf_kick_and_synchrotron_tune[両]`、`test_proton_radiation_loss_v2`、`test_proton_radiation_loss_k64`、`test_electron_radiation_same_in_both_modes` (電子の放射は両モードでビット単位で同じ) |
| 使える範囲 | `MASS`・`CHARGE` を FFS の中で与えた陽子・重い粒子・多価イオンの TrackParticles・Emittance[]・CALC6D・RF の位相・ソレノイド。放射の U0 は γ ≫ 1 で正しい (v2) |
| 使えない・未確認 | (1) `CHARGE` は FFS に入るたびに 1 に戻る (gs と同じ。FFS の中で与える)。(2) 放射の β の因子 (K64 の超相対論的な式): β < 0.99 の粒子の放射は U0 より 1/β² ほど大きく出る。放射が意味を持つのは γ ≫ 1 なので影響は小さいが、v2 は直していない。(3) **TouschekLifetime・Intrabeam (`touschek.py` の `rclassic = charge²·elradi`) と RADINT は K64 の電子の古典半径のまま** (v2 で直していない、陽子での確認なし)。(4) MASS・CHARGE を MAIN で与える経路は `maintrack` が既に渡していた |

## 項目 3: `emit6d._tquade` の解析ヤコビアン

| 欄 | 内容 |
|---|---|
| 状況 | docs/46 §4 の記述 (x のずれた軌道で δ の微分を落とし、`emit_extra`・`calc6d`・`touschek` で T の δ 列が 1e-5 ずれる) を調べた。(a) 四極 1 個 (`_jac_quad_full`、軌道 x, y, δ ≈ 1e-3、DX・DY・ROTATE あり) の解析ヤコビアンの δ の列は数値微分と 1e-12 以内。(b) 要素単位の T (`_tquade`、`State.piece` が z の行を組み直す) は Python 版・numba 版 (`_quade_fast`)・写像の数値微分で 1e-12〜1e-9 以内で一致 (K1 = ±0.6〜±3、L = 0.3〜1.5、DX・DY・ROTATE、軌道 5 mm までの 24 例)。(c) 大きなリング A (DX・DY をランダムに 3e-5 m 入れた 440 台の四極、軌道 0.9 mm) の 1 周の T: Python 版 (`emit_extra._turn_record`)・numba 版 (`tturne`) が 1e-13 以内で一致 (旧記述の 1e-5 は出ない)。(d) 全 3563 要素を 1 個ずつ数値微分と比べると、差が出るのは BEND だけ (T[z, δ] で 2.3e-7/個。他の要素は 5e-12 以下) |
| 結論 | **旧記述の差は再現しない**。記述の後に直っていた可能性が高い (原因の `_jac_tsolque_body` の δ 依存の抜けは確認できず)。直す対象は無かった |
| 見つけたこと | (i) `_jac_tsolque_body` が返す 6×6 の z の行 (列 x px y py) は 0 のままで、`State.piece` がシンプレクティック条件で組み直していた (K64 の `trans(5,j)` と同じ) ので結果に影響しない。v2 では関数単体として完全にした (`J[4,0..3] = −Z_{x,p,y,q}/4`、数値微分と 1e-10)。k64 は今まで通り 0 の行 (結果は同じ)。(ii) 解析 T と 1 周の写像の数値微分の差 (大きなリング A で 2e-5) は BEND の T[z, δ] の K64 の近似 (`tbendebody0` の (5,6) 成分) で、gs と同じ (K64 の仕様) |
| 試験 | `test_quad_jacobian_full_v2`、`test_quad_jacobian_k64_keeps_zero_z_row`、`test_tquade_T_matches_numeric_with_orbit[k64、v2]` (24 例) |
| 使えない・未確認 | `RADCOD` ON の 1 周の T の解析微分と写像の数値微分の差は 2e-3 (T[x, px]、小リング。k64 でも v2 でも同じで、K64 の `tradke` の微分が近似) |

## 項目 4: BEND の F1・FRINGE + RAD の追跡

| 欄 | 内容 |
|---|---|
| 状況 | 最小の例は、`F1=0.1 FRINGE=3` の BEND 1 個の輸送路 (INS) と、隣り合う BEND の対のリング。NORAD は全桁一致、RAD で差が出る。ずれの量は F1 に比例し (F1 = 0.01 で 0.3%、0.1 で 2.8%)、FB1 だけ・FB2 だけならその半分 |
| K64 の動作 | `tbend.f`・`tbendi.f`: `iprev(l_track) == 0 .and. fb1 /= 0` のとき入口側に F1 の断片 (n = −1、0。長さ f1r = F1/2 の 2 つに放射 `tradk(…, alr = f1r)`)、`inext == 0 .and. fb2 /= 0` のとき出口側に断片 (n = ndiv+1、ndiv+2) を入れる。`inext`・`iprev` は `tfgeo.f` が作り、「同じ型 (BEND) の次の要素。間は長さ 0 の要素だけ」を繋ぐ。リングの終わりは先頭に回って探すので、**先頭が長さ 0 の MARK の輸送路の最後の BEND は自分自身に繋がり (iprev = inext = 自分)、F1 の放射の断片が無い** |
| 原因 | pysad は `first`・`last` (iprev・inext が 0 か) を常に True にしていた (`track6n` の `k_bend0`・`k_bendi` への引数、`track6` の `_tbend_ens` の引数が常に True)。隣り合う BEND (B と B2、または B MARK B2) と、1 個しかない BEND の輸送路で F1 の放射の断片を余計に入れていた。DX・DY・ROTATE だけなら断片が無いので 1e-13 で一致していた |
| v2 の動作 | k64 と同じ (gate なし。K64 の動作)。`track6n.bend_links` (tfgeo の inext・iprev の再現)、`Prog`・`TrackRun`・`track6.track_particles`・`ffscmd.fn_track` に渡す |
| 根拠 | gs: 1 個の BEND (INS) の 1 周 (RAD): δ = −2.17142054679e-5、9.782447092820e-4、−1.0021281980e-2 (旧 pysad −2.11131e-5、3%)。隣り合う BEND の 12 セルのリング: 3.6% → 1e-6 (8 桁表示の丸め)。K1≠0 の BEND の対 (`tbendi`)・POL (Python の追跡 `track6`) も 1e-6 で一致 |
| 大きなリング | 大きなリング A・PF リング・ASTRID2 (MPW 入り) の BEND は F1 なし (または隣り合う対なし) なので結果は変わらず、gs と一致 (RAD の 5 周・10 周で 4e-10 以内) |
| 試験 | `test_bend_f1_single_bend_line_matches_gs[k64、v2]`、`test_bend_f1_adjacent_bends_match_gs[k64、v2]`、`test_bend_links_function` |
| **見つけた別の差 (直していない)** | **`Emittance[]`・EMIT・CALC6D は BEND の F1・FRINGE を無視する** (`emit6d._tbende` は F1 = 0 を前提。pysad が「未対応」と言わず黙って F1 なしの結果を返す)。F1=0.1 FRINGE=3 の 12 セルのリングで gs は U0 739674 eV、νy = −0.13957 (Emittance[] の Tunes) だが pysad は 781633 eV、−0.133128 (F1 なしの値と同じ)。CALC (Twiss) は F1 を入れていて gs と 1e-9 で一致するので、CALC と Emittance[] のチューンが食い違う。K64 の `tbende` (tblfre・tbendef1・f1r・f2r・bradprev) の直訳が必要。v2 は黙って違う値を返さないよう、Emittance[]・EMIT・CALC6D で F1 のある BEND を見つけると 1 度だけ警告する (`???-Emittance[]-Warning-BEND …`。試験 `test_emittance_warns_about_ignored_bend_f1`)。**要判断・残した項目** (下) |

## 項目 5: L≠0 の CAVI + DX・DY

| 欄 | 内容 |
|---|---|
| 状況 | 例: L = 0.5 m の CAVI に DX = 1e-3、DY = −2e-3 のリングの RAD の追跡。gs では x = y = 0 から始めた粒子に y = 6.128e-9、py = −1.6e-13 が入る (pysad は 0)。x 軸の差は 1e-3 |
| 原因 | K64 の `tcav` は入口で `include 'inc/TENT.inc'` (ずらして θ 回す)、出口で `TEXIT.inc` を実行するが、pysad の追跡の CAVI (`track6n.k_cav`、`track6._tcav`) には DX・DY・ROTATE が無かった (K64 の `tcave` (Emittance[]) は `tchge` で入れていた)。さらに、**パーサが CAVI の V1・V20・V11・V02 を捨てていた** (`parser.py` は PHI・DPHI・DVOLT・AUTOPHI だけ属性にした)。`k_cav` は V1・V20・V11・V02 の式 (tcav.f の v1a・v20a・v11a・v02a) を持たなかった (TRPT の `cav_full` だけ持つ) |
| v2 の動作 | k64 と同じ (gate なし。K64 の動作)。`track6n.tent_all`・`texit_all` (TENT・TEXIT の書き方どおり)、`k_cav_q` (V1・V20・V11・V02 つき。全部 0 のときは今まで通り `k_cav` なので既存の結果はビット単位で同じ)、`parser` が V1 V20 V11 V02 を CAVI の属性にする、`track6._tcav` にも |
| 根拠 | gs と全桁一致 (表示 9 桁、NORAD・RAD、POL の Python の追跡も): DX・DY のみ、DX・DY・ROTATE=0.3、L=0 の CAVI、V1=1e5 V11=0.5e6 V20=1e6 V02=−3e5 と DX・DY・ROTATE=0.2 (L=0.5) |
| 大きなリング | 大きなリング A (603 m): 空洞に DX=2e-4・DY=−3e-4、QDA2 を 0.03 rad 回して RAD の 5 周追跡: gs と差 0 (表示 11 桁)。EMIT のあと 5 周: 1e-11。PF リング: 差 0 (DX=2e-4、DY=−3e-4) |
| 試験 | `test_cavity_offset_and_transverse_voltage_match_gs[2 つ × k64、v2]`、`test_cavity_offset_changes_orbit` |
| 使えない・未確認 | WakeFunction のある CAVI の V1 V20 V11 V02 (`NotImplementedError`)、SOL 領域内の CAVI の V1 V20 V11 V02 (既存の `NotImplementedError`)、Emittance[] の CAVI の V1 V20 V11 V02 (K64 の `tcave` は V=0 を前提とした pysad の直訳、docs/history の記録どおり) |

## 項目 6: EMIT のあとの USE

| 欄 | 内容 |
|---|---|
| 状況 | gs で確かめた: `EMIT` (または `Emittance[]`、`CALC6D`) のあとに `USE RING;` しても、追跡の RF の位相は EMIT の値のまま。`CELL; CALC;` のあとも変わらない。pysad は `USE` で `self.trf0 = 0.0` にしていた |
| K64 の動作 | `trf0` は `tmacro` の大域変数。`temit.f` が 0 にしてから `tturne` で計算し、`tffsa.f` の MEASURE (`trf0 = 0`) だけが他に触る。`USE` (tffsa.f の `chguse`) は触らない。ほかの EMIT の跡 (`emit` の結果、Twiss) は USE で作り直される |
| v2 の動作 | k64 と同じ (gate なし。K64 の仕様)。`USE` で `trf0` を持ち越す (最初の値は 0) |
| 根拠 | gs: EMIT の前 δ = (−3.0771610e-4、6.82902209e-4、−1.023444e-2)、EMIT 後 (−4.1421124e-5、9.50018237e-4、−9.97733e-3)、USE 後 (同じ)、さらに CELL CALC 後 (同じ)。pysad (旧) は USE 後に EMIT 前の値に戻った (相対 0.87)。新: EMIT・Emittance[]・CALC6D のどれでも gs と全桁一致 |
| 大きなリング | 大きなリング A: CELL CALC・RAD・EMIT の 5 周追跡 → USE RING → 5 周追跡: gs と 1.3e-11 (表示の丸め)。PF リング: 差 0 |
| 試験 | `test_use_after_emit_keeps_trf0[EMIT、Emittance[] × k64、v2]` |
| 使えない・未確認 | `emit`・Touschek の結果 (`self.emit`、`self.touschek`) は USE で捨てる (K64 では emit の結果が何処に残るかは確かめていない) |

## 項目 7: RADTAPER の追跡の閉軌道

| 欄 | 内容 |
|---|---|
| 状況 | 再現: 3 GeV、周長 105.6 m の FODO 24 セル、RF 4 MV: Emittance[] (RFSW・RAD・RADCOD) の閉軌道 (OrbitAtExit、z = 0 から) から TrackParticles (RAD、NOFLUC) を 1 周すると、x が −4.349e-7 m、z が 2.375e-6 m 動く (gs も同じ: −4.3492e-7、2.3753e-6)。RADTAPER あり: x 1.967e-6 m、z −1.077e-5 m (gs・pysad k64 とも)。要素ごとに `tturne([el])` と追跡 1 要素を比べると、差が出るのは BEND だけ (x 1.5e-7、px 2.6e-7 m/個、RAD ON だけ。K0 に依らず) |
| 原因 | K64 の `tradke` (Emittance[]・閉軌道、`photons.f:1417`) と `tradkn` (追跡、`photons.f:1247`) で、放射の反作用の角度の項の符号が逆: tradke は `ddpx = +0.5 dpx dg`、tradkn は `ddpx = −0.5 dpx dg` (dg < 0)。エネルギーを失う間は同じ磁場が運動量の小さい粒子をより強く曲げる (BEND の中の平均 δ = dg/2 → 偏向角 +|dg|/2 倍) ので、**追跡の符号が物理的に正しく、tradke の符号が誤り** (K64 の誤りと判断した。作者の意図の裏づけは無い) |
| v2 の動作 | tradke の `ddpx`・`ddpy` (`cod` の更新と、行列 D・B の微分 `tradke_db` の `dddpx`) の符号を反転 (numba の `elkern._tradke`、Python の `State.tradke`、`perfkern.tradke_db`)。RADTAPER の追跡の BEND の係数も tturne と同じ (3δ_l + δ_l+1)/4 + dptaper にした (K64 の追跡は (δ_l + δ_l+1)/2 + dptaper。この重みの違いも閉じない原因の 1 つで、符号だけ直すと RADTAPER あり x 4e-9) |
| 根拠 | 閉軌道の 1 周の追跡の動き (v2 / k64): 小リング (RADTAPER なし) x 1.4e-9 / −4.35e-7 m、z −3.7e-9 / 2.4e-6 m。小リング (RADTAPER) x 8e-10 / 1.96e-6 m、z −1.6e-9 / −1.08e-5 m。大きなリング A (RADCOD): x −4.5e-11 / −2.0e-8 m、z −1.2e-11 / 5.6e-10 m。大きなリング A (RADTAPER): −3.0e-11 / −1.6e-8 m。PF リング: 2.7e-11 / −6.2e-9 m。ASTRID2 (MPW 入り): −5.4e-9 / −1.17e-7 m。20 周の追跡も閉軌道 (−2.19e-7 m) の周りに止まる |
| 他の量への影響 | 閉軌道そのもの、U0、チューン、ダンピング、エミッタンスが相対 1e-5〜1e-3 変わる。大きなリング A (RADCOD): ε_x 1.54937e-10 → 1.54966e-10 (+1.9e-4)、ダンピング率 +2e-4、U0 296074.8485 → 296074.8816 eV (+1e-7)、チューン +8e-6、閉軌道 x −1.8232e-6 → −1.8257e-6 m。小リング (RADTAPER なし) U0 1563263.2244 (gs) → 1563263.2192 eV。**v2 の Emittance[] (RADCOD) は gs と相対 1e-4 〜 1e-3 程度ずれる (k64 は gs と全桁一致)** |
| 他の試験の変更 | v2 の Emittance[] (RADCOD) が gs と 1e-6〜3e-3 ずれるので、v2 で gs の値と比べていた 5 件の許容誤差を v2 だけ緩めた (k64 は同じ): `test_v2_C_track.py::test_radtaper_matches_gs` (U0 1e-7 → 5e-6、チューン 2e-6 → 2e-5、2 周の追跡の gs との比較は k64 だけ)、`test_v2_C2_optics.py::test_radtaper_in_calc6d` (軌道 1e-4 → 5e-3、DDP 1e-6 → 1e-5)、`test_v2_C2_spch.py::test_wspac_emittance_beam_matrix_ring`・`test_example_wspac_script` (チューンシフト 2e-3 → 6e-3、NOSELFCOD 1e-3 → 3e-3)、`test_v2_A_spin.py::test_tradke_px_py_diffusion_uses_b_y_py_v2` (v2 と k64 の B の他の成分の一致を `==` から相対 1e-10 へ) |
| 試験 | `test_radcod_closed_orbit_closes_in_tracking_v2[RADTAPER なし・あり]`、`test_radcod_closed_orbit_does_not_close_k64`、`test_tradke_python_and_numba_agree_v2` |
| 使える範囲 | RFSW・RAD・RADCOD (・RADTAPER) の閉軌道から追跡を始めると 1 周で閉じる (x・px 1e-9 m、z 1e-8 m。3 GeV、周長 100〜600 m の電子リング、U0/E ~ 1e-4〜5e-4)。RADCOD 付きの Emittance[] を gs と比べるときは k64 を使う |
| 使えない・未確認 | (1) tradke の符号反転の根拠は解析的な導出と閉軌道の閉じ方で、gs の Emittance[] は K64 のまま (v2 の値は gs と合わない)。(2) 1 周の解析 T (RADCOD) と数値微分の差 2e-3 は残る (項目 3)。(3) RADTAPER の重みの変更は「閉じる」ことを根拠にした選択で、どちらの重みが物理的に正しいかは判断していない (**要判断**)。(4) SOL 領域内の `tradke` (bzh ≠ 0) も同じ式なので符号が変わるが、追跡 (`tsolqur` など) との整合は確かめていない。(5) U0/E が大きいリング (8 GeV、周長 1.4 km: 大きなリング B 相当) は確かめていない |

## 項目 8: PSPAC の x と y が gs と 16% 違う

| 欄 | 内容 |
|---|---|
| 状況 | 元の例 (docs/history/v2_changes_C_spch.md §5 末尾: gs の rms x 1.461 mm、y 1.259 mm、v2 1.228/1.226 mm、解析式 1.248 mm。台本が残っておらず、ビームの作り方が不明) を、同じ記述 (5 MeV・1e11・半径 1 mm の一様な円、`D(0.5) SPCH D(0.5) SPCH`、`PSPACNZ=32 PSPACDZ=0.005 PSPACDX=1e-4`) から作り直した (KV 分布 2048 個、準乱数 4096 個)。**gs の x と y の 16% の非対称は再現しなかった**: gs は x 1.2338、y 1.2356 mm (対称)、v2 は 1.2270、1.2282 mm (解析式 1.2481 mm)。gs の x が解析式より 17% 大きかった元の値は、別の乱数のビームか SPC8 (下) の状態の持ち越しによる可能性が高い |
| 調べたこと | gs と v2 を、解析式 (薄いキックの 2 回: rms x = 0.5 mm (1 + s/2)、s = 2 r_e N L/(lb β² γ³ a²)) と、メッシュ (縦 lb/DZ、横 a/DX)・粒子数・エネルギーを変えて比べた (各 gs は別のセッションで取る: SPC8) |
| 結果 (縦のメッシュ DZ = 0.005 m、a/DX = 10) | rms x の解析式からの差 (v2 / gs): lb/DZ = 20: −1.7% / −1.1%、10: −4.0% / −2.8%、4: −12% / −7.8%、2: −27% / −17%、1: −48% / −20%。DZ = 0.001 m (lb = 0.02 m、20 セル): −2.4% / −1.5%、(lb = 0.005 m、5 セル): −9.7% / −2.8%。縦のメッシュを自動 (`PSPACNZ<4`) にすると v2 は 1% 以内 (a/DX ≥ 4)。横: a/DX = 10・4・2・1 (自動の縦): −0.6%・−1.0%・−2.8%・−13% (lb = 0.1 m) |
| 粒子数・γ | 粒子数 500〜4096: 解析式との差 0.3% (gs・v2 とも。128 個は統計で +2〜4%)。5〜50 MeV (s を揃えて): 0.3% (gs・v2 とも)。1/γ には依らない |
| 原因 | **メッシュの分解能** (縦と横)。gs (K64 の spch.c) は電荷を横 CIC (双線形)、縦「ほぼ最近接」(`rhocal` の SPC6: セルの前半は両方の重みを同じノードに足す) で置き、電位を双 3 次補間、λ を 3 次スプラインで補間する。v2 は電荷を横・縦とも 2 次の B スプライン (TSC) で置き、電位を TSC で補間し (自己力の補正・横のぼかしの補正つき)、λ を TSC で補間する。TSC の縦のぼかしは、一様な長さ lb のバンチの両端 (±1.5 セル) の λ を下げるので、lb が数セルのバンチでは λ が小さく出て、キックが小さい。gs の「ほぼ最近接」はぼかしが小さい。両方とも 20 セル以上 (a は 10 セル以上) では 2% 以内 |
| v2 の動作 | 結果は変えない (v2 のモデルが K64 と別の形状関数なので、粗いメッシュの gs との差は残る。v2 は 1 セル以下の粗いメッシュで gs より解析式から遠い)。**粗いメッシュの警告を足した**: 横の rms が PSPACDX・DY の 1.5 セル未満 (キックが 10% 以上小さい)、縦の rms が PSPACDZ の 1.5 セル未満 (PSPACNZ ≥ 4 のとき)、粒子の 1% 以上がメッシュの外 (`???-PSPAC-Warning-…`) |
| 根拠 | 解析式 (上)、粒子数・γ の依存 (上の表)、gs の別セッションの値 |
| 試験 | `test_pspac_resolution_vs_analytic_and_gs` (v2、解析式と gs の値、縦・横の粗いメッシュの警告) |
| 使える範囲 | メッシュが十分細かいとき (a ≥ 4 DX、一様なバンチの lb ≥ 10 DZ または `PSPACNZ<4`) は、解析式と 1〜2%、gs と 1% で合う。使えない: 粗いメッシュの定量値 (警告が出る)、キックが強い (s ≳ 30、ビームがメッシュの外に広がる) 場合 |
| 使えない・未確認 | gs の x と y の 16% の非対称の原因は未確認 (再現しない)。gs は PSPAC の状態を C の static 変数に持つので、1 つのセッションで TrackParticles を繰り返すと 2 回目から PSPAC が働かない (SPC8。例 `A1_pspac_resolution.gs.txt` の 2 行目以降が該当)。gs と直接比べるときは 1 セッション 1 回。K64 の形状関数 (CIC・双 3 次・「ほぼ最近接」) を v2 に入れていないので、粗いメッシュでも gs と同じ値が要るときは別途 (k64 モードに K64 の形状関数の PSPAC を入れる案) |

## 要判断・残した項目

1. **`Emittance[]`・EMIT・CALC6D が BEND の F1・FRINGE を無視する** (項目 4 の発見)。gs と U0 で 5%、チューンで 0.006 ずれる。`tbende` の F1 の部分 (tblfre、tbendef1、f1r・f2r、bradprev) の直訳が必要。K64 のバグではなく pysad の未対応。理由: 項目 4 は追跡の範囲で、`tbende` は別の大きな直訳 (隣り合う BEND の bradprev も含む) なので、この項目では扱わなかった。
2. **tradke の ddpx の符号 (項目 7) は K64 の誤りと判断して v2 で反転した。RADTAPER の追跡の BEND の重みも変えた。** どちらも「閉軌道が追跡で閉じる」を根拠にした選択。K64 の作者の意図の裏づけは無い。gs の Emittance[] (RADCOD) と v2 の値は 1e-4 程度ずれる。
3. 大きなリング B (8 GeV、周長 1.4 km 級) では確認していない (項目 7: U0/E が大きいとき。項目 2: 陽子の大きなリング)。
4. 放射の β の因子 (項目 2)、Touschek・RADINT の古典半径 (項目 2)。
5. DETR ≥ 1 の INS と MEASURE (項目 1)。
6. PSPAC の gs 互換の形状関数 (項目 8)。

## 全試験・回帰

全試験 (`--ignore` の 3 本を除く、`PYSAD_SAD_PACKAGES` つき、他の担当の追加分を含む): 1015 通過・4 skip (A1 の追加 40 件を含む)。回帰 (`QUICK=1`、k64): 50 ベンチ FAIL 0 (145 s)。
