# v2.0.0 の変更記録: 担当 C (追跡・集団効果・要素)

作成: 2026-11-17。この文書は v2.0.0 の作業のうち担当 C の分の記録です (まとめは docs/47 に入る予定)。
MEASURE・TRCOD (B41・B42) は [v2_changes_C_measure.md](v2_changes_C_measure.md)、空間電荷 (SPAC・PSPAC・SPCH) は [v2_changes_C_spch.md](v2_changes_C_spch.md) に分けて書いた。

- 互換モード: `pysad.compat` (`PYSAD_COMPAT=v2|k64`)。v2 (既定) が修正した動作、k64 が K64 (gs) の動作。
  表の「k64 モード」が「同じ」のものは、モードで切り替えていない (gs に合わせた未対応機能の追加、pysad 自身の不具合の修正)。
- 根拠: v2 の値は、解析式・独立な数値計算 (scipy の solve_ivp によるローレンツ力の積分など)・保存則・極限で確かめた。gs (SAD 1.4.4.2k64、WSL) は
  「K64 の動作の確認」と「K64 が正しく動く機能 (LOSSMAP、RADTAPER、ResponseMatrix の TRPT、MULT のスキュー、スピン追跡の損失) の参照」に使った。
- 試験: `pysad/tests/test_v2_C_*.py` (11 ファイル。件数は §4)。走る例: `example/v2_C/` (6 本と出力)。
- 使った条件: 電子。3 GeV・周長 31〜106 m の FODO リング (空洞 1.5 MV か 4 MV)、1 GeV・3 GeV の短い輸送路、RADLIGHT は 10.2 MeV (γ = 20)、
  ソレノイドは 0.3〜3 T/(Bρ) 相当 (`BZ` 2〜3 T @3 GeV、試験のカーネルは bz = 0.3・3 1/m)。実機の大きなリングでは確かめていない。

## 1. 一覧

| ID | K64 の動作 | v2 の動作 | 根拠 (式・数値) | 試験 | k64 モード |
|---|---|---|---|---|---|
| B10 | `tdusrn` は期待光子数 anp が NaN だと npmax = 10000 回まわり乱数を 20001 個引く | anp が NaN・±Inf なら光子を出さず、乱数は 1 個 (通常の「光子なし」と同じ) | 光子数の期待値が定義できない (発散した) 粒子は放射させない。有限の anp では K64 と同じ乱数列 (同じ種で同じ戻り値と種) | `test_v2_C_track.py::test_b10_*` (3) | 10000 回のループを再現 (種が 20001 回進む) |
| B17 | gfortran の MAX/MIN が NaN を無視する前提 (asinz・gfloor・sqrt1) | 乱数の消費は B10 で NaN を明示的に扱う。asinz・gfloor・sqrt1 の NaN の返し方は変えていない | §3 要判断 1 | (B10 の試験) | 同じ |
| B33 | セッション最初の TrackParticles の最初の tapert は l_track = 0 で、LOSSAMPL の外の粒子に落とす印が付かず追跡が続く | 最初の tapert で落ちる (座標は 3·LOSSAMPL に頭打ち、落ちた周 1、落ちた要素は始点の 1 つ前) | 「絞りの外から始まる粒子は落ちる」。2 回目以降の K64 の動作と同じ | `test_b33_*` (2) | 再現 (従来どおり) |
| B34 | BEAMBEAM の傾き回転が cost = sin(傾き角)、sint = 傾き角 (回転にならず座標が潰れる) | (cos t, sin t) の回転 | 回転に対する共変性: 包絡を t 回した強ビームの蹴り = 粒子を −t 回し、回していない包絡で蹴り、+t 回す。t = 0.3 で 1e-9 (蹴りに対する相対) | `test_v2_C_beambeam.py` (2) | 再現 |
| B35 | SOL 領域に DRIFT BEND QUAD SOL MULT CAVI MAP APERT MARK MON 以外があると SIGSEGV | `SolRegionError` (K64 の tsol.f が出すはずの文 `Only DRIFT, BEND, … are supported within SOL: 名前 (型)`)。長さ 0 の SEXT なども黙って通さない | 仕様 (tsol.f の select case の case default) | `test_b35_*` (1) | 従来どおり (NotImplementedError、L=0 は素通し) |
| B37 | UND に磁場のキーワードが無く、磁場なしのドリフト。さらに undulator.f は陰的な刻みの右辺に dpds を二重に掛けていて (新規。下の N1)、磁場があっても収束力が刻み幅に比例して消える | `FBX`・`FBY` [T] を受ける。追跡・CALC の Twiss・Emittance[]・EMIT に同じ写像が入る。RAD で放射 (FLUC も) | 平面アンジュレータの解析式: 軸上の軌道 (出口で x = x' = 0)、行路長 −(F0/kz)² L/4 (1e-6)、垂直の収束 k = F0²/2 (行列要素 5e-4、1/(1+δ)²)、KX≠0 の kx = −F0²KX²/(2kz²)・ky = F0²(kz²+KX²)/(2kz²) (5e-3)、シンプレクティック性 1e-7、放射損失 (2/3) r_e γ³ (F0²/2) L (1e-4)。リング: Δνy = ∫βk ds/4π (1%)、ΔU0 (2e-3)、減衰率の和 = 2U0/E (1e-3)、量子励起は同じ ⟨1/ρ²⟩ の BEND の列との比が ⟨\|1/ρ\|³⟩ の比 1.2004 (2%)、収束は BEND の列と Δν 2e-5 | `test_v2_C_und.py` (9) | 従来どおり (FBX は NotImplementedError、光学は恒等) |
| B38 | PHSROT は行列がすべて 0 になり粒子の座標が 0 になる。PX PY PZ は指定できない | 意図された写像 M = N·Rot(2πPX, 2πPY, 2πPZ)·(1−D)·N⁻¹ + 拡散 (RAD・FLUC)。追跡・CALC・Emittance[] | M0 のシンプレクティック性 1e-13、固有値 = exp(±2πi P)、Courant-Snyder の 2×2 (1e-14)、リング (ST PHSROT) の Twiss・チューン = 入力 (1e-9)、平衡エミッタンス b/(2d−d²) (Emittance[] 1e-6、追跡の統計 8%)、減衰率 = D | `test_v2_C_phsrot.py` (5) | 従来どおり (追跡は NotImplementedError、光学は恒等) |
| B40 | DELCOR が MARK には DETR・DPX、CAVI には HARM・PHI に (符号も逆に) 書く | MARK・CAVI も DX・DY に、ほかの磁石と同じ符号で入れる | 相関長 ≫ 周長なら全要素のずれが同じ値・同じ符号になる (QF・RFC・ST が 5% 以内。相関長 5 km では全桁同じ) | `test_v2_C_delcor.py::test_b40_*` (2) | 再現 |
| B62 | tvcorr の corlength≠0 の枝は波数 −1 の成分 cv(n) が未設定 (2 回目以降は古い記憶域の値)、負の波数の重みが 1 つずれる | 正負の波数で対称に置く (乱数を引く順は同じ) | 相関関数 ⟨c(s) c*(s+Δ)⟩ = exp(−Δ²/(2 x0²)) (実数)。300 標本の平均で 5% 以内、虚部 2% 以内 | `test_b62_*` (1) | 再現 (cv(n) = 0) |
| B41・B42 | MEASURE は無出力、TRCOD は SIGSEGV。ttinit の添字 | MEASURE・TRCOD・LOGMEASURE・MAIN の `TRACK TURNS=-1` を実装 | [v2_changes_C_measure.md](v2_changes_C_measure.md) | `test_v2_C_measure.py` (19) | 従来どおり (無視) |
| B59 | CSRSetupOY の dpjn が zjn と同じ記憶域 (zj に運動量が入り、dpj は 0)。しきい値以上の点が無いと i1 が未初期化 | zj は位置、dpj は運動量。点が無ければ失敗 (CSR::irreg) | 蛙跳びの関係 z_{k+1} − z_k = p_k r56 dt1 (1e-9)、半周期の間 p の符号が同じで両端 (折り返し点) で \|p\| が小さい、z がバンチの範囲 | `test_v2_C_csr.py` (2) | 再現 |
| B65 | RADLIGHT の軌道に横方向の座標が入らず、BEND で向きが変わらず、BEND の中を記録しない | 横の軸を幾何から取り、BEND は円弧 (tlstore の dphi≠0 の式) で進め、中も 1/γ の角度あたり `RADLIGHTDIV` 点 (既定 16) 記録する。曲げのある MULT は出口だけ、CAVI も記録 | 振幅 0 の粒子の点 = 幾何 `LINE["GX"/"GY"/"GZ"]` (ROTATE つき・負の角を含む。1e-11)、時刻 = 行路長/(βc) (1e-9)。偏向磁石の放射場のピーク E = e β²/(4πε0 ρ R (1−β)²) (Liénard-Wiechert) と RadiationField が 0.25% (RADLIGHTDIV 16)、5e-5 (64) | `test_v2_C_radlight.py` (3) | 再現 |
| B66 | 旧式の `WAKE` コマンドは配布に無い WAKECOMMAND.n を読もうとして失敗する (tfwake.f は死んだコード) | tfwake.f の構文 `WAKE L\|T 要素 s1 w1 … ENDWAKE` を WakeFunction の表として登録する (CAVI だけ) | 同じ表を `WakeFunction[…] = …` で与えた追跡 (gs と一致済み) と同じ結果 (2e-9) | `test_v2_C_wakecmd.py` (3) | 従来どおり |
| B69 (docs/46) | SOL 領域内の BEND の偏向を gs がほぼ無視する。pysad v1 は領域内の BEND・MULT・CAVI・APERT が未対応 | BEND = 直線の座標のまま一様な双極磁場 (tsol.f の `tdrift(al, bzs, phiy, phix)`)、MULT = tmulti1、CAVI = tmultiacc を bz つきで。追跡・CALC (INS/CELL) の Twiss と軌道・Emittance[] に入る | ローレンツ力 dp/ds = (p × B)/pz の数値積分 (solve_ivp、rtol 1e-13) と力学的運動量で一致: BEND・K0 だけの MULT 1e-14、K1 つき 1e-10、6 極・八極・スキューは bz = 0 の MULT と同じ分割誤差 (1e-8、EPS で減る)。放射損失 (2/3) r_e γ³ (K0/L)² L (bz = 0・0.4・2 で 2e-4)。bz → 0 の領域内の要素 = 領域の外の要素 (リングの Twiss・Emittance[] が 2e-7)。CALC の軌道 = 追跡 (1e-12)。bz = 3 T のリングで減衰率の和 = 2U0/E (2e-3) | `test_v2_C_sol.py` (30) | 従来どおり (BEND・MULT・CAVI は NotImplementedError) |
| APERT (領域内) | tsol.f は tapert1 を呼ぶ。pysad v1 は未対応 | 領域内の APERT で落ちる (la は変えない) | 12 mm の粒子が ±10 mm の絞りで落ち、LOSSMAP の要素が APERT | `test_sol_region_calc_orbit_equals_tracking_and_apert_drops` | 同じ (追加) |
| LOSSMAP | gs は TrackParticles の戻りに 2 行 (落ちた要素 l_track、落ちた周。生存は 0 と nend+1) を足す。pysad v1 は足さない | 足す (POL のときは 10・11 行目)。l_track はセッションを通じて持つ | gs の例 (`example/05_tracking/flag_LOSSMAP`): 要素 {0,0,242,242,242}、周 {31,31,4,1,1} と一致 | `test_lossmap_rows_match_gs` (2) | 同じ (追加)。v2 は l_track = 0 で落ちた粒子にも落ちた周を返す |
| RADTAPER | 磁石の強さを閉軌道の局所の運動量に合わせる。pysad v1 はフラグを持つだけ | Emittance[]・EMIT・TrackParticles に入る。`Twiss["DDP"]` に閉軌道の δ | gs: 軌道 x −7.887e-5 → −2.198e-7 (pysad −2.173e-7)、U0 1564904.0 (8 桁一致)、εx 1.7372e-7、εz 1.9077e-5、チューン 6 桁、DDP {0, 3.9113e-4, 2.6068e-4} (1e-4)、2 周の追跡 1〜2% (gs の閉軌道の収束判定が緩いため) | `test_radtaper_matches_gs` (2) | 同じ (追加) |
| ResponseMatrix (TRPT) | 応答に sqrt(gammab(stt+1)/gammab(mon)) を掛ける。pysad v1 は未対応 | 掛ける | gs と全桁一致 (1 → 1.2 GeV の加速)。比 = sqrt(1/1.2) | `test_response_matrix_trpt_matches_gs` (2) | 同じ (追加) |
| スピン追跡の損失 | POL + RAD の TrackParticles も tapert で落ちる。pysad v1 は LOSSAMPL・LOSSDZ を見ない (発散は例外) | track6 の経路でも最初の要素・15 要素ごと・区間の終わりに tapert。LOSSMAP の 2 行も付く | gs: 落ちた粒子の要素 74・周 1・座標 .022037636775、生存粒子のスピンが一致 (k64 で 1e-12)。2 回目の呼び出しは全項目一致 | `test_spin_tracking_loss_matches_gs` (2) | 同じ (追加)。最初の呼び出しの B33 の癖 (復活) は近似 (gs は要素 60、pysad は 15) |
| MULT のスキュー (新規 N2) | (K64 は正しい) pysad v1 は MULT の SK0 SK1 SK2 … (ANGLE = 0) を黙って無視していた | `akr(n) = cmplx(K_n, SK_n)·cr1^(n+1)`、cr1 は複素 K1 の主軸 (tpara) | gs: 追跡 4 成分が 2e-13、Twiss (BX BY R1 R2 DX DY) が 1e-9 | `test_mult_skew_components_match_gs` (2) | 同じ (pysad の不具合の修正) |
| SPAC・PSPAC・SPCH | [v2_changes_C_spch.md](v2_changes_C_spch.md) | 同左 | 同左 | `test_v2_C_spch.py` (22) | 従来どおり (空間電荷なし) |

新しく見つけたもの (docs/21 への採番は未):

- **N1 (K64、undulator.f)**: `fgx=dpds*(ff*fx+gg*gx)`・`fgy=dpds*(ff*fy+gg*gy)` は、すでに dpds を掛けた fx・gx・fy・gy にもう一度 dpds を掛けている。
  陰的な刻み p' = p − ds ∂H/∂q(q, p') の右辺は p − (ff·fx + gg·gx) (fx は ds/(1+δ)·∂ff/∂x) で、dpds は 1 回。K64 の式では収束力が刻み幅に比例して消える
  (200 刻みと 2000 刻みで py が 10 倍違う)。v2 の磁場つきのカーネル (`und.und_seg1`) は 1 回だけ掛ける。磁場なし (k64 の経路 `und.k_und`) は 0 を掛けるだけなので影響なし。
- **N2 (pysad v1)**: 上の表の MULT のスキュー。`lattice._mult_params` が `KsL` を読んでいなかった (CALC・Emittance[]・追跡のすべて)。
- **N3 (pysad、未修正)**: `CalculateOptics` は TRPT の加速を入れない (CALC は入れる)。ResponseMatrix の TRPT の試験は、CALC の Twiss から `opt` を作って gs と比べた。
- **N4 (K64、phsrot.f)**: 拡散の乱数は `v = U·√λ` (ベクトル) を成分ごとに正規乱数へ掛けてから N を掛ける。U が対角でないと共分散が N B Nᵀ にならない。v2 は N U (√λ ⊙ ξ)。
- **N5 (K64、tdrift.f)**: ソレノイド + 双極の中の放射 (tradk) は 1 刻みの前後の向きの変化 (弦) から損失を出すので、らせんの回転角 bz·L が大きいと小さく出る
  (相対 (bz·L)²/12。bz·L = 1 で 8%)。v2 の領域内の BEND・MULT は bz·L/ndiv ≤ 0.02 に刻む。

## 2. 項目ごとの使い方・使える範囲・使えない範囲

### 2.1 UND (アンジュレータ) の磁場 (B37)
```
UND U1 = (L=1 FBX=1.8 POLE=20 SLICE=2000);      ! FBX: 垂直磁場のピーク [T] (x 方向に蛇行)、FBY: 水平磁場のピーク [T]
```
- 軸上の磁場は By = (FBX·c/p) cos(kz s)、Bx = (FBY·c/p) cos(kz s + PHI)、kz = 2π POLE/L。`KX`・`QY` は横方向の波数 (K64 と同じ)。
- 追跡 (`TrackParticles`、MAIN の TRACK、DA: track6n)、CALC の Twiss、Emittance[]・EMIT・CALC6D・Touschek (写像は追跡のカーネルの数値微分、放射は 1 周期を 64 個以上に分けて tradke)。
- **v2 では磁場がなくても UND は長さ L のドリフトとして光学に入る** (K64 と v1 は光学では恒等で、追跡だけがドリフトさせる)。`SLICE=0` は追跡も光学も恒等。
- 条件: `SLICE` は `POLE` の 20 倍以上 (刻みは 1 次の陰的な蛙跳び。刻みが周期の整数倍だと sin(kz s) = 0 で磁場が見えない)。試験は 3 GeV、1 T・1.8 T、λ = 5 cm、20 周期。
- 使えない: track6 (POL + RAD のスピン追跡) の経路、`BX`・`BY` という名前 (K64 にも無い。ValueError)、端の磁極の整形 (軌道は x' = −(F0/kz) sin(kz s) で、中心が −F0/kz² ずれる)。
  PHI≠0 の入口では正準運動量をそのままつなぐので、ベクトルポテンシャルの飛びが蹴りになる (K64 の式のまま)。

### 2.2 PHSROT (B38)
```
PHSROT PR = (AX=0.5 BX=4 PX=0.31 AY=-0.3 BY=2 PY=0.22 AZ=0 BZ=10 PZ=0.01 D11=1e-4 … B11=2e-12 …);
```
- PX PY PZ はチューンの単位の位相、Dij は正規座標での減衰行列 (RAD ON で効く)、Bij (i ≤ j) は 1 回あたりの拡散行列 (RAD かつ FLUC)。R1..R4、EX EPX EY EPY ZX ZPX ZY ZPY も使える。
- 追跡 (track6n)、CALC (シンプレクティックな M0)、Emittance[]・EMIT・CALC6D・Touschek (T ← M0 T、D・B に減衰・拡散)。
- 注意: (z, δ) の回転は向きが逆 (matrix.f m_prot_set_6d。縦のチューンは −PZ と出る)。分散つきの Twiss を CALC で読むと、4 次元の分解なので BX・EX は入力と 1e-3 ほど違う (H 行列の定義の差)。
- 使えない: track6 の経路 (POL + RAD)、k64 モード。拡散行列に負の固有値があると ValueError。

### 2.3 ソレノイド領域内の BEND・MULT・CAVI・APERT (B69、B35)
- BEND: `L`、`ANGLE + K0` (蹴りの大きさ)、`ROTATE + DROTATE` (向き)、`RAD`。**直線の座標のままの一様な双極磁場** (幾何は曲がらない)。`K1`・`DX`・`DY` は NotImplementedError (四極成分は QUAD か MULT で)。
- MULT: K0..Kn と SK0..SKn、`DX DY ROTATE`、端のフリンジ (bz = 0 の式)、RAD、`VOLT` (RFSW か TRPT)。本体は「ソレノイド + 四極」か「ソレノイド + 双極」の厳密解 (K0 と K1 が両方あるときは四極の中心のずれとして厳密)、ほかの多極は薄い蹴り。
- CAVI: `L`、`VOLT`、`PHI`、`DX DY ROTATE`、端のフリンジ。Emittance[] では [ソレノイドのドリフト半分 → 薄い空洞 → 半分] に分ける。
- APERT: 領域の外と同じ判定。
- CALC の Twiss と閉軌道 (solmap: 追跡の写像の数値微分)、Emittance[]・EMIT (solemit) に入る。CALC6D・RADCOD つきの CALC は SOL があると従来どおり SADError。
- 使えない (明確なエラー): TRPT の加速 (VOLT·sin(PHI)≠0 の CAVI・MULT。基準運動量が領域内で変わる)、wake のある CAVI・MULT、`AUTOPHI`、`V1 V20 V11 V02`、ANGLE≠0 の MULT、RF つき MULT の Emittance[]、
  SEXT・COORD・BEAMBEAM など (SolRegionError)、SOL の `DX DY DZ CHI1..3`・`F1` (従来どおり未対応)、track6 の経路 (POL + RAD のスピン追跡では領域内は DRIFT・QUAD だけ)。
- 未確認: 領域内の MULT・CAVI の POL、RF つき MULT の RAD、bz ≥ 3 1/m を超える場。

### 2.4 LOSSMAP・B33・スピン追跡の損失
- `LOSSMAP;` で TrackParticles の戻りの flag の後に 2 行 (落ちた要素の番号、落ちた周) が付く。落ちた要素は判定 (tapert) のときの l_track なので、15 要素ごとの判定点か APERT の番号。
- POL + RAD の TrackParticles (track6 の経路) も LOSSAMPL・LOSSDZ で落ちる。発散した粒子は例外にせず落とす。
- 使えない: MAP 要素 (ExternalMap) が落とした粒子の要素番号 (0 になる)、空間電荷 (spch) が落とした粒子の要素番号。POL なしで高速追跡に乗らない要素の並びは `LOSSMAP: … 記録できない` の SADError。

### 2.5 RADTAPER
- `RFSW; RAD; RADCOD; RADTAPER;` のすべてが ON のとき、Emittance[]・EMIT が磁石の強さ (BEND の K0+ANGLE・K1、QUAD、SEXT〜DODECA、MULT) を閉軌道の局所の運動量に合わせる。続く TrackParticles も同じ係数を使う (直前の Emittance[] の結果)。
- 実装は K64 の tcod の反復の直訳ではなく、同じ不動点 (δ(1) = DP0、一様な係数 dptaper を未知数にして 6 次元で閉じる) を外側の反復で求める (`radtaper.py`)。gs は収束判定が緩いので、軌道の残り (0.2 μm) は 1% 違う。
- 使えない: CALC (4 次元の光学) のテーパー (K64 の qtwiss の radtaper は未翻訳)、SOL 領域内の要素、RFSW OFF。見つからないときは `???-RADTAPER-closed orbit with taper not found` を出してテーパーなしで続ける。

### 2.6 RADLIGHT (B65)
- `RADLIGHT;` の TrackParticles の軌道が、v2 では実際の 3 次元の軌道になる。`RADLIGHTDIV = n;` で BEND の中の点数 (1/γ の角度あたり) を変える (既定 16)。
- 中の点は、状態の写しを刻んだ BEND で進めて作る (本当の追跡は要素全体を 1 回)。DX・DY・F1・FRINGE のある BEND は中の点だけ 1e-6 m 程度ずれる。
- 使えない: SOL 領域を含む区間 (従来どおり SADError)、COORD による幾何の変化、曲げのある MULT の中の点 (出口だけ)。高いエネルギーでは点数が γ に比例して増える (3 GeV・0.26 rad の BEND で 2.5 万点)。

### 2.7 そのほか
- DELCOR (B40・B62): 使い方は同じ。`corlength≠0` の rms は実部・虚部それぞれ 1/√2 (K64 の規格化 1/sqrt(4Σw²+2) のまま。§3 要判断 2)。
- 旧式の WAKE (B66): `WAKE LONGITUDINAL|TRANSVERSE 名前 s1 w1 s2 w2 … ENDWAKE;`。**WAKE から ENDWAKE までを 1 行に書く** (pysad は文を行ごとに読む)。CAVI だけ。追跡は WakeFunction と同じ (TRPT と LWAKE/TWAKE)。
- CSRSetupOY (B59): 戻りの 4 番目 (zj) が位置、5 番目 (dpj) が運動量になる。`CSROYMatrix` はそれを使う (gs が wrongtype を返すので行列自体は未照合のまま)。
- BEAMBEAM (B34): 強ビームに x-y 結合 (R1..R4 か R13 など) があるときだけ結果が変わる。

## 3. 要判断 (直していない・決めきれないもの)

1. **B17 (NaN と MAX/MIN)**: v2 で NaN を明示的に扱ったのは乱数の消費 (B10) だけ。`asinz(NaN) = −π/2`、`gfloor`、`sqrt1` の下限は、発散した粒子が落ちるまでの座標の埋め方にだけ効き、
   有限の粒子の結果は変わらない。v2 で「NaN になった粒子をその場で落とす」に変えると、落ちる要素・周が K64 と変わる (LOSSMAP の値に出る)。物理としてはどちらも「落ちた粒子」なので変えなかった。
   マッチングの fuzz (tffsmatch) は担当外。
2. **DELCOR の規格化**: corlength≠0 のとき DX の rms は delx/√2 (複素数の絶対値の rms が 1)。corlength = 0 の枝 (ランダムウォーク、1 m あたり分散 1) とは意味が違う。意図が読み取れないので K64 のまま。
3. **DELCOR の CAVI の符号**: K64 は (設計の DX) − delx を書く (SOL と同じ形の式)。v2 は他の磁石と同じ +delx にした (相関のある地盤の動きとして、隣の四極と同じ向きにずれるのが自然)。SOL の式を写した意図があるなら逆。
4. **B33 の落ちた要素**: v2 は「始点の 1 つ前の要素 (周回で数える)」とした。K64 の 2 回目以降は「前の追跡の最後の要素」。
5. **UND の光学**: v2 は磁場なしでも光学にドリフトとして入れた (追跡と揃える)。K64 の「光学では恒等」を仕様とみるなら、磁場つきのときだけにする選択もある。
6. **スピン追跡の損失の k64 の最初の呼び出し**: K64 の B33 の癖 (落とす印 0 の粒子が後で入れ替えで戻る) は track6 の経路では近似 (gs 要素 60、pysad 15)。高速追跡 (track6n) は従来どおり全桁一致。

## 4. 試験

| ファイル | 件数 | 内容 |
|---|---|---|
| `test_v2_C_beambeam.py` | 2 | B34 |
| `test_v2_C_und.py` | 9 | B37 (カーネルの解析式 6、スクリプト 1、リング 2) |
| `test_v2_C_phsrot.py` | 5 | B38 |
| `test_v2_C_delcor.py` | 3 | B40・B62 |
| `test_v2_C_csr.py` | 2 | B59 |
| `test_v2_C_wakecmd.py` | 3 | B66 |
| `test_v2_C_radlight.py` | 3 | B65 |
| `test_v2_C_track.py` | 14 | B10 (3)、B33 (2)、LOSSMAP (2)、B35 (1)、RADTAPER (2)、ResponseMatrix の TRPT (2)、スピン追跡の損失 (2) |
| `test_v2_C_sol.py` | 32 | SOL 領域 (ODE 18、放射 3、空洞 1、リング 6、輸送路 1、k64 1)、MULT のスキュー (2) |
| `test_v2_C_measure.py` | 19 | MEASURE・TRCOD (別記) |
| `test_v2_C_spch.py` | 22 | 空間電荷 (別記) |
| 計 | 114 | |

## 5. 変更したファイル

- 新規: `pysad/phsrot.py`、`pysad/radtaper.py`、`pysad/spch.py` (別記)、試験 11 本、`example/v2_C/`、この文書と別記 2 本。
- `pysad/beambeam.py` (B34: `apply`)、`pysad/track6n.py` (`make_rs`・`tdusrn` の B10、`turn` の LOSSMAP `lmap`、領域内の APERT、PHSROT・UND の分岐、`TrackRun.lossmap`)、
  `pysad/track6n_ext.py` (`tdrift_bz1`・`k_sol_bend`・`sol_mult_full`・`_compile_sol_v2`、`step_kind`)、`pysad/track6.py` (`SolRegionError`、`_tapert_loss`、`_tturn1`・`track_particles` の損失)、
  `pysad/sadrandom.py` (B10)、`pysad/und.py` (B37)、`pysad/solmap.py`・`pysad/solemit.py` (UND・PHSROT・領域内の BEND/MULT/CAVI)、`pysad/csr.py` (B59)、
  `pysad/lattice.py` (`_mult_params` のスキュー)、`pysad/calc6d.py`・`pysad/touschek.py` (UND・PHSROT の呼び出し)、`pysad/maintrack.py`・`pysad/dapert.py`・`pysad/track6n_ext.py` (`make_rs`)、
  `pysad/measure.py`・`pysad/parser.py` (別記)、
  `pysad/ffscmd.py` (B33、LOSSMAP、DELCOR、WAKE、`_RadLight`、RADTAPER、ResponseMatrix の TRPT、スピン追跡の損失、UND・PHSROT の annotate、MEASURE・空間電荷のつなぎ)。
- 触っていない: `pysad/sadscript.py`、既存の docs・README・CHANGELOG、他の担当の試験。

## 6. 残っているもの (担当 C の範囲で未対応)

- SOL の `DX DY DZ CHI1..3`・`F1` (trots・tradks)、領域の `tdgeo`。調べた範囲: tsol.f の trots の呼び出しと引数は読んだが、翻訳・検証はしていない。
- 領域内の BEND・MULT・CAVI の track6 の経路 (POL + RAD のスピン追跡) と、領域内の TRPT の加速。
- CALC (4 次元) の RADTAPER、`CalculateOptics` の TRPT の加速 (N3)。
- `TRCOD FIXCOD/MVCOD/STNCOD/ALIGN`、`PUTMEASURE`・`GETMEASURE` (別記)。`MULT` の `RADIUS` の空間電荷、`WSPAC` (別記)。
- UND・PHSROT の track6 の経路。UND の端の磁極の整形。
