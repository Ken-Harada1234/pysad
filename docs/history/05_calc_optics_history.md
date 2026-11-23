# 05_calc_optics.md から取り除いた履歴の記述

05_calc_optics.md (現在の仕様の解説) から、版・経過に結びついた記述を取り除いた。下は取り除く前の元の行を、節の見出しつきでそのまま並べたもの (行番号は整理前のもの)。

## # 05. CELL/INS の光学計算 (CALC) — SAD の qcell/qcod/qtwiss と pysad の lattice.py

整理前の 3 行目:

````text
> 最終更新: 2026-11-18 (v2.1.0 の後) (現状のスナップショット)。SAD ソース (`<SAD>/src`) と [`pysad/lattice.py`](../pysad/lattice.py), [`pysad/ffscmd.py`](../pysad/ffscmd.py), [`pysad/match6.py`](../pysad/match6.py) を照合して記述。行番号は 2026-10 時点。
````

## ## 2. 要素ごとの CALC 用の写像 (`qtwiss.f:184-400` の分岐)

整理前の 81 行目:

````text
**SAD の癖 (k64 モード)**: tchge の入口で th≠0 かつ回した座標系のずれ≠0 (または dz≠0, dtheta≠0) のとき、tbrote が行列を上書きするので行列に回転が入らない (COD は回る)。一周行列が R^T·M になり結合が物理とずれる。`_tchge_jac` が再現 (docs/11 A25)。**v2 モードでは直した** (B09): 入口の回転を行列に掛ける (行列 = 軌道の写像の微分)。`_tchge_jac` が `_tchge_cod` の数値微分と 1e-8 で一致し、ロール + ずれのある6 極・BEND を含むリングの CALC のチューンが追跡の一周写像の固有値と一致する (6 極 1e-9、BEND 4e-6。K64 は 0.02〜0.18 ずれる。[`test_v2_A_tchge.py`](../pysad/tests/test_v2_A_tchge.py))。同じ tbrote のスピン (srot) の更新も v2 は (I+rr)·srot の行列積にした (B56。直交性 1e-14)。
````

## ## 4. 安定性判定

整理前の 136 行目:

````text
| qmdiag | `(TrA−TrB)² < 0` → `***qmdiag---> Sum resonance: …` (lfno=0 で次の表示まで 1 回)、`stab=.false.` (twiss(lb) は以前の値のまま); NaN も不安定。trab2 の絶対値が 1e-14 以下は特殊 6 分岐 | `lattice._qmdiag` (全分岐の直訳)。通常の場合は従来の `periodic_optics` (ビット単位で不変) |
````

## ## 4b. CALC6D と新しい要素の Twiss の範囲

整理前の 147 行目:

````text
- **CALC6D** (`qcell6d`、列 22–28): `calc6d.py`。`CELL`/`INS` に関わらず `qcell6d` に入る (1 節)。SOL を含む CALC6D と、RADCOD つきの SOL の CALC は v2.1.0 で対応した (両モード。k64 は gs と 1e-9、docs/29・44)。CELL の CALC6D のあと `DP0` は先頭の `DDP` (閉軌道の δ) の別名になり、RADCOD の CALC6D のあとの 4D の CALC はそのエネルギーで光学を計算する (SOL なしの RING でも gs と一致: チューン 2.70642542976106、DP0 なしの値は 2.70675785)。
````

## ### v2 モードで変わった CALC・光学の動作 (k64 モードは K64 のまま)

整理前の 166 行目:

````text
モードによらない修正 (pysad の不具合): **`CalculateOptics` が TRPT の加速を入れなかった** (CALC は入れる) のを直した。TRPT かつ RFSW の加速を区間に入れる (区間の途中からでも)。gs と終端の BX AX BY AY EX EPX NX NY が全桁、CalculateOptics − CALC = 0 (`test_calculate_optics_includes_trpt_acceleration`)。また、加速区間の分散の qmat2twiss の gr (= gammab(出口)/gammab(入口) を入口の分散に掛ける) の向きが逆だった (入口に分散のある加速区間の EX EPX が gs と合わなかった: 0.1109 / 0.0734) のを直した。gs と全桁 (薄い空洞・厚い空洞。`test_dispersion_through_accelerating_cavity_matches_gs`)。
````

整理前の 168 行目:

````text
**ソレノイド領域の CALC** (docs/29・44): v2 モードは領域内の BEND・MULT・CAVI・K1=0 の四極、ずれて傾いたソレノイド (`DX DY DZ CHI1..3`・`F1`) を Twiss と閉軌道に入れる (solmap)。領域内の CAVI・RF つき MULT の TRPT の加速は CALC に入る (追跡と 2e-9)。**CALC6D・RADCOD つきの CALC は SOL 領域があっても走る** (v2.1.0。両モード。領域の写像は `solmap`・`solemit`。SOL 領域のあるリング、3 GeV、BZ = 2 T、回転四極、QFE の DX・DY: CALC 4D 2.73475252978800・1.87531224259842 (pysad k64 2.73475252978794)、CALC6D 2.73495437251036・1.87554377922481・NZ −.01173531080122、RADCOD つき CALC6D 2.73497629185364 (1e-10)、β も 1e-9)。**RADTAPER つきの CALC6D** は両モードで入る (Emittance[] と同じ閉軌道・係数。k64 も入れる: gs の CALC6D + RADTAPER はテーパーを入れ、チューンが Emittance[] の値と同じ (gs と 8e-8)。**4 次元の CALC は K64 と同じくテーパーなし**: gs で RFSW・RAD・RADCOD・RADTAPER ON の CELL CALC (4D) が RADTAPER なしの CALC と全桁同じ。実装しない)。`test_radtaper_in_calc6d`・`test_v2_B1_optics.py::test_calc6d_radtaper_matches_gs_k64`。
````

整理前の 169 行目:

````text
- **`Twiss["DETR"]`** (N12。v2.1.0 で修正。両モード): INS の CALC の後で 0 になっていた (`ffscmd._rows_from_tws` が `detr` を渡していなかった)。小さな結合リング (45 度回した四極 K1=0.8)、大きなリング A (QDA2 を 0.03 rad 回し空洞に DX=2e-4・DY=−3e-4: `Twiss["DETR", {1, 500, 1000, 2000}]` と R1 が gs と 10 桁)、PF リングで gs と全桁。DETR = R1R4 − R2R3 の恒等式が 1e-12。`MEASURE` の初期分布が読むのは先頭 MARK の生の値 (CELL CALC が書く) で、元から正しかった。**DETR ≥ 1 (モードが入れ替わる) の INS と MEASURE は未確認** (DETR < 1 のリングは R1 = 0.3 まで入れ替わらなかった。入口を DETR > 1 の生の値で与える輸送路は gs と pysad で入口の規格化が違い比較できなかった)。
````

整理前の 170 行目:

````text
- **ずれた BEND の CALC (B114)・K1≠0 の BEND の端 (B117)**: v2.1.0 で v2 は追跡と整合させた (下の「v2 モードで変わった CALC」の表)。k64 は K64 のまま。K1≠0 の BEND の CALC の解析行列は y ≠ 0 の軌道で自身の写像の数値微分と y = 2 mm で 1.5e-6 違う (tbfrie の解析行列の近似)。CALC の解析行列と `tbende` は DISFRIN=1 でも Maxwell のフリンジを入れる (追跡の `tbendi` は DISFRIN に従う。B134)。
````

整理前の 171 行目:

````text
**WSPAC (空間電荷) つきの CALC・CALC6D** (v2.1.0。**v2 だけ**。K64 は入口以外の Twiss が NaN になる: B133): `WSPAC` ON で CELL の CALC6D と、RFSW ON の 4D の CALC (6D の経路 `calc6d(sc=…)` で、反復の最終パスの T_{1→l} と正規形から Twiss) に空間電荷のチューンシフトが出る。0.2 GeV の 12 セル FODO で、一次のラスレットの積分と 0.15%、Emittance[] のチューンと 2e-6、BX($$$) = BX(1) (周期的)。使えない・未確認: RF なしの CALC (σz が決まらない)、強い空間電荷 (0.2 GeV で PBUNCH ≳ 1e9 は不安定)、K64 の `qwspac` (中点間のレンズ) の直訳ではない (gs が NaN なので検証できない)。k64 モードは従来どおり空間電荷を入れない。詳細は docs/07。
````

## ## 6. 対応状況と既知の未対応

整理前の 200 行目:

````text
- **SOL の CALC の残り**: k64 モードでは領域内の BEND・MULT・CAVI (NotImplementedError → 警告して直進扱い)、領域内の L=0/K1=0 の四極、ずれ・傾きのある SOL が未対応 (v2 モードは対応。上の「v2 モードで変わった CALC」)。両モード共通で未対応: BOUND=0 で BZ≠0 の SOL (NotImplementedError)。(CALC6D・RADCOD つきの CALC は v2.1.0 で対応: 4b 節)。RADTAPER + SOL は未確認。v2 の領域内の TRPT 加速では、`V1 V20 V11 V02`、ANGLE≠0 の MULT が使えない (docs/44)。
````

整理前の 205 行目:

````text
- **解析行列がない要素** (数値微分): ずれ・ロール付き BEND (DROTATE≠0)、MULT、F1 付き BEND、ANGLE=0 の BEND (四極と薄肉多極は解析行列。2 節)。F1 付き BEND は CALC に元から入る (gs と 1e-9)。ANGLE=0 の BEND (ステアリング)・MULT の F1 は v2.1.0 で CALC・Emittance[] に入った (docs/07)。
````

