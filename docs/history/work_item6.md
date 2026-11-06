# 作業ログ: 項目 6 — disp / emit の出力差、emit の高速化、Emittance[] の使い方 (2026-10-16)

担当範囲: `pysad/ffscmd.py`, `pysad/emit6d.py`, `pysad/lattice.py`。参照: work_script_compat.md (disp の差 3 点)、07_emittance.md、03_ffs_commands.md、work_match_fixes.md (emit 8〜12 s)。
SAD ソース: `tfdisp.f`, `temit.f` (temitprint, tecalc)、`tcod.f`、`tturne.f`、`photons.f` (tradke)、`tbende.f`、`tqfrie.f` ほか。
参照実行は WSL の `gs` (ラティス `pysad/tests/kek-ls_v03_12.lat` を**連結した 1 本のファイル**で渡す。`gs lat script` と 2 引数で渡すと SAD は標準入力待ちで止まる)。

## 0. 結果の要約

| 項目 | 結果 |
|---|---|
| A1 反転配置の要素名の `-` | `disp` が `-BPI7.2` 等 380 要素に `-` を付ける。名前から s(m) までの桁位置まで SAD と全 780 行一致 (ヘッダは 66 行ごと 12 回、末尾に余計な `$$$` 行は出さない)。`LINE["NAME"]`・`LINE["ELEMENT"]` は SAD も `-` を付けない (pysad も付かない)。`LINE["DIR",k]` は反転で −1 (以前は常に +1) |
| A2 `Beam tilt` 行 | `Beam tilt` と、同じ入口ビーム行列から決まる `Beam size xi/eta` の行を出す。xi・eta は SAD と一致 (.04093289 / .00000000)。tilt は 0 同然の量 (SAD 9.9e-13、pysad 3.9e-13) で値は一致しない |
| A3 `Equilibrium dz` | SAD 20.618669 mm、pysad **20.618669 mm** (旧 20.619466)。原因は pysad の閉軌道 (trf0 固定の Newton) と SAD の tcod (z を trf0 に移しながら反復) で z と trf0 の分け方が違っていたこと。tcod.f を直訳して trf0 が SAD と 12 桁一致。Bunch Length も 3.77333112 mm で SAD と一致 (旧 3.77333234) |
| B 高速化 | ESRF (3563 要素) の `emittance6d`: RFSW+RADCOD **22.7 s → 1.6 s**、RFSW のみ 10.3 → 1.6〜1.9 s、RF なし 8.9 → 1.3〜1.8 s。`disp` スクリプト全体 25 s → 3.4 s。16b は全桁 ≤1.3e-11 のまま。RADCOD は SAD と 2e-11 以内に一致 (旧 1e-7) |
| C Emittance[] の使い方 | esrf スクリプトで実行される 9 通りの使い方 (`Print[Emittance[]]`, `a=Emittance[]`, `x/.Emittance[]`, `(Emittances/.Emittance[])[1]`, `Emittance[OneTurnTransferMatrix -> True]` …) を SAD と並べて一致を確認。SaveEMIT の既定でも EMITX.. の扱いは SAD と同じ。影響なし |

回帰: pytest 85 本 (新規 test_emit_jac.py 22 本) 全 pass、`run_regression.py` (全 28 ベンチ) FAIL 0 (§5)。

## 1. 項目 A: disp の出力差

### A1. 反転配置の `-`

- SAD: `tfdisp.f:438-448`。`direlc(l) <= 0` (`cmp%ori=.false.` = LINE の `-` で反転) の要素は名前を `' -'//name` にし `max(10, lname+2)` 桁に、そうでなければ `lname<=8` で `'  '//name` を 10 桁、それより長い名前は `' '//name`。
  `-` を付けるのは DISP の表だけ。`tftwiss.f:858` の `LINE["NAME"]` は `elname`(`elname.f`) が名前だけを返し、`LINE["ELEMENT"]` も `pnamec`。`Twiss["NAME",k]` は SAD ではエラー (`General::wrongval`)。
- pysad (`ffscmd.disp_cmd`): 反転は `el._rev` (parser が `-` で逆向きに通すとき立てる。以前の実装は存在しない `_reversed` を見ていたので一度も出なかった)。
  名前欄の幅を SAD の式に合わせた。副産物として DISP の他の差も直した: (i) 以前は名前の後ろが 1 桁広かった、(ii) s(m) 欄が `10.6` で 100 m 以上の位置が `100.21800` と 1 桁欠けていた (SAD は `11.6` で `100.218000`)、
  (iii) Value 欄は `10.7` (負の値も 7 桁。以前は 9 桁で `-.300000` と 1 桁欠けた)、(iv) 値が 0 の欄は `tdtrimz` で `.0` (s=0 の MARK)、(v) 一致した行だけを出す (以前は末尾に `$$$` を常に足していた)、
  (vi) ヘッダは 66 行ごとに繰り返す (`nlc=66`、SAD の出力は 12 回)、(vii) 名前が 10 桁を超えるときの空白の詰め (`tfsqueezespace`、後ろの欄から lname−10 個、足りなければ前の欄から。`ffscmd._squeeze_space`)、
  (viii) `disp *;` (後ろに何もない `*`) が SADScript の式 (積) と判定されて何も出さなかった → `sadscript._is_ffs_command` で `*` 1 字はワイルドカードに (`free *` なども同じ)。
- `LINE["DIR",k]` (tftwiss.f:942 `direlc`): 反転で −1。
- 照合: `<work-esrf>/disp` (`cd <work-esrf>; python -u ../run_sad.py ../pysad/tests/kek-ls_v03_12.lat disp`) を `disp.sad.out` と比較 (`probes_item6/compare_disp_columns.py`)。要素名から s(m) までの文字列が **780 行すべて一致** (修正前 654 行が不一致)。
  AX…DetR の全値を 5 桁の表示のまま比べると最大相対差 0 (11700 値)。

### A2. `Beam tilt`

- SAD: `temit.f:1126-1141` (`temitprint`、書式 9102)。`xxs=emitp(1)-emitp(6)`, `yys=-2 emitp(4)`、`btilt=atan2(yys,xxs)/2` (どちらかが 0 なら 0)、`sig1=|emitp(1)+emitp(6)|/2`, `sig2=hypot(emitp(1)-emitp(6), 2 emitp(4))/2`、
  `sigx=max(√(sig1+sig2), √|sig1−sig2|)`, `sigy=min(…)`。`emitp` は入口の平衡ビーム行列 (パックの 1=σxx, 4=σxy, 6=σyy)。
- 同じ報告の `Closed orbit:` も SAD と同じく Entrance (tcod の cod) と Exit (最後の一周の出口) の 2 行にした (`tput.f`: 各値は `' '//autofg(v,'9.6')`、`emittance6d` が `orbit_entrance`/`orbit_exit` を返す)。値の書式は SAD と一致 (`-1.82E-6 -2.41E-8 … -1.14E-6 9.799E-6` / Exit は `9.798E-6`)。
- pysad: `emit6d.beam_tilt_sizes(emitp)` (tecalc が返す `beam_matrix` から)。`print_emit` が `Bunch Length … mm  Beam tilt … rad` と `Beam size xi … mm  Beam size eta … mm` の 2 行を出す。閉軌道が見つからない (NaN の報告) ときは NaN。
- 照合: `Bunch Length = 3.77333112 mm`、`Beam size xi = .04093289 mm`、`Beam size eta = .00000000 mm` は SAD と一致。`Beam tilt` は 0 同然の量 (SAD 9.9e-13、pysad 3.9e-13) なので桁が違う。
- **εy / eta の雑音の原因**: ヤコビアン・D・B は y が厳密に 0 (結合なし) だが、`tecalc` の平衡ビーム行列を解く `_tsolva` が numpy の SVD (特異値の硬い打ち切り) で、y 方向に εy~1e-21 の雑音を作っていた (SAD は `tsolva` → `tsvd` で ~1e-38)。
  SAD の `tsvd` の直訳 (`pysad/tsvd.py`、マッチングで移植済み) に差し替えた (`_tsvd_tsolva`) → εy = -1.8e-30、eta = 0、他の量は 1e-13 以内で同じ・時間も同じ。

### A3. `Equilibrium dz` (20.618669 mm vs 20.619466 mm)

- 報告の `dz` は `trf0*1e3` (`temit.f:1034`)。`trf0` は空洞の同期位置 (位相 `phis = trf0·w`)。RADCOD ON では `tcod.f` が決める。
- **SAD の tcod (RFSW ON, RADCOD ON)** (`tcod.f`): まず RFSW OFF で自分自身を呼び (4D の閉軌道 + trf0)、`cod(5)=asin((U0 pgev − vcacc)/vceff)/wrfeff − trf0` とおいてから RFSW ON の本体に入る。本体は Newton 法 (`tsolvg` = 特異値 1e-8 の SVD、ヤコビアン T−1)。
  RADCOD のときは**一周ごとに z の平均 `dz=(codi(5)+codf(5))/2` を `codi(5)` から引き、`trf0 += dz/v0`** (RF の位相の基準 z+trf0 は保たれる)。収束判定は `r=Σ(dcod1/codw)² < 1e-4` (codw = 1e-6,1e-7,1e-7,1e-8,1e-5,1e-6)。
  判定が緩い (z は 1e-7 まで) ので、z と trf0 の最終値は反復の履歴 (どこで止まったか) で決まる。
- **旧 pysad**: trf0 を RADCOD OFF の手順 (z=0 の閉軌道から) で固定し、6D の固定点を厳密に Newton で解いた。物理量は z+trf0 の組だけで決まるので emittance などは近いが、z と trf0 の分け方が違い、報告の dz も SAD と 8e-4 mm ずれた。
- **修正**: `emit6d._tcod_radcod` が tcod.f を直訳 (入れ子の RFSW OFF → 本体、`fact` による減速・`fmin`、`dpthre` の判定、`sign(1,trf0)·mod(…)` の折り返しまで)。初期値は SAD と同じく CALC の閉軌道 (4D、z=0、dp=0)。
  結果 (ESRF, RFSW RADCOD): trf0 = 0.02061866961998 (SAD `Emittance[]` の EquilibriumPosition と 12 桁以上一致)、z = -1.14207e-6 (SAD の出口 cod -1.14198e-6)。反復の履歴: RFSW OFF の本体 3 回 (r=1.9e5, 2.7e-3, 5e-9) → RFSW ON の本体 2 回 (r=3.6, 1.9e-6) で収束。
- **同時に直した表示**: `Revolution freq. f0` は `c·β/C` (β=p/h) で、以前は β を落としていて 4e-8 ずれ (497069.34 → SAD 497069.32)。`dz` も `f0` も `autofg('10.7')` の桁 (切り捨てに見える) に合わせて `fmt_real(x,'10.7')`。
- Bunch Length の 2e-5 (SAD 3.77333112 vs 旧 pysad 3.77333234、実際は相対 3e-7) の主因は dz ではなく、数値微分の雑音 (§2) だった。dz の修正 (3.7733306) だけでは 1.4e-7 残り、解析ヤコビアン導入後に 3.77333112 と一致した。
- 全桁照合 (`bench_item6_emit.sad`、SAD の `Emittance[Matrix->True]`、RFSW+RADCOD):

| 量 | pysad | SAD | 相対差 |
|---|---|---|---|
| εx | 1.5493685877683e-10 | 1.54936858777e-10 | 1.1e-12 |
| εz | 2.6213050634472e-6 | 2.621305063471e-6 | 9e-12 |
| σE | 6.9469480483149e-4 | 6.946948048315e-4 | 2e-14 |
| σz | 3.7733311172975e-3 | 3.77333111733e-3 | 9e-12 |
| α_c | 2.3062643821984e-4 | 2.306264382157e-4 | 2e-11 |
| U0 [eV] | 296074.84850396 | 296074.84850396 | 0 |
| trf0 [m] | 0.02061866961998 | 0.02061866961998 | 1e-14 |
| νx, νy, νz | | | 1.5e-12, 1.6e-12, 9e-12 |
| 減衰率 | | | ≤8e-13 |
| ビーム行列 σxx, σxδ, σzz, σδδ | | | 4e-11, 2e-11, 2e-11, 1e-13 |

  (旧実装は εx 5e-7、σxx 5e-6 など 1e-7〜1e-5。docs の「RADCOD ON は 1e-7 (SAD の tcod 収束誤差)」は誤りで、原因は pysad の数値微分の雑音だった。)

## 2. 項目 B: emit の高速化

### B-1. プロファイル (旧、`emittance6d` を RFSW+RADCOD で)

- 1 回の `cell calc emit` は 22.7 s。内訳: ① RADCOD の経路でも先に走らせていた z 閉じの割線法 (`_qcod_periodic` 付きで 3 回 ≈ 5 s、結果は捨てていた)、② 閉軌道の Newton 反復で毎回「D・B 付きの完全な一周」(1.9 s) を最大 30 回、③ 最後の一周。
- 1 周 (3563 要素、piece 5611 個) のうち**数値ヤコビアン `_num_jacobian` が 60%** (1.9 s 中 1.29 s)。写像評価 12 回/piece × 5611 ≈ 67000 回。残りは `tradke` 0.48 s、`np.eye`・行列積。
- 写像の種類別 (1 周の piece 時間): 四極本体 `f_body` 720 個 0.143 s、薄い多極 840 個 0.122 s、BEND 本体 520 個 0.108 s、BEND エッジ 1200 個 0.136 s、ドリフト 1020 個 0.054 s、四極フリンジ・半ドリフト 1320 個 0.109 s。
- cod だけ追う一周 (ヤコビアン・D・B なし) は 0.18 s。つまり「何周するか」と「1 周でヤコビアンを作るか」が効く。

### B-2. 変更 (結果は変えない: 1e-11 を守る)

1. **解析ヤコビアン** (`emit6d` の `st.piece(f, jac=…)`): 数値微分の置き換え。`_jac_drift` (厳密ドリフト)、`_jac_kick`/`_jac_thin` (tthine = ドリフト L/6・キック・2L/3・キック・L/6)、`_jac_tqente`、`_jac_tqfrie`、
   `_jac_tsolque_body`、`_jac_quad_full` (四極全体)、BEND は lattice の `_bend_edge_chain` (tbfrie.f の trans) と `_tbendebody0_mat` (tbende.f tbendebody0 の trans1) に `dv` 項 (`h0/h1³·L`) を足す。
   `_dvz`: `z −= dv(δ)L` の δ 微分は `+h0/h1³·L` で、`-dv'(δ)` と代数的に一致する (SAD の `h0/h1emit**3`、数値で 1e-12 まで確認)。
   z 行 (5,1..4) は従来どおり `piece` がシンプレクティック条件から組み直す (解析式で要るのは x,px,y,py の行と (5,6) だけ)。クリップがかかる等で解析式が使えない場合は None を返し数値微分に戻る。
   分割数 (nrad, ndiv, nq の `int()`) は写像側のままなので、SAD と同じ分割。`emit6d._USE_ANALYTIC_JAC = False` で全部数値微分に戻せる (テスト・デバッグ用)。
2. **cod だけ追う一周** (`tturne(…, want_jac=False)`): ヤコビアン・D・B・(RADCOD OFF では tradke) を省く。z 閉じの割線法 (`z_slip`) と、収束判定だけの反復に使う。`want_db=False` は SAD の irad=6 (D・B なし)。
3. **z 閉じの qcod 省略** (`_qcod_accepts_start`): 軌道 0 から cod だけ追った 1 周の残りが `qcod` の受理条件 (`r ≤ conv`、qtwiss.f) を満たすなら `_qcod_periodic` (ヤコビアン付きの反復、0.5 s) を呼ばない。
   誤差のない理想ラティスでは qcod も最初の反復で受理して初期値 0 を返すので同じ値。条件は線形の一周行列 (解析式、0.016 s) の Courant-Snyder パラメータで評価し、境界 (conv の 1/10 以上) は本物の qcod に任せる。
4. **RADCOD の経路の整理**: z 閉じを RADCOD の分岐の後ろへ (以前は先に走らせて捨てていた)。tcod の反復は、前回の残差が小さい (r<0.1) ときはまず cod だけの一周で判定し、収束していなければ同じ cod でヤコビアン付きにやり直す。
   さらに前回の残差が r<100 の反復は最初から D・B 付きの完全な一周にし、収束したらその状態を temit の最後の tturne として使う (やり直さない)。
5. 閉軌道の求め方を `closed_orbit6d` に独立 (`emittance6d` と `calc6d._closed_orbit` (= `Emittance[]` / CALC6D) が共用)。
6. `tecalc` の `tsolva` を numpy の SVD から SAD の `tsvd` (`pysad/tsvd.py`) に差し替え (上記 A2)。
7. 検討して**採用しなかった**: Newton 法のヤコビアンの使い回し (chord 法)。軌道が 1.7e-14 変わるだけで εx 等が 2e-8 変わる (量子励起の分割数の `int()` の境界。同様に dp の 1e-9 の変化で σxx が 1.6e-5 変わった)。`_REUSE_STEP = 0` で無効 (コードは残してある)。

### B-3. 時間 (ESRF v03_12、3563 要素、`emittance6d`、秒。測るたびに ±10% ばらつく)

| ケース | 変更前 (旧 emit6d) | A3 後 (tcod 直訳のみ) | B 後 |
|---|---|---|---|
| RFSW + RADCOD (`disp` の設定) | 22.7 (測り直し 14.1 は A3 後) | 14.1 | **1.6** |
| RFSW のみ (16b の設定) | 10.3 | 10.3 | **1.6〜1.9** |
| RF なし | 8.9 | 8.9 | **1.3〜1.8** |
| RFSW + RADCOD, z 閉じなし | 5.5 | 5.5 | 1.7 |
| `<work-esrf>/disp` スクリプト全体 (読み込み+CALC+EMIT+DISP B*) | 25 | 17 | **3.4** |
| `Emittance[]` (emit_extra 経由) RADCOD ON / OFF | 約 2 / 6〜7 | — | 2.4 / 2.3 (BeamMatrices 付きで +1〜3) |

RFSW+RADCOD の内訳 (B 後、ESRF): 解析ヤコビアン付きで 1 周 0.45 s、cod だけ 0.19 s、D・B 付き 0.9 s。tcod は RFSW OFF の本体 3 反復 (ヤコビアン付き 2 + cod だけ 1) + RFSW ON の 2 反復 (ヤコビアン付き 1 + 完全 1) = 約 1.5 s。

### B-4. 精度

- 非 RADCOD (旧 numeric vs 解析): εx 1.4e-11、εz 1.8e-12、σE 2.4e-13、σz 2.0e-12、α_c 4.9e-12、νx,νy 6.6e-11 (旧 pysad と新 pysad の差。SAD との差は変わらない)。
  **bench_level16b** (ESRF、RFSW NORAD NOFLUC、`$FORM="18.14"`) の全桁: 修正前と修正後で SAD との最大相対差が同じ (εx 1.29e-11、εz 1.1e-12、σE 2.9e-13、σz 2.7e-12、νx 6.7e-11 (旧 6.6e-11)、νy 3.3e-11、νs 7.4e-12、α_c 1.0e-11、bh 5.3e-12)。U0・Vc・trf0・減衰率は一致。
- 旧数値微分の雑音: K1=0 の BEND 本体のヤコビアンは刻み 1e-7 の中心差分で相対 1e-8〜1e-7 の丸め雑音 (`rhob`≈20 倍の桁落ち)。刻みを 1e-5 にすると解析式と 2e-10 で一致 (刻みを小さくすると悪くなる)。
  これが軌道のあるとき (RADCOD) に量子励起の分割 (`nq`, `ndiv` の `int()`) を跨ぐ形で効き、εx 5e-7 の雑音を作っていた。解析式で雑音が消え SAD と 2e-11 以内になった (§1 A3 の表)。
- **見つけたバグ**: `lattice._tbendebody0_mat` (4D の CALC は z 行を使わないので影響なし) の `xsn = sn − phi0` は符号が逆 (mathfun の `xsincos` は xs = x − sin x)。BEND 本体の (5,6) が `rhob·2(φ−sinφ)` ずれていた。直した。
- 解析ヤコビアンのユニットテスト (`pysad/tests/test_emit_jac.py`): 各写像を数値微分 (刻み 1e-6) と比べて相対 1e-8 以内、ESRF 1 セルの全 piece で 1e-7 以内、ESRF 1 セルのリング (`bench_item6_tiny.sad`、RFSW OFF/ON/RADCOD) の SAD の値と照合。

## 3. 項目 C: Emittance[] の使い方への影響

- esrf スクリプト (読み取り専用、`grep -rho "Emittance\[[^]]*\]"`) の `Emittance[...]` は 90 ファイルに計 150 回 (`Emittance[]` 138、`OneTurnInformation->True` 6 (全部コメント行)、`OneTurnTransferMatrix -> True` 6)。
  コメントでない (実行される) 行は 50 行で、形は 9 通り: `Print[Emittance[]]` (規則のリスト全体)、`a = Emittance[]; …(Emittances/.a)[1]…`、`H = HarmonicNumber/.Emittance[]`、`Print[j, BucketHeight/.Emittance[], BunchLength/.Emittance[]]`、
  `(Emittances/.Emittance[])[1]*1e9` (FitFunction の中)、`Print[(Emittances/.Emittance[])[1]*1000000000, (EnergyLossU0/.Emittance[])/1000000]` など。`EMITX/EMITY/EMITZ/SIGE/SIGZ` を `Emittance[]` の直後に読むスクリプトは無い
  (`touscheklife` は `emit;` の後、`touscheck` は `cell calc emit;` の後に `EMITX`/`SIGE` を読む)。
- 照合: `bench_item6_emitusage.sad` (10 通りの使い方 + `SaveEMIT->True` + `Length[Emittance[OneTurnInformation->True]]=30`、`Length[Emittance[OneTurnTransferMatrix->True]]=22`) を SAD と pysad で実行して比べた。
  `Print[Emittance[]]` の規則名・順序・値 (Tunes, U0, Vc, EquilibriumPosition, α_c, OrbitDilation, BucketHeight, HarmonicNumber=992, OrbitAtExit, DampingRate, Emittances, MomentumSpread, BunchLength, TuneShiftByRadiation …)、
  `a=Emittance[]`、`/.Emittance[]`、EMITX の保存 (既定は元に戻る、`SaveEMIT->True` で更新) が一致。εy は 0 同然 (SAD -4.7e-38、pysad -1.8e-30)。
- OrbitAtExit の z が SAD と 9 桁一致 (−1.14197668e-6)、`EquilibriumPosition` も一致 (以前の `Emittance[]` は閉軌道 (`calc6d._closed_orbit`) が旧手順で、trf0 が 20.619466 mm だった。A3 の共通化で解決)。
- 対応: 影響なし。ただし RFSW OFF の `Emittance[]`/EMIT は次の点が SAD と違っていたので直した (§4)。

## 4. 追加で直したもの

- **RFSW OFF でも SAD は RF の量を報告する**: tturneg/tecalc は RFSW に関係なく空洞の積算量 (vc0, dvcacc …) から Vc・wrfeff を作り、`omegaz=√|αp 2π heff Vc/pgev cos φrf| ω0/2π`、`sigz=|αp| σE c p0/h0/ωz`、`bh`、`trf0` (tturne の update) を出す。
  pysad は RFSW OFF で Vc=0, Equilibrium dz=0, Bucket height=0, εz=0, σz=0 を返していた。`tecalc` を直し (`bench_item6_tiny` の RFSW OFF: Vc 一致、bh 1.6e-12、εz 1.6e-12、σz 2.7e-12、trf0 一致)、
  RFSW OFF の trf0 は T(5,6) の符号と U0・Vc から (`emittance6d` は最後の一周から、`Emittance[]` は `closed_orbit6d(norf_trf0=True)` から)。`self.trf0` は SAD のグローバル trf0 と同じく以後の追跡に使われる。
- `calc6d._closed_orbit` (Emittance[] / CALC6D の閉軌道) を `emit6d.closed_orbit6d` に任せる (放射あり)。放射なし (`radiation=False`) は従来のまま。

## 5. 回帰

- pytest (`--ignore` 3 本のまま): **85 本 pass** (作業前 63 + 新規 `test_emit_jac.py` 22)。
- `run_regression.py` (QUICK でない全部、7・12t・13a・16a を含む 28 ベンチ): **FAIL 0**。KNOWN は 3C の光学 EX 2.8e-7 / EPX 7.8e-8 (既知の「10/13 からの結合時の分散」)と 13A-2 の TWB (既知) のみ。
  L12 (ESRF emit εx) 2.9e-15、L13B 5.7e-13、L14A 1.0e-10、L16A は DA 地図・Score・EMITS 全 PASS、L24A〜F は全 PASS。`QUICK=1` も FAIL 0。
- 作業前に保存した 13a/13b/14a/16a/18a/22a/23a/16b の出力との diff: 13a, 13b, 22a, 23a は**完全に同一**。16b・18a は雑音 (εy、TW の 1e-12) だけ。14a/16a は表示の変更 (f0・dz・Closed orbit の 2 行・Beam tilt/size の行) と εy だけ。
  16a の DA 地図は作業前 (旧 tar `pysad_20261016.tar.gz` の pysad でも) SAD と 1 文字違い (境界の 1 粒子、`1181` vs SAD `1171`、損失周の差 2 周で RES が FAIL) だったが、作業後は SAD と完全に一致 (PASS 6/6)。
  DA の初期粒子は EMIT の σE と正規モードで決まるので、境界の粒子が 1e-13 の差で変わる。tsvd への差し替え (εy・σ の雑音が SAD に近づいた) で SAD と同じ側に落ちたと思われる (因果の確認はしていない)。
- Level 17 (INTRA、`bench_level17_pysad.py`): 17a (RFSW) 15 ケース ≤6.4e-12 (作業前と同じ ≤7e-12)、17b (RADCOD ON) ≤1.7e-8 (旧 1e-7)。
- 作業中の全比較は `probes_item6/` (emit_modes_timing.py などの調査スクリプトと、`emit_base.pkl` = 作業前、`emit_new5.pkl` = tsvd 差し替え前の作業後の emittance6d 4 モードの結果)。

## 6. 未解決・気付いた点 (次の作業)

1. EMIT の報告の欠けている行: `Extended Twiss Parameters` のブロック、`Orbit dilation dl`、`Effective harmonic # h`、`Synchrotron frequency`、`Imag.tune` 行、
   `Damping time (sec)`、`Tune shift due to radiation`、`Damping partition number`、`Nominal spin tune`、`Polarization time`。値は `emit_extra`/`tecalc` に揃っているものが多い (TuneShiftByRadiation, HarmonicNumber, NominalSpinTune)。
   Beam tilt・εy は 0 同然の桁で SAD と違う (§1 A2)。
2. `DISP [loc [loc1]]` の位置番号の範囲 (`disp 94 98`) は pysad が数字をパターンとして扱い 1 行も出さない (SAD は 94〜98 の行)。オプション (GEO, ORBIT, …) も未対応 (既知)。
3. `calc6d.py`・`touschek.py`・`emit_extra.py` は薄い多極の分岐を各自のコピーで持ち、`st.piece(f_thin)` が数値微分のまま (emit6d.tturne の `jthin` を使っていない)。`Emittance[]` の 1 周あたり 0.1 s 程度の差で、共通化は未了。
4. `_qcod_periodic` (4D の閉軌道) は K1≠0 の BEND・四極で数値微分 (1 反復 0.5 s)。軌道のある場合の z 閉じ (誤差入り) は従来どおり遅い。
5. 解析ヤコビアンを入れていない piece: 四極の線形フリンジ (FRINGE=1,2,3 = tqlfre)、ずれ・ロールのある四極・薄い多極、MULT、ANGLE=0 の BEND (tsteee)、空洞、K1≠0 の BEND の `_tbendecorr` 以外 (これは元から解析行列)。それらは数値微分のまま (正しいが遅い)。
6. RFSW ON, RADCOD ON の tcod は、空洞がない・収束しないとき旧 Newton 法 (z と trf0 の分け方が SAD と違う) に落ちる。収束しない場合の SAD の挙動 (`closed orbit not found` の後も続行) は未検証。

## 7. 変更・追加したファイル

- `pysad/emit6d.py`: `_tcod_radcod`, `closed_orbit6d`, `emittance6d` の整理、解析ヤコビアン (`_jac_*`, `_dvz`, `_ROT_IN/_OUT`, `_USE_ANALYTIC_JAC`)、`_State.piece(jac=)`・`want_jac`・`want_db`、`tturne(want_db=, want_jac=)`、`_qcod_accepts_start`、`beam_tilt_sizes`、`tecalc` の RFSW OFF の RF 量。
- `pysad/ffscmd.py`: `disp_cmd` の名前・桁・ヘッダ・`$$$`、`_squeeze_space`、`print_emit` (f0, dz, Beam tilt/size)、`LINE["DIR"]`。
- `pysad/sadscript.py`: `_is_ffs_command` に 1 行 (`disp *;` を FFS コマンドとして扱う)。範囲外のファイルだが 1 か所のみ。
- `pysad/lattice.py`: `_tbendebody0_mat` の `xsn` の符号。
- `pysad/calc6d.py`: `_closed_orbit` が放射ありのとき `closed_orbit6d` を使う (範囲外のファイルだが 1 か所のみ)。
- `pysad/tests/test_emit_jac.py` (新規)、`bench_item6_*.sad` と `.result` (SAD の参照値: emit, emitusage, tiny, names)、`probes_item6/`。


## 7. 2026-10-16 追記: §6 の解決状況 (docs/history/work_emit_report.md)

- §6-1 EMIT 報告の欠けている行 → **解決** (Extended Twiss、dl、h、fs、Imag.tune、Damping time、Tune shift、分配数、Nominal spin tune、Polarization time。ESRF の RFSW × RADCOD 4 通り 148 行が SAD と一致。EPX の欄に EY が出る SAD の癖も再現)。偏極 (POL ON) の行は未実装。
  副産物: RFSW OFF + RADCOD ON の閉軌道 (放射込み) を SAD の tcod の rfsw=F 枝として実装、FSHIFT (dvfs) を emit の閉軌道に入れた。Beam tilt と εy の 0 同然の桁は SAD の雑音と一致しない (再現しない)。
- §6-3 薄い多極の各自のコピー → 既に `emit6d._tthine` を共有していた (記述が古い)。
- §6-5 解析ヤコビアンを入れていない piece → 四極の線形フリンジ (tqlfre)、ずれ・ロール付き四極・薄い多極は **解決**。MULT、ANGLE=0 の BEND (tsteee)、空洞、DROTATE 付き BEND は未了。
- §6-4 `_qcod_periodic` の数値微分 → 四極・薄い多極は解析化 (BEND K1≠0 は元から解析)。MULT・F1 付き BEND は数値微分。
- §5 の「KNOWN は 3C の光学 EX 2.8e-7 / EPX 7.8e-8」は誤り (再測定で 3.2e-11 / 6e-12。`bench_level3c_pysad.result` が古かった)。
- §6-2 (`DISP loc loc1` の範囲) と §6-6 (tcod が収束しないときの続行) は未了 (後者は NaN の報告の書式までは SAD に合わせた)。
