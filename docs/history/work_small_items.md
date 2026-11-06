# 小項目 A〜E の作業ログ (2026-10-16): RADINT・DISP の位置範囲・CALC6D 単独・FFS 内の ON/OFF・FSHIFT の追跡

担当: pysad/ffscmd.py, lattice.py, calc6d.py, track6.py, track6n.py, sadscript.py (新規 radint.py、周辺に track6n_ext.py / maintrack.py / dapert.py の 1〜2 行)。
docs/history/00_RESUME.md の未解決リストの 5 項目を、SAD (K64, WSL) の実測とソースの読みで 1 項目ずつ合わせた。

## 0. 結果の一覧

| 項目 | 結果 | ベンチ (SAD 参照 `.result` 付き) | 一致 |
|---|---|---|---|
| A. RADINT | **実装** (INTGRL.f の直訳 `pysad/radint.py`)。SAD の癖 (直前の rhob/rho0 の使い回し、単精度定数 ほか) も再現 | 29a (ESRF 6 ブロック 246 行)、29a2 (小リング 3 ブロック 120 行) | 29a 最大相対差 1.8e-7 (SAD の表示 7 桁)、29a2 3e-5 (V 面の小さい量だけ) |
| B. DISP の位置範囲 | **実装** (tfdisp.f の解析の直訳: 数値・名前・NAME.n±k・パターン・ALL・周回・R/O/E モード) | 29b (36 ブロック) | 行単位で全一致 (0 同然の雑音の桁を除く) |
| C. `CALC6D;` 単独 | **SAD は再計算しない** (想定していた挙動は無かった)。代わりに別の癖を見つけて合わせた: 4D の CALC は縦の列 (AZ..ZPY) を作り直さず直前の 6D の先頭行の値を全行に残す | 29c (14 ケース) | 全一致 (3.4e-12) |
| D. FFS 内の `ON/OFF flag;` | **ON/OFF の語は無効** (`Out[1]:= USE` を出すだけ) で、後ろの語は普通にフラグとして読まれる。`OFF INTRA;` は INTRA を**立てる**。pysad を合わせた | 29d (36 ブロック) | 全一致 |
| E. FSHIFT の追跡 | **実装** (dvfs を tconvm / tradk / tcav / tcavfrin の dv に足す。track6, track6n, track6n_ext, MAIN TRACK, DA) | 29e (NORAD / RAD / RAD+FLUC、乱数の最終 seed も) | 座標 ≤1.9e-11 |

回帰: §6。

## A. RADINT (INTGRL.f)

### A.1 SAD のソース (src/INTGRL.f 693 行、tffsa.f:617 `RADINT` → `tfsetparam; intgrl(lfno)`)

- 入口は直前の CALC の `twiss(j,0,*)` (DP0 の行、結合 R1..R4・閉軌道 DX..DPY を含む)。BEND/QUAD/SEXT を 2·IR 個のスライスに分け、スライスごとに FUNI(i,1..12) を作って SYMPS (複合 Simpson) で積分する。
  BEND: `IR=min(14,max(1,int(|ANG|/0.002)))`、QUAD: `int(|K1L|/0.02)` (1〜14)、SEXT: `IR=0` (点 1 つ)。ANGLE=0 の BEND と K1=0 の四極は読み飛ばし (長さだけ RLT に足す)。OCT/MULT などは入らない。
- スライスの行列: BEND は `tbendebody` の 1 次の行列 (K1 と分散ベクトル VM)、QUAD は `QUADRU` の閉じた式 (K1<0 は双曲線)、いずれも `TROT` でロール θ を掛ける。Twiss の伝搬は R 行列で結合を基準化する `TMATRS` (R0CAL/R0I、DETR≤1 と >1 の 2 枝) + `TWTRANS` (x/y のβ変換と分散・閉軌道)。
  エッジ: `E1/E2 ≠ 0` で `BEDGE` (tan 項) + `TEDGE` (I4 と I9, I10 の補正)。
- 被積分 (TCAL): I1=CKL0X·(R0I の η)、I2=CURV²DL、I3=|CURV|³DL、I4=(CKL0X+2CKL1X/CURV)CURV²(R0I(1,1)EU+R0I(1,2)EPU)+(CKL0Y−2CKL1Y/CURV)CURV²(R0I(3,1)EU+R0I(3,2)EPU)、I5=|CURV|³(GU EU²+2AU EPU EU+BU EPU²)DL、I6〜I8 は V 面、**I9=γx·DL、I10=γy·DL (色収差。K1 の項は入らない)**、I11=ΣCKL0X (全曲げ角)、I12=ΣCKL0Y。
- 出力: 見出し + `SYNCHROTRON RADIATION INTEGRALS` の ENERGY/DP/長さ/曲げ角/RF/I1X I1Y I2 I3 I4U I4V I5U I5V、`NU-X-`/`NU-Y-` (最終行の位相/2π)、運動量圧縮率 X/Y (I1/L)、ELOSS (kV)、DU/DV、減衰時間 (msec)、エネルギー幅、エミッタンス U/V、NU-S・バケット高・自然バンチ長・量子寿命 (ELOSS≤RFV のとき。そうでなければ `NOT ENOUGH RF VOLTS;  ELOSS,RFV = ` + リスト入力出力)、線形色収差 GX/GY (`-I9/4π`, `-I10/4π`)。書式は F9.3 / F8.5 / F12.8 / `1pg15.7` / `1pg12.4` / `1pD12.4` (gfortran と照合済みの `maintrack.fmt_g/fmt_f`)。

### A.2 SAD の癖 (docs/11 A34)

1. **tbendebody は bendib モジュールの rhob/rho0 を使うが INTGRL は設定しない** → 「直前に tbend 系 (tbend/tbende/tbendi/tsteee) が設定した値」= **ビーム順で最後に処理した曲率のある BEND の L/ANGLE** が**全 BEND**に使われる。
   ESRF ではリングの最後が反転した ARCAH2 なので最後の曲がりは BID5 (L/ANGLE = 20.757 m)。この値で I1X…I9 が SAD と 7 桁一致し、別の値 (BPI7 の 37.6 や BID1 の 77.5) だと I9 が 8e-4〜1e-3 ずれる (感度が高いので、この仮定は一致で裏づけられた)。
   EMIT の後 (tbende が rhob を設定) でも最後の曲がりは同じなので結果は変わらない。pysad は `radint.default_stale_rho(els)` で同じ規則 (最後の BEND の L/ANGLE)、`compute(..., stale_rho=)` で上書き可。
   (docs/16 §11 の「RADINT と Emittance[] の 13% 差」の候補の 1 つ。因果は未検証。)
2. 六極のブロックは DL を設定しない: TCAL は**直前の BEND/QUAD の最後のスライス長 DL** を使う (I9, I10 = γ·DL が六極でも足される)。DKL2 は代入されないので六極の CKL0X/CKL1X は常に 0。RLE=0 の六極は TCAL の**後で** `DL=0.2` (単精度の 0.2) にして ` L(SEXTUPOLE) EQ 0. SET L=0.2` を出す。DL, TETA, R0/R0I, FUNI などモジュール変数はプロセス中ずっと残る (2 回目の RADINT の DL の初期値は前回の最後の値)。pysad も `radint._STATE` に残す。
3. **単精度定数**: `COEFF=(…)*1e27` の `1e27` は単精度 (9.99999988e26) なので ELOSS が相対 1.15e-8 小さい。`FREQ*1.E-6`、`DL=0.2`、`EXP(2.3026*TAUQ)` も単精度。SAD の ELOSS (678.50075327596767) が 17 桁で一致して確認した。
4. エネルギーは `amass/1e9*H0` で amass は大域 `MASS` (既定 `elmass`=0.51099895e6、tfstk.f:232)。**pysad の MASS の既定値は 0.51099906e6 で SAD と違っていた** → 0.51099895e6 に修正 (ffscmd.py。影響は gamma を使う 1 箇所だけ)。
5. CAVI は VOLT/HARM/FREQ を**最後のもので上書き**する (合計しない)。HARM=0 なら FREQ·RLT/c、FREQ=0 なら HARM·c/RLT。pysad の parser は HARM から FREQ を補うので、SAD の定義値 (raw の FREQ が 0 なら 0) に戻して判定する。
6. CKL1X=DKL1·DX/DR (DR=|(DX,DY)|) は軌道の**向き**だけで決まる (大きさは効かない)。COD が 0 に近い四極は CURV が 1e-18 程度になり、I2〜I5 に寄与しない。

### A.3 pysad

- 新規 `pysad/radint.py` (INTGRL.f の直訳: `RadInt` クラスが TCAL/TMATRS/TWTRANS/TEDGE/BEDGE/BMAG/QUADRU/SEXTU/TROT/R0CAL を持ち、`compute`・`report` が本体と出力)。BMAG の `tbendebody` は `emit6d._BendIBody` (tbendiinit/tbendeinit の直訳) の 1 次の行列を使う。
  `FFS.radint_cmd` (ffscmd.py) が `RADINT` を呼ぶ (以前は `[pysad] FFS command RADINT は未対応`)。
- ゼロ除算 (L=0 の K1≠0 の四極) は SAD と同じく inf/NaN (numpy)、`  DL=0. IN TCAL ` を出す (SAD との照合は未了)。

### A.4 照合 (SAD の `.result` との行単位比較 `cmp_radint.py`、0 同然の V 面の量 <1e-10 は雑音として一致扱い)

- `bench_level29a.sad` (ESRF 連結): R1 CELL CALC 後、R2 EMIT 後、R3 NORFSW、R4 結合 (`LINE["DY","QDA2.1"]=2e-6`、`LINE["ROTATE","QFA4.1"]=5e-4`、`LINE["ROTATE","BID3H.1"]=5e-4`)、R5 `DP0=0.003`、R6 `DX QFA4 1e-4` (閉軌道≠0 の四極キック)。6 ブロック 246 行が最大相対差 1.8e-7 (SAD の表示が 7 桁) で一致。
  (R4 は最初 `DY 1e-4` を全 QDA2 に入れて R≈0.3 の強結合にしたら I4U が 5e-6 ずれた。ずれは BEND の I4/I9 の 1e-6 で、Twiss 入力 1e-9 の差が打ち消し合う項で拡大されたもの。要素ごとの DI を SAD の write 入りデバッグ版と比べて、四極・六極は ≤1e-8、BEND の大半も ≤1e-8 と確認し、ベンチは中程度の結合に変えた。)
- `bench_level29a2.sad` (13a の小リング: E1=E2=0.5·ANGLE の BEND、RF なし → `NOT ENOUGH RF VOLTS` の行が 17 桁・空白まで一致、RF あり、ロール付き四極 + DY): 3 ブロック 120 行が ≤3.2e-5 (V 面の EY ~1e-9 の絶対精度が I5V/I4V で拡大される量だけ)。
- 道具: `cmp_radint.py` (RADINT ブロックの比較)、`<work-scratch>/` (SAD のデバッグ版 (WSL のホームの作業フォルダ) の INTGRL.f に `DBGI` の write を入れる `patch_intgrl.py`、要素ごとの DI を比べる `t_dbg.py` `t_pat.py` `t_scan.py`)。

## B. DISP の位置範囲 (tfdisp.f)

### B.1 SAD の解析 (tfdisp.f:132-203、tffsa.f:595、ielm.f、temat.f)

- 語を 1 つずつ読む (`getwdl2`): (1) `| < { * %` のどれかを含む語 (`***` を除く) は**パターン** (seldis): 表示する行を要素名で絞る (`temat`)。(2) `ALL` は範囲を 1..nlat に。(3) モード語 (`A_CCELERATION G_EOMETRY GA_MMA OG_EOMETRY O_RBIT R_MATRIX B_EAM P_HYSICAL D_UMPOPTICS Z E_XTREMUM RE_FERENCE DRE_FERENCE`)。(4) それ以外は位置 (`ielmex`): **1 つ目の位置 = loc1**、続く**1 語が位置ならそれが loc2**、そうでなければ loc2=loc1 (その語は次のコマンドとして読まれる)。位置でない非空の語に出会ったら**何も表示せずに戻る** (その語が次のコマンド)。語がなければ前回の範囲 (id1, id2)。
- **範囲 (id1, id2) は FFS が保持** (tffsa.f:108 `id1=1; id2=nlat`、USE のたびに戻る)。`DISP 94 98` の後の `DISP;` や `DISP QF*;` は前回の範囲を使う。
- 範囲の走査: `loc1==loc2` は**その要素から $$$ まで** (刻み `nlat-loc2`、つまり loc1 と $$$ の 2 行)、`loc2<loc1` は nlat を超えて周回 (`l=mod(l1-1,nlat)+1`)、それ以外は loc1..loc2。
- 位置 (ielmex): 式を評価して実数ならその値 (`v>=0: int(v+0.499)`、負は `nlat+1+v+0.5`、1..nlat に収める)。`^` で始まる語 (`^^^`=1) は評価しない。名前は `NAME`・`NAME.n`・`NAME.n+k`・`$$$`/`***` (=nlat)・`^^^` (=1)。`NAME` は 1 台だけの要素かファミリーの**代表 (最初)**、`NAME.n` は同名が複数あるときの n 番目 (**1 台だけの要素の `RFC.1` は未定義**)。`+k/-k` の後は `ielma` で 1..nlat に収める (周回しない)。
  未定義の `NAME.n` は ` ???-FFS-Error-?Undefined location NAME.n` + `Out[1]:= USE` (以降のその文は捨てる)。`.` のない未定義の語は何も言わずコマンドとして読まれ ` ?Undefined command or element`。
- パターン (temat): 標準名 (`NAME` / `NAME.n`) に合う、または 1 台だけの要素なら `NAME.1` に合う、複数ある要素の代表ならファミリー名に合う。`2*47` のように `*` を含めば**パターン**になる (数式として評価されない)。
- 行の書式は既存どおり (`tdtrimz`: 値が厳密に 0 なら `.0` + 空白。以前は s(m) の列だけで、他の列は `.00000` だった)。ヘッダは表示する行が 1 行以上あるときだけ、66 行ごとに繰り返す (パターンで 1 行も一致しなければ何も出ない)。
- モード: R (R1〜R4 の列)、O (DX DPX DY DPY の列)、E (β の極値: 前後の β と比べて x か y が極大/極小の行だけ。パターン一致の行だけで状態が進む) を実装。**未対応** (`[pysad] DISP のモード … は未対応 (表示なし)`): A, G, GA, OG, B, P, D, Z, REFERENCE, DREFERENCE。

### B.2 pysad

`ffscmd.py`: `disp_cmd` を tfdisp の解析の直訳に書き直し (`_ielmex/_ielmf/_ielmh/_ielma/_temat/_eval_real`、行の出力は `_disp_show`)。`_disp_id1/_disp_id2` (USE で 1/nlat に戻す)。
`sadscript._is_ffs_command`: `disp ^^^ 5` の `^^^` をべき乗の演算子と見て SADScript の式にしていたのを、FFS 語の後に `^^^` が続くときは FFS コマンドにした。

### B.3 照合 `bench_level29b.sad` (`cmp_blocks.py`)

D1 `disp 94 98`、D2 名前 (`QFA1H.1 QDA2.1`)、D3 単独の数値、D4 `disp;` (前回の範囲の繰り返し)、D5 周回 (`3550 20`)、D6/D18 パターン (前回の範囲・指定範囲)、D7 `ALL QFA4.2`、D8 負の位置、D9 小数 (`94.4 98.6`)、D10/D31 `NAME.n±k`、D11 `$$$`、D12 `^^^ 5`、D13 `QD* 94 400`、D14〜D17、D20 `ALL` (3618 行)、D21 R、D22 O、D23 E、D24 未定義 (`ZZ.1`)、D25/D26 1 台だけの要素 (`RFC`, `RFC.1`)、D27/D28 範囲外 (3563 3570, 0 3)、D29 `QFA4.2 QFA4.3`、D30、D32 `disp 94 CALC` (DISP の後に CALC)、D33 `2*47 100` (パターン扱い)、D34 `{QB}FA1H.%`、D35 `QFA1H*|BID5*`、D36 USE のあと範囲が戻る。**36 ブロック全部が行単位で一致** (SAD の 1e-17 などの丸め雑音の桁と pysad の 0 は一致扱い)。

## C. `CALC6D;` 単独 (CALC なし) の直後

### C.1 SAD の実測 (`bench_level29c.sad`)

- `tflag` (tflag.f) は `calc6d` フラグを立てる/下ろすだけ。tffsa.f にも CALC6D/CALC4D の切替で再計算する処理はない。
- 実測: `CELL; CALC;` (4D) → `CALC6D;` → `Twiss[...]` は**直前の 4D の値のまま** (C2 = C1)。`CALC4D;` 単独、`Emittance[]`/`EMIT` を挟んでも同じ (C3, C5, C7, C11〜C14)。**SAD は再計算しない**。docs/history/work_emit_report.md の想定は 26d の最初の試行の読み違いと思われる (そこでは `CALC6D; CELL; CALC;` で、6D になるのは当然)。pysad も再計算しないので元から一致していた。
- 別の癖 (A36): **4D の CALC は縦方向の列 (AZ BZ NZ ZX ZPX ZY ZPY = 22〜28) を作り直さない**。tffscalc.f:87-90 は 4D (orbitcal=偽) では `ntfun=mfitdetr` (15) 列だけ先頭行に入れ直し、各行は qtwiss1 (qtwiss.f:536) で直前の行から縦の列を**写す**。結果、直前の 6D の CALC の**先頭行**の AZ/BZ/… が全行に残る (C6: AZ=−.0025607 が全位置で同じ、BZ=5.4317 が全位置で同じ)。
  DZ/DDP は 4D の閉軌道が書き直す (雑音の桁)。

### C.2 pysad

`ffscmd.py`: `FFS._z_stale` (USE で None)。6D の CALC の後に先頭行の `alfz betz muz zx zpx zy zpy` を保存し、4D の CALC (CELL も INS も) の rows の全行に入れる。

### C.3 照合 `bench_level29c.sad`: C1〜C14 (CALC6D/CALC4D の単独、`CALC6D; CALC;`、4D 後の Emittance[]/EMIT) が 3.4e-12 で一致。

## D. FFS 内の `ON flag;` / `OFF flag;`

### D.1 SAD の実測 (`bench_level29d.sad`、フラグは `FLAG[{"RAD","FLUC","INTRA","RFSW","COD","INTRES","CELL","CALC6D"}]` で確認)

| 文 | 結果 |
|---|---|
| `ON INTRA;` | ` Out[1]:= USE` と出る。INTRA は 0→**1** |
| `OFF RAD;` | `Out[1]:= USE`。RAD は 1 のまま (切れない) |
| `OFF INTRA;` (INTRA=0 のとき) | `Out[1]:= USE`。INTRA は 0→**1** (OFF なのに立つ) |
| `OFF NOINTRA;` | `Out[1]:= USE`。INTRA=0 (NOINTRA が読まれる) |
| `ON INTRA CALC6D;` | `Out[1]:= USE` が **1 回**。INTRA, CALC6D とも立つ |
| `ON NORAD NOFLUC;` | `Out[1]:= USE`。RAD, FLUC とも切れる |
| `ON;` | `Out[1]:= USE` |

→ **FFS 内の ON/OFF の語は無効** (`Out[1]:= USE` を出すだけ) で、後ろの語は普通に次の文として読まれる (フラグ名なら立てる、`NO`+名なら下ろす)。つまり `ON INTRA;` ≡ `INTRA;`、`OFF INTRA;` も ≡ `INTRA;` (立てる)。INTRA を切るには `NOINTRA;`。
MAIN 部 (FFS; の前) の `ON/OFF` は別で有効 (`parser.parse_main_commands` → `FFS.main_exec`)。ユーザースクリプト 57 本の FFS 内 `ON INTRA` は、メッセージが出て INTRA が立つ点で期待どおりに動く (OFF を FFS 内に書いたものは INTRA が立ってしまう)。

### D.2 pysad

`ffscmd._exec_word`: ON/OFF は `Out[1]:= USE` を書いて 1 語だけ消費 (旧: 後ろのフラグをまとめて ON/OFF していた)。`main_exec` に MAIN 用の ON/OFF 処理を移した。

### D.3 照合 `bench_level29d.sad` (E0〜E16、36 ブロック、フラグの値と `Out[1]` の行) 全一致。

## E. FSHIFT の追跡

### E.1 SAD で dvfs が入る経路 (grep dvfs src/*.f)

`tsetdvfs` (tturne.f:859) = `merge(rgetgl1('FSHIFT'), 0, rfsw)` が tftrack.f:219 / tracka.f:45 / trackd.f:86 / track.f:111 で呼ばれる。dv に足される場所 (追跡に関わるもの):
- 入力の変換: `tconvm.f:18` / `tconv.f:28` / `trackd.f tinip:521` (`dv = −g(1+pr)/h1/(h1+pr h0) + dvfs`)
- 放射で g が変わった後: `photons.f` tradkn (1064, 1146, 1202, 1256)、`trad.f`、`tsdrad.f`、`tserad.f`、`tqrad.f`
- 空洞: `tcav.f` の tcav (140: `h1=p1r h0/(1−dv+dvfs)`、158: `dv=…+dvfs`)、tcavfrin (230, 246)、`tmulti.f` (571, 585, 672)
- ドリフトなどは保持している dv を使うだけ (`z += dpz·al1 − dv·al`)。tcav の `dvn` (dzn の計算) は dvfs を含まない。

### E.2 pysad

`emit6d._Phys.dvfs` (既存) を追跡にも使う: `track6.track_particles(…, fshift=)` / `track6n.track_particles_fast(…, fshift=)` / `track6n_ext.track_particles_ext(…, fshift=)` が `ph.dvfs = fshift if rfsw else 0` を設定。
- track6.py (純 Python): 入力の dv は `ph.dv` が dvfs を含むのでそのまま。`_tradk` の dv、`_tcav`/`_tcavfrin` の h1 と dv に `+dvfs`。
- track6n.py (numba): 粒子状態 S を (n, 10) → **(n, 11)** にして最後の列に dvfs を持たせた (カーネルの引数を増やさずに済む。`make_state` が `ph.dvfs` を入れる)。`tradk_core1`、`tcavfrin1`、`k_cav` の dv に `+S[i,10]`。
- 呼び出し側: `FFS.fn_track` (TrackParticles)、`FFS.fn_dasurvey` → `dapert.dynamic_aperture(fshift=)`、`maintrack.run_track` (MAIN の TRACK) が `FSHIFT` (RFSW ON のとき) を渡す。

### E.3 照合 `bench_level29e.sad` (ESRF、`ON RAD FLUC COD RFSW RADCOD`、`FSHIFT=1e-6`、5 粒子 (閉軌道のまわり、δ≠0・z≠0 を含む))

N1/N10 (NORAD 1 周・10 周)、R1/R10 (RAD NOFLUC)、F10 (RAD+FLUC、`SeedRandom[17]`、乱数の最終 seed も) の 40 行が**座標 ≤1.9e-11** で一致 (dvfs を 0 にすると 10 周後に x が 6e-6、z が 5.9e-3 ずれるので、影響の大きい項目)。
- 課題の `FSHIFT=1e-4` は **SAD 自身が成り立たない** (周長 603 m で dl=−60 mm、閉軌道の δ=−0.094、`???-Emittance[]-closed orbit not found.`、粒子はすべて落ちる)。1e-6 (δ=−4.3e-3) を使った。
- 別件 (未解決): `RFSW ON + RADCOD OFF` で FSHIFT≠0 の Emittance[] の閉軌道は SAD の tcod と一致しない (x が 5e-5 相対、px 6e-9 の閉じ残し。emit6d の z を閉じる δ の割線法 `close_z` が FSHIFT≠0 で精度不足)。RADCOD ON (tcod の直訳) は 26c で ≤1.5e-10 一致なのでベンチは RADCOD ON にした。

## 6. 回帰

- 作業前に 13a/13b/14a/16a/22a/23a/28a の出力を `base29/` に保存して、作業後と diff した: **13a・13b・22a は同一**。14a・16a・23a・28a は **「条件表」の改行だけ**が違う。
  tfshow.f:70,418-433 は行の長さを `lw=max(79,min(256,PageWidth))` で決め、セルを ' ' でつないで超えるなら改行する。pysad は常に 3 セルで改行していたので、`PageWidth=1999` のスクリプトで SAD (6 セルが 1 行) と違っていた。これを SAD の規則にした (`print_fit_table`)。SAD の `bench_level14a.result` などと一致する側に直った。
- `QUICK=1` なしの `run_regression.py` 全 43 ベンチ (29a, 29a2, 29b, 29c, 29d, 29e を追加): **FAIL 0** (13a の TWB が既知 KNOWN のまま)。29a/29a2 は RADINT ブロック (SAD の表示 7 桁で一致)、29b/29d は全ブロック一致、29c は 3.4e-12、29e は 1.9e-11。
- pytest (`--ignore` test_bend_feeddown / test_coord / test_sext_radiation): **159 passed** (新規 `test_small_items.py` 6 本: RADINT の小リング照合と単精度定数、DISP の位置記法と範囲の保持、DISP の 29b 照合、FFS 内の ON/OFF、追跡の FSHIFT (track6 と track6n が一致し、FSHIFT で値が変わり、RFSW OFF では効かない))。
- 配布: `pysad_20261016.tar.gz` を更新 (radint.py, 29a〜e ベンチ, cmp_radint.py, cmp_blocks.py, docs/ 含む)。`python make_html.py` で docs/html を再生成 (43 ファイル)。

## 7. 未解決・次の候補

1. RADINT: ゼロ長の K1≠0 四極 (DL=0 の NaN 表示) の SAD との照合。`stale_rho` の規則が「最後の曲がりのある BEND の L/ANGLE」でなくなる場合 (INS で一部だけ CALC するなど) は未確認。
2. DISP: モード A/G/GA/OG/B/P/D/Z と REFERENCE/DREFERENCE (SAD の列の定義は tfdisp.f にある。tdrefheader/trefheader の見出しの書き換えも要る)。
3. FFS の `RFSW ON + RADCOD OFF + FSHIFT≠0` の Emittance[] の閉軌道 (上記 E.3)。
4. `CALC6D;` 単独に関する想定は消えた (SAD は再計算しない)。
