# 03. FFS のコマンド解釈 (src/tffsa.f) と pysad の ffscmd.py

最終更新: 2026-11-06 (現状のスナップショット)

対象: SAD K64 (<SAD>\src) と pysad (<repo>\pysad)。
すべて Fortran ソースと pysad のコードで確認した内容。行番号は 2026-10 時点のもの。

## 1. 1 文の解釈順序 (tffsa.f: subroutine tffsa)

FFS は `getwrd` で語 (word) を 1 つ取り、次の順で「誰かが受け取るまで」試す。受け取った (exist) ら次の語へ (`go to 10`)。

| 順 | 呼び出し | 受け取る語 | 備考 |
|---|---|---|---|
| 0 | `tfprint` (tfprint.f) | SADScript 式 | `itfpeeko` で先読みし、定義済みシンボルや `name[`、`name=` など「式」なら評価して値を表示。`;` で終わる式は表示しない。FFS コマンド名と一致しない語はここで式として評価される |
| 1 | `UNTIL` / `REP_EAT` | ループ | `nrpt(levelr)`、`irptp` で入力位置を巻き戻す (maxrpt=32) |
| 2 | `tfif` | `IF/ELSE/ENDIF` | |
| 3 | `tflag` (tflag.f) | フラグ | §2 |
| 4 | `tgetfv` (tgetfv.f) | FIT 条件 (`NX 7.5 2`, `BXM 10`, `BXI`, `BXSCALE`) | §3 |
| 5 | `terror` (terror.f) | `DELX/DELY/DTHETA/DL/DK/DDK/DBEAM/DUMP` | §4 |
| 6 | `tffile` | `INC_LUDE`, `OUT`, `LOG` などファイル系 | |
| 7 | `tfgetlineps(word,…,mode=1)` | 要素名・パターン (`QF`, `Q*`, `QF.2`) | 一致したら label 7000 の `tfgetv` (§4) |
| 8 | コマンド群 | `STOP/EXIT/QUIT/ABORT/USE/VISIT/BYE/SPLIT/MAXI_TERATION/ATT_RIBUTE/SAVE/EXPAND/RESET/FREE/FIX/FIT/FITP_OINTS/CAL_CULATE/GO/REC_OVER/T_YPE/DISP_LAY/SCALE/AP_ERTURE/MAT_RIX/CHRO_MATICITY/RADINT/REF_ERENCE/DIMAD/ZAP/TRC_OD/PUSHS_EED…/DELW_AVE` | 第 1 群 |
| 9 | | `VARY/REJ_ECT/COUP_LE/DECOUP_LE/INDEP_ENDENT/ORI_GIN(ORG)/SHOW/STAT_US/VAR_IABLES(VARS)/RENUM_BER/REV_ERSE/DRAW/GEO/PRI_NT` | 第 2 群 |
| 10 | | `MEA_SURE/ORBITJ_ITTER/PUTM_EASURE/GET_MEASURE/LOGM_EASURE` | 第 3 群 |
| 11 | | `WAKE/DELC_OR/LTR_ACK/EMIT_TANCE/SYNCHROB_ETA/BEAM_SIZE/ALI_GN` | 第 4 群 |
| 12 | label 7000 `tfgetv` | 要素値 | 失敗なら `?Undefined command or element:` |

**abbrev の記法** (abbrev.f): `abbrev(word,'CAL_CULATE','_')` は `_` より前 (`CAL`) が必須部分、以後は省略可能な接尾辞。`CAL`, `CALC`, `CALCU`, … `CALCULATE` のどれでも真、`CALCULATED` のように長い語や `CAX` は偽。`_` のない参照は完全一致。pysad では `_abbrev(word, 'CAL_CULATE')` が同じ規則。

**REPEAT n … UNTIL cond** (tffsa.f:255-290): `levelr`, `nrpt`, `irptp` で読み位置を巻き戻す。n 省略は 65535、UNTIL に条件が無ければ回数だけ、条件が真なら抜ける。多重可。pysad は `Interp.run` が最上位 FFS の文列を添字つきのループで回し、`_control_word`/`_control` が処理する (`FFS["…"]` の文字列の中では働かない)。

**IF / ELSEIF / ELSE / ENDIF** (tfif.f): 偽なら同じ深さの ELSE/ELSEIF まで読み飛ばし、ELSE に来たら ENDIF まで読み飛ばす。語の直後が `[ * / ) } ] > < = . & ^ |` なら制御語ではない。`IF 条件 コマンド` の同じ行の続きも実行する。pysad は `_if_cond`/`_skip_if`。条件式は空白で終わる (演算子の前後の空白は式の続き)。エラーは ` ???-FFS-Error-?ELSE without IF ` など SAD と同じ書式。

`CAL`/`GO` は label 1000 に飛び、`tfsetparam` → (CELL なら) `tfsetupcell` (区間 1..nlat の AX, AY = 0 の条件を追加、mfitp は既存条件の最大 |mfitp|) → 有効条件の抽出 (`kfit<=mfit` かつ `mfitp/=0`) → `tffssetupcouple` → `nfr = max((|mfitp|-1)/2)` と dp 格子 (`DP`, `DPM`) → `tffsmatch` → `tfshow` → CELL で不安定なら `Horizontal/Vertical Unstable.` を出す。`busy` フラグで再帰を禁止 (`?Recursive CAL or GO`)。

## 2. フラグ (tffs.f fname/sino, nflag=53)

`tflag` は `NAME` → 真、`NONAME` → 偽、`SINO` (反対語) → 偽、`NOSINO` → 真。既定は initb1.f の `defflg` と tffsa.f の初期化 (cellstab, simulate, absweit, jitter, bipol, fseed = T; cell, gauss, ideal, codplt = F)。pysad は `_FLAGS_DEFAULT` と `_SINO`/`_FLAG_ALT` に同じ表を持つ。

| fname | sino | 既定 | 意味 |
|---|---|---|---|
| TRPT | RING | F (RING) | 輸送路。Emittance[] の平衡計算を省く |
| CELL | INS | F (INS) | 周期解 (閉軌道・周期 Twiss) か、先頭 MARK の Twiss から伝搬か |
| RFSW | | F (FFS 進入時に track.f:148 で ON に強制。TRPT=OFF, CALC6D=OFF も同時) | 空洞を ON (6D) |
| RAD | | F | 放射減衰 (追跡・EMIT) |
| FLUC | DAMPONLY | T | 量子励起 (乱数) |
| RADCOD | | F | 閉軌道に放射損失を入れる |
| RADTAPER | | F | 局所運動量に合わせて磁石をスケール |
| CALC6D | CALC4D | F | 6D の CAL/GO |
| KEEPEXP / CALEXP | | | 要素値展開 (EXPAND/NOEXPAND) の既定 |
| EMIOUT | | F | EMIT の詳細出力 |
| CODPLOT | | F | EMIT で 6D ビームサイズも |
| COD | | T | 閉軌道計算 |
| INTRA | | F | 粒子内散乱 |
| POL / RADPOL | | F / T | 偏極 |
| LWAKE / TWAKE / WSPAC / SPAC / PSPAC | | F | ウェーク・空間電荷 |
| SELFCOD | | F | |
| CONV | | | GO の収束 (内部) |
| STABLE | UNSTABLE | T | 不安定でも条件に入れるか (cellstab) |
| GEOCAL | GEOFIX | T | 幾何の再計算 |
| RADLIGHT / PHOTONS / LOSSMAP / SORG | | F | |
| INTRES / HALFRES / SUMRES / DIFFRES | | T / T / T / F | GO で整数・半整数・和・差共鳴をまたぐ解を不安定扱い |
| FFSPRMPT / CONVCASE / PRSVCASE / SUS | | | 入力の大文字化など |
| K64 | LEGACY | T | K64 仕様 (TYPE で DETR/BZ/NZ/DZ を出す等) |
| GAUSS | UNIFORM | F | 追跡の初期分布 |
| FIXSEED | MOVESEED | T | 乱数種の固定 |
| BIPOL | UNIPOL | T | |
| ORBITCAL / CALOPT (ORBONLY) | | T | |
| DAPERT | | F | |
| IDEAL | REAL | F | 誤差なしの理想値で計算 |
| FOURIER / TRACKSIZ | | F | |
| SIMULATE | OPERATE | T | |
| ABSW | RELW | T | FIT 重みを絶対/相対 |
| JITTER | QUIET | T | |
| TRGAUSS | TRUNI | | 追跡の運動量分布 |
| BARYCOD / WAKEOPT | | F | |

**FFS 内の `ON flag;` / `OFF flag;` は無効** (SAD K64、docs/11 A35): `ON`/`OFF` の語は `Out[1]:= USE` を出すだけで何もせず、後ろの語が普通に次の文として読まれる。つまり `ON INTRA;` ≡ `INTRA;` (立てる)、`OFF INTRA;` も INTRA を**立てる** (切るのは `NOINTRA;`)、`ON NORAD NOFLUC;` は RAD と FLUC を切る。MAIN 部 (FFS; の前) の ON/OFF は有効。(pysad: `_exec_word`、docs/history/work_small_items.md D)

**EMIT コマンドは GAUSS を ON にする** (temitf の終了後 `gauss=.true.`、docs/11 A39)。pysad も同じ。`?FLAG` でフラグの値を読める。

`STAT_US` で全フラグ (偽は sino か `NO`+名) と CHARGE/MASS/MOMENTUM/DP/DP0/EMITX…/NPARA、`twelm` による `FIT`, `DISP_LAY`, `ORG`, `MEA_SURE` の位置が表示される。

## 3. FIT 条件 (tgetfv.f)、FREE/FIX (tffsfreefix.f)、GO の解法とフック

### 3.1 位置: `FIT [loc [loc1]]` (tffsa.f)

- `FIT ALL` → `mfpnt=1, mfpnt1=nlat`。
- `FIT loc` → `mfpnt=mfpnt1=ielme(loc)`。`FIT loc loc1` → 区間 (ifitp=min, ifitp1=max は tgetfv で正規化)。
- 位置が解決できない語 (`FIT NX 7.5` の `NX`、`FIT;` の空語) は `ielmf` が `nlat` ($$$) を返すので **`FIT;` は位置を末尾 $$$ に戻すだけで条件は消さない** (消すのは `REJECT TOTAL`)。位置は次の `FIT` まで保持され、以後の `NX …` はその位置に付く。
- `ielmf` (ielm.f): `$$$`/`***` → nlat、`^^^` → 1、`NAME.n` (n は式も可)、`NAME+k`/`NAME-k` (k は式、小数部は要素内の割合)。`fit QF` のようにファミリー名だけなら `ielmh(name,0)` = 単独要素か多重要素の **先頭 (代表)**。

### 3.2 関数: `FUN value [npoints]`

- `FUN` は nlist (tffs.f:1622) の 64 個: `AX BX NX AY BY NY EX EPX EY EPY R1 R2 R3 R4 DETR DX DPX DY DPY DZ DDP AZ BZ NZ ZX ZPX ZY ZPY PEX PEPX PEY PEPY PZX PZPX PZY PZPY GMX GMY GMZ BMAGX BMAGY BMAGZ TRX TRY TRZ LENG GX GY GZ CHI1 CHI2 CHI3 DEX DEPX DEY DEPY DDX DDPX DDY DDPY PDEX PDEPX PDEY PDEPY`。`mfit` までが FIT 可能 (TRX/TRY は末尾のみ、区間 FIT は mfit 以下で TRX/TRY 以外)。GX/GY/GZ/CHI1-3 の fit は、pysad では値を geo.py の幾何から取り、微分は数値微分、CHI の scale は π/180 (docs/history/work_geo_sig.md)。
- value: 実数、`{min,max}` (不等式。両端が同じなら実数扱い)、`{v}`、`*` (登録済みの値を保つ: 無ければ `No default value of`)、`@`/`@-`/`@*` (位置の MARK 要素の保存値、`@-` は符号反転、`@*` は進行方向 direlc で AX/EPX/R2 等を反転)。`:` を前置すると scale を掛けない (NX/NY/NZ は 2π、CHI1-3 は π/180 が scale)。
- npoints: `mfitp = n+1`。省略 (または 0) は既存値を保ち未設定なら 2。`-1` で無効 (mfitp=0)。`FUNM` (末尾 M) は最大/最小フィット (mfitp を負)。BX/BY に 0 以下を与えると `Non-positive value for`。
- `FUNI v` (`BXI 1.0`) は先頭要素 (MARK) の初期 Twiss `rlist(latt(1)+i)` を設定。`FUNSCALE v` は scale(i) を変更。
- `FITP_OINTS n` (tffsa.f): **|mfitp|>2 の条件だけ** `sign(|n|+1, mfitp)` に一括変更。オフモメンタム点数 `nfr=max((|mfitp|-1)/2)` が 1 以上なら `DP` (dpmax) が必要 (`?No momentum bandwidth.`)。
- FITP>1 のとき `MatchingAmplitude = {{dp, xamp, yamp}, …}` と `NFAMP` で有限振幅ファミリー (tffamsetup) を足す: 各 dp 点の軌道 + 正規モード振幅の軌道のまわりの周期解 (閉軌道は探さない) で条件を評価する。FitValue の {id, dp} の id は kfam (x: 1..NFAMP、y: −1..−NFAMP)。pysad に実装はあるが SAD との照合はない (自己検証のみ。§9 のとおり FITP>1 の参照値が取れない)。
- `REJ_ECT` (tfrej.f): `TOTAL` は全条件 mfitp=0 + `tfinitcalc`、`TOTALFIT` は mfitp=0 のみ、関数名 (パターン可、`…M` も) は現在位置 (mfpnt,mfpnt1) の該当条件だけ mfitp=0。
- `CAL_CULATE [FUN|FUN-]…` (tfcalc.f): 表示に関数を追加/削除 (`txcalc`) してから計算。

### 3.3 FREE / FIX (tffsfreefix.f)

`FREE pat [KEY] pat …`: `tfgetlineps(mode 0)` で要素位置を集め、`.` を含まない名前は **ファミリー先頭 (klp==k) だけ** を変数にする。`KEY` (次の語) が要素型のキーワード (`tfkwrd`) なら ivvar にそれを、無ければ既定キーワード `nelvx%ival` (tfinit.f: DRIFT は L (ival=1, 下限 vlim(1)=0.1)、BEND/QUAD/SEXT/OCT/DECA/DODECA/SOL/CAVI/TCAV は 2 番目のキーワード (ANGLE/K1/K2/K3/…/BZ/VOLT)、MULT は K1、MARK/MONI は 0 → `Can't use as variable`)。`BXI` 等の `nlistI` は先頭 MARK の初期値を変数にする (it=41)。`QF.2` のように成分を指定すると `ivcomp` 付き (`Info-Component … is coupled to`)。`valvar = 値/errk` で保持。`FIX pat [KEY]` は一致する変数を削除。ワイルドカード (`%*<{|`) があれば複数語を続けて読む。

### 3.4 GO の解法・微分・フック (tffsmatch.f)

- **tsvd (tfsolv → tsolva)**: 2 重対角化の後、対角の最小値 > `anorm*eps` なら特異値を落とさず直接解く (quick path)。そうでなければ暗黙 QR + 滑らかな減衰。pysad は `pysad/tsvd.py` (直訳)。J の条件数が 1e6 以上のとき (INS のフィット) numpy の硬い打ち切りと解が大きく違う。`PYSAD_TSVD_NUMPY=1` で numpy 版に切り替え。
- **INS のフィットの微分**: SAD は INS でも解析微分 (tffsqu→qdcell→qdtwis、cell=False)。`FFS$NumericalDerivative` のときだけ数値微分。pysad も match6d の解析微分を INS で使う (`PYSAD_INS_NUMERIC=1` で数値微分)。未対応の条件・状況は中心差分に自動フォールバック (`Matcher._analytic_unsupported`)。
- **FitValue / FitWeight のフック**: `FitValue[comp, fn, {id_,dp_}, goal_, now_]` (2 点条件は `FitValue[c1,c2,fn,{id,dp},goal,now1,now2]`)、`FitWeight[comp, fn, {id_,dp_}, w_]`。pysad は {id, dp} を `{0, dp}` のリストで渡す (裸の実数では `{_,dp_}` に一致しない)。
- **VariableRange / VariableWeight** (tffsvarfun、tffsmatch.f:717): VariableRange が実数を返せば dlim、0 なら変更前の値、`{min,max}` で範囲を狭める。VariableWeight は変数の重み gw を置換。pysad は `Matcher.var_range/var_weight` (`_vlimit`/`_tweigh`)。81 本のユーザースクリプトが使う。
- **OpticsProlog / OpticsEpilog**: `OpticsProlog := …` / `OpticsEpilog := …` を定義すると CALC (と fit の各評価) の前後で呼ばれる。`=.` で外す。

## 4. 要素値コマンド (tfgetv.f) と誤差コマンド (terror.f)

### 4.1 `NAME [KEY|R|MAX|MIN|MAXABS] value …`

`tfgetlineps(mode 1)` でファミリーを列挙し、次の語が `R_ELATIVE` (相対: `v*(1+x)`)、`MAX`/`MIN`/`MAXABS|ABSMAX|MAXMIN|MINMAX` (変数の上下限 `nelvx%vlim` を設定し、現在値をその範囲に丸める)、要素型のキーワード (`itftypekey`)、または値。キーワードなしは既定キーワード (`?No default keyword for` は MARK など)。パターンで型が混ざると `Info-Different types of elements match`。既定キーワードの代入は `cmpd%value(iv)=vx*errk(1,ii)` (誤差係数を掛けた実効値)。それ以外は `tftouch`。`S* 0`, `SDA MAX 0`, `QF K1 0.5`, `QF R 0.01` が典型。`QF.3 K1 v` (`.` 入り) は成分だけ。

`LINE[key,name]=` で触ったキーは USE まで残り、CALC のたびに代表の値に揃える (tftouch)。

### 4.2 誤差: `DELX/DELY/DTHETA/DL/DK/DDK/DBEAM sigma [opts] pat …`、`DUMP [pat]`

- `sigma` の前後に `U_NIFORM`/`R_ANDOM` (uni)、`C_OHERENT`/`INC_OHERENT` (cohe)、`A_DD`/`P_UT` (add) を置ける。既定は RANDOM・INCOHERENT・PUT。
- 各要素 (パターン一致、DRFT 以外は `master(i)` の先頭から成分列) に `errg = sigma*tgauss()` (uni なら sigma、coherent なら 1 回だけ引いた値を共用) を加える。**乱数は要素ごとに 1 回**引く (tgauss は GCUT で切る正規乱数、追跡と同じ乱数列)。
- master グループ (tmast.f): 同じ要素の分割片 (間にゼロ長要素のみ) はグループの先頭ごとに乱数を 1 回引き、同じ誤差を共有する。pysad は `_master_array`。
- 代入先: DELX→DX、DELY→DY、DTHETA→DROT (無ければ ROT)、DL→L (`max(0, …)`)、DDK→K0、DK→`errk(1,j)=1+errg` (ADD なら `errk+errg`) で、既定キーワードの実効値 `= 設計値/er*errk*couple(j)`。PUT の基準は設計値 `idval` (SAVE で更新)、ADD は現在値。DK は主強さだけに掛ける。`LINE["K1"/"ANGLE", …]` は DK 込みの値を返し、`Element[]` は設計値を返す。
- `DBEAM` は先頭の Twiss を揺らす。`DUMP pat` は `element DL(mm) DK DDK DTHETA(mrad) DELX(mm) DELY(mm)` を `(1x,a,3p,f10.5,1p,2g12.5,3p,3f10.5)` で、ずれのある要素だけ表示。
- uni でなければ最後に `tfadjst` (tfadjst.f): 隙間なく並ぶ磁石の鎖は先頭磁石の DX/DY のずれを共有する。pysad は `_tfadjst`。

詳細は docs/history/work_codcorr.md。

## 5. USE / CELL / CALC / GO / EMIT / TYPE / DISP / SAVE / RESET / COUPLE / CHRO

- **USE line** (`VISIT` は保存付き、`BYE` で戻る): `[NO]EXP_AND` → `tfbeamline` でビームライン決定 → `expnln` で展開 → `nlat=nlat0+1` → `chguse=.true.` で label 101 に戻り再初期化。LINE 名は書かれたとおりに引く (大文字小文字混じりの `LinacPFBT` も可。pysad `Env.resolve_beamline`)。`USE BeamLine[ORG,..]` / `VISIT BeamLine[..]` は新しい LINE `L000001$` を作る (先頭は MARK 必須)。VISIT/BYE は直前のラインに戻る (pysad `_use_beamline_expr`, `_visit_stack`)。
  最初の USE では `{EMITX,EMITY,EMITZ,SIGZ,SIGE}=LINE[{"EMITX","EMITY","EMITZ","SIGMAZ","SIGE"},1]; DP0=LINE["DDP",1]` で **先頭 MARK の値をグローバルに取り込む** (tffsa.f:135 `if(.not. chguse)`: FFS に入った最初の 1 回だけ)。以後の USE では直前の値 (emit 失敗時は NaN も) が残る。scale(NX,NY,NZ)=2π、ExponentOfResidual=2、NP=1000 も設定 (tffsa.f:141-177)。CONVERGENCE=1e-9 は FFS 進入時に track.f:35 の `tfevals` で設定される。
- **CELL/INS** はフラグ。CELL では周期解を求め、INS では MARK の初期値 (twiss(1,0,:)) から伝搬する。
- **CALC の閉軌道 (qcod)**: 一周行列の跡が NaN なら over で打ち切り、twiss は前回の値のまま (pysad `_qcod_periodic`、`_last_rows`)。
- **CALC/GO 後の表示** (tffsmatch.f:170-191): 収束で ` Matched. (res)` 、未収束で ` Residual =res`、不安定で `Unstable = n Residual =res`。書式 `(a,1pG11.4,a,' DP =',0pf8.5,'  DP0 =',f8.5,'  ExponentOfResidual =',f4.1,a,' OffMomentumWeight =',f8.3,a)`、末尾は `  CALEXP`/`  NOCALEXP`。続いて `tfshow` (tfshow.f) が条件表: nfam>0 なら ` DP `/`Res.` 行、各条件は `名前(11 桁, FIT 位置なら 11 桁目 'f')`+`関数(5)`+`目標値 (voutfv, S 書式)`+`点数 (|mfitp|-1, I3)`+`現在値 (autofg, 幅 lf=6..10)`。未登録の表示関数は ` #######  #`。区間条件は `名/名1` と差 (BX 等は比、BMAG は tgetbmagu2)。
- **EMIT_TANCE [pspan]** → `temitf` (EMIT は自動的に RING)。閉軌道、P0/f0、U0/Vc、dz/alpha、Bucket height、Real tune、Damping、Emittance X/Y/Z、Energy spread、Bunch Length、**Beam tilt、Beam size xi/eta** を表示し EMITX… を更新。pysad `print_emit` はこの書式を再現 (Extended Twiss・Orbit dilation・h・fs・Imag.tune・Damping time・Tune shift・partition・spin 関係の行は未出力)。
  7 桁目は `autofg('10.7')(1:10)` で切り捨てに見える (dz, f0)。f0=c·β/C。
  EMIT の後に DTSYNCH / PHICAV / EFFVCRATIO / EFFVC / EFFRFFREQ も書き出す。OMEGA0 は USE と EMIT で更新。EMITXE/EMITYE/EMITZE はエミッタンスの下限。EMIT は GAUSS を ON にする (§2)。
- **CHRO_MATICITY** (tffsa.f:609 → tfchro.f): 四極・六極・MULT の寄与を `1x,a12,1x,2f17.3` で 1 要素ずつ + `TOTAL`。BEND は数えない。六極単独の ξy は代入なので K2=0 で `-0.000`。pysad は `FFS.chro_cmd`。ESRF 960 要素 + MULT 入りの小リングで全行一致 (bench_level24b/24b2)。
- **T_YPE [pat …]** (tftype.f): 型 (ltyp の順: DRIFT, BEND, QUAD, SEXT, OCT, DECA, DODECA, MULT, SOL, CAVI, TCAVI, MAP, BEAM, SPCH, INS, COORD, APRT, MONI, MARK) ごとに ` ;` 行→型名→`NAME    =(KW =v KW =v )` を `twbuf` (タブ 7/5、行幅 lpw-2) で折り返す。値は `autos1` の最短表現、0 の項目は既定キーワード以外省略、`ROTATE` は DEG。MARK は twiss(l,0,:) (l=1 は保存値) と EMITX/EMITY/EMITZ/SIGZ/SIGE/DP。
- **DISP_LAY [opts] [loc [loc1]] [pat]** (tfdisp.f): オプション `ALL, A_CCELERATION, G_EOMETRY, GA_MMA, OG_EOMETRY, O_RBIT, R_MATRIX, B_EAM, P_HYSICAL, D_UMPOPTICS, Z, E_XTREMUM, RE_FERENCE, DRE_FERENCE`。既定ヘッダ `   AX      BX      NX      EX      EPX   Element   Length   Value      s(m)       AY      BY      NY      EY      EPY    DetR     #`、各行は 8.5 幅の tdispv (NX/NY は /2π)、Length 10.6、Value は既定キーワードの値、s(m) は 11.6 (値 0 は `tdtrimz` で `.0`)、Value は 10.7、`(a,a,a,1x,i5)` で位置番号。ヘッダは 66 行ごとに繰り返す (`nlc=66`)。パターンに一致した行だけを出す (`$$$` は一致したときだけ)。名前が 10 桁を超えると後ろの欄の連続空白を詰めて桁を保つ (`tfsqueezespace`)。
  s(m) は tfgeo の pos (F1/FB1/FB2 付き BEND/MULT のフリンジの長さ補正、OFFSET MARK の位置、§6)。Value 欄は既定キーのある要素なら値が 0 でも ` .0`、MARK や既定キーの無い要素は ` 0`。
  **位置の指定**: 語を 1 つずつ読み、`| < { * %` を含む語はパターン (表示する要素名の絞り込み `temat`)、`ALL` は範囲 1..nlat、モード語、それ以外は位置 (数式・`NAME`・`NAME.n`・`NAME.n±k`・`$$$`・`^^^`)。1 つ目が loc1、続く 1 語が位置なら loc2 (そうでなければ loc2=loc1 で、その語は次のコマンド)。位置でない語に出会うと何も表示せずに戻る。**範囲 (id1, id2) は FFS が保持** (USE で 1..nlat に戻る) し、位置を省いた `DISP;` や `DISP QF*;` は前回の範囲を使う。`loc1==loc2` はその要素から `$$$` まで (2 行)、`loc2<loc1` は nlat を超えて周回する。数式は `int(v+0.499)` (負は `nlat+1+v+0.5`) で 1..nlat に収める。`NAME` は 1 台だけの要素かファミリーの代表、`NAME.n` は同名が複数あるときの n 番目 (1 台だけの `RFC.1` は未定義)。未定義の `NAME.n` は ` ???-FFS-Error-?Undefined location NAME.n` + `Out[1]:= USE`。パターン (temat) は標準名 (`NAME`/`NAME.n`)、1 台だけの要素の `NAME.1`、複数ある要素の代表のファミリー名に合うもの (`2*47` のように `*` を含む語は数式でなくパターン)。モード R (R1〜R4)・O (DX DPX DY DPY)・E (β の極値のみ)・**G/OG (幾何: Gx Gy Gz, s, Length, Value, Chi1-3 (度)。OG は閉軌道を乗せた tforbitgeo。docs/history/work_geo_sig.md)** は実装、A/GA/B/P/D/Z と REFERENCE/DREFERENCE は未対応 (pysad の解析は docs/history/work_small_items.md B)。
  **反転配置 (LINE の `-`、`cmp%ori=.false.`) の要素名は `' -'//name`** (`max(10,lname+2)` 桁。通常は `lname<=8` で `'  '//name` を 10 桁、それより長ければ `' '//name`、tfdisp.f:438-448)。`-` が付くのは DISP の表だけで、
  `LINE["NAME"]`・`LINE["ELEMENT"]`・`Twiss[…]` の返値には付かない (`Twiss["NAME",k]` は SAD ではエラー)。`LINE["DIR",k]` が −1。pysad は parser が立てる `el._rev` を見る。66 行ごとに端末なら `(c_ontinue, q_uit, a_ll)?`。
- **SAVE [ALL] [pat…]** (tfsave.f): 現在値を設計値 `idval` に書き戻す (既定キーワードは `/errk`)。MARK は twiss(l,0,:) (l=1 は現在の初期値) と EMITX 等を保存。`EXPAND` も同じ。**RESET [ALL|TOTAL] [pat]** (tfrst) は設計値を現在値に戻す (`TOTAL` は `tffsresetall`)。
- **COUP_LE slave master coef** (tfcoup.f): 同型・同キーワードのみ。`icomp(slave)=master`, `couple(slave)=coef` で以後 slave の値は `master*coef`。`COUPLE *` は `tffssetupcouple`。`DECOUP_LE`/`INDEP_ENDENT` で解除。
- **REC_OVER**: 変数の `valvar` と `valvar2` (GO 直前の値) を交換して `tfsetv`。**MAXI_TERATION n**: `flv%itmax` (既定 40)。このビルドでの実際の効き方は §9。

## 6. 線形加速器 (TRPT) の加速と MARK の OFFSET

- TRPT かつ RFSW のとき、CAVI と VOLT≠0 の MULT は加速する (tfgeo1 の gammab = tfacc、qcav/qmult1 → tcave/tmulte の acc、qcopymatg の sqrt(gammab(l+1)/gammab(l))、qmat2twiss の gr)。RING では空洞は RFSW に依らずドリフト (docs/11 A4)。pysad は `pysad/trptacc.py`。
- `LINE["GAMMABETA", k]` = gammab(k) (p/mc、TRPT の加速分布込み)、`LINE["GAMMA", k]` = sqrt(1+gammab²)。
- OFFSET≠0 の MARK は位置 OFFSET+l (前後の要素の中の分数位置) の Twiss と s を自分の行にする (qcell.f:318-345)。
- CAVI の DX/DY/ROTATE は emit / CALC6D / TRPT で tcave の tchge (回してからずらす) として効く。追跡 (tcav.f) では SAD も無視する。

詳細は docs/history/work_trpt_acc.md。

## 7. SADScript 側の組み込み (tftwiss.f, tffsa.f tfffs, tftrack.f, tfoptics.f)

- **Twiss[key, loc]** (`tftwiss`): key は nlist の名前、`"*"/"ALL"` (全関数のリスト)、`"DALL"/"RALL"`、先頭 `D`/`R` 付き (`"DBX"` = 参照との差、`"RBX"` = REFERENCE 列)、`"LENGTH"` (=nlat)。`SIG…/SIZE…/GAMMA/GAMMABETA/S` は tfline に委譲。loc 省略で全位置のリスト、loc は実数 (端数は要素内の割合、`qtwissfrac` で補間)、文字列、文字列のリスト。`Twiss[...]=v` は参照 (ref) で代入可。
- **位置の解決 `tflinestk` → `tflinenamestk`**: 実数なら `floor` + 端数 (`ielma` で負数は末尾から)。文字列は `@` 前置で名前そのまま、`NAME+expr`/`NAME-expr` で位置オフセット (端数可)、`NAME.*` でファミリー全成分、`*%{|` を含めばパターン (temat)、`$$$`/`***`/`^^^`/`NAME.n` は `ielmf`。一致が 1 件ならパターンの形に関わらず実数 (`***`、`RFC*` も)、0 件は `{}`、複数はリスト (bench_level31b)。名前は文字列でも裸のシンボル (`{Q1A, Q2A}` の要素) でもよい (tfgetstrns)。要素の型名と同じ名前の LINE (`TEST = (…)`) も読める。
- **名前パターン** (tmatch.f): `* % {abc} {^abc} |`。`?` は普通の文字。pysad は `_wild_match` (`?` はワイルドカードではなく、`{^..}` に対応)。
- **LINE[key, loc]** (`tfline`/`tfline1`): `LENGTH`, `EXPAND` (値展開), `SPOSITION` (s→位置), `S|LENG` (端数で補間), `GAMMABETA`, `GAMMA`, `MULT`, `TYPE` (idtype 番号), `TYPENAME`, `NAME` (elname: 複数なら NAME.n), `ELEMENT` (ファミリー名, 末尾は `$$$`), `POSITION`, `GEO/GX/GY/GZ/GCHI1-3`, `OGEO/OGX…`, `DIR`, `BZS`, `BZT`, `UPDATE`, `DK` (errk), `NDIV`; その他は要素キーワード (`tflinek`: `LINE["K1","QF"]`)。`LINE[key,loc]=v` は代入 (ref)。K1.. は積分強度なので `LINE["L","QF"]=x` / `Element["L","QF"]=x` で K1 は保たれる (pysad は `set_value('L')` が `_PER_LENGTH` の積分強度を保存)。
- **FFS["cmd" [, unit]]** (`tfffs`): 文字列をファイルとして開き `tffsa(lfnp+1, …)` で実行。unit 省略時 `outfl=0` で **出力は抑制**、`-1` は現在の出力、OpenWrite 等の番号へ向ける。`CAL/GO/VAR/SHOW` は結果のリスト (`{dp, kind, {residual,xstab,ystab}, function-values}`) を返す。エラーは `FFS::error`。IF/REPEAT は 1 回の FFS[] 内で完結させる。`FFS$SHOW[]` = `FFS["SHOW"]` (fit 条件の平らなリスト)。
- **Element[key, name [, Saved->True]]** (`tfelement`): `LENGTH` (=nele), `EXPAND` (`tffsadjust`), `NAME`, `VALUE` (既定キーワード), その他キーワード。name はファミリー名。`Element / SetElement / FLAG / ExtractBeamLine / BeamLine / PrintBeamLine / WriteBeamLine / GetMAIN` (Packages/element.n, beamline.n, src/tftwiss.f) は pysad では `pysad/sadelem.py` (bench_level24f で 164 タグ一致)。
- **CalculateOptics[lb,le,inicond,cell,idim]** (tfoptics.f → qcell1): `{各位置の 28 個の Twiss, {stab}, {trace}}` を返し、大域の Twiss も書き換える。pysad は `FFS.fn_calculate_optics` (cell は周期解、cell=0 は `ins_twiss`、idim=3 は全周の CALC6D)。stabz と idim=3 の stab/trace は SAD で未初期化。
- **Emittance[]**: Packages (`CalculateEmittance` 内部関数 1004) で、`Stable, Tunes, EnergyLossU0, RfVoltageVc, EquilibriumPosition, MomentumCompaction, BucketHeight, OrbitAtExit, DampingRate, Emittances, MomentumSpread, BunchLength, EquilibriumBeamMatrix, OneTurnTransferMatrix` の規則リスト。
- **光学関数ほか** (SADHelp の仕様に合わせたもの): `TransferMatrix[p1, p2, Fold->False, Symplectic->False]`、`BeamMatrix[p]`、`NormalCoordinates[p]` (4D は R 行列+Twiss、CALC6D なら 6×6)、`OrbitGeo[p]`、`GeoBase[{χ1,χ2,χ3}]`、`SymplecticJ[n]`、`SymplecticInverse[m]`、`SurvivedParticles[beam]`、`AccelerateParticles[beam, p(n), {n, nmax}]`、`ParabolaRandom[n]`。位置 p は要素名か番号。(docs/history/work_trpt_acc.md 11 節)
- **TrackParticles[beam, dest, nbegin, nend]** (tftrack.f:35-100): beam は `{loc, {{x..},{px..},{y..},{py..},{z..},{dp..}[,{sz},{psi}],{flag..}[,{lost..}]}}`。dest は到達位置 (省略 nlat)、nbegin/nend で周回数 (`#4 >= #3`)。`loc` と dest は nlat-1 を法として巻く (ls∈[1,nlat]、ld∈[2,nlat]; ld<=ls なら 1 周足す)。mt=nend−nbegin+1 周: mt>1 なら ls→nlat、(mt−2) 周、1→ld、mt==1 なら ls→ld。戻りは `{正規化した ld, 座標}` (`merge(ld,ls,normal)`)。RFSW OFF では tturn (区間) ごとに z=0。NPARA≥2 で粒子 9 個以上なら fork。tapert (LOSSAMPL/LOSSDZ) で落ちた粒子は alive=0 になり、座標は落ちた時点の値 (limittrack で抑えたもの) で返る。初期フラグ 0 は追跡しない。ANGLE=0 で K1≠0 の BEND の追跡は tbend の tquad/tmulti 呼び出しに従う (K1<0 の符号を失う SAD の癖、docs/11 A38)。pysad `fn_track` は区間ごとに `track_particles_fast(els[a-1:b-1], …, circ=全周長)` (docs/history/work_zdp.md)。
- **SeedRandom[n]** (tfRandom_.c): 種を設定し `{"SAD", seed, …}` 形式の状態リストを返す。`tfaddseed` は第 2 要素に kseed を足して再設定。
- **Random[]** は RandomUniform0 (tfRandom_.c `REG8("Random", RandomUniform0)`) = random_plugin_sad.c の RANDOM_UNIFORM0: (seed + 2³¹)/2³² の [0,1)。RANDOM_UNIFORM (0,1) は +0.5 の式で別物 (こちらを使うと seed 列が稀に 2 ずれる)。GaussRandom[] は tgauss (Box–Muller、GCUT 付き)。
- **単位付き数値**: SADScript の `1GEV` `30DEG` `3MM` は SAD では `1*GEV` = 1 になる (docs/11 C24)。pysad は単位を掛けて読む (仕様差)。
- **FFS[] の中のエラー**: `FFS["…"]` の中のコマンドがエラーを出すと `???FFS::error: Error in FFS command n in FFS["…"]` で最上位の文 (Do など) を打ち切る。`LINE[key,"NAME"]=v` の未定義要素は `???General::invset`。qcod の溢れは止まらず NaN (直前の twiss があればそれ) で続行 (docs/11 C25, C26)。

## 8. pysad ffscmd.py の対応表

| SAD | pysad (`class FFS`) |
|---|---|
| 文の振り分け | `sadscript.Interp.run(src,'ffs')` → `_is_ffs_command` が先頭語を `FFS.is_word` で判定 (フラグ/`_COMMANDS`/`_ERRCMDS`/`_FUNS`(+M/I/SCALE)/要素名・パターン)。式 (`x=1`, `f[`) は SADScript へ |
| 1 文の解釈 | `execute(stmt)` → `_split_words` ({} [] "" 内は分けない) → `_exec_word(words,i)` が消費語数を返す。順序: フラグ → FIT 条件 (`_fit_condition`) → `FUNI`/`FUNSCALE` → 誤差 (`_error_cmd`) → コマンド → 要素値 (`_element_value`) → `?Undefined command or element:` |
| REPEAT/UNTIL, IF/ELSEIF/ELSE/ENDIF | `Interp.run` の最上位 FFS 文列の添字ループ + `_control_word`/`_control`、`_if_cond`/`_skip_if` (§1) |
| abbrev | `_abbrev(word,'CAL_CULATE')` |
| tflag | `self.flags` + `_FLAG_ALT` (sino/NOsino) |
| tgetfv | `_fit_condition` → `Matcher.fit(fun[+'M'], value, npoints)` (match6.py)。`{min,max}`, `*` 対応、`@` は未対応 (`SADError`)。`_fit_location` が `FIT [loc [loc1]]` (`FIT;` は位置を $$$ に戻すだけ) |
| FITP | `Matcher.fitp(n)` (|mfitp|>2 のみ変更)。FITP>1 の MatchingAmplitude/NFAMP も実装 (§3.2) |
| FREE/FIX | `_freefix` → `_free_family` (既定キーワードは `default_key`: DRIFT は L、他は `_KEYS[型][1]`、MARK/MONI は不可、vlim は DRIFT L [0.1,1e10]) / `Matcher.fix` |
| tfgetv | `_element_value` (R/MAX/MIN/MAXABS、キーワード、`set_family` で couple も反映、`_vlim`) |
| terror | `_error_cmd`/`_apply_error` (`errk` は `el._errk`、PUT の基準は `env._elements` の設計値、乱数は `self.rng.gauss` = track6 と同じ `SadRandom`、master グループは `_master_array`、`_tfadjst`)、`_dump` |
| USE | `use(name)`: `build_lattice`、**最初の USE だけ** 先頭 MARK の EMITX/EMITY/EMITZ/SIGZ/SIGE/DDP(→DP0)/DP を取り込む (tffsa.f:135 の `if(.not. chguse)`)。`USE BeamLine[..]` は `_use_beamline_expr`、`VISIT`/`BYE` は `_visit_stack` (§5) |
| CAL/GO | `go(fit)`: MAXITERATION/CONVERGENCE/ExponentOfResidual/OffMomentumWeight/DP/DP0/DPM/INTRES… を `Matcher` に渡し、`Matcher.go` (`Matched.`/`Residual =` 行は match6.py:781-792、`fmt_g` = 1PG11.4) → `print_fit_table` (tfshow 風)。解法は `pysad/tsvd.py`、VariableRange/VariableWeight は `Matcher.var_range/var_weight` (§3.4) |
| 微分 (GO) | `go()` が `FFS$NumericalDerivative` を `Matcher.numerical_derivative` に渡す。既定は `match6d.analytic_qu` (INS でも解析微分。未対応の条件・状況は中心差分に自動フォールバック: `Matcher._analytic_unsupported`) |
| EMIT | `run_emit` (emit6d) + `print_emit` |
| CHRO | `chro_cmd` (bench_level24b/24b2 で全行一致) |
| TYPE / DISP | `type_cmd` (型順・タブ 26,31,…) / `disp_cmd` (tfdisp の解析の直訳: 位置範囲・パターン・ALL・R/O/E/G/OG モード、bench_level29b 36 ブロック一致。A/GA/B/P/D/Z と REF は未対応。反転配置の `-` と桁は SAD と一致、docs/history/work_item6.md / docs/history/work_small_items.md B) |
| CALC6D | `calc()` が CELL + CALC6D で `calc6d.calc6d` を呼び `self.rows` (拡張 Twiss 28 成分) と `self.calc6d_result` を作る。分数位置は `row_frac` → `calc6d.twiss_frac` |
| TRPT の加速 | `pysad/trptacc.py` (§6) |
| INTRA / TouschekLifetime | `run_emit` が INTRA ON で `touschek.touschek_table` を呼び EMITX.. / EMITXC.. / `TouschekTable` をグローバルに入れる。`fn_touschek` = `touschek.touschek_lifetime` |
| GEO / ORG | `ensure_geo()` (= `geo.compute_geo`, GEOCAL/GEOFIX)、`line_value` の `GX GY GZ GCHI1-3 GEO OGX.. OGEO` と `S/LENG` (geo の pos)、`_org` (ORG [要素 gx gy gz chi1 chi2 chi3])。`Twiss["GX"]` は SAD と同じくエラー、γ は `GMX`/`GMY` |
| Emittance[opt] | `fn_emittance` = `emit_extra.emittance_rules` (mode: Matrix 3 ⊃ Orbit 2 ⊃ OneTurnInformation 1)。SaveEMIT->False (既定) では EMITX.. を元に戻す。`LINE["SIG"+key, l]` / `Twiss["SIG"+key, l]` は `beam_size_value`: CALC6D なら `beam_matrix_at` (6D 平衡ビーム行列)、そうでなければ `beam_matrix_twiss` (tfbeam.f の直訳: 基準座標の Twiss + EMITX/EMITY (MINCOUP 下限) + SIGE² を R1-R4 で物理座標へ。**pysad の仕様: EMIT の事前実行が必要**、docs/history/work_geo_sig.md)。`BEAM_SIZE` コマンドは EMIT 済みの確認だけ |
| CalculateOptics | `fn_calculate_optics` (§7) |
| SAVE/RESET | `save_cmd` (設計値 `env._elements` と `el._saved`) / `reset_cmd` |
| COUPLE | `_couple` → `self.couples[slave]=(master,coef)`; DECOUPLE/INDEPENDENT で全解除 |
| Twiss/LINE/Element/FFS/Emittance/TrackParticles/SeedRandom | `fn_twiss`/`fn_line`/`fn_element`/`fn_ffs`/`fn_emittance`/`fn_track`/`fn_seedrandom`、代入は `set_line`/`set_element` (Set をラップ)。位置解決は `_resolve_loc` → `positions` (tflinenamestk 相当: `$$$`, `NAME.*`, パターン, `NAME.n`, `NAME±k`) |
| Element/SetElement/BeamLine 系 | `pysad/sadelem.py` (§7) |

**既知の未対応** (`[pysad] FFS command … は未対応 (無視)` と出すもの): DRAW, GEO, SHOW, STAT, VARS, VARY, BEAM, MATRIX, ALIGN, ATT, INC, REC, SCALE, RENUM, REV, ZAP, DIMAD, SYNCHROB, REF, MEASURE, TRCOD, PRINT。そのほか:

- MCOD/TCOD/MON (この K64 ビルドでは FFS コマンドではなく式として評価されて何もしないので pysad も無視)。
- INS での fit (Matcher は `cell=False` で安定判定を省くだけで、SAD の INS 条件照合は未了)。
- FITP>1 のオフモメンタム条件と MatchingAmplitude (match6 に実装済みだが SAD との照合はない: §9)。
- FIT の `@`。
- DISP のモード A/GA/B/P/D/Z と REFERENCE/DREFERENCE。
- TYPE の PROF 等リスト値。
- `FUNSCALE` (読み飛ばし)。
- DBEAM。
- FFS の termes エラーでの外側 Do の中断、`FFS["…"]` 文字列内の REPEAT/UNTIL、一致なしパターンの無言処理 (docs/history/work_codcorr.md)。

## 9. この SAD ビルド (WSL, K64) で分かった注意

- **NPARA の既定は 4** (Packages/init.n:219 `NPARA=4;`。initb1.f の `RsetGL('NPARA',1)` を init.n が上書き)。TrackParticles は粒子 9 個以上で fork するため乱数列が分岐し、pysad と比較するには `NPARA=1;` をスクリプトに入れる。
- **FITP>1 の GO で SAD が再起動する** (`*** Run time Environment` バナーが再び出る)。NPARA=1 でも起きるので、オフモメンタム条件の参照値は得られていない。
- **NPARA>1 では出力が子プロセスで重複する** (tffsmatch の fork)。参照出力は NPARA=1 で取る。
- **FFS["…"] の unit 省略は非表示** (`outfl=0`)。画面に出すには `FFS["…", 6]` か `-1`。
- **MAXITERATION はこのビルドでは GO の反復数に効かない** (`MAXITERATION 1` でも 7 反復)。pysad はコマンドを受理する。
- EMIT は `codin=twiss(1,0,DX:DDP)` を初期軌道にし、終了後 `gauss=.true.`。DISP/TYPE の前には `tfsetparam` が呼ばれる。
