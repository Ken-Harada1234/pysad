# 02 SAD の MAIN レベル (ラティス定義) と pysad parser.py

最終更新: 2026-11-06 (現状のスナップショット)。照合対象: `<SAD>\src` (K64 版 SAD, Fortran) と `pysad/parser.py`。SAD ソースの行番号は 2026-10 時点。
表記: 無印 = コードで確認した事実 / 「解釈」= 推論 / 「未検証」= 未確認。parser.py は改訂が続くため行番号ではなく関数名で示す。

---

## 1. SAD 入力ファイルの構造 (MAIN レベル)

### 1.1 トップレベルの読み取り (`toplvl.f: toplvl`)

`gs file.sad` で起動した SAD は `MAIN.f` → `initbl.f` (コマンド表) → `initb1.f` (要素型・キーワード・単位・グローバル変数の表) を経て `toplvl` に入り、トークン単位で次のように振り分ける (toplvl.f:1700-1760)。

| トークン | 種別 (idtype) | 処理 |
|---|---|---|
| `DRIFT`, `BEND`, `QUAD`, … | icDEF (< icMXEL) | `doelem` — 要素定義ブロック (doelem.f) |
| `LINE` / `CELL` | icDEF (icLINE=100) | `doline` — ライン定義ブロック (doline.f, lread.f) |
| `FFS`, `TRACK`, `EMIT`, `QUICK`, `LIE`, `PLOT` | icACT | `doACT` → `ActTra` (track.f) など。`FFS;` 以降は FFS コマンド層 (tffsa.f) |
| `ON` / `OFF` / `READ` / `PRINT` / `STOP` / `EXPAND` / `REVERSE` / `LIST` | sethtbr 登録のサブルーチン (initbl.f) | `doonfl`/`dooffl` (doflag.f), `doread`, `dostop` など |
| 既存の要素名 / LINE 名 | 1..icMXEL / icLINE | `doelm2` / `dolin2` — 型キーワードなしの再定義 |
| それ以外の識別子 | — | `dAssgn` (dAssgn.f): `NAME = 式;` でグローバル実数 (icGLR) を定義・代入 |
| `;`, `}`, CR | — | 読み飛ばし |

- `!` から行末まではコメント (gettok.f:55 `if (curc == '!') call myfflush`)。文の区切りは `;`、ブロックは `{ }` でも囲める (`doelem`/`doline` は `{` を任意で受ける)。
- `FFS;` で ActTra に入り、FFS 内の `STOP`/`END` で MAIN に戻る。MAIN の `STOP;` (`dostop`) で SAD が終了する (bench_level12.result 末尾の `*** dostop *** Good Bye!`)。
- **gs の慣習**: ラティスファイル (MAIN 部のみ) と FFS スクリプトを連結して 1 本の `.sad` として `bin/gs` に渡す。本プロジェクトでは `bench_level12.sad` = `pysad/tests/kek-ls_v03_12.lat` + FFS 部であり、pysad 側は `load_sad_files(lat, script)` / `run_sad.py lattice script` が同じ連結を行う (docs/history/00_RESUME.md Level 14)。

### 1.2 `MOMENTUM=`, グローバル変数, 式と単位

- `MOMENTUM = 3GEV;` は `dAssgn` → `rdterm` → `yyparse` (calc.y) で評価される。文法は `NUM UNIT → $1*$2` (calc.y:65)、演算子 `+ - * / ^`、関数 `SIN COS TAN ATAN ATAN2 SINH COSH TANH SQRT EXP LOG LN` (yylex_.c:25-38)。
- 単位表 (initb1.f:1618-1670, `icUNIT`): `%`=0.01, `RAD`=1, `MRAD`=1e-3, `DEG`=π/180, `M`=1, `CM`=1e-2, `MM`=1e-3, `T`=1, `GAUSS`=1e-4, `JOULE`=1/e, `EV/KEV/MEV/GEV`=1/1e3/1e6/1e9, `V/KV/MV/GV`, `HZ/KHZ/MHZ/GHZ`。大文字小文字の両方が登録されている。
- 既定のグローバル (initb1.f:1697-1790, `defglb`+`RsetGL`): `MOMENTUM`=30e9 (**既定 30 GeV**), `MASS`=電子質量, `CHARGE`=0, `PBUNCH`=1e10, `NBUNCH`=1, `EMITX/Y/Z`=0, `SIGE`, `SIGZ`, `FSHIFT`=0, `DP`… など。
- 式の中の識別子が icGLR のグローバルなら値に置き換わる (yylex_.c:93-96)。したがって `QFK1=1.2; QUAD QF=(K1=QFK1);` のように要素パラメータに変数を使える。

### 1.3 `ON`/`OFF` 行 (doflag.f)

`ON ECHO; OFF LOG;` は `$ECHO$`, `$LOG$` という整数フラグ (icFLAG) に FLAGON/FLAGOF を入れる。未定義の名前は新規フラグとして定義され "New flag … is defined." と表示される (doflag.f:36-48)。MAIN の `$RAD$ $FLUC$ $COD$ $RADCOD$ $INTRA$ $POL$ $DAPERT$ $EMIOUT$ …` は FFS 入口 (track.f:80-107、`IgetGL('$RAD$')` 等) で FFS のフラグに写される (<work-item>/mflag*.sad での実測: `OFF RAD FLUC COD RFSW RADCOD;` で RAD/FLUC/COD/RADCOD が F)。**`$RFSW$` だけは届かない**: `isynch=IgetGL('$RFSW$')` を読んだ後、FFS に入るとき `rfsw=.true.` に固定される (track.f:149)。FFS 入口の表 `RFSW RADCOD RAD …` の RFSW 列は `$RFSW$` の値 (既定 F) で、FFS の `FLAG[{"RFSW"}]` は 1。pysad: `run_sad_string`/`main_exec` は RFSW を `main_flags` (MAIN の TRACK 用) に入れ、FFS の RFSW は ON のまま (docs/history/work_zdp.md)。既定で ON の MAIN フラグは FLUC と COD。

### 1.4 要素定義 (doelem.f, rdkwdl.f, rdexpr.f)

```
QUAD   QFA1 =(L =.25  K1 =1.0302859090451 )
       QFA1H=(L =.125 K1 =1.0302859090451/2 )
 ;
```

- `TYPE name=(kw=式 kw=式 …) name2=(…) … ;` — 名前はいくつでも続けられ、区切りは空白または `,`。
- 新しい名前: `idtype=型`, `idval=ktcaloc(kwMAX+1)` (doelem.f:43-48)。`ktcaloc` はゼロクリアするので (italoc.f:86-94) **書かなかったキーワードは 0**。
- キーワードは `CAPITA` で大文字化 (rdkwdl.f:40) し、icKWRD かつその型で有効 (`kytbl(kw,型)≠0`) でなければ "invalid key word name" / "Undefined keyword" でエラー (rdkwdl.f:43-66)。値は `rdexpr` → `yyparse` なので `L=.413/2`, `K3=-6000*0.0125`, `ROTATE=-90 DEG` のような式が書ける。
- **再定義 (後勝ち)**: 同じ型で同じ名前を再定義すると同じ格納領域に `rdkwdl` を再実行する (doelem.f:57-60)。つまり **書いたキーワードだけが上書きされ、書かなかったキーワードは前の値を保持**する。別の型で再定義すると "already used as …" エラー (doelem.f:61-66)。
- 要素の値が LINE 展開にコピーされるのは `USE` 時 (§1.5) なので、LINE を書いた後で要素を再定義しても反映される (pysad/tests/kek-ls_v03_12.lat は LINE の後に BEND/QUAD/SEXT を 3 回再定義しており、最後の値が有効)。

### 1.5 `LINE` 定義と展開 (doline.f, lread.f, expln.f)

```
LINE  ARCA = (ARCAH -ARCAH2)
      RING = (MKS ARCA, ARCA, … RFC, ARCA … MKE)
 ;
```

- `lread` が項目ごとに多重度 `mult` と名前を積む: `-` で `mult=-mult`、`N*` で `mult=mult*N` (lread.f:36-71)。したがって `-ARCAH2` は逆順、`2*SDA` は 2 回、`-2*X` は逆順 2 回。`+` は無視。項目は要素名でも LINE 名でもよく (入れ子)、定義時点で未定義の名前も許される (展開時に "is not defined yet" 警告, expln.f:110-116)。
- LINE の再定義は "… is redefined" の通知 (レベル −1) とともに旧定義を解放して **全体を置き換える** (doline.f:49-78)。要素の部分上書きとは挙動が違う。
- 展開は FFS の `USE` で `expnln` (expln.f:11-77) が行う: `explnstk(idxl, direct)` が再帰的に平坦化し、LINE 項目は `|mult|` 回・向き `dir=sign(direct*mult)` で再帰、要素項目は `|mult|` 回積む。**反転は要素順を逆にし、各要素インスタンスに `ori=.false.`** を付ける (expln.f:96-123, 51)。要素の値は展開時に定義からコピーされる (`tfvcopycmpall`)。
- `ori=.false.` の BEND では `tpara` が E1↔E2, AE1↔AE2, FB1↔FB2 を入れ替える (tffs.f:1259-1270)。QUAD/CAVI では FRINGE モードの入れ替え (tffs.f:1322-1325, 1447-1450)。


### 1.6 MAIN の `TRACK` / `EMIT` コマンド (詳細は docs/history/work_main_track.md)

`TRACK USE=RING NPART=21 TURNS=1000 NX=(1,10,20) NY=(0) NZ=(1);` は `initbl.f` が icACT として登録 (`ActTra.f` → `track.f` → `tracka.f`)。
キーワードは USE TURNS NPART XEMI YEMI SIGS SIGE ELEC POSI PROTON RADI SPAN CENTER AFTER AXIS SYNC NX NY NZ (ユーザースクリプトが使うのは USE NPART TURNS NX NY NZ だけ)。
粒子 j の規格化振幅 = (x,px,y,py,z,dp) = (NX[j], 0, NY[j], 0, 0, NZ[j]) (リストは繰り返し、無ければ x=y=1, dp=0)、初期座標 = r·(振幅×(emx,emx,emx,emx,emz,emz)) + codin
(r, codin, EMITX/Y/Z, SIGE は直前の `EMIT USE=..;` が作る古典 SAD の意味)。出力は 追跡前のフラグ一覧、損失の報告、smear の表 (SMEAR ON 既定)、tune (FFT ピーク 8 本) の表、`Tracking end:`、NX=/NY=/NZ= の確認。
`PHSPLOTS=n;` の TopDrawer 描画は無視してよい。
**K64 版では MAIN の EMIT は何もせず (`Use EMIT_TANCE or Emittance[] within FFS.`)、TRACK は `tinip` で SIGSEGV する**ので SAD で出力を取れない。
pysad は `parser.parse_main_commands` (ON/OFF, EMIT, `NAME=値`, TRACK, STOP を順に読む) + `FFS.main_exec` + `maintrack.run_track` で古典の意味に実装している (SADScript で同じ手順を書いた版と一致)。

---

## 2. 要素型とキーワード

型コード (tfstk.f:277-285): DRIFT=1, BEND=2 (`ST` も BEND, initb1.f:1004), QUAD=4, SEXT=6, OCT=8, DECA=10, DODECA=12, SOL=20, MULT=22, CAVI=31, TCAV=32, MARK=41, MONI=42, APERT=43, LINE=100。
キーワード名の別名 (initb1.f:729-994): `ANGLE`, `ROT`=`TILT`=`ROTATE`, `DISRAD`(RAD), `DISFRIN`(FRIN), `FRINGE`(FRMD), `ACHROMA`(CHRO), `RANKICK`, `ALPHAX`=`AX`, `BETAX`=`BX`, `PSIX`, `EMIX`=`EMITX`, `SIGZ`=`SIGMAZ`, `DROT`=`DROTATE`, `K0`…`K21`, `SK0`…`SK21`。

| 型 (initb1.f) | キーワード (格納順 `ky_*`) | 単位・意味 (根拠) |
|---|---|---|
| DRIFT (995-1004) | L, DISKIN, COUPLE, RADIUS, AX, AY, DISRAD | L [m] |
| BEND (1008-1042) | L, ANGLE, E1, E2, ROTATE, DISRAD, K1, DX, DY, K0, DISFRIN, EPS, RANKICK, F1, FRINGE, COUPLE, DROTATE, AE1, AE2, FB1, FB2, INDEX, AX, AY, LDEV, LRAD, CHI2, DZ | ANGLE [rad] = 基準座標系の偏向角 (phi0)。**K0 [rad] は追加の双極場で、磁場の偏向角は phib = ANGLE+K0** (tbend.f:565, qtwiss.f:281)。**E1/E2 は ANGLE に対する比**: 実エッジ角 psi1 = E1·ANGLE + AE1 (tffs.f:1261)。**K1 は積分値 K1·L [1/m]** (tbend.f:566-573 で tquad に渡す)。ROTATE [rad] = 設計ロール角 (tffs.f:1294 で cos/sin を取る)。DX/DY [m] 据え付けずれ。DISFRIN=0 でフリンジ有効 (tturn.f:342)。FRINGE ≤0 のとき FB1/FB2 を 0 にする (−1 で入口のみ, −2 で出口のみ残す; tffs.f:1271-1278)。F1+FB1/FB2 は有効長補正 `al −= ((φ·fb1)²+(φ·fb2)²)/al/48·sin(w/2)/sin(φ/2)` (tffs.f:1280-1284)。 |
| QUAD (1044-1073) | L, K1, ROTATE, DX, DY, DISRAD, ACHROMA, DISFRIN, F1, F2, FRINGE, EPS, DISKIN, COUPLE, INDEX, AX, AY, LDEV, LRAD, F1K1F… | **K1 は積分値 K1·L [1/m]**: `akk = K1/al` (tffs.f:1316)、L=0 なら薄レンズ強度 (tquad.f:22-24)。ROTATE [rad] (tffs.f:1312)。 |
| SEXT/OCT/DECA/DODECA (1076-1106, THIN 共通) | L, K2/K3/K4/K5 (同じスロット `ky_K_THIN`), ROTATE, DX, DY, DISRAD, DISFRIN, COUPLE, INDEX, AX, AY, LDEV, LRAD | K は積分値 (tthine で K/fact(n)/2 のキック; Level 3A で SAD の色収差と一致, docs/history/00_RESUME.md)。 |
| MULT (1110-) | L, DX, DY, DZ, CHI1, CHI2, CHI3=ROTATE, EPS, DISRAD, DISFRIN, F1, F2, FRINGE, VOLT, DVOLT, HARM, PHI, FREQ, COUPLE, ANGLE(23), K0..K21 (30,32,…), SK0..SK21 (31,33,…) | K_n の単位は未検証 (解釈: 積分値。tffs.f:1440 の分割数推定は K/al を使う。tmulti 本体は本文書では未照合)。 |
| SOL (1239-1258) | L, BZ, DX, DY, DZ, DPX, DPY, BOUND, CHI1-3, GEO, F1, DBZ, COUPLE, DISFRIN, DISRAD, FLAT | BZ [T] |
| CAVI (1268-1295) | L, VOLT, DVOLT, HARM, PHI, FREQ, RANVOLT, RANPHASE, LWAKE, TWAKE, DX, DY, ROTATE, V1, V20, V11, V02, COUPLE, DPHI, DISFRIN, FRINGE, AX, AY, LDEV, AUTOPHI | VOLT [V], FREQ [Hz], PHI [rad] (**既定 0**: ゼロクリア)。波数 `w = 2π·FREQ/c` (HARM=0) または `ω0·HARM/c`、`ω0 = 2π c p0/h0/C` (tffs.f:1452-1454, tfsetparam.f:124) → FREQ = HARM·βc/C。 |
| MARK (1512-1556) | AX, BX, PSIX, AY, BY, PSIY, EX, EPX, EY, EPY, R1-R4, DETR, DX, DPX, DY, DPY, DZ, DDP, AZ, BZ, PSIZ, ZX, ZPX, ZY, ZPY, DP(29), OFFSET, SIGZ, SIGE, GEO, JDX…, EMITX(40), EMITY, EMITZ, COUPLE | 先頭 1..28 は FFS の mfit 順と同じ (ky_AX_Mark=1; tffsmatch.f:133 `rlist(latt(1)+mfitddp)`)。`tffsinitparam` (tffs.f:2904-2965): 先頭 MARK の EMITX/EMITY/EMITZ/SIGE/SIGZ が >0 ならグローバルを上書き、DP ≤0 なら `DP` 変数 → それも 0 なら 0.01、BX/BY/BZ ≤0 なら 1。 |
| MONI (1576-1582) | DX, DY, OFFSET, COUPLE, ROTATE | **L なし** (L を書くとエラー)。 |

---

## 3. pysad parser.py の対応

### 3.1 構成

| 部品 | 役割 |
|---|---|
| `tokenize` / `_TOKEN_RE` | 正規表現トークナイザ。`(* *)` と `!…行末` をコメントとして捨てる。NUMBER は `.05` 形式も可、IDENT は `$` と `` ` `` を含む。 |
| `_Parser.parse_program` → `parse_statement` → `_try_special_form` | 文ごとに特別形を判定: `SetElement[…]`、要素型キーワード (`_SAD_ELEMENT_TYPES`、DECA/DODECA を含む) → `_parse_native_element_block`、`LINE` (直後が `[` でない) → `_parse_line_block`、`x = BeamLine[…]`、`Twiss/Emittance/RadiationIntegrals[…]`。それ以外は式/代入 (`MOMENTUM = 3GEV;` → `env.set('MOMENTUM', 3e9)`)。エラーは警告を出して次の `;` まで読み飛ばす。 |
| `_parse_native_element_block` / `_parse_sad_params` | `TYPE name=(kw=式 [単位] …) …` を読む。`;`、次の型キーワード、`LINE` でブロック終了。型名と同じ要素名 (`CAVI CAVI=(…)`) も可。`kw` だけで `=` がなければ `True`。単位は `_SAD_UNITS` (§1.2 の initb1.f の単位表: DEG/MRAD/RAD, M/CM/MM, EV/KEV/MEV/GEV, V/KV/MV/GV, HZ/…/GHZ, %, T, GAUSS, JOULE)。`env.define_element(name, type, dict)`。同じ型で同じ名前を再定義したときは書いたキーワードだけを上書きする (doelem.f と同じ、§1.4)。 |
| `_parse_line_block` / `_parse_beamline_refs_sad` | `LINE name=(…) …;`。項目は `ID`、`N*ID` → `('REPEAT', N, ID)`、`-ID` → `('REVERSE', ID)`、`-N*ID` は反転 + 繰り返し。 |
| `_parse_primary` | 数値直後の単位 (`3GEV`, `-90 DEG`) を掛ける。未知の識別子は文字列として返す。 |
| `SADEnvironment` | `_vars` (組込み定数+グローバル), `_elements` {名前: {TYPE, kw: 値}}, `_beamlines` {名前: 参照列}, `_ffs_state` (Twiss[]/LINE[] のディスパッチ先)。 |
| `SADEnvironment.build_lattice` | `_expand_refs` で平坦化 → `_build_element` → `Lattice`。最後に周長 C を取り、`_harm` を持つ空洞に `FREQ = HARM·c/C` を入れる。 |
| `_expand_refs` | 再帰展開 (深さ 50 まで)。`REVERSE` は参照列を逆順にし `is_reversed` を反転して再帰、要素は `('REV_EL', name)` にする (二重反転で正順に戻る)。 |
| `_build_element` / `_build_element0` / `_apply_extra_params` | 下表の変換で pysad 要素を作り、`_apply_extra_params` が DX, DY, K0, F1, FB1, FB2, FRINGE, DROTATE, AE1/AE2 などを属性に渡す (BEND の K0 はファイルで定義しても読まれる = bench_level23h)。数値パラメータ全部を `el._raw` に残す (FFS の `LINE["…"]` 表示用, ffscmd.get_value)。QUAD/SEXT/OCT の ROTATE, DX, DY, K0, F1, F2, FRINGE, DISFRIN、CAVI の PHI, DX, DY のうち属性に無いものは `_raw` のみで、読み込み時の計算には効かず、FFS 側 (`ffscmd.set_value`, `ffs.set_element`) で設定した値だけが効く。SAD の既定 (PHI=0, FRINGE=0, DISFRIN=0 など) と一致する限り結果は変わらない。 |
| `SADScript` / `load_sad_file` / `load_sad_string` / `load_sad_files` | 公開 API。`load_sad_files(lat, script)` は複数ファイルを同じ環境に順に読む (= gs への連結)。`.ffs` / `run_ffs` で FFSState に接続。 |

### 3.2 `_build_element0` の変換 (SAD 値 → pysad 属性)

| SAD | pysad | 変換 |
|---|---|---|
| DRIFT L | `Drift(L)` | そのまま |
| QUAD L, K1 | `Quadrupole(L, K1)` | **K1 → K1/L** (L>0; pysad の K1 は長さあたり [m⁻²])。`K1L` キーワードがあればそれを /L。 |
| BEND L, ANGLE, K1, E1, E2, AE1, AE2, ROTATE | `Bend(L, ANGLE, K1, E1, E2, ROTATE)` | K1 → K1/L。**E1 → E1·ANGLE + AE1, E2 → E2·ANGLE + AE2** (SAD の比→ラジアン、tbende.f:448)。逆向き通過 (`REV_EL`) では E1↔E2。ROTATE は Bend では設計ロール角 (elements.py Bend.transfer_matrix)。**ANGLE=0 でも ROTATE があれば Drift にせず `Bend(ANGLE=0, ROTATE=…)` のまま保持する** (esrf の `VD=(ROTATE=-90 DEG)` や pf2 の VD/PV/VS に K0 を与えたとき効かせるため)。SAD では ANGLE=0, K1=0 の BEND は tsteee/tsteer (K0=0 ならドリフト、L=0 なら薄いキック `px −= K0`) で、ROT は tchge で回す (docs/history/work_rotate.md)。BEND の ROT は `Bend.ROTATE` (設計ロール)、DROT は `DROTATE`。FFS の `LINE["ROTATE",…]` もそれぞれの属性に入る。 |
| SEXT K2 / OCT K3 / DECA K4 | `Sextupole(K2/L)`, `Octupole(K3/L)`, `Decapole(K4/L)` | 積分値 → 長さあたり (L=0 ならそのまま) |
| CAVI VOLT, FREQ, HARM | `Cavity(VOLT, FREQ)` + `_harm` | FREQ=0 かつ HARM>0 なら `_harm` を記録 → build_lattice で FREQ = HARM·c/C。物理計算 (emit6d.py, track6.py) は `_harm` と ω0 から SAD と同じ `w` を作るので β=1 近似の FREQ は使われない。 |
| MONI / MARK | `Monitor(L)`, `Marker()` | MARK の Twiss キーワード (BX, AX, EMITX, DP, …) は `_raw` と `env._elements` に残り、ffscmd/ffs が参照する。 |
| MULT | `Multipole(L)` | K0..K21 は属性に渡らない (`_raw` のみ)。 |
| SOL / 不明な型 | `Drift(L)` | BZ 等は無視 |

### 3.3 既知の制限・SAD との差

1. **型が違う再定義**: SAD は "already used as …" エラー (§1.4)、pysad は黙って型を変える。
2. MAIN 式の関数 (`SIN( )`, `SQRT( )`, `ATAN2( )` …) は未対応 (四則と `^` のみ)。未定義の要素名は SAD が警告するのに対し pysad は黙って飛ばす。
3. `FFS;` 以降を同じファイルに含めても `USE RING; CELL; CALC;` などは読み飛ばされるが、`Print[…]`, `LINE["…"]`, `FFS["…"]` はその場で評価される。FFS 部は `run_sad.py` / `sadscript.py` + `ffscmd.py` で扱う (docs/history/00_RESUME.md Level 14)。

---

## 4. 実例: ESRF 型 20 セル (Level 12)

- ラティス: `pysad/tests/kek-ls_v03_12.lat` (209 行)。`MOMENTUM = 3GEV;` の後に DRIFT/BEND/QUAD/SEXT/OCT/CAVI/MARK/MONI を定義 (式 `L=.413/2`, `K3=-6000*0.0125`, `VD=(ROTATE=-90 DEG)`, MKS/MKE に Twiss・EMITX・DP)、`LINE ARCAH … ARCAH2 … ARCA=(ARCAH -ARCAH2) RING=(MKS 20×ARCA + RFC MKE);`、その後 BEND/QUAD/SEXT (+DRCT) を 3 回再定義 (最後が有効。再定義は全キーワードを書いている)。
- SAD 側: `bench_level12.sad` = 上記 + FFS 部 (`USE RING; LINE["NAME"|"TYPE"|"L"|"K1"|"ANGLE"|"K2"|"K3"|"E1"|"ROTATE","*"]` を出力し `CELL; CALC;` → Twiss、`RFSW; emit;`)。結果 `bench_level12.result`。
- pysad 側: `bench_level12_pysad.py` が `parser.load_sad_file(pysad/tests/kek-ls_v03_12.lat)` → `build_lattice('RING')` し、SAD の `LINE[]` 出力と要素ごとに比較する (pysad の K1·L, K2·L, K3·L を SAD の積分値と、Bend.E1 を SAD の E1×ANGLE と比較)。結果 `bench_level12_pysad.result`:

| 項目 | pysad | SAD | 差 |
|---|---|---|---|
| 要素数 / 周長 | 3563 / 603.120000 m | NLAT−1 = 3563 / 603.120000 m | 一致 |
| 要素名 (順序込み) | 不一致 0 | | |
| L, ANGLE, K1, K2, K3, E1, ROTATE | 最大差 K1 4.9e-11, ANGLE 4.0e-11, K2 4.5e-11 | (`$FORM="15.10"` の 10 桁出力) | 表示桁の丸め |
| Twiss 末尾 `$$$` (BX, AX, BY, AY, EX, EPX, NX, NY) | | | ≤ 4e-11 |
| Twiss `MKC` | | | ≤ 4e-11 |
| emit (εx, σE, σz, α, U0, Qx/Qy/Qz, 減衰) | | | 相対差 ≤ 2e-5 (SAD 5 桁出力) |

ANGLE=0 の `VD`/`HD` の扱いは §3.2 (Bend のまま保持) のとおりで、要素数・名前・周長は SAD と一致する。反転ライン `-ARCAH2` の中の要素も (`REV_EL` 経由で) 同じ順序になることをこの比較が保証している (ARCAH2 の BEND に E1/E2 は無いので入替の効果は Level 12 では見えない)。

**BEND の K0 / E1 の FFS 側の意味 (bench_level23)**: `LINE["K0","BPI6.1"]=v` は読み出しでは 1 台だけ変わるが、次の CALC / Emittance[] / TrackParticles でファミリー全体が最初の 1 台 (.1) の値に揃えられる (`.3` だけ設定した値は CALC 後に消える)。`LINE["K0","BPI6"]` (ドット無し) や FFS の要素値コマンド `BPI6 K0 v` はファミリー全体に即座に効く。`BPI6.1 K0 v;` と書くと SAD は式として評価して何もしない。`LINE["E1",name]` は ANGLE 比 (esrf の BEND なら 0、pf2 の BMS なら .5) を返す (ANGLE=0 の BEND でもファイルの値)。pysad は `ffscmd.FFS.spread_pending`、`_e_ratio` で対応。
