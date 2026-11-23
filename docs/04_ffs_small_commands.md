# 04. FFS の小さなコマンドと DELCOR

**この文書で分かること**: FFS の小さなコマンド (VARIABLES・MATRIX・SHOW・PRINT・OUTPUT など、1 行で使えるもの) の働きと pysad での扱い (§1)、リングに沿った相関のあるずれを入れる DELCOR の使い方・作り方・使える範囲 (§2)、限界 (§3)。
**前提**: [03_ffs_commands.md](03_ffs_commands.md) の FFS の基本 (USE・CALC・FIT・FREE) を読んでいること。
**読む順**: §1 の表で目的のコマンドを探す。誤差入りのリングを作りたい人は §2。

用語: **FFS** = SAD の対話コマンド層。**K64** = 翻訳の元にした SAD (Fortran)、**gs** = その実行ファイル。**k64 モード** = `PYSAD_COMPAT=k64` で K64 と同じ動作、**v2 モード** = pysad の既定 (K64 の誤りを直した動作)。

小さなコマンドは K64 の tfvars.f・tfattr.f・tfmat.f・tshow.f・tffsa.f の直訳で、実装は `pysad/ffsextra.py`。表の書式は autofg (`'12.9'` など) まで gs の出力と同じ (`tests/test_ffsextra.py`)。

## 1. 小さなコマンドの表

下の表は、FFS の小さなコマンドと、その入力・出力を一覧にしたもの。

| コマンド | 内容 |
|---|---|
| `VARIABLES` (`VARS`、`VAR`) | FREE した変数の表 (Now、Previous、Saved、Minimum、Maximum、結合先、係数)。`FFS["VAR"]` は同じ値のリスト |
| `ATTRIBUTE` (`ATTR`) | 要素の既定キーの値と限界値 |
| `MATRIX [SYMPLECTIC\|PHYSICAL] [from to]` | 転送行列 (CALC4D は 4×5、CALC6D は 6×6) |
| `SHOW` | 現在のマッチング条件 (`FFS["SHOW"]` はリスト) |
| `VARY keyword pattern` | 要素の既定キーを変える |
| `RECOVER` | FREE した要素の値を、直前の GO の前の値と交換 |
| `REVERSE` | 入口の AX AY EPX EPY R2 R3 DPX DPY の符号を反転 |
| `PRINT expr`、`EXECUTE "文字列"` | 式の値を `Out[n]:= 値` で出す、文字列を FFS コマンドとして実行 |
| `OUTPUT(OUT) file` / `APPEND(APP) file` / `CLOSE OUTPUT` | 出力 (FFS の表示と `Print`) をファイルへ切り替え・追記・戻す。ファイル名は引用符つき、またはファイル番号 |
| `INPUT(IN) file` / `READ file` | ファイルの文を続けて実行 (末尾で元の入力に戻る) |
| `STATUS` (`STAT`) | フラグ・特別変数の一覧 (twbuf と同じ 8 桁の桁揃え)。CTIME・DT は実行時間なので gs と一致しない |
| `CLOSE`・`TERMINATE`・`SUSPEND`・`RESUME` | 入力側は何もしない (入力はファイルを丸ごと実行するため) |

### 1.1 動作の補足

- `MATRIX` の名前だけの引数 (例 `MATRIX 要素名`) は `NAME.1` として扱い、ラティスに無い名前は読み飛ばす (to=`$$$`)。
- `DUMP pat`: 語の変数は `character*20` で、パターンは 20 文字で切れる。合う要素がなければ全要素。DK・DDK 列は Fortran の `1p g12.5` (0 は `0.0000`、小さい値は `4.22829E-04`) で、`ffscmd._gwd` が同じ書式にする。DDK 列は K0 のキーワードのある要素 (BEND・MULT・TCAVI) だけ (QUAD にはない)。
- FIT のフラグ UNIPOL・RELW: UNIPOL は変数の符号を変えず (tffsmatch.f:669)、RELW は重みを相対にする。pysad は `ffscmd.matcher_` で match6 のフラグとつなぐ。
- `DISP` のモード R・O・P・GA・Z・A・B (BEAM、EMIT 後)・D・G・OG と REFERENCE・DREFERENCE (見出しの文字の置換は tfdisp.f の trefheader/tdrefheader のとおり)、`REFERENCE` コマンド、`Twiss["R…"/"D…"/"RALL"/"DALL"]` を実装 (詳細は [03_ffs_commands.md](03_ffs_commands.md) §5・§8)。
- `RENUMBER comp` と `LINE["MULT"]` (tfinimult の直訳。同名 1 台だけは 0)、`USE [NO]EXPAND`・`CALC|GO [NO]EXPAND` を実装 ([03_ffs_commands.md](03_ffs_commands.md))。`NOKEEPEXP` (tffsa.f の `calexp = .not. tfvcomp()`) と RENUMBER の INDEX 属性 (2 以上なら同名の要素の番号をその値から振る。`LINE["INDEX",…]` も読める) を実装 (両モード)。`tfvcomp` は、既定でないキーの成分の変数 (成分ごとの FREE) があれば真 ([03_ffs_commands.md](03_ffs_commands.md) §3.3)。
- `QF.3 K1 0.9;` のように `名前.番号` で始まる文は、K64 では式として評価されて何も設定しない。pysad も何も設定しない (`Out[n]:=` 行は再現しない)。
- K64 の tffsa.f に無いコマンド (MON・MCOD・TCOD): FFS では式として評価されて何も起きない。v2 モードは「K64 に無いコマンド」と警告して文全体を読み飛ばす (機能は実装しない)。
- フラグの効果: GAUSS は MAIN の TRACK・MEASURE の初期分布、CTIME は CPU 時間の表示、CONV は収束の判定・表示、EMIOUT は EMIT 報告の追加行 (pysad は spin の行だけ対応)、BIPOL は tffsmatch.f だけ。

## 2. DELCOR (相関のあるずれ)

K64 の `tcorr.f` と `tvcorr.f` (FFT は `tcftr.f`)。実装は `ffscmd._delcor`、`_tvcorr`。

```
DELCOR [ADD] corlength delx dely element [element1 ...];
```

リングに沿った相関のあるランダムなずれを、指定した要素の DX・DY に入れる。`corlength` は相関長 [m]、`delx` `dely` は大きさ。ADD が付けば足し、なければ置き換える。


### 2.1 k64 モードと v2 モードの違い

以下の 3 点が違う (それ以外は両モードで同じ)。k64 モードは K64 (gs) と同じ動作、v2 モード (既定) は修正した動作。

| | k64 モード | v2 モード |
|---|---|---|
| MARK・CAVI の書き込み先 (B40) | MARK は DETR・DPX、CAVI は HARM・PHI (符号も逆) | MARK・CAVI も DX・DY に、ほかの磁石と同じ符号 (+delx) で入れる |
| corlength≠0 の波数 −1 の成分 (B62) | `cv(n)` が未設定 (2 回目以降は古い記憶域の値)。負の波数の重みが 1 つずれる | 正負の波数で対称に置く (乱数を引く順は同じ) |
| corlength≠0 の規格化 (B99) | `s = 1/sqrt(4Σw²+2)`。複素数 c の絶対値の 2 乗平均が 1 なので、DX の rms は delx/√2 | `s = 1/sqrt(2Σw²+1)`。実部・虚部それぞれの分散が 1。DX の rms = delx、DY の rms = dely |

corlength = 0 の枝 (ランダムウォーク、1 m あたり分散 1) は両モードでビット単位で同じ。


### 2.2 作り方

- 周長 `al = pos(nlat)` (INS は 2 倍)。点数 `nc = 128` から `nc/2 < al` になるまで 2 倍 (最大 65536)。
- `tvcorr`: `corlength=0` は 2 つのガウス乱数を足していくランダムウォークで、両端が合うように直線を引く。≠0 は周波数空間に `exp(-(k·2π cl/(2 al))²/4)` の重みのガウス乱数 (実部・虚部) を置き、
  `tcftr` (exp(+2πi jk/n) の和 = `n·ifft`) で戻して規格化の定数を掛ける (k64 は `1/sqrt(4Σw²+2)`、v2 は `1/sqrt(2Σw²+1)`)。乱数は `tgauss` と同じ列 (`SadRandom`、GCUT 適用)。呼び出しの最初に必要な分をまとめて引く。
- 要素ごとに、中心位置 s の周期的な線形補間の値 (実部 × delx、虚部 × dely) を入れる。同じ要素の連なり (master) には同じ値、ワイルドカード以外は最初の 1 つだけ。最後に `tfadjst` (隙間なく並んだ磁石は先頭と同じずれ)。
- 対象は BEND QUAD SEXT OCT (と MARK、CAVI)。DECA・DODECA・MULT・MONI・ドリフトは `Undefined element` 扱い (K64 の型の判定どおり)。


### 2.3 v2 モードの根拠と使える範囲 (統計・解析式。`test_v2_C_delcor.py`、`test_v2_C2_delcor.py`)

- 規格化: 場の rms (500 m、512 点、400 標本) は x0 = 3、15、5000 m で実部・虚部とも 1 ± (2、4、8)% 以内 (k64 は 0.708、0.703)。DELCOR コマンド (12 セルのリング 62.4 m、四極 24 台、corlength 5 m、150 回) で rms(DX)/delx = 1.021、rms(DY)/dely = 1.043 (統計誤差 約 2.6%。k64 は 0.713、0.746)。
- 相関関数 ⟨c(s) c*(s+Δ)⟩/⟨|c|²⟩ = exp(−Δ²/(2x0²)) が 5% 以内 (虚部 2%)。相関長 1000 km ではリング全体が同じ値 (差 < 1e-9) で、rms 1.0 ± 5% のガウス分布。
- B40: 相関長 ≫ 周長なら全要素のずれが同じ値・同じ符号になる (QF・RFC・ST が 5% 以内。相関長 5 km では全桁同じ)。
- 使える範囲: v2 モード、corlength > 0。rms は多数回の平均としての値で、GCUT を掛ければその分小さい。相関長が格子間隔 (周長/nc、nc ≥ 128、INS は周長の 2 倍) より短いと線形補間で rms が下がる (間隔 = 相関長/4 で約 1%)。corlength = 0 の枝は別の量 (分散 = delx²·s(1−s/周長) [m]) で、delx の意味は「1 m あたりの rms」。
- 未確認・判断: 規格化 (delx を DX の rms と読む) は作者の文書による裏づけが無い推定 (SADHelp に DELCOR の項は無い)。K64 の値が必要なら k64 モード、または delx に 1/√2 を掛ける。CAVI の符号 (+delx) は決定済みだが、K64 が SOL と同じ形の式 ((設計の DX) − delx) を写した意図があるなら逆。
- 例: [`example/02_optics_cod/flag_DELCOR_rms.sad`](../example/02_optics_cod/flag_DELCOR_rms.sad)。


### 2.4 k64 モードの動作と限界

- **MARK には DX・DY ではなく DETR と DPX が入る (k64 モード。v2 は DX・DY)。** K64 は `latt(j)+15` と `+17` に書くが、値の添字では `ky_DETR_Mark=15`、`ky_DPX_Mark=17` (`ky_DX_Mark=16`、`ky_DY_Mark=18`)。
  一つずれた書き方に見える (docs/30 B40)。gs と合わせるためそのまま写した。INS の始点 (MARK の DETR・DPX) が変わる。
- SOL は master が 0 なので K64 でも対象にならない (何も起きない)。**CAVI** は (k64 モードでは) K64 の添字のずれで DX・DY ではなく **HARM と PHI** に書く (tcorr.f の `le=latt+3`、`le1=le+1`)。add でなければ設計の DX・DY (`idvalc+13`, `+14`) から delx・dely を引いた値、add なら -delx・-dely を足す。gs で確かめた (`CAVI RFC=(VOLT=1.5e6 HARM=64 DX=1e-5 DY=2e-5)`、`DELCOR 3 1e-4 1e-4 RFC` で HARM=-.00015145432812411818、PHI=.00011480572629164664。pysad は同じ値)。HARM は実数のまま持つ。
- 各要素の DX・DY の添字 (QUAD は 5、BEND は 9 など) は gs の結果で合うことを確認 (QUAD・SEXT・BEND)。OCT は対象に含めた DELCOR で、チューンと QF の DX が gs と 12 桁で一致した (2026-11-16 夜)。ただし OC 自体の DX は gs も pysad も 0 を返すので、OCT への効果そのものを見る強い比較ではない。
- 'DELCOR' 以外の Rule 系の誤差コマンド (ERROR の DELX など) は従来どおり (docs/29)。
- **2 回目以降の `DELCOR` で corlength≠0 の場合、gs と値が一致しない (乱数の列は同じ)。k64 モードのこと。v2 は cv(n) を対称に置くので、この不一致は起きない (gs とは違う値になる)。** K64 の `tvcorr` の corlength≠0 の枝は `cv(n)` (最後の要素) を設定しない
  (`cv(2..n/2)` と `cv(n/2+1..n-1)` と `cv(1)` だけ) ので、`tcftr` は `ktaloc` で得た記憶の古い内容を読む。1 回目は 0 で pysad (0) と一致するが、2 回目以降は前の確保の中身が入って
  数 % 違う (corlength=0 の枝と ADD の回は一致。確かめた 4 通りの組み合わせ: 0→3 と 3→2 と 3→3 は不一致、3→0 と 0→0 は一致)。再現できないので pysad は 0 のままにし、docs/30 に B62 として残した。

**検証**: CELL・INS のリングで、DELCOR のあとのチューンと閉軌道が gs と一致する (CAVI の書き込みも)。詳細は [history/04_ffs_small_commands_validation.md](history/04_ffs_small_commands_validation.md)。

## 3. 限界・K64 との差

- `PRINT` の `Out[n]` は、USE が `Out[1]:= USE` を出したものとして 2 から数える (pysad は USE でこの行を出さない)。
- `MATRIX` の 5 列目 (分散の項) は η(l2) − M·η(l1)。要素が 1e-16 程度の雑音を持つ箇所は gs が `-4.798E-16` などと出すのに対し、pysad は 0 を出す。
- `MEASURE`・`TRCOD`・`LOGMEASURE`・`ORBITJITTER`・`PUTMEASURE`・`GETMEASURE`: **v2 モードだけ**実装した ([`pysad/measure.py`](../pysad/measure.py)。使い方は docs/03 の「v2 モードの MEASURE・TRCOD」)。k64 モードは未対応 (警告して無視)。`MEASURE` は gs が何も出力せず (B42、`TRCOD` は SIGSEGV)、**gs との数値比較は無い** (K64 のソースの読みと解析式・pysad の独立な経路で確かめた)。`SPLIT` は K64 でも obsolete。`SEED` も obsolete。
- v2 モードの台本の読み方 (MAIN の `OFF RFSW` が FFS に届く、`ORG` の負の値、MOMENTUM を設定しない計算の警告、`FFS;` の書き方、MAIN の代入 `NAME = 式;`) は docs/03 の同じ節に書いた。
- gs は `TERMINATE INPUT` の後で入力が尽きて止まるが、pysad は続ける。

**検証** (小さなコマンド): MATRIX・DUMP・DISP・REFERENCE・FIT の UNIPOL/RELW・RENUMBER・NOKEEPEXP は gs の出力と一致 (K64 に無いコマンドと MEASURE 系は除く)。詳細は [history/04_ffs_small_commands_validation.md](history/04_ffs_small_commands_validation.md)。
