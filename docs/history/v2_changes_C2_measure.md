# v2 の変更 (担当 C2): TRCOD FIXCOD / MVCOD / STNCOD / ALIGN、ORBITJITTER、PUTMEASURE / GETMEASURE と MEASURE の残り (2026-11-17)

担当 C の記録 `v2_changes_C_measure.md` の「使えない範囲・未確認」(§5) と「要判断」(§6) の続き。
実装は `pysad/measure.py`、コマンドのつなぎは `pysad/ffscmd.py` (2 か所)、試験は `pysad/tests/test_v2_C2_measure.py` (17 関数・20 件)。
例は `example/v2_C/C2m_trcod_reference.sad`・`C2m_putget_measure.sad` と出力 `*.pysad.txt` (v2 モードの出力)。

**参照について**: K64 の実行ファイル (gs) は MEASURE が無出力、TRCOD が SIGSEGV (B42) なので、ここの項目も gs との比較は 1 件も無い。
動作は K64 のソース (`tffsa.f`、`trcoda.f`、`nalign.f`、`tfojit.f`、`tfgetm.f`、`nfgetm.c`、`tracke.f`、`tdrwdt.f`、`twsdrw.f`、`mbmp.f`、`tftwiss.f`) の読み、
値は解析式と pysad の独立な経路 (CALC の Twiss と軌道、1 粒子の TrackParticles) で確かめた。ヘルプ (`SADHelp.html`) にこれらのコマンドの説明は無い。

すべて **v2 モードだけ**で動く。k64 モードは従来どおり (`TRCOD FIXCOD` などは「未対応 (無視)」、PUTMEASURE・GETMEASURE・ORBITJITTER は語として認識しない)。
k64 モードの経路に入れた変更は無い (`ffscmd.py` の 2 か所は `_compat.v2()` の中、`measure.py` の新しい関数は v2 からだけ呼ばれる)。

## 1. 一覧

| 項目 | K64 の動作 (ソースの読み) | v2 の動作 | 根拠 | 試験 |
|---|---|---|---|---|
| `TRCOD FIXCOD [位置]` | 位置 1..measp の今の Twiss の EX EPX EY EPY と「15..18 番」を `twiss(*,ndim,*)` (基準の軌道) と私用の領域 itw に写し、` trcodf 位置 名前` を表示。位置を書くと measp になる。15..18 は DETR を入れる前の番号で、今は DETR DX DPX DY (**DPY が入らない**)。`twiss(*,ndim,*)` を読むのは旧補正パッケージ (`mbmp.f`) と DRAW のモニターの読み (`tdrwdt.f`) だけで、どちらも K64 の FFS からは呼ばれない → **保存した基準を利用者が読む手段が無い** | EX EPX EY EPY DX DPX DY DPY を、REFERENCE と同じ参照光学 (K64 の `twiss(*,-1,*)`) の 8 関数と私用の写しに書く。`Twiss["RDX", k]` で基準、`Twiss["DDX", k]` で 現在 − 基準 (EX.. も同じ)。`Twiss["DX"]` などは絶対値のまま。他の光学関数の参照は触らない | §3 | `test_fixcod_*` 4 件、`test_stncod_without_fixcod_and_fixcod_after_trcod` |
| `TRCOD MVCOD [位置]` | 1 点 (省略すると measp。measp は変えない) の基準の DX DPX DY DPY を今の値にする (私用の写しがあればそれも) | 同じ | 式どおり | `test_fixcod_position_mvcod_stncod` |
| `TRCOD STNCOD` (`STN`) | 私用の写しを `twiss(*,ndim,*)` へ戻す。範囲は値の定まらない局所変数 meas0。FIXCOD がまだなら何もしない | 私用の写しを参照光学の 8 関数へ戻す。範囲は FIXCOD で保存した範囲。FIXCOD がまだなら注意を出して何もしない。用途: REFERENCE で参照光学を上書きした後に基準を戻す | 式どおり | 同上 |
| `TRCOD ALIGN 要素1 要素2` | 要素1 (最初に合う位置) の DX・DY を、その次から 要素2 (最初に合う位置。見つからなければ終端) までの BEND QUAD SEXT OCT の DX・DY に書く。逆順なら入れ替える。要素1 が磁石でなければ 0 を書く | 同じ。対象に DECA DODECA MULT を足した。分割された磁石の扱いと SOL は §2 | §3 (平行移動との一致) | `test_align_*` 3 件 |
| `ORBITJITTER xj yj` (`ORBITJ`) | 位置 1..nlat の Twiss の DX に xj·tgauss()、DY に yj·tgauss() を足す (各位置で x、y の順)。値が足りなければ `?Missing parameter for ORBITJ_ITTER x y` | 同じ | §3 | `test_orbitjitter_*` 2 件 |
| `PUTMEASURE` (`PUTM`) | 直前の MEASURE / TRCOD の xa ya xxa xya yya を、環境変数 SADPOUT のファイル (FIFO) に `%15.7lg` 5 個のコンマ区切りで追記。標準エラーに `Putting data:…`。開けなければ exit(1) | 同じ形式で追記 (普通のファイルでよい)。開けない・SADPOUT が無いときは理由を表示して続ける | 往復で 7 桁 | `test_putmeasure_getmeasure_*` 2 件 |
| `GETMEASURE [alpha]` (`GET`…) | SADPIN のファイルから 5 個を読み、measp の 1 点の fit 条件のうち関数番号 32 と 34 の目標値を (xa1 − xa)/2·α + 今の値·(1 − α) (y も同じ) にする。32・34 は今の番号では PEPY・PZPX | DX と DY の条件を同じ式で書き換える。`FFS["GETMEASURE"]` は読んだ 5 個を返す。読めなければ理由を表示して目標値は変えない | 式どおり | `test_getmeasure_sets_fit_goals` |
| DETR > 1 の分布 (C の未確認) | `ttinit.f`: DETR > 1 なら x・y のモードを入れ替える | 変更なし (確認だけ) | §3 | `test_ttinit_detr_gt1_consistent_with_calc_coupled_twiss` 2 件 |
| LWAKE / TWAKE ON の MEASURE (C の未確認) | 全粒子 z = 0 の分布 (wake 用の分布はコメントアウト) | 変更なし (z の分布を作る。確認だけ) | §3 | `test_measure_lwake_*` 3 件、`test_measure_twake_uniform_wake_centroid` |
| DRAW の中の粒子追跡 (`tdrwdt.f`) | **K64 でも動かない死んだコード** (§5) | 実装しない | §5 | — |

## 2. 仕様

### TRCOD FIXCOD / MVCOD / STNCOD (基準の軌道)
- 何のためか: TRCOD (追跡の重心) か CALC が Twiss に入れた軌道 DX DPX DY DPY と分散 EX EPX EY EPY を「基準」として取っておき、
  後で誤差やステアリングを変えたときの軌道を基準との差で見る。
- `TRCOD FIXCOD [位置]`: 位置 1 から位置 (省略すると前回の MEASURE / TRCOD / FIXCOD の位置、最初は終端) までを基準にする。位置より後ろの基準は変えない (初めてなら 0)。
- 読み方: `Twiss["RDX", k]` (基準)、`Twiss["DDX", k]` (現在 − 基準)。R・D を付ける関数は EX EPX EY EPY DX DPX DY DPY。`DISP` の `R` / `D` の表でも同じ。
  Python からは `measure.trcod_reference(ffs)` ((nlat, 8) の配列。列は `measure.COD8_KEYS`)。
- **REFERENCE との関係**: K64 の基準は `twiss(*,ndim,*)`、REFERENCE の参照光学は `twiss(*,-1,*)` で別の場所。K64 には前者を読む手段が無いので、
  v2 は同じ参照光学に書く。`REF;` は参照光学を今の光学で全部上書きする (基準も上書きされる)。戻すには `TRCOD STNCOD;`。
  逆に FIXCOD は 8 関数だけを書くので、REFERENCE で取った BX などの参照は残る。
- 次の CALC で Twiss は計算し直されるが、基準は残る (`USE` で格子の長さが変わると、次の FIXCOD で作り直す)。

### TRCOD ALIGN
- `TRCOD ALIGN 要素1 要素2`: 要素1 の最初の位置 is の DX・DY を、is+1 から 要素2 の最初の位置 ie までの磁石に書く (絶対値を写す。設計値との差ではない)。
  同じ架台に載った磁石の並びを、先頭の磁石のずれに揃える使い方。名前はパターンでもよい (最初に合う位置)。
- 分割された磁石 (同じ名前が隙間なく並ぶ): is はその先頭へ、ie はその末尾へ広げる。
- 要素2 を書かない・見つからないとき: 終端までを対象にし、その旨を表示する (K64 は黙って終端まで。K64 と同じく第 2 の語は消費する)。
- 要素1 が見つからないとき: 何もせず `?TRCOD ALIGN: undefined element …` (K64 は黙って戻る)。
- 要素1 が磁石でない (MARK、DRIFT など): 区間の磁石の DX・DY を 0 にする (K64 と同じ)。
- 実行後は光学を無効にする (次の Twiss[] / CALC で計算し直す)。

### ORBITJITTER
- `ORBITJITTER xj yj`: 今の Twiss の DX・DY に、位置ごとに独立な正規乱数 (GCUT で切る。FFS の乱数列を 2·nlat 個進める) を足す。DPX・DPY と要素は変えない。
  モニターの読みの揺らぎを Twiss の軌道に重ねる使い方。次の CALC で元に戻る。

### PUTMEASURE / GETMEASURE
- `PUTMEASURE`: 直前の MEASURE / TRCOD の重心 x、y と <xx> <xy> <yy> (LOGMEASURE と同じ 5 個) を SADPOUT のファイルに 1 行追記する。
  ファイルは最初の PUTMEASURE で開き、開いたままにする (行ごとに flush)。SADPOUT の値を変えると開き直す。
- `GETMEASURE [alpha]`: SADPIN のファイルから次の 1 行を読む (呼ぶたびに次の行)。`FIT 位置; DX …; DY …;` (位置は MEASURE の位置、1 点の条件) があれば
  目標値を DX: (xa1 − xa)/2·alpha + 今の目標·(1 − alpha)、DY: (ya1 − ya)/2·alpha + 今の目標·(1 − alpha) にする (xa、ya は自分の直前の MEASURE の重心、alpha の既定は 1)。
  二つのビームを別々の SAD で追い、互いの重心の差の半分ずつ軌道を寄せる (衝突点の軌道フィードバック) ための道具と読める。
- 環境変数は SADScript から `Environment["SADPOUT", "ファイル名"];` で設定できる。Python からは `measure.close_measure_files(ffs)` で閉じる
  (次の GETMEASURE はファイルの先頭から)。
- **OS による違い**: K64 は名前付きパイプ (FIFO) を想定している (相手の SAD が書くまで読む側が待つ)。Windows には FIFO が無いので、v2 は普通のファイルで同じ形式を
  読み書きする (待たない。データが無ければ `?GETMEASURE: no more data in …`)。Linux・WSL では FIFO を指定すれば K64 と同じ待ち合わせになる (未確認: 試験は普通のファイルだけ)。
- 標準エラーの表示 (` opening pipe:…`、`Putting data:…`、`data gotten from :…`) は K64 と同じ文。

## 3. 根拠 (数値)。すべて v2 モード

| 確かめたこと | 独立な基準 | 結果 |
|---|---|---|
| FIXCOD の基準 (輸送路 14 要素、1 GeV、入口の軌道 1 mm) | 取った時点の Twiss の 8 関数 | 全位置でビット一致。`Twiss["DX"]` は絶対値のまま |
| FIXCOD 後の差 (輸送路。薄いステアリング Δθ = 1e-6) | Δx = −Δθ·√(β1β2)·sin(ψ2−ψ1)、Δx' = −Δθ·√(β1/β2)·(cos − α2 sin) (β α ψ は CALC) | 4 点で相対 1e-5 以内 (例: 終端 Δθ = 1e-4 で DDX = −4.68118e-4、解析値 −4.68118e-4) |
| FIXCOD 後の差 (リング: 24 セルの FODO + CAVI、2 GeV、周長 139.2 m、νx の小数部は共鳴から離れている。Δθ = 1e-5) | 閉軌道の式 −Δθ·√(β1βk)·cos(\|ψk−ψ1\|−πν)/(2 sin πν) | 6 点で相対 1e-4 以内 (例: 位置 1 で 7.364202e-5、式 7.364176e-5) |
| TRCOD の後の FIXCOD (NP = 400、ε = 1e-16、DP = 1e-5) | CALC の軌道・分散 | 基準 − CALC: DX 1.7e-11 m、EX 4e-4 m (統計誤差) |
| ALIGN (QUAD → QUAD SEXT OCT MULT、dx = 0.1 mm、dy = −0.2 mm) | ずれのない格子に入口の軌道 (−dx, −dy) を入れた CALC の軌道 + (dx, dy) | 出口の DX DPX DY DPY が 1e-13 以内で一致 |
| ALIGN (区間に BEND を含む) | 同上 | 4e-10 m の差 (ずれの 2 次 dx²·θ ≈ 5e-10。pysad の BEND の DX の扱いによる。ALIGN の問題ではない) |
| ORBITJITTER (xj = 10 µm、yj = 30 µm) | 同じ種の状態から独立に作った乱数列 × xj、yj | 全位置で 1e-18 以内。DPX・DPY は不変 |
| ORBITJITTER と GCUT = 1.5 (4001 位置) | 切った正規分布の rms = σ·√(1 − 2aφ(a)/erf(a/√2)) = 0.7437σ、最大 1.5σ | rms 4% 以内、最大値は 1.3σ〜1.5σ |
| PUTMEASURE → GETMEASURE | 書いた値 | 7 桁で一致 (`%15.7g`)。2 行書いて 2 回読むと順に戻る |
| GETMEASURE の目標値 (alpha = 0.25 と 1) | tfgetm.f の式 | 相対 1e-12 |
| DETR > 1 の分布 (結合のある輸送路: スキュー四極 3 台。入口 R = 0.1, 0.3, −0.05, 0.2、DETR = 1.5 と 0) | ttinit の線形写像 A から Σ0 = A·Aᵀ、1 粒子の追跡の転送行列 M で M Σ0 Mᵀ。これを CALC がその位置まで運んだ Twiss から作った Σ と比べる | 4 点とも「通常」か「入れ替え」の片方の形で 1e-9 以内 (実際は 1e-11)、もう片方は 0.1 以上ずれる。固有エミッタンスは EMITX・EMITY (相対 1e-10) |
| LWAKE ON の MEASURE (一様な縦ウェイク 1e12 V/C、PBUNCH = 1e10、1 GeV、NP = 600 と 1500) | 重心の Δδ = −e·N·W/(p0c) = −1.602177e-6 | 相対 1e-9。NP に依らない。広がりは不変 |
| TWAKE ON の MEASURE (一様な横ウェイク 1e14 V/C/m、重心 x ≈ 1 mm、y ≈ −0.5 mm) | Δpx = e·N·W1·<x>/(p0c)、Δpy も同じ | 相対 1e-9 (Δpx = 1.5978e-7、Δpy = −8.076e-8) |
| LWAKE ON の MEASURE (階段状の縦ウェイク、σz = 3 mm、NP = 2000) | 平均 −e·N·W/(2 p0c)、δ の分散の増加 (e·N·W/p0c)²/12 | 平均 相対 5e-4 以内、分散の増加 3% 以内 (z の分布が効いている。K64 の z = 0 では増えない) |

## 4. v1 (k64 モード) と v2 の違い (親の決定: v2 の実装のまま)

| 項目 | K64 のソース | gs の実動作 | pysad k64 モード (v1 と同じ) | pysad v2 モード |
|---|---|---|---|---|
| MEASURE と TRPT | `trpt=.true.` (と `trf0=0`、`vcalpha=1`) にして**戻さない**。MEASURE の後の CALC・EMIT・TrackParticles は TRPT ON で走る | MEASURE は何も出力しない (B42) | MEASURE は「未対応 (無視)」。TRPT は変わらない | 追跡の間だけ TRPT ON、終わったら MEASURE の前の値に戻す。後の計算に影響しない。K64 のソースどおりにしたければ MEASURE の後に `TRPT;` と書く |
| TRCOD と TRPT | TRPT を変えない (今の TRPT のまま追う) | SIGSEGV (B42) | 「未対応 (無視)」 | TRPT を変えない (K64 のソースと同じ) |
| MAIN の `TRACK … TURNS=-1` | `trpt=.true.`、`rfsw=.true.`、1 回通過。`ttinit` の分布を作り `Statistics at Entrance` / `Exit` を表示 (輸送路の 1 回通過) | MAIN の TRACK が落ちる | TURNS=-1 を 1 周のリングの追跡として扱う (|TURNS| = 1、RFSW ON。Entrance / Exit の統計は出ない) | K64 のソースどおり輸送路の 1 回通過 (`measure.track_trpt`)。TRPT・RFSW は終わったら戻す |
| `TRCOD FIXCOD` など、PUTMEASURE、GETMEASURE、ORBITJITTER | §1 | (TRCOD は SIGSEGV。他は未確認) | 従来どおり動かない | §1 |

## 5. DRAW の中の粒子追跡 (`tdrwdt.f`): 実装しない理由
- `tdrwdt.f` は、ビームサイズ (kv = 28..32) を描くとき `trsize` が真なら `tracke` の STANDBY / TRACK / SAVE … RESET を呼ぶ。
- `trsize` は K64 では**常に偽**: フラグ名 `TRACKSIZ` は `tffs.f:298` でコメントアウトされ、`tffsa.f:151` で `trsize=.false.` にしている。利用者が ON にする手段が無い。
- 仮に真でも、追跡の結果を使う行 (`call tracke(…'TRACK'…)` と `beam(1)=sv(3)` など、tdrwdt.f:121〜130) はコメントアウトされていて、
  描く値 `beam` は更新されない (前の要素の値のまま)。
- したがって K64 に「意図された動作」として翻訳すべきものが残っていない。DRAW のビームサイズは `tfbeam` (EMIT のビーム行列) の経路だけ。
- 同じ `tdrwdt.f` の kv = 33、34 (モニターの読み = 軌道 − 基準 `twiss(*,ndim,*)` − 四極のずれ − モニターの誤差) は、モニターの表 (`imon`、`nmon`) を作る
  旧補正パッケージ (`mbmp.f` の `mcmon`) が K64 の FFS から呼ばれない (`tffsa.f:499` の `mcmess` はコメントアウト、MON コマンドは無い) ので、これも到達しない。
  添字も `2*kv-51` = 15、17 (DETR を入れる前の DX、DY の番号) のまま。

## 6. 使える範囲と条件
- v2 モード (既定) の FFS。`USE` で格子を選んだ後。輸送路 (INS) とリング (CELL) のどちらでも。
- FIXCOD / MVCOD / STNCOD: Twiss に軌道が入っていること (CALC か TRCOD の後。まだなら FIXCOD が CALC を走らせる)。
- ALIGN: 対象は BEND QUAD SEXT OCT DECA DODECA MULT。DX・DY だけ (回転 ROTATE、縦のずれ DZ は揃えない)。
- PUTMEASURE / GETMEASURE: 環境変数 SADPOUT / SADPIN。GETMEASURE が目標値を書き換えるのは、MEASURE の位置の 1 点の DX・DY の条件だけ。
- MEASURE と LWAKE / TWAKE: `WakeFunction` の表は 3 行以上 (TrackParticles と同じ条件)、TRPT は MEASURE が自分で ON にする。PBUNCH を設定すること。
- 速さ: 上の試験 20 件で 1.5 秒。

## 7. 使えない範囲・未確認
- **gs との比較は無い** (B42)。ここの動作は全部ソースの読み。
- **GETMEASURE の対象の関数 (DX・DY) は判断**: K64 の番号 32・34 が指していた関数名は、今のソースからは確定できない (今の番号では PEPY・PZPX)。
  式が x・y の重心の差の半分なので軌道の関数と判断した。旧い関数表 (`nlist` の末尾に名前だけ残る DDX DDPX DDY DDPY = 基準との差) だった可能性もある。
  その場合の意味は「基準 (FIXCOD) からの軌道の変化の目標」で、v2 の DX の目標とは基準の分だけ違う。
- FIFO での待ち合わせ (二つの SAD を同時に走らせる使い方) は試していない。
- ALIGN の SOL・MONI: K64 の id = 20 の枝 (今の番号では SOL。DX DY を「設計値 − ずれ」にする) は旧い型番号の名残と判断して、v2 は SOL・MONI を動かさない。
- `Twiss["DETR"]` が pysad の CALC の後で常に 0 になる (§8 の 1)。リングの周期解が「入れ替え」の形 (DETR > 1) になる強い結合のリングでは、
  MEASURE の初期分布が「通常」の形として作られる。**強い結合のリングの MEASURE は未確認**。輸送路で先頭 MARK に DETR を書く場合は正しい (§3)。
- 区間に BEND を含む ALIGN は、平行移動との一致が 4e-10 m まで (§3)。
- MEASURE の wake: 確かめたのは一様・階段状のウェイクでの重心と広がり。実際の形のウェイクでの値は TrackParticles の wake (gs と 10 桁一致、`docs/38_wake_guide.md`) に依る。
  K64 の MEASURE (全粒子 z = 0) とは結果が違う。

## 8. 見つけたこと (親が判断)
1. **`Twiss["DETR"]` が常に 0** (k64・v2 とも): `ffscmd.py` の CALC の行を作る箇所 (`rows.append(TwissRow(… r1=…, r2=…, r3=…, r4=…, x_cod=…` の所。1121 行付近) が
   `detr=tw.get('detr', 0.0)` を渡していない (`_transport_twiss` は detr を計算している。小数位置の行を作る 1502 行付近は渡している)。
   gs の `Twiss["DETR"]` は R の行列式 (入れ替えの形では 1 + xyth − 行列式)。共有ファイルの CALC の中心なので直していない。
2. K64 のバグ (ソースの読みのみ。docs/21 は未採番):
   (a) `trcoda.f` の trcodf / trcods が添字 15..18 を直書き (DETR を入れる前の番号。今は DETR DX DPX DY で、DPY が入らない)。trcodm は名前で書いてあり正しい。
   (b) `nalign.f:41` `if (master(is) .lt. 0) then is=-master(is)`: 固定形式では変数 `thenis` への代入になり、何もしない。
   (c) `nalign.f` の id = 20 の枝が今は SOL に当たる。
   (d) `tfgetm.f` の関数番号 32・34 が今は PEPY・PZPX。
   (e) `tffsa.f` の STN_COD が、値を入れていない meas0 を範囲に渡す。
   (f) `tdrwdt.f` の `2*kv-51` が旧い番号 (到達しないコード)。
3. **前任 C の試験の 1 行が古くなった**: `pysad/tests/test_v2_C_measure.py::test_measure_particle_loss_and_errors` の最後の
   `assert '?TRCOD FIXCOD: not supported' in out`。FIXCOD を実装したので、この文は出ない (` trcodf` の行が出る)。
   他人の試験は触らない規則なので直していない。`assert '  trcodf' in out` に替えれば通る (他の 18 件は通る)。
4. FIXCOD の基準を REFERENCE と同じ参照光学に置いたこと (§2)。別の置き場にして専用の読み方を足す方がよければ、`measure._ref_rows` の 1 か所を替える。
5. ALIGN の対象に DECA DODECA MULT を足したこと (K64 は BEND QUAD SEXT OCT だけ)。K64 のとおりにするなら `measure._ALIGN_V2` を `_ALIGN_K64` にする。

## 9. 使い方の例 (`example/v2_C/C2m_*.sad`。`py -3.12 ../../run_sad.py C2m_trcod_reference.sad`)

```
FFS; USE TL; INS; CALC;
TRCOD FIXCOD;                  ! 今の軌道・分散を基準にする (  trcodf 16 $$$ と表示)
ST K0 1E-4; CALC;              ! ステアリングを変える
Print[Twiss["DDX","$$$"]];     ! 基準との差 −4.68118e-4 (解析値 −Δθ·√(β1β2)·sin Δψ = −4.68118e-4)
TRCOD MVCOD MM;                ! MM の 1 点だけ基準を今の値に
REF; TRCOD STNCOD;             ! REFERENCE で上書きされた基準を戻す
NP=500; GAUSS; TRCOD; TRCOD FIXCOD;   ! 追跡の重心を基準にする
ORBITJITTER 1E-5 5E-6;         ! Twiss の DX DY に 10 µm、5 µm の揺らぎ
TRCOD ALIGN QA MM;             ! QA から MM までの磁石の DX DY を QA の値に揃える
```

```
Environment["SADPOUT","C2m_measure.dat"]; Environment["SADPIN","C2m_measure.dat"];
r=FFS["MEASURE IP"]; PUTMEASURE; LOGMEASURE;     ! 相手のビームの測定値を書く
ST K0 5E-5; r=FFS["MEASURE IP"];                 ! 自分のビーム
FIT IP; DX 0; DY 0;
r=FFS["GETMEASURE"];           ! {5.330096e-5, −1.439175e-4, …}。DX の目標は (5.3301e-5 − (−7.8317e-5))/2 = 6.580878e-5
SHOW;
```

`example/v2_C/README.md` に足す行の案 (親がまとめる):
- `C2m_trcod_reference.sad`: TRCOD FIXCOD / MVCOD / STNCOD (基準の軌道と差 `Twiss["DDX"]`)、ORBITJITTER、TRCOD ALIGN。v2 だけ (k64 では「未対応」)。
- `C2m_putget_measure.sad`: PUTMEASURE / GETMEASURE (測定値をファイルで受け渡し、fit の目標値にする) と LWAKE ON の MEASURE。v2 だけ。実行すると `C2m_measure.dat` ができる (消してよい)。

## 10. 変更したファイル
- `pysad/measure.py`: `compat` の import を先頭へ (関数内の局所 import をやめた)。`cmd_trcod` の副コマンドの分岐、`ffs_command` に 3 コマンド。
  新規: `trcod_fix`、`trcod_move`、`trcod_restore`、`trcod_reference`、`cmd_align`、`cmd_orbitjitter`、`cmd_putmeasure`、
  `read_measure_record`、`cmd_getmeasure`、`close_measure_files`。
- `pysad/ffscmd.py` (2 か所、どちらも `_compat.v2()` の中): `is_word` に ORBITJ_ITTER / PUTM_EASURE / GET_MEASURE の省略形、
  `_exec_word` の MEASURE の分岐の条件に同じ 3 語。
- 新規: `pysad/tests/test_v2_C2_measure.py`、`docs/history/v2_changes_C2_measure.md`、`example/v2_C/C2m_trcod_reference.sad`・`.pysad.txt`、
  `example/v2_C/C2m_putget_measure.sad`・`.pysad.txt`。
