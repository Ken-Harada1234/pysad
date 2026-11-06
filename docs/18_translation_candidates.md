# 18. pysad 翻訳候補の調査 (ユーザー全スクリプト × SAD 機能)

最終更新: 2026-11-06 (現状のスナップショット)

目的: `ユーザーの SAD スクリプトのディレクトリ` と `ユーザーの旧スクリプトのディレクトリ` の全スクリプトで実際に使われている FFS コマンド・SADScript 組み込み関数・SAD 標準パッケージ関数を集計し、
pysad に未翻訳のものを「翻訳候補」として難易度・検証可能性つきで並べる。あわせて、ユーザーのスクリプトに現れないが SAD の重要機能で未翻訳のものも載せる。
集計はユーザーのディレクトリの grep・一覧のみで行い (読み取り専用)、スクリプトは実行していない。SAD の実行による参照値取得もしていない。
そのため「ベンチ検証可能性」の列はソース読解と既存ドキュメント (11/12/14, work_*.md) からの推定で、SAD での実走で確かめたものではない (「要確認」と書いたもの)。

## 0. 集計方法と読み方 (必ず先に読むこと)

| 項目 | 内容 |
|---|---|
| 対象ファイル | `sad` と `sad_archive` 配下で拡張子が `.sad` / `.n` / なし のもの。15603 本のうち、バイナリ・3 MB 超 (1 本) を除き 15513 本を読んだ |
| 「スクリプト」の定義 | コメント (`!` 以降、`(* *)`) を除いた後に、FFS コマンド・フラグ・`Name[` 形式の関数呼出しのいずれかを 1 つ以上含むファイル。13583 本。ラティスファイルでも末尾に `CELL; CALC;` 等を持つものは含まれる (ラティスだけのファイルは除外される) |
| 重複 | 13583 本のうち、内容が完全に同一のファイルを 1 つにまとめると 7048 本 (sad_archive の `memoryNN` 系ディレクトリ、modules_*.n の複製が多い)。表では「全/ユニーク」で書く。**頻度の判断はユニーク数を主にする** |
| 内訳 (全) | sad_archive 10858、3GSR 1699、sad のその他 (PF, elegant, RitsSOR, pf3, aurora …) 778、pf2 248。ユニークでは 5238 / 1036 / 583 / 191 |
| 数え方 | 1 ファイルに何回出ても 1 と数える (ファイル数)。FFS コマンドは (a) `FFS["…"]` 文字列内、(b) ベタ書きの行頭語、の文頭語を数える。`;` と改行で文を分ける。ベタ書きは、コマンド語 (下記語彙) の直後が `[ ( = :` でないものに限る |
| 語彙 | コマンド: tffsa.f の `abbrev(word,…)` 一覧 + pysad `_COMMANDS`。フラグ: tffs.f `fname` 53 個 + sino 反対語。FFS の `ON x y z` は x,y,z を各フラグとして数える |
| 組み込み/パッケージ/ユーザー定義の判定 | 組み込み = `src/*.f` の `itfunaloc('名前',…)` 467 個 (tfinitn.f)。パッケージ = `Packages/*.n` 内の定義 (`Name[...]:=`、`Name=Class@…`) 537 個。ユーザー定義 = スクリプト内の `Name[…]:=` / `Name=Function[…]`。Emittance は Packages/emit.n (CalculateEmittance のラッパー)、Element / SetElement は Packages/element.n (Element$ / SetElement$ のラッパー) |
| 誤差 | 正規表現ベースなので ±数 % の誤差。`MOMENTUM`, `NPARA`, `DZ`, `END`, `OUT` などは FFS コマンドではなく変数代入・ラティス語が混ざるため表から外した |
| 除外 | Tcl/Tk・GUI (KBF*, Canvas, Window, Frame, Button, ListPlot, Show, Plot, Graphics, Rectangle, Text, OpticsPlot, OpticsMatchPanel, TkWait ほか)、EPICS (CaGet/CaPut/CaMonitor/CaOpen/CaClose/CaRead/CaWrite)、純粋な IO・文字列 (Open*, Close, Read*, Write*, StringJoin, ToString, DateString …)、ユーザー定義関数本体 (WriteGraphs, DAsv, FitValue 等。ただしそれらが呼ぶ組み込みは対象) |

主なディレクトリの偏り (ファイル数、全): VSX_Ring (CHROMATICITY 572 / TRACK 492 / UNTIL 506)、PF_before_STRUG、Super_SOR、DBA、DR/DRBT、demin (CALC6D/WSPAC/SPAC/SELFCOD/GAUSS)、AR、cERL、3GSR (INTRA 48、CalculateOptics/Monitor/Steer 各 125)。

## 1. FFS コマンドの頻度と pysad の状況 (sad_archive 込み)

pysad の判定は `pysad/ffscmd.py` の dispatch (1000〜1135 行) と `_COMMANDS` / 警告出力の grep による。
「対応」= 12_user_scripts_status.md の判定を流用 (Level 13/14 で SAD と照合済み)。

### 1.1 対応済みのコマンド (候補ではない)

| コマンド | ファイル数 全/ユニーク | pysad |
|---|---|---|
| STOP | 6692 / 3164 | 対応 |
| USE | 5197 / 2820 | 対応 (VISIT 162/107 も同じ経路。`USE BeamLine[..]` も可) |
| DISP / TYPE | 3987 / 2018、3240 / 1617 | 対応 (反転要素名 `-NAME` の差は 12 の §4-7) |
| SAVE | 3295 / 1718 | 対応 |
| FIT / FREE / GO | 2911 / 1528、2901 / 1520、2878 / 1492 | 対応 (既定は解析微分。多解問題の 2 段目以降は SAD と別解のことがある)。VariableRange / VariableWeight も対応 |
| CALC / CAL | 2753 / 1321、154 / 132 | 対応 |
| FIX | 1546 / 767 | 対応 |
| FITP | 1468 / 697 | 対応 (FITP>1 は SAD 側が再起動を繰り返して落ちるので参照値なし) |
| DP / DP0 / GCUT | 1821 / 924、688 / 239、363 / 192 | 対応 |
| COUP / COUPLE | 639 / 292、141 / 59 | 対応 (`COUP` は ffscmd の `_couple`) |
| DELX / DELY / DTHETA / DK / DDK / DL | 834 / 422、747 / 363、731 / 369、786 / 397、221 / 79、47 / 19 | 対応 (terror.f の直訳、照合済み) |
| EXPAND / REJ / RESET / EXIT / ABORT / BYE | 77 / 33、93 / 84、39 / 36、70 / 51、214 / 124、138 / 97 | 対応 |
| TRPT / CODPLOT | 407 / 239、464 / 256 | 対応 (フラグ。CODPLOT は表示なし) |
| DUMP | 411 / 174 | 受理 (ffscmd.py:1465。出力内容の SAD との一致は未確認) |
| EMIT (コマンド/フラグ) | 2800 / 1266 | 対応 |
| 要素値設定 (`SDA K2 …`, `MAX`/`MIN`)・FITFUNCTION・入口 Twiss (`AXI` `BXI` … `EPXI`)・fit 条件 (`NX` `NY` `BX` …) | NX/NY 1029、AXI 283、BXI 282 (FFS 文字列内) | 対応 |
| UNTIL / REPEAT | 1159 / 569、1067 / 501 (VSX_Ring 506、PF_before_STRUG 356、DBA 70) | 対応 (最上位の FFS 文列。候補 C-04) |
| IF / ELSEIF / ELSE / ENDIF | 932 / 632、342、190 / 88、512 / 247 (Super_SOR 320、PF_before_STRUG 297、VSX_Ring 104) | 対応 (同上) |
| CHROMATICITY (CHRO) | 817 / 396 (VSX_Ring 572、PF_before_STRUG 81、DBA 47、pf2 17) | 対応 (候補 C-03) |
| TRACK (MAIN の `TRACK USE= NPART= TURNS= NX= NY= NZ=`) | 727 / 174 (再集計: 725 ファイル、TRACK 文 1727。VSX_Ring 492、PF_before_STRUG 145) | 対応 (`pysad/maintrack.py`。SAD K64 の MAIN の TRACK は SIGSEGV で動かず、SADScript で同じ手順を書いた版と一致)。PHSPLOTS の描画は無視、FOURIE ON は未対応 (docs/history/work_main_track.md、候補 C-14) |

### 1.2 ON/OFF のフラグ (全 / ユニーク)

| フラグ | ファイル数 | pysad | 備考 |
|---|---|---|---|
| LOG / COD / CELL / RADCOD / RFSW / FLUC / RAD / RING / INS | 6849 / 3249、6809 / 3246、6598 / 3031、6181 / 2900、6055 / 3012、5866 / 2786、5848 / 2912、4475 / 2692、2372 / 1199 | 対応 | INS の初期軌道オフセットも一致 |
| DAPERT | 1872 / 844 | 対応 (dapert.py、DA の地図一致) | 実寸 DAsv の最終 seed はカオスで不一致 (Score 1〜3 %) |
| CTIME / ECHO / BIPOL / EMIOUT / GAUSS / UNIPOL / CONV | 854 / 470、826 / 468、553 / 305、71 / 37、59 / 57、19 / 13、23 / 10 | フラグとして受理。**効果の有無は未確認** (GAUSS = 初期分布、BIPOL = 磁場の両極性、EMIOUT = emit 出力) | 低優先 |
| **INTRA** | **95 / 57** (3GSR 48、VSX_Ring 14、Super_SOR 13、pf2 4) | 対応 (INTRA で EMIT が IBS 付き、TouschekLifetime の表) | C-05 |
| CALC6D | 56 / 49 (demin 55) | 対応 (CELL と INS) | U-01 |
| WSPAC / SPAC / SELFCOD | 56 / 48、33 / 33、42 / 38 (すべて demin) | 受理のみ。空間電荷 (tspac.f) 未翻訳 | 低優先 |

### 1.3 未対応・一部のコマンド (候補)

| コマンド | ファイル数 全/ユニーク | 主な出所 | pysad | 候補 ID |
|---|---|---|---|---|
| MON / MCOD / TCOD | 868 / 372、384 / 195、316 / 169 | 3GSR 約 100、archive | 受理のみ (K64 ビルドでは SAD も何もしない) | 対象外 |
| PRINT | 295 / 71 | archive | 未対応 (警告して無視) | C-15 |
| REF / RENUM / ORG / SHOW / STAT / BEAMSIZE / MATRIX / GEO | 112 / 57、118 / 56、69 / 60、80 / 43、22 / 17、34 / 22、26 / 20、31 / 7 | PF (REF 103)、VSX_Ring (RENUM 63)、DR (ORG 20) | 未対応 (警告して無視)。SAD 側は RENUM = 要素番号の振り直し (HelpMessages.n:332)、MATRIX = 転送行列の表示 (tfmat.f 106 行)、ORG = 座標原点 (tffsa.f:779) | C-15 |
| DRAW | 601 / 410 | 全般 | 警告して無視 | 除外 (描画) |

## 2. SADScript 関数の頻度と pysad の状況

「組」は 組み込み (B) / パッケージ (P) / ラッパー (W: パッケージ経由で組み込みを呼ぶ)。ファイル数は 全 / ユニーク。

### 2.1 光学・追跡・ラティスに関わる関数 (翻訳対象)

| 関数 | 組 | ファイル数 | 出所 | pysad | 候補 ID |
|---|---|---|---|---|---|
| Twiss | B | 4781 / 2182 | 3GSR+pf2 940 | 対応。キー NX NY BX BY EX AX AY DX DY EY EPX DPX S `*` は対応。**SIGX/SIGY/SIGPX/SIGPY/SIGZ/SIGDP/SIGE/SX/SY (計 約 15 ファイル、ビームサイズ系) は未対応**。GAMMABETA/GAMMA も要確認 | C-09 の一部 |
| LINE | B | 4295 / 1958 | 3GSR+pf2 758 | 対応。キー NAME K1 LENGTH K2 POSITION K0 L ANGLE S ELEMENT DIR DX DY TYPE TYPENAME DK VOLT ROTATE E1 E2 は対応または値あり。GX GY GZ (66 ファイル)、CHI1-3 (28)、OFFSET (57) は対応。**RHO/RHOX/RHOY/KX/SK0 は要確認** | C-12 |
| FFS | B | 4259 / 1904 | | 対応 | |
| Emittance | W (emit.n) | 650 / 419 | | 対応。オプション `OneTurnInformation->True` (約 100)、`Matrix->True` (約 36)、`OneTurnTransferMatrix` (7)、`SaveEMIT->True` (54)、`Output->6` (17) は emit_extra.py で対応 | C-09 |
| DynamicApertureSurvey | P (dapert.n) | 1681 / 868 | | 対応 (遅い) | |
| TrackParticles | B | 436 / 302 | | 対応 | |
| CalculateOptics | B (tfoptics.f) | 574 / 296 | 3GSR+pf2 143 | 対応 (4 次元の cell/INS・区間・結合・分散、idim=3 は全周) | C-02 |
| **Monitor / Steer** | P (corutility.n) | 491 / 259、618 / 324 | 3GSR 125 | **未対応** (関数として未登録。FFS の MON とは別物) | C-06 |
| ExtractBeamLine | B (tfbeamline.f) | 448 / 281 | 3GSR 59、cERL ほか | 対応 (`USE BeamLine[..]` / `VISIT` / `BYE` も) | C-07 |
| BeamLine / PrintBeamLine / WriteBeamLine / ExpandBeamLine | P (beamline.n) / B | 104 / 94、12 / 11、7、1 | | 対応 | C-07 |
| SetElement | W (element.n → SetElement$) | 508 / 326 | 3GSR 56 | 対応 (`SetElement["名前", {"K2"->…}]` の値設定、`SetElement[名前, "MULT"/"bend", {…}]` の型指定・新規定義とも。PROFILE, Saved->True, DK の連動を除く) | C-08 |
| Element | W (element.n → Element$) | 540 / 285 | archive 509 | 対応 (`Element[キー, 名前]` のスカラー値、`Element["KEYWORDS"]` 10、`"DEFAULT"` 6、`"VALUE"` 61、`"POSITION"` 11 などの形とも。PROFILE, Saved->True, DK の連動を除く) | C-08 |
| FLAG | B (pgflag.f) | 256 / 124 | archive 239 | 対応 (`FLAG["CELL"]` が大半) | C-10 |
| GetMAIN | B (tffs.f tfmain) | 172 / 162 | archive 168 | 対応 (MAIN 形式ラティスを読み込む。スクリプト中の絶対パスは /users/… で古い) | C-10 |
| FFSFlags | B | 8 | archive | 登録名あり。内容は未対応 | C-10 |
| CorrectOrbit | P (correction.n) | 155 / 92 | archive 146 | 未対応 | C-11 |
| MakeBump | P (bump.n) | 66 / 21 | archive 58 | 未対応 | C-11 |
| OptimizeOptics | P (optimize.n) | 49 / 47 | archive | 未対応 (Amoeba/DownhillSimplex 10 を含む SADScript) | C-06 |
| TransferMatrix / BeamMatrix / NormalCoordinates / SymplecticInverse | P (Optics.n, emit.n) | 57 / 52、8、7、39 / 34 | archive のみ | 対応 (optics_fn.py、SAD と 1e-9) | |
| PolynomialFit | P (Fit.n) | 394 / 171 | 3GSR+pf2 20 | 対応 (Fit.n をそのまま評価) | |
| Fit / Spline / NIntegrate / Integrate / FindRoot / D | P (Fit.n, D.n) | 47 / 32、48 / 40 (Spline は対応)、56 / 47、28 / 24、18 / 15、13 / 12 | archive | Spline のみ対応。他は未対応 (D.n は SADScript 実装) | C-06 |
| TouschekLifetime | P (emit.n) | 12 / 11 | 3GSR 9 | 対応 (15 ケース ≤7e-12) | C-05 |
| SeedRandom / GaussRandom / Random | B | 640 / 325、298 / 160、152 / 117 | | 対応 (SAD と同じ乱数列) | |
| PolynomialFit, CalculateBeta 等のユーザー定義 | — | — | — | 対象外 | |

### 2.2 SADScript 言語系の組み込み (頻度)

頻度は全 / ユニーク。`pysad/sadlist.py` に MapThread, ScanThread, MapAt, Complement, DeleteCases, Cases, MatchQ, Extract, Intersection, Union, Position, Count, MemberQ, FreeQ, Select, Sort, Take/Drop, Insert/Delete/ReplacePart, Apply, Map, Fold, Nest, FixedPoint, Through ほかを実装済み (候補 C-01)。Total/Cross/Norm/Riffle/StringSplit は SAD にも無い (未評価のまま返る)。表のそれ以外の関数の対応状況はこの文書では追っていない。

| 関数 | ファイル数 | 備考 |
|---|---|---|
| MapThread | 150 / 139 | 3GSR+pf2 57 |
| ScanThread | 66 / 59 | |
| MapAt | 54 / 50 | |
| Complement | 51 / 49 | |
| DeleteCases | 45 / 41 | |
| Cases / SelectCases / SwitchCases | 25 / 24、17 / 17 | パターンマッチ (Cases は MatchQ と同じ基盤) |
| MatchQ | 18 / 18 | |
| Extract | 23 / 23 | |
| Depth / Restrict / MinMax / Second / Order / Intersection | 23、23、13、13、10、12 | |
| Operation / StandardForm / StringFill / FromCharacterCode / ToCharacterCode | 76 / 33、38 / 28、26 / 17、33 / 33、9 | |
| Sleep | 159 / 57 | no-op でよい |
| Fourier / InverseFourier / BesselJ / BesselK | 10 / 7、18 / 11、10、8 | tfinitn の実装は FFT / Bessel |
| ReleaseHold / MemoryCheck | 14 / 14 | |

### 2.3 除外した関数 (集計のみ。翻訳しない)

GUI/描画 (ユニーク数): KBMainFrame 160 + KBF* 数十、Canvas 260、OpticsPlot 276、Show 220、Graphics 156、Rectangle 121、Window 129、ListPlot、Plot、TkWait 479 (no-op 済み)、TextLabel、DrawBeamLine2D 30 など。
EPICS: CaGet 104、CaPut 60、CaMonitor 32、CaOpen/CaClose/CaRead/CaWrite。
IO/文字列: Close 904、OpenWrite 891、OpenAppend 1031、Read 542、OpenRead 465、Write 1253、ToString 258、StringPosition 128 など (すべて対応済み)。
ユーザー定義: WriteGraphs 1450 (全)、DAsv 994、FitValue 963、RMSWrite、ChromCorrect、AmplitudeDependentTuneShift、MatchingRingTune 等 (modules_*.n は 12 §2.2 で整理済み)。

## 3. 翻訳候補表 (ユーザーのスクリプトで使われているもの)

難易度: 低 = SADScript の組み込みや単純なコマンド。中 = 既存モジュールの拡張・接続。高 = 新しい物理計算の直訳。
ベンチ検証: 取れる = SAD の WSL ビルド (`bin/gs`) で参照値が取れると判断。要確認 = 見込みはあるが SAD での実走は未確認。不可 = SAD が落ちる/無効/GUI 依存。
ファイル数は 全/ユニーク。「pysad の現状」列に対応状況を書く。

| ID | 候補 | 使用 | pysad の現状 | SAD 側の実装箇所 | 難易度と根拠 | ベンチ検証 |
|---|---|---|---|---|---|---|
| C-01 | SADScript 言語系の欠落組み込み (MapThread, ScanThread, MapAt, Complement, DeleteCases, Cases, MatchQ, Extract, Intersection, Depth, SwitchCases, StringFill, Sleep, MinMax …、§2.2) | MapThread 139、ScanThread 59、Complement 49、DeleteCases 41、MapAt 50 ほか。3GSR+pf2 でも MapThread 57 | **対応済み**: `pysad/sadlist.py` (MapThread, ScanThread, MapAt, Complement, DeleteCases, Cases, MatchQ, Extract, Intersection, Union, Position, Count, MemberQ, FreeQ, Select, Sort, Take/Drop, Insert/Delete/ReplacePart, Apply, Map, Fold, Nest, FixedPoint, Through, …)。Total/Cross/Norm/Riffle/StringSplit は SAD の実測で **SAD にも無い** (未評価のまま) | tfinitn.f の `itfunaloc` 登録と tfefun*.f / tfeeval.f のディスパッチ。パターン系は `Pattern`/`Blank` が sadscript.py にある | 低: 既存インタプリタ (`B[...]`) に関数を足すだけ。Cases/MatchQ/DeleteCases はパターンマッチの再利用が要るので低〜中 | 取れる (SADScript 単体でリスト演算を Print して比較。ラティス不要) |
| C-02 | CalculateOptics[lb, le, inicond(28), cell, idim, geo] | 574 / 296 (3GSR+pf2 143) | **対応済み** (4 次元の cell/INS・区間・結合・分散、idim=3 は全周) | `src/tfoptics.f` (117 行): 初期条件 28 個を行 lb に書き込み、qcell1 を実行して 28 関数の表と安定性・トレースを返す。idim=2 → CALC6D off、3 → on | 低〜中: 本体は qcell1 の薄いラッパー。pysad は lattice の INS 計算 (`ins_twiss`) があるので 2 次元版は直結できる。idim=3 は calc6d.py | 取れる (`Print[CalculateOptics[1,LINE["LENGTH"],Twiss["*",1],0,2]]` を SAD と比較) |
| C-03 | FFS コマンド CHROMATICITY (CHRO) | 817 / 396 (VSX_Ring 572 ほか) | **対応済み** | `tffsa.f:609` → `tfchro.f` (51 行): QUAD / SEXT / MULT ごとの薄肉近似クロマティシティ (X, Y) と合計の表を出力 | 低: 閉形式の式 51 行。α, β, η の表だけで計算できる | 取れる (表を diff。薄肉式なので精度は ≤1e-10 が期待できる) |
| C-04 | FFS の REPEAT/UNTIL、IF/ELSEIF/ELSE/ENDIF | UNTIL 569、REPEAT 501、IF 632、ENDIF 247 (ユニーク) | **対応済み** (最上位の FFS 文列) | `tffsa.f:255-290` (UNTIL/REPEAT。入力位置 `ipoint` を保存して巻き戻す。`REPEAT n` 回数 / `UNTIL 式`)、`tfif.f` (95 行) | 低〜中: 式評価は `_getva` がある。難所は「FFS の入力ストリームを巻き戻す」ことで、pysad が文を 1 行ずつ実行する構造 (`run_text` / 行リスト) に依存する | 取れる (FFS 単体で `repeat 3; … until` を使うテストを SAD と比較) |
| C-05 | `ON INTRA` (IBS 付き emit / Touschek) の ffscmd 接続と TouschekLifetime | INTRA 57 ユニーク (3GSR 48、VSX_Ring 14、Super_SOR 13、pf2 4)。TouschekLifetime 11 | **対応済み** (INTRA で EMIT が touschek_table を呼び、TouschekLifetime の表を返す) | `temit.f` の INTRA 経路 (tintra 系)、`Packages/emit.n:129-` の TouschekLifetime / TouschekTable | 中: 計算部はあり、接続 (EMIT 時に touschek_table を呼ぶ、Emittance[] に IBS 後の値を載せる、表を SADScript に渡す) だけ | 取れる (Level 17A の 15 ケースを流用。history/work_touschek.md に接続案) |
| C-06 | Packages/*.n の自動ロード基盤と、それで動くパッケージ関数 (Monitor, Steer, CorrectOrbit 以外の bump/optimize/D/Optics/beamline) | Monitor 259、Steer 324 (ユニーク)、OptimizeOptics 47、TransferMatrix 52、NIntegrate 47、Integrate 24、Fit 32、FindRoot 15、DownhillSimplex 10 | 未対応。PolynomialFit だけ Fit.n を直接評価。汎用の AutoLoad/Require はない | `Packages/Library.n` (AutoLoad, Require)、`init.n` の AutoLoad 表、corutility.n (376 行、Monitor/Steer)、optimize.n、D.n、Optics.n。Monitor/Steer は FFS の MON とは別で、`Master[]` (pgmast.f 83 行) と `LINE["type"/"element",…]` を呼ぶだけ | 低〜中: 自前翻訳せず SAD のパッケージ (SADScript) をそのまま sadscript.py に読ませる方針。必要なのは `Master[]`、小文字キー (`LINE["element"]`)、シングルクォート文字列、`Module`/`With` の細部。12 §4-5 の単一引用符 (`'…'`) と `Twiss['*',k]` が先に要る | Monitor/Steer は取れる (戻り値のリストを Print)。FFS の `MON` 自体は SAD で無効 (12 §3-7) |
| C-07 | ExtractBeamLine / BeamLine / ExpandBeamLine / WriteBeamLine / PrintBeamLine | ExtractBeamLine 281 (3GSR 59)、BeamLine 94、PrintBeamLine 11 | **対応済み** (`USE BeamLine[..]` / `VISIT` / `BYE` も) | `tfbeamline.f:362-510` (ExtractBeamLine は `BeamLine[名前…]` を返す。逆向き要素は `-1*名前`)、`beamline.n` (BeamLine の演算子定義、WriteBeamLine) | 低〜中: lattice.py の `els` (名前列) から BeamLine を作り、beamline.n の規則 (入れ子展開・負符号反転) を SADScript で評価する。`RING = ExtractBeamLine["RING"]` は USE 先の再定義として使われる | 取れる (`Print[ExtractBeamLine["RING"]]` と WriteBeamLine の出力) |
| C-08 | Element / SetElement の全形 | Element 285、SetElement 326 (3GSR 56) | **対応済み** (PROFILE, Saved->True, DK の連動を除く) | `tfbeamline.f:10-148` (Element$: ELEMENT/POSITION/KEYWORDS/DEFAULT/VALUE)、`149-340` (SetElement$: 型指定・新規定義・キー設定)、`Packages/element.n` | 中: 型指定での要素の新規作成・型変更は parser/lattice の要素表を動的に変える必要がある (`Do[SetElement[s//i,"MULT",{"L"->…}]]`)。キーごとの値設定は `set_family` で足りる | 取れる (SetElement 後に `LINE["K1",…]`/Twiss を比較) |
| C-09 | Emittance[] のオプション接続 (OneTurnInformation, Matrix, OneTurnTransferMatrix, SaveEMIT, Output, Orbit, InitialOrbit, InitialBeamMatrix) と Twiss[SIGX…] | OneTurnInformation 約 100、Matrix 約 36、SaveEMIT 54、Output->6 17、Twiss SIGX 系 15 | **対応済み** (emit_extra.py、2.1e-10。Twiss["SIGX"..] も)。偏極・SpinTune は 0 固定 | `tfemit.f` (226 行)、`temit.f` (mode 1〜3)、`Packages/emit.n:20-100` | 中: 接続のみ。`fn_emittance` の戻りに規則リストを足し、`Twiss["SIGX"]` を beam_size_at に結ぶ | 取れる (Level 21A の 2.1e-10 を流用) |
| C-10 | FLAG[] / GetMAIN[] / FFSFlags[] | FLAG 124、GetMAIN 162 (ほとんど古い絶対パス)、FFSFlags 8 | **FLAG / GetMAIN は対応済み**。FFSFlags は未 | `pgflag.f` (104 行)、`tffs.f` の tfmain (GetMAIN)、`tfetok.f` の tfflags (FFSFlags) | 低: FLAG は `self.flags` を引くだけ。GetMAIN は Get と同じ MAIN 形式ラティスの読込み (parser.py) | FLAG は取れる。GetMAIN は絶対パスのため参照を取りにくい (ファイルを手元にコピーすれば取れる) |
| C-11 | ResponseMatrix + LinearSolveConditioned (+ CorrectOrbit, MakeBump, correction.n/bump.n) | CorrectOrbit 92、MakeBump 21、ResponseMatrix は correction.n/bump.n 経由 | 未対応 | `pgrmat.f` (503 行、解析的な COD 応答行列。Twiss から作る)、`pgsolvcond.f` (251 行、条件付き SVD 解)、correction.n (490 行)、bump.n (370 行) | 中〜高: 物理は解析式で素直だが 750 行ある。パッケージ (SADScript) 側は C-06 の基盤で動く | 取れる見込み (SAD の `CorrectOrbit` はユーザー側で MON を使わない形の例が要る。要確認) |
| C-12 | LINE["GX"/"GY"/"GZ"/"CHI1".."CHI3"] と fit の GX (幾何) | GX/GY 66、GZ 28、CHI 28 | **対応済み** (geo.py、4.8e-13。fit の GX/GY/GZ/CHI は数値微分、tdgeo は未翻訳) | `tfgeo.f` (411 行)、`tdgeo.f` | 中: LINE キーの接続は低。`FIT GX` は微分 (tdgeo) が必要で中〜高 | 取れる (Level 19 のケース流用) |
| C-13 | 解析微分 (match6d.py) を FIT/GO の既定にする | fit/go 約 1500 ユニーク、FITP 697 | **対応済み** (match6d.py が Matcher の既定。esrf matching を 12 桁再現) | `qdcell.f` の tffsqu、`tffsmatch.f` | 中: 接続と既定化。FITP>1 のオフモメンタム fit は微分の対応が別途必要 | FITP=1 は取れる。FITP>1 は **SAD が落ちる (不可)** |
| C-14 | MAIN の TRACK コマンド (`TRACK USE= NPART= TURNS= NX NY NZ PHSPLOTS`) | 727 / 174 (VSX_Ring 492、PF_before_STRUG 145) | **対応済み** (`pysad/maintrack.py`、docs/history/work_main_track.md)。SAD K64 の MAIN の EMIT/TRACK は動かない (EMIT は何もせず TRACK は tinip で SIGSEGV) ので、SADScript で同じ手順を書いた版 (bench_level27a/b) と照合。PHSPLOTS の描画は無視、FOURIE ON・ResultsOfTracking フックは未対応 | `track.f`, `tracka.f`, `ttinit.f` (初期分布 NX NY NZ はσ単位)、TopDrawer 出力 | 中: 追跡コアは既存。初期分布生成と出力形式が必要。`STO M A`/`RECALL` などの旧サブコマンドもある | 描画 (PHSPLOTS) が出力なので、粒子座標の dump 比較に限り可。基本は **GUI/描画依存 (不可)** |
| C-15 | 小コマンド: PRINT、REF、RENUM、ORG、MATRIX、BEAMSIZE、GEO、SHOW、STAT、BEAM | 各 17〜71 ユニーク | 未対応 (警告して無視) | `tffsa.f` 内の各分岐、tfmat.f (106 行)、tffsa.f の tfsetref | 低: 表示系が主。RENUM は要素番号の再割当て (要素参照 `NAME.n` に影響) | 取れる (表示 diff) |
| C-16 | BEND の K0 キック (steering) が COD/分散に効かない | pf2 の bump/disp_response 系 26 本が応答 0 (work_pf2_survey)。LINE["K0"] 994 ファイル | **対応済み** (CALC の tsteee、追跡の tsteer。bump3_r/calc_bump_r/kick_r が SAD と一致) | `tsteer.f`、`tbende.f` の K0 キック、`tbend.f` | 中: 14 の優先度「高」。複合の bendib とは別に K0≠0 の経路が要る | 取れる (esrf/pf2 の `DELX`+`K0` で比較) |

表外: 2 段目のフィット (重み・tsvd・INS の解析微分・フック) は対応済み。VariableRange / VariableWeight は対応済み。

## 4. ユーザー未使用だが SAD の重要機能で未翻訳のもの

「ユーザー未使用」= `sad`/`sad_archive` の全スクリプトで 0 ファイル (または極少数)。MULT/SOL/beam-beam/wake/CSR は低優先。

| ID | 機能 | ユーザー使用 | pysad の現状 | SAD 側の実装箇所 | 難易度と根拠 | ベンチ検証 |
|---|---|---|---|---|---|---|
| U-01 | CALC6D (INS を含む 6D 光学計算) の ffscmd 接続 | CALC6D フラグ 49 (demin のみ。ユーザーの主要スクリプトでは未使用) | **対応済み** (CELL と INS、1.7e-10) | `qcell.f` の qcell6d、`tfetwiss` ほか | 中: 接続のみ。INS 初期条件付きは未照合 | 取れる (Level 18A/B/C) |
| U-02 | POL / 偏極 (SPIN、Polarization, SpinTune, PolarizationTime) | 0 | 未翻訳 (emit_extra は 0 固定) | `temit.f:491,1188` の calpol、`polpar.f` (111 行)、`photons.f`・`tturn.f` の spin 部、`tftrack.f` (calpol で mcf=9) | 高: スピン行列・放射偏極の新規計算。追跡 (mcf) も 9 次元になる | 要確認 (K64 ビルドで POL が動くか) |
| U-03 | SynchroBetaEmittance / SYNCHROB | 0 | 語を受理するだけ。未翻訳 | `tfemit.f:144` (tfemits)、`temits.f` (1325 行)、`emit.n:91` | 高: シンクロ・ベータ結合の放射減衰 (大きなルーチン) | 要確認 |
| U-04 | RADINT コマンド (放射積分 I1〜I5 の表) | 0 (ユーザーは calc_RadInt.sad で Twiss[] から自前計算) | **対応済み**: `radint.py` (INTGRL.f の直訳) を FFS RADINT に接続。ESRF 6 通り・小リング 3 通りが SAD と 7 桁一致 | `INTGRL.f` (693 行) | 低: SAD の表形式で出力する | 取れる (出力 diff) |
| U-05 | TrackParticles の PHOTONS (放射光子リスト `PhotonList`) | PHOTONS フラグ 0 | 未翻訳 | `tftrack.f:144,315` (photons 分岐)、`photons.f` の kphlist | 中: 追跡中の放射光子を記録するだけ。乱数列は一致済み | 取れる見込み |
| U-06 | OneTurnInformation の出力詳細 (`OneTurnTransferMatrix`, `Orbit`) | 約 100 (C-09 に含む) | C-09 | | | |
| U-07 | Apert 系 (APERTURE コマンド、ZAP) | APERTURE 4、ZAP 0 | tapert (追跡側の損失判定) は翻訳済み。コマンドは未対応 | `tfaprt.f` (151 行)、`tfzap.f` (16 行) | 低: 表示系 | 取れる |
| U-08 | LMA (局所運動量アクセプタンス) | calc_LMA.sad は Twiss[] の表示だけ | SAD に LMA の組み込みはない。DA の追跡 (DynamicAperture) の変種で計算するのが SAD の使い方 | 該当なし | — | — |
| U-09 | TouschekTable の IBS 部 (intra の積分) | INTRA 57 → C-05 | touschek.py に実装済み (接続のみ) | `temit.f` intra | C-05 | C-05 |
| U-10 | Element 系の全キーワード (`Element["KEYWORDS"]`, `"DEFAULT"`, `"VALUE"`, ElementValues) | 10 / 6 / 61 ファイル | C-08 (対応済み) | `tfbeamline.f`、`tffsa.f:1632` ElementValues | 中 | 取れる |
| U-11 | SaveBeamLine / WriteBeamLine の FFS/MAIN 書出し | 7 (WriteBeamLine) | C-07 (対応済み) | `beamline.n` | C-07 に含む | 取れる |
| U-12 | OpticsCalc の Transport (TRPT 付き光学計算、INS で Twiss 伝搬) | TRPT 239、INS 1199 | INS: 照合済み (≤2.6e-7)。TRPT の Emittance 無視は未確認 | `qcell.f` | 低: 確認のみ | 取れる |
| U-13 | tfgeo の SOL / MAP (ソレノイド内の軌道追跡で幾何を決める) | 0 | geo.py は SOL と MAP を直進扱い (注意を出す) | `tsgeo.f`、`tgmap` | 中 | 取れる見込み |
| U-14 | Monitor/Steer 系の COD 補正関数 (SADScript 側) | Monitor 259、Steer 324 (C-06、C-11) | C-06/C-11 | corutility.n、correction.n、bump.n | C-06/C-11 | C-06/C-11 |
| U-15 | MULT ANGLE≠0 (tmulta.f 433 行)、emit の SOL (tsole.f/tsolque.f)、TCAVI (ttcav.f 72 行)、COORD (tcoord.f 82 行)、MAP (temap.f 625 行)、UND (undulator.f 165 行) | 0〜数本 (MULT は PF-II の ID で pf2 に少数) | 未翻訳 (NotImplementedError) | 左記 | 中〜高 | 取れる (低優先) |
| U-16 | beam-beam (bbstrhl.f)、wake (tffswake.f, txwake.f)、CSR (csroy.f ほか)、BBBrem (tbbbrem.f)、空間電荷 (tspac.f)、ルミノシティ (tlum.f) | WSPAC/SPAC/SELFCOD は demin のみ 33〜56 | 未翻訳 | 左記 | 高 | 低優先 (衝突型・FEL 用) |
| U-17 | DELCOR / tcorr.f / tvcorr.f (相関付き誤差・補正) | 0 (FFS コマンド DELCOR) | 未翻訳 | `tcorr.f` 151 行、`tvcorr.f`、tffsa.f:982 | 中 | 要確認 |
| U-18 | FFS の MEASURE / PUTMEASURE / GETMEASURE / LOGMEASURE / ORBITJITTER / TRCOD / FIXCOD / MVCOD / STNCOD / LTRACK / PKSEED / PUSHSEED / POPSEED / EXCGSEED / DELWAVE | 0 | 受理のみ (TRCOD, MEASURE) または未定義 | `tffsa.f:639-985`、`trcoda.f` (158 行)、`trackb.f` | 中: 測定と補正の組 (MON を使うので K64 では参照が取れない) | 不可に近い (MON 無効) |
| U-19 | OneTurnInformation 以外の SADScript パッケージ: LDSolve, ExpMatrix, rational.n, extmap.n (ExternalMap), Misc.n | 0 | 未翻訳 | Packages/*.n | 低 (C-06 の基盤で動く) | 取れる |

## 5. 推奨順位 (頻度 × 物理的重要度 ÷ 難易度)

頻度はユニークファイル数 (全体)、物理的重要度は「ユーザーの最終目的 (esrf/pf2 の設計スクリプトを SAD と同じ数値で再現)」への寄与で判断する。「対応状況」列は §3 の現状。

| 順位 | 候補 | 頻度 | 難易度 | 順位の理由 | 対応状況 |
|---|---|---|---|---|---|
| 1 | C-01 SADScript 言語系の欠落組み込み (MapThread, ScanThread, MapAt, Complement, DeleteCases, Cases, MatchQ, Extract, Intersection …) | 最大 139 (MapThread) | 低 | 名前未登録だとスクリプトが途中で止まる。ほぼ機械的に足せる。ラティス不要で検証できる | 対応済み (sadlist.py) |
| 2 | C-03 FFS CHROMATICITY | 396 (全 817) | 低 | 51 行の閉形式。old スクリプトの大半 (VSX_Ring) が使う。単独で検証できる | 対応済み |
| 3 | C-04 REPEAT/UNTIL、IF/ELSEIF/ELSE/ENDIF | 約 570 / 630 | 低〜中 | 1 つ欠けるとスクリプトの制御フローが壊れる。SAD の入力ストリーム処理 (tffsa.f 255-290) の直訳 | 対応済み (最上位の FFS 文列) |
| 4 | C-02 CalculateOptics | 296 (3GSR+pf2 143) | 低〜中 | 117 行のラッパー。COD 補正系 (correct_cod*) の前提で、Monitor/Steer と組で使われる | 対応済み |
| 5 | C-05 ON INTRA 接続 (Touschek/IBS) | 57 (3GSR 48) | 中 | 3GSR の touscheck 系に直結 | 対応済み |
| 6 | C-13 解析微分の既定化 | fit/go 約 1500、FITP 697 | 中 | 最大の頻度。FITP>1 は SAD が落ちるので照合は FITP=1 に限る | 対応済み |
| 7 | C-06 Packages ローダ (+ `Master[]`、単一引用符、小文字キー) → Monitor/Steer/OptimizeOptics/NIntegrate/Fit/TransferMatrix | Monitor/Steer 約 300 | 低〜中 | 自前翻訳せず SAD の SADScript を読む。1 つの基盤で多数の関数が解決する | 未対応 |
| 8 | C-16 BEND の K0 キック | pf2 26 本、LINE["K0"] 994 | 中 | bump/応答系の前提 | 対応済み |
| 9 | C-09 Emittance[] オプション + Twiss[SIGX…] の接続 | 約 100 | 中 | | 対応済み |
| 10 | C-07 ExtractBeamLine / BeamLine | 281 | 低〜中 | `RING = ExtractBeamLine["RING"]` が最初の行で使われる | 対応済み |
| 11 | C-08 Element / SetElement の全形 | 285 / 326 | 中 | SetElement の型指定・新規定義 (3GSR 56) が必要なスクリプトで効く | 対応済み (PROFILE, Saved->True, DK の連動を除く) |
| 12 | C-10 FLAG / GetMAIN / FFSFlags | 124 / 162 | 低 | FLAG は 1 行。GetMAIN は参照パスの問題が別にある | FLAG / GetMAIN 対応済み、FFSFlags 未 |
| 13 | C-12 LINE["GX"] と geo の接続 | 66 | 中 | FIT GX の微分は数値微分 (tdgeo は未翻訳) | 対応済み |
| 14 | C-11 ResponseMatrix + CorrectOrbit/MakeBump | 92 / 21 | 中〜高 | 物理は素直だが 750 行。archive の COD 補正・バンプ系 | 未対応 |
| 15 | C-14 MAIN TRACK、C-15 小コマンド (PRINT/REF/RENUM/ORG/MATRIX …)、U-01 CALC6D 接続 | 174 / 17〜71 / 49 | 中 / 低 / 中 | 頻度は中程度。TRACK は出力が描画中心なので数値の検証が弱い | C-14 対応済み (maintrack.py)、C-15 未対応、U-01 対応済み |

## 6. 補足

- **K64 SAD の性質**: `MON` / `MCOD` / `TCOD` は SAD でも無効。`FITP>1` の `GO` は SAD が再起動を繰り返す (11/12)。これらは「SAD 側で落ちる/無効」で参照値が取れない。
- **K64 でも使える見込みで未確認のもの**: POL (U-02)、SYNCHROB (U-03)、DELCOR (U-17)。いずれもユーザー未使用なので、着手時に小さい SAD スクリプトで動作を確認してから決める。
- **パッケージ経由のものの注意**: Emittance / Element / SetElement / Monitor / Steer / DynamicApertureSurvey / TouschekLifetime / PolynomialFit は SAD では SADScript パッケージ (`Packages/*.n`) の関数で、Fortran の組み込みではない。pysad が Python で個別に書き直している場合は、SAD のパッケージの既定値・オプションとの差が出やすい (例: DynamicApertureSurvey の既定 Output / ExpandElementValues は ffscmd.py:880 で再現)。
- **sad_archive の古い MAIN コマンド**: `TRACK … PHSPLOTS`、`STO M A` / `RECALL M A`、`SET WINDOW X … Y …`、`OUT 'a' draw …`、`cor x y`、`exec c` は旧 SAD の対話型コマンド。pysad の対象は FFS と SADScript で、MAIN の TRACK だけ maintrack.py で手順を再現している。
- **集計スクリプト**: scratchpad の `survey.py` (ファイル走査)、`agg.py` / `dedup.py` / `keys.py` (集計)。再集計時は `survey.json` を作り直す。
