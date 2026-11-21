# 21. SAD K64 のバグ報告 (pysad 移植中に見つかったもの)

最終更新: 2026-11-21 (現状のスナップショット。B01〜B162 の 162 件。B72〜B110 は v2.0.0 の作業、B111〜B134 は v2.1.0 の作業、B135〜B159 は v2.2.0 の作業、B160・B161・B162 は v2.3.0 の作業 (T1・U1・T2) で見つけたもの (B160 は U1 で「誤りではなく K1 の暗黙の定義」に再分類、B161 は U2 で v2 で直した、B162 は Touschek の損失率の式の適用範囲で v2 は Piwinski 型の式を既定にした)。B40・B42・B60・B81・B83・B99・B106〜B109 は v2.1.0 で、B42・B108・B111・B112・B116・B118〜B120・B128・B129 は v2.2.0 で更新)

> 対象: SAD K64 (github KatsOide/SAD、WSL Ubuntu でビルドした `bin/gs`)。ソースは `<SAD>/src`、`Packages/`。
> 材料: docs/11_sad_quirks.md、docs/history/work_*.md、docs/history/00_RESUME.md。行番号は 2026-10-13〜16 の照合時点のソース。
> 「仕様・癖」(E1/E2 が ANGLE 比、K1 が積分値、`FIT;` の意味、NPARA の既定と fork、LINE[] のファミリー展開など) は含めない。
> 明らかな誤り (綴り・引数の並び・未初期化・クラッシュ・同種ルーチン間の不整合) だけを選んだ。
> 「再現手順」の `bench_levelXX.sad` は `<repo>\` にある。「pysad での扱い」の「再現」は SAD に数値を合わせるため誤りもそのまま写したもの、「正しい動作」は pysad が正しい挙動にしたもの。
> **v2.0.0 から**: pysad は互換モードを持つ。既定の `v2` は K64 のバグを直した動作、`PYSAD_COMPAT=k64` (または `pysad.compat.set_mode('k64')`) は v1 と同じ K64 (gs) の動作。
> 表の「pysad」の欄は、v2 で直したものを「v2 で修正 (k64 モードは再現)」の形で書く。「再現」と書いたものは k64 モードの動作で、v2 では直してあるものがある (欄の文を参照)。
> 一覧は [`47_v2_changes.md`](47_v2_changes.md)。B72〜B110 は v2.0.0 の作業の記録 (`history/v2_changes_*.md`)、B111〜B134 は v2.1.0、B135〜B159 は v2.2.0 の作業の記録 (`history/v2next_*.md`) の「新しく見つけたもの」を採番したもの。B160・B161 は v2.3.0 の T1・U1 (`history/v23_T1.md` §6、`history/v23_U1.md`)。
> pysad 自身の不具合は K64 のバグではないのでここには入れない (47 の (d))。
> 確証の度合い: 「確定」はソースの読みと SAD の実行の両方で確認、「ソース」はソースの読みのみ、「実測」は実行のみ (原因未調査)。

## 0. 一覧 (重要度順)

重要度: 高 = 結果が変わる・落ちる、中 = 特定条件で結果が変わる、低 = 表示・性能。

| ID | 題名 | 場所 | 影響 | 重要度 | 確度 | pysad |
|---|---|---|---|---|---|---|
| B01 | MAIN の `TRACK` が SIGSEGV (tinip で r 未設定、MAIN の `EMIT` は何もしない) | trackdlib の tinip、ActTra.f | 落ちる | 高 | 確定 | 古典の意味で実装 |
| B02 | ANGLE=0 かつ K1≠0 かつ K0≠0 の BEND で FFS 初期化が SIGSEGV | tbend 系 (原因行は未調査) | 落ちる | 高 | 確定 (実行) | 照合不可 |
| B03 | `SetElement` を含むスクリプトの出力をファイルへリダイレクトすると SIGSEGV (pty なら動く) | 未調査 | 落ちる | 高 | 確定 (実行) | pty で参照を取得 |
| B04 | `draw` で core dump | 未調査 (Tk/描画) | 落ちる | 高 | 確定 (実行) | 描画は対象外 |
| B05 | FITP>1 の `GO` で「Run time Environment」から再起動を繰り返す | 未調査 | 落ちる (最終解が取れない) | 高 | 確定 (実行) | 未照合 |
| B06 | RADINT: 全 BEND に「最後に処理した BEND」の L/ANGLE を使う | INTGRL.f (tbendebody の bendib) | I1..I10 などが変わる | 高 | 確定 | v2 で修正 (その BEND 自身の L/ANGLE。k64 モードは再現: `stale_rho`) |
| B07 | RADINT: 6 極ブロックが直前の BEND/QUAD の DL を使う (L=0 の6 極は DL=0.2) | INTGRL.f TCAL | I9, I10 が変わる | 高 | 確定 | v2 で修正 (6 極は自分の L と K2。L=0 は寄与なし。k64 モードは再現) |
| B08 | tsolqur (RAD ON, bz≠0) が最後の粒子の bzp を全粒子に使う | tsolqu.f tsolqur :268-273 | 粒子の並びで結果が変わる | 高 | 確定 | v2 で修正 (粒子ごとの BZ/(1+δ)。k64 モードは再現) |
| B09 | tchge の入口で tbrote 等が回転行列を上書きし、一周行列に回転が入らない | tchge.f :45, :64-68 | 結合・チューンが物理とずれる | 高 | 確定 | v2 で修正 (入口の回転を行列に掛ける。k64 モードは再現) |
| B10 | tdusrn が NaN の anp で npmax=10000 回ループし乱数を大量に消費 | drndsr.f :71-111 | 損失粒子の周以後の乱数列がずれる | 高 | 確定 | v2 で修正 (anp が NaN・±Inf なら光子を出さず乱数は 1 個。k64 モードは 10000 回のループを再現) |
| B11 | 4D の CALC が縦方向の Twiss 列 (AZ BZ NZ ZX ZPX ZY ZPY) を作り直さない | tffscalc.f :87-90、qtwiss1 | Twiss[] に古い値が残る | 中 | 確定 | v2 で修正 (縦の列を既定値 AZ=0・BZ=1・他 0 に戻す。k64 モードは再現) |
| B12 | RADINT: 定数が単精度 (`1e27`, `1.E-6`, `0.2`, `2.3026`) で ELOSS が相対 1.15e-8 小さい | INTGRL.f :35, :172, :212, :266 | 数値が 8 桁で変わる | 中 | 確定 | v2 で修正 (倍精度。k64 モードは再現) |
| B13 | RADINT: CAVI の VOLT/HARM/FREQ を合計せず最後のもので上書き | INTGRL.f | 空洞が複数あると誤り | 中 | 確定 | v2 で修正 (VOLT は総和、FREQ・HARM は電圧のある最初の空洞。k64 モードは再現) |
| B14 | tsolqum の z の dv 項が tsolqu と違う (dv に 0.25 が掛かる) | tsolqu.f :86-88 と :250 | δ=1e-3 で 1 周 5e-10 | 中 | 確定 (実測 6.7e-10) | v2 で修正 (tsolqu と同じ −dv·aln0。k64 モードは再現) |
| B15 | `UnstableLevel` が綴りの誤り (文字数 14) で常に 0 | tffscalc.f :58 | 不安定レベルが読めない | 中 | 確定 | v2 で修正 (不安定な点の数 nstab。k64 モードは 0) |
| B16 | CalculateOptics の stabz と idim=3 の stab/trace が未初期化 | tfoptics.f → qcell1 の `optstat` | 戻り値が不定 | 中 | 確定 (実行) | v2 で修正 (idim=3 は tracez = 2cos(μz) と stabz、idim=2 は stabz=1・tracez=0。k64 モードは 1 と 0) |
| B17 | gfortran の MAX/MIN が NaN を無視する前提のコード (tradkfn、asinz、fuzz) | photons.f、tbend.f、tffsmatch.f :255-289 | NaN の扱いが箇所ごとに違う | 中 | 確定 | 再現 (乱数の消費だけ v2 で NaN を明示的に扱う: B10。asinz・gfloor・sqrt1 は変更なし) |
| B18 | PolynomialFit が平均 0 (xa=0) のとき xl=1 になり次数補正が壊れる | Packages/Fit.n :37-56 | 係数が誤る | 中 | 確定 | v2 で修正 (平均 0 のデータの係数。k64 モードは K64 の値。平均が 0 でないデータは両モード同じ) |
| B19 | qdmdia の label 1300 で dr4 が未設定 | qdmdia.f | 縮退分岐の微分が不定 | 中 | ソース | dr4=0 にした |
| B20 | `Level[{1,{2,{3}}},{-1}]` が `{1,3}` を返す | SADScript の Level | リスト操作の結果が誤り | 低 | 実測 | 再現せず |
| B21 | tcod の rfsw=.false. 枝の `trs(:,5:6)=0` と実測 (z を閉じる δ) が合わない | tcod.f | RFSW/RADCOD OFF の δ | 低〜中 | 実測 (原因未調査) | 近似 (v2 でも変更なし。z が閉じる δ を 1e-15 まで解く) |
| B22 | MON/MCOD/TCOD が何もせず、`MON 0 0 0 0 PM*` が "Missing value for PM*" | FFS コマンド表 (K64 の tffsa.f に MON・MCOD・TCOD のコマンドは無い。ソースで確認 2026-11-16) | 補正スクリプトが止まる | 低〜中 | 実測 + ソース | v2 は「K64 に無いコマンド」と警告して文全体を読み飛ばす (定義はしない。k64 モードは K64 と同じ) |
| B23 | Extended Twiss の EPX の欄に EY が出る | temit.f :1010-1017 | 表示のみ | 低 | 確定 | 再現 |
| B24 | `LINE["LENGTH","$$$"]` が narg エラー | LINE[] の引数処理 | 位置が取れない | 低 | 実測 | v2 でも値は返さず、エラーに位置の取り方 (`LINE["POSITION",loc]`・`LINE["LENG",loc]`・`LINE["S",loc]`) の案内を足した (k64 モードは再現) |
| B25 | TrackParticles の POL ON 列順がソース内で自己矛盾 | tftrack.f :344 と :502 | 文書・メッセージの不整合 | 低 | ソース | 実コードの列順 (コメント側) で対応 (POL ON の TrackParticles は 9〜11 列) |
| B26 | `ORG` の引数が式として読まれ `ORG B3 1.5 -2.5 ...` が `1.5-2.5` になる | tffsa.f ORG、getvad | 負数の書き方で値が変わる | 低 | 確定 | v2 で修正 (空白で区切った `-2.5` は別の値。k64 モードは K64 と同じ) |
| B27 | SeedRandom[] の表示と seed の末尾桁の丸め | 表示ルーチン | 表示のみ | 低 | 実測 (調査浅い) | 符号なし表示に修正 |
| B28 | 未初期化変数の値が出力に出る (EMIT 報告の Closed orbit Entrance、MAIN TRACK の l_track) | temit.f、tracka.f | 表示のゴミ | 低 | 実測 | 再現せず |
| B29 | QUAD に RAD=1 を書くと COD・エミッタンスが 0 になる | tturn.f:360、tsol.f:677、tturne.f:572 (`krad = (RAD_QUAD == 0) .and. al /= 0`) | 結果が変わる。**バグではなく仕様** (RAD=1 はその要素の放射を止める。tturn.f:360 ほか) | 仕様 | 確定 (ソース) | 仕様どおり (DampingRate などで切り分ける) |
| B30 | MAIN の `OFF RFSW` (`$RFSW$`) が FFS の RFSW に届かない (他の RAD/FLUC/COD/RADCOD/INTRA は届く) | track.f :107, :149 (`rfsw=.true.` で tffs を呼ぶ) | `OFF RFSW` を書いても空洞が効く。FFS 入口の表の RFSW 列は MAIN の値で FFS の値と違う | 中 | 確定 (実測 + ソース) | v2 で修正 (MAIN の `OFF RFSW`・`ON RFSW` が FFS に届く。FFS の中の RFSW・NORFSW が後から優先。k64 モードは再現: `main_flags`) |
| B34 | BEAMBEAM の傾き回転 (tilt) が cost=sin(傾き角)、sint=傾き角を使う (bbinit は xangle(6)=sin、(7)=cos を書く) | bbstrhl.f:beambeam (sint=colb%xangle(5)、cost=colb%xangle(6)) | 強ビームに x-y 結合があると座標がほぼ 0 に潰れる | 高 | 確定 (gs で再現) | v2 で修正 ((cos t, sin t) の回転。k64 モードは再現) |
| B35 | TRPT で SOL 領域 (SOL BOUND=1 の間) に DRIFT/BEND/QUAD/SOL/MULT/CAVI/MAP/APERT/MARK/MON 以外の要素 (SEXT、COORD、BEAMBEAM) があると SIGSEGV (本来は 'Only DRIFT, ... are supported within SOL' を出して無視する分岐) | tsol.f:798 (case default 付近) | 落ちる | 中 | 確定 (デバッグビルドで特定、メッセージの write を消しても落ちる。原因は未特定) | v2 は `SolRegionError` (K64 の tsol.f が出すはずの文)。k64 モードは NotImplementedError |
| B36 | L≠0 の TCAVI を含むリングの Emittance[] で EnergyLossU0 が NaN、縦チューン 0、エミッタンスが桁違い (K0=0 でも。ドリフトなら正常) | ttcave.f (tdrife0 に krad を渡す第 2 ドリフト付近。RAD=1 でも NaN なので原因未特定) | 結果が壊れる | 中 | 確定 (gs で再現)、原因は未特定 | 変更なし (pysad は自然な値。L=0 は gs と一致) |
| B37 | UND に磁場のキーワードが無い (FBX FBY BX BY はどれも Undefined keyword)。undulator は磁場なしのドリフトになる | initb1.f:1208-1209 (`kytbl(kwFBX,icUND)=…`。kwFBX が未定義の変数で 0 になり、kytbl(0,icUND) を上書きするだけ) | 機能が使えない | 低 | 確定 (gs で再現) | v2 で `FBX`・`FBY` [T] を受ける (追跡・CALC・Emittance[]・EMIT。RAD の放射も)。k64 モードは SAD のまま (FBX 等は NotImplementedError、光学は恒等) |
| B38 | PHSROT 要素は動作しない: phsinit の行列 (PtoN など) がどんな入力でもすべて 0 になり、追跡で粒子の座標が 0 になる。さらに PX PY PZ (位相) のキーワードは名前表に未登録で指定できず、`i_p0` が未定義の整数変数 (=0) でその dv の更新は blist(1)・blist(2) を使う | phsrot.f (phsinit、phsrot)、initb1.f (kwPX 等) | 機能が使えない | 低 | 確定 (gs で再現、行列の原因は未特定) | v2 で意図された写像を実装 (追跡・CALC・Emittance[]。PX PY PZ も指定できる)。k64 モードは追跡が NotImplementedError、Twiss・Emittance[] は恒等 |
| B39 | tradke (ソレノイド内、bzh≠0) の dbeam(8) (px・py の共分散) が `beam(2)-beam(10)+b*beam(4)` と書かれている。力学的運動量 (px+b·y, py−b·x) の共分散なら `beam(2)-beam(9)+b*beam(4)` (B12 − B34 + b B13) のはずで、beam(10)=B44 は誤記に見える | photons.f tradke | SOL 領域内の放射のビーム行列 (px-py 結合) | 低 (bzh=0 では効かない) | 疑い (gs と合わせるため pysad はそのまま写した。beam(9) に直しても Emittance[] の ε は強い場の試験 (BZ=5、軌道 2.7 mm) でも全桁同じで、試験では区別できない) | v2 で修正 (B12 − B34 + b·B13。k64 モードは K64 の式: `perfkern.tradke_db`) |
| B40 | DELCOR (tcorr.f) が MARK に `latt(j)+15`、`+17` を書く。値の添字では ky_DETR_Mark=15、ky_DPX_Mark=17 (DX=16、DY=18) で、DX・DY の一つ手前を書いている。CAVI も `le=latt+3` に書いて元の値は `idvalc+13` から読む食い違いがある | tcorr.f | MARK への DELCOR は初期条件の DETR・DPX を変える (DX・DY は変わらない)。CAVI は不定 | 低 | 確定 (MARK は gs で再現、CAVI は読みのみ)。CAVI の `tcorr.f` の `id=20` (SOL)・`id=31` (CAVI) の枝は `dx1 = (設計の値) − delx` で、他の磁石 (+delx) と符号が逆。SOL の枝は境界の SOL の組に同じ値を書く (tsol は入口で x −= dx、出口で x += dx なので −delx は領域全体をずらす符号と読める) ので、CAVI の枝はその写しとみられる。しかし K64 の CAVI は追跡で DX・DY を使わない (wake の `txwake` に dx = dy = 0 を渡す) ので、**符号は観測できない** (`LINE["DX","RFC"]` の読みだけ) | v2 で修正 (MARK・CAVI も DX・DY に書く。符号は他の磁石と同じ +delx のまま変更なし: 実害なし)。k64 モードは再現: MARK は DETR・DPX の欄、CAVI は HARM・PHI |
| B41 | ttinit.f (MEASURE の初期分布) が `rlist(kx+26..29)` を軌道の揺らぎ (JDX..JDPY) として読むが、MARK の値の添字では 26..29 は ky_ZPX, ky_ZY, ky_ZPY, ky_DP。JDX.. は 34..37 | ttinit.f | MEASURE の初期の py の揺らぎが DP·tgauss になり、乱数も 4 回余分に消費する | 低 | 確定 (読解。MEASURE が出力しないため gs では未確認) | v2 で MEASURE を実装し、JDX JDPX JDY JDPY を読む (ソースの読みのみ。gs では確認できない)。k64 モードは MEASURE を無視 (`measure.ttinit(b41=True)` で K64 の読み方) |
| B42 | `TRCOD` が gs (この環境の K64 ビルド) で SIGSEGV。`MEASURE` は何も出力しない | trcoda.f / trackb.f | TRCOD・MEASURE の参照が取れない | 低 | 確定 (gs で再現) | v2 で MEASURE・TRCOD・LOGMEASURE・MAIN の `TRACK TURNS=-1` を実装 (ソースの読み + 解析式・pysad の独立な経路で確認。gs との比較は無い)。k64 モードは未対応 (無視) **v2.2.0 の更新**: gs の 10-19 ビルドは MEASURE を出力する (先頭の要素 = 入口の MARK の位置だけ。それ以外は無出力、TRCOD は SIGSEGV)。k64 は先頭の位置の MEASURE・TRCOD・LOGMEASURE を gs と同じ統計 (4 桁) で実装した (B151・B152・B154)。先頭以外は無出力・警告。v2 の追跡を通した値は gs と比べられない |
| B43 | `Hypergeometric2F1[.5,.5,1,.99]` が 2.348414 (真値 2.352716)。`BesselJ[.5,30.]` が -.143924709693 (真値 -.143929653370)。級数・漸近展開の打ち切りが粗い | gamma.f / hg.f (hgrr, cbesj = confhg0) | 引数が 1 に近い 2F1、x が大きい Bessel で相対 1e-3〜3e-5 の誤差 | 低 | 確定 (gs で再現) | pysad は scipy の値 (真値) を返す (以前から。モードによらない) |
| B44 | `ExpIntegralEi[1.5]` が複素数 `.164-2.985 I` (真値 3.3013)、`Ei[-1.5]` が -2.0089 (真値 -.1000)、`Ei[.3]` が `-.166-1.417 I` (真値 -.3027)。D.n の式が HypergeometricU の負引数 (NaN 経路) に依存 | Packages/D.n:ExpIntegralEi | Ei が全域で誤り。LogIntegral は別式で正しい | 中 | 確定 (gs で再現) | pysad は真値 (実数。モードによらない) |
| B45 | `HermiteH[3,.7]` が NaN (整数 n で Gamma[(1-n)/2] が極)。`LegendreP[3,1,.4]`・`LegendreP[2,2,.3]` が NaN (Hypergeometric2F1Regularized の 1-m が非正整数)。`Binomial[-3,2]` が NaN | Packages/D.n | 整数次数で NaN | 中 | 確定 (gs で再現) | pysad は真値 (モードによらない)。LegendreP の奇数 m の符号の誤り (pysad 自身) も v2 で直した |
| B46 | `HypergeometricU[1,2,.5]` が INF (真値 2)、`BesselY[0,.1]` が NaN (真値 -1.5342)、`BesselK[1.5,.7]` が INF。`BesselJ[-1,1.5]`・`BesselY[1,2.5]`・`BesselJZero[0,1]` は gs が応答しない (無限ループ) | gamma.f (chgu, cbesy, cbesk, confhg0) | 整数次数 (b が整数の U、cbesy の n 整数での 0/0) | 中 | 確定 (gs で再現。応答しないのは 20 秒で打ち切り) | pysad は真値 (モードによらない) |
| B47 | `InverseErf` は Newton の係数に sqrt(pi/2) を使い (正しくは sqrt(pi)/2)、`df^2 <= no*f1^2*eps` で止めるため、真値から相対 1e-8 ずれる (`InverseErf[.3]` = .272462711811、真値 .272462714727)。FitGaussian の出力 (σ の 8 桁目) に出る | gamma.f:inverseerf | 8 桁以降の誤差 | 低 | 確定 (gs で再現) | v2 は真値 (`InverseErf[.3]` = .2724627147)。k64 モードは K64 の値 (`specfn.inverseerf`) |
| B48 | `FixedPointList` は gs で未評価のまま返る (autoload の表に `FixedPoinstList` と綴り違い (Misc.n の SetAttributes と init.n)) | Packages/init.n:597、Misc.n:23 | FixedPointList が使えない | 低 | 確定 (gs で再現) | pysad は Misc.n の定義どおり動かす (以前から。モードによらない) |
| B49 | `RadiationField` が gs (この K64 ビルド) で SIGSEGV (どの引数でも。tlcalc の `loc_sad(ktfaddr(ka)-3,dl)` が不正な番地を指すとみられる)。RadiationSpectrum は動く | tlinit.f:tlcalc | RadiationField の参照が取れない | 低 | 確定 (gs で再現) | pysad はソースどおりに翻訳 (`radfield.py`)。等速運動の電荷の速度場と相対 6e-3 で一致を確認 |
| B50 | `FileByteCount` は `wc -l` (行数) を返す (14 バイトの 2 行のファイルで 2)。FileSystem.n のファイル操作 (CopyFile RenameFile DeleteFile …) は gs のこの環境 (WSL の /mnt/c) で失敗が多い | Packages/FileSystem.n | バイト数が取れない | 低 | 確定 (gs で再現) | pysad は Python の os・shutil で動かし、FileByteCount はバイト数 (`fileops.py`) |
| B51 | `RenumberElement[name,{n1,..}]` が gs (この K64 ビルド) で SIGSEGV (tfrenumber_element)。参照が取れない | tfrenumber_element.f | RenumberElement が使えない | 低 | 確定 (gs で再現) | pysad はソースどおりに翻訳 (`ffscmd.fn_renumber_element`。NAME.n の番号の付け替え、足りない分は最大から続ける) |
| B52 | `Fit[data,fun,{x1,x2},...]` (ベクトルの独立変数。data の行は {{x1,x2},y}) が、式を `Function[{v},fun]@@x$[[k]]` で呼ぶ際に General::narg を出し、最小二乗解にならない値を返す (例: ChiSquare 1.1519 に対し真の最小は .9320) | Packages/Fit.n | ベクトル独立変数の Fit が使えない | 低 | 確定 (gs で再現) | pysad は正しい最小二乗解 (モードによらない) |
| B53 | `BesselJ`・`BesselK`・`BesselY` が x が大きいと不正確 (x=33 で BesselJ[0,33.] が 3e-3 相対、BesselK[0,12.5] が 1e-5)。`FresnelS[6.1]` = 1.3e7、`FresnelC[6.1]` = -2.1e7、`AiryAi[15.5]` = -16 など、引数が大きい Fresnel・Airy は桁落ちでデタラメ (超幾何級数の直接和) | gamma.f (cbesj = confhg0 の級数)、D.n (FresnelS/C は HypergeometricPFQ) | 引数が x≳6 (Fresnel・Airy)、x≳12 (Bessel) の値が使えない | 中 | 確定 (gs で再現) | pysad は scipy の値 (真値。モードによらない) |
| B54 | `ExpIntegralEi[x]` は x<0 でも誤り (Ei[-0.1] = -2.2265、真値 -1.8229。B44 の続き)。`PolyLog[n,z]` は gs が複素数 (虚部 1e-15 の雑音つき) を返す | Packages/D.n | | 低 | 確定 (gs で再現) | pysad は真値 (実数。モードによらない)。非整数次の PolyLog の大きい |z| の精度も v2 で改善 |
| B55 | POL の平衡偏極度 `Polarization` / `Polarizations` が NaN になる (試した 8 つのリング、平面でも非平面でも)。`spdepol` → `spintrm` の `exp(cb**2-cphi)*cerfc(cb)` が \|cb\|>18 (cezlim) で漸近式 `exp(-cb²)` → 0 と `exp(cb²)` → ∞ の積になる。cb = −(e^{−cφ}−1)/a、a = √2·ssprd(スピン軌道結合 ~0.05 など) で、cb は通常 20 以上。erfcx 相当 (`erfc·exp` を桁あふれしない形) に書き直すと有限 (4e-7〜8e-7、つまり 4e-5〜8e-5 %) になるが、これも 92% 級の理論値からは遠く、spdepol の式か単位の誤りの疑いがある | photons.f:spintrm (cexp1・cerfc)、gamma.f:cerfc | 偏極度が使えない (PolarizationTime=1/ΣT 式は有限で正しい) | 中 | NaN は確定 (gs で再現)。原因は推定 (pysad で同じ式を写すと同じ NaN、erfcx 版で有限)。有限値の妥当性は未確認 | v2 で修正 (exp(cb²)·erfc(cb) を Faddeeva 関数で。平衡偏極度は線形のスピン軌道結合の脱偏極率から求め、平面リングで 92.376%)。k64 モードは K64 の式のまま (NaN も写す) |
| B56 | `tbrote` (tchge.f:197-202) の srot 更新は `do concurrent (i=1:9)` の中で行を逐次更新するので、2 行目は更新後の 1 行目を、3 行目は更新後の 1・2 行目を使う (rr·srot の行列積にならない)。`serot` (photons.f:850-857) も `c1` を更新してから `c2` の更新に使う。`tradke` の末尾 (photons.f:1614-1615) の `srot(3,:)=(-sp*srot(1,:)+srot(3,:))/cp` も更新後の 1 行目を使う (これは偏向角 sp≠0 の bend で効く。平面リングの gs と一致したので意図かもしれない) | tchge.f:197、photons.f:850、1614 | 偏向磁石のずれ (DX・DZ・DROTATE・CHI2) で POL の結果がずれる可能性 | 低 | ソースの読み (gs での確認は serot と tradke 末尾だけ。tbrote は未確認) | v2 で修正 ((I+rr)·srot。k64 モードは再現) |
| B57 | `tradke` の `call tadddsg(dsg,srot,dgx,dgy,dgx,dbeam)` は dgz の引数に dgx を渡す (photons.f:1589)。dsg は `spequ` (使われない旧ルーチン) でしか読まれないので結果に影響はない | photons.f:1589 | なし | 低 | ソースの読み | pysad は dsg を持たない |
| B58 | `tbendebody`・`tbendebody0` の `bsir0` の更新で `bsi2` の項は enarad が偽の最後のスライスにしか付かず、常に捨てられる (`bsi2=1` は n==n2 のときだけで、そのとき `enarad .and. n /= n2` は偽)。意図は偏向磁石の出口の bsi の補正とみられる | tbende.f:215・313、535-537 | 縦のオフセットがある偏向磁石の POL で出口側の bsi が入らない | 低 | ソースの読み | v2 で修正 (入口で +(y/ρ + k x y)、出口で −。スピン軌道結合にも入れる。追跡の tbendi も同じ)。k64 モードは K64 と同じ (出口の項を捨てる) |
| B59 | `CSRSetupOY` (tfcsroysetup) は dpjn を zjn と同じ記憶域に割り当てる (`c_f_pointer(c_loc(rlist(kazjn)),dpjn,…)` が 2 回、dpjn 用の `kadpjn` を使わない)。`zjn(k,i)=z` の直後に `dpjn(k,i)=p` が同じ場所を上書きするので、返す 4 番目の表 (zj) は z ではなく dp の値になり、5 番目 (dpj) は書かれない (0 のまま)。続く `CSROYMatrix` は zj を位置として使うので行列が狂う。`csrfindc` は見つけた値を返した後の `csrfindc=0.d0` で常に 0 を返す (未使用)。`csrjasign` は `rho` がしきい値以上の点が無いと `i1` が未初期化 | csroy.f:tfcsroysetup (kazjn, dpjn)、csrfindc | CSROY の行列が使えない | 中 | ソースの読みと gs の出力 (dpj の和が 0、zj の値が p のスケール) | v2 で修正 (zj は位置、dpj は運動量。点が無ければ失敗)。k64 モードは K64 と同じ (zj に p を入れ、dpj は 0) |
| B60 | `CSROYMatrix[{wi,J,omega,zj,dpj,rhoj},{dz,inten,nus,sigz,sige}]` が gs (この K64 ビルド) で、`CSRSetupOY` の出力を渡すと General::wrongtype を返す (手で作った正しい形の実数リストなら gs でも動く: v2.0.0 の記述「手で作った実数リストでも」は誤り)。原因は D2 の調査で特定: `tfcsroysetup` が 4 番目 (zj) の実数リストへ `c_f_pointer(c_loc(rlist(kazjn)), zjn, [mphi,nj])` で書くとき、先頭の語 `rlist(kazjn)` (リストの頭 = head) から書き始める (要素は `rlist(kazjn+1)` から書くべきで、ほかの 4 つの配列は `rlist(k+1:…)` で正しい)。head が zjn(1,1) で上書きされ、`tfcsroymat` の `tfreallistq(las%dbody(4), lazj)` が偽になる (B59 の `dpjn` が `zjn` と同じ記憶域も同じ箇所)。登録名は `CSRSetupOY`・`CSROYMatrix` (tfinitn.f:1044-1045。ソース内の関数名 tfcsroysetup・tfcsroymat と綴りが違う) | csroy.f:tfcsroysetup・tfcsroymat | `CSRSetupOY` の出力から `CSROYMatrix` を作れない。`OYStability`・`OYScan` は gs で完走しない | 低 | 確定 (gs で再現し、`CSRSetupOY` の出力の各要素の `Head` を調べて原因を特定: zj の Head が実数 1.2952814481588496e-5) | pysad の `CSROYMatrix` は両モードで動く (`CSRSetupOY` の出力も受ける: gs の wrongtype は再現しない。直していない: 要判断)。手で作った入力 (nr=64、nj=3、nphi=2) の 6x6 の行列 (36 要素) は gs と相対 4e-20 で一致 |
| B61 | `FindRoot` に複素の初期値を与えると壊れた結果になる (`FindRoot[z^2+1==0,{z,1.+1.*I}]` が `{z->Restrict[(.25+I),-Infinity,Infinity],Residual->0}`、`Exp[z]==2` も同様) | tffindroot.f:tfsetupvars (352-435。初期値は実数だけ受け付ける作りで、:414-418 が複素に `General::wrongval` を出す経路。`Restrict[…]` が返る理由は未調査) | 複素の根が求まらない (実数の初期値なら正常: `FindRoot[x^2-2==0,{x,1.}]` は gs と pysad が 1e-14 で一致) | 低 | 確定 (gs で再現、2026-11-17。原因の特定はソースの読みのみ) | v2 は複素の初期値を複素の Newton 法で解く。k64 モードは従来どおり TypeError (gs の壊れた値は再現しない。実数の初期値は両モードとも gs と一致) |
| B62 | `DELCOR` (tvcorr) の corlength≠0 の枝は `cv(n)` を設定せず、`ktaloc` の古い内容が入る (2 回目以降の DELCOR で値が数 % 変わる)。CAVI への DELCOR は DX・DY ではなく HARM・PHI に書く (tcorr.f の `latt+3`、添字のずれ) | tvcorr.f:22-29 (`cv(2..n/2)`・`cv(n/2+1..n-1)`・`cv(1)` だけ設定、`cv(n)` なし)。tcorr.f:83-102 (CAVI の `le=latt(j)+3` 系の添字) | 2 回目以降の corlength≠0 の DELCOR が再現できない (乱数の列は同じ)。CAVI は DX・DY が変わらず HARM・PHI が変わる | 中 (CAVI)・低 (`cv(n)`) | 確定 (gs で再現。確かめた組み合わせは docs/31。CAVI は `CAVI RFC=(VOLT=1.5e6 HARM=64 DX=1e-5 DY=2e-5)` で HARM・PHI の値が gs と一致。`cv(n)` の原因は未初期化の読みで、ソースの読みと gs の実測) | v2 で修正 (波数 −1 の成分 cv(n) を対称に置く。CAVI は B40 のとおり DX・DY)。k64 モードは再現 (cv(n)=0。CAVI は HARM・PHI: `test_delcor_cavity_writes_harm_and_phi_like_gs`。2 回目以降の corlength≠0 は gs と数 % 違う) |
| B31 | 光学の wake (WAKEOPT) で、縦と横の表の片方しかない要素があると SIGSEGV | tffswake.f:102 | 落ちる | 高 | 確定 (gdb 相当: addr2line、修正して確認) | 直した SAD と同じ動作 |
| B32 | TrackParticles (TRPT、wake) で SOL 領域の始点の SOL にウェイクがあると SIGSEGV | tturn.f:567 (twxiwp) | 落ちる | 高 | 確定 | 直した SAD と同じ動作 |
| B33 | TRPT で、最初の `tapert` の時点で LOSSAMPL の外にいる粒子があると SIGSEGV | tapert.f:184 (`idtypec(l_track)`、l_track=0) | 落ちる | 高 | 確定 | v2 で修正 (最初の tapert で落とす。k64 モードは再現: `TrackRun.l_track0`) |
| B63 | SOL (BOUND≠0) の長さが周長に入らない (omega0・MomentumCompaction・空洞の位相・偏極時間に効く)。SOL に GEO=1 がないと S 座標が途切れて周長がさらに狂う | 周長の計算 (場所は未調査) | 周回周波数・運動量圧縮・偏極時間が SOL の長さの分だけずれる。GEO なしでは偏極時間が約半分、MomentumCompaction が約 100 倍 | 中 | 確定 (gs で確認: PLANAR + SOL で、SOL を除いた周長 45.6 のときと偏極時間が 15 桁一致)。意図か誤りかは未判定 | バグではないと判断 (SOL に L のキーワードは無い: SADHelp)。pysad は SOL の L を軌道・光学・S 座標・周長に入れない (gs と同じ。モードによらない: `test_sol_has_no_length`。以前の pysad は S 座標だけに入れていた)。GEO なしの S 座標の途切れは再現しない |
| B64 | `tmulke` (CALC・Emittance[] の MULT) は m=0 の枝 (L=0 の薄いキックと端) で `bsir0` を更新せず、m≠0 の枝 (スライス) では `bsir0=merge(bsir0+imag(cx)/aln,0.d0,m == 1)` と書く。`merge` は m==1 のときだけ加算し、m≥2 では 0 に戻すので、意図 (全スライスの積算) と合わないように見える | tmulte.f:6-43 (tmulke) | 縦オフセットのある MULT の POL で bsi の更新が抜ける可能性 | 低 | 疑い (ソースの読みのみ。gs の試験リングでは効かず、差が出る例が無い) | v2 で修正 (ψ = Im Σ k_n z^(n+1)/(n+1)! を入口で足し出口で引く。四極は k x y。POL のときは分割を細かくする)。k64 モードは K64 のとおり |
| B65 | RADLIGHT: `TrackParticles` の軌道 {t,x,y,z} に横方向の座標が入らない。tlinit が geo0 の軸 1・2 を設定せず (初期値 0)、位置は straight (dphi=0) の積み上げで BEND で向きが変わらず、BEND と曲げのある MULT の中は記録しない | tlinit.f (tlinit / tlstore / tlresult) | RadiationField の入力が進行方向の点だけになる。BEND を含む実リングでは意味のある場にならない (`example.oide.sad` は B を ANGLE=0 の STEER にして避けている) | 中 | 確定 (gs と差 0 で再現、[`test_radlight.py`](../pysad/tests/test_radlight.py))。意図か誤りかは未判定。なお RadiationField 自体は gs が SIGSEGV (B49) | v2 で修正 (軌道に横の座標が入り、BEND は円弧、BEND の中も `RADLIGHTDIV` 点記録)。k64 モードは再現 |
| B66 | 旧式の `WAKE` コマンドは `tfwake.f` を呼ばず `WAKECOMMAND[]` (init.n の AutoLoad) を評価する。`Packages/WAKECOMMAND.n` は配布に無く、fileopen のエラーを出して行の残りを捨てる | tffsa.f:976 (docs の旧記述は 972)、init.n:696 | WAKE コマンドが使えない (tfwake.f は死んだコード) | 低 | 確定 (gs で再現) | v2 は旧式の `WAKE L|T 要素 s1 w1 … ENDWAKE` を WakeFunction の表として登録 (CAVI だけ。行をまたいで書ける)。k64 モードは従来どおり |
| B67 | `MOMENTUM` の初期値が 30 GeV (`RsetGL('MOMENTUM',30.d9)`)。MOMENTUM を設定しないスクリプトの追跡・Emittance[] が 30 GeV で計算される | initb1.f:1698 | MOMENTUM の書き忘れに気づきにくい | 低 | 確定 (ソースと gs)。仕様の可能性が高く「癖」 (docs/11 A43) | 再現 (既定 30e9)。v2 は MOMENTUM を設定しないまま運動量を使う計算に入ると 1 度だけ警告 (値は変えない。k64 モードは警告なし) |
| B68 | SOL (BZ≠0) 領域内の四極 + RAD の追跡で、現在の gs (2026-10-19 ビルド) が、保存した参照 (`bench_level9`、10-14 の gs) と別の値を返す (S R T1 P1 の x が -.0027867898 から -.0014539780)。NORAD は一致 | tsol 系 (ソースのどこかは未特定) | 同じ入力で gs のビルドにより RAD の追跡結果が変わる。現在の gs が RAD ON・BZ≠0 の領域内の四極の力を入れていない疑いが強い (local/wake/verify_log.md の bench_level9 の項: `SB D1 QD SE` が `SB D1 SE` と 8 桁一致。10-19 と 11-14 -O0 の 2 ビルドは同値で、最適化の差ではない) | 中 | 疑いが強い (gs 側。ソースの箇所は未特定) | pysad は 10-14 の参照・物理と一致。現在の gs とは合わない |
| B69 | SOL 領域 (BOUND≠0 の間) 内の直線の BEND (ANGLE か K0 付き) の偏向を gs がほぼ無視する (BZ=0.5 T の領域で ANGLE=0.1 でも位置・角度が 1e-5 しか変わらない)。K64 の tsol.f は `tdrift(al,bzs,phiy,phix,…)` を呼ぶので、その直訳は ANGLE·距離だけ曲がる | tsol.f、tpara (領域内の BEND への K0/ANGLE の入り方は未調査) | SOL 領域内の BEND の追跡 | 低 | 疑い (原因は未確認。K64 側の癖かもしれない)。生の gs の出力は `local/wake/ps23.sad`〜`ps28.sad` | v2 で SOL 領域内の BEND・MULT・CAVI・APERT に対応 (直線の座標のままの一様な双極磁場。ローレンツ力の数値積分と一致)。k64 モードは未対応 (NotImplementedError) |
| B70 | `MessageName[f,"usage"]` (関数形) は f を無視して `Undefined MessageName General::"usage"` になる。`f::usage` の形は SAD の MessageName.n の表を引いて本文を返す | MessageName の処理 (未調査) | 関数形の MessageName が使えない | 低 | 実測 (gs で確認、原因は未調査) | v2 は関数形が f の本文を返す (無ければ `General::tag`)。`f::usage` の形は一致 (モードによらない) |
| B71 | 組み込みの `For[ini,cond,incr,body]` は 4 番目の引数 (本体) を実行せず、3 番目 (incr) だけを繰り返す (`For[i=0,i<3,i++,Print[i]]` は何も印字せず i=3)。`Packages/functions.n` の定義 `For[`ini_,`cond_,`incr_,`body_]:=(ini;While[cond,body;incr];)` は本体を実行するはずで、同じ式を別名 (For2) で定義すると gs も実行する。3 引数の `For` は未評価のまま | `Packages/functions.n:34` (読み込み時の名前解決と思われるが原因は未特定) | `For` の本体が動かない (`Do`・`While` は正常) | 低 | 実測 (gs で確認。`local/wake/` の f1〜f5 相当) | v2 は定義どおり本体を実行する (Break・Continue も効く)。3 引数は本体なしのループ。k64 モードは gs と同じ (本体を実行せず、3 引数は未評価。以前の pysad は本体を実行していた) |
| B72 | RADINT の `QUANTUM LIFETIME(hours)` が exp(2.3026·τ[秒]) で、時間への換算になっていない | INTGRL.f (RADINT の出力) | 量子寿命の表示値が桁違い | 中 | ソースの読み + 解析式 (τq = τE/2·exp(ξ)/ξ。gs での照合は記録に無い) | v2 で修正 (τ[秒]/3600。k64 モードは再現) |
| B73 | RADINT の `LINEAR CHROMATICITY` が −I9/4π (磁石の中だけの ∮γ ds で、色収差にならない。例: −2.36、実際は −7.30) | INTGRL.f (RADINT の出力) | GX・GY (クロマティシティ) の値 | 中 | ソースの読み + 運動量を ±1e-4 ずらした CALC のチューンの差分と 2e-5 で一致 (gs での照合は記録に無い) | v2 で修正 (運動量をずらした光学のチューンの差分 dν/dδ。6 極・偏向の寄与を含む。k64 モードは再現) |
| B74 | RADINT は MULT・SOL を黙って無視する。BEND の K0 も無視する (ANGLE だけ。仕様かもしれない) | INTGRL.f | MULT・SOL のあるリングの積分値に入らない (表示も出ない) | 低 | ソースの読み (K0 は仕様の可能性が残る) | v2 は場のある MULT・SOL があれば「入っていない」と表示する (積分には入れない。k64 モードは表示なし)。BEND の K0 は K64 のまま入れない |
| B75 | Emittance[]・EMIT (tradke) と追跡 (tradk) は、偏向磁石の端の角 e の蹴りを偏向角 θ に入れるが、場の中の行路 al1 を変えない。放射 ∝ θ²/al1 なので I4 の端の項 −η tan(e)/ρ² が 2 倍になる | tradke・tradk (放射の刻み)。RADINT (INTGRL.f) は正しい | 端の角のある偏向磁石のリングで Jx・εx・σδ が誤る (矩形磁石の例: Jx = 1.097 (正しくは 1.006)、εx 8% 小さく σδ 2% 大きい)。K64 の中で EMIT と RADINT が食い違う (docs/46 の Jx の 2.4% の差) | 高 | 確定 (gs で再現: Emittance[] の Jx = 1.0974、同じ gs の RADINT は 1.00575。厳密な幾何 + Robinson の式と RADINT が一致) | v2 で修正 (端の角のある面で行路 al1 を x·tan(e) だけ短くする。k64 モードは再現) |
| B76 | 追跡の tbendi (K1≠0 の BEND) は、放射で運動量が変わっても曲率半径 rhoe を入口の値のままにする | tbendi (追跡) | 減衰の和 Jx+Jy+Je が 4 にならない (K1=0.1 で 4.05、Jx=0.7385。正しくは Jx=0.6875) | 中 | ソースの読み + 解析式 (Jx = 0.68762) と RADINT・Emittance[] (gs での照合は記録に無い) | v2 で修正 (スライスごとの放射のあとに rhoe = rhob·(1+δ) を更新。k64 モードは再現) |
| B77 | Emittance[] の MomentumCompaction: RFSW OFF は −R56/C だけで入口の分散の項 Σ R5j η_j が抜ける (始点で値が変わる。FODO の例で 3.4% 大きい)。RFSW ON はシンクロトロンチューンから μ² で求めるので μ²/12 だけ大きい (0.25%) | Emittance[] (temit 系) | MomentumCompaction・BucketHeight・BunchLength・νs (√α に比例)。docs/46 の「BUCKET HIGHT の差」(0.010045 と 0.010035) の原因 | 中 | 確定 (独立な I1/C と 2e-4 で一致。K64 の値は pysad の k64 モードで gs と合う) | v2 で修正 (RFSW OFF は −(R56 + Σ R5j η_j)/C、RFSW ON は 4 sin²(μ/2)。k64 モードは再現) |
| B78 | undulator.f は陰的な刻みの右辺 `fgx=dpds*(ff*fx+gg*gx)`・`fgy=dpds*(ff*fy+gg*gy)` で、すでに dpds を掛けた fx・gx・fy・gy にもう一度 dpds を掛ける | undulator.f | 磁場があっても収束力が刻み幅に比例して消える (200 刻みと 2000 刻みで py が 10 倍違う)。K64 では UND に磁場が無い (B37) ので現状は効かない | 低 | ソースの読み + 刻み数を変えた数値実験 (gs では磁場キーワードが無く確認できない) | v2 で修正 (磁場つきのカーネル `und.und_seg1` は 1 回だけ掛ける。k64 モードの磁場なしの経路は影響なし) |
| B79 | phsrot.f の拡散の乱数は v = U·√λ (ベクトル) を成分ごとに正規乱数へ掛けてから N を掛ける。U が対角でないと共分散が N B Nᵀ にならない | phsrot.f | PHSROT の拡散 (RAD かつ FLUC) の共分散。K64 では PHSROT が動かない (B38) ので現状は効かない | 低 | ソースの読み (gs では PHSROT が動かず確認できない) | v2 で修正 (N U (√λ ⊙ ξ))。k64 モードは PHSROT が従来どおり動かない |
| B80 | ソレノイド + 双極の中の放射 (tradk) は 1 刻みの前後の向きの変化 (弦) から損失を出すので、らせんの回転角 bz·L が大きいと小さく出る (相対 (bz·L)²/12。bz·L = 1 で 8%) | tdrift.f (tradk) | SOL 領域内の BEND・MULT の放射損失。K64 では領域内の BEND の偏向が効かない (B69) ので現状は見えにくい | 低 | ソースの読み + 解析式 (gs での照合は記録に無い) | v2 で修正 (領域内の BEND・MULT は bz·L/ndiv ≤ 0.02 に刻む)。k64 モードは領域内の BEND・MULT が未対応 (NotImplementedError) |
| B81 | スピン追跡 sprot の縦磁場: `bz = bz0/pzm + bsi` で bsi = bz·al を p0 で規格化したまま z 軸の向きに足す。(1) δ≠0 で (1+a)B∥ の歳差が 1/(1+δ) 倍にならない (δ = −5e-3、0.12 rad で 6e-4 rad)、(2) 軌道が傾くと bsi·ẑ の運動量に垂直な成分が B⊥ に二重に入る、(3) 双極 + ソレノイド (tdrift) の `bsi0 − ak0x·y − ak0y·x` が T-BMT と 3.6e-3 rad 合わない | sprot (SOL 領域・UND の中のスピン追跡) | SOL 領域・UND の中のスピンの歳差 (POL + RAD の追跡) | 中 | ソースの読み + Thomas-BMT の数値積分 (v2 は 1e-7〜1e-5、K64 の式は 6e-4〜4e-3 ずれる)。gs での照合は記録に無い | v2 で修正 (SOL 領域と UND では運動量の向きの磁場の積分 ∫B∥dl/(Bρ) を運動量の向きに置く。k64 モードは再現)。Emittance[]・EMIT の POL の量 (同じ縦磁場の足し方) は B130 で v2.1.0 が修正 (k64 モードは K64 の式のまま) |
| B82 | 領域内の四極 (RAD) の tsolqur は、bsi を最初のスライスは半分、最後のスライスは半分で上書きするので 1 スライス分足りない。放射の基準の向き pxr0 が正準運動量のままで、最初のスライスに (bz/2)(y, −x) の見かけの偏向が入る | tsolqur (tquad → tsolqur) | 軸上の粒子のスピンの回転が半分 ((1+a) bz L の 1/2)、軸外の粒子の放射損失が 2.2 倍 | 中 | ソースの読み + 解析式・MULT の放射の式との比較 (gs での照合は記録に無い) | v2 で修正 (スライスごとに bz·aln/(1+δ)、放射の基準は力学的運動量。高速追跡も。k64 モードは再現) |
| B83 | スピン追跡 sprot の歳差の角が (1+aγ)·2 sin(θ/2) (外積 sin θ を中点の pz ≈ cos(θ/2) で割る)。刻みが粗いと θ²/24 小さい (0.26 rad の BEND を 3 刻みで 3e-4) | sprot (すべての要素の追跡のスピン) | 平面リングの 1 周のスピンの回転角が 2π aγ から 4e-3 rad ずれる | 低 | ソースの読み + 解析式。v2.1.0 の確認 (A2): 24 セル FODO (3 GeV、RF 4 MV) の 1 周のスピン角が k64 −1.1940907 rad、v2 −1.1935299 rad (差 5.6e-4 rad/周)。偏向磁石に `EPS=1E-5` を付けると gs −1.1935395510、v2 −1.1935395636 (1e-8 で一致)。k64 の誤差は刻みを細かくすると消え、v2 は 57 倍速く収束する。example/07 の 500 周後の保存値 0.8794 rad (gs・k64) は v2 の 1.1516 rad と違うが、**この K64 の刻みの誤差が 500 周に積算した値で物理の値ではない** (EPS=1E-5 で gs 0.3969、v2 0.3970) | v2 で修正 (厳密な θ。k64 モードは再現) |
| B84 | tcav.f・tmulti.f の `fw=fw0*(wsn+ws(n-1)*.5d0)` で、空洞が 2 つ以上に分割されると wake の蹴りの合計が (1.5n − 0.5)/n 倍になる | tcav.f、tmulti.f | 空洞・RF つき MULT の wake の大きさ (分割数に依る。20 MV (1 分割) と 200 MV (4 分割以上) で 1.25〜1.5 倍の差) | 中 | ソースの読み + 位置によらない wake の表でのエネルギー変化 (v2 は分割数によらず 5e-3 で一致)。k64 モードは gs と一致する追跡を再現 | v2 で修正 (台形の重み fw0·(ws(n) + ws(n−1))/2。領域の外の空洞・RF つき MULT も。k64 モードは再現) |
| B85 | `phis=merge(0.d0,merge(phic,w*trf0,autophi),trpt)` は TRPT でないとき phis = phic にするので dphis = 0 になり、足したばかりの DZ·w も PHI も効かない。光学 (tcave.f: phii = phic) とも合わない。tcave.f は autophi のとき dvcacc・ddvcacc を足さないので、空洞がすべて AUTOPHI だと実効の RF の傾きが 0 | tcav.f、tmulti.f、tcave.f | AUTOPHI のリングの追跡 (位相の差が消える)・Emittance[] の RF の傾き。SADHelp に AUTOPHI の説明は無い (意図は未確認) | 低 | ソースの読み。TRPT の場合は gs と全桁一致 (v2 の定義 = K64 の TRPT の動作) | v2 は TRPT のときの形 (phis = 0) にそろえ、追跡・CALC (TRPT の加速)・Emittance[] に AUTOPHI を実装 (従来は NotImplementedError。k64 モードは従来どおり) |
| B86 | tsol.f の trots の入口は z に DZ だけを足し、出口は DZ·(1 − dv) を足す | tsol.f (trots) | 入口の DZ≠0 (GEO が出口側の SOL にあるとき) の z の差 dv·DZ (3 GeV で 1e-13 m 以下) | 低 | ソースの読み (参考。gs では未確認) | v2 は入口・出口とも (1 − dv)。k64 モードは SOL のずれ・回転が従来どおり未対応 |
| B87 | `ttinit.f` の JITTER OFF (QUIET) の重心合わせの基準が、分布の作り方と合わない。px の基準は `epx·dp0 − ax/bx·ex·dp0 + dpx` だが、分布の px の平均は `epx·dp0 + dpx` (ax/bx の項は分散を含まない xi に掛かる)。結合 R1..R4 も基準に入っていない | ttinit.f (MEASURE・TRCOD・MAIN の TRACK TURNS=-1 の初期分布) | DP0≠0 で EX≠0 かつ AX≠0 (または結合あり) のとき、QUIET の重心が設計値からずれる (例: px の差 2.5e-5、統計誤差 3e-7) | 低 | ソースの読み + 統計 (n=200000 の JITTER ON の平均が v2 の基準と 4σ 以内、K64 の基準から 20σ 以上)。gs は MEASURE が無出力 (B42) で確認できない | v2 で修正 (分布の作り方と同じ写像で求めた設計の重心に合わせる。k64 モードでは MEASURE を無視。`measure.ttinit(center_fix=False)` で K64 の式) |
| B88 | TWAKE / LWAKE が ON のとき ttinit.f は pk = 0 (全粒子 z = 0) にする。wake 用の分布はコメントアウトされて残っていない (名残とみる) | ttinit.f | wake ON の MEASURE の分布が縦に広がらず、wake の分散への寄与が出ない | 低 | 疑い (ソースの読み。意図は未確認) | v2 は z の分布を普通に作る (`measure.ttinit(wake=True)` で K64 の z = 0)。k64 モードでは MEASURE を無視 |
| B89 | SPAC の力が 1/γ 倍に弱い。spfield は静止系の電位 φ' (= Σ 個数/r') の勾配をそのまま実験室系の場として使う (`ex = −∂φ'/∂x`、`ez = −∂φ'/∂ζ / γ`)。正しくは E⊥ = γ E'⊥、E_z = E'_ζ で、どちらも γ が 1 つ足りない | spkick.f (spfield) | SPAC の空間電荷の力が 5 MeV で 1/9.8、20 MeV で 1/39 | 高 | 確定 (gs で実測: 5 MeV・1e11 個の rms x 0.6226 mm は K/γ の式 0.645 mm に近く、K の式 1.639 mm と合わない) | v2 は正しい係数で SPAC を追跡に入れる (新機能)。k64 モードは K64 の SPAC (SPC1・SPC2・SPC4) をビット単位で再現していないので、SPAC は追跡に入らない |
| B90 | SPAC のキックの後、`px(i)=px(i)/pr` と、すでに p で割ってある px をもう一度 (1+δ) で割る | spkick.f:199-200 | δ≠0 の粒子の横の運動量が 5 cm ごとに 1/(1+δ) 倍になる | 中 | ソースの読み (spdeflect の入口の px は tturn の内部変数 px/p) | v2 は正準運動量に戻して足し、新しい δ で割る |
| B91 | MULT の `RADIUS` の空間電荷は、`spkick` の呼び出しの引数が合っていない (`sx,sy,sz` を余分に渡す)。さらに空間電荷のキックの長さ `al*wm` は最後の回だけ wm = wi/2 で、合計が L の (ndiv − 1/2)/ndiv 倍 | tmulti.f:222、541 | MULT の `RADIUS≠0` の SPAC が使えない | 中 | ソースの読み (呼び出しが壊れているので gs と比べる相手が無い) | v2 は MULT の RADIUS の SPAC を実装 (ANGLE=0・VOLT=0・SOL 領域の外。空間電荷のキックは ndiv 回とも L/ndiv)。入れられない MULT はメッセージを出して空間電荷なしで追う。k64 モードは空間電荷なし |
| B92 | SPAC のメッシュ: `nint(5.d0*dble(np/1000)**0.67d0)` が整数の割り算で、np<1000 では縦 10 枚に固定される。横は `max(4, nint(npz^(1/4)/2))` で、np=2000 でも 8×8 | spkick.f (spmesh) | SPAC の結果が粗いメッシュで決まる (gs の結果は x と y が 16% 違う) | 低 | ソースの読み + gs の結果 | v2 は 24〜64 ノードの 2 次の B スプライン |
| B93 | PSPAC の既定 (`PSPACNZ=1`) では rhocal が SIGSEGV で落ちる。縦のメッシュが 1 枚で、λ(z) の配列の外を読む | spch.c (rhocal) | PSPAC を既定のまま使うと落ちる (`PSPACNZ=64` を与えると走る) | 高 | 確定 (gs で再現) | v2 は `PSPACNZ<4` のときバンチに合わせた 64 ノードを自動で張る (PSPAC を追跡に入れる新機能。k64 モードは空間電荷なし) |
| B94 | PSPAC の λ(z) の配分 (rhocal): セルの前半では `rhoz[ii]+=rzi; rhoz[ii]+=rzin;` と両方を同じノードに足す (後半は ii と ii+1)。重みがセルの中央で不連続 | spch.c (rhocal) | 線密度の配分の精度 | 低 | ソースの読み | v2 は 2 次の B スプライン |
| B95 | PSPAC の周の終わりのキックは `pos(lend)-sspac0` で、`sspac0` は最後の要素の中点 (WSPAC 用の変数)。最後の SPCH から終点までの分が入らない。始めは `sspac1=0` で、区間を途中から始めても s=0 からの長さになる | tturn.f:201、270-277、584-597 | PSPAC のキックの長さの合計が周長と合わない | 低 | ソースの読み | v2 は区間の始めから数え、終わりに最後の SPCH からの残りを入れる |
| B96 | PSPAC の状態 (周の数、メッシュ) が C の static 変数で、`TrackParticles` を呼び直しても続きから数える。`PSPACNTURN` 周で解放した後も使う | tspch_.c | 2 回目以降の TrackParticles の PSPAC が 1 回目の状態に依る | 低 | ソースの読み | v2 は呼び出しごとに作り直す |
| B97 | `tspac` (粒子の rms からのガウス分布のキック) の係数が `1/(p0 (1+p²))` = 1/(βγ·γ²) で、β が 1 つ足りない。K64 のどこからも呼ばれない死んだコード | tspac.f | なし (呼ばれない) | 低 | ソースの読み | 使わない (WSPAC の twspac は係数 1/(β²γ³) で正しい) |
| B98 | SPAC の中の粒子数が少ないと落ちる・NaN になる。spmesh は縦を nza = 10 枚 (np < 4000) のスライスに npz = ceil(np/nza) 個ずつ分け、`dz2 = zz(itab(np)) − zz(itab((nza−1)·npz+1))` を計算する。np < 9·ceil(np/10)+1 では添字が np を超えて配列の外を読み、等しいと最後のスライスが 1 粒子で dz2 = 0 | spkick.f (spmesh) | np = 2・15・33 は SIGSEGV、5 は Bus error、9・10・11・19・25・28・81・82 は NaN (全粒子が失われる)。RADIUS が小さい絞りで粒子が減る途中でも当たる (docs/46 の「小さい RADIUS>0 で gs が落ちる」の原因) | 高 | 確定 (gs で再現し、原因まで特定) | v2 は粒子数によらず走る (2・5・15・33・81 個で試験。SPAC を追跡に入れる新機能。k64 モードは空間電荷なし) |
| B99 | tvcorr.f の corlength≠0 の枝は `s = 1/sqrt(4Σw²+2)` で、複素数 c の絶対値の 2 乗平均を 1 に規格化する。tcorr.f は DX = Re(c)·delx、DY = Im(c)·dely と使うので DX の rms は delx/√2。同じサブルーチンの corlength=0 の枝は成分ごとの規格化 (1/2 なし) で整合しない | tvcorr.f | corlength≠0 の DELCOR のずれの rms が delx の 1/√2 になる | 低 | 規格化の根拠は確定 (v2.1.0 の C: K64 は複素数の ⟨\|c\|²⟩ = 1 にする。corlength = 0 の枝は実部・虚部それぞれ分散 al/n で、2 つの枝で実部・虚部の分散の約束が揃っていない。試験: k64 は ⟨\|c\|²⟩ = 1.0、v2 は 2.0)。「DX の rms = delx」が作者の意図かは推定のまま (作者の文書による裏づけは無い) | v2 で修正 (s = 1/sqrt(2Σw²+1)。DX の rms = delx)。k64 モードは K64 のまま (delx に 1/√2 を掛けても同じ) |
| B100 | `trcodf`・`trcods` が添字 15..18 を直書きしている (DETR を入れる前の番号。今は DETR DX DPX DY で DPY が入らない)。`trcodm` は名前で書いてあり正しい。FIXCOD が保存する基準 `twiss(*,ndim,*)` を読む手段が K64 の FFS に無い (読むのは旧補正パッケージ mbmp.f と DRAW のモニターの読みだけで、K64 の FFS から呼ばれない) | trcoda.f (trcodf、trcods) | TRCOD FIXCOD の基準が使えない | 低 | ソースの読み (gs は TRCOD が SIGSEGV: B42) | v2 は EX EPX EY EPY DX DPX DY DPY の 8 関数を保存し、`Twiss["RDX",k]` (基準)・`Twiss["DDX",k]` (現在 − 基準) で読める (REFERENCE と同じ参照光学に置く)。k64 モードは未対応 |
| B101 | `if (master(is) .lt. 0) then is=-master(is)` は、固定形式では変数 `thenis` への代入になり、何もしない | nalign.f:41 | TRCOD ALIGN で分割された磁石の master の扱いが効かない | 低 | ソースの読み (gs は TRCOD が SIGSEGV: B42) | v2 の TRCOD ALIGN は分割された磁石を、先頭 (is) と末尾 (ie) へ広げる。k64 モードは未対応 |
| B102 | id = 20 の枝が今の型番号では SOL に当たる (DX DY を「設計値 − ずれ」にする。旧い型番号の名残と判断) | nalign.f | TRCOD ALIGN が SOL に作用する | 低 | ソースの読み | v2 は SOL・MONI を動かさない |
| B103 | `GETMEASURE` が目標値を書き換える fit 条件の関数番号 32・34 が、今の番号では PEPY・PZPX (旧い関数表の番号) | tfgetm.f | GETMEASURE の目標値の対象が確定できない | 低 | ソースの読み | v2 は DX・DY の条件に書く (式が x・y の重心の差の半分なので軌道の関数と判断。旧い DDX.. だった可能性は残る)。k64 モードは認識しない |
| B104 | `tffsa.f` の STN_COD が値を入れていない局所変数 meas0 を範囲に渡す | tffsa.f | TRCOD STNCOD の範囲が不定 | 低 | ソースの読み | v2 は FIXCOD で保存した範囲 (FIXCOD がまだなら注意を出して何もしない) |
| B105 | DRAW のビームサイズ (kv = 28..32) の追跡は `trsize` が常に偽 (TRACKSIZ は tffs.f:298 でコメントアウト、tffsa.f:151 で `trsize=.false.`)。追跡の結果を使う行もコメントアウト (tdrwdt.f:121〜130)。kv = 33・34 のモニターの読みは `2*kv-51` = 15・17 (DETR を入れる前の番号) で、呼び出し元 (mbmp.f の mcmon) が K64 の FFS から呼ばれない | tdrwdt.f | なし (到達しないコード) | 低 | ソースの読み | 実装しない |
| B106 | tradke の dbeam の b = bzh·0.5。力学的運動量は px + bzh·y なので b = bzh のはず (B11 + b(2B13 + b B33) の b が半分)。bzh の定義の別の約束がある可能性が残る | photons.f | SOL 内の放射のビーム行列の「ビームの広がり」の項だけ (ε は全桁変わらない) | 低 | 式の一貫性 (数値の効果は 1e-13 以下: 結合のある SOL つきリングの Emittance[] は 13 桁変わらない。拡散の項は ⟨u²⟩⟨P⊥²⟩ で、横方向の減衰 ⟨u⟩ に対し ⟨u²⟩/⟨u⟩ ≈ 1e-6 だけ小さい) | v2 で修正 (b = bzh。B39 と同じ `fixb39` の v2 だけ)。k64 モードは b = bzh·0.5 のまま |
| B107 | RADINT の BMAG の `CKL1Y = −cos(θ)·DKL1` は sin の誤記。θ=0 (回転なし) の K1≠0 の BEND でも勾配ベクトルの y 成分が 0 にならず、Y 方向の分散 η_y に 2 k h η_y の項が余計に入る。(QUADRU の `−DKL1` と軌道のずれの扱いは B126) | INTGRL.f (BMAG) | 結合のあるリングで K1≠0 の BEND に Y 方向の分散があるときの I4V: 24 セル FODO + skew 四極 + K1 = ±0.1 の BEND で I4U + I4V が物理的な I4 (結合のある 5x5 行列) と +1.6% (K1 = +0.1)・+2.7% (−0.1) 違い、DU・DV・Je が 1 次の摂動論と 0.3 違う | 低 | 確定 (y 成分を +sin(θ)·K1 にすると I4U + I4V が独立な I4 と 1e-7、DU・DV・Je が 1 次の摂動論と 3e-3 以内。K1=0 または結合なしでは同じ。θ ≠ 0 の数値の確認は無い: 傾いた BEND は閉軌道が大きく動く) | v2 で修正 (`CKL1Y = +sin(θ)·DKL1`)。k64 モードは gs と全桁一致 (再現) |
| B108 | ANGLE=0 の BEND (tsteee・tsteer、K0 だけの補正磁石) と K0 つき MULT の端の角は、B75 と同じ「端の蹴りを偏向角に入れるが行路は変えない」作りのはず | tsteee・tsteer 等 | 端の角つきの補正磁石 (AE1・AE2 つきの挿入光源の極) の放射 | 中 | 確定 (硬い端の一様磁場の極の厳密な幾何 (円弧と直線) とロビンソンの定理 Je = 2 + dln(Σ B_i² l_i(δ))/dδ。24 セルの FODO に極 3 個の組 (K0・−2K0・K0) を入れ、AE = (0.3, 0.3) で ΔJx = 厳密 +0.01440、v2 +0.01487 (1.032)、K64 +0.02992 (2.08 倍)。(−0.2, 0.4) で 1.037 と 2.16 倍。K64 は端の項が約 2 倍 = B75 と同じ誤り) | v2 で修正 (`_tsteee` の `tradke` に tin = tan(AE1)・tout = tan(AE2) を渡す。ANGLE≠0 の `tbende` と同じ)。k64 モードは K64 のまま。残る 3〜4% は EMITDIV = 20 の有限の分割と 2 次の項 (傾き 30 mrad の極では 2〜3%)。K1≠0・K0≠0 の ANGLE=0 BEND と MULT の端の角は pysad が AE を読まず未確認 **v2.2.0 の更新**: MULT・K1≠0・K0≠0 の ANGLE=0 BEND の端の角 AE1・AE2 を実装し (CALC・Emittance[]・追跡。k64 は gs と全桁一致)、放射の行路も v2 で修正した (B143)。追跡 (`tsteer`) の端の放射も v2 で修正 (ΔJx は厳密な幾何の 1.03〜1.04 倍)。傾き 30 mrad の極の 0.02〜0.03 は K0 に依らない 2.7% の系統差 (AE に比例、極の長さ・間隔に依る。原因未特定) |
| B109 | tbgrote の trans(5,6) が、写像の微分 zf/pzf·(1 − pr/pzf·r33·dpzidp) (コメントアウトされた式) ではなく r33·(dpzf·pr/pzf + dpzi/pzi) − dr33 になっている。dr = 0 でも dpzi·(pr+1)/pzi が残る | tchge.f (tbgrote) | 2 次の小さい量 | 低 | 確定 (tchge の解析的な 6x6 行列と写像 `_tchge_cod` の有限差分 (3 つの条件、入口・出口): K64 は (5,6) が 1.3e-5 (= px² + py²) ずれ、v2 は全成分 4e-12 以内) | v2 で修正 (コメントアウトされた式。dr = 0 では 0)。k64 モードは再現 |
| B110 | tthine (長さのある6 極など) の bsir0 が `∓aki·imag(cx)/al` で、次元が合わない | tthine | なし (放射なしで呼ばれ、入口と出口で打ち消す) | 低 | ソースの読み | 実装していない (影響なし) |
| B111 | 放射の反作用の角度の項の符号が tradke (Emittance[]・閉軌道: `ddpx = +0.5·dpx·dg`) と tradkn (追跡: `−0.5·dpx·dg`、dg < 0) で逆。エネルギーを失う間は同じ磁場が運動量の小さい粒子をより強く曲げる (BEND の中の平均 δ = dg/2 → 偏向角 +\|dg\|/2 倍) ので追跡の符号が物理的に正しい | photons.f (tradke:1417、tradkn:1247) | RADCOD ON の閉軌道が追跡で閉じない (3 GeV・周長 105.6 m の FODO で 1 周 x −4.3e-7 m、z 2.4e-6 m。RADTAPER ありで 2e-6・1e-5 m)。閉軌道・U0・チューン・ダンピング・エミッタンスが相対 1e-5〜1e-3 動く (大きなリング A (3 GeV、603 m): εx +1.9e-4、ダンピング率 +2e-4) | 中 | 確定 (gs で再現。符号は解析的な導出と、反転後に閉軌道が 1 周で 1e-9 m で閉じることによる。作者の意図の裏づけは無い) | v2 で修正 (tradke の ddpx・ddpy と行列 D・B の微分の符号を反転。RADTAPER の追跡の BEND の重みも tturne と同じ (3δ_l + δ_l+1)/4 にした: これは「閉じる」ことを根拠にした選択で、どちらの重みが物理的かは未判断)。k64 モードは再現。v2 の Emittance[] (RADCOD) は gs と相対 1e-4〜1e-3 ずれる **v2.2.0 の更新**: 符号は確定 (放射反作用つきのローレンツ力の数値積分と v2 の tradke が 6 桁一致。k64 は 1 台の偏向磁石・既定の EMITDIV = 1 で Δpx が 9% 違う)。RADTAPER の BEND の係数は (δ_l + δ_l+1)/2 が物理的なので、v2 は Emittance[] も追跡もこれにした (3/4 は K64 の粗い離散化に合わせた係数)。v2 の Emittance[] が gs とずれるのは K64 の誤差。微分の誤りは B138 |
| B112 | 放射の古典半径が粒子の質量によらない (`re0 = elradi`、電子の値。`rclassic = charge²·re0`)。追跡の `tradkn` は光子数 (5/(2√3)) α γ θ と臨界エネルギーから δ を下げる超相対論的な式で、β の因子が無い (U0_K64 = U0/β³) | tfsetparam.f、photons.f (tradkn) | MASS≠電子の粒子の放射: 30 GeV/c の陽子 (γ = 32) の U0 が古典式 6.88e-4 eV に対し K64 は 1.2651 eV (mp/me = 1836 倍)。β の因子は γ = 32 で +0.15% (β < 0.99 では 1/β² ほど大きい) | 中 | 確定 (K64 の値は gs と全桁一致。古典式 U0 = (2/3) r_p mc² β³γ⁴ Σθ²/L と v2 の追跡の 1 周の δ の変化が 1e-6 で一致) | v2 で修正 (`re0 = elradi·elmass/amass`)。β の因子は直していない (γ ≫ 1 の粒子だけ意味を持つ)。k64 モードは再現。TouschekLifetime・Intrabeam・RADINT の `rclassic` は K64 のまま (陽子での確認なし) **v2.2.0 の更新**: β の因子は直していない (要判断)。k64 は gs と全桁 (1.5・3・30 GeV/c の陽子の U0)。v2 の U0 は古典式のちょうど 1/β³ 倍 (1.5・3・30 GeV/c で 64%・15%・0.15%。γ ≫ 1 では 3/(2γ²)) |
| B113 | tfgeo の `inext`・`iprev` (「同じ型 (BEND) の次の要素。間は長さ 0 の要素だけ」を繋ぐ。リングの終わりは先頭に回って探す) により、先頭が長さ 0 の MARK の輸送路の最後の BEND は自分自身に繋がる (iprev = inext = 自分)。tbend・tbendi は `iprev == 0 .and. fb1 /= 0` で入口側、`inext == 0 .and. fb2 /= 0` で出口側に F1 の放射の断片を入れるので、この BEND は F1 の断片の放射が無い。隣り合う BEND の間も同様 | tfgeo.f、tbend.f、tbendi.f | F1・FRINGE のある BEND の RAD の追跡 (1 個の BEND の輸送路で 3%、隣り合う BEND の 12 セルのリングで 3.6%)。意図かどうかは未確認 | 低 | 確定 (gs で再現: 1 個の BEND (INS) の 1 周 δ = −2.17142054679e-5)。意図は疑い | pysad は K64 の動作を再現 (`track6n.bend_links`、両モード、gate なし。v2.0.0 は断片を常に入れていた pysad の不具合で、K64 の動作に合わせた) |
| B114 | ずれた BEND の CALC (tbende の tchge) の出口で、tbrote の `dz1 = dz·c0 − dx·s0` (入口は `dz·c0 + dx·s0`)。出口の (dx1, dz1) の両方の s0 の符号が逆で、出口のドリフトが逆向きになる | tchge.f (tbrote) | 矩形磁石 1 台の DX = 2 mm でチューンが (2.7e-5, 2.0e-5) 動く (追跡は動かない)。K1=0 の BEND の CALC と追跡の写像が行列で 5.2e-4 ずれる。大きなリング A に DX を入れた CALC のチューンが追跡と 3e-7 (v2 は 6e-8) | 中 | 確定 (gs で再現。厳密な幾何 (弦に垂直なずれの帯の中の円弧) で矩形は光学が不変、扇形は x-px の行列が 8.6e-6。v2 の CALC は 2e-8 で一致) | v2 で修正 (出口も入口と同じ回し方 R(phig)(dx, dz))。k64 モードは再現。K1≠0 の BEND は B117 |
| B115 | RADINT は ANGLE=0・K0≠0 の BEND (挿入光源の極: 基準軌道は曲がらず磁場だけ) を積分に入れない | INTGRL.f (RADINT) | ASTRID2 + MPW (0.58 GeV、周長 45.7 m) で U0 6.17 keV (Emittance[] 7.80 keV)、εx 1.268e-8 (Emittance[] 1.056e-8) | 中 | 確定 (gs で再現。I2 は Σ(ANGLE+K0)²/L と 1e-6) | v2 で修正 (K0 を曲げ角として扱う。端の角なし、TOTAL BENDING ANGLE には足さない)。U0・I2・I3・σδ は使える (U0 1.1e-3、σδ 2.2e-3)。Jx は 0.7%、εx は 3.5% Emittance[] とずれる (極の中の分散の伝搬を扇形の行列で近似)。k64 モードは再現 |
| B116 | `Polarization`・`PolarizationTime` は、逆向きの磁場 (挿入光源の交互の極) があっても符号つきの率 sdamp = n·gintd = Σ r_i (b_i·n) を緩和率に使う。逆向きの極が相殺する | temit の偏極の量 | ASTRID2 + MPW で P = 92.376% (一定。独立な式は 70.16%)、τ = 28842 s (独立な式 21886 s。実際より 31% 長い)。K64 の RADPOL 追跡の緩和は Σ r (b·n)² の側 | 中 | 確定 (独立な式 P = (8/5√3)⟨ρ⁻³⟩/⟨\|ρ\|⁻³⟩、1/τ = (5√3/8) r_e ƛ_C c γ⁵⟨\|ρ\|⁻³⟩ と、K64 の追跡の緩和率の比 (MPW あり/なし 1.311) で確認) | v2 で修正 (緩和率に符号に依らない Σ r_i (b_i·n)² (`sdamp_rel`) を使う。P = 0.70190、τ = 21915 s (3.9e-4、1.3e-3 の差))。k64 モードは再現 (P は NaN: B55)。n が傾く・非平面のリングで Σ r (b·n)² が最適かは未確認 **v2.2.0 の更新**: n が傾く・非平面のリングは厳密式 Σ r (1 − 2/9 (n·v)²) を v2 に入れた (`spin.ST_EXACT`。垂直の偏向・傾いた n・水平磁場の MPW で独立な Derbenev-Kondratenko の数値積分と τ が 1e-5〜2e-3)。分子の運動論的偏極 (`spin.KINETIC`)、脱偏極率の放射の刻みごとの積分 (`spin.STEP_RATE`) も入れ、P は 0.1〜0.6% (脱偏極率が分母の 15% 以上では 2〜3%) で合う。A2 の Σ r (b·n)² は n が b に垂直で 0 になる (水平磁場の MPW で τ が 32% 短い) |
| B117 | `tbendi` (K1≠0 の BEND の追跡) の端は薄い写像 (`pxi = px + tanp1·x/ρe` など) で、端面までのドリフトと端面から扇形の面までの弧 (x tan e の項) を入れない。CALC の `tbedge` (tbedgedrift → tbedgemaxwell → tbedgebody) は厳密 | tbendi.f | K1≠0 の BEND の CALC と追跡の x-px の行列が軌道のずれ (x0 = 1 mm) で 1.1e-5 違う (E = 0.25 で 2.7e-6)。K1=0 の追跡 (tbendcore) は厳密で CALC と合う | 低 | 確定 (ローレンツ力の数値積分の厳密な写像と CALC が 1e-11、追跡は 1e-10〜9e-9) | v2 で修正 (`DISFRIN=0` かつ端の角 ≠ 0 のとき、追跡の端を `tbedge` と同じ厳密な写像にする。追跡は厳密な写像と 1e-10)。k64 モードは再現。四極の磁場が端面でなく半径方向の面で始まる模型は K64 の仕様で、CALC・追跡とも同じ (直さない) |
| B118 | `HurwitzZeta[s,a]` の s が複素で、a を 1 以下に畳んだ値が 0.5 以上 1 未満 (a = 0.5・1.5・3.5 など。0.25・1・2・10 は正しい) のとき誤る。`Zeta[s]` の臨界線上の零点は 1e-16 の雑音 | gamma.f | 複素の s の HurwitzZeta: 相対 1e-6〜1e9 の誤り (例 `HurwitzZeta[2+I,.5]` = 3.9549+2.1791 I、真値 3.5075+2.0314 I)。Zeta 40 点・HurwitzZeta の他の点は 1e-13 以内で正しい | 低 | 確定 (別に書いた Euler–Maclaurin と 1e-10、Zeta は Borwein のアルゴリズム・関数等式と 1e-11) | v2 は真値を返す (`specfn_cx.py`)。**k64 モードも真値を返す** (gs の誤った値は再現しない: 再現には複素の不完全ガンマなどの直訳が要る。要判断)。gs の値が必要なら gs を使う **v2.2.0 の更新**: 複素の特殊関数の gs の誤りの全体は B158 (732 式の走査。追加した 15 関数で 641 式が gs と 1e-9) |
| B119 | 非整数次の `PolyLog[s,z]` が複素の z で、Re z < 0・Im z ≤ 0・\|z\| > 1 のとき誤る。整数次と第 3 象限以外は正しい | gamma.f (`cpolylogl`) | 例 `PolyLog[1.5,-3.28-3.54 I]` = 3.47+9.43 I (真値 −2.09−1.06 I) | 低 | 確定 (冪級数・spence・積分表示・整数次の反転公式と 2e-11 で比較) | v2 は真値 (冪級数・μ = ln z の展開・倍角の公式の再帰)。k64 モードも真値 (B118 と同じ扱い。要判断)。負の実数の z は上側 (arg z = +π)、実軸上の z > 1 は従来どおり切断線の下側 |
| B120 | `HypergeometricU[a,b,z]` で、b が整数だと複素の z でも NaN (実数の z は B46 と同じ)。\|z\| = 20 で相対 1e-8 の誤差 | gamma.f | 複素の z の U (b が整数) が使えない | 低 | 確定 (積分表示・接続公式・多項式・漸近展開と 1e-10。gs が正しい 18 点は gs と 2e-13) | v2 は真値 (Re a > 0 は積分表示、a が負の整数は多項式、ほかは a の漸化式と Kummer の変換)。k64 モードも真値 (要判断)。Re a ≤ 0 の複素の a・b は SADError、z = 0 の U は未対応 |
| B121 | `Union[list, SameTest->f]` は、並べ替えも比較もせず、結合した先頭の要素だけを返す | Packages/functions.n | `Union[{1,2,3,4,5,6},SameTest->(Abs[#-#2]<1.5&)]` が {1,3,5} でない | 低 | 確定 (gs で再現) | v2 は真値 ({1,3,5})。k64 モードは gs と同じ |
| B122 | `WriteBeamLine[..., Name->"L2"]` (文字列) が `Insufficient number of Names` (`Name->{"L2"}` とリストで書かないと通らない) | Packages/beamline.n | 文字列の Name が使えない | 低 | 確定 (gs で再現) | v2 は文字列も受ける。k64 モードは gs と同じ |
| B123 | `tphrec` の `c1 = (cost·ppx1 + sint·ppy1)/pp` は、すでに規格化した `ppx1 = ppx/pp`・`ppy1 = ppy/pp` をさらに pp (≈ 曲げ角、数 mrad) で割る。c1・s1 は偏光の回転角 χ の cos・sin のはずなので \|c1\| ≤ 1 | photons.f (tphrec) | 光子の偏光の Stokes ξ1 が ±1 を超え 1e3〜1e6 (gs の PhotonList の 8 列目が −3735.37、−1.0015e6 など) | 低 | 確定 (gs で再現。radangle1 の ξ2² + ξ3² = 1 と合わない) | v2 は pp で割らない (\|ξ1\|・\|ξ3\| ≤ 1、ξ1² + ξ2² + ξ3² ≈ 1)。k64 モードは gs と同じ (全列一致) |
| B124 | tdgeo の CHI3 の微分の式は、領域の終端で GCHI3 が 0 に戻される (tsgeo の `tfderotgeo`) ので幾何の変化と合わない | tdgeo (tsgeo) | 9e-4 を返すが幾何は動かない。CHI3 を使う幾何の fit は K64 でも収束しない | 低 | ソースの読み (領域の終端のリセットと式の不整合) | 領域内の四極の DX・DY の fit は gs と同じ 2 回の Newton で収束 (変数が 1e-12 以内)。CHI3 の fit は未対応のまま |
| B125 | `spdepol` の脱偏極率が emittance/減衰率に比例する。2 次の零位相の項 (`spintrm` の `cm0 = −mat/(e^{cφ} − 1)`、位相が打ち消し合う 2 次の積で cφ = −2λ のとき分母 ≈ 2λ) が、平衡のエミッタンス (すでに 1/(2λ) を含む) にさらに 1/(1 − e^{−2λ}) を掛ける。正しい率は光子の跳びの拡散 D だけで決まり、減衰率に依らない。`spintrm` のガウスの項の 2π も余分 (幅 0 の極限で 2π·cm0。正しくは √π/a) | photons.f (`spdepol`・`spintrm`) | `Polarization` が使えない (B55 の NaN は桁あふれで、直しても 4.2 GeV・12 MV の 12 セル FODO で P = 1.27e-4 (率 4.3e-5 /周)。独立な率 2.3e-8〜4.5e-8 の約 1000 倍。追跡は (10.4 ± 1.6)e-8) | 高 | 位置の特定は数値の依存性 (e に比例、λ に反比例、ssprd に依らない) による。K64 の式の導出 (Mathematica 由来の 6 次までの展開) には戻っていない。独立な値は導出・線形の軌道 + 厳密なスピン回転の MC・追跡 | v2 は独立な式で置き換え (SLIM の 1 次の摂動にスピンの歳差の位相の変調 (シンクロトロン側帯) を厳密に取った項を入れる。`spin.sideband_depolarization`)。4.2 GeV・12 MV で追跡の 0.43 倍 (過小評価。軌道の非線形が残る)。k64 モードは K64 の式のまま (B55) |
| B126 | RADINT の QUADRU は、四極の中の軌道の位置 (dx, dy) に Twiss の閉軌道 DX・DY だけを渡し (四極自身の据え付けのずれ DX・DY を引かない)、K1<0 の四極では勾配ベクトル CKL1 の符号が逆 (DKL1 の符号つき)。(B107 の 2 つ目) | INTGRL.f (QUADRU) | 軌道がずれた四極の RADINT の Jx・Jy: 24 セル FODO で ΔJx = +0.0028 (正しくは +0.0078。追跡 +0.00781、Emittance[] +0.00786)。Y 方向のずれ + skew 四極の ΔJy は符号が逆 | 中 | 確定 (平面リングの 3x3 行列の放射積分と RAD つき 6 次元追跡の 1 周写像の両方。経路長が変わらない QF : QD のずれの比 (DX = 1e-4・2.2974e-4) で比較) | v2 で修正 ((dx, dy) = 閉軌道 − 要素の DX・DY、勾配ベクトルは \|K1\|·(dx, dy)/r)。k64 モードは gs と全桁一致 (再現)。SEXT (フィードダウン)・BEND (DX・K1 の勾配) の軌道のずれは RADINT に入れていない (K64 が軌道を無視する作り)。経路長が変わるずれでは RADINT (δ = 0 の約束) と Emittance[] の J が δ の分だけ違う |
| B127 | WSPAC の `twspace` は `p1 = (1+δ)·(gammab(l) + gammab(l+1))/2` で、`tturne.f` の最後の `tintrb(…, iend+1)` (l = nlat、一周の最後のスライス) では `gammab(nlat+1)` が配列の外で 0 になり、p1 が半分、レンズ ∝ 1/(p1² h1) が 8 倍強い | tintrb.f (twspace)、tturne.f | WSPAC の Emittance[] のチューンシフトが 3〜4% 強い (0.2 GeV の 12 セル FODO で x 4.3%・y 3.0%。大きなリング A でラスレットの積分に対し K64 は 1.0085・1.0109) | 中 | 確定 (gs の `tintrb.f` の複製に出力を足して再ビルド: l = nlat で `GAMB 75 391.3902 0.000000 195.6951 391.3902`。108 個のレンズは一致、最後の 1 個だけ 8.0 倍) | v2 で修正 (p1 = (1+δ) p0)。v2 のシフトは一次のラスレットの積分と 0.2〜1.8% (大きなリング A・PF・ASTRID2)。k64 モードは再現 (gs と全桁一致) |
| B128 | `SynchroBetaEmittance` が大きなリング A (3 GeV、周長 603 m、νs = 0.0041、RF 1.4 MV) で 6.18e-9 を返し、Emittance[] の 1.55e-10 と 40 倍違う。AzimuthalModes 3〜8 で変わらず、12・20・30・32 以上で 3.19e-9・1.46e-9・8.50e-10・7.87e-10 と動き続ける。νs の入力・エネルギー・6 極の強さにも依らない | SynchroBeta (temits) | 大きなリングで使えない (独立な 6x6 の平衡 1.5173e-10 と、RAD + FLUC の 6 次元追跡 (60 粒子、4000 周) の 1.6e-10 と合わない)。PF・ASTRID2・SP8 (大きなリング B、8 GeV) は Emittance[] と 0.5〜1% で合う | 中 | 確定 (v2.2.0 で原因を特定: 右端の欄の更新。v2.1.0 までは「疑い」で、gs も pysad k64 も 6.18e-9 だった) | v2 は値を変えず、AzimuthalModes m と 1.5 倍の m の emitx が 20% 以上違えば警告を出す。k64 モードは gs と同じ (警告なし)。使うときは AzimuthalModes を変えて収束を確認し Emittance[] と比べる。`AzimuthalModes->1` は pysad が IndexError (gs は 8.6e-7 を返すが意味がない) **v2.2.0 の更新 (原因を特定して v2 で修正)**: K64 の `tefsetup` が 1 周の写像 X' = T X + h1 の源泉 h1 に閉軌道 c(δ) (1 次の分散 η δ を含む) をそのまま使う (固定点が c にならず、2 次元の回転では平均が c の半分になる)、`tesumb` が x_b の 2 次のモーメントに η η σδ² だけを足し相関 η σδ² (`beam(16:19)`) を入れない。リングの始点の分散 η ≠ 0 で εx が始点に依り数倍〜40 倍過大 (24 セル FODO の始点を QF・BEND・BEND の途中にして 6.54・6.54・4.74 倍、gs でも同じ。η が小さい PF・ASTRID2 は 1%)。確度は確定。v2 は源泉から 1 次の分散を引き、軌道の δ 依存の効果を解析的な 2 項に分けて入れる (`tesumb_v2`): Emittance[] の 1.0003〜1.0005 (FODO)・1.0084 (大きなリング A、AzimuthalModes 3。32 で 1.0011)・1.0045 (SP8 (大きなリング B))・1.005 (PF)・1.006 (ASTRID2)。k64 は K64 のまま (gs と 2e-7)。警告 (上記) は v2 の大きなリング A では出なくなった |
| B129 | `CorrectOrbit` の `Solver->{"Micado",n}` (`Condition` なし) が、選んだ nx 本の補正量でなく、全補正量の最小二乗解の先頭 nx 個を使う。ハウスホルダー変換のあとの `call tsvd(a,b,x,eps,svd)` で、新しい tsvd は配列の大きさを `size()` で取る (コメントアウトされた旧 `tsolva(a,b,x,nx,nx,mdim,eps)` は nx を渡していた) ので全 n 列を解き、x(nx+1:n) = 0 にする | pmicad.f | 24 セル FODO (BPM 48 台、垂直ステアラ 24 台) の補正後の垂直 rms が補正前 9.41e-5 m に対し Micado n=1 で 1.46e-4 (悪化)、n=3 8.2e-5、n=6 1.19e-4、n=8 1.75e-4 (本数を増やしても単調に下がらない) | 中 | 確定 (gs で再現。選んだ nx 列の最小二乗 `lstsq` と v2 が 1e-9、選び方は独立な貪欲法と同じ) | v2 で修正 (上三角の nx x nx だけを解く。n=1 4.74e-5 … 24 本 2.06e-5 と単調に下がる)。k64 モードは gs と 1e-6 以内で再現。`Condition` つき (CorrectLocal) の Micado は選んだ列だけで `psvdcond` を解く元の分岐で K64 でも正しい。水平・2 平面同時・`Bump`・`SVDmode->1` つきは未確認 **v2.2.0**: 水平・2 平面同時の Micado も gs の k64 と 5e-6 で一致し、v2 は単調に下がる。`SVDmode->1` は B157 |
| B130 | Emittance[]・EMIT の `tradke` は、運動量の向きの磁場の積分 (縦磁場 bz·al、端の ψ) を z 軸の向きの磁場として足す。運動量が z 軸から傾くと (傾いたソレノイド・軸を外れた粒子) (1) 向きの変化から作る B⊥ (すでに運動量に垂直) に垂直成分が二重に入り、(2) δ≠0 で (1+a)B∥ の歳差が 1/(1+δ) 倍にならない (追跡側は B81) | photons.f (tradke)、tsolqu 系 | ずれて傾いたソレノイド (GEO=1 の SOL に DX DY DPX DPY) の POL: SpinTune の追跡との差 1.5e-4〜1.1e-3、n0 の向き 2.5e-3〜1.5e-2 rad (小さな FODO・大きなリング A・PF・ASTRID2・ASTRID2 + MPW で同じ) | 中 | 確定 (回転した縦磁場の歳差の閉じた式 (傾き 1〜50 mrad・1〜3 T・0.3〜1 GeV の 5 通り) と追跡が 5e-7、Emittance[] の srot が 1e-8。Thomas-BMT の数値積分と 3e-5) | v2 で修正 (B⊥ = 向きの変化 (asin)、B∥ = [bz·al + (bx px + by py)·al/pz + ψ]/(1+δ) を運動量の向きに分けて足す。SpinTune の差は 1e-13〜5e-7、n0 の差は 0〜1e-5 rad)。k64 モードは K64 の式 (`spin.SOL_BPAR = False` でも旧式)。1〜3 T・0.3〜3 GeV で確認、bz ≥ 3 1/m の QUAD・MULT は 1e-4 台 |
| B131 | SOL 領域内の CAVI (VOLT=0 でも) を `tsol` が `krad = .false.` で呼ぶので、追跡 (POL) でスピンが歳差しない (リング全体で追跡と Emittance[] の R が 0.12 違った) | tsol.f | 領域内の CAVI のスピンの歳差が追跡に入らない | 低 | ソースの読み + Thomas-BMT (v2 は 1e-9 台で一致)。意図は未確認 | v2 は CAVI も放射・スピンの刻みに通す (物理的に歳差させる)。k64 モードは K64 のまま |
| B132 | PhotonList のエネルギー列 E は、γ でなく gammab = p/(mc) = βγ を掛ける | photons.f (tphotonlist) | 光子のエネルギーが相対 1/(2γ²) 小さい (3 GeV で 1.45e-8) | 低 | 確定 (gs と pysad の E 列の比が、ほぼ全光子で一定の 1.4e-8。βγ を掛けると全列一致) | v2 は γ (h0)、k64 は p0 (gs と全列一致) |
| B133 | WSPAC ON の CALC・CALC6D は入口以外の Twiss が NaN になる (0.2 GeV のリングで PBUNCH = 1e4 でも `Horizontal/Vertical Unstable` と NaN) | qtwiss.f (qwspac)・qtwissfrac | WSPAC の光学の CALC が使えない (Emittance[]・EMIT は使える: B127) | 中 | 確定 (gs で再現。原因は未調査) | v2 は CALC6D と RFSW ON の CALC (4D) に空間電荷を入れる (6D の経路。Emittance[] のチューンと 2e-6、一次のラスレットと 0.15%)。K64 の `qwspac` の直訳ではないので gs との比較はできない。k64 モードは従来どおり空間電荷を入れない (gs が NaN) |
| B134 | CALC の解析行列と `tbende` は、DISFRIN=1 でも Maxwell のフリンジを入れる (追跡の `tbendi` は DISFRIN に従う) | tbfrie.f・tbende.f | DISFRIN=1 の BEND の CALC と追跡が合わない | 低 | ソースの読み (gs での確認は記録に無い) | 直していない (要判断。k64・v2 とも CALC は DISFRIN に依らず、追跡は DISFRIN=1 で薄い写像のまま) |
| B135 | 追跡 (track6: `tbshift`・`tent`) は、BEND (K1=0・K1≠0 の `tbshift`) と ANGLE=0 の BEND (挿入光源の極の `tsteer`) の座標の回転 (ROTATE) をスピンの系に掛けない (四極の `tsolrot` は回す)。スピンが回さない系の鉛直軸のまわりに歳差する | tbshift・tent (追跡のスピン) | ROTATE≠0 の BEND の POL の追跡。ROTATE = 90°・K0 = 0.064・L = 0.0875 m の BEND 1 台で、追跡のスピン回転行列は y のまわりの回転 (R[0,2] = −0.1479)、厳密な Thomas-BMT の積分と Emittance[] の srot は x のまわり (R[1,2] = +0.1479)。水平磁場の挿入光源 (ROTATE = 90°) のリングの追跡のスピン | 中 | 確定 (pysad の k64 = K64 の翻訳が厳密な Thomas-BMT の積分と 0.1 以上違い、Emittance[] は一致。gs の追跡のスピンで同じかは未確認で、K64 の翻訳からの推定) | v2 で修正 (`track6._spin_rot_cs` が座標と同じ回転をスピンに掛ける。厳密な積分と 1e-9)。k64 モードは K64 のまま |
| B136 | 偏向磁石の本体で、鉛直方向の傾き py ≠ 0 の粒子のスピン回転が厳密な値と (2+a) θ py だけ違う。運動量の向きの磁場の歳差 (1+a) B∥ ds (B∥ = B·v̂) を入れない分 (1+a) θ py と、さらに θ py (フレームの回転と運動量の向きの外積の系の混在と推定。未確定)。K64・追跡・Emittance[] の 3 つとも同じ | tradke・tradk・tbend の本体 | 扇形 (E1=E2=0)・K1=0・L = 1.2 m・15°・0.5 GeV・py = 5e-3 で R[0,1] が厳密 −1.462e-3、追跡・Emittance[] −4.05e-3・−4.07e-3 (差 −2.58e-3 ≈ −(2+a) θ py)。D[z,py] の差は −(2+a)θ。n0 の傾き・SpinTune・脱偏極の D に効く。軸上・x・δ は 1e-12〜1e-7 で厳密に一致。Emittance[] の D と追跡の D が鉛直の軌道・分散が大きいリング (0.5 GeV・y 軌道 6 mm) で食い違う (追跡の δ キックへの応答 0.0963 に対し、Emittance[] の N の δ 列 0.0025) 原因の一部と推定 | 中 | 確定 (厳密な Thomas-BMT の積分 [`example/v2_next/P1_bend_vertical_angle_spin.py`](../example/v2_next/P1_bend_vertical_angle_spin.py)。原因の分解のうち θ py の項は推定) | **直していない** (要判断。直すには tradke・tradk・tbend の全ての要素の本体に B∥ を足し、フレームの系の扱いを揃える必要があり、B58・B64・B130 のスピンの結果が全部動く)。k64・v2 とも K64 のまま |
| B137 | 追跡の K1=0 の BEND のスピン軌道結合 D が、放射・スピンの刻みに強く依存する (A2 の「K1=0 の BEND の追跡は分割が粗い」の延長) | tbend + tradk | 同じ BEND の D の行 0 の y・py 列が、刻みなし −0.075・0.054、0.05 m の刻み −0.032・0.099 (0.005 m まで収束)、Emittance[] −0.094・0.068。収束値は厳密な積分に近い (B135 を除く) | 低 | 確定 (刻みの収束と厳密な積分) | 直していない (K64 の分割の粗さ。B136 と合わせて要判断)。P の脱偏極率の部分の 2〜3% の残差の原因の一つ |
| B138 | `tradke` の微分 (`irad > 6` の節) に 4 点の誤り・近似がある: (1) 区間の始点の座標に関する微分 `tr1` を、終点の座標に関する微分を要する D・B の更新 `D += tr1·(T + D)`・`B = (I+tr1) B (I+tr1)ᵀ` にそのまま使う (正しくは tr1 ← tr1·tr2⁻¹)、(2) `dddpx` の 2 項目が ddpx·ddg (正しくは dpx·ddg)、(3) `tr1(1,:)` の x' = x + ddpx·c1 の c1 = al1/(3 pr) の微分が無い、(4) `dxpa = (xp·dxp)/xpa/pr²` (正しくは /xpa/pr⁴) で δ≠0 の軌道に 2δ の相対誤差 | photons.f (tradke) | 放射の減衰行列 D・拡散行列 B (Emittance[]・EMIT の減衰率・エミッタンス・チューンシフト) が区間が粗いとき 1e-3 程度ずれる。RADCOD の 1 周の解析 T + D と写像の数値微分 (4 次の中心差分) の差: BEND 1 台 5.1e-7、QF 1 台 1.1e-8、1 セル 1.9e-6、24 セルのリング 4.3e-6 (3 GeV の FODO)。WSPAC の Emittance[] のチューンシフトで 2e-3 (相対) | 低 | 確定 (写像の 4 次の中心差分と v2 が 1e-10 以内) | v2 で修正 (`perfkern.tradke_db` の `exactd`)。BEND 1 台 1.7e-10、QF 1 台 6.7e-16、1 セル 6.4e-9、24 セル 1.9e-7。k64 モードは K64 のまま。残る 6.4e-9・1.9e-7 は (z, δ) 成分だけで B139 による |
| B139 | BEND の T[z, δ] (`tbendebody0` の (5,6)) の近似: px² に比例する項 | tbendebody0 | 1 周で 1e-7 (大きなリング A で 2e-5 とされた RADCOD の T の誤差の主因。RADCOD OFF の T だけでも 1.2e-7)。px = −5.9e-5 で 6.5e-9、px = 1e-5 で 1.5e-10 | 低 | 確定 (写像の数値微分との差) | 直していない (K64 の仕様の近似。v2 で直す案は要判断)。k64・v2 とも同じ |
| B140 | 追跡の `tsteer` (ANGLE=0 の BEND) の出口の端の角 AE2 は、`if(n .eq. n2 .and. fb2 .ne. 0.d0)` の中でだけ tanp2 を設定し、**fb2 ≠ 0 のときだけ**入る。CALC・Emittance[] の `tsteee` は常に入れる | tsteer.f | AE2 だけのステアリングの追跡の一周写像のチューン (gs) は AE なしと同じ 0.1340952 (FB2 = 1e-6 を付けると 0.1338835。CALC は AE2 を入れた値)。AE2 = 0.3・K0 = 0.0015 のステアリングで追跡と CALC のチューンが 1.2e-4、大きなリング A (AE2 = 0.3 の補正磁石) で 1.4e-4 違う | 低 | 確定 (gs で再現: `P3_steer_ae2_tune.sad`) | v2 で修正 (AE2 を常に入れる。追跡の一周写像の固有値と CALC のチューンが 6e-6 (K0 の 2 次))。k64 モードは K64 のまま |
| B141 | `tmulte` (CALC・Emittance[] の MULT) は出口の端の蹴り (`tmulke` の apsi2) を最後の半スライス (`tsolque`) の前に置く。追跡の `tmulti` は後 (`tmulk`) | tmulte.f・tmulti.f | K0 = ±0.01・AE = 0.3 の MULT で CALC と追跡のチューンが 3e-5〜1.2e-4 (AE1 + AE2) 違う | 低 | 確定 (gs の CALC (pysad k64 が全桁一致) と追跡 (全桁一致) の差) | v2 で修正 (CALC・Emittance[] の出口の蹴りを最後の半スライスの後に。差 3e-8)。k64 モードは K64 のまま |
| B142 | K1≠0・K0≠0 の ANGLE=0 BEND の追跡 (`tbend`) が `tmulti` に端の角として psi = 0 を渡す (`0.d0,0.d0,0.d0,0.d0` = bz・phia・psi1・psi2)。CALC・Emittance[] の `tbende` は `tmulte` に AE (apsi) を渡す | tbend.f・tbende.f | この BEND の追跡が AE を無視し、CALC のチューンと 5e-4 以上違う。**gs はこの BEND で USE の中で SIGSEGV する (B02) ので実行で確認できない** | 低 | ソースの読み (gs で再現できていない。要確認: 読みが正しいか) | v2 で修正 (追跡も AE を使う。CALC のチューンと 2e-6)。k64 モードは K64 のまま (AE を無視) |
| B143 | `tmulte`・`tmulti` (MULT) の端の角 AE1・AE2 の蹴りは、放射の場の中の行路を変えない (B75・B108 と同じ型) | tmulte.f・tmulti.f | 端の角つきの MULT (勾配つきの極を含む) の放射: 端の項が約 2 倍 (厳密な幾何との ΔJx の比 K64 : v2 = 2.0 : 1.03) | 中 | 確定 (硬い端の磁石の厳密な幾何 + ロビンソンの定理。表は `history/v2next_P3.md` §3.3) | v2 で修正 (追跡 `tsteer`・`tmulti` の `zr0` に x tan(e)、Emittance[] の `_tmulte` の `tradke` に tin・tout)。厳密な幾何と 3〜4%。残る 2.7% (AE に比例、極の長さ・間隔に依る) は原因未特定。k64 モードは K64 のまま。B108 の MULT への拡張 |
| B144 | 補正磁石に ROTATE と DX・DY を同時に付けたとき、gs の CALC のチューンが gs の追跡 (pysad k64 と全桁一致) の一周写像のチューンと合わない (B09 と同じ型: ロール + ずれの tchge の行列) | tchge 周辺 (CALC) | 12 セルの FODO で ROTATE = 0.1・DX = 1e-4 のステアリング: gs の CALC の νy 0.1258111、追跡 0.13371 (νx は 0.2668346 対 0.2668366 で合う)。pysad の k64 の CALC は gs と 1e-4〜7e-4 違う (K64 を再現していない) | 低 | 確定 (gs の CALC と追跡の不整合) | v2 で修正 (v2 の CALC は追跡と 5e-7 で一致。B09 の修正が入る)。k64 の CALC が K64 を再現しない点は未調査のまま |
| B145 | `tsteer` の F1 の端の模型は、fb1 ≠ fb2 のとき x が厳密な Maxwell 場の端の運動と (K0/al)(fb1² − fb2²)/24 だけ違う (fb1 = fb2 では 1e-10 で一致) | tsteer.f | K0 = 0.003・al = 0.2・fb = (0.3, 0.0001) で 5.6e-5 m (`P3_maxwell_ramp.out.txt`) | 低 | 確定 (Maxwell 場のベクトルポテンシャルの正準方程式の数値積分) | 直していない (K64 の模型の近似。入口の面の取り方の判断が要る)。k64・v2 とも同じ |
| B146 | SOL 領域内の端の縦磁場の積分 ψ (QUAD の k x y、MULT の k0 y + k1 x y + …) が DISFRIN に関係なく入る | tsol 系 (追跡)・tradke (Emittance[]) | フリンジなしの模型 (DISFRIN=1) の領域内の MULT・QUAD のスピンが Thomas-BMT の積分と 1e-4 ずれる (MULT は bz に依らず 2e-5〜1.6e-4) | 低 | 確定 (Thomas-BMT の積分。K64 の領域内の MULT は未対応なので、v2 の式 (E) が受け継いだ形。意図は未確認) | v2 で修正 (ψ は DISFRIN=0 のときだけ。DISFRIN=1 で 9e-7 以下)。DISFRIN=0 の ψ が物理的に正しいかは未確認 |
| B147 | PSPAC のキックは平均のエネルギーを保存しない: 縦の 3 次スプライン λ' と形状関数の不整合 | spch.c・psn.c | 一様な円・長さ 0.1 m のバンチ (4096 個) の平均 δ が +1.8e-5 (rms δ 3.6e-4 の 5%、v2 は 1.9e-7)。粗い縦メッシュ (4 セル) では平均 δ が rms の 12% | 低 | 確定 (gs と pysad k64 が全桁一致する値。物理的な基準はキックの合計の保存) | k64 モードは K64 の形状関数を再現 (`spch_k64.py`、gs と rms 12 桁以上)。v2 は別の形状関数 (2 次の B スプライン) のまま |
| B148 | PSPAC の C コード (spch.c・tspch_.c) の不整合: (1) 電位・λ のメモリを PSPACNTURN 周で解放するのにその後の呼び出しが解放済みのメモリを読む (値は凍結した電位と同じだが未定義の動作)、(2) メッシュの端のセルの `makebcu` が配列の外を読む、(3) 縦のスプライン係数 `rhoz_sp[nz]` が未初期化、(4) 1 セッションの 2 回目の `TrackParticles` で PSPAC が働かない (SPC8、既出)、(5) `PSPACNZ = 1` で SIGSEGV (B93、既出)、(6) NPARA > 1 の `TrackParticles` は粒子を子プロセスに分け、PSPAC は各プロセスが自分の粒子だけで電位を作る (粒子数が 1/NPARA で、キックの係数 Kfac も子ごとの np で変わる。NPARA = 1 に替えると rms x 1.23379 → 1.23401 mm) | spch.c・tspch_.c | PSPAC の結果が NPARA・セッションの履歴に依る。A1 の gs の 16% の x・y の非対称は (6) に由来する可能性がある (再現は未確認) | 低 | ソースの読み + gs ((4)(5)(6) は gs で確認) | k64 は TrackParticles ごとに新しく始め、PSPACNZ < 3 をエラーにする (gs の SIGSEGV を再現しない)。端のセルの配列の外は 0 とする。NPARA = 1 でだけ gs と比べた |
| B149 | gs の `TrackParticles` は WSPAC ON でも `Emittance[]` のあとでは入力をそのまま返す (`CODPLOT; WSPAC; EMIT;` が要る。`beamsize` が EMIT の plot でだけ作られるため)。追跡の WSPAC のビーム行列には MINCOUP の下限が入らず、垂直のエミッタンスが 0 の平らなリングは NaN (ビームの高さ 0 のまま追跡する) | tturn の `twspac`・tfsetparam | 追跡の WSPAC を使う台本が `CODPLOT; WSPAC; EMIT;` を要する。結合のある 0.5 GeV のリング (2e10 個) の試験粒子のチューンシフトは gs (−2.938e-5、−1.210e-3) で、MINCOUP の下限つきの空間電荷つきの行列 (−2.776e-5、−3.807e-4) と合わず、ゼロ電流の行列 (−2.982e-5、−1.201e-3) と 1〜2% で合う | 低 | 確定 (gs で再現: `P4_wspac_tracking_beam.gs.txt`) | 既定は K64 に合わせてゼロ電流の行列 (v2)。v2 のオプション `WSPACSC = 1;` が空間電荷つきの行列 (MINCOUP の下限つき) を使う (K64 にはない)。k64 は追跡に WSPAC を入れない |
| B150 | RAD ON・bz≠0 の SOL 領域内の QUAD は光子を記録しない (B68 の延長の疑い。領域内の MULT・DRIFT は記録する) | tsol 系 | 3 GeV・BZ = 1 T の領域の QUAD の光子が 0 個 (MULT は 955 個、gx は要素の範囲、エネルギーの合計 312.5 eV/粒子 = 古典値 285 eV の 1.10 倍) | 低 | 確定 (gs で再現: `P4_sol_photons.gs.txt`。原因は未特定) | v2 で、領域内の DRIFT・QUAD・MULT (ANGLE=0) の光子を実装 (位置が軌道上、エネルギーの合計が古典値と 8% 以内)。BEND は未対応。k64 は K64 のまま |
| B151 | `MEASURE`・`TRCOD` の追跡 (`tracker` の `tturn`) が、先頭の要素より後ろを通ると全粒子を落とす: 出力が出るのは位置が先頭の要素 (入口の MARK) のときだけ。それ以外の MEASURE は無出力 (`FFS["MEASURE"]` は np = 1)、TRCOD は SIGSEGV (B42 の更新: 先頭の MARK の位置は動く) | trcoda.f・trackb.f (`tturn` の入力) | MEASURE・TRCOD が先頭の位置でしか使えない | 中 | 確定 (gs の 10-19 ビルドで、小さな FODO・大きなリング A・PF・ASTRID2、エネルギー・エミッタンス・フラグを変えて再現。落ちる原因は未特定) | k64 は先頭の位置の MEASURE・TRCOD を gs と同じ統計 (4 桁) で実装し、それ以外は無出力・警告 (`SIGSEGV (B42)`)。v2 は追跡を通せる (値は gs と比べられない) |
| B152 | `NISTACK$FIXSEED@Push/Pop` が MEASURE で効かず、FIXSEED でも NOFIXSEED でも種が戻らない (MEASURE を続けて呼ぶと分布が変わる。位置が先頭でない MEASURE (無出力) も種を進める)。MEASURE の初期分布の重心 (B41・B87) とは別 | NISTACK.n (`GetSeed[] := SeedRandom[]` の戻り値が種になっていない可能性) | MEASURE・TRCOD の再現性 | 低 | 確定 (gs で再現。原因は推定) | k64 は種を戻さない (gs と同じ)。v2 は戻す (K64 の設計の意図) |
| B153 | 入れ替わり (DETR ≥ 1) をまたぐ純粋な 90 度回転 (COORD CHI3 = ±π/2) の要素の位相の跳び ±π が、`atan2(y, x)` の x < 0 で y が丸め誤差 (0 か ±1e-31) になるため丸め誤差の符号で決まり、一周の NX・NY の整数部が不定 (`qtwiss.f` の `dpsix=atan2(u32,aa)` など) | qtwiss.f | pysad の k64 は 3 つのリング (ASTRID2・PF の INS・大きなリング A) で整数部が gs と +1・+1・−1 ずれる (小数部は 8 桁で一致)。gs は 3 つとも 2 つの跳びが打ち消し、本体の物理的なチューンを返す | 低 | 確定 (gs と pysad の 3 リング。gs が 3 つともたまたま打ち消した可能性: 丸め誤差の符号は gs の計算順で決まる) | v2 で修正 (入れ替わりに入る要素の跳びを +π、出る要素を −π と決めて打ち消す。`lattice._exch_phase`)。k64 は K64 のまま (丸め誤差の符号) |
| B154 | `TRCOD` の粒子数が TRCOD 時の `NP` でなく、最後の `tfsetparam` (CALC・EMIT・MEASURE など) の時の値 (なければ 1000)。`NP` の既定は 1000 | trcoda.f・tfsetparam | `NP=200; CALC; TRCOD ORIGIN` は 200、`CALC; NP=200; TRCOD` は 1000 | 低 | 確定 (gs で再現。他の `tfsetparam` を呼ぶ命令の全リストは未調査) | k64 は gs と同じ (`ffscmd._k64_np0`)。v2 は NP の指定を要求する |
| B155 | `CorrectOrbit[…, Bump->b]` は `CorrectOrbitBump[kind,mon,str,opt,op]` と 5 引数で呼ぶが、定義は `CorrectOrbitBump[orbitkind_,bumpset_,bpm_,str_,opt_,options___]` (bumpset が抜けている) | Packages/correction.n | `Index 15 out of range` で補正は何も起きない。`CorrectOrbitBump` の直接呼びは gs が 300 s で終わらない (比較していない) | 低 | 確定 (gs と pysad が同じエラー) | 直していない (k64・v2 とも K64 と同じエラー。`MakeBump` は gs と 1e-12) |
| B156 | MAIN の要素定義の式の末尾の演算子 (`1+`・`2*3+`) が構文エラーにならず値 1 になる (`2*3+` も 6 でなく 1) | MAIN の式の評価 | 入力ミスが黙って通る | 低 | 確定 (gs で再現) | k64 は gs と同じ (値 1)。v2 は構文エラー。ほかに gs は `SIN[0.1]`・`45DEG` を構文エラーにする (pysad の `SIN[0.1]` は値 0 を入れていた誤りで、両モードで構文エラーにした。`45DEG` は k64 だけ構文エラー) |
| B157 | `CorrectOrbit` の `SVDmode->1` が Micado の n を無視する (n = 1・2・3・6・12 の結果が互いに同じ。n = 0 とも 1e-8) | Packages (`LinearSolveConditionedS` が `nx` を使わない?) | `SVDmode->1` で Micado が働かない | 低 | 確定 (gs と pysad k64 が一致。`LinearSolveConditionedS` のソースは未確認) | 直していない (k64・v2 とも同じ)。B129 の `SVDmode` なしの経路は v2 で修正 |
| B158 | 複素引数の特殊関数の gs の誤り・不正確 (732 式の走査で、pysad が独立な基準と合い gs が違う 91 式): `LogGamma` (Re z < 0 で虚部の分枝が 4π 違う)、`Erf`・`Erfc` (Re z < 0: `Erfc[−1.2+.7I]` gs −.0664+.1202I、真値 2.0664−.1202I)、`ExpIntegralEi` (複素引数 12 点すべて: `Ei[.3+.4I]` gs −.629+2.409I、級数 .1595+1.389I。D.n の式で、実引数の B44 と同じ系統)、`ArcCosh` (6 点で符号が逆)、`FresnelS`・`FresnelC` (\|z\| ≥ 12 で 1e38 や符号違い)、`BesselY` (複素引数は NaN)、`BesselK`・`AiryAi/Bi`・`ExpIntegralE`・`Factorial` (精度 1e-9〜1e-1) | gamma.f・D.n | 複素引数の該当の関数 | 低 | 確定 (scipy・級数・恒等式と比較。B118〜B120 は既出の別の誤り) | v2 は真値。**k64 も真値** (gs の誤りを再現しない: 規則が単純な LogGamma の分枝・ArcCosh の符号・Erf/Erfc・ExpIntegralEi は再現できる可能性がある。要判断) |
| B159 | MARK の `DETR` キーワードが INS の入口で期待どおりに効かない: gs は `DETR=1` を書いても入口の DETR を R から計算し直した値 (R1 の値が出る) にし、伝搬も DETR を書かないときと途中まで同じ (QS の後で違う)。`DX` を書くと B40 の取り違えで 1e200 になる | MARK の値の添字の取り違え (B40 の周辺か) | 生の DETR ≥ 1 を入れ替わった入口にする INS (A1 の「入口の規格化が違う」の原因) | 低 | 確定 (gs の 6 通り。原因は未特定) | 直していない (pysad は生の DETR ≥ 1 を入れ替わった入口として使う。k64 を gs に合わせる修正は途中で止めた: gs の挙動を説明できない) |
| B160 | **誤りではなく K1 の暗黙の定義 (U1 で再分類)**: K1≠0 の BEND の本体は、中央面の磁場 B_y = B0 + (K1/h) ln(r/ρ0)、B_r = K1 y ρ0/r (= (1 + hx)∂b/∂x が一様、B''(ρ0) = −h B'(ρ0)、x y² の項なし) の中の運動を厳密に実装している。elegant (CSBEND K2 = 0) など B' が x に一様な定義とは K1 が同じ記号でも別の磁石で、本体に K2 = −h K1 の6 極成分を足したのと同じ。SADHelp の BEND の K1 は B'L/Bρ とだけ書き B'' を指定しない | BEND の本体の写像 (`tbend.f` tbendiinit・tbendicorr・tbendibody、`tbendi.f`、CALC の `tbende.f` tbendecorr) | オフモメンタムのチューン (ξ・∂²ν/∂δ²)、LMA・DA・クロマティック ADTS を B' 一様の定義と比べるとき。δ = 0 の光学・エミッタンスには効かない。n = 0.4 の扇形磁石 16 台 (3 GeV) の dν/dδ は gs・pysad・独立な積分 'log' が (−0.17213、+0.21082)、B' 一様 (独立な積分・解析式・elegant) は (−0.60246、+0.73786)。ASTRID2 の自然な ξy は SAD −11.8、elegant −8.2 | 低 (仕様・癖。SAD → 他コードの変換は K2 = −h K1 を足す。docs/11 A59) | 実装の同定は確定 (ソース、導出、独立なローレンツ力の積分と 16 台の非線形の写像が 1e-10、MULT (ANGLE≠0) の別経路が 1e-12 で同じ)。「定義であって誤りではない」は中〜高 (BEND の K1 の B'' は文書にない)。T1 の「疑い〜確定の間」を更新 | 直さない (k64・v2 とも既定は K64 と同じ 'log')。**v2 のオプション (v2.3.0 U2): `compat.set_bend_k1_def('lin')` (環境変数 `PYSAD_BEND_K1_DEF=lin`) で B' 一様の定義 ('log' + 6 極 K2 = +hK1) にできる (k64 では無視)。n = 0.4 の 16 台の dν/dδ が (−0.60246、+0.73786) になり、独立な積分・解析式・elegant (K2 = 0) と一致。ASTRID2 の自然な ξ も elegant (K2 = 0、EDGE_ORDER = 2) と 4 桁で一致。`history/v23_U2.md`** CALC の tbendecorr の dk の不整合は B161 |
| B161 | CALC・Emittance[] の BEND (K1≠0) の解析行列 `tbendecorr` の dk = −K1 xr (1 − xr (1 − 2 xr)) (xr = x/ρ0) が、同じ関数の kick の導関数 −K1 xr (1 − xr/2 + xr²/3) と xr の 2 次の項が違う (xr = 0.001 で −0.05%、0.1 で −3.5%、0.3 でたまたま一致)。写像 (kick) と行列の不整合 | `tbende.f` tbendecorr の 86 行 (kick)・87 行 (dk) | CALC・Emittance[] の転送行列だけ (追跡は kick を使うので影響を受けない)。スライス 1 cm・K1 = −0.06 /m²・xr = 1e-3 で t[1,0] に 3e-10、リングのチューンで 1e-8 程度以下と見積もった (全体の量は見ていない) | 低 | 確定 (U2: CALC の行列が kick の数値微分と合わないことを `bendecorr` で確認、K1 の強い扇形磁石のリングで CALC のチューンが独立な積分と合わない: δ = 1% で 3e-5、2% で 1e-4、4% で 5e-4。K64 の作者の意図は不明。実際のリングへの影響は小さい: ASTRID2・大きなリング A で ξ に 4e-6、∂²ν/∂δ² に 0.01、ν(δ = 2%) に 5e-7、Emittance[] の EMITX・SIGE に 1e-9 (相対) 以内) | **v2 で直した (v2.3.0 U2、`compat.bend_dk_fix()`、`perfkern_bend.py` `bendecorr` の dk を kick の導関数 −K1 xr (1 − xr/2 + xr²/3) に)。k64 は K64 のまま。** `history/v23_U1.md` §8、`history/v23_U2.md` |
| B162 | **K64 の Touschek の損失率 (`TouschekTable`・`TouschekLifetime[INF, INF, nz]`) の式の適用範囲の限界 (誤りというより近似)**: `tintrb.f` の `touckf` (emit.n の `TouschekF`) は、平坦ビームの Bruck・Piwinski の関数 D(ε) の当てはめ (C(x)/x の有理式 × e^{−x} × log。厳密な関数との差は x = 1e-6〜10 で ±1%) に、運動量の広がりを ptrans² = (位置を固定した運動量の共分散 P の 3 固有値の和) として入れたもの。x-y 面の運動量の広がりが等方でない (結合が小さくない、丸いビーム) と、2 次元ガウス分布の厳密な積分 (Bessel I0) と合わない。ASTRID2 (0.58 GeV、IBS 込み PBUNCH 1e9、一様な受容 0.4〜1.6%)・KEK-LS (3 GeV、0.4〜2%) で、独立な式 (T2) に対する K64 の寿命の差は、結合 1% で +6%〜+0.1%、10% で +1〜+14%、100% (丸いビーム) で −14%〜+21% (受容が小さいほど平坦な式に近く、大きいほど過大)。 | `tintrb.f` touckf・`ptrans` (tintrb の `ptrans = sqrt(Σ eig)`)、emit.n `TouschekF` | 結合が 1% 以上のリングの Touschek 寿命・局所の受容の寿命の係数 (最大 ±20%)。TouschekTable の横の振幅の受容 (nx・ny が有限) も同じ式に依存 | 中 | 確定 (独立な式・3 次元の一般の積分・実験室系の対のモンテカルロ (厳密な相対論) が互いに 0.3% 以内で一致し、K64 だけ違う。elegant (標本化を細かくしたとき) も独立な式と −1.6%〜+0.3%。T2、`history/v23_T2.md`) | k64 は K64 の式のまま。**v2 は `TouschekLifetime[INF, INF, nz]`・`TouschekLifetimeLMA[]` の既定を Piwinski 型の式 (`touschek.touschek_piwinski`) にし、`Formula->"K64"` で表に戻せる**。`TouschekTable` (EMIT の報告・横の振幅の受容) は K64 のまま |

## 1. 落ちる・動かない (重要度 高)

### B01 MAIN の TRACK が SIGSEGV で落ちる

- 現象: MAIN 部で `EMIT USE=RING;` に続けて `TRACK USE=RING NPART=5 TURNS=10;` を実行すると、フラグ一覧 (RFSW RADCOD RAD ... SMEAR) を表示した後に `Program received signal SIGSEGV`。`NPART=1 TURNS=-1` でも同じ。`NX=(1,2)` のリストは `*** $TRACK *** list is not allowed` で拒否される (`NX=1` と書けば通るが、その後同じ場所で落ちる)。MAIN の `EMIT` は `Use EMIT_TANCE or Emittance[] within FFS.` と言うだけで何も計算しない。
- 再現: [`docs/history/work_main_track.md`](history/work_main_track.md) §3 の最小スクリプト。
- 期待される動作: 古典 SAD と同様に `EMIT` が正規モード行列 r・閉軌道・EMITX/Y/Z・SIGE を作り、`TRACK` が追跡する。少なくとも未設定のとき SIGSEGV ではなくエラーメッセージを出す。`NX=(1,2)` の形式も古典 SAD では許される。
- 原因: バックトレースを `nm` (obj/AMD64-Linux5/sad1.exe) で引くと `toplvl_ → doact_ → acttra_ → track_ → tracka_ → __trackdlib_MOD_tinip` (+0x1db)。`tinip` の `tmap(xa,-1)` が `use temw, only: rn=>r` の行列を使うが、MAIN の EMIT が r を作らないので未設定。FFS 内で `EMIT` すれば r はできるが、FFS から MAIN に戻って `TRACK` を実行する経路がない (`END;` の後は SADScript の対話、`GetMAIN["file"]` はファイルを SADScript として読む)。
- 影響: MAIN の `TRACK` の参照出力が K64 では取れない。
- 重要度: 高 (落ちる)。
- pysad での扱い: ActTra.f → track.f → tracka.f を直訳して [`pysad/maintrack.py`](../pysad/maintrack.py) に実装 (EMIT は古典の意味)。参照は「同じ手順を SADScript で書き追跡だけ TrackParticles に任せた版」(`bench_level27a.sad`, `27b`) との照合で代用 (座標 ≤5e-11)。

### B02 ANGLE=0 かつ K1≠0 かつ K0≠0 の BEND で FFS 初期化が SIGSEGV

- 現象: `BEND` に `ANGLE=0` かつ `K1≠0` かつ `K0≠0` を同時に与えると、FFS の初期化 (`USE RING` 後の最初の CALC の前) で SIGSEGV (exit 139)。L の値、ROTATE の有無に依らない。K1 のみ (K0=0)、K0 のみ (K1=0) は動く。
- 再現: `bench_level25e` の初版 (tmp25/bis3.py)。[`docs/history/work_rotate.md`](history/work_rotate.md) §4、[`docs/history/work_bend_k0.md`](history/work_bend_k0.md) §7。`BEND` に `ANGLE=0 K0=… K1=…` を同時に与えた LINE を `USE` する (最小スクリプトは未保存)。
- 期待される動作: tbend.f:565-611 の分岐 (ANGLE=0 で K1=0 なら tsteer、K0=0 なら tquad、他は tmulti) に従って計算される。
- 原因: 未調査 (初期化時に落ちるので、tmulti 系へ分岐する要素の設定に不備があると推測)。
- 影響: この組み合わせの要素を含むラティスは読み込めない。参照値が取れない。
- 重要度: 高 (落ちる)。
- pysad での扱い: CALC/emit は tmulte (`_bend_as_mult`) で計算、追跡は NotImplementedError。SAD との照合は不可。

### B03 SetElement を含むスクリプトの出力をリダイレクトすると SIGSEGV

- 現象: `SetElement[QF,QUAD,K1->0.3]` などを含むスクリプトを、標準出力をファイルにリダイレクトして実行すると SIGSEGV。pty 経由 (`script -qec`) で実行すると同じスクリプトが動く。
- 再現: `bench_level24f.sad` (164 タグ)。[`docs/history/work_item4.md`](history/work_item4.md) §F。
- 期待される動作: 出力先が端末かファイルかで挙動が変わらない。
- 原因: 未調査 (出力先の種別に依存するので、端末判定を使う箇所の未初期化と推測)。
- 影響: バッチ実行 (出力をファイルに取る使い方) で落ちる。
- 重要度: 高 (落ちる)。
- pysad での扱い: 参照値は pty 経由で取得 (`.result` はこの方式)。

### B04 draw で core dump

- 現象: `draw` が WSL の SAD で core dump (`dumped core`)。Tk 初期化失敗の後に SIGSEGV。
- 再現: `draw` を含む台本。
- 期待される動作: Tk/描画が使えない環境では描画をスキップするかエラーで続行する。
- 原因: 未調査 (WSL に表示環境がないための可能性が高い。環境依存)。
- 影響: 描画を含むスクリプトが最後まで走らない。
- 重要度: 高 (落ちる。ただし環境依存)。
- pysad での扱い: 描画は対象外。参照の取得時は `FFS["draw ..."]` を `Print["nodraw ..."]` に置換。

### B05 FITP>1 の GO で再起動を繰り返す

- 現象: `FITP` が 1 より大きい条件 (オフモーメンタムのマッチング、クロマティシティ補正) で `GO` すると、`*** Run time Environment` からスクリプトを最初から再実行する。`NPARA=1` でも起きる。
- 再現: `bench_level15a`。
- 期待される動作: `GO` が最後まで走って結果を返す。
- 原因: 未調査 (プロセスが異常終了して親が再起動していると見られる。実行環境の再起動機構の影響の可能性もある)。
- 影響: 最終解が取れず、比較は「再起動の途中で出た値」との照合になる。
- 重要度: 高 (落ちる)。
- pysad での扱い: match6/ffscmd の FITP>1 は直訳のまま未照合。`FITP 1` に落とした版で解の一致を確認する。

## 2. 結果が変わる誤り (重要度 高)

### B06 RADINT が全 BEND に最後に処理した BEND の L/ANGLE を使う

- 現象: RADINT の BEND の rhob/rho0 が、その BEND 自身の L/ANGLE ではなく、直前に tbend 系が設定した値 (ビーム順で最後に処理した曲がりのある BEND の L/ANGLE) になる。
- 再現: `bench_level29a.sad`、`29a2`。別の値 (BPI7 の 37.6、BID1 の 77.5) を仮定すると I9 が 8e-4〜1e-3 ずれる。20.757 m で SAD と 7 桁一致するので裏づけられた。
- 期待される動作: 各 BEND の L/ANGLE を使う。
- 原因: INTGRL.f が呼ぶ tbendebody は bendib モジュールの rhob/rho0 を使うが、INTGRL は設定しない (tbendiinit/tbendeinit を呼ばない)。
- 影響: I1X, I2, I3, I4, I5, I9, I10 など。RADINT の値が BEND ごとの正しい値からずれる (docs/18 §11 の RADINT と Emittance[] の 13% 差の候補の 1 つ。因果は未検証)。
- 重要度: 高。
- pysad での扱い: **v2 (既定)**: その BEND 自身の L/ANGLE を使う ([`test_v2_A_radint.py`](../pysad/tests/test_v2_A_radint.py))。**k64 モード (= v1)**: 再現 (`radint.default_stale_rho(els)`、`compute(stale_rho=)` で上書き可)。

### B07 RADINT の6 極ブロックが直前の DL を使う

- 現象: 6 極のブロックが DL を設定せず、直前の BEND/QUAD の最後のスライス長 DL を使うので、I9, I10 (=γ·DL) に6 極でも DL が足される。DKL2 は代入されず6 極の CKL0X/CKL1X は 0。L=0 の6 極は TCAL の後で `DL=0.2` (単精度の 0.2) にして ` L(SEXTUPOLE) EQ 0. SET L=0.2` を出す。DL, TETA, R0/R0I, FUNI などモジュール変数がプロセス中残るので、2 回目の RADINT の DL の初期値も前回の値。
- 再現: `bench_level29a.sad`、`29a2`。
- 期待される動作: 6 極のスライス長 (L) を DL に設定する。L=0 は薄レンズとして扱い、定数 0.2 を足さない。
- 原因: INTGRL.f の6 極の分岐 (SEXTU 呼び出し) が DL を代入しない。
- 影響: I9, I10 (色収差の GX/GY)。6 極が多いリングで値が変わる。
- 重要度: 高。
- pysad での扱い: **v2 (既定)**: 6 極は自分の L と K2 を使い、L=0 は寄与なし。モジュール変数を毎回初期化する。**k64 モード (= v1)**: 再現 (docs/history/work_small_items.md A.2-2)。

### B08 tsolqur が最後の粒子の bzp を全粒子に使う

- 現象: 放射あり (RAD ON) でソレノイド内 (bz≠0) の四極は、最後の半ステップで `tz%tz0%bzp` (最後に tzsetparams を呼んだ粒子の値) を全粒子に使う。結果が粒子の並び順・粒子集合に依存する (差 ∝ δ_i − δ_last)。放射なしの tsolqu は粒子ループ内で完結するので起きない。
- 再現: Level 9 (ベンチ名は docs/history/00_RESUME.md の Level 9 の項) で同じ座標の粒子セットを変えると SAD 自身が 1.15e-6 変わる。
- 期待される動作: 各粒子が自分の bzp を使う。
- 原因: tsolqu.f tsolqur の tzsetparams (:268) と bzp の使用 (:273-)。粒子ループの外で tzsetparams を呼んだ後の値を使っている。
- 影響: 粒子の並び・数で追跡結果が変わる。DA 計算の粒子スロット順にも影響する。
- 重要度: 高。
- pysad での扱い: **v2 (既定)**: 粒子ごとの BZ/(1+δ) を使う。**k64 モード (= v1)**: 再現 (track6 を要素外側・粒子内側に組み直し、`_tsolqu_bz_final` に最後の粒子の bzp を渡す)。SAD との比較は同じ粒子セットで行う。

### B09 tchge が回転行列を上書きし一周行列に回転が入らない

- 現象: tchge の入口で th≠0 かつ回した座標系のずれ≠0 (または dz, dtheta≠0) のとき、tbrote/tbgrote/tdrife1 が回転行列 `trans1` を「上書き」する。COD は回るが一周行列に回転が入らない。出口は tmultr5 で掛けるので回転が入る。一周行列が R^T·M になり、結合 (R1..R4)・チューンが物理とずれる。phig≠0 では tdrife1 (BEND 中点のドリフト) が行列に入る。四極・MULT の tsolrote にはない。
- 対象: BEND (ROTATE と DX/DY/DROTATE)、薄い多極 (SEXT/OCT/DECA と L=0 の四極)。
- 再現: `bench_level25f.sad` (25g は docs/11 A25 の記載のみ)。
- 期待される動作: 入口でも行列に回転が掛かり、一周行列が R_out·M·R_in になる。
- 原因: tchge.f の `transa=>trans1` (:45) と tbrote への渡し方 (:64-68, :83)。入口の呼び出しが trans を代入で上書きし、掛け合わせない (出口は :113-114 で tmultr5 により掛ける)。
- 影響: ロールした BEND を含むリングのチューン・結合が変わる。この癖のある要素を含むリングの Emittance[] は SAD の tcod の収束で決まる。
- 重要度: 高。
- pysad での扱い: **v2 (既定)**: 入口の回転を行列に掛ける (`_tchge_jac` が `_tchge_cod` の数値微分と 1e-8 で一致)。**k64 モード (= v1)**: 再現 (`lattice._tchge_jac`、`_thin_cod_matrix`、emit6d `_tbende`/`_tthine_gen`)。ただしこの要素を含むリングの Emittance[] (δ、チューン) は pysad と一致しない。

### B10 tdusrn が NaN の anp で 10000 回ループして乱数を消費する

- 現象: 粒子が落ちる周は発散中の値が NaN/inf のまま tradkfn が光子の乱数を引く。anp=NaN (g=NaN など) だと tdusrn が npmax=10000 回ループして 20001 個の乱数を引く (anp<0 なら 1 個)。ftest で tdusrn(NaN) と tdusrn(+Inf) は 10000 回、tdusrn(−Inf) は 0 回を確認。
- 再現: `bench_level16f`、`16l`、`20` (small)。[`docs/history/work_loss_rng.md`](history/work_loss_rng.md)、`history/work_loss_rng2.md`。
- 期待される動作: 非有限の anp では光子を引かない (あるいは損失粒子には放射を呼ばない)。
- 原因: drndsr.f :71-111 の tdusrn。`p .le. pn` と `pn .gt. cutoff` が NaN で偽になるので exit せず、npmax (drndsr.f:3、10000) まで回る。呼び出し側 (photons.f tradkfn) が損失判定前の粒子にも呼ぶ。
- 影響: 損失粒子が出た周以降の乱数列 (seed) が大きくずれ、他の粒子の放射も変わる。DA 計算 (DynamicApertureSurvey) の再現性にも関わる。
- 重要度: 高 (結果が変わる)。ただし NaN の扱いは B17 (gfortran -Ofast の MAX/MIN) と連動する未定義動作的な面がある。
- pysad での扱い: **v2 (既定)**: anp が NaN・±Inf なら光子を出さず、乱数は 1 個だけ引く。**k64 モード (= v1)**: 再現 (`Prog.la_count`、NaN の特別扱い)。再同期テストで 30 周 seed 一致。発散中の分岐は丸めの 1 ulp で逆になるので、自由走行の最終 seed は SAD と一致しない (SAD 自身も 1e-16 の摂動で変わる)。

## 3. 特定条件で結果が変わる誤り (重要度 中)

### B11 4D の CALC が縦方向の Twiss 列を作り直さない

- 現象: `CALC4D;` の後の CALC で、縦方向の列 (AZ BZ NZ ZX ZPX ZY ZPY) が作り直されず、直前の 6D の CALC の先頭行の値が全行に残る。`CALC6D;`/`CALC4D;` 自体はフラグを変えるだけで再計算しない (これは仕様)。
- 再現: `bench_level29c.sad`。6D の CALC の後に 4D の CALC をして `Twiss["AZ",...]` などを全行で読む。
- 期待される動作: 4D では縦方向の列を 0 にする (または NaN など明示的な値にする)。
- 原因: tffscalc.f:87-90 は 4D では 15 列だけ入れ直し、qtwiss1 が縦の列を直前の行から写す。
- 影響: Twiss[] の縦の列を読むスクリプトが古い値を読む。
- 重要度: 中。
- pysad での扱い: **v2 (既定)**: 縦方向の Twiss 列を既定値 (AZ=0、BZ=1、他 0) に戻す。**k64 モード (= v1)**: 再現 (`FFS._z_stale`)。

### B12 RADINT の定数が単精度

- 現象: `COEFF=(elradi/(elmass**3)/3)*1e27` の `1e27` は単精度 (9.99999988e26) なので ELOSS が相対 1.15e-8 小さい。`FREQ*1.E-6`、`DL=0.2`、`EXP(2.3026*TAUQ)` も単精度。
- 再現: `bench_level29a.sad`。SAD の ELOSS (678.50075327596767) が 17 桁で一致して確認。
- 期待される動作: `1d27`、`1.d-6`、`0.2d0` と倍精度で書く。
- 原因: INTGRL.f :35 (1e27)、:172 (DL=0.2)、:212 (FREQ*1.E-6)、:266 (2.3026)。
- 影響: ELOSS、減衰時間、エネルギー幅など ELOSS から作る量が 8 桁目で変わる。
- 重要度: 中。
- pysad での扱い: **v2 (既定)**: 倍精度の定数を使う。**k64 モード (= v1)**: 再現。

### B13 RADINT の CAVI が最後のもので上書き

- 現象: CAVI の VOLT/HARM/FREQ を合計せず最後に読んだ空洞の値で上書きする。HARM=0 なら FREQ·RLT/c、FREQ=0 なら HARM·c/RLT。
- 再現: `bench_level29a.sad`。空洞が複数あるラティスで RF 電圧・シンクロトロンチューン・バケット高さが総和にならない。
- 期待される動作: 電圧を総和し、周波数/ハーモニック数を整合させる。
- 原因: INTGRL.f の CAVI の処理 (代入であって `+=` でない)。
- 影響: 複数空洞の RADINT の RF 関連量。
- 重要度: 中。
- pysad での扱い: **v2 (既定)**: VOLT は総和、FREQ・HARM は電圧のある最初の空洞 (違う周波数が混ざると注意を出す)。**k64 モード (= v1)**: 再現 (parser が HARM から FREQ を補う分は SAD の定義値に戻して判定)。

### B14 tsolqum の z の dv 項が tsolqu と違う

- 現象: tsolqu は `−0.25·(… + ap·aln0 + …) − dv·aln0`、tsolqum は `−0.25·(… + (ap − dv)·aln0 + …)` で dv に 0.25 が掛かる。δ=1e-3 で 1 周 5e-10 の z のずれ。どちらが正しいかは物理から決める必要があるが、少なくとも両者は一致すべき。
- 再現: `bench_level11.sad` (観測 6.7e-10)。
- 期待される動作: 同じ量を計算する 2 つのルーチンで dv 項が同じ。
- 原因: tsolqu.f :86-88 と :250。
- 影響: 放射なし MULT 追跡の z。
- 重要度: 中。
- pysad での扱い: **v2 (既定)**: tsolqu と同じ −dv·aln0 にする。**k64 モード (= v1)**: 再現 (`_tsolqu(qum=True)`)。

### B15 UnstableLevel の綴り誤りで常に 0

- 現象: マッチング終了後に `UnstableLevel` が更新されず常に 0。`MatchingResidual=10^n·nstab+r` には nstab が入る。
- 再現: `bench_level28*.sad` (不安定のまま終わるマッチング)。
- 期待される動作: nstab (不安定な点の数) が `UnstableLevel` に入る。
- 原因: tffscalc.f:58 `isl = ktfsymbolz('UnstableLevel',14)-4`。名前は 13 文字なのに 14 を渡しているので、別のスロットに書かれる (`rlist(isl)=nstab` が UnstableLevel に届かない)。`MatchingResidual` (16) と `NetResidual` (11) は文字数が合っている。
- 影響: スクリプトから不安定レベルが読めない。
- 重要度: 中。
- pysad での扱い: **v2 (既定)**: 不安定な点の数 nstab を返す (MatchingResidual = 10·nstab + NetResidual)。**k64 モード (= v1)**: 再現 (`ffscmd.py:2141-2143`、`UnstableLevel=0` 固定。詳細は docs/history/work_unstable.md)。

### B16 CalculateOptics の stabz と idim=3 の stab/trace が未初期化

- 現象: `CalculateOptics[lb,le,inicond,cell,idim]` の戻りの `{stabx,staby,stabz}` の stabz が呼び出しごとに 0/1 が不定 (部分区間 cell=0 では 0 になった)。idim=3 の stab/trace も未初期化 (`{0,0,0} {1,0,0}` が返った)。
- 再現: `bench_level24d.sad`、`24d2` (docs/history/work_item4.md §D)。
- 期待される動作: 6 次元でないときは stabz を 1 か 0 に決めて返す。
- 原因: tfoptics.f → qcell1。`ffs_stat optstat` が intent(out) の局所変数で、縦方向 (z) の安定判定が代入されない経路がある。
- 影響: 戻り値を使う判定が不定。
- 重要度: 中。
- pysad での扱い: **v2 (既定)**: idim=3 は tracez = 2cos(μz) と stabz、idim=2 は stabz=1・tracez=0 を返す。**k64 モード (= v1)**: 1 を返す。ベンチでは stab は先頭 2 要素だけ、idim=3 は stab/trace を除いて比較。

### B17 gfortran の MAX/MIN が NaN を無視する前提のコード

- 現象: `-Ofast` の gfortran では `max`/`min` の NaN の返り値が箇所ごとに違う (tradkfn の `max(gmin, g+dg)` は NaN で gmin、asinz の `max(-1,x)` は −1、sqrt1 の `max(1e-200,·)` は後者)。tffsmatch の fuzz は NaN でモードを変えない (残差が悪化すると crate<0 → log10 が NaN → min/max が NaN を無視して badcnv=0 → DESCEND に切り替えず fmincube の直線探索に行く)。比較演算も NaN で真になるものがある (`v==v` が T、`v/=0` が F)。
- 再現: `bench_level16f/16l/20`。
- 期待される動作: NaN で振る舞いが言語実装に依存しない (NaN を明示的に判定する)。
- 原因: photons.f tradkfn、tbend.f tbendcore の asinz、tffsmatch.f:255-289 の fuzz。コンパイラの未定義動作への依存。
- 影響: 損失粒子の乱数消費 (B10)、マッチングの途中経過 (再現しないと早く "Residual =" で止まる)。
- 重要度: 中 (バグというより未定義動作依存)。
- pysad での扱い: **v2 (既定)**: 再現 (乱数の消費だけ v2 で NaN を明示的に扱う: B10。asinz・gfloor・sqrt1 は変更なし)。**k64 モード (= v1)**: 再現 (match6 は `_log10` が NaN を返し fuzz が 1。track6n は `sqrt1`/`pxy2dpz`/`asinz`/`gfloor`、NAN_MODE)。

### B18 PolynomialFit が平均 0 のとき xl=1

- 現象: 平均 0 (xa=0) の x に対して `PolynomialFit` が xl=1 にして係数の次数補正が壊れる。クロマティシティ補正の DP0 スキャンで 1e10 級の係数になる。
- 再現: DP0 を ±対称に振って tune を `PolynomialFit` に渡す。
- 期待される動作: xa=0 でも次数補正が正しい (`xl` は xa の大きさのスケールで、0 にならないように別に決める)。
- 原因: Packages/Fit.n :37-56 `If[xa, Table[...], 1]`。xa=0 だと偽 (0) として else 側の 1 が使われる。
- 影響: 係数が誤る。
- 重要度: 中。
- pysad での扱い: **v2 (既定)**: 平均が 0 のとき係数を平行移動なしの値にする (numpy.polyfit と 1e-13)。**k64 モード (= v1)**: sadscript.py が Fit.n をそのまま評価するので同じ挙動 (癖も再現)。

### B19 qdmdia の dr4 の未設定

- 現象: qdmdia.f の label 1300 に至る経路で dr4 が代入されない。
- 再現: 解析微分 (qdcell.f から呼ばれる縮退分岐) を使うマッチング。docs/history/work_analytic_deriv.md の対応表。
- 期待される動作: dr4 を設定する。
- 原因: qdmdia.f label 1300 (行番号は未記録)。ソースの読みのみで SAD の実行では未確認。
- 影響: 縮退分岐の微分が不定になりうる。
- 重要度: 中。
- pysad での扱い: dr4=0 にした。

## 4. 軽微・調査中 (重要度 低〜要調査)

### B20 Level[{1,{2,{3}}},{-1}] が {1,3} を返す

- 現象: `Level[{1,{2,{3}}},{-1}]` が `{1,3}` を返す (正しくは `{1,2,3}`)。同様に `DeleteCases[..,n]` の第 4 引数、非正方・ゼロ行列の `Inverse` (別の疑似逆/NaN) も Mathematica と異なる (後者は仕様の違いの可能性)。
- 再現: docs/history/work_item4.md (bench_level24 系)。
- 原因: 未調査 (SADScript の Level 実装)。
- 影響: リスト処理の結果。
- 重要度: 低。
- pysad での扱い: 再現せず。

### B21 tcod の trs(:,5:6)=0 と実測の不整合

- 現象: rfsw が偽 (再帰呼び) では tcod は `dcod1(5)=0`、`trs(:,5:6)=0` なので z も δ も動かないはずだが、RFSW OFF かつ RADCOD OFF の実測は z を閉じる δ になる (3C: δ=−4.85e-8)。
- 再現: Level 3C (docs/history/00_RESUME.md)、history/work_emit_unstable.md の 28 系 (`bench_level28a.sad` など)。
- 原因: 未調査 (temit の呼び出し側が cod の δ を決めている可能性)。バグと断定できないので要調査。
- 影響: RF なしのリングの Emittance[] の δ。
- 重要度: 低〜中。
- pysad での扱い: **v2 (既定)**: 近似 (v2 でも変更なし。z が閉じる δ を 1e-15 まで解く)。**k64 モード (= v1)**: tcod の直訳ではなく close_z (割線法) の近似。

### B22 MON/MCOD/TCOD が何もせず "Missing value for PM*"

- 現象: このビルドでは `MON`/`MCOD`/`TCOD` コマンドが何もせず式として評価される。`MON 0 0 0 0 PM*` は `Missing value for PM*` で停止する。
- 再現: `MON 0 0 0 0 PM*` を含む FFS 文。
- 期待される動作: モニタを定義する (古典 SAD の MON コマンド)。
- 原因: K64 の tffsa.f に MON・MCOD・TCOD のコマンドは無い (`grep` で 0 件。2026-11-16 にソースで確認)。FFS では式として評価されて何もしない。
- 影響: モニタを使う軌道補正スクリプトが動かない。
- 重要度: 低〜中。
- pysad での扱い: **v2 (既定)**: 「K64 に無いコマンド」と警告して文全体を読み飛ばす。**k64 モード (= v1)**: 無視 (ffscmd.py:1060)。

### B23 Extended Twiss の EPX の欄に EY が出る

- 現象: EMIT の報告 (temitprint) の Extended Twiss で、EPX の欄に EY の値が出る。`ZX/ZPX/ZY/ZPY/AZ/BZ/PSIZ` は正しい。`Emittance[]` の ExtendedTwissParameters は正しい (EPX は EPX)。
- 再現: `bench_level26a.sad` の報告 (R1 など)。
- 期待される動作: EPX の欄に EPX (`vout9(mfitepx)`)。
- 原因: temit.f:1010-1017 の `write(lfno,9001)` の引数並びで、`vout9(mfitzpx)` の次が `vout9(mfitey)` (本来は mfitepx)。format 9001 (:1021) の `EPX:` に対応する位置に EY が渡る。EY は後ろの `ZY: EY:` の位置にも出るので 2 回出る。
- 影響: 表示のみ。
- 重要度: 低。
- pysad での扱い: 再現 (`ffscmd._print_emit_exttwiss`)。報告の行を SAD と一致させるため。

### B24 `LINE["LENGTH","$$$"]` が narg エラー

- 現象: `LINE["LENGTH","$$$"]` が narg エラーになり、値が返らない (`LINE["LENGTH", loc]` 全般)。
- 再現: docs/history/work_calc6d.md (POS 行)、`bench_level19*.sad`。
- 期待される動作: 他の LINE[key, loc] と同様に末尾要素の位置 (または全長) を返す。
- 原因: 未調査 (LINE["LENGTH",...] の引数処理)。
- 影響: 位置が取れない。位置は要素名の n 番目から求めて回避。
- 重要度: 低。
- pysad での扱い: **v2 (既定)**: エラーのまま。メッセージに位置の取り方の案内を足した (K64 のソースが `LENGTH` を引数 1 個に限っており、仕様と判断)。**k64 モード (= v1)**: SAD と同じ narg エラーを再現。

### B25 TrackParticles の POL ON の列順がソース内で自己矛盾

- 現象: tftrack.f が列順を 2 通りに書いている。
- 再現: ソースの読みのみ。
- 原因: エラーメッセージ (:344) は `{z},{dp}[,{sz},{psis}],{flg}`、コメント (:502) は `(x, x', y, y', z, delta, phis, sz, alive)`。実コード (:455-460) は kaj(7)=atan(...)=psis、kaj(8)=zx(8)/|zx(7:9)|=sz の順で、コメントと一致する。メッセージが誤り (または逆)。
- 影響: 列順の文書不整合 (POL ON を使う場合)。
- 重要度: 低。
- pysad での扱い: POL ON の TrackParticles は実コード・コメントの列順 `(x, x', y, y', z, delta, phis, sz, alive[, lost])` で読む (`ffscmd` の TrackParticles、スピン追跡は `track6.py`。docs/36・43)。メッセージ (:344) の `{sz},{psis}` の順は使わない。

### B26 ORG の引数が式として読まれる

- 現象: `ORG B3 1.5 -2.5 0.5 10 20 30` が `1.5 -2.5` を式 1.5−2.5=−1.0 と読み、gx=−1.0 として 1 引数分欠落する。欠落分は 0 になり直前の ORG の値は残らない。
- 再現: `bench_level19c.sad` (6 通り)。
- 期待される動作: スペース区切りの引数は別の値として読まれる (または仕様として明記される)。
- 原因: tffsa.f の ORG が getvad (式として読む) で引数を取る。FFS の他のコマンドと同じ仕様なので、バグというより仕様の落とし穴。
- 影響: 負数を含む ORG の値が意図と変わる。`0-2.5` のように続けて書くか別行にする。
- 重要度: 低 (仕様の可能性が高い)。
- pysad での扱い: **v2 (既定)**: 空白で区切った `-2.5` は別の値として読む。**k64 モード (= v1)**: 再現。

### B27 SeedRandom[] の表示丸め

- 現象: SAD の seed の表示の末尾の桁が丸められ、seed の値が 1 違うように見える (単粒子 12 件のうち 2 件)。
- 再現: docs/history/work_loss_rng2.md の `l16l_check.py` (16l の単粒子)。
- 原因: 未調査 (表示の桁数)。値自体は 32 bit 整数。
- 影響: 表示のみ。
- 重要度: 低 (調査が浅く、バグでない可能性がある)。
- pysad での扱い: `SeedRandom[]` は符号なし (seed mod 2³²) で表示するよう修正済み。

### B28 未初期化変数の値が出力に出る

- 現象: (a) EMIT の報告の Closed orbit の Entrance の値 (閉軌道が見つからないとき) が未初期化のゴミ。(b) MAIN の TRACK の 1 周目の最初の判定の `l_track` が未初期化 (0 のはず) で、そこで落ちる粒子は印が付かない (初期振幅が LOSSAMPL=1 m 超の場合だけ)。
- 再現: (a) history/work_emit_report.md、閉軌道が収束しないラティスの EMIT。(b) history/work_main_track.md §1 (B01 のため SAD では確認不可。ソースの読みのみ)。
- 重要度: 低 (表示のみ)。
- pysad での扱い: 再現せず。

### B29 QUAD に RAD=1 を書くと COD・エミッタンスが 0 になる

- 現象: QUAD の定義に `RAD=1` を書くと COD・エミッタンスが 0 になる。原因不明。
- 再現: `sadrun/c3ye_nq` (Level 3C)。
- 原因: バグではなく仕様。K64 は `krad = (RAD_QUAD == 0) .and. al /= 0` (tturn.f:360、tsol.f:677、tturne.f:572) で、要素の RAD が 0 以外ならその要素の放射を止める (RAD=1 は「この要素は放射しない」。2026-11-16 にソースで確認)。
- 影響: 切り分けには DampingRate を使う。
- 重要度: 仕様 (バグではない)。
- pysad での扱い: 仕様どおり。

### B30 MAIN の `OFF RFSW` が FFS の RFSW に届かない

- 現象: MAIN 部 (FFS; の前) に `OFF RAD FLUC COD RFSW RADCOD;` と書くと、FFS 入口の表は `RFSW … F` と出るが、FFS 内の `FLAG[{"RFSW"}]` は 1 で、TrackParticles は空洞ありの結果 (z≠0、dp が位相で変わる) になる。RAD/FLUC/COD/RADCOD/INTRA/POL は MAIN の ON/OFF が届く。
- 再現: `<repo>/<work-item>/mflag1〜6.sad` (MAIN の ON/OFF の 6 通り → FFS で FLAG[] を表示)、`zdp_chk2_full.sad` (FFS 内の `NORFSW;`/`RFSW;` で TrackParticles の z/dp が変わる)。
- 原因: track.f:80-107 で `isynch=IgetGL('$RFSW$')`、`rfsw=isynch/=0` とした後、`nturn<0 .and. np0==0` (FFS 起動) の枝で `rfsw=.true.` に固定して `tffs` を呼ぶ (track.f:149)。表 (track.f:138) は固定前の値を出す。
- 影響: RFSW を切りたいスクリプトは FFS 内で `NORFSW;` と書く必要がある。
- 重要度: 中 (仕様かもしれないが、表と実体が食い違う)。
- pysad での扱い: **v2 (既定)**: MAIN の `OFF RFSW`・`ON RFSW` が FFS に届く。FFS の中の `RFSW;`・`NORFSW;` が後から優先。**k64 モード (= v1)**: 再現 (MAIN の RFSW は `main_flags` に入れ、FFS の RFSW は ON で始める)。

### B31-B33 wake まわりの SIGSEGV (原因と修正の差分: `docs/history/sad_patch/wake_sigsegv_fixes.patch`)

- B31: `tffswake` は `iwl=abs(kwaketbl(1,iwp))`、`iwt=abs(kwaketbl(2,iwp))` の後に、表の有無を見ずに `nlw=(ilist(1,iwl-1)-1)/2`、`ntw=(ilist(1,iwt-1)-1)/2` を計算し `c_f_pointer(c_loc(rlist(iwl)),…)` する。表がない側は 0 なので `ilist(1,-1)` を読む。後段の `tffswakekick` は `iwl /= 0` で守られている。
- B32: `tturn` は始点の SOL の前半キックを `tsol` の前、後半キックを `tsol` の後で呼ぶ。`tsol` は領域全体を進めて戻る。領域の各要素で `if(sol)` 分岐が `go to 1020` するとき `aewak` が古い真のままで、1020 の後半キックが領域の要素ごとに呼ばれ、`nwak` が表の数を超えて `twxiwp` が表の外を読む。始点のウェイクが `nextwake` を始点の位置に残すので、`tsol` の中のウェイクは (元でも) 効かない。
  修正は、`tsol` の前に始点の後半キックと `nwak` の更新を済ませ、`aewak` を偽にすること。直した SAD は pysad の追跡と、始点だけのウェイクでも領域内との併用でも 10 桁以上一致した。
- B33: `tapert` の「落ちた粒子を報告する」条件が `idtypec(l_track) == icAprt` を Fortran の非短絡の `.or.` で評価する。最初の `tapert` は `tturn` の最初の要素の前に呼ばれ、`l_track` はまだ 0。最初から LOSSAMPL の外にいる粒子があると落ちる (wake とも SOL とも無関係。w3 の LOSSAMPL=0.0009 で再現)。
  直した SAD の出力座標が pysad と合わなかった件は 2026-11-16 に原因が分かり解消した: `tapert` は落とす印に `kptbl(i,4)=l_track` を書くので、セッション最初の TrackParticles (l_track=0) の最初の `tapert` では印が 0 (生きている) になり、最初から LOSSAMPL の外の粒子は落ちず (座標は 3·LOSSAMPL に頭打ちされて追跡を続け)、後の `tapert` で落ちる。2 回目の TrackParticles からは最初の要素で落ちる。pysad は `TrackRun.l_track0` でこれを写し、gs (1・2 回目) と直した SAD に全桁一致 (test_small_items.py)。

### B63-B71 後から見つかった癖・疑い

B63〜B70 は [`docs/history/work_unimpl.md`](history/work_unimpl.md) と docs/11_sad_quirks.md (A42〜A48) に経緯がある。B71 は gs の実測 (下)。「確定」は gs の実測、「疑い」は原因・意図が未確認のもの。

- B63: 確かめた条件は、PLANAR のリングに SOL (BOUND≠0、L=0.3・0.6、BZ=0) を足した周長の比較 (偏極時間が 15 桁一致。docs/29・docs/44)。ほかの磁場・エネルギーでの確認はしていない。GEO のない SOL の S 座標の途切れは、pysad が GEO を見ないので再現しない (docs/29)。
- B64: 疑いのまま。pysad は K64 と同じに写したが、gs と差の出る例が無く (試験のリングでは効かない)、この `merge` が意図かどうかと、結果に効くかは未確認。
- B65: RADLIGHT の「記録しない」のは tlstore の仕様 (型 < CAVI の要素の後だけを記録)。意図なら癖、そうでなければ誤り。
- B66: init.n の AutoLoad の行が指す `WAKECOMMAND.n` が配布に無いこと自体は確認済み (ファイルが無い)。
- B67: 仕様と考えられるが、1 GeV と思い込んだ文書が過去にあったため載せた。
- B68: 現在の gs が RAD ON・BZ≠0 の領域内の四極の力を入れていない疑いが強い (証拠は local/wake/verify_log.md の bench_level9 の項。単一粒子 1 周追跡で、領域内がドリフトだけなら gs=pysad、四極ありは BZ≠0 で gs が四極なしの結果と一致)。ソースのどこかは未特定、保存した参照 (10-14 の gs) を出した版は今は動かず確認できない。
- B69: どちらが正しいか決められない。pysad 側の誤りの可能性も残る。pysad と gs が合わない未解決の差のうち、B68・B69 以外は、SAD の誤りか pysad の誤りかも未判定なのでこの表に載せていない (docs/history/work_unimpl.md)。v2.0.0 のときの 3 件のうち、BEND の F1=0.1 FRINGE=3 の RAD ON 追跡の差 (B113) と L≠0 の CAVI の DX・DY の差 (CAVI の TENT・TEXIT が pysad に無かった) は v2.1.0 で解決した (gs と全桁一致、`history/v2next_A1.md`)。K3 が強い MULT リングの CALC (qcod) の閉軌道は未解決のまま。
- B70: gs で `MessageName[Sin,"usage"]` が `Undefined MessageName General::"usage"` を返すことを 2026-11-17 に再確認 (原因は未調査)。
- B71: gs で `For[i=0,i<3,i++,Print[i]]` が何も印字せず `i` は 3 になることを 2026-11-17 に再確認 (pysad は 0 1 2 を印字して i=3)。`Packages/functions.n:34` の定義は本体を実行する形なので、組み込みの `For` が先に評価されると思われるが、原因は未特定。

### B72-B162 v2.0.0・v2.1.0・v2.2.0・v2.3.0 の作業で新しく見つけたもの

B72〜B162 は表 (§0) に症状・場所・影響・重要度・確度・pysad の扱いを書いた。詳細 (根拠の数値、試験、使った条件) は次の記録にある (B111〜B134 は v2.1.0。記録の項目名との対応は `local/wake/v21_bnum.md`)。

- B72〜B77 (RADINT・減衰・運動量圧縮): `history/v2_changes_A.md` §1、B75 は §2。B106〜B110 は同 §4 の「要判断」。
- B78〜B80 (UND・PHSROT・ソレノイド + 双極の放射): `history/v2_changes_C.md` §1 の「新しく見つけたもの」(N1・N4・N5)。
- B81〜B86 (スピン追跡の縦磁場・刻みの角・wake の重み・AUTOPHI・trots): `history/v2_changes_C2.md` §1・§3 (N6〜N8、N10、N11、N13)。
- B87・B88 (MEASURE の初期分布): `history/v2_changes_C_measure.md` §1・§6。B100〜B105 (TRCOD・ALIGN・GETMEASURE・DRAW の旧いコード): `history/v2_changes_C2_measure.md` §8。
- B89〜B98 (空間電荷 SPAC・PSPAC): `history/v2_changes_C_spch.md` §2 (SPC1〜SPC9)、`history/v2_changes_C2_spch.md` §1・§7 (SPC3 の補足、SPC10)。
- B99 (DELCOR の規格化): `history/v2_changes_C2_misc.md` §2。作者の意図の裏づけは無い (推定)。

- B111〜B134 (v2.1.0): `history/v2next_*.md` の各担当の「K64 の誤り・新しく見つけたもの」。B111・B112・B113 は A1 (項目 7・2・4)、B114・B115・B116 は A2、B117 は B1 (項目 2)、B118〜B122 は B2、B123・B124 は B3、B125 は B3 と D1、B126 は C (B107 の 2 つ目)、B127・B128・B132 は D1、B129 は D2、B130・B131 は E、B133・B134 は B1。B40・B60・B81・B83・B99・B106〜B109 の更新は C・D2・E・A2。
- B135〜B159 (v2.2.0): `history/v2next_P1.md`〜`P5.md` の各担当の「K64 の誤りの候補」。B135・B136・B137 は P1、B138・B139 は P2、B140〜B145 は P3、B146〜B150 は P4、B151〜B159 は P5。B42・B108・B111・B112・B116・B128・B129 の更新は P5・P3・P2・P1。B128 は P2 が原因を特定した (SB の式の誤り) ので新しい番号は付けていない。P4 の「DAPWIDTH の窓を pysad が追っていなかった」は pysad 自身の不具合で、B 番号を付けていない (docs/47 (c))。
- B160 (v2.3.0): `history/v23_T1.md` §3.2・§6 (K1≠0 の BEND のオフモメンタムの光学。gs・解析式・elegant との比較、[`example/v2_3/T1_gs_bend_k1.py`](../example/v2_3/T1_gs_bend_k1.py))。U1 で K1 の暗黙の定義 (誤りではない) に再分類: `history/v23_U1.md`。
- B161 (v2.3.0): `history/v23_U1.md` §8 (CALC の `tbendecorr` の dk と kick の導関数の xr² の項の不整合。疑い)。
- B162 (v2.3.0): `history/v23_T2.md` (Touschek の損失率の式 `touckf` の適用範囲。独立な式 `touschek_piwinski` との比較、elegant との比較、`example/v2_3/T2_*`)。

確度が「疑い」のもの (B88、B99 の意図、B113 の意図、B131) と、ソースの読みだけのもの (gs で再現できていないもの) は、K64 の作者の意図を確かめていない。v2 で直したものも、直したのは pysad の v2 モードだけで、K64 のソースは変えていない。

## 5. 開発者への依頼の要点

優先して直してほしい順に:

1. B01, B02, B03 (落ちる。再現が最小で、原因がソースから追える)。
2. B06, B07, B12, B13 (RADINT。INTGRL.f の中で閉じていて、修正が小さい)。
3. B08, B09, B10, B14 (追跡・行列の結果が変わる)。
4. B15, B16, B23 (綴り・未初期化・引数の並び。修正は 1 行)。

なお、v1 系の pysad は「SAD の数値を再現する」ことが目的なので、上記のうち結果が変わるものは修正されても再現を外す必要がある (対応箇所は各項の「pysad での扱い」)。v2 系の pysad は互換モードで両方を持つ (`PYSAD_COMPAT=k64` が K64 の動作、既定の `v2` が修正した動作)。

## 6. 本レポートに含めなかったもの (参考)

- pysad 側の誤りだったもの: fmincube の `m_euler` (SAD はオイラー定数 γ、pysad がネイピア数 e を入れていた)、MASS の既定値、`MINCOUP` の既定値、`Cavity.PHI` の既定値、`lattice._tbendebody0_mat` の z 行の符号。
  v2.0.0 の作業で見つかった pysad 自身の不具合 (MULT のスキュー成分、`CalculateOptics` の TRPT、加速区間の分散、`LegendreP` の符号、EMITDIV ほか) は [`47_v2_changes.md`](47_v2_changes.md) の (d)。
- 仕様・癖: docs/11_sad_quirks.md の A1-A48 のうち上記 (B63〜B71) 以外、B1-B5、C1-C23、D1-D5。
