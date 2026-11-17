# 21. SAD K64 のバグ報告 (pysad 移植中に見つかったもの)

最終更新: 2026-11-17 (現状のスナップショット。B63〜B71 は 2026-11-17 までに見つかった癖・疑いの追加。B01〜B62 は 11/06〜11/15 のまま)

> 対象: SAD K64 (github KatsOide/SAD、WSL Ubuntu でビルドした `bin/gs`)。ソースは `<SAD>/src`、`Packages/`。
> 材料: docs/11_sad_quirks.md、docs/history/work_*.md、docs/history/00_RESUME.md。行番号は 2026-10-13〜16 の照合時点のソース。
> 「仕様・癖」(E1/E2 が ANGLE 比、K1 が積分値、`FIT;` の意味、NPARA の既定と fork、LINE[] のファミリー展開など) は含めない。
> 明らかな誤り (綴り・引数の並び・未初期化・クラッシュ・同種ルーチン間の不整合) だけを選んだ。
> 「再現手順」の `bench_levelXX.sad` は `<repo>\` にある。「pysad での扱い」の「再現」は SAD に数値を合わせるため誤りもそのまま写したもの、「正しい動作」は pysad が正しい挙動にしたもの。
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
| B06 | RADINT: 全 BEND に「最後に処理した BEND」の L/ANGLE を使う | INTGRL.f (tbendebody の bendib) | I1..I10 などが変わる | 高 | 確定 | 再現 |
| B07 | RADINT: 六極ブロックが直前の BEND/QUAD の DL を使う (L=0 の六極は DL=0.2) | INTGRL.f TCAL | I9, I10 が変わる | 高 | 確定 | 再現 |
| B08 | tsolqur (RAD ON, bz≠0) が最後の粒子の bzp を全粒子に使う | tsolqu.f tsolqur :268-273 | 粒子の並びで結果が変わる | 高 | 確定 | 再現 |
| B09 | tchge の入口で tbrote 等が回転行列を上書きし、一周行列に回転が入らない | tchge.f :45, :64-68 | 結合・チューンが物理とずれる | 高 | 確定 | 再現 |
| B10 | tdusrn が NaN の anp で npmax=10000 回ループし乱数を大量に消費 | drndsr.f :71-111 | 損失粒子の周以後の乱数列がずれる | 高 | 確定 | 再現 |
| B11 | 4D の CALC が縦方向の Twiss 列 (AZ BZ NZ ZX ZPX ZY ZPY) を作り直さない | tffscalc.f :87-90、qtwiss1 | Twiss[] に古い値が残る | 中 | 確定 | 再現 |
| B12 | RADINT: 定数が単精度 (`1e27`, `1.E-6`, `0.2`, `2.3026`) で ELOSS が相対 1.15e-8 小さい | INTGRL.f :35, :172, :212, :266 | 数値が 8 桁で変わる | 中 | 確定 | 再現 |
| B13 | RADINT: CAVI の VOLT/HARM/FREQ を合計せず最後のもので上書き | INTGRL.f | 空洞が複数あると誤り | 中 | 確定 | 再現 |
| B14 | tsolqum の z の dv 項が tsolqu と違う (dv に 0.25 が掛かる) | tsolqu.f :86-88 と :250 | δ=1e-3 で 1 周 5e-10 | 中 | 確定 (実測 6.7e-10) | 再現 |
| B15 | `UnstableLevel` が綴りの誤り (文字数 14) で常に 0 | tffscalc.f :58 | 不安定レベルが読めない | 中 | 確定 | 再現 |
| B16 | CalculateOptics の stabz と idim=3 の stab/trace が未初期化 | tfoptics.f → qcell1 の `optstat` | 戻り値が不定 | 中 | 確定 (実行) | 1 を返す |
| B17 | gfortran の MAX/MIN が NaN を無視する前提のコード (tradkfn、asinz、fuzz) | photons.f、tbend.f、tffsmatch.f :255-289 | NaN の扱いが箇所ごとに違う | 中 | 確定 | 再現 |
| B18 | PolynomialFit が平均 0 (xa=0) のとき xl=1 になり次数補正が壊れる | Packages/Fit.n :37-56 | 係数が誤る | 中 | 確定 | 再現 |
| B19 | qdmdia の label 1300 で dr4 が未設定 | qdmdia.f | 縮退分岐の微分が不定 | 中 | ソース | dr4=0 にした |
| B20 | `Level[{1,{2,{3}}},{-1}]` が `{1,3}` を返す | SADScript の Level | リスト操作の結果が誤り | 低 | 実測 | 再現せず |
| B21 | tcod の rfsw=.false. 枝の `trs(:,5:6)=0` と実測 (z を閉じる δ) が合わない | tcod.f | RFSW/RADCOD OFF の δ | 低〜中 | 実測 (原因未調査) | 近似 |
| B22 | MON/MCOD/TCOD が何もせず、`MON 0 0 0 0 PM*` が "Missing value for PM*" | FFS コマンド表 (K64 の tffsa.f に MON・MCOD・TCOD のコマンドは無い。ソースで確認 2026-11-16) | 補正スクリプトが止まる | 低〜中 | 実測 + ソース | 無視 |
| B23 | Extended Twiss の EPX の欄に EY が出る | temit.f :1010-1017 | 表示のみ | 低 | 確定 | 再現 |
| B24 | `LINE["LENGTH","$$$"]` が narg エラー | LINE[] の引数処理 | 位置が取れない | 低 | 実測 | 再現 |
| B25 | TrackParticles の POL ON 列順がソース内で自己矛盾 | tftrack.f :344 と :502 | 文書・メッセージの不整合 | 低 | ソース | 未対応 (POL) |
| B26 | `ORG` の引数が式として読まれ `ORG B3 1.5 -2.5 ...` が `1.5-2.5` になる | tffsa.f ORG、getvad | 負数の書き方で値が変わる | 低 | 確定 | 再現 |
| B27 | SeedRandom[] の表示と seed の末尾桁の丸め | 表示ルーチン | 表示のみ | 低 | 実測 (調査浅い) | 符号なし表示に修正 |
| B28 | 未初期化変数の値が出力に出る (EMIT 報告の Closed orbit Entrance、MAIN TRACK の l_track) | temit.f、tracka.f | 表示のゴミ | 低 | 実測 | 再現せず |
| B29 | QUAD に RAD=1 を書くと COD・エミッタンスが 0 になる | tturn.f:360、tsol.f:677、tturne.f:572 (`krad = (RAD_QUAD == 0) .and. al /= 0`) | 結果が変わる。**バグではなく仕様** (RAD=1 はその要素の放射を止める。tturn.f:360 ほか) | 仕様 | 確定 (ソース) | 仕様どおり (DampingRate などで切り分ける) |
| B30 | MAIN の `OFF RFSW` (`$RFSW$`) が FFS の RFSW に届かない (他の RAD/FLUC/COD/RADCOD/INTRA は届く) | track.f :107, :149 (`rfsw=.true.` で tffs を呼ぶ) | `OFF RFSW` を書いても空洞が効く。FFS 入口の表の RFSW 列は MAIN の値で FFS の値と違う | 中 | 確定 (実測 + ソース) | 再現 (`main_flags`) |
| B34 | BEAMBEAM の傾き回転 (tilt) が cost=sin(傾き角)、sint=傾き角を使う (bbinit は xangle(6)=sin、(7)=cos を書く) | bbstrhl.f:beambeam (sint=colb%xangle(5)、cost=colb%xangle(6)) | 強ビームに x-y 結合があると座標がほぼ 0 に潰れる | 高 | 確定 (gs で再現) | SAD のまま (直訳) |
| B35 | TRPT で SOL 領域 (SOL BOUND=1 の間) に DRIFT/BEND/QUAD/SOL/MULT/CAVI/MAP/APERT/MARK/MON 以外の要素 (SEXT、COORD、BEAMBEAM) があると SIGSEGV (本来は 'Only DRIFT, ... are supported within SOL' を出して無視する分岐) | tsol.f:798 (case default 付近) | 落ちる | 中 | 確定 (デバッグビルドで特定、メッセージの write を消しても落ちる。原因は未特定) | NotImplementedError |
| B36 | L≠0 の TCAVI を含むリングの Emittance[] で EnergyLossU0 が NaN、縦チューン 0、エミッタンスが桁違い (K0=0 でも。ドリフトなら正常) | ttcave.f (tdrife0 に krad を渡す第 2 ドリフト付近。RAD=1 でも NaN なので原因未特定) | 結果が壊れる | 中 | 確定 (gs で再現)、原因は未特定 | pysad は自然な値 (L=0 は gs と一致) |
| B37 | UND に磁場のキーワードが無い (FBX FBY BX BY はどれも Undefined keyword)。undulator は磁場なしのドリフトになる | initb1.f:1208-1209 (`kytbl(kwFBX,icUND)=…`。kwFBX が未定義の変数で 0 になり、kytbl(0,icUND) を上書きするだけ) | 機能が使えない | 低 | 確定 (gs で再現) | SAD のまま (FBX 等は NotImplementedError) |
| B38 | PHSROT 要素は動作しない: phsinit の行列 (PtoN など) がどんな入力でもすべて 0 になり、追跡で粒子の座標が 0 になる。さらに PX PY PZ (位相) のキーワードは名前表に未登録で指定できず、`i_p0` が未定義の整数変数 (=0) でその dv の更新は blist(1)・blist(2) を使う | phsrot.f (phsinit、phsrot)、initb1.f (kwPX 等) | 機能が使えない | 低 | 確定 (gs で再現、行列の原因は未特定) | 追跡は NotImplementedError。Twiss・Emittance[] は恒等 |
| B39 | tradke (ソレノイド内、bzh≠0) の dbeam(8) (px・py の共分散) が `beam(2)-beam(10)+b*beam(4)` と書かれている。力学的運動量 (px+b·y, py−b·x) の共分散なら `beam(2)-beam(9)+b*beam(4)` (B12 − B34 + b B13) のはずで、beam(10)=B44 は誤記に見える | photons.f tradke | SOL 領域内の放射のビーム行列 (px-py 結合) | 低 (bzh=0 では効かない) | 疑い (gs と合わせるため pysad はそのまま写した。beam(9) に直しても Emittance[] の ε は強い場の試験 (BZ=5、軌道 2.7 mm) でも全桁同じで、試験では区別できない) | `perfkern.tradke_db` にコメント |
| B40 | DELCOR (tcorr.f) が MARK に `latt(j)+15`、`+17` を書く。値の添字では ky_DETR_Mark=15、ky_DPX_Mark=17 (DX=16、DY=18) で、DX・DY の一つ手前を書いている。CAVI も `le=latt+3` に書いて元の値は `idvalc+13` から読む食い違いがある | tcorr.f | MARK への DELCOR は初期条件の DETR・DPX を変える (DX・DY は変わらない)。CAVI は不定 | 低 | 確定 (MARK は gs で再現、CAVI は読みのみ) | pysad は MARK を同じ欄に書く。CAVI は NotImplementedError |
| B41 | ttinit.f (MEASURE の初期分布) が `rlist(kx+26..29)` を軌道の揺らぎ (JDX..JDPY) として読むが、MARK の値の添字では 26..29 は ky_ZPX, ky_ZY, ky_ZPY, ky_DP。JDX.. は 34..37 | ttinit.f | MEASURE の初期の py の揺らぎが DP·tgauss になり、乱数も 4 回余分に消費する | 低 | 確定 (読解。MEASURE が出力しないため gs では未確認) | 未翻訳 (MEASURE 全体) |
| B42 | `TRCOD` が gs (この環境の K64 ビルド) で SIGSEGV。`MEASURE` は何も出力しない | trcoda.f / trackb.f | TRCOD・MEASURE の参照が取れない | 低 | 確定 (gs で再現) | 未翻訳 |
| B43 | `Hypergeometric2F1[.5,.5,1,.99]` が 2.348414 (真値 2.352716)。`BesselJ[.5,30.]` が -.143924709693 (真値 -.143929653370)。級数・漸近展開の打ち切りが粗い | gamma.f / hg.f (hgrr, cbesj = confhg0) | 引数が 1 に近い 2F1、x が大きい Bessel で相対 1e-3〜3e-5 の誤差 | 低 | 確定 (gs で再現) | pysad は scipy の値 (真値) を返す |
| B44 | `ExpIntegralEi[1.5]` が複素数 `.164-2.985 I` (真値 3.3013)、`Ei[-1.5]` が -2.0089 (真値 -.1000)、`Ei[.3]` が `-.166-1.417 I` (真値 -.3027)。D.n の式が HypergeometricU の負引数 (NaN 経路) に依存 | Packages/D.n:ExpIntegralEi | Ei が全域で誤り。LogIntegral は別式で正しい | 中 | 確定 (gs で再現) | pysad は真値 |
| B45 | `HermiteH[3,.7]` が NaN (整数 n で Gamma[(1-n)/2] が極)。`LegendreP[3,1,.4]`・`LegendreP[2,2,.3]` が NaN (Hypergeometric2F1Regularized の 1-m が非正整数)。`Binomial[-3,2]` が NaN | Packages/D.n | 整数次数で NaN | 中 | 確定 (gs で再現) | pysad は真値 (eval_hermite、lpmv、積) |
| B46 | `HypergeometricU[1,2,.5]` が INF (真値 2)、`BesselY[0,.1]` が NaN (真値 -1.5342)、`BesselK[1.5,.7]` が INF。`BesselJ[-1,1.5]`・`BesselY[1,2.5]`・`BesselJZero[0,1]` は gs が応答しない (無限ループ) | gamma.f (chgu, cbesy, cbesk, confhg0) | 整数次数 (b が整数の U、cbesy の n 整数での 0/0) | 中 | 確定 (gs で再現。応答しないのは 20 秒で打ち切り) | pysad は真値 |
| B47 | `InverseErf` は Newton の係数に sqrt(pi/2) を使い (正しくは sqrt(pi)/2)、`df^2 <= no*f1^2*eps` で止めるため、真値から相対 1e-8 ずれる (`InverseErf[.3]` = .272462711811、真値 .272462714727)。FitGaussian の出力 (σ の 8 桁目) に出る | gamma.f:inverseerf | 8 桁以降の誤差 | 低 | 確定 (gs で再現) | pysad は K64 と同じ値を返す (`specfn.inverseerf`) |
| B48 | `FixedPointList` は gs で未評価のまま返る (autoload の表に `FixedPoinstList` と綴り違い (Misc.n の SetAttributes と init.n)) | Packages/init.n:597、Misc.n:23 | FixedPointList が使えない | 低 | 確定 (gs で再現) | pysad は Misc.n の定義どおり動かす |
| B49 | `RadiationField` が gs (この K64 ビルド) で SIGSEGV (どの引数でも。tlcalc の `loc_sad(ktfaddr(ka)-3,dl)` が不正な番地を指すとみられる)。RadiationSpectrum は動く | tlinit.f:tlcalc | RadiationField の参照が取れない | 低 | 確定 (gs で再現) | pysad はソースどおりに翻訳 (`radfield.py`)。等速運動の電荷の速度場と相対 6e-3 で一致を確認 |
| B50 | `FileByteCount` は `wc -l` (行数) を返す (14 バイトの 2 行のファイルで 2)。FileSystem.n のファイル操作 (CopyFile RenameFile DeleteFile …) は gs のこの環境 (WSL の /mnt/c) で失敗が多い | Packages/FileSystem.n | バイト数が取れない | 低 | 確定 (gs で再現) | pysad は Python の os・shutil で動かし、FileByteCount はバイト数 (`fileops.py`) |
| B51 | `RenumberElement[name,{n1,..}]` が gs (この K64 ビルド) で SIGSEGV (tfrenumber_element)。参照が取れない | tfrenumber_element.f | RenumberElement が使えない | 低 | 確定 (gs で再現) | pysad はソースどおりに翻訳 (`ffscmd.fn_renumber_element`。NAME.n の番号の付け替え、足りない分は最大から続ける) |
| B52 | `Fit[data,fun,{x1,x2},...]` (ベクトルの独立変数。data の行は {{x1,x2},y}) が、式を `Function[{v},fun]@@x$[[k]]` で呼ぶ際に General::narg を出し、最小二乗解にならない値を返す (例: ChiSquare 1.1519 に対し真の最小は .9320) | Packages/Fit.n | ベクトル独立変数の Fit が使えない | 低 | 確定 (gs で再現) | pysad は正しい最小二乗解 (`fitroot.py`) |
| B53 | `BesselJ`・`BesselK`・`BesselY` が x が大きいと不正確 (x=33 で BesselJ[0,33.] が 3e-3 相対、BesselK[0,12.5] が 1e-5)。`FresnelS[6.1]` = 1.3e7、`FresnelC[6.1]` = -2.1e7、`AiryAi[15.5]` = -16 など、引数が大きい Fresnel・Airy は桁落ちでデタラメ (超幾何級数の直接和) | gamma.f (cbesj = confhg0 の級数)、D.n (FresnelS/C は HypergeometricPFQ) | 引数が x≳6 (Fresnel・Airy)、x≳12 (Bessel) の値が使えない | 中 | 確定 (gs で再現) | pysad は scipy の値 (真値) |
| B54 | `ExpIntegralEi[x]` は x<0 でも誤り (Ei[-0.1] = -2.2265、真値 -1.8229。B44 の続き)。`PolyLog[n,z]` は gs が複素数 (虚部 1e-15 の雑音つき) を返す | Packages/D.n | | 低 | 確定 (gs で再現) | pysad は真値 (実数) |
| B55 | POL の平衡偏極度 `Polarization` / `Polarizations` が NaN になる (試した 8 つのリング、平面でも非平面でも)。`spdepol` → `spintrm` の `exp(cb**2-cphi)*cerfc(cb)` が \|cb\|>18 (cezlim) で漸近式 `exp(-cb²)` → 0 と `exp(cb²)` → ∞ の積になる。cb = −(e^{−cφ}−1)/a、a = √2·ssprd(スピン軌道結合 ~0.05 など) で、cb は通常 20 以上。erfcx 相当 (`erfc·exp` を桁あふれしない形) に書き直すと有限 (4e-7〜8e-7、つまり 4e-5〜8e-5 %) になるが、これも 92% 級の理論値からは遠く、spdepol の式か単位の誤りの疑いがある | photons.f:spintrm (cexp1・cerfc)、gamma.f:cerfc | 偏極度が使えない (PolarizationTime=1/ΣT 式は有限で正しい) | 中 | NaN は確定 (gs で再現)。原因は推定 (pysad で同じ式を写すと同じ NaN、erfcx 版で有限)。有限値の妥当性は未確認 | K64 の式のまま (NaN も写す) |
| B56 | `tbrote` (tchge.f:197-202) の srot 更新は `do concurrent (i=1:9)` の中で行を逐次更新するので、2 行目は更新後の 1 行目を、3 行目は更新後の 1・2 行目を使う (rr·srot の行列積にならない)。`serot` (photons.f:850-857) も `c1` を更新してから `c2` の更新に使う。`tradke` の末尾 (photons.f:1614-1615) の `srot(3,:)=(-sp*srot(1,:)+srot(3,:))/cp` も更新後の 1 行目を使う (これは偏向角 sp≠0 の bend で効く。平面リングの gs と一致したので意図かもしれない) | tchge.f:197、photons.f:850、1614 | 偏向磁石のずれ (DX・DZ・DROTATE・CHI2) で POL の結果がずれる可能性 | 低 | ソースの読み (gs での確認は serot と tradke 末尾だけ。tbrote は未確認) | 写している |
| B57 | `tradke` の `call tadddsg(dsg,srot,dgx,dgy,dgx,dbeam)` は dgz の引数に dgx を渡す (photons.f:1589)。dsg は `spequ` (使われない旧ルーチン) でしか読まれないので結果に影響はない | photons.f:1589 | なし | 低 | ソースの読み | pysad は dsg を持たない |
| B58 | `tbendebody`・`tbendebody0` の `bsir0` の更新で `bsi2` の項は enarad が偽の最後のスライスにしか付かず、常に捨てられる (`bsi2=1` は n==n2 のときだけで、そのとき `enarad .and. n /= n2` は偽)。意図は偏向磁石の出口の bsi の補正とみられる | tbende.f:215・313、535-537 | 縦のオフセットがある偏向磁石の POL で出口側の bsi が入らない | 低 | ソースの読み | 捨てる (K64 と同じ) |
| B59 | `CSRSetupOY` (tfcsroysetup) は dpjn を zjn と同じ記憶域に割り当てる (`c_f_pointer(c_loc(rlist(kazjn)),dpjn,…)` が 2 回、dpjn 用の `kadpjn` を使わない)。`zjn(k,i)=z` の直後に `dpjn(k,i)=p` が同じ場所を上書きするので、返す 4 番目の表 (zj) は z ではなく dp の値になり、5 番目 (dpj) は書かれない (0 のまま)。続く `CSROYMatrix` は zj を位置として使うので行列が狂う。`csrfindc` は見つけた値を返した後の `csrfindc=0.d0` で常に 0 を返す (未使用)。`csrjasign` は `rho` がしきい値以上の点が無いと `i1` が未初期化 | csroy.f:tfcsroysetup (kazjn, dpjn)、csrfindc | CSROY の行列が使えない | 中 | ソースの読みと gs の出力 (dpj の和が 0、zj の値が p のスケール) | K64 と同じ (zj に p を入れ、dpj は 0) |
| B60 | `CSROYMatrix[{wi,J,omega,zj,dpj,rhoj},{dz,inten,nus,sigz,sige}]` が gs (この K64 ビルド) で、`CSRSetupOY` の出力でも手で作った実数リストでも General::wrongtype を返す。登録名は `CSRSetupOY`・`CSROYMatrix` (tfinitn.f:1044-1045。ソース内の関数名 tfcsroysetup・tfcsroymat と綴りが違う) | csroy.f:tfcsroymat | CSROYMatrix の参照が取れない | 低 | 確定 (gs で再現)、原因は未調査 | ソースどおりに翻訳 (gs 未照合) |
| B61 | `FindRoot` に複素の初期値を与えると壊れた結果になる (`FindRoot[z^2+1==0,{z,1.+1.*I}]` が `{z->Restrict[(.25+I),-Infinity,Infinity],Residual->0}`、`Exp[z]==2` も同様)。pysad は実数の `FindRoot` だけ gs と一致させ、複素は TypeError のまま (gs の壊れた値は再現しない) | K64 | — | — | — | — |
| B62 | `DELCOR` (tvcorr) の corlength≠0 の枝は `cv(n)` を設定せず、`ktaloc` の古い内容が入る (2 回目以降の DELCOR で値が数 % 変わる)。CAVI への DELCOR は DX・DY ではなく HARM・PHI に書く (tcorr.f の `latt+3`、添字のずれ) | K64 | — | — | — | — |
| B31 | 光学の wake (WAKEOPT) で、縦と横の表の片方しかない要素があると SIGSEGV | tffswake.f:102 | 落ちる | 高 | 確定 (gdb 相当: addr2line、修正して確認) | 直した SAD と同じ動作 |
| B32 | TrackParticles (TRPT、wake) で SOL 領域の始点の SOL にウェイクがあると SIGSEGV | tturn.f:567 (twxiwp) | 落ちる | 高 | 確定 | 直した SAD と同じ動作 |
| B33 | TRPT で、最初の `tapert` の時点で LOSSAMPL の外にいる粒子があると SIGSEGV | tapert.f:184 (`idtypec(l_track)`、l_track=0) | 落ちる | 高 | 確定 | 再現 (最初の tapert は l_track=0 で落とす印が 0 になる癖を `TrackRun.l_track0` で写した。直した SAD と全桁一致。出力座標の差は解消) |
| B63 | SOL (BOUND≠0) の長さが周長に入らない (omega0・MomentumCompaction・空洞の位相・偏極時間に効く)。SOL に GEO=1 がないと S 座標が途切れて周長がさらに狂う | 周長の計算 (場所は未調査) | 周回周波数・運動量圧縮・偏極時間が SOL の長さの分だけずれる。GEO なしでは偏極時間が約半分、MomentumCompaction が約 100 倍 | 中 | 確定 (gs で確認: PLANAR + SOL で、SOL を除いた周長 45.6 のときと偏極時間が 15 桁一致)。意図か誤りかは未判定 | 再現 (`emit6d._emit_circ`、追跡の `circ`)。GEO なしの S 座標の途切れは再現しない |
| B64 | `tmulke` (CALC・Emittance[] の MULT) は m=0 の枝 (L=0 の薄いキックと端) で `bsir0` を更新せず、m≠0 の枝 (スライス) では `bsir0=merge(bsir0+imag(cx)/aln,0.d0,m == 1)` と書く。`merge` は m==1 のときだけ加算し、m≥2 では 0 に戻すので、意図 (全スライスの積算) と合わないように見える | tmulte.f:6-43 (tmulke) | 縦オフセットのある MULT の POL で bsi の更新が抜ける可能性 | 低 | 疑い (ソースの読みのみ。gs の試験リングでは効かず、差が出る例が無い) | K64 のとおり (`lattice.kick_cx`、`emit6d._tmulte`)。`merge` の振る舞いも同じに写した (`emit6d.py` の tmulke 相当、m≠1 で 0 に戻す)。m=0 の枝は更新なし |
| B65 | RADLIGHT: `TrackParticles` の軌道 {t,x,y,z} に横方向の座標が入らない。tlinit が geo0 の軸 1・2 を設定せず (初期値 0)、位置は straight (dphi=0) の積み上げで BEND で向きが変わらず、BEND と曲げのある MULT の中は記録しない | tlinit.f (tlinit / tlstore / tlresult) | RadiationField の入力が進行方向の点だけになる。BEND を含む実リングでは意味のある場にならない (`example.oide.sad` は B を ANGLE=0 の STEER にして避けている) | 中 | 確定 (gs と差 0 で再現、`test_radlight.py`)。意図か誤りかは未判定。なお RadiationField 自体は gs が SIGSEGV (B49) | 再現 (`ffscmd._RadLight`) |
| B66 | 旧式の `WAKE` コマンドは `tfwake.f` を呼ばず `WAKECOMMAND[]` (init.n の AutoLoad) を評価する。`Packages/WAKECOMMAND.n` は配布に無く、fileopen のエラーを出して行の残りを捨てる | tffsa.f:976 (docs の旧記述は 972)、init.n:696 | WAKE コマンドが使えない (tfwake.f は死んだコード) | 低 | 確定 (gs で再現) | 同じ経路 (エラー文の表示だけ未再現)。WAKE は `WakeFunction` + LWAKE/TWAKE (docs/23) |
| B67 | `MOMENTUM` の初期値が 30 GeV (`RsetGL('MOMENTUM',30.d9)`)。MOMENTUM を設定しないスクリプトの追跡・Emittance[] が 30 GeV で計算される | initb1.f:1698 | MOMENTUM の書き忘れに気づきにくい | 低 | 確定 (ソースと gs)。仕様の可能性が高く「癖」 (docs/11 A43) | 再現 (既定 30e9)。以前の pysad の文書は 1 GeV と誤記していた |
| B68 | SOL (BZ≠0) 領域内の四極 + RAD の追跡で、現在の gs (2026-10-19 ビルド) が、保存した参照 (`bench_level9`、10-14 の gs) と別の値を返す (S R T1 P1 の x が -.0027867898 から -.0014539780)。NORAD は一致 | tsol 系 (ソースのどこかは未特定) | 同じ入力で gs のビルドにより RAD の追跡結果が変わる。現在の gs が RAD ON・BZ≠0 の領域内の四極の力を入れていない疑いが強い (local/wake/verify_log.md の bench_level9 の項: `SB D1 QD SE` が `SB D1 SE` と 8 桁一致。10-19 と 11-14 -O0 の 2 ビルドは同値で、最適化の差ではない) | 中 | 疑いが強い (gs 側。ソースの箇所は未特定) | pysad は 10-14 の参照・物理と一致。現在の gs とは合わない |
| B69 | SOL 領域 (BOUND≠0 の間) 内の直線の BEND (ANGLE か K0 付き) の偏向を gs がほぼ無視する (BZ=0.5 T の領域で ANGLE=0.1 でも位置・角度が 1e-5 しか変わらない)。K64 の tsol.f は `tdrift(al,bzs,phiy,phix,…)` を呼ぶので、その直訳は ANGLE·距離だけ曲がる | tsol.f、tpara (領域内の BEND への K0/ANGLE の入り方は未調査) | SOL 領域内の BEND の追跡 | 低 | 疑い (原因は未確認。K64 側の癖かもしれない)。生の gs の出力は `local/wake/ps23.sad`〜`ps28.sad` | 未対応 (NotImplementedError) |
| B70 | `MessageName[f,"usage"]` (関数形) は f を無視して `Undefined MessageName General::"usage"` になる。`f::usage` の形は SAD の MessageName.n の表を引いて本文を返す | MessageName の処理 (未調査) | 関数形の MessageName が使えない | 低 | 実測 (gs で確認、原因は未調査) | `f::usage` の形は一致。関数形は見送り (docs/20) |
| B71 | 組み込みの `For[ini,cond,incr,body]` は 4 番目の引数 (本体) を実行せず、3 番目 (incr) だけを繰り返す (`For[i=0,i<3,i++,Print[i]]` は何も印字せず i=3)。`Packages/functions.n` の定義 `For[`ini_,`cond_,`incr_,`body_]:=(ini;While[cond,body;incr];)` は本体を実行するはずで、同じ式を別名 (For2) で定義すると gs も実行する。3 引数の `For` は未評価のまま | `Packages/functions.n:34` (読み込み時の名前解決と思われるが原因は未特定) | `For` の本体が動かない (`Do`・`While` は正常) | 低 | 実測 (gs で確認。`local/wake/` の f1〜f5 相当) | pysad は定義どおり本体を実行する。3 引数は gs と同じく未評価 |

## 1. 落ちる・動かない (重要度 高)

### B01 MAIN の TRACK が SIGSEGV で落ちる

- 現象: MAIN 部で `EMIT USE=RING;` に続けて `TRACK USE=RING NPART=5 TURNS=10;` を実行すると、フラグ一覧 (RFSW RADCOD RAD ... SMEAR) を表示した後に `Program received signal SIGSEGV`。`NPART=1 TURNS=-1` でも同じ。`NX=(1,2)` のリストは `*** $TRACK *** list is not allowed` で拒否される (`NX=1` と書けば通るが、その後同じ場所で落ちる)。MAIN の `EMIT` は `Use EMIT_TANCE or Emittance[] within FFS.` と言うだけで何も計算しない。
- 再現: `docs/history/work_main_track.md` §3 の最小スクリプト。
- 期待される動作: 古典 SAD と同様に `EMIT` が正規モード行列 r・閉軌道・EMITX/Y/Z・SIGE を作り、`TRACK` が追跡する。少なくとも未設定のとき SIGSEGV ではなくエラーメッセージを出す。`NX=(1,2)` の形式も古典 SAD では許される。
- 原因: バックトレースを `nm` (obj/AMD64-Linux5/sad1.exe) で引くと `toplvl_ → doact_ → acttra_ → track_ → tracka_ → __trackdlib_MOD_tinip` (+0x1db)。`tinip` の `tmap(xa,-1)` が `use temw, only: rn=>r` の行列を使うが、MAIN の EMIT が r を作らないので未設定。FFS 内で `EMIT` すれば r はできるが、FFS から MAIN に戻って `TRACK` を実行する経路がない (`END;` の後は SADScript の対話、`GetMAIN["file"]` はファイルを SADScript として読む)。
- 影響: MAIN の `TRACK` の参照出力が K64 では取れない。
- 重要度: 高 (落ちる)。
- pysad での扱い: ActTra.f → track.f → tracka.f を直訳して `pysad/maintrack.py` に実装 (EMIT は古典の意味)。参照は「同じ手順を SADScript で書き追跡だけ TrackParticles に任せた版」(`bench_level27a.sad`, `27b`) との照合で代用 (座標 ≤5e-11)。

### B02 ANGLE=0 かつ K1≠0 かつ K0≠0 の BEND で FFS 初期化が SIGSEGV

- 現象: `BEND` に `ANGLE=0` かつ `K1≠0` かつ `K0≠0` を同時に与えると、FFS の初期化 (`USE RING` 後の最初の CALC の前) で SIGSEGV (exit 139)。L の値、ROTATE の有無に依らない。K1 のみ (K0=0)、K0 のみ (K1=0) は動く。
- 再現: `bench_level25e` の初版 (tmp25/bis3.py)。`docs/history/work_rotate.md` §4、`docs/history/work_bend_k0.md` §7。`BEND` に `ANGLE=0 K0=… K1=…` を同時に与えた LINE を `USE` する (最小スクリプトは未保存)。
- 期待される動作: tbend.f:565-611 の分岐 (ANGLE=0 で K1=0 なら tsteer、K0=0 なら tquad、他は tmulti) に従って計算される。
- 原因: 未調査 (初期化時に落ちるので、tmulti 系へ分岐する要素の設定に不備があると推測)。
- 影響: この組み合わせの要素を含むラティスは読み込めない。参照値が取れない。
- 重要度: 高 (落ちる)。
- pysad での扱い: CALC/emit は tmulte (`_bend_as_mult`) で計算、追跡は NotImplementedError。SAD との照合は不可。

### B03 SetElement を含むスクリプトの出力をリダイレクトすると SIGSEGV

- 現象: `SetElement[QF,QUAD,K1->0.3]` などを含むスクリプトを、標準出力をファイルにリダイレクトして実行すると SIGSEGV。pty 経由 (`script -qec`) で実行すると同じスクリプトが動く。
- 再現: `bench_level24f.sad` (164 タグ)。`docs/history/work_item4.md` §F。
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
- pysad での扱い: 再現 (`radint.default_stale_rho(els)`、`compute(stale_rho=)` で上書き可)。

### B07 RADINT の六極ブロックが直前の DL を使う

- 現象: 六極のブロックが DL を設定せず、直前の BEND/QUAD の最後のスライス長 DL を使うので、I9, I10 (=γ·DL) に六極でも DL が足される。DKL2 は代入されず六極の CKL0X/CKL1X は 0。L=0 の六極は TCAL の後で `DL=0.2` (単精度の 0.2) にして ` L(SEXTUPOLE) EQ 0. SET L=0.2` を出す。DL, TETA, R0/R0I, FUNI などモジュール変数がプロセス中残るので、2 回目の RADINT の DL の初期値も前回の値。
- 再現: `bench_level29a.sad`、`29a2`。
- 期待される動作: 六極のスライス長 (L) を DL に設定する。L=0 は薄レンズとして扱い、定数 0.2 を足さない。
- 原因: INTGRL.f の六極の分岐 (SEXTU 呼び出し) が DL を代入しない。
- 影響: I9, I10 (色収差の GX/GY)。六極が多いリングで値が変わる。
- 重要度: 高。
- pysad での扱い: 再現 (docs/history/work_small_items.md A.2-2)。

### B08 tsolqur が最後の粒子の bzp を全粒子に使う

- 現象: 放射あり (RAD ON) でソレノイド内 (bz≠0) の四極は、最後の半ステップで `tz%tz0%bzp` (最後に tzsetparams を呼んだ粒子の値) を全粒子に使う。結果が粒子の並び順・粒子集合に依存する (差 ∝ δ_i − δ_last)。放射なしの tsolqu は粒子ループ内で完結するので起きない。
- 再現: Level 9 (ベンチ名は docs/history/00_RESUME.md の Level 9 の項) で同じ座標の粒子セットを変えると SAD 自身が 1.15e-6 変わる。
- 期待される動作: 各粒子が自分の bzp を使う。
- 原因: tsolqu.f tsolqur の tzsetparams (:268) と bzp の使用 (:273-)。粒子ループの外で tzsetparams を呼んだ後の値を使っている。
- 影響: 粒子の並び・数で追跡結果が変わる。DA 計算の粒子スロット順にも影響する。
- 重要度: 高。
- pysad での扱い: 再現 (track6 を要素外側・粒子内側に組み直し、`_tsolqu_bz_final` に最後の粒子の bzp を渡す)。SAD との比較は同じ粒子セットで行う。

### B09 tchge が回転行列を上書きし一周行列に回転が入らない

- 現象: tchge の入口で th≠0 かつ回した座標系のずれ≠0 (または dz, dtheta≠0) のとき、tbrote/tbgrote/tdrife1 が回転行列 `trans1` を「上書き」する。COD は回るが一周行列に回転が入らない。出口は tmultr5 で掛けるので回転が入る。一周行列が R^T·M になり、結合 (R1..R4)・チューンが物理とずれる。phig≠0 では tdrife1 (BEND 中点のドリフト) が行列に入る。四極・MULT の tsolrote にはない。
- 対象: BEND (ROTATE と DX/DY/DROTATE)、薄い多極 (SEXT/OCT/DECA と L=0 の四極)。
- 再現: `bench_level25f.sad` (25g は docs/11 A25 の記載のみ)。
- 期待される動作: 入口でも行列に回転が掛かり、一周行列が R_out·M·R_in になる。
- 原因: tchge.f の `transa=>trans1` (:45) と tbrote への渡し方 (:64-68, :83)。入口の呼び出しが trans を代入で上書きし、掛け合わせない (出口は :113-114 で tmultr5 により掛ける)。
- 影響: ロールした BEND を含むリングのチューン・結合が変わる。この癖のある要素を含むリングの Emittance[] は SAD の tcod の収束で決まる。
- 重要度: 高。
- pysad での扱い: 再現 (`lattice._tchge_jac`、`_thin_cod_matrix`、emit6d `_tbende`/`_tthine_gen`)。ただしこの要素を含むリングの Emittance[] (δ、チューン) は pysad と一致しない。

### B10 tdusrn が NaN の anp で 10000 回ループして乱数を消費する

- 現象: 粒子が落ちる周は発散中の値が NaN/inf のまま tradkfn が光子の乱数を引く。anp=NaN (g=NaN など) だと tdusrn が npmax=10000 回ループして 20001 個の乱数を引く (anp<0 なら 1 個)。ftest で tdusrn(NaN) と tdusrn(+Inf) は 10000 回、tdusrn(−Inf) は 0 回を確認。
- 再現: `bench_level16f`、`16l`、`20` (small)。`docs/history/work_loss_rng.md`、`history/work_loss_rng2.md`。
- 期待される動作: 非有限の anp では光子を引かない (あるいは損失粒子には放射を呼ばない)。
- 原因: drndsr.f :71-111 の tdusrn。`p .le. pn` と `pn .gt. cutoff` が NaN で偽になるので exit せず、npmax (drndsr.f:3、10000) まで回る。呼び出し側 (photons.f tradkfn) が損失判定前の粒子にも呼ぶ。
- 影響: 損失粒子が出た周以降の乱数列 (seed) が大きくずれ、他の粒子の放射も変わる。DA 計算 (DynamicApertureSurvey) の再現性にも関わる。
- 重要度: 高 (結果が変わる)。ただし NaN の扱いは B17 (gfortran -Ofast の MAX/MIN) と連動する未定義動作的な面がある。
- pysad での扱い: 再現 (`Prog.la_count`、NaN の特別扱い)。再同期テストで 30 周 seed 一致。発散中の分岐は丸めの 1 ulp で逆になるので、自由走行の最終 seed は SAD と一致しない (SAD 自身も 1e-16 の摂動で変わる)。

## 3. 特定条件で結果が変わる誤り (重要度 中)

### B11 4D の CALC が縦方向の Twiss 列を作り直さない

- 現象: `CALC4D;` の後の CALC で、縦方向の列 (AZ BZ NZ ZX ZPX ZY ZPY) が作り直されず、直前の 6D の CALC の先頭行の値が全行に残る。`CALC6D;`/`CALC4D;` 自体はフラグを変えるだけで再計算しない (これは仕様)。
- 再現: `bench_level29c.sad`。6D の CALC の後に 4D の CALC をして `Twiss["AZ",...]` などを全行で読む。
- 期待される動作: 4D では縦方向の列を 0 にする (または NaN など明示的な値にする)。
- 原因: tffscalc.f:87-90 は 4D では 15 列だけ入れ直し、qtwiss1 が縦の列を直前の行から写す。
- 影響: Twiss[] の縦の列を読むスクリプトが古い値を読む。
- 重要度: 中。
- pysad での扱い: 再現 (`FFS._z_stale`)。

### B12 RADINT の定数が単精度

- 現象: `COEFF=(elradi/(elmass**3)/3)*1e27` の `1e27` は単精度 (9.99999988e26) なので ELOSS が相対 1.15e-8 小さい。`FREQ*1.E-6`、`DL=0.2`、`EXP(2.3026*TAUQ)` も単精度。
- 再現: `bench_level29a.sad`。SAD の ELOSS (678.50075327596767) が 17 桁で一致して確認。
- 期待される動作: `1d27`、`1.d-6`、`0.2d0` と倍精度で書く。
- 原因: INTGRL.f :35 (1e27)、:172 (DL=0.2)、:212 (FREQ*1.E-6)、:266 (2.3026)。
- 影響: ELOSS、減衰時間、エネルギー幅など ELOSS から作る量が 8 桁目で変わる。
- 重要度: 中。
- pysad での扱い: 再現。

### B13 RADINT の CAVI が最後のもので上書き

- 現象: CAVI の VOLT/HARM/FREQ を合計せず最後に読んだ空洞の値で上書きする。HARM=0 なら FREQ·RLT/c、FREQ=0 なら HARM·c/RLT。
- 再現: `bench_level29a.sad`。空洞が複数あるラティスで RF 電圧・シンクロトロンチューン・バケット高さが総和にならない。
- 期待される動作: 電圧を総和し、周波数/ハーモニック数を整合させる。
- 原因: INTGRL.f の CAVI の処理 (代入であって `+=` でない)。
- 影響: 複数空洞の RADINT の RF 関連量。
- 重要度: 中。
- pysad での扱い: 再現 (parser が HARM から FREQ を補う分は SAD の定義値に戻して判定)。

### B14 tsolqum の z の dv 項が tsolqu と違う

- 現象: tsolqu は `−0.25·(… + ap·aln0 + …) − dv·aln0`、tsolqum は `−0.25·(… + (ap − dv)·aln0 + …)` で dv に 0.25 が掛かる。δ=1e-3 で 1 周 5e-10 の z のずれ。どちらが正しいかは物理から決める必要があるが、少なくとも両者は一致すべき。
- 再現: `bench_level11.sad` (観測 6.7e-10)。
- 期待される動作: 同じ量を計算する 2 つのルーチンで dv 項が同じ。
- 原因: tsolqu.f :86-88 と :250。
- 影響: 放射なし MULT 追跡の z。
- 重要度: 中。
- pysad での扱い: 再現 (`_tsolqu(qum=True)`)。

### B15 UnstableLevel の綴り誤りで常に 0

- 現象: マッチング終了後に `UnstableLevel` が更新されず常に 0。`MatchingResidual=10^n·nstab+r` には nstab が入る。
- 再現: `bench_level28*.sad` (不安定のまま終わるマッチング)。
- 期待される動作: nstab (不安定な点の数) が `UnstableLevel` に入る。
- 原因: tffscalc.f:58 `isl = ktfsymbolz('UnstableLevel',14)-4`。名前は 13 文字なのに 14 を渡しているので、別のスロットに書かれる (`rlist(isl)=nstab` が UnstableLevel に届かない)。`MatchingResidual` (16) と `NetResidual` (11) は文字数が合っている。
- 影響: スクリプトから不安定レベルが読めない。
- 重要度: 中。
- pysad での扱い: 再現 (`ffscmd.py:2141-2143`、`UnstableLevel=0` 固定。詳細は docs/history/work_unstable.md)。

### B16 CalculateOptics の stabz と idim=3 の stab/trace が未初期化

- 現象: `CalculateOptics[lb,le,inicond,cell,idim]` の戻りの `{stabx,staby,stabz}` の stabz が呼び出しごとに 0/1 が不定 (部分区間 cell=0 では 0 になった)。idim=3 の stab/trace も未初期化 (`{0,0,0} {1,0,0}` が返った)。
- 再現: `bench_level24d.sad`、`24d2` (docs/history/work_item4.md §D)。
- 期待される動作: 6 次元でないときは stabz を 1 か 0 に決めて返す。
- 原因: tfoptics.f → qcell1。`ffs_stat optstat` が intent(out) の局所変数で、縦方向 (z) の安定判定が代入されない経路がある。
- 影響: 戻り値を使う判定が不定。
- 重要度: 中。
- pysad での扱い: 1 を返す。ベンチでは stab は先頭 2 要素だけ、idim=3 は stab/trace を除いて比較。

### B17 gfortran の MAX/MIN が NaN を無視する前提のコード

- 現象: `-Ofast` の gfortran では `max`/`min` の NaN の返り値が箇所ごとに違う (tradkfn の `max(gmin, g+dg)` は NaN で gmin、asinz の `max(-1,x)` は −1、sqrt1 の `max(1e-200,·)` は後者)。tffsmatch の fuzz は NaN でモードを変えない (残差が悪化すると crate<0 → log10 が NaN → min/max が NaN を無視して badcnv=0 → DESCEND に切り替えず fmincube の直線探索に行く)。比較演算も NaN で真になるものがある (`v==v` が T、`v/=0` が F)。
- 再現: `bench_level16f/16l/20`。
- 期待される動作: NaN で振る舞いが言語実装に依存しない (NaN を明示的に判定する)。
- 原因: photons.f tradkfn、tbend.f tbendcore の asinz、tffsmatch.f:255-289 の fuzz。コンパイラの未定義動作への依存。
- 影響: 損失粒子の乱数消費 (B10)、マッチングの途中経過 (再現しないと早く "Residual =" で止まる)。
- 重要度: 中 (バグというより未定義動作依存)。
- pysad での扱い: 再現 (match6 は `_log10` が NaN を返し fuzz が 1。track6n は `sqrt1`/`pxy2dpz`/`asinz`/`gfloor`、NAN_MODE)。

### B18 PolynomialFit が平均 0 のとき xl=1

- 現象: 平均 0 (xa=0) の x に対して `PolynomialFit` が xl=1 にして係数の次数補正が壊れる。クロマティシティ補正の DP0 スキャンで 1e10 級の係数になる。
- 再現: DP0 を ±対称に振って tune を `PolynomialFit` に渡す。
- 期待される動作: xa=0 でも次数補正が正しい (`xl` は xa の大きさのスケールで、0 にならないように別に決める)。
- 原因: Packages/Fit.n :37-56 `If[xa, Table[...], 1]`。xa=0 だと偽 (0) として else 側の 1 が使われる。
- 影響: 係数が誤る。
- 重要度: 中。
- pysad での扱い: sadscript.py が Fit.n をそのまま評価するので同じ挙動 (癖も再現)。

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
- pysad での扱い: tcod の直訳ではなく close_z (割線法) の近似。

### B22 MON/MCOD/TCOD が何もせず "Missing value for PM*"

- 現象: このビルドでは `MON`/`MCOD`/`TCOD` コマンドが何もせず式として評価される。`MON 0 0 0 0 PM*` は `Missing value for PM*` で停止する。
- 再現: `MON 0 0 0 0 PM*` を含む FFS 文。
- 期待される動作: モニタを定義する (古典 SAD の MON コマンド)。
- 原因: K64 の tffsa.f に MON・MCOD・TCOD のコマンドは無い (`grep` で 0 件。2026-11-16 にソースで確認)。FFS では式として評価されて何もしない。
- 影響: モニタを使う軌道補正スクリプトが動かない。
- 重要度: 低〜中。
- pysad での扱い: 無視 (ffscmd.py:1060)。

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
- pysad での扱い: SAD と同じ narg エラーを再現。

### B25 TrackParticles の POL ON の列順がソース内で自己矛盾

- 現象: tftrack.f が列順を 2 通りに書いている。
- 再現: ソースの読みのみ。
- 原因: エラーメッセージ (:344) は `{z},{dp}[,{sz},{psis}],{flg}`、コメント (:502) は `(x, x', y, y', z, delta, phis, sz, alive)`。実コード (:455-460) は kaj(7)=atan(...)=psis、kaj(8)=zx(8)/|zx(7:9)|=sz の順で、コメントと一致する。メッセージが誤り (または逆)。
- 影響: 列順の文書不整合 (POL ON を使う場合)。
- 重要度: 低。
- pysad での扱い: POL は未対応 (docs/history/review_notes.md に保留)。

### B26 ORG の引数が式として読まれる

- 現象: `ORG B3 1.5 -2.5 0.5 10 20 30` が `1.5 -2.5` を式 1.5−2.5=−1.0 と読み、gx=−1.0 として 1 引数分欠落する。欠落分は 0 になり直前の ORG の値は残らない。
- 再現: `bench_level19c.sad` (6 通り)。
- 期待される動作: スペース区切りの引数は別の値として読まれる (または仕様として明記される)。
- 原因: tffsa.f の ORG が getvad (式として読む) で引数を取る。FFS の他のコマンドと同じ仕様なので、バグというより仕様の落とし穴。
- 影響: 負数を含む ORG の値が意図と変わる。`0-2.5` のように続けて書くか別行にする。
- 重要度: 低 (仕様の可能性が高い)。
- pysad での扱い: 再現。

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
- pysad での扱い: 再現 (MAIN の RFSW は `main_flags` に入れ、FFS の RFSW は ON で始める)。

### B31-B33 wake まわりの SIGSEGV (原因と修正の差分: `docs/history/sad_patch/wake_sigsegv_fixes.patch`)

- B31: `tffswake` は `iwl=abs(kwaketbl(1,iwp))`、`iwt=abs(kwaketbl(2,iwp))` の後に、表の有無を見ずに `nlw=(ilist(1,iwl-1)-1)/2`、`ntw=(ilist(1,iwt-1)-1)/2` を計算し `c_f_pointer(c_loc(rlist(iwl)),…)` する。表がない側は 0 なので `ilist(1,-1)` を読む。後段の `tffswakekick` は `iwl /= 0` で守られている。
- B32: `tturn` は始点の SOL の前半キックを `tsol` の前、後半キックを `tsol` の後で呼ぶ。`tsol` は領域全体を進めて戻る。領域の各要素で `if(sol)` 分岐が `go to 1020` するとき `aewak` が古い真のままで、1020 の後半キックが領域の要素ごとに呼ばれ、`nwak` が表の数を超えて `twxiwp` が表の外を読む。始点のウェイクが `nextwake` を始点の位置に残すので、`tsol` の中のウェイクは (元でも) 効かない。
  修正は、`tsol` の前に始点の後半キックと `nwak` の更新を済ませ、`aewak` を偽にすること。直した SAD は pysad の追跡と、始点だけのウェイクでも領域内との併用でも 10 桁以上一致した。
- B33: `tapert` の「落ちた粒子を報告する」条件が `idtypec(l_track) == icAprt` を Fortran の非短絡の `.or.` で評価する。最初の `tapert` は `tturn` の最初の要素の前に呼ばれ、`l_track` はまだ 0。最初から LOSSAMPL の外にいる粒子があると落ちる (wake とも SOL とも無関係。w3 の LOSSAMPL=0.0009 で再現)。
  直した SAD の出力座標が pysad と合わなかった件は 2026-11-16 に原因が分かり解消した: `tapert` は落とす印に `kptbl(i,4)=l_track` を書くので、セッション最初の TrackParticles (l_track=0) の最初の `tapert` では印が 0 (生きている) になり、最初から LOSSAMPL の外の粒子は落ちず (座標は 3·LOSSAMPL に頭打ちされて追跡を続け)、後の `tapert` で落ちる。2 回目の TrackParticles からは最初の要素で落ちる。pysad は `TrackRun.l_track0` でこれを写し、gs (1・2 回目) と直した SAD に全桁一致 (test_small_items.py)。

### B63-B70 2026-11-16 までに見つかった癖・疑い

どれも `docs/history/work_unimpl.md` と docs/11_sad_quirks.md (A42〜A48) に経緯がある。「確定」は gs の実測、「疑い」は原因・意図が未確認のもの。

- B63: 確かめた条件は、PLANAR のリングに SOL (BOUND≠0、L=0.3・0.6、BZ=0) を足した周長の比較 (偏極時間が 15 桁一致。docs/29・docs/44)。ほかの磁場・エネルギーでの確認はしていない。GEO のない SOL の S 座標の途切れは、pysad が GEO を見ないので再現しない (docs/29)。
- B64: 疑いのまま。pysad は K64 と同じに写したが、gs と差の出る例が無く (試験のリングでは効かない)、この `merge` が意図かどうかと、結果に効くかは未確認。
- B65: RADLIGHT の「記録しない」のは tlstore の仕様 (型 < CAVI の要素の後だけを記録)。意図なら癖、そうでなければ誤り。
- B66: init.n の AutoLoad の行が指す `WAKECOMMAND.n` が配布に無いこと自体は確認済み (ファイルが無い)。
- B67: 仕様と考えられるが、1 GeV と思い込んだ文書が過去にあったため載せた。
- B68: 現在の gs が RAD ON・BZ≠0 の領域内の四極の力を入れていない疑いが強い (証拠は local/wake/verify_log.md の bench_level9 の項。単一粒子 1 周追跡で、領域内がドリフトだけなら gs=pysad、四極ありは BZ≠0 で gs が四極なしの結果と一致)。ソースのどこかは未特定、保存した参照 (10-14 の gs) を出した版は今は動かず確認できない。
- B69: どちらが正しいか決められない。pysad 側の誤りの可能性も残る。pysad と gs が合わない未解決の差のうち、B68・B69 以外 (BEND の F1=0.1 FRINGE=3 の RAD ON 追跡で相対 0.7、L≠0 の CAVI の DX・DY で 1e-3、K3 が強い MULT リングの CALC (qcod) の閉軌道) は、SAD の誤りか pysad の誤りかも未判定なのでこの表に載せていない (docs/history/work_unimpl.md)。

## 5. 開発者への依頼の要点

優先して直してほしい順に:

1. B01, B02, B03 (落ちる。再現が最小で、原因がソースから追える)。
2. B06, B07, B12, B13 (RADINT。INTGRL.f の中で閉じていて、修正が小さい)。
3. B08, B09, B10, B14 (追跡・行列の結果が変わる)。
4. B15, B16, B23 (綴り・未初期化・引数の並び。修正は 1 行)。

なお、pysad は「SAD の数値を再現する」ことが目的なので、上記のうち結果が変わるものは修正されても pysad 側の再現を外す必要がある (対応箇所は各項の「pysad での扱い」)。

## 6. 本レポートに含めなかったもの (参考)

- pysad 側の誤りだったもの: fmincube の `m_euler` (SAD はオイラー定数 γ、pysad がネイピア数 e を入れていた)、MASS の既定値、`MINCOUP` の既定値、`Cavity.PHI` の既定値、`lattice._tbendebody0_mat` の z 行の符号。
- 仕様・癖: docs/11_sad_quirks.md の A1-A48 のうち上記 (B63〜B71) 以外、B1-B5、C1-C23、D1-D5。
