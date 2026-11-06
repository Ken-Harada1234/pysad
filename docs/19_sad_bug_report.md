# 19. SAD K64 のバグ報告 (pysad 移植中に見つかったもの)

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
| B04 | `draw` (WriteGraphs) で core dump | 未調査 (Tk/描画) | 落ちる | 高 | 確定 (実行) | 描画は対象外 |
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
| B22 | MON/MCOD/TCOD が何もせず、`MON 0 0 0 0 PM*` が "Missing value for PM*" | FFS コマンド表 | 補正スクリプトが止まる | 低〜中 | 実測 | 無視 |
| B23 | Extended Twiss の EPX の欄に EY が出る | temit.f :1010-1017 | 表示のみ | 低 | 確定 | 再現 |
| B24 | `LINE["LENGTH","$$$"]` が narg エラー | LINE[] の引数処理 | 位置が取れない | 低 | 実測 | 再現 |
| B25 | TrackParticles の POL ON 列順がソース内で自己矛盾 | tftrack.f :344 と :502 | 文書・メッセージの不整合 | 低 | ソース | 未対応 (POL) |
| B26 | `ORG` の引数が式として読まれ `ORG B3 1.5 -2.5 ...` が `1.5-2.5` になる | tffsa.f ORG、getvad | 負数の書き方で値が変わる | 低 | 確定 | 再現 |
| B27 | SeedRandom[] の表示と seed の末尾桁の丸め | 表示ルーチン | 表示のみ | 低 | 実測 (調査浅い) | 符号なし表示に修正 |
| B28 | 未初期化変数の値が出力に出る (EMIT 報告の Closed orbit Entrance、MAIN TRACK の l_track) | temit.f、tracka.f | 表示のゴミ | 低 | 実測 | 再現せず |
| B29 | QUAD に RAD=1 を書くと COD・エミッタンスが 0 になる | 原因未調査 | 結果が変わる (要調査) | 要調査 | 実測 | 未対応 |
| B30 | MAIN の `OFF RFSW` (`$RFSW$`) が FFS の RFSW に届かない (他の RAD/FLUC/COD/RADCOD/INTRA は届く) | track.f :107, :149 (`rfsw=.true.` で tffs を呼ぶ) | `OFF RFSW` を書いても空洞が効く。FFS 入口の表の RFSW 列は MAIN の値で FFS の値と違う | 中 | 確定 (実測 + ソース) | 再現 (`main_flags`) |

## 1. 落ちる・動かない (重要度 高)

### B01 MAIN の TRACK が SIGSEGV で落ちる

- 現象: MAIN 部で `EMIT USE=RING;` に続けて `TRACK USE=RING NPART=5 TURNS=10;` を実行すると、フラグ一覧 (RFSW RADCOD RAD ... SMEAR) を表示した後に `Program received signal SIGSEGV`。`NPART=1 TURNS=-1` でも同じ。`NX=(1,2)` のリストは `*** $TRACK *** list is not allowed` で拒否される (`NX=1` と書けば通るが、その後同じ場所で落ちる)。MAIN の `EMIT` は `Use EMIT_TANCE or Emittance[] within FFS.` と言うだけで何も計算しない。
- 再現: `docs/history/work_main_track.md` §3 の最小スクリプト。実スクリプト `VSX_Ring/memory17/v18_67` + `phase` でも同じ。
- 期待される動作: 古典 SAD と同様に `EMIT` が正規モード行列 r・閉軌道・EMITX/Y/Z・SIGE を作り、`TRACK` が追跡する。少なくとも未設定のとき SIGSEGV ではなくエラーメッセージを出す。`NX=(1,2)` の形式も古典 SAD では許される。
- 原因: バックトレースを `nm` (obj/AMD64-Linux5/sad1.exe) で引くと `toplvl_ → doact_ → acttra_ → track_ → tracka_ → __trackdlib_MOD_tinip` (+0x1db)。`tinip` の `tmap(xa,-1)` が `use temw, only: rn=>r` の行列を使うが、MAIN の EMIT が r を作らないので未設定。FFS 内で `EMIT` すれば r はできるが、FFS から MAIN に戻って `TRACK` を実行する経路がない (`END;` の後は SADScript の対話、`GetMAIN["file"]` はファイルを SADScript として読む)。
- 影響: MAIN の `TRACK` (ユーザーの実スクリプトの 725 本が使う書き方) の参照出力が K64 では取れない。
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

- 現象: `WriteGraphs` 内の `draw` が WSL の SAD で core dump (`dumped core`)。plot_inj は Tk 初期化失敗の後に SIGSEGV。
- 再現: esrf の `disp_b`、`disp_arc_optics`、pf2 の `plot_inj_r` (docs/history/work_recheck.md、history/work_pf2_survey.md)。
- 期待される動作: Tk/描画が使えない環境では描画をスキップするかエラーで続行する。
- 原因: 未調査 (WSL に表示環境がないための可能性が高い。環境依存)。
- 影響: 描画を含むスクリプトが最後まで走らない。
- 重要度: 高 (落ちる。ただし環境依存)。
- pysad での扱い: 描画は対象外。参照の取得時は `FFS["draw ..."]` を `Print["nodraw ..."]` に置換 (`<work-recheck>/nodraw/`)。

### B05 FITP>1 の GO で再起動を繰り返す

- 現象: `FITP` が 1 より大きい条件 (オフモーメンタムのマッチング、クロマティシティ補正) で `GO` すると、`*** Run time Environment` からスクリプトを最初から再実行する。`NPARA=1` でも起きる。
- 再現: `bench_level15a`、esrf の chrom。再起動の回数は `chrom1.sad.out` に 61 回、`chrom_b.sad.out` に 3 回、`<work-pf2>/chrom_r.sad.out` に 47 回。
- 期待される動作: `GO` が最後まで走って結果を返す。
- 原因: 未調査 (プロセスが異常終了して親が再起動していると見られる。実行環境の再起動機構の影響の可能性もある)。
- 影響: 最終解が取れず、比較は「再起動の途中で出た値」との照合になる。
- 重要度: 高 (落ちる)。
- pysad での扱い: match6/ffscmd の FITP>1 は直訳のまま未照合。`FITP 1` に落とした版で解の一致を確認する。

## 2. 結果が変わる誤り (重要度 高)

### B06 RADINT が全 BEND に最後に処理した BEND の L/ANGLE を使う

- 現象: RADINT の BEND の rhob/rho0 が、その BEND 自身の L/ANGLE ではなく、直前に tbend 系が設定した値 (ビーム順で最後に処理した曲がりのある BEND の L/ANGLE) になる。ESRF では最後が反転した ARCAH2 の BID5 (L/ANGLE=20.757 m) が全 BEND に使われる。
- 再現: `bench_level29a.sad`、`29a2`。別の値 (BPI7 の 37.6、BID1 の 77.5) を仮定すると I9 が 8e-4〜1e-3 ずれる。20.757 m で SAD と 7 桁一致するので裏づけられた。
- 期待される動作: 各 BEND の L/ANGLE を使う。
- 原因: INTGRL.f が呼ぶ tbendebody は bendib モジュールの rhob/rho0 を使うが、INTGRL は設定しない (tbendiinit/tbendeinit を呼ばない)。
- 影響: I1X, I2, I3, I4, I5, I9, I10 など。RADINT の値が BEND ごとの正しい値からずれる (docs/16 §11 の RADINT と Emittance[] の 13% 差の候補の 1 つ。因果は未検証)。
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
- 再現: `bench_level24d.sad`、`24d2`、esrf の `CalculateOptics` を使う補正スクリプト (docs/history/work_item4.md §D)。
- 期待される動作: 6 次元でないときは stabz を 1 か 0 に決めて返す。
- 原因: tfoptics.f → qcell1。`ffs_stat optstat` が intent(out) の局所変数で、縦方向 (z) の安定判定が代入されない経路がある。
- 影響: 戻り値を使う判定が不定。
- 重要度: 中。
- pysad での扱い: 1 を返す。ベンチでは stab は先頭 2 要素だけ、idim=3 は stab/trace を除いて比較。

### B17 gfortran の MAX/MIN が NaN を無視する前提のコード

- 現象: `-Ofast` の gfortran では `max`/`min` の NaN の返り値が箇所ごとに違う (tradkfn の `max(gmin, g+dg)` は NaN で gmin、asinz の `max(-1,x)` は −1、sqrt1 の `max(1e-200,·)` は後者)。tffsmatch の fuzz は NaN でモードを変えない (残差が悪化すると crate<0 → log10 が NaN → min/max が NaN を無視して badcnv=0 → DESCEND に切り替えず fmincube の直線探索に行く)。比較演算も NaN で真になるものがある (`v==v` が T、`v/=0` が F)。
- 再現: `bench_level16f/16l/20`、esrf/matching。
- 期待される動作: NaN で振る舞いが言語実装に依存しない (NaN を明示的に判定する)。
- 原因: photons.f tradkfn、tbend.f tbendcore の asinz、tffsmatch.f:255-289 の fuzz。コンパイラの未定義動作への依存。
- 影響: 損失粒子の乱数消費 (B10)、マッチングの途中経過 (再現しないと早く "Residual =" で止まる)。
- 重要度: 中 (バグというより未定義動作依存)。
- pysad での扱い: 再現 (match6 は `_log10` が NaN を返し fuzz が 1。track6n は `sqrt1`/`pxy2dpz`/`asinz`/`gfloor`、NAN_MODE)。

### B18 PolynomialFit が平均 0 のとき xl=1

- 現象: 平均 0 (xa=0) の x に対して `PolynomialFit` が xl=1 にして係数の次数補正が壊れる。クロマティシティ補正の DP0 スキャンで 1e10 級の係数になる。
- 再現: DP0 を ±対称に振って tune を `PolynomialFit` に渡す (esrf の calc_chromaticity)。
- 期待される動作: xa=0 でも次数補正が正しい (`xl` は xa の大きさのスケールで、0 にならないように別に決める)。
- 原因: Packages/Fit.n :37-56 `If[xa, Table[...], 1]`。xa=0 だと偽 (0) として else 側の 1 が使われる。
- 影響: 係数が誤る。ユーザーは自作の MyPolyFit で回避している。
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
- 影響: リスト処理の結果。実スクリプトでの使用は確認していない。
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

- 現象: このビルドでは `MON`/`MCOD`/`TCOD` コマンドが何もせず式として評価される。`MON 0 0 0 0 PM*` は `Missing value for PM*` で停止する (esrf の `correct_cod` v05 版)。nomon 版は完走する。
- 再現: esrf の `correct_cod` (docs/history/work_recheck.md 項目 1)。
- 期待される動作: モニタを定義する (古典 SAD の MON コマンド)。
- 原因: K64 の FFS コマンド表に MON が登録されていない可能性 (未確認)。
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

### B24 LINE["LENGTH","$$$"] が narg エラー

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
- 原因: 未調査。バグかどうかも未確認 (RAD が QUAD で別の意味の可能性)。
- 影響: 切り分けには DampingRate を使う。
- 重要度: 要調査。
- pysad での扱い: 未対応。

### B30 MAIN の `OFF RFSW` が FFS の RFSW に届かない

- 現象: MAIN 部 (FFS; の前) に `OFF RAD FLUC COD RFSW RADCOD;` と書くと、FFS 入口の表は `RFSW … F` と出るが、FFS 内の `FLAG[{"RFSW"}]` は 1 で、TrackParticles は空洞ありの結果 (z≠0、dp が位相で変わる) になる。RAD/FLUC/COD/RADCOD/INTRA/POL は MAIN の ON/OFF が届く。
- 再現: `<repo>/<work-item>/mflag1〜6.sad` (MAIN の ON/OFF の 6 通り → FFS で FLAG[] を表示)、`zdp_chk2_full.sad` (FFS 内の `NORFSW;`/`RFSW;` で TrackParticles の z/dp が変わる)。ユーザーの esrf `disp_R56`、`matching_test`。
- 原因: track.f:80-107 で `isynch=IgetGL('$RFSW$')`、`rfsw=isynch/=0` とした後、`nturn<0 .and. np0==0` (FFS 起動) の枝で `rfsw=.true.` に固定して `tffs` を呼ぶ (track.f:149)。表 (track.f:138) は固定前の値を出す。
- 影響: RFSW を切りたいスクリプトは FFS 内で `NORFSW;` と書く必要がある。`disp_R56` は結果的に空洞ありの R56/R66 を測っている。
- 重要度: 中 (仕様かもしれないが、表と実体が食い違う)。
- pysad での扱い: 再現 (MAIN の RFSW は `main_flags` に入れ、FFS の RFSW は ON で始める)。

## 5. 開発者への依頼の要点

優先して直してほしい順に:

1. B01, B02, B03 (落ちる。再現が最小で、原因がソースから追える)。
2. B06, B07, B12, B13 (RADINT。INTGRL.f の中で閉じていて、修正が小さい)。
3. B08, B09, B10, B14 (追跡・行列の結果が変わる)。
4. B15, B16, B23 (綴り・未初期化・引数の並び。修正は 1 行)。

なお、pysad は「SAD の数値を再現する」ことが目的なので、上記のうち結果が変わるものは修正されても pysad 側の再現を外す必要がある (対応箇所は各項の「pysad での扱い」)。

## 6. 本レポートに含めなかったもの (参考)

- pysad 側の誤りだったもの: fmincube の `m_euler` (SAD はオイラー定数 γ、pysad がネイピア数 e を入れていた)、MASS の既定値、`MINCOUP` の既定値、`Cavity.PHI` の既定値、`lattice._tbendebody0_mat` の z 行の符号。
- 仕様・癖: docs/11_sad_quirks.md の A1-A33 のうち上記以外、B1-B5、C1-C23、D1-D5。
