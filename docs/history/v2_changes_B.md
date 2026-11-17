# v2.0.0 の変更記録: 担当 B (特殊関数・数値、SADScript、FFS の小項目)

作成: 2026-11-17。この文書は v2.0.0 の作業のうち担当 B の分の記録です (まとめは docs/47 に入る予定)。

- 互換モード: `pysad.compat` (`PYSAD_COMPAT=v2|k64`)。v2 (既定) が修正した動作、k64 が K64 (gs) の動作。
  下の表の「k64 モード」の欄が「同じ」のものは、モードで切り替えていない (pysad 自身の不具合の修正・未対応だった機能の追加)。
- 根拠: 数値は gs ではなく、独立な計算 (閉じた式、級数、AGM、積分表示の数値積分、Wronskian などの恒等式、numpy の最小二乗) と比べた。
  mpmath はこの環境に無いので使っていない (scipy.special と、scipy に依らない独立な計算の 2 つを突き合わせた)。
  gs (SAD 1.4.4.2k64、2026-10-19 ビルド、WSL) は「K64 の動作の確認」と「K64 が正しい仕様 (NOKEEPEXP、INDEX、文脈) の参照」に使った。
- 試験: `pysad/tests/test_v2_B_specfn.py` (13 件)、`test_v2_B_script.py` (25 件、うち 1 件は POSIX だけ)、`test_v2_B_ffs.py` (14 件)。計 52 件。
- 走る例: `example/v2_B/` (`B1_special_functions.sad`、`B2_package_context.sad`、`B3_shared_system.sad`、`B4_ffs_v2.sad` と pysad の出力)。
- 使った条件: 特殊関数は実数の倍精度。FFS の項目は小さな試験用リング (1 GeV・周長 6〜7 m の四極だけのリング、3 GeV・周長 105.6 m の FODO 24 セル + 空洞 1.5 MV・500 MHz)。
  実機の大きなリングでは確かめていない。

## 1. 一覧

| ID | K64 の動作 | v2 の動作 | 根拠 (式・数値) | 試験 | k64 モード |
|---|---|---|---|---|---|
| B43 | `Hypergeometric2F1[.5,.5,1,.99]` = 2.348414、`BesselJ[.5,30.]` = -.143924709693 (級数の打ち切り) | 2.352715816780、-.143929653370 (scipy) | 2F1(½,½;1;m) = 2K(m)/π を AGM で計算、J_{1/2}(x) = sqrt(2/πx) sin x。相対 1e-13 で一致 | `test_b43_hyp2f1_besselj_precision` | 同じ (pysad は以前から真値) |
| B44・B54 (Ei) | `ExpIntegralEi[1.5]` = .164-2.985 I、`Ei[-.1]` = -2.2265 | 3.301285449130、-1.822923958419 (実数) | 級数 γ + ln|x| + Σ x^k/(k k!)、負の側は -E1(x) = -∫ e^-t/t dt、大きい x は漸近級数 e^x/x Σ k!/x^k。1e-10〜2e-12 | `test_b44_b54_expintegralei` | 同じ |
| B45 | `HermiteH[3,.7]`・`LegendreP[3,1,.4]`・`LegendreP[2,2,.3]`・`Binomial[-3,2]` が NaN | -5.656、+.274954541697、2.73、6 | H3 = 8x³-12x、P_3^1 = -(1-z²)^½ (15z²-3)/2、P_2^2 = 3(1-z²)、P_4^3 = -105 z (1-z²)^{3/2}、二項係数の積 | `test_b45_hermite_legendre_binomial`、`test_b45_legendre_integer_m_is_limit_of_formula` | 同じ。**pysad の誤りを 1 件直した** (下の §2.1) |
| B46 | `HypergeometricU[1,2,.5]` = INF、`BesselY[0,.1]` = NaN、`BesselK[1.5,.7]` = INF、`BesselJ[-1,1.5]`・`BesselY[1,2.5]`・`BesselJZero[0,1]` は応答なし | 2、-1.534238651350、1.806573612779、-.557936507910、.145918137967、2.404825557696 | U(1,2,z) = 1/z、Y_0 の級数、K_{3/2} の閉じた式、J_n の積分表示、Wronskian J_1 Y_0 - J_0 Y_1 = 2/(πx)、J_0(零点) = 0 | `test_b46_integer_order_u_bessely_besselk` | 同じ |
| B47 | `InverseErf[.3]` = .272462711811 (Newton の係数が sqrt(π/2)、早い打ち切り。相対 1e-8 の誤差) | .2724627147267543 (scipy の erfinv)。`InverseErfc` は小さい z で桁落ちしない erfcinv | erf(InverseErf(x)) = x が 4e-16 以内、正しい係数 (sqrt(π)/2 e^{y²}) の Newton 法と 1e-13 で一致 | `test_b47_inverseerf_v2_true_value`、`test_fitgaussian_uses_true_inverseerf_in_v2` | K64 の値 (`test_b47_inverseerf_k64_unchanged`、既存の `test_k64_inverseerf_is_reproduced`) |
| B53 | 大きい引数の `BesselJ/K/Y` が不正確 (`BesselJ[0,33.]` 相対 3e-3)、`FresnelS[6.1]` = 1.3e7、`AiryAi[15.5]` = -16 | 真値 (scipy)。x = 33・100・1000 の Bessel、x = 6.1・20 の Fresnel、Airy | J_n・I_0 の積分表示 (台形則)、K_ν = ∫ exp(-x cosh t) cosh νt dt、Fresnel は Simpson 積分、Ai = (1/π) sqrt(x/3) K_{1/3}(2/3 x^{3/2})、Ai Bi' - Ai' Bi = 1/π | `test_b53_large_argument_bessel`、`test_b53_large_argument_fresnel_airy` | 同じ |
| B54 (PolyLog) | 複素数 (虚部 1e-15 の雑音つき) | 実数 (z ≤ 1)。z > 1 は切断線の下側 (Li_2(2) = π²/4 - iπ ln 2)。**非整数次は |z| の上限なし** (§2.2) | Li_2(½) = π²/12 - ln²2/2、Li_3(½) の閉じた式、-∫ ln(1-t)/t dt、Fermi–Dirac 積分、Im = -π ln^{s-1} z/Γ(s)、Li_{s-1} = z dLi_s/dz | `test_b54_polylog_real`、`test_b54_polylog_large_argument_noninteger_order` | 同じ |
| (追加) 複素引数 | — | Bessel J/Y/I/K (実数の次数)、Hypergeometric0F1・1F1・2F1 (実数のパラメータ)、Airy が複素引数で使える (以前は TypeError) | 積分表示・級数・Wronskian (複素でも成り立つ形) と 1e-12〜1e-13 | `test_complex_argument_bessel_hypergeometric` | 同じ |
| B18 | `PolynomialFit` は x の平均が 0 のとき (`xl=1` のまま二項係数の和を取り) 係数が壊れる: `{{-1,1},{0,0},{1,1},{2,4.1},{-2,4}}` の 2 次で {1.0257, 2.0486, 1.0143} | {-.0085714, .02, 1.0142857} (平行移動が無いので係数は sl そのもの) | numpy.polyfit と 1e-13 で一致。3 次の多項式を対称な点に当てはめて係数を回復 | `test_b18_polynomialfit_mean_zero_v2`、`test_chromaticity_scan_in_user_package` | K64 の値 (`test_b18_polynomialfit_k64_and_nonzero_mean`)。平均が 0 でないデータは両モード同じ |
| B48 | `FixedPointList` は未評価 (autoload の綴り違い) | Misc.n の定義どおり動く | Newton の反復 x → (x + 2/x)/2 を Python で再計算 | `test_b48_fixedpointlist` | 同じ (pysad は以前から動く) |
| B52 | `Fit[data,fun,{x1,x2},…]` が narg を出して最小にならない値を返す | 最小二乗解 | numpy.linalg.lstsq の解・残差と一致 (係数 1e-6、ChiSquare 1e-8) | `test_b52_fit_vector_independent_variable` | 同じ |
| B61 | 複素の初期値の `FindRoot` が壊れた値 (`{z->Restrict[(.25+I),…]}`) | 複素の Newton 法で解く (§2.3) | z² = -1 → i、e^z = 2 → ln 2、z³ = 1 → (-1+i√3)/2、連立 zw = i・z + w = 2 を代入して残差 1e-12 以下 | `test_b61_findroot_complex_v2` | 従来どおり TypeError (`test_b61_findroot_complex_k64_is_unchanged`) |
| B70 | `MessageName[f,"usage"]` (関数形) は f を無視して `Undefined MessageName General::"usage"` | f の本文を返す (`f::usage` と同じ)。無ければ `General::tag` の本文 | 定義した文字列がそのまま返ること。`General::mytag` の継承は gs の `h::mytag` の実測と同じ | `test_b70_messagename` | 同じ (関数形は pysad が以前から f を見る。`General::tag` への継承を追加: gs の `::` 形と同じ) |
| B71 | `For[ini,cond,incr,body]` は本体を実行せず incr だけ繰り返す。3 引数は未評価 | 本体を実行 (Break・Continue も効く)。3 引数は本体なしのループ | t = Σ i² = 30、Break・Continue の動作 | `test_b71_for_body_v2` | **gs と同じにした** (本体を実行しない。3 引数は未評価): `test_b71_for_body_k64`。以前の pysad は本体を実行していた |
| B22 | `MON`・`MCOD`・`TCOD` はコマンドでなく、後ろの語が式として読まれて `Missing value for PM*` で文が止まる | 「K64 に無いコマンド」と警告して文全体を読み飛ばす (定義はしない) | K64 の tffsa.f にこれらのコマンドは無い (docs/21 B22)。モニタの定義は実装していない (§3) | `test_b22_mon_mcod_tcod` | K64 と同じ |
| B24 | `LINE["LENGTH","$$$"]` が narg エラー | エラーのまま。メッセージに位置の取り方 (`LINE["POSITION",loc]`、`LINE["LENG",loc]`・`LINE["S",loc]`) を足した | K64 のソース (tftwiss.f:508) が `LENGTH` を引数 1 個に限っている (要素数を返す)。バグではなく仕様と判断 (§3) | `test_b24_line_length_with_location` | K64 と同じ |
| B26 | `ORG B3 1.5 -2.5 0.5 10 20 30` が `1.5-2.5` と読まれ値が 1 つ足りなくなる | 空白で区切った `-2.5` は別の値 (gx=1.5, gy=-2.5)。`3/2`、`2 *5` などの式は書ける | K64 で曖昧さのない書き方 (`1.5 0-2.5 …`) と同じ ORG の値・同じ GX/GY/GZ | `test_b26_org_negative_numbers` | K64 と同じ |
| B27 | `SeedRandom[]` の seed が符号つきで表示される | 符号なし (seed mod 2³²) | -5 → 4294967291 | `test_b27_seedrandom_display_unsigned` | 同じ (以前から) |
| B30 | MAIN の `OFF RFSW` が FFS に届かない (FFS は RFSW ON で始まる) | 届く (`ON RFSW` も)。FFS の中の `RFSW;`・`NORFSW;` が後から優先 | FFS の `NORFSW;` と全桁同じ追跡結果 (δ が保存)、`RFSW;` と同じ結果 | `test_b30_main_off_rfsw_reaches_ffs` | K64 と同じ |
| B50 | `FileByteCount` が行数 (`wc -l`) | バイト数 | 14 バイト・2 行のファイルで 14 | `test_b50_filebytecount_is_bytes` | 同じ (以前から) |
| B67 | MOMENTUM の既定が 30 GeV で、書き忘れても何も出ない | 既定は同じ 30 GeV。MOMENTUM を MAIN でも FFS でも設定しないまま、運動量を使う計算 (追跡・EMIT・Emittance・RADINT・CALC6D・DA ほか) に入ると 1 度だけ警告 | 値は変えていない (30e9)。明示の `MOMENTUM=30e9` では警告しない | `test_b67_momentum_warning` | 警告なし |
| NOKEEPEXP | `NOKEEPEXP` のとき、EXPAND/NOEXPAND を付けない CALC・GO は `calexp = .not. tfvcomp()` | 同じ (tffsa.f:553, 569 を写した) | gs の実測 10 行 (K1・BX・フラグ) と一致 | `test_nokeepexp_matches_gs` | 同じ |
| RENUMBER の INDEX | BEND・QUAD・SEXT… の `INDEX` 属性が 2 以上なら、同名の要素の番号はその値から (QF の INDEX=5 → QF.5 QF.6 …)。USE でも RENUMBER でも | 同じ (tfinit.f tfinimult を写した)。`LINE["INDEX",…]` も読める | gs の `LINE["NAME"]`・`LINE["MULT"]`・`LINE["POSITION",…]`・RENUMBER 後の名前と一致 | `test_renumber_with_index_attribute_matches_gs` | 同じ |
| 文脈 | `BeginPackage`・`Begin`・`End`・`EndPackage` で名前空間を作る | 同じ規則で名前を解決する (§2.4)。`Names["Stat`*"]` は v2 がその文脈の記号の一覧、k64 が gs と同じ `{"Stat`"}` | gs の実測 (15 行の調査スクリプトと example の出力) | `test_context_package_file`、`test_context_resolution_rules`、`test_context_package_using_ffs_and_builtins`、`test_chromaticity_scan_in_user_package` | 同じ (Names だけ切り替え) |
| Shared | fork した子と共有するメモリ | 無名の mmap に値を符号化して書く。大きさの判定 (`Shared::toolarge`)・`SharedSize` は K64 と同じ式 | gs: `SharedSize[{1,2,3}]` = 32 | `test_shared_memory_in_process`、`test_fork_and_shared_on_posix` | 同じ |
| Fork | プロセスを分ける | POSIX では本当に分ける (`Wait`・`Wait4`・`Kill`・`Sleep` も)。Windows は v2 でエラー (§2.5) | gs の example と同じ結果 (子が書いた {10,20,30} を親が読む) を WSL で確認 | `test_fork_and_shared_on_posix` (POSIX)、`test_fork_on_windows` | Windows: 従来どおり pid を返す (子は作らない) |
| System | 任意のシェルコマンドを実行。戻り値は起動できれば常に 1 | 既定では実行しない。許可したときだけ実行し、終了コードを返す (§2.6) | SADHelp「returns the system return code」、Packages が `System[…]==0` で成功を判定している | `test_system_is_disabled_by_default`、`test_system_runs_only_with_permission` | 許可なし: 従来どおり 0。許可あり: K64 と同じ 1 |
| Out[n] | `;` なしの式文の値を `Out[n]:=` で表示。未定義の記号 1 つの文・`QF.3 K1 0.9;` などでは前の Out の値を繰り返し表示する | 変更なし (以前からの表示: `;` なしの最上位の式文だけ) | gs の繰り返し表示は K64 の癖と判断し、写さない (§3) | (既存の `test_scriptgs.py`) | 同じ |

## 2. 補足

### 2.1 LegendreP の奇数 m の符号 (pysad の誤り)
`LegendreP[n,m,z]` の整数 m の経路は scipy の `lpmv` に `(-1)^m` を掛けていたが、`lpmv` は Condon–Shortley の位相つきで、D.n の式
`((1+z)/(1-z))^(m/2) 2F1reg(-n, n+1; 1-m; (1-z)/2)` の m → 整数の極限と同じ符号である。奇数 m で符号が逆になっていた
(`LegendreP[3,1,.4]` が -.2749…、正しくは +.2749545416973507)。gs はここで NaN (B45) なので、参照が無く気づけなかった。
式の m = 1 ± 1e-6 の値 (.27495157…、.27495750…) と閉じた式で確かめて直した。モードによらない。

### 2.2 PolyLog の非整数次・大きい |z|
従来の μ = ln z の展開は |μ| → 2π で収束が遅く、z < -12・z > 60 で精度が落ちていた (z = -140、s = -1.5 で相対誤差 1。|z| ≥ 150 は NaN)。
`specfn._polylog_far` を足した: s < 2.2 は反転公式 (Hurwitz ζ を Euler–Maclaurin で複素の a に対して計算)、s ≥ 2.2 は積分表示
(z < -1 は Fermi–Dirac 積分、z > 1 は極を避けた半直線 t = r e^{iθ}、tan θ = π/ln z)。z < -1 は実数を返す。
- 使える範囲: 実数の z、実数の非整数 s。確かめた精度は s ≤ 7.9、|z| ≤ 1e12 で相対 1e-11 以内 (z = 40・80 で従来の展開と 1e-15)。
- 使えない範囲: 複素数の z (SADError)。負の非整数 s が大きい (s < -10 程度) ときの精度は確かめていない。

### 2.3 複素の FindRoot (v2)
初期値に複素数があると、式を変数について正則と仮定して、前進差分のヤコビアンと最小二乗 (numpy lstsq) の Newton 法 (歩幅を半分にする減速つき) で解く。
- 使い方: `FindRoot[z^2+1==0,{z,1.+1.*I}]` → `{z->(…+I),Residual->0}`。連立も可 (`{z*w==I,z+w==2}`)。`MaxIterations`・`AccuracyGoal` が効く。
- 使えない範囲: 変数の範囲 `{var, ini, {min, max}}` は複素では無視する。正則でない式 (Abs、Conjugate、Re、Im を含む) は収束を保証しない。
  実数の初期値で複素の値になる式は、従来どおり `General::wrongval`。

### 2.4 文脈 (BeginPackage / Begin)
K64 (init.n の定義と gs の実測) の規則を写した: 名前は「現在の文脈 ($Context) → $ContextPath の順」に探し、無ければ現在の文脈に作る。
`BeginPackage[c1,c2,…]` は $ContextPath を `{c1,c2,…,System`}` にし、`Begin[c]` は $Context を c にする。`` `x `` は大域 (`` ` ``) の x。
文脈の中の記号は `Stat`名` という名前になり、外からは `Stat`名` で呼ぶ。`$Context`・`$ContextPath`・`Context[sym]` が読める。
- 使い方: `example/v2_B/B2_package_context.sad`。大域に公開する関数は `` `Name[…] := … `` と定義する。大域の変数を使うときは `` `var `` と書くか、
  `BeginPackage[My`, `]` のように $ContextPath に大域を入れる。組み込み関数・FFS の関数 (Twiss、LINE、FFS[…]、DP0 など) はそのまま使える。
- SAD 自身の Packages (`.n`) は従来どおり文脈を無視して大域で読む (自動読み込み `load_package` と、`PackagesDirectory` の下のファイルの `Get`)。
  GUI・CorrectOrbit などの読み込みは変わらない (`PYSAD_SAD_PACKAGES` を設定して、GUI・描画・パッケージ系の試験 test_tclbridge・test_plot・test_sadclass・test_gui_lang・test_corr ほかが通ることを確認)。
- K64 との差 (意図的): (1) pysad の System` は「組み込み関数 + pysad が記号として扱う名前 (`AccuracyGoal`・`Plot` などのオプション名を含む)」で、K64 より広い。
  K64 ではパッケージの中の `AccuracyGoal` は `Stat`AccuracyGoal` になり、`` `AccuracyGoal `` と書かないとオプションとして効かない。pysad はどちらの書き方でも効く。
  (2) 表示: K64 は $ContextPath に無い文脈の記号に文脈名を付けて表示する (大域の記号が `` `sq `` になるなど)。pysad は文脈つきの名前だけ、$ContextPath にある文脈名を省く。
  (3) `Symbol["Stat`x"]` は K64 が値を返し、pysad は記号を返す (pysad の `Symbol` の従来の動作)。
- 使えない範囲: `$Context = …`・`$ContextPath = …` への代入、`WithContext`、文字列を `ToExpression` するときの文脈は「その時点の現在の文脈」。
  入れ子の文脈 (`A`B``) は名前としては通るが、確かめていない。

### 2.5 Fork・Shared・Wait
- 共有メモリ: `s = OpenShared[n]` (n バイト)、`Shared[s] = 値`、`Shared[s]`、`Shared[s]++`、`Shared[s] += v`、`Close[s]`。値は数・文字列・記号・リスト・式・規則。
  `SharedSize[値]` が n を超えると `Shared::toolarge` (K64 と同じ式: リストは (1 + 要素数) × 8 + 中身)。
- `Fork[]`: Linux・WSL・macOS では `os.fork` で本当に分ける。親には子の pid、子には 0。子は `Exit[]` か最上位の文の終わりで `os._exit` する
  (呼び出し元の Python に戻らない)。`Wait[]` は `{pid, status}`、`Wait4[pid, flags…]` は `{pid, 終了コード, Exited->…, …}` (K64 と同じ並び)、
  `Kill[pid]` はこの pysad が作った子だけに送れる (安全のための制限)、`Sleep[t]` は本当に待って 1 を返す。`Pause[t]` は子がいる間だけ待つ (従来は待たない)。
- Windows には fork が無い。v2 は `System::error: Fork[] is not available on this platform …` を出して文を打ち切る (従来は pid を返して子の分岐が黙って実行されなかった)。
  並列の計算は `Parallelize[f, n, lshare, npara]` (pysad は逐次で `{f[1],…,f[n]}` を返す) を使う。
- 確かめた範囲: WSL (Ubuntu、Python 3.14、numba は代用の空モジュール) で example と 3 個の子の試験。ラティスを使う重い計算を子で走らせる使い方 (numba のスレッドと fork の併用) は確かめていない。
  子の `Print` は親と同じ出力先に出る (出力先が Python の StringIO のときは親に届かない)。

### 2.6 System
`System[cmd]` は任意のシェルコマンドを実行するので、既定では実行しない。
- 許可: 環境変数 `PYSAD_ALLOW_SYSTEM=1` (pysad を起動する前に設定する)、または Python から `pysad.miscfn.allow_system(True)`。SADScript の中からは変えられない
  (環境変数は起動時の値だけを見るので、スクリプトの `SetEnv`・`Environment` で書き換えても許可にならない)。
- 許可があるとき: シェルで実行する。作業ディレクトリは `Directory[]`、標準出力・標準エラーは出力先に出る、標準入力は空。
  シェルは OS の既定 (Windows は cmd.exe、ほかは /bin/sh)。`PYSAD_SYSTEM_SHELL` に sh・bash の実行ファイルを書くと `<shell> -c cmd`。
  戻り値は v2 が終了コード (成功 0)、k64 が 1。
- 許可が無いとき: v2 は 1 行の案内を出して `$Failed` を返す。`[ -e path ]`・`[ -d path ]`・`[ -f path ]`・`[ -r path ]`・`[ -s path ]`・`test …` の形だけは
  シェルを使わずに判定して 0 (真)・1 (偽) を返す (KBFrame などの `System["[ -d dir ]"]==0` のため)。k64 は従来どおり「実行しません」を出して 0。
- 確かめた範囲: Windows (cmd.exe) と WSL (sh) で `echo … > file`、`exit 3`、標準出力。SAD の Packages が呼ぶ個々のコマンド (ImageMagick の `import`、`latex` など) は確かめていない。

### 2.7 MOMENTUM の警告 (B67)
既定値のまま (MAIN にも FFS にも MOMENTUM が無い) かどうかは、既定値を印つきの実数 (`ffscmd._DefaultMomentum`) にして見分ける。
CALC (4D の光学) と幾何では出さず、`FFS.momentum()` を呼ぶ計算 (追跡、EMIT、Emittance、RADINT、CALC6D、DA、Touschek、SynchroBeta) で 1 度だけ出す。
SADScript の式の中で `MOMENTUM` を直接使う計算 (例: `MOMENTUM/MASS`) では出ない。`OpticsPlot` など MOMENTUM を一時的に書き換える関数のあとは出なくなることがある。

### 2.8 評価器に足した補助
- `Compat$V2[]`: 互換モードが v2 なら 1 (prelude の PolynomialFit が使う)。
- `Interp._ctx_*` (文脈)、`Interp._pkg_loading` (SAD の Packages を読んでいる間)、`sadscript._exit_child` (Fork の子の終了)。
- `ffscmd.py`: `_exec_word` の中の `from . import compat as _compat` (ほかの担当の MEASURE の行) を外した。関数の中で import し直すと、同じ関数の前の方の
  `_compat` の参照 (B22) が UnboundLocalError になるため。モジュールの先頭で import してある。

## 3. 要判断・直さなかったもの

| 項目 | 判断 | 理由・調べたこと |
|---|---|---|
| B24 `LINE["LENGTH",loc]` | 値は返さず、案内だけ足した | K64 のソース (tftwiss.f:166, 508, 749) が `LENGTH` の引数を 1 個に限っている。`LINE["LENGTH"]` は要素数で、位置を付けたときの意味 (位置か、軌道長か) がソースからもヘルプからも決まらない。位置は `LINE["POSITION",loc]`、軌道長は `LINE["LENG",loc]` で取れる |
| B22 MON・MCOD・TCOD の機能 | 定義・補正は実装しない (警告して読み飛ばすだけ) | K64 にコマンドが無く、古い SAD の MON の書式・意味を確かめられる資料 (K64 のソース・SADHelp) が無い。軌道補正は CorrectOrbit パッケージ (docs/34) で動く |
| `Out[n]:=` の gs の表示 | 写さない | gs は、未定義の記号 1 つの文で前の `Out[n]` を繰り返し、`QF.3 K1 0.9;` で `Out[3]:= 3`・`Out[4]:= .3` を 2 回ずつ出すなど、値と関係のない行を出す (2026-11-17 に gs で確認)。K64 の表示の癖で、利用者の役に立たない |
| NOKEEPEXP の `tfvcomp` | 常に偽として実装 | `tfvcomp` は「成分ごとのフィット変数 (ivcomp≠0) があれば真」。pysad の `FREE` はファミリー単位で、成分ごとの変数を持たないので真になる場合が無い。成分ごとの FREE を実装するときに `FFS._tfvcomp` を直す |
| RENUMBER 後の名前の引き方 | pysad は `LINE["NAME"]` の表示と一致する位置を返す | gs は RENUMBER のあと `LINE["POSITION","QF.1"]` が `LINE["NAME"]` の表示と食い違う (INDEX なしで QF.1 が 2 番目を指すなど)。gs の内部の表の更新の癖と見て、写していない (INDEX の有無によらず以前からの差) |
| 複素引数の PolyLog・HypergeometricU・Zeta | 未対応のまま (SADError・TypeError) | scipy に複素の実装が無い (hyperu) か、分岐の向きを確かめる独立な参照が取れていない |
| `Parallelize` の本当の並列 | 逐次のまま | POSIX では Fork で実装できるが、numba のスレッドと fork の併用を確かめていない。結果は逐次と同じ |
| `Symbol["x"]` の値 | そのまま | gs は値を返し pysad は記号を返す (文脈の調査で見つけた差)。担当の項目の外 |
| `Message[f::usage]` の表示 | そのまま | gs は `???f::usage:  本文` を表示し、pysad は何も表示しない (B70 の調査で見つけた差)。担当の項目の外 |

## 4. k64 モードで以前と変わった点 (意図したもの)

k64 モードは「K64 (gs) の動作」なので、次の 2 つは gs に合わせて変えた (ほかは以前と同じ値):
- `For` の本体を実行しない (B71。以前の pysad は実行していた。SAD の Packages は `For` を使っていない)。
- 文脈: 利用者のコードの `BeginPackage`・`Begin` が名前空間を作る (以前は無視して全部大域。SAD の Packages の読み込みは変えていない)。

pysad 自身の誤り・精度の修正で、モードによらず値が変わったもの: `LegendreP` の奇数 m の符号、非整数次の `PolyLog` の z < -12・z > 60
(以前は誤差が大きいか NaN)、`SharedSize` (以前は 1、gs と同じ 32 に)、`LINE["INDEX",…]` (以前は 0)、INDEX 属性つきの要素の番号、`NOKEEPEXP`、
`h::tag` の `General::tag` への継承、複素引数の Bessel・超幾何 (以前は TypeError)、POSIX の `Fork` (以前は分けない)、`Sleep` (以前は未定義)。
