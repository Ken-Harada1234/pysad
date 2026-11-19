# 46. オリジナル SAD (K64) に対する pysad の未対応・未確認の一覧

最終更新: 2026-11-20 (v2.2.0 の現状のスナップショット)

この文書は、オリジナルの SAD (K64 版、`src/*.f`・`Packages/*.n`・SADHelp) に対して、pysad が**まだ対応していないもの**と**対応したが gs で確かめていないもの**を、1 か所にまとめたものです。

- **互換モード (v2.0.0 から)**: pysad は 2 つのモードを持つ。`k64` は K64 (gs) と同じ動作 (v1 と同じ。gs との数値比較はこのモード)、`v2` (既定) は K64 のバグを直し機能を足した動作。環境変数 `PYSAD_COMPAT=k64|v2` か `pysad.compat.set_mode()` で切り替える ([47_v2_changes.md](47_v2_changes.md))。**各項目の「モード」の欄が `k64` ならそのモードだけの話、`v2` なら v2 モードだけの話、`両` ならモードによらない。**
- v2 の値の根拠は、解析式・独立な数値計算 (ローレンツ力・Thomas-BMT の数値積分など)・保存則・別のコードとの一致で、gs ではない。**v2 で実装・修正した機能の多くは、gs が動かない・K64 が誤っているので gs と比較できない。** 「未確認」にはその旨を書いた。使った条件は電子、0.2〜4.2 GeV のリングと輸送路 (空間電荷は 5〜50 MeV、RADLIGHT は 10.2 MeV)、ソレノイド bz = 0.3〜3 1/m (一部 20 1/m まで) で、**実機の強い結合・挿入光源のずれのあるリングでは確かめていない**。大きなリング A (3 GeV、周長 603 m、`pysad/tests/kek-ls_v03_12.lat`)・PF リング・ASTRID2 (MPW 入りを含む) では多くの項目を gs と比べた。大きなリング B (8 GeV、周長 1436 m) は数値だけ。
- 出所: [16_fortran_status.md](16_fortran_status.md) §9、[history/work_unimpl.md](history/work_unimpl.md)、[21_sad_bug_report.md](21_sad_bug_report.md)、v2 の作業の記録 `history/v2_changes_*.md` (v2.0.0) と `history/v2next_*.md` (v2.1.0: A1・A2・B1・B2・B3・C・D1・D2・E、v2.2.0: P1〜P5)、[47_v2_changes.md](47_v2_changes.md) の v2.2.0 の (d) 要判断・未解決、`example/` の各例の「gs と pysad の比較」・「未確認・限界」節、[34_sad_packages.md](34_sad_packages.md) ほか各解説文書。
- 基準は gs (SAD K64) の実動作です。gs が落ちる・誤った値を返すものは「gs 側の問題」(§5) に分けています。K64 のバグの番号は [21_sad_bug_report.md](21_sad_bug_report.md) の B01〜B159 (B72〜B110 は v2.0.0、B111〜B134 は v2.1.0、B135〜B159 は v2.2.0 の作業で見つけたもの)。
- `例` の欄は、その項目を示す `example/` の例です (パスは `example/` からの相対)。v2 の機能の例は `example/v2_B/`・`example/v2_C/`・`example/v2_next/`。

凡例: 「未翻訳」は pysad に物理が入っていない。「一部未対応」は一部の入力・経路だけ未対応。「未確認」はコードはあるが gs と数値で比べていない。「モード」は上のとおり。

## 1. 未翻訳

| 項目 | モード | 理由・状況 | 例 |
|---|---|---|---|
| `wjfit`、`wfres`、`wfbin`、`wiord`、`wjdraw` (MAIN の TRACK の Fourier 解析、FOURIER フラグ) | 両 | FFS から使えない。gs も SIGSEGV (`tracka` → `tinip`) で出力が取れず、直訳しても検証できない。pysad は FOURIER を無視し `wtunem` の表だけ出す。チューンの抽出は NAFF (`pysad/naff.py`) を足した | v2_next/B3_naff_tunes.py |
| 光子 (`tphrec`) のうち SOL 領域内の BEND・SEXT・OCT、COORD・TCAVI・MAP のあるラティス | 両 | BEND・QUAD・MULT・SEXT・OCT・磁場つきドリフト・BEAMBEAM + BEND は記録する (v2.1.0)。**SOL 領域内の DRIFT・QUAD・MULT (ANGLE=0) は v2.2.0 の v2 モードで記録する** (位置が軌道上、エネルギーの合計が古典値と統計の範囲。3 GeV、bz = 0.1 1/m、傾き・ずれのない領域)。領域内の BEND・SEXT・OCT、傾いた・ずれた領域の光子の位置 (`geo_frac` が領域の軸の傾きで別の系になりうる)、偏光 ξ の値は未対応・未確認。k64 は K64 のまま (領域内の QUAD は記録しない: B150。gs と比べられない)。COORD・TCAVI・MAP のあるラティスは従来の高速追跡で BEND があっても光子を記録しない | v2_next/B3_photons_bend.sad、P4_sol_photons.sad |
| TrackParticles の F1 の残り (`tsteer`・`tmulti`) | 両 | v2.2.0 で `tsteer` の F1・FB1・FB2・FRINGE と `tmulti` の F1・F2・F1K1F/B・FB1・FB2・FRINGE を実装した (§6)。残り: 逆向きの要素 (`direlc < 0`) の F1 の入れ替え (BEND の F1・FB1・FB2、MULT の入口・出口。Emittance[] も同じ)、RF つき MULT (`tmultiacc`) の F1 (NotImplementedError)、POL (スピン) つきの MULT の F1 の端の項、F1 つきの MULT の追跡の高速化 (numba の MULT カーネルが F1・FRINGE・FB・AE に未対応で Python の経路に落ちる: 1 周 1 粒子 0.4 s 程度)。K64 の F1 の模型は fb1 ≠ fb2 で x が厳密な Maxwell 場の端の運動と (K0/al)(fb1² − fb2²)/24 違う (B145。直さない) | v2_next/P3_mult_ae_emit.sad |
| CALC (4 次元) の RADTAPER | 両 | **実装しない**: K64 の 4D の CALC は RADTAPER を見ない (gs で確認)。CALC6D は RADTAPER を入れる | 04_emittance_radiation/flag_RADTAPER |
| MON・MCOD・TCOD の機能 | 両 | K64 にコマンドが無く、古い SAD の MON の書式・意味を確かめられる資料が無い。v2 は「K64 に無いコマンド」と警告して文全体を読み飛ばす (定義はしない)。軌道補正は CorrectOrbit パッケージ ([34_sad_packages.md](34_sad_packages.md)) | |
| `DRAW` の中の粒子追跡 (`tdrwdt.f`) | 両 | K64 でも動かない死んだコード (`trsize` が常に偽。B105) | |
| `Out[n]:=` の gs の表示 | 両 | gs は未定義の記号 1 つの文で前の `Out[n]` を繰り返すなど値と関係のない行を出す。K64 の表示の癖として写さない | |
| 複素の特殊関数で gs が誤る点の再現 (B118〜B120・B158) | k64 | k64 モードも真値を返す (gs の誤った値を再現しない)。再現には K64 の `gamma.f` (複素の不完全ガンマ `cgamma2`・`czconvb`・`cpolylogl`) と D.n の式の直訳が要る。要判断。v2.2.0 の走査 (732 式) で gs が誤る・不正確な式は 91 (LogGamma の分枝・Erf/Erfc の Re z < 0・ExpIntegralEi・ArcCosh の符号・FresnelS/C・BesselY・BesselK ほか)。規則が単純な LogGamma・ArcCosh・Erf/Erfc・ExpIntegralEi は k64 で再現できる可能性がある | v2_next/B2_special_functions_complex.sad、P5_specfn_cx_scan.sad |
| `DynamicApertureSurvey` の NPARA > 1・`Damping->True` | 両 | v2.2.0 で DAPWIDTH の窓 (既定 7) は `trackd.f` の直訳にした (§6)。残り: **NPARA > 1** (gs は粒子を子プロセスに分け、子ごとに乱数の種をずらす。再現しない。RAD + FLUC の最終 seed の比較は NPARA = 1)、`Damping->True`・`DampEnough`。v2 は放射の式が K64 と違うので RAD + FLUC の seed は gs と合わない (NORAD の地図は合う)。損失の境界の粒子は gs と 1〜数周違うことがある | v2_next/P4_dapwidth_window.sad |
| `TRCOD` の FIFO での待ち合わせ (`PUTMEASURE`・`GETMEASURE`)・`ORBITJITTER` のパイプ | v2 | 普通のファイルで同じ形式を読み書きする (待たない)。Linux・WSL の FIFO は未確認 (試験は普通のファイルだけ) | |
| K64 の `qwspac` (4 次元の CALC の空間電荷) の直訳 | v2 | v2 は CALC6D と同じ 6D の経路で入れる。gs が NaN を返す機能 (B133) なので直訳の検証ができない | v2_next/B1_wspac_emittance.sad |

## 1a. モードで状態が違うもの (k64 モードでは未対応・v2 モードで実装済み、または結果が違う)

v2 モードの実装は、特記がなければ**gs と比べていない** (gs が動かない・K64 が誤っているため)。根拠は解析式・数値積分。v2.1.0 で k64 と v2 の結果が変わった項目の一覧は [47_v2_changes.md](47_v2_changes.md) の v2.1.0 の (a) (B111・B112・B114〜B117・B121〜B123・B125〜B127・B129〜B133 ほか)。

| 項目 | k64 モード | v2 モード | 例 |
|---|---|---|---|
| MEASURE・TRCOD (FIXCOD・MVCOD・STNCOD・ALIGN)・LOGMEASURE・ORBITJITTER・PUTMEASURE・GETMEASURE、MAIN の `TRACK TURNS=-1` (`trackb.f`・`tracke.f`・`trcoda.f`・`ttinit.f`・`ttstat.f`) | k64: 先頭の MARK の位置の MEASURE・TRCOD・LOGMEASURE だけ実装 (v2.2.0。gs と統計が 4 桁一致)。それ以外は警告して無視。gs は先頭以外の位置で全粒子を落とし MEASURE は無出力、TRCOD は SIGSEGV (B42・B151) | 実装 (K64 のソースの読み。**v2 の追跡を通した値は gs と比べられない** (gs が追跡で落ちる)。解析式と pysad の独立な経路で確認: docs/03、C の記録 `history/v2next_C.md` §3.1) | example/v2_C/C2m_trcod_reference.sad、v2_next/C3_trcod_orbit.sad、P5_measure_*.sad |
| SOL 領域内の BEND・MULT・CAVI・APERT (`tsole1`・`tmulte1`)、SOL の `DX DY DZ CHI1..3`・`F1`、領域の幾何 | 未対応 (NotImplementedError)。BEND は K64 の tsol.f どおりに実装したが gs が偏向をほぼ無視して合わず取り下げ (B69) | 実装 (ローレンツ力の数値積分との一致。docs/44 §7a) | example/v2_C/C3_sol_region.sad、C2a_sol_tilt.sad |
| 領域内の CAVI・RF つき MULT の TRPT の加速・wake・`AUTOPHI` | 未対応 | 実装 (加速は運動方程式の積分と一致。AUTOPHI の意味は K64 のソースの読みで、gs と一致するのは TRPT の場合だけ) | example/v2_C/C2b_sol_rf.sad |
| 空間電荷 (SPAC・PSPAC・SELFCOD) の追跡 | 追跡に入らない (K64 の SPAC は力が 1/γ 倍に弱い: B89 など。ビット単位で再現する翻訳はしていない)。WSPAC は Emittance[]・EMIT に入る (gs と全桁)。**PSPAC は v2.2.0 から K64 の形状関数で入る** (gs と rms が 12 桁以上。NPARA = 1。B147・B148) | 実装 (SPAC は解析式・Laslett の式、PSPAC は gs と 3〜4 桁、WSPAC は gs と 1%。粗いメッシュの PSPAC は警告。WSPAC の CALC・CALC6D にも入る) | example/v2_C/C2s_mult_radius.sad、C2s_wspac.sad、05_tracking/flag_SPAC |
| `PHSROT` | 追跡は NotImplementedError、光学は恒等 (K64 で動作しない: B38) | 実装 (意図された写像。固有値・シンプレクティック性・平衡エミッタンスで確認) | example/v2_C/C2_phsrot.sad |
| UND の磁場 (`FBX`・`FBY`) | 磁場なしのドリフトのみ (K64 が設定できない: B37)。光学は恒等 | 実装 (解析式で確認)。磁場がなくても光学にドリフトとして入る (要判断: K64 の「光学では恒等」を仕様と見るなら磁場つきのときだけにする選択もある) | example/v2_C/C1_undulator.sad |
| `Polarization` (平衡偏極度) の有限値 | K64 が NaN (B55)・値が大きすぎる (B125)。逆向きの BEND (挿入光源) で緩和率が相殺する (B116) | 有限値 (平面リング 92.376%。誤差のあるリングは側帯の項とクロマティックな変調つきで FLUC の追跡と 4.2 GeV・12 MV で 1.05 倍: §2。逆向きの BEND は符号に依らない率。n が傾く・非平面のリングは厳密な緩和率と運動論的偏極)。docs/43 | 07_spin_solenoid_special/flag_POL_RADPOL、v2_next/B3_pol_sideband.sad、A2_pol_mpw.sad、P1_pol_*.sad |
| Emittance[]・EMIT の POL の量 (ずれて傾いたソレノイド) | K64 の式 (縦磁場を z 軸の向きに足す: B130)。SpinTune が追跡と 1.5e-4〜1.1e-3 違う | 運動量の向きに分けて足す (追跡と 1e-13〜5e-7) | v2_next/E_sol_tilt_spin.py |
| スピン追跡 (POL + RAD) の損失の処理、領域内の BEND・MULT・CAVI・ずれ・F1、UND、PHSROT、TRPT の加速 | 損失の処理なし。領域内は DRIFT・QUAD だけ。UND・PHSROT は NotImplementedError。TRPT の加速は明示のエラー | 損失を処理 (LOSSMAP の 2 行も付く)、領域内・UND を通す (PHSROT はスピンを動かさない)。**POL + TRPT + 加速は SOL 領域の外の CAVI だけ** | example/v2_C/C2c_spin.sad、v2_next/D1_spin_accel.py |
| RADLIGHT の軌道の横方向・BEND の中 | 横方向の座標が入らず BEND の中を記録しない (K64 の癖。gs と差 0)。放射場の計算は意味を持たない | 3 次元の軌道、BEND の中も記録 (B65)。複数周の軌道でも `RadiationField` が動く | example/v2_C/C5_radlight.sad |
| 旧式の `WAKE … ENDWAKE` の表 | `WAKECOMMAND[]` を評価して失敗 (K64 と同じ) | `WakeFunction` の表として登録 (CAVI だけ。行をまたいで書ける) | example/v2_C/C2w_wake_table.sad |
| `System`・`Fork` | `System` は実行せず 0、`Fork` は pid を返すだけ | `System` は許可 (`PYSAD_ALLOW_SYSTEM=1`) があるときだけ実行し終了コードを返す。`Fork` は POSIX で本当に分ける (Windows はエラー) | example/v2_B/B3_shared_system.sad |
| 複素初期値の `FindRoot` | TypeError (gs は壊れた値: B61) | 複素の Newton 法 (変数の範囲は無視、正則でない式は収束を保証しない) | 10_math_data/fn_Gamma_Erf |
| `MessageName[f,"usage"]` の関数形 | pysad は以前から f を見る (gs は f を無視する: B70) | f の本文を返す (無ければ `General::tag`) | 13_packages/fn_OptimizeOptics ほか |
| `Union[..., SameTest->f]`・`Names["組み込み"]`・`WriteBeamLine[..., Name->"L2"]` (文字列)・`Begin` に文字列の文脈 | gs と同じ (B121・B122 の誤りを含む)。`Begin["A`"]` は拒否 | 真値の `Union`、組み込みも含める `Names`、文字列の Name・文脈も受ける | v2_next/B2_packages_compare.sad、B2_context.sad |

## 2. 一部未対応

| 機能 | モード | 未対応の部分・使えない範囲 | 例 |
|---|---|---|---|
| MULT `ANGLE`≠0 | 両 | Twiss の `BZ`・`RF`・`F1`・`SK0`・`FB1/FB2` (追跡は `FB1/FB2` まで対応) | 07_spin_solenoid_special/elem_MULT_ANGLE |
| MULT の端の `FB1/FB2` の `tblfre` (`multae.py`)・逆向きの要素の F1・AE | 両 | **AE1・AE2 は v2.2.0 で実装** (CALC・Emittance[]・追跡。k64 は gs と全桁一致。§6)。残り: Emittance[]・CALC の MULT の `FB1/FB2` の `tblfre` は NotImplementedError、逆向きの要素 (`direlc < 0`) の F1・AE の入れ替えは未対応。K1≠0・K0≠0 の ANGLE=0 の BEND の追跡の AE の K64 の値は gs が SIGSEGV (B02) で確認できない (ソースの読み: B142)。実験室系で y に傾く端の面 (K1 < 0 の MULT・ROTATE のある MULT) の厳密な幾何との比較は無い (K64 の約束は面が要素と一緒に回る)。POL (スピン) つきの MULT の AE の bsi の端の項は未確認。AE のある極の放射の残る 2.7% (B143) | v2_next/E_mult_f1.sad、P3_mult_ae_emit.sad |
| MAP | 両 | `GEO`・`CompiledMap`。MAIN の TRACK・DA・Touschek は ExternalMap を呼ばない。MAP のスピンは未対応 | 06_collective/fn_ExternalMap |
| POL (SOL 領域内) | k64 | 領域内の MULT・BEND・CAVI のスピン (Emittance[]・追跡とも)。スピン追跡は損失の処理なし | 07_spin_solenoid_special/flag_POL_RADPOL |
| POL (SOL 領域内) | v2 | 領域内の APERT のスピン追跡 (NotImplementedError)。領域内の CAVI の POL (Emittance[]) は VOLT = 0 と 0.2 MV で追跡と 1.6e-7 (位相 0 なので弱い試験。加速・位相が 0 でないときは未確認)。**領域内の QUAD・MULT のスピンは v2.2.0 で Thomas-BMT の積分と 9e-7 以下** (bz = 0.3〜9 1/m、0.4〜0.5 m、DISFRIN=1。QUAD は bz·al/ndiv ≤ 0.02 に刻み、端の縦磁場の積分 ψ は DISFRIN=0 のときだけ: B146)。**DISFRIN=0 (フリンジあり) の ψ を、硬い端の B_s と軌道のフリンジ写像を合わせた厳密な基準で確かめていない**。ずれ・ROTATE を持つ MULT の ψ、大きなリング A・PF・8 GeV は未確認 (0.1〜1 GeV・ASTRID2 までのリング)。領域内の四極の DX・DY のある POL の追跡 (track6) は NotImplementedError。領域内の光子は DRIFT・QUAD・MULT だけ記録する (上)。大きなリング A への SOL の確認は小リングのみ。領域の外の B∥ (垂直な軌道の傾きの中の偏向磁石の B_y y') は K64 も追跡も入れていない (B136) | |
| POL + TRPT + 加速 (track6) | v2 | SOL 領域の外の VOLT≠0 の CAVI だけ。**SOL 領域内の加速は実装しない (v2.2.0 の判断)**: 領域内では刻みごとに基準運動量 p0(s) を進め、ソレノイドの bz を p0/p(s) で直し、スピンの歳差の aγ(s)・運動量の向きの磁場の積分を p(s) で作り直す必要があるが、RF の電場を含む運動方程式の積分 (E_z(t) cos ωt) の独立な基準が無い。実装するなら `sol_mult_full` の刻みで `_RadHook` が基準運動量を更新し `_tradk_sol` に p(s) を渡す。RF つき MULT・TCAVI の加速・リングの何周もの TRPT (K64 の TRPT は輸送路用)・スピン追跡の損失の処理は SADError。現状は明示のエラー (`test_pol_with_trpt_acceleration_is_an_explicit_error`) | v2_next/D1_spin_accel.py |
| `Polarization` | v2 | 軌道の振幅について 1 次の SLIM に、スピンの歳差の位相の変調 (シンクロトロン側帯。B125) と、**ベータトロンのチューンの δ 依存 (クロマティシティ) による側帯の重みの変化** (v2.2.0、`spin.CHROMATIC`) を入れた。4.2 GeV・12 MV の 12 セル FODO で式 / 追跡 (FLUC) = 1.05 ± 0.16 (v2.1.0 は 0.43)。2.903 GeV・1.5 MV 1.22 (統計 ±24%)、ASTRID2 0.99、3 GeV・1.5 MV 1.8 (追跡の統計 ±22%。式が高い側 = P を低く出す安全側)。**追跡の測定は 4.2 GeV・12 MV の 1 点が主**: 6 MV・24 MV・PF・大きなリング A は式の値だけ (大きなリング A・PF はクロマティックな変調を入れない: 変調指数 0.008〜0.088 か、要素数 > 1200・K > 6e4)。式のモード分解は決定論的な応答を追跡の 1.5 倍に見積もり、結合項 b(u) の 2 次の欠落 (−10%) と打ち消している可能性がある。六極の振幅依存のチューンシフトと、六極で補正した小さな ξ の模型は未検証。目安: 変調指数 ≥ 0.03 かつスピンチューンがベータトロンの側帯 νx + kνs の近く (離調 ≲ 2νs) のリングで B3 の式の 2 倍前後。そうでないリングは B3 の式と ±5%。限界の表は docs/43。不安定なリングは **P = NaN** (v2.2.0: 1 周の写像の固有値の大きさが 1 を超える。以前は線形の式を不安定な軌道の上で解いて 0.92 などを返した)。n が傾く・非平面のリングの緩和率は厳密式 (`spin.ST_EXACT`): 傾いた n の大きなリングの独立な計算は未実施。P の脱偏極率の部分は Emittance[] の D (追跡の D と違う: B136・B137) を使うので、脱偏極率が分母の 15% を超えると P は 2〜3% 違う。運動論的偏極は位置・角度を固定した ∂n/∂δ (SLIM の δ 列。追跡との 13% の差を含む)。`Polarization2/4/6` は実装していない。POL のときは四極・MULT を細かく分割するので Emittance[] が遅い | v2_next/B3_pol_sideband.sad、P1_pol_chromatic.py |
| TCAVI | 両 | `RANKICK`・`RANPHASE`。`Emittance[]` は L=0 のみ gs と一致 (L≠0 は gs が壊れる: B36。pysad は有限値を出す) | |
| UND | k64 | 磁場なしのドリフトのみ (K64 が磁場を設定できない) | |
| UND | v2 | 端の磁極の整形なし (軌道は x' = −(F0/kz) sin(kz s) で中心が −F0/kz² ずれる)、`BX`・`BY` という名前は K64 にも無い、PHI≠0 の入口はベクトルポテンシャルの飛びが蹴りになる (K64 の式のまま)。`SLICE` は `POLE` の 20 倍以上が条件。挿入光源のあるリングでは未確認 | example/v2_C/C1_undulator.sad |
| PHSROT | v2 | スピンを動かさない (リング全体を PHSROT 1 個で表すとスピンチューンは 0 になり偏極の計算に使えない)。拡散行列に負の固有値があると ValueError。分散つきの Twiss を CALC で読むと BX・EX は入力と 1e-3 ほど違う (H 行列の定義の差)。(z, δ) の回転は向きが逆 (縦のチューンは −PZ と出る) | example/v2_C/C2_phsrot.sad |
| `APERT` | k64 | SOL 領域内は NotImplementedError | 05_tracking/elem_APERT_LOSSAMPL |
| `APERT` | 両 | 再入射 (`JDPX`・`JDPY`) は粒子 12 個以上で gs と乱数列が違う (NPARA の分け方は `NPARA=n;` で再現できる)。`LOSSAMPL`・RAD・FLUC つきの TRPT、MAIN の TRACK の損失報告、DAPERT での APERT は gs と未照合 | 05_tracking/elem_APERT_LOSSAMPL |
| LOSSMAP | 両 | MAP 要素 (ExternalMap) が落とした粒子の要素番号は 0。空間電荷が落とした粒子の要素番号も。POL なしで高速追跡に乗らない要素の並びは `LOSSMAP: … 記録できない` の SADError。落ちた要素は判定 (tapert) のときの l_track なので 15 要素ごとの判定点か APERT の番号 | 05_tracking/flag_LOSSMAP |
| RADTAPER | 両 | CALC (4 次元) は K64 の仕様で見ない (§1)。SOL 領域 + RADTAPER は未確認、RFSW OFF (見つからないときは `???-RADTAPER-closed orbit with taper not found` を出してテーパーなしで続ける)。v2 は追跡の BEND の重みを (3δ_l + δ_l+1)/4 にして閉軌道が追跡で閉じる (B111。どちらの重みが物理的かは未判断)。gs と比べるときは gs 側で `CODCONV = 1e-4;` が要る。k64 の CALC6D は RADTAPER を入れる (gs と 8e-8)。`DP0` の別名は k64 が gs と同じ使い方で、v2 は CALC6D が書いた DP0 を RADTAPER の探索に使わない | 04_emittance_radiation/flag_RADTAPER |
| RF のある MULT (`tmultiacc`) | 両 | `F1`・`FRINGE`。`AUTOPHI` は k64 で NotImplementedError、v2 で実装。CAVI の `V1 V20 V11 V02` は追跡で対応 (TENT・TEXIT つき)。WakeFunction のある CAVI・SOL 領域内の CAVI の V1 V20 V11 V02 (NotImplementedError)、Emittance[] の CAVI の V1 V20 V11 V02 (`tcave` は V = 0 を前提) は未対応。領域内の CAVI・MULT: ANGLE≠0 の MULT・領域の外と bz = 0 の領域の RF つき MULT の Emittance[] は未対応 (v2 の C2 の記録に「bz = 0 の領域の中の CAVI・MULT は領域の外と同じ写像」という記述もあり、**記録内の表現に食い違いがある: 要確認**) | |
| SOL 領域 (v2) | v2 | CALC6D・RADCOD つきの CALC・INTRA・Emittance[]・EMIT は対応 (v2.1.0。k64 も gs と 1e-9)。使えない: SEXT・COORD・BEAMBEAM など (`SolRegionError`)、領域内の BEND の `K1`・`DX`・`DY`、追跡を領域の途中から始める・途中で終える (`ソレノイドの終端がない` のエラー)、領域より上流の幾何に軸まわりの回転 (GCHI3≠0) があるラティスで追跡を先頭以外から始めた場合 (CHI3 が入らない)。未確認: TouschekLifetime (gs が NaN)、ソレノイドの軸が設計軌道と 10 mrad 以上ずれる場合、領域内の SOL (BOUND=0) にずれ・回転を書いた場合、GEO のない領域に手で書いた CHI の符号の約束、領域内の空洞の wake の横の蹴り (gs と未比較)、領域内の RF つき MULT の RAD、RADTAPER + SOL、WSPAC・INTRA + SOL の併用の値。ずれた SOL (k64 は SOL 自体のずれが未対応) の幾何の fit。CHI3 を使う幾何の fit は K64 でも収束しない (B124)。放射は bz·L/ndiv ≤ 0.02 に刻んで N5 (B80) の誤差を避けている | |
| SOL (両モード) | 両 | BOUND=0 で BZ≠0 の SOL (領域の外のソレノイド) の Twiss・CALC は NotImplementedError。QUAD の FRMD の向きの入れ替えは未実装。行列は数値微分で解析行列の gs と 1e-10 程度の差。RADCOD の閉軌道探索中の分割の固定は未確認 | 07_spin_solenoid_special/elem_SOL |
| 空間電荷 (v2・k64) | 両 | **PSPAC は v2.2.0 から k64 モードにも入る** (K64 の形状関数: 横は CIC、縦は「ほぼ最近接」、電位は 5 点差分 + 双 3 次補間、線密度は 3 次スプライン。`spch_k64.py`。gs と rms が 12 桁以上一致 (5 条件・多周回 4 条件・ASTRID2 に 3 SPCH で 10 周)。ただし NPARA = 1 で比べた)。k64 は PSPAC だけ入れる (SPAC は B89 の 1/γ・WSPAC の追跡は入れない)。**使えない**: gs の NPARA > 1 の再現 (子プロセスごとの電位: B148)、`PSPACNZ < 3` (k64 はエラー。gs は SIGSEGV)、SPAC・WSPAC との併用、ウェーク・放射との組み合わせ、PF・大きなリング A での PSPAC (ASTRID2 までの確認)。電荷の 2 乗は入れない (K64 は q を使わない)。v2 の PSPAC は粗いメッシュで gs・解析式から離れる (rms x が縦 20・4・2 セルで gs と 0.6%・4.5%・12%)。v2 の SPAC・PSPAC のほか: `ANGLE≠0`・`VOLT≠0` の MULT と SOL 領域内の MULT の SPAC (メッセージを出して空間電荷なしで追跡) など v2.0.0・v2.1.0 の未確認はそのまま | example/v2_C/C2s_mult_radius.sad、v2_next/P4_pspac_k64_shape.sad、P4_pspac_astrid2.sad |
| WSPAC | v2 | **追跡の WSPAC** (`LINE["SIG…"]`) のビーム行列は既定でゼロ電流の平衡値・MINCOUP なし (gs の追跡のビーム行列にも MINCOUP の下限が入らず、ゼロ電流の行列と 1〜2% で合う: B149)。**v2.2.0 のオプション `WSPACSC = 1;`** は空間電荷つきの自己無撞着な行列 (MINCOUP の下限つき) を追跡に使う (K64 にない。チューンシフトが同じ行列の一次の式と 2e-4。0.2 GeV の 12 セル FODO、PBUNCH = 1e6〜1e8。gs の追跡とは合わない)。強い空間電荷 (0.2 GeV で PBUNCH ≳ 1e9) は反復が収束しない、RF なし、大きなリングでの収束・時間は未確認。既定を空間電荷つきの行列に変えるかは要判断 (docs/47 v2.2.0 (d))。ビームの大きさは追跡の間じゅう固定 (強-弱。エミッタンスの増大・コヒーレントな振動は扱えない)。垂直のエミッタンスが 0 の平らなリングは既定では Emittance の経路で使えない (`WSPACBEAM` を与える。`WSPACSC = 1` なら MINCOUP の下限で使える)。`WSPACBEAM` は σz 一定・x-z の相関なし・MINCOUP なし。TRPT の輸送路ではエミッタンスの断熱減衰を入れていない。縦の力は既定で「軸上で 0」(K64 と同じ)。**Emittance[]・EMIT・CALC6D・CALC (RFSW ON) の WSPAC は実装済み (§6)**: 強い空間電荷、RF なしの CALC、TouschekLifetime との併用、要素ごとのレンズの位置 (K64 の tintrb は磁石の中で厳密でない近似)、実機の大きなリングは未確認。SOL 領域内・`SELFCOD` + 閉軌道のあるリング・放射つきの追跡は走るが値の確認なし。gs で `EMIT` を `NOWSPAC` で行ってから `WSPAC` で追跡すると NaN。gs の追跡の WSPAC は `CODPLOT; WSPAC; EMIT;` のあとでだけ働く (`Emittance[]` のあとでは入力をそのまま返す) | v2_next/B1_wspac_emittance.sad、P4_wspac_tracking_beam.sad |
| MEASURE・TRCOD (v2) | v2 | **gs と比べられない** (gs は先頭の MARK 以外の位置で追跡が全粒子を落とす: B151。k64 は先頭の位置の統計が gs と 4 桁一致で、v2 は同じ初期分布を追跡に通す)。C の独立な基準 (閉軌道の式、固有エミッタンス、ビームサイズ) で誤りは見つからなかった。`TRCOD FIXCOD`・`MVCOD`・`STNCOD` の基準は REFERENCE と同じ参照光学に置いた (`measure._ref_rows` の 1 か所で替えられる)。`ALIGN` の対象に DECA DODECA MULT を足した (K64 は BEND QUAD SEXT OCT だけ。`measure._ALIGN_V2` を `_ALIGN_K64` にすれば K64 のとおり)。SOL・MONI は動かさない。区間に BEND を含む ALIGN は平行移動との一致が 4e-10 m まで。GETMEASURE の対象の関数 (DX・DY) は K64 の関数番号 32・34 から確定できず式から判断 (B103)。**MEASURE の入口が DETR ≥ 1 の場合は未確認** (gs は先頭の MARK の位置だけ動き、MARK の生の DETR が gs で効かない: B159)。TWAKE / LWAKE ON は z の分布を作る (K64 は全粒子 z = 0)。MEASURE の TRPT は追跡の間だけ ON にして戻す (K64 は戻さない)。MAIN の `TRACK TURNS=-1` は粒子種の語 (PROTON など) を見ない。TRGAUSS OFF の分布と表示の桁 (autofg) は gs の実出力と比べていない。TRCOD の途中の統計は出さない。周期解を matched 分布として 1 周追う試験は δ が大きいと崩れる (DP = 1e-3 で 40%)。QUIET の重心合わせ (B87) と TWAKE/LWAKE の z = 0 (B88) は K64 の誤りと判断して直した推定。k64 の `FFS["MEASURE"]` の追加の規則 (先頭以外は無出力) は先頭の位置の実験 (4 つのリング) から決めた。TRCOD FIXCOD などの k64 の実装は無い。GETMEASURE・FIFO は未確認 | example/v2_C/C2m_putget_measure.sad、v2_next/C3_trcod_orbit.sad |
| RADLIGHT (v2) | v2 | SOL 領域を含む区間 (SADError)、COORD による幾何の変化、曲げのある MULT の中の点 (出口だけ)。DX・DY・F1・FRINGE のある BEND は中の点だけ 1e-6 m 程度ずれる。高エネルギーでは点数が γ に比例して増える (3 GeV・0.26 rad の BEND で 2.5 万点)。`RadiationField` の振幅・位相・スペクトル・全エネルギーは Liénard-Wiechert の厳密な数値と合う (1 GeV の BEND、点間隔 1/γ の角度あたり 16 点以上で振幅 3e-3、64 点で 2e-4。観測点がスペクトルの波長で Fraunhofer の距離に無いと光源の大きさ ρ/γ の分だけ物理的に違う)。RAD が ON のため軌道が放射で変わる点に注意 (同じ軌道を別の計算に渡すとき)。gs との比較は無い (gs は SIGSEGV: B49) | 07_spin_solenoid_special/flag_RADLIGHT_fn_RadiationField_fn_RadiationSpectrum、v2_next/C3_radfield_bend.sad |
| DELCOR (v2) | v2 | corlength≠0 の規格化 (delx を DX の rms と読む) の根拠は確定 (K64 は複素数の ⟨\|c\|²⟩ = 1。corlength = 0 の枝は成分ごとで約束が揃っていない) が、作者の意図かは推定 (B99。K64 の値が必要なら k64 モードか delx に 1/√2 を掛ける)。CAVI の符号 (B40) は追跡に効かず観測できない (+delx のまま)。rms は多数回の平均で、GCUT を掛ければその分小さい。相関長が格子間隔 (周長/nc) より短いと線形補間で rms が下がる | |
| RADINT・Emittance[]・EMIT (v2) | v2 | RADINT は MULT・SOL が積分に入らない (表示が出る)、ANGLE≠0 の BEND の K0 は K64 のまま、ANGLE=0 の K0 の極は入るが εx は Emittance[] と 3.5%・Jx は 0.7% 違う (極の中の分散の伝搬を扇形の行列で近似し、K0 の場の端は無視。B115)。SEXT のフィードダウン・BEND の DX と K1 の勾配の軌道のずれは入れていない (B126)。軌道のずれで経路長が変わるとき RADINT (δ = 0 の約束) と Emittance[] の J が違う。結合リングで νy ≈ νs のときの J は RADINT と追跡・Emittance[] が 0.015 違う (モードが混ざり定義できない条件)。K1≠0 の BEND の傾き θ ≠ 0 の CKL1Y の数値の確認が無い。空洞は VOLT を単純に足す (PHI を見ない。ベクトル和にするかは仕様の判断)。既定の分割 (EMITDIV=1) では Emittance[] の εx が低エミッタンスのリング (大きなリング A) で 2.1% 大きく Jx が 1.5e-3 大きい (`EMITDIV=5` で εx 8e-5、`EMITDIV=20` で 6e-6・Jx 1e-4)。Jx・Je・τ・σδ は EMITDIV を上げても U0/E のオーダーの差が残る (RADINT は U0/E の 1 次の理論)。**v2 の Emittance[] と RADINT が矩形磁石で 1.4e-3 違う** (大きなリング B: v2 6.9651e-9、RADINT 6.9554e-9、k64 は 1e-4 で整合。v2 の修正の側の問題かもしれない: 未調査)。**v2 の Emittance[] (RADCOD) は gs と 1e-4〜1e-3 ずれる** (B111。k64 は全桁一致)。K1≠0 の BEND の CALC の解析行列は y ≠ 0 の軌道で自身の写像の数値微分と y = 2 mm で 1.5e-6 違う。K1=0 の BEND の追跡のスピン (tbend + tradk) が Emittance[]・tbendi (v2) と 2.2e-5 rad 違う (放射の分割が粗い K64 の精度)。RADCOD の 1 周の解析 T と写像の数値微分が 2e-3 違う (K64 の tradke の微分が近似、k64・v2 とも)。tthine の bsir0 (B110) は直していない (影響なし)。K1≠0 の BEND の四極が半径方向の面で始まる模型は K64 の仕様 (直さない)。ANGLE=0 の BEND の AE1・AE2 の残る 3〜4% (傾き 5 mrad の極。30 mrad では 2〜3%)  **v2.2.0 (P2)**: 上の「矩形磁石で 1.4e-3」は**v2 の誤りではない**: EMITDIV の既定 (1) の離散化の誤差 (+1.4e-3。扇形でも k64 でも同じ。EMITDIV = 2 で 9e-5、4 で 6e-6) で、v2 は EMITDIV = 16 の収束値が独立な平面リングの放射積分 (I1〜I5・Jx・Je を 3x3 の行列で計算) と 1.3e-6 で一致する (大きなリング B)。k64 の「1e-4 の一致」は離散化 +1.3e-3 と端の項が 2 倍の誤り (B75) −1.3e-3 の偶然の相殺。**RADINT との比較には EMITDIV ≥ 2 (できれば 4) を使う**。v2 の既定は変えていない (要判断)。RF が分散のある場所にある矩形リングの 6D の εx は独立な 6x6 の平衡 (`_c3_indep.coupled_equilibrium`) と扇形で 4e-5、矩形 (E = 0.5) で +1.7e-3 違い、原因を特定できていない (v2 と独立な 6x6 のどちらが誤りか決められない)。RADCOD の 1 周の解析 T の誤差は B138 を v2 で直して BEND 1 台 1.7e-10・1 セル 6.4e-9・24 セル 1.9e-7 に減り (k64 は 5.1e-7・1.9e-6・4.3e-6)、残りは BEND の T[z, δ] の近似 (B139。K64 の仕様)。**K1≠0 の BEND の端の模型**: v2 のオプション `compat.set_bend_k1_edge('face')` (環境変数 `PYSAD_BEND_K1_EDGE=face`) は四極の磁場も端面から始まる硬い端の磁石にする (CALC・追跡・Emittance[] の Python の経路で、ローレンツ力の厳密な写像と CALC 3e-10・追跡 1e-9)。既定は K64 の `radial` (k64 は常に `radial`)。使えない: DISFRIN=1、ANGLE=0 の K1≠0 の BEND・MULT の端、SOL 領域内、RADINT、F1・FB の端。軌道が中心 (x = 0) のリングでは 2 つの模型は完全に同じ | |
| `SynchroBetaEmittance` | 両 | **v2.2.0 で原因を特定して v2 で修正** (B128): K64 は源泉に閉軌道 c(δ) をそのまま使い (固定点が c にならず平均が半分)、x_b と δ の相関 η σδ² を入れず、リングの始点の分散 η ≠ 0 で εx が始点に依って数倍〜40 倍過大 (k64 は gs と 2e-7 で一致するので大きなリング A は 6.18e-9 のまま)。v2 は Emittance[] の 1.0003〜1.0005 (24 セル FODO の 3 通りの始点)・1.0084 (大きなリング A、AzimuthalModes 3。32 で 1.0011)・1.0045 (大きなリング B)・1.005 (PF)・1.006 (ASTRID2)・1.019 (ASTRID2 + MPW)。残る 0.3〜0.8% は物理的な効果 (色収差・β(δ)) か離散化かを区別していない。**使えない・未確認**: SB が示すシンクロ・ベータ共鳴の効果 (νx ± mνs = 整数の近く) の独立な基準 (長時間の 6 次元追跡) との比較、平均の非断熱な部分 (νs が大きい・ベータトロンチューンが整数に近い)、結合リング (SQ・SOL) で始点の分散 η_y ≠ 0 の確認は FODO + skew 四極の走行だけ。v2 の警告 (AzimuthalModes m と 1.5 倍の m が 20% 以上違うとき) は残し、v2 は大きなリング A で出さなくなった。結合のあるリング (24 セル FODO + skew 四極) では gs も落ちず k64 は 7 桁一致。`AzimuthalModes->1` は pysad が IndexError (gs は 8.6e-7 を返すが意味がない) | 04_emittance_radiation/fn_SynchroBetaEmittance、v2_next/C3_synchrobeta_coupled.sad、P2_sb_dispersion_start.sad、P2_sb_rings.py |
| 追跡 (RAD) の値 | v2 | FLUC の乱数列は K64 と同じ (`NPARA=n;` で gs の子プロセスの分け方も再現)。v2 では放射の式 (B111・B112)・端の角のある BEND (B75) が変わるので値は K64 と違う (k64 モードは gs と一致) | |
| BEAMBEAM の傾き (v2) | v2 | 傾き回転 (B34) の修正は共変性の確認のみ (t = 0.3 の 1 例。gs は座標が潰れるので比較できない)。実機の衝突点の規模では未確認。強ビームに x-y 結合があるときだけ結果が変わる | |
| `tsteer` (ANGLE=0 の BEND) | 両 | v2.2.0 で F1・FB1・FB2・FRINGE と端の放射 (B108・B140・B143) を実装・修正した (§6)。残り: 追跡と CALC の差は K0 の 2 次で残る (AE = 0.3 のステアリングで K0 = 0.0015・0.003・0.006 に対し 4e-6・1.6e-5・6e-5。K64 の端の蹴りの 2 次の違い、直さない)。ステアリングに `ROTATE` と `DX`/`DY` を同時に付けると gs の CALC のチューンが gs の追跡と違う (B144。v2 は追跡と 5e-7 で一致。k64 の CALC は gs と 1e-4〜7e-4 違い、K64 を再現していない: 未調査) | |
| `FindRoot` (v2) | v2 | 変数の範囲 `{var, ini, {min, max}}` は複素では無視。Abs・Conjugate・Re・Im を含む式は収束を保証しない。実数の初期値で複素の値になる式は `General::wrongval` | |
| 特殊関数 | 両 | 複素の PolyLog・HypergeometricU・Zeta・HurwitzZeta は実装済み (§6)。**v2.2.0 で ExpIntegralE・PolyGamma・WhittakerM/W・ProductLog・CosIntegral・GammaRegularized・Beta・Factorial・`Gamma[a,z]`・複素の次数の BesselJ など 15 関数を複素に拡張** (`specfn_cx2.py`。732 式の走査で gs と 1e-9 が 436 → 641 式。残り 91 式は gs が誤る・不正確: B158)。使えない: 複素の次数の PolyLog・ExpIntegralE・`HarmonicNumber[z,r]`・`LegendreP[ν,m,z]` (m ≠ 0)・複素の a の GammaRegularized・BesselY/I/K の複素の次数、Re a ≤ 0 の複素の a・b の HypergeometricU (SADError)、z = 0 の U、\|s\| > 500 (Zeta・HurwitzZeta)、\|z\| > 1e8 (PolyLog)・\|z\| > 30 (U・W・追加した関数。Legendre・WhittakerW は ≲ 20、複素の次数の BesselJ は \|x\| ≤ 25)。PolyLog の負の非整数 s が大きい (s < −10 程度) ときの精度は未確認。HurwitzZeta の q が非正の整数 (`HurwitzZeta[2,-1]` は gs が 2.6449、pysad は NaN) | v2_next/P5_specfn_cx_scan.sad |
| `ElementValues` (`Couple.n`)・`FREE` の成分 | 両 | フィット変数との連動 (`CouplingMatrix`・`CoupledKeys`)、COUPLE で結合した要素の成分、`USE NOEXPAND` の展開済みの成分の値の保存は未確認。**先頭でない成分の既定のキーの変数は作らない** (gs は値を持たない変数として数える: 反復の行の変数の数が 1 つ違い、その変数だけで GO すると pysad は `?No variable.`。`VISIT`・`BYE` の退避・復元は実装済み。`BYE ALL` は gs が SIGSEGV で比べられず、VISIT の未定義の LINE のメッセージ (gs は `Missing beamline in VISIT.`) は違う)  **v2.2.0: 先頭でない成分の既定のキーの変数は作るようにした** (gs と全桁一致、§6)。`COUPLE` で結合した要素の成分、`ElementValues` との連動、`BYE ALL` は未確認。反復の途中の `qcod`・`qtwiss` のオーバーフロー行の数が gs と違う (値には影響しない)。上の「作らない」は v2.1.0 までの記述 | 13_packages/fn_Couple、v2_next/B2_free_component.sad |
| `RENUMBER` | 両 | 実装済み (INDEX 属性を含む)。ただし RENUMBER のあと `LINE["POSITION","QF.1"]` が `LINE["NAME"]` の表示と食い違う gs の癖は写さない (pysad は `LINE["NAME"]` の表示と一致する位置を返す) | |
| SADScript の細部 | 両 | `LINE["POSITION",{{1,2},{3,4}}]` (入れ子のリスト。gs は入れ子のまま返し pysad は平らにする)、`LINE["POSITION",0.5]` (gs は 244、pysad は 1)。`f[x_]:=Hold[x]` の仮引数 x や `Hold[1+2]` の算術の評価 (gs は値を埋め込む / 評価する)。`LINE["LENGTH",loc]` (引数 1 個に限る K64 の仕様で、位置の取り方の案内だけ足した: B24)。`BeginPackage` の戻り値 (`{$ContextPath,"`"}`) の gs との一致は確かめていない (要確認)。除算 `a/b` は gs が `a*(1/b)` で最後の桁が 1 ulp 違うことがある (`Log[100]/Log[10]` は gs 1.9999999999999998、pysad 2。評価器全体は変えていない: 全試験の基準値に影響)。`Through[{f,g}[1,2]]` はどちらもエラーだが文言が違う。`PROFILE` の代入は gs も `General::invset`。MARK の Twiss のキー (AX など) は CALC のあとの SetElement で gs は MAIN の 0、pysad は CALC が書き戻した値。`Element["K1","QF"]=0.9` の読み戻しが 1 ulp (0.8999999999999999) | 13_packages/fn_OptimizeOptics ほか |
| 文脈 (`BeginPackage`・`Begin`・`$Context`・`WithContext`) | 両 | 利用者のコードでは実装済み (代入・`WithContext`・入れ子は gs と全行一致)。差: pysad の ``System` `` が K64 より広い (`AccuracyGoal`・`Plot` などのオプション名を含む)。文脈つきの `Names` の規則が決められない (gs は状況で `{}` と `{"A`"}` を返す)。init.n が大域に作る名前は 51 個しか調べていない (gs の `Names` が列挙しない)。`Off[f::tag]` による抑止 (pysad の `Off` は何もしない)。`Message` の `Tkinter::undef` は出さない (pysad の Tcl 橋が gs と違う経路で通るため)。SAD 自身の Packages の読み込みは従来どおり文脈を無視 | 08_scripting/fn_Get_Needs_BeginPackage、v2_next/B2_context.sad |
| 台本の読み方 | 両 | `Ffs;`・`FFs;` は gs と同じ構文エラー、LINE 名は書かれたとおり (`USE ring` は ring が無いとエラー。gs の「USE の失敗のあとの仮の 1 要素の LINE」は作らない)。`FFS` (`;` なし) の次の行が USE 以外の文のとき gs はエラーにして行の残りを捨てるが pysad は FFS 部として実行する。MAIN の `A = 3 B = 4;` は gs では A = 3 (エラーつき) だが pysad は代入しない。MAIN の要素定義の式は gs の文法で読む (v2.1.0)。gs が構文エラーにするが pysad が受ける形: `45DEG`・`SIN[0.1]`・未定義の名前 (`L=A`)。gs が受けるが pysad が構文エラーにする形: 末尾の `1+`。MAIN の変数を `E`・`Pi`・`Degree` で代入しても要素定義の式の `E` は SADScript の定数。`LINE[{"EMITY",…},1]` は gs が MAIN の値、pysad が MARK に書いた値 (0) | v2_next/D2_main_expr.sad |
| `Parallelize` | 両 | **Windows は逐次** (`multiprocessing` の spawn は Interp・FFS・ラティスを pickle できない)。POSIX は fork で並列。numba と fork の併用は実際には試せていない (WSL に numba が無く Windows に fork が無い)。(i) numba を `parallel=True` にする、(ii) 初回のコンパイルを子で行ってキャッシュを複数の子が同時に書く、(iii) OpenMP 版の BLAS の numpy、のどれかなら未確認。親でカーネルを一度走らせてから呼ぶ | v2_next/B2_parallelize.sad |
| `Packages/*.n` の関数 | 両 | `functions.n`・`Misc.n`・`element.n`・`beamline.n`・`emit.n`・`dapert.n` は小さな台本を gs と行ごとに突き合わせた (v2.1.0)。残った差: `DAPWIDTH` の窓 (§1)、`DynamicApertureSurvey` の `Axes->"XZ"` で z が 3 点のときの gs のエラー、`SymplecticInverse` の `-0` の表示、`Names`、不正な Axes のメッセージ | v2_next/B2_packages_compare.sad |
| `CorrectOrbit` の Micado・`CorrectLocal`・`PeriodicR` | 両 | gs と比較済み (PeriodicR は 12 成分が 3e-16、CorrectLocal の補正量が 1e-14、Micado は水平・2 平面同時・`SVDmode->1` も含め k64 が gs と 5e-6、v2 は B129 の修正で n = 1〜12 が単調に下がる、`MakeBump` は gs と 1e-12)。未確認: `CorrectOrbit` の `Bump->` は K64 の呼び出しの誤りで gs も pysad も動かない (B155)、`CorrectOrbitBump` の直接呼びは gs が 300 s で終わらず比べていない、`SVDmode->1` は Micado の n を無視する (B157。k64・v2 とも)、BPM の誤差 (`Weight`・`Orbit->` の乱数) つきの補正 | v2_next/D2_micado.sad、P5_README.md |
| 描画 | 両 | `CanvasSymbol`、`Canvas3D*`、`UseBG->True` の文字背景。`OpticsPlot` の text 項目数 (gs 69、pysad 74。Tk のフォントによる `RemoveOverlap` の差)。`BeamPlot` (gs は `SIGXX=0` で楕円を描かない) | 11_graphics/fn_OpticsPlot、fn_BeamPlot |
| GUI | 両 | 最初の `Window[]` の親 (gs は Toplevel、pysad はルート窓)。起動時の窓の再現が要る。Tk 9.0 (Python 3.14) では `cget` が長さを pt で返す (`-padx` 4p など) ので、SAD の `.n` が `cget` の値を数として使うと壊れる (v2.2.0 の確認でも該当する `.n` は見つからなかったが、網羅ではない。見つかったら `winfo fpixels` で直す) | 12_gui_tkinter/README |
| `Packages/CSR.n` | 両 | 全体は読み込まない (必要な部分だけ) | 06_collective/fn_CSRInit_CSRTrack |

## 3. 未確認 (コードはあるが gs と数値で比べていない)

| 項目 | モード | 状況 | 例 |
|---|---|---|---|
| v2 モードの機能全般 | v2 | 値の根拠は解析式・独立な数値計算で、gs との比較は無い (上の §1a・§2)。特に MEASURE・TRCOD (K64 のソースの読みのみ)、SOL 領域内の BEND・MULT・CAVI・ずれ・F1 (ローレンツ力の積分)、UND・PHSROT、`Polarization` の有限値 (線形の軌道 + 厳密なスピン回転の MC と FLUC つきのスピン追跡)、領域内の AUTOPHI (gs と一致するのは TRPT の場合だけ)、RADLIGHT の 3 次元の軌道、v2.1.0 の修正の値 (B111・B114〜B117・B125〜B127・B130 ほか: 独立な式・厳密な幾何で確認し、gs とは k64 モードの再現だけ比べた) | example/v2_B、example/v2_C、example/v2_next |
| 大きなリング B (8 GeV、周長 1436 m) | 両 | 数値だけ確認 (CALC・RADINT・Emittance[]・COD は k64 が gs と 1e-10、v2 は εx +0.13%)。U0/E が大きいリングの RADTAPER・tradke の符号 (B111)、陽子の大きなリング、POL (四極 60 台に 0.1 mm の誤差を入れると光学が不安定で NaN) は未確認 | |
| `FITP>1` の最終解、`tffamsetup` (fam 点) | 両 | FITP>1: 12 セルの FODO の 3 例 (FITP 5・3・1) で最終の残差が gs と 2e-4 以内 (`beta` は 9.0772e-4 で全桁一致、K1 は 6 桁一致、K2 は 1e-4 ずれる: 六極の谷が浅く反復の経路で決まる)。経路が違う原因は gs が試行点の `Overflow & closed orbit not found` で 1〜2 反復早く止まること。24 セルでは独立な最小化と 5e-8 で一致 (gs は最小値の手前で止まる)。fam 点: 初期条件の Courant-Snyder の不変量と、八極の薄いキックのヤコビアンによるチューン (1e-4) は確認。有限振幅の周期解の条件そのもの、振幅依存チューンシフトを目標にした FIT の収束、結合・六極のあるラティスの fam 点 (`MatchingAmplitude`・`NFAMP`)、大きなリング・PF・ASTRID2 の FITP>1 (オフモーメンタムの閉軌道が不安定で gs が長くかかる・落ちる可能性があり未実施) は未確認 | 03_matching/FFS_FITPOINTS |
| `CSROYMatrix`・`OYStability`・`OYScan` | 両 | `CSROYMatrix` の行列そのものは、手で作った入力 (nr=64、nj=3、nphi=2) で gs と 4e-20 で一致。gs で `CSRSetupOY` の出力から作ると wrongtype (原因は特定: `tfcsroysetup` が zj の実数リストを head から書く。B60)。**v2.2.0: `OYStability` の後半 (CSROYMatrix → Eigensystem → 成長率) を pysad の `CSRSetupOY` の出力を台本に書き込んで gs に通し、18 × 18 の行列・固有値・成長率が 1e-12 で一致** (安定 gr = −1e-4、不安定 0.027279622602。nr = 256、nj = 6、nphi = 3)。`OYStability`・`OYScan` の通し (CSRStability の構築・ZL の用意。MakeZL は 7 分) は gs が完走せず比較なし。実用の物理 (閾値の強度) の独立な基準 (Vlasov 方程式の既知の解など) は作っていない。pysad は両モードで `CSRSetupOY` の出力を受ける (k64 で gs の wrongtype を再現するかは要判断: 再現すると OYStability・OYScan が使えなくなる)。nj・nphi が大きい手作りの入力は未確認 | 06_collective/fn_CSRInit_CSRTrack、v2_next/D2_csroymatrix.sad |
| `tbrote` の `srot` 更新、`tbendebody` の `bsi2` | k64 | ソースの読みだけで、gs の個別確認なし。v2 は行列積・入口と出口の縦磁場に直した (Emittance[] と追跡で 1e-5 以下、Thomas-BMT との比較) | |
| `measure.py` の `tlum`・`tluma` | 両 | gs で MEASURE が動かない。解析値のみ。K64 の MEASURE・TRACK からは呼ばれず、`ttstat_lines(…, lum=True)` で出せる | |
| `tfloor.f` の `cfunc` (`xsincos` など) | 両 | 個別実装 | |
| `BeamPlot`・`ColumnPlot`・`LaTeX`・`PointSymbol` の見た目 | 両 | 動くことだけ確認。LaTeX は `$HaveLaTeX` が空の状態での一致 | 11_graphics |
| Tkinter の見た目・操作 | 両 | 色・大きさ・フォントは SAD と比べていない。Python 3.14 (Tk 9.0.4) では、ウィジェットの木・幾何・フォント・色を 3.12 (Tk 8.6.13) と全ウィジェットで比べ、標準出力は Tk の版の行以外同じ。Checkbutton・Radiobutton の幅が 2〜3 px 大きい。実際の画面の撮影とマウス操作は未確認  v2.2.0: Python 3.14.9 + Tk 9.0 で Window・Button・CheckButton・TextLabel・Entry・ListBox・Scale・Canvas・Menu・KBMainFrame・KBFComponentFrame が動き、窓の大きさは Tk 8.6 と同じ (白以外の画素数の差 0.1% 以内)。PNG 4 枚を `example/v2_next/P5_tk90_*.png` に撮った。gs の窓 (WSLg の X11) の撮影、高 DPI は未確認 | 12_gui_tkinter |
| Python 3.14 (Tcl/Tk 9.0) | 両 | v2.1.0 で `test_tclbridge.py`・`test_plot.py`・`test_gui_lang.py` の 25 件と、`test_v2_D2_*`・`test_spin.py`・`test_corr.py`・`test_csr.py`・`test_fam_points.py` の 90 件が通る (1 skip)。全試験は 3.14 では走らせていない。3.10・3.13 は未確認 (3.10 は同じ `nan+nanj` の分岐なので通る見込み) | |
| 乱数の消費の NaN の扱い (B17) | 両 | 決定済み: 現状のまま。`asinz(NaN)`・`gfloor`・`sqrt1` の下限は、発散した粒子が落ちるまでの座標の埋め方にだけ効く。「NaN になった粒子をその場で落とす」に変えると落ちる要素・周が K64 と変わる (LOSSMAP の値)。マッチングの fuzz (tffsmatch) は未確認 | |
| k64 モードの最初の呼び出しのスピン追跡の損失 | k64 | K64 の B33 の癖 (落とす印 0 の粒子が後で入れ替えで戻る) は track6 の経路では近似 (gs は要素 60、pysad は 15)。高速追跡 (track6n) は全桁一致 | |
| MASS・CHARGE を与えた粒子の Touschek・Intrabeam・RADINT | 両 | `rclassic = charge²·elradi` は電子の値のまま (B112)。陽子での確認なし。放射の β の因子も直していない (β < 0.99 の粒子の放射は U0 の 1/β² ほど大きく出る)  v2.2.0: 放射の β の因子 (B112) は直さないと判断した (v2 の U0 は古典式のちょうど 1/β³ 倍。陽子 1.5・3・30 GeV/c で 64%・15%・0.15%。直すには光子数とエネルギー分布の β 依存を厳密な放射のスペクトルで入れ、独立な基準と比べる必要があり、`crad` に β³ を掛ける近道は量子励起を無視するので入れなかった。k64 は gs と全桁)。重イオン (電荷 > 1) の放射の β の因子、β < 0.8 の粒子は未確認 | v2_next/A1_proton_ring.sad |
| `Twiss["DETR"]` ≥ 1 | 両 | **v2.2.0**: 90 度回転 (COORD CHI3 = ±π/2) の対で DETR = 1.375 の区間を作り、CELL CALC・INS・Emittance[] を 3 つのリング (大きなリング A・PF・ASTRID2) で gs と比べ、β・DETR・Emittance[]・U0 が一致。一周の位相の整数部は K64 では丸め誤差の符号で決まり、pysad の k64 は gs と ±1 ずれる (B153。v2 は入れ替わりの入口を +π・出口を −π と決めて gs と同じ整数部にした)。**未確認**: MARK に生の DETR ≥ 1 を書いて入れ替わった入口にする INS (gs は入口の DETR を R から計算し直す: B159。k64 を gs に合わせる修正は途中で止めた)、MEASURE の入口が DETR ≥ 1 の場合、非常に強い結合・回転以外のモード入れ替わり (ソレノイドを含む位相のゆっくりした入れ替わり) の整数部 | v2_next/A1_n12_ins_detr.sad、P5_detr_*.sad |

## 4. 既知の差 (未解決・または意図した差)

| 差 | モード | 状況 |
|---|---|---|
| SOL 内の四極 + RAD | k64 | gs が RAD ON・BZ≠0 の領域内で四極の力を入れていない疑いが強い (gs 側の欠陥。B68)。v2 の領域内の四極は物理 (Thomas-BMT・運動方程式の積分) に合わせた。現在の gs と bench_level9 の保存値が違う |
| `pm1.sad` の NaN (K3 が強くずれの大きい MULT リングの CALC の閉軌道が gs と別の解) | 両 | 原因未調査 |
| `LOSSAMPL`・RAD・FLUC つきの TRPT、MAIN の TRACK の損失報告、DAPERT での APERT | 両 | gs と比べていない |
| v2 の Emittance[] (RADCOD) と gs | v2 | 1e-4〜1e-3 ずれる。**意図した差** (tradke の符号の修正。B111。k64 は全桁一致)。大きなリング A (RADCOD) の εx +1.9e-4、ダンピング率 +2e-4、U0 +1e-7 v2.2.0: **v2 側が正しいと確定した** (放射反作用つきのローレンツ力の積分と v2 の tradke が 6 桁一致。1 台の偏向磁石・既定の EMITDIV = 1 で k64 は Δpx が 9% 違う。B111) |
| v2 の Emittance[] と RADINT (矩形磁石) | v2 | 1.4e-3 (大きなリング B)。v2 の修正の側の問題かもしれない (未調査) |
| PSPAC の gs との x と y の 16% の非対称 | v2 | 再現しない (v2.1.0。同じ記述で gs は x 1.2338・y 1.2356 mm で対称、v2 は 1.2270・1.2282 mm、解析式 1.2481 mm)。元の例の乱数のビームか PSPAC の状態の持ち越し (B96) による可能性が高い。粗いメッシュの差の原因はメッシュの分解能で、警告を足した |
| ASTRID2 の MPW・K0 の極の RADINT | v2 | εx が Emittance[] と 3.5%、Jx が 0.7% 違う (B115) |
| 追跡の FLUC の脱偏極 (v2 の Polarization の過小評価) | v2 | 解決 (v2.2.0): 原因はベータトロンのチューンの δ 依存 (クロマティシティ) がスピンの側帯の共鳴の重みを変えること。クロマティックな変調を側帯の式に入れ、4.2 GeV・12 MV で式 / 追跡 = 0.43 → 1.05 ± 0.16。ただし追跡の測定が 1 点 (§2) |
| 偏向磁石の本体の B∥ の欠落 (B136)・追跡の K1=0 の BEND の D の刻み依存 (B137) | 両 | 追跡・Emittance[]・K64 の 3 つとも、鉛直方向の傾き py ≠ 0 の粒子のスピン回転が厳密な値と (2+a) θ py 違う (0.5 GeV の扇形 15°・py = 5e-3 で R[0,1] が −1.462e-3 に対し −4.05e-3)。n0 の傾き・SpinTune・脱偏極の D に効き、v2 の脱偏極率が鉛直の軌道・分散の大きいリングで追跡と合わない成分を持つ (低エネルギーで目立つ)。**直していない (要判断)**: 直すには tradke・tradk・tbend の全ての要素の本体に B∥ を足しフレームの扱いを揃える必要があり、B58・B64・B130 のスピンの結果が全部動く |
| TrackParticles の F1 の fb1 ≠ fb2 (B145)・AE のある極の放射の 2.7% (B143)・AE の和則の破れ | 両 | F1: x が厳密な Maxwell 場と (K0/al)(fb1² − fb2²)/24 違う (K64 の模型の近似)。AE: v2 の ΔJx は厳密な幾何の 1.03〜1.04 倍で、**残る 2.7% は AE の大きさ・K0・分割に依らず極の長さ・間隔に依る** (傾き 30 mrad の 0.02〜0.03 はこの 2.7% と同じ)。追跡の ΔJx + ΔJy + ΔJe は 0 でなく +2.7e-4 (AE = 0.3)。原因は未特定 (要判断) |

## 5. gs 側の問題で比べられないもの

| 項目 | 状況 |
|---|---|
| `For` | gs は 4 番目の引数 (本体) を実行しない (B71)。k64 モードの pysad は gs と同じ (本体を実行しない。以前の pysad は実行していた)、v2 モードは本体を実行する |
| Bessel・HermiteH・Ei などの一部 | gs が誤った値を返す (`BesselJ`(100) など。B43〜B47、B53、B54)。pysad は scipy の値 (真値)。`InverseErf` だけ k64 モードが K64 の値、v2 モードが真値 (B47) |
| 複素の `PolyLog`・`HypergeometricU`・`HurwitzZeta` | gs が誤る点がある (B118〜B120: 非整数次の PolyLog の Re z < 0・Im z ≤ 0・\|z\| > 1、HypergeometricU の b が整数、HurwitzZeta の複素の s で a の畳み方によって誤る)。pysad は両モードで真値 |
| gs が落ちる | `RadiationField` の SIGSEGV (B49。ほかの引数の渡し方でも落ちる)、`TRCOD` の SIGSEGV と MEASURE の無出力・MAIN の TRACK (B42。先頭の MARK の位置だけ動く: B151)・`ON FOURIER` の TRACK の SIGSEGV (`tinip`)、`TCAVI` の L≠0 (B36)、PSPAC の既定 (`PSPACNZ=1`、B93)、`SPAC` で粒子数が少ない (2〜33 個などで SIGSEGV・NaN。**原因は粒子数が少ないことで RADIUS の大きさではない**: spkick.f の spmesh の添字。B98)、BEAMBEAM が SOL 領域内にある配置 (B35)、`RenumberElement` (B51)、`BYE ALL` (SIGSEGV)。**`SynchroBeta` の結合ラティスで gs が落ちるという以前の記述は再現しなかった** (24 セル FODO + skew 四極で gs は落ちず k64 と 7 桁一致。SOL・強い結合で落ちるかは未確認) |
| gs が誤った値を返す | SPAC の力が 1/γ 倍に弱い (B89。解析式の K/γ に近い)、領域内の BEND の偏向をほぼ無視する (B69)、SOL (BZ≠0) 領域内の四極 + RAD (B68)、WSPAC のチューンシフトが 3〜4% 強い (B127)、WSPAC ON の CALC・CALC6D が NaN (B133)、`Polarization` (B55・B125)、光子の偏光 ξ1 が ±1 を超える (B123)、`Union` の `SameTest` (B121)、`CSRSetupOY` の出力で `CSROYMatrix` が wrongtype (B60)、複素引数の特殊関数 (B158)、`MARK` の `DETR` (B159)、`CorrectOrbit` の `Bump->` (B155)、MAIN の末尾の演算子 (B156) |
| `draw` コマンド | gs は表示環境がないと出力せず、`FFS["draw ...",12]` は内部エラー。pysad は無視する |

## 6. 実装・解消して一覧から外したもの (現状)

記録 (`history/v2_changes_*.md`・`history/v2next_*.md`) と試験で確かめた。「両」はモードによらない。

| 項目 | モード | 現状 |
|---|---|---|
| `LOSSMAP` | 両 | `TrackParticles` の戻りに 2 行 (落ちた要素・落ちた周) を足す。gs の例と一致 |
| `RADTAPER` | 両 | Emittance[]・EMIT・追跡・CALC6D に入る。gs (CODCONV=1e-4) と x 6 桁、U0 10 桁、チューン全桁。v2 は追跡の閉軌道が 1 周で閉じる (B111) |
| MULT の `SK0 SK1 SK2 …` (スキュー成分) | 両 | 以前は黙って無視していた (pysad の不具合)。gs と追跡 2e-13、Twiss 1e-9 |
| `ResponseMatrix` の TRPT、`CalculateOptics` の TRPT の加速、加速区間の分散 | 両 | gs と全桁 |
| `NOKEEPEXP`、`RENUMBER` の INDEX 属性、`LINE["INDEX",…]` | 両 | gs の実測と一致。成分ごとの `FREE` で `tfvcomp` も gs と同じ |
| `BeginPackage`・`Begin` の文脈、`$Context`・`$ContextPath` への代入、`WithContext`、入れ子 | 両 | 利用者のコードで実装 (名前空間を作る)。gs の 19 行の台本と全行一致 |
| `Shared`・`SharedSize`・`OpenShared`、`Parallelize` (POSIX) | 両 | fork した子と共有 (無名の mmap)。`SharedSize[{1,2,3}]` = 32 (gs と同じ)。Parallelize は fork で並列 (Windows は逐次) |
| 複素引数の Bessel・超幾何・Airy | 両 | 実数のパラメータで実装 (1e-12〜1e-13) |
| 複素引数の `PolyLog`・`HypergeometricU`・`Zeta`・`HurwitzZeta` | 両 | v2.1.0 で実装 (`specfn_cx.py`。独立な参照と 1e-10〜1e-11。gs の誤り B118〜B120 の点でも真値) |
| PolyLog の非整数次の大きい \|z\| (z < -12・z > 60 の精度、\|z\| ≥ 150 の NaN) | 両 | `_polylog_far` で相対 1e-11 (s ≤ 7.9、\|z\| ≤ 1e12) |
| `LegendreP` の奇数 m の符号 | 両 | 以前は符号が逆 (pysad の不具合)。+.2749545416973507 |
| `EMITDIV`・`EMITDIVB`・`EMITDIVQ` | 両 | 以前は無視していた。Emittance[] に効く (MAIN の `EMITDIV = 20;` も) |
| MAIN の `NAME = 式;` の代入、`FFS;` の書き方、要素定義の式 (`ANGLE=PI/8`) | 両 | gs の動作に合わせた (代入 32 個の式を gs と 2e-14。要素定義は 37 通りの値と 19 通りの構文エラー) |
| `For` の本体 | k64 / v2 | 上の §5 |
| RADINT と EMIT の Jx の 2.4% の差 (0.9415 と 0.9651) | v2 | 原因は K64 の Emittance[]・追跡が偏向磁石の端の角で放射の行路を変えないこと (B75)。**RADINT が正しく K64 の EMIT が誤り**。v2 は Emittance[] を直し、E1=E2=0.5 で解析値 1.005751 に対し 1.00586 (`EMITDIV=20` で 1e-4)。k64 モードは 0.9651 のまま |
| `BUCKET HIGHT` (`RADINT`) と `BucketHeight` (`Emittance[]`) の 0.010045 と 0.010035 | v2 | 原因は K64 の MomentumCompaction の誤り (B77)。v2 は 0.010047 と 0.010045。k64 モードは K64 のまま |
| RADINT の量子寿命・線形クロマティシティ・BEND ごとの半径・六極・単精度の定数・CAVI | v2 | 直した (B06・B07・B12・B13・B72・B73)。k64 モードは K64 のまま |
| DELCOR の MARK・CAVI の書き込み先、波数 −1 の成分 (B40・B62) | v2 | 直した。k64 モードは K64 のまま (CAVI は HARM・PHI に入る) |
| `PolynomialFit` の平均 0 のデータ (B18)、`InverseErf` (B47) | v2 | 直した。k64 モードは K64 の値 |
| MOMENTUM を設定しない計算の警告 (B67) | v2 | 1 度だけ警告 (値は変えない) |
| **以下は v2.1.0** | | |
| BEND の `F1=0.1 FRINGE=3` + RAD の追跡の gs との差 (相対 0.7) | 両 | 解決: tfgeo の inext・iprev の再現 (B113)。gs と 1e-6 (8 桁表示の丸め)。大きなリング A・PF・ASTRID2 も一致 |
| L≠0 の CAVI + `DX`・`DY` の差 (1e-3) | 両 | 解決: CAVI の TENT・TEXIT と V1 V20 V11 V02 を追跡に足した。gs と全桁 |
| EMIT のあとに `USE` すると追跡の RF の位相が変わる差 | 両 | 解決: `trf0` を持ち越す (K64 の仕様)。gs と全桁 |
| `TrackParticles` が `MASS`・`CHARGE` を追跡に渡さない | 両 | 解決: 渡す。陽子・イオンの追跡・Emittance[]・CALC6D・RF・ソレノイドが gs と全桁一致。放射の古典半径は v2 だけ質量に比例 (B112) |
| 乱数列 (FLUC つきの追跡) の δ の rms が gs と違う (gs 4.21e-4、pysad 4.52e-4) | 両 | 解決: 原因は gs の既定 NPARA=4 (子プロセスに分ける)。NPARA=1 で pysad (k64) は gs と全桁一致、`NPARA=n;` で分け方を再現 (FODO 200 粒子 1・10・100 周、PF リング 40 粒子) |
| 光学・Emittance への空間電荷 (WSPAC) | Emittance[]・EMIT は両、CALC・CALC6D は v2 | 実装 (反復して空間電荷つきの一周の写像から)。k64 は gs と全桁一致 (K64 の最後のスライスの 8 倍のレンズ B127 を再現)、v2 は一次のラスレットの積分と 0.2〜1.8% |
| Emittance[]・EMIT・CALC6D の BEND の F1・FRINGE、ANGLE=0 の BEND (補正磁石)・MULT の F1 | 両 | 実装。k64 は gs と全桁 (FODO と大きなリング A・PF・ASTRID2・MPW 入り)。CALC のチューン = Emittance[] のチューン |
| 光子 (`tphrec`) の BEND・QUAD・MULT・SEXT・OCT・磁場つきドリフト・BEAMBEAM + BEND | 両 | 実装。k64 は gs と光子 358・325・71 個が全列一致 (E 列は 1.4e-8 違うのは B132 の仕様) |
| tdgeo の領域内の四極の DX・DY の解析微分と幾何の fit | 両 (領域内の四極のずれ) / v2 (解析微分) | 実装 (幾何が gs と全桁、fit が 2 回の Newton で収束) |
| CALC6D・RADCOD つきの CALC・Touschek (INTRA) と SOL 領域 | 両 | 実装。k64 は gs と 1e-9 |
| 傾いてずれたソレノイドの Emittance[] の POL の量 | v2 | 修正 (B130)。SpinTune が追跡と 1e-13〜5e-7 |
| POL + TRPT + 加速 (SOL 領域の外の CAVI) | v2 | 実装 (断熱加速の解析式と 5e-7) |
| ステアリング・MULT の K0・端の角 (B108) | v2 | 修正 (厳密な幾何と 3〜4%) |
| NORFSW; RAD の Emittance[] で閉軌道がずれたときの gs との差 (1.8〜4%) | 両 | 解決: pysad の不具合 (z を δ で閉じていた)。大きなリング B (QF5 に DX = 5e-5) で gs 7.1529048e-9、k64 7.1529045e-9 |
| 不安定なリングの POL の Emittance[] が終わらない (大きなリング B で 20 分以上) | 両 | 解決: NaN で打ち切る (2 秒、`depol_error` つき) |
| `emit6d._tquade` の解析ヤコビアンが δ の微分を落とす | 両 | **再現しなかった** (Python 版・numba 版・数値微分が 1e-12 で一致)。z の行の関数単体の欠落を補った |
| 複数セグメント・複数周の `RadiationField` の失敗 | 両 | 解決 (重なった点を除く) |
| `Twiss["DETR"]` が INS の CALC の後で 0 | 両 | 解決 (gs と全桁。DETR ≥ 1 は未確認: §3) |
| `Parallelize`・複素の特殊関数・`$Context` の代入・`Symbol`・`Message`・VISIT/BYE・成分ごとの FREE | 両 | 実装・gs に合わせた (上) |
| `Optics.n` の `PeriodicR`・`CorrectLocal`・Micado | 両 | gs と比較 (Micado の誤りは B129 で v2 が修正) |
| `dapert.n`・`element.n`・`emit.n`・`beamline.n`・`functions.n`・`Misc.n` の gs との比較 | 両 | 小さな台本で行ごとに比べた。差の一部は直した (`Log2$`・`FunctionDefinition`・`Fill$String`・`Remove$Paren`・`Emittance->False`・TURNS の既定・`Element[..]+=` のフック) |
| Python 3.14 (Tk 9.0) の試験の失敗 2 件 | 両 | 解決 (`test_tk_version_and_font_system`、不安定なリングの `_spintrm`) |
| TrackParticles の F1・FB1・FB2・FRINGE (`tsteer` (ANGLE=0 の BEND)・`tmulti` (MULT: F1・F2・F1K1F/B)) | 両 | v2.2.0 で実装 (K64 の動作の追加)。k64 は gs と 10 桁表示の全桁一致 (12 セルの FODO・ステアリング 6 通り・MULT 9 通り、大きなリング A の補正磁石 2 台)。numba の高速追跡と Python は 1e-13。独立な基準 (軸上の偏向 1e-14、放射の積分 1e-5、Maxwell 場の端の運動 1e-3〜1e-10) |
| MULT の AE1・AE2 (CALC・Emittance[]・追跡)・K1≠0・K0≠0 の ANGLE=0 BEND の AE (CALC・Emittance[]) | 両 | v2.2.0 で実装 (これまで黙って無視)。k64 は gs と全桁一致 (5 通り + 大きなリング A)。v2 は B140・B141・B142・B143 を直す。大きなリング A で追跡と CALC のチューンが 2e-7 以内 |
| TrackParticles の F1・AE の放射の行路 (B108 の MULT・追跡への拡張) | v2 | v2.2.0 で修正 (厳密な幾何 + ロビンソンの定理と ΔJx が 1.03〜1.04 倍。K64 は約 2 倍) |
| 要素ごとの FREE の何もしない変数 | 両 | v2.2.0 で gs と全桁一致 (`FREE QF.3 K1` は結合先 QF.1 を持つ何もしない変数として VAR に出る、変数の数・競合のエラー・FIX・Saved の欄・RESET) |
| `DynamicApertureSurvey` の DAPWIDTH の窓 (既定 7) | 両 | v2.2.0 で `trackd.f` のスロット管理を直訳。NORAD の地図・Score が gs と一致 (24 セル FODO 4 通り、大きなリング A 2 通り)。k64 は RAD + FLUC の最終 seed まで一致。既定を 51 から 7 に変えた (51 を前提にした台本は `DAPWIDTH = 51;`) |
| SynchroBetaEmittance の大きなリング A での 40 倍 (B128) | v2 | v2.2.0 で原因を特定して修正 (上の §2) |
| 追跡のスピンが BEND の ROTATE で回らない (B135) | v2 | v2.2.0 で修正 (厳密な Thomas-BMT の積分と 1e-9) |
| tradke の微分 (B138)・tradke の符号 (B111)・RADTAPER の重み | v2 | v2.2.0 で確定・修正 (B111 の符号は放射反作用つきのローレンツ力の積分と 6 桁。RADTAPER の BEND の係数は (δ_l + δ_l+1)/2 に統一) |
| n が傾く・非平面のリングの緩和率 (B116) | v2 | v2.2.0 で厳密式 Σ r (1 − 2/9 (n·v)²) と運動論的偏極を入れた (独立な DK の数値積分と τ が 1e-5〜2e-3) |
| 領域内の QUAD・MULT のスピン (bz ≥ 3)・領域内の光子 | v2 | v2.2.0 で QUAD・MULT のスピンは 9e-7 以下 (B146)、DRIFT・QUAD・MULT の光子を実装 (上の §1・§2) |
| PSPAC の K64 の形状関数 (k64 モード) | k64 | v2.2.0 で実装 (gs と rms が 12 桁以上。NPARA = 1) |
| WSPAC の追跡のビーム行列 | v2 | v2.2.0 でオプション `WSPACSC = 1;` を実装 (既定はゼロ電流の行列) |
| MEASURE・TRCOD の k64 (先頭の位置)・DETR ≥ 1 の INS・複素の特殊関数 15 関数・`OYStability` の後半・Micado の変種・FITP>1・Tk 9.0 の見た目・文法の細部 (`SIN[0.1]`・`45DEG`・MAIN の `E`・`Pi`) | 両 | v2.2.0 で gs と比較・実装 (上の §2・§3 と docs/47 v2.2.0 の (b)) |

## 7. 優先度 ([16_fortran_status.md](16_fortran_status.md) §9.4 から)

| 優先 | 項目 | 理由 |
|---|---|---|
| 中 | 誤差入り 6D emit と DA の反復 (性能) | 反復が 600 s を超える |
| 中 | `Polarization` の追跡との差の確認 (脱偏極率) | 追跡の測定は 4.2 GeV・12 MV の 1 点。6 MV・24 MV・PF・大きなリング A の FLUC の追跡 (1 条件数時間)。偏向磁石の本体の B∥ (B136) |
| 中 | v2 の EMITDIV の既定 (1)・RF が分散のある場所にある矩形リングの εx の +1.7e-3 | EMITDIV の既定を 2〜4 にするか。+1.7e-3 は端の角のある BEND の 6D の厳密な線形化 (ローレンツ力の数値積分の写像の 6x6 の数値微分 + 放射) で v2 と独立な 6x6 のどちらが誤りかを決める |
| 中 | `SynchroBetaEmittance` のシンクロ・ベータ共鳴の効果 | 長時間の 6 次元追跡 (FLUC) との比較 (PF・ASTRID2 で 0.5% 以内。CPU が必要) |
| 中 | 逆向きの要素の F1・AE の入れ替え、F1 つきの MULT の追跡の高速化、AE のある極の放射の 2.7% (B143)、SOL 領域内の加速 (POL + TRPT)・BEND の光子、DISFRIN=0 の ψ の独立な基準 | §2 |
| 中 | v2 の要判断 (docs/47 v2.2.0 の (d)) | 複素の特殊関数の k64 の値、DAPWIDTH 既定 7・WSPACSC・EMITDIV の既定・β の因子・`face` の既定・ψ の DISFRIN・B∥ ほか |
| 低 | `wfres` 系 (MAIN の FOURIER) | 必要になったら。NAFF で代用 |
| 不要 | MON/MCOD/TCOD の機能、Obsolete・非リンクのファイル、`tdrwdt.f` の粒子追跡 | K64 で無効 |

## 8. 注意

- 条件: 例の数値は、サンプルのラティス (`example/00_lattice/FODO_ring.sad`: 3 GeV、周長 105.6 m) と、PF・ASTRID2 のラティス (`example/pf`、`example/astrid`) で得たものです。v2 の機能の数値は、各記録の条件 (上の冒頭) の範囲です。
- 表示や数字の桁の微妙な差は、ここには入れていません。完走することと、意味のある値が合うことを基準にしています。
- 一覧は変わります。項目を直したら、この文書と [16_fortran_status.md](16_fortran_status.md) §9 を一緒に更新してください。
