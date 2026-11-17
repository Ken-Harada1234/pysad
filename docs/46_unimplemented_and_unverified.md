# 46. オリジナル SAD (K64) に対する pysad の未対応・未確認の一覧

最終更新: 2026-11-17 (現状のスナップショット)

この文書は、オリジナルの SAD (K64 版、`src/*.f`・`Packages/*.n`・SADHelp) に対して、pysad が**まだ対応していないもの**と**対応したが gs で確かめていないもの**を、1 か所にまとめたものです。

- **互換モード (v2.0.0 から)**: pysad は 2 つのモードを持つ。`k64` は K64 (gs) と同じ動作 (v1 と同じ。gs との数値比較はこのモード)、`v2` (既定) は K64 のバグを直し機能を足した動作。環境変数 `PYSAD_COMPAT=k64|v2` か `pysad.compat.set_mode()` で切り替える ([47_v2_changes.md](47_v2_changes.md))。**各項目の「モード」の欄が `k64` ならそのモードだけの話、`v2` なら v2 モードだけの話、`両` ならモードによらない。**
- v2 の値の根拠は、解析式・独立な数値計算 (ローレンツ力・Thomas-BMT の数値積分など)・保存則・別のコードとの一致で、gs ではない。**v2 で実装・修正した機能の多くは、gs が動かない・K64 が誤っているので gs と比較できない。** 「未確認」にはその旨を書いた。使った条件は電子、0.5〜3 GeV のリングと輸送路 (空間電荷は 5〜50 MeV、RADLIGHT は 10.2 MeV)、ソレノイド bz = 0.3〜3 1/m で、**実機の大きなリング・挿入光源・強い結合では確かめていない**。
- 出所: [16_fortran_status.md](16_fortran_status.md) §9、[history/work_unimpl.md](history/work_unimpl.md)、[history/work_example_fixes_20261117.md](history/work_example_fixes_20261117.md)、[history/work_apert.md](history/work_apert.md)、[history/work_ffs_gs_match.md](history/work_ffs_gs_match.md)、[21_sad_bug_report.md](21_sad_bug_report.md)、v2 の作業の記録 `history/v2_changes_A.md`・`_B`・`_C`・`_C_measure`・`_C_spch`・`_C2`・`_C2_measure`・`_C2_misc`・`_C2_spch`、[47_v2_changes.md](47_v2_changes.md) の (e) 要判断・未解決、`example/` の各例の「gs と pysad の比較」・「未確認・限界」節、[34_sad_packages.md](34_sad_packages.md) ほか各解説文書。
- 2026-11-17 に、v2 の記録と試験 (`pysad/tests/test_v2_*.py`) で、旧版の一覧の項目 (MEASURE・TRCOD、LOSSMAP、空間電荷、RADTAPER、SOL 領域内の BEND・MULT・CAVI、`System`・`Fork`、`BeginPackage` の文脈、UND の磁場、PHSROT、`Polarization`、DELCOR、APERT の SOL 領域、複素 `FindRoot`、複素引数、RADINT と EMIT の Jx の差、`BUCKET HIGHT` の差など) を確かめ直して反映した。EMIT のあとに `USE` する差など、記録で解消を確かめられなかったものは未解決のまま残した。ほかは出所の記載のままです。
- 基準は gs (SAD K64) の実動作です。gs が落ちる・誤った値を返すものは「gs 側の問題」(§5) に分けています。K64 のバグの番号は [21_sad_bug_report.md](21_sad_bug_report.md) の B01〜B110 (B72 以降は v2 の作業で見つけたもの)。
- `例` の欄は、その項目を示す `example/` の例です (パスは `example/` からの相対)。v2 の機能の例は `example/v2_B/`・`example/v2_C/`。

凡例: 「未翻訳」は pysad に物理が入っていない。「一部未対応」は一部の入力・経路だけ未対応。「未確認」はコードはあるが gs と数値で比べていない。「モード」は上のとおり。

## 1. 未翻訳

| 項目 | モード | 理由・状況 | 例 |
|---|---|---|---|
| `wjfit`、`wfres`、`wfbin`、`wiord`、`wjdraw` | 両 | MAIN の TRACK の Fourier 解析で、FFS から使えない | |
| 光学 (CALC)・Emittance への空間電荷 (K64 の `twspace`・`qwspac`、Emittance の反復) | 両 | v2 の WSPAC は追跡だけ。WSPAC のビーム行列は**空間電荷なし (ゼロ電流) の平衡値**で、空間電荷によるビームの大きさの変化は入らず、Twiss のチューンにもシフトが出ない (追跡でだけ見える)。入れるにはビーム行列を自己無撞着に解く `emit6d`・`calc6d` の変更が要る (要判断: docs/47 (e)) | 06_collective/flag_SPAC |
| `tphrec` (光子) | 両 | `BEAMBEAM` の分だけ対応。ほかの要素の光子は未対応 | |
| `tdgeo.f` の四極などの DX・DY の解析微分 (`qdtrans`・`tftmatu`) | 両 | 数値微分で代用。v2 で領域内の変数の幾何 fit は中心差分 | |
| CALC (4 次元) の RADTAPER | 両 | K64 の `qtwiss` にテーパーの式が無い (CALC6D は v2 で RADTAPER が入る) | 04_emittance_radiation/flag_RADTAPER |
| CALC6D・RADCOD つきの CALC と SOL、Touschek (INTRA) と SOL | 両 | SOL があると SADError | |
| MON・MCOD・TCOD の機能 | 両 | K64 にコマンドが無く、古い SAD の MON の書式・意味を確かめられる資料が無い。v2 は「K64 に無いコマンド」と警告して文全体を読み飛ばす (定義はしない)。軌道補正は CorrectOrbit パッケージ ([34_sad_packages.md](34_sad_packages.md)) | |
| `DRAW` の中の粒子追跡 (`tdrwdt.f`) | 両 | K64 でも動かない死んだコード (`trsize` が常に偽。B105) | |
| `Out[n]:=` の gs の表示 | 両 | gs は未定義の記号 1 つの文で前の `Out[n]` を繰り返すなど値と関係のない行を出す。K64 の表示の癖として写さない | |
| `Parallelize` の本当の並列 | 両 | 逐次のまま (`{f[1],…,f[n]}` を返す)。POSIX では Fork で実装できるが、numba のスレッドと fork の併用を確かめていない | |
| 複素引数の `PolyLog`・`HypergeometricU`・`Zeta` | 両 | scipy に複素の実装が無い (hyperu)・分岐の向きを確かめる独立な参照が無い。SADError・TypeError。複素引数の Bessel・超幾何 (実数のパラメータ)・Airy は実装した | 10_math_data/fn_Gamma_Erf |
| `$Context = …`・`$ContextPath = …` への代入、`WithContext`、入れ子の文脈 (``A`B``) | v2 | 文脈 (`BeginPackage`・`Begin`) の実装 (両モード) の範囲外。文字列を `ToExpression` するときの文脈は「その時点の現在の文脈」 | 08_scripting/fn_Get_Needs_BeginPackage |
| `TRCOD` の FIFO での待ち合わせ (`PUTMEASURE`・`GETMEASURE`)・`ORBITJITTER` のパイプ | v2 | 普通のファイルで同じ形式を読み書きする (待たない)。Linux・WSL の FIFO は未確認 (試験は普通のファイルだけ) | |

## 1a. モードで状態が違うもの (k64 モードでは未対応・v2 モードで実装済み)

v2 モードの実装は、特記がなければ**gs と比べていない** (gs が動かない・K64 が誤っているため)。根拠は解析式・数値積分。

| 項目 | k64 モード | v2 モード | 例 |
|---|---|---|---|
| MEASURE・TRCOD (FIXCOD・MVCOD・STNCOD・ALIGN)・LOGMEASURE・ORBITJITTER・PUTMEASURE・GETMEASURE、MAIN の `TRACK TURNS=-1` (`trackb.f`・`tracke.f`・`trcoda.f`・`ttinit.f`・`ttstat.f`) | 未対応 (警告して無視)。gs で TRCOD が SIGSEGV、MEASURE は出力なし (B42) | 実装 (K64 のソースの読み。**gs との比較は 1 件も無い**。解析式と pysad の独立な経路で確認: docs/03) | example/v2_C/C2m_trcod_reference.sad |
| SOL 領域内の BEND・MULT・CAVI・APERT (`tsole1`・`tmulte1`)、SOL の `DX DY DZ CHI1..3`・`F1`、領域の幾何 | 未対応 (NotImplementedError)。BEND は K64 の tsol.f どおりに実装したが gs が偏向をほぼ無視して合わず取り下げ (B69) | 実装 (ローレンツ力の数値積分との一致。docs/44 §7a) | example/v2_C/C3_sol_region.sad、C2a_sol_tilt.sad |
| 領域内の CAVI・RF つき MULT の TRPT の加速・wake・`AUTOPHI` | 未対応 | 実装 (加速は運動方程式の積分と一致。AUTOPHI の意味は K64 のソースの読みで、gs と一致するのは TRPT の場合だけ) | example/v2_C/C2b_sol_rf.sad |
| 空間電荷 (SPAC・PSPAC・WSPAC・SELFCOD) の追跡 | 追跡に入らない (K64 の SPAC は力が 1/γ 倍に弱い: B89 など。ビット単位で再現する翻訳はしていない) | 実装 (SPAC は解析式・Laslett の式、PSPAC は gs と 3〜4 桁、WSPAC は gs と 1%) | example/v2_C/C2s_mult_radius.sad、C2s_wspac.sad、05_tracking/flag_SPAC |
| `PHSROT` | 追跡は NotImplementedError、光学は恒等 (K64 で動作しない: B38) | 実装 (意図された写像。固有値・シンプレクティック性・平衡エミッタンスで確認) | example/v2_C/C2_phsrot.sad |
| UND の磁場 (`FBX`・`FBY`) | 磁場なしのドリフトのみ (K64 が設定できない: B37)。光学は恒等 | 実装 (解析式で確認)。磁場がなくても光学にドリフトとして入る (要判断: K64 の「光学では恒等」を仕様と見るなら磁場つきのときだけにする選択もある) | example/v2_C/C1_undulator.sad |
| `Polarization` (平衡偏極度) の有限値 | K64 が NaN (B55)。NaN のまま写す | 有限値 (平面リング 92.376%、誤差のあるリングで追跡の脱偏極率と 0.5〜1.6 倍。docs/43 §8a) | 07_spin_solenoid_special/flag_POL_RADPOL |
| スピン追跡 (POL + RAD) の損失の処理、領域内の BEND・MULT・CAVI・ずれ・F1、UND、PHSROT | 損失の処理なし。領域内は DRIFT・QUAD だけ。UND・PHSROT は NotImplementedError | 損失を処理 (LOSSMAP の 2 行も付く)、領域内・UND を通す (PHSROT はスピンを動かさない) | example/v2_C/C2c_spin.sad |
| RADLIGHT の軌道の横方向・BEND の中 | 横方向の座標が入らず BEND の中を記録しない (K64 の癖。gs と差 0) | 3 次元の軌道、BEND の中も記録 (B65) | example/v2_C/C5_radlight.sad |
| 旧式の `WAKE … ENDWAKE` の表 | `WAKECOMMAND[]` を評価して失敗 (K64 と同じ) | `WakeFunction` の表として登録 (CAVI だけ。行をまたいで書ける) | example/v2_C/C2w_wake_table.sad |
| `System`・`Fork` | `System` は実行せず 0、`Fork` は pid を返すだけ | `System` は許可 (`PYSAD_ALLOW_SYSTEM=1`) があるときだけ実行し終了コードを返す (許可なしは `$Failed`、`[ -d path ]` などの形だけ判定)。`Fork` は POSIX で本当に分ける (Windows はエラー) | example/v2_B/B3_shared_system.sad、09_io_system/fn_System_TemporaryName、fn_Fork_Shared |
| 複素初期値の `FindRoot` | TypeError (gs は壊れた値: B61) | 複素の Newton 法 (変数の範囲は無視、正則でない式は収束を保証しない) | 10_math_data/fn_Gamma_Erf |
| `CSRSetupOY` の zj・dpj | K64 の取り違え (zj に dp が入り dpj は 0: B59) | 位置・運動量が入る (`CSROYMatrix` の行列そのものは gs が wrongtype を返すので未照合: B60) | 06_collective/fn_CSRInit_CSRTrack |
| `MessageName[f,"usage"]` の関数形 | pysad は以前から f を見る (gs は f を無視する: B70) | f の本文を返す (無ければ `General::tag`) | 13_packages/fn_OptimizeOptics ほか |

## 2. 一部未対応

| 機能 | モード | 未対応の部分・使えない範囲 | 例 |
|---|---|---|---|
| MULT `ANGLE`≠0 | 両 | Twiss の `BZ`・`RF`・`F1`・`SK0`・`FB1/FB2` (追跡は `FB1/FB2` まで対応) | 07_spin_solenoid_special/elem_MULT_ANGLE |
| MAP | 両 | `GEO`・`CompiledMap`。MAIN の TRACK・DA・Touschek は ExternalMap を呼ばない | 06_collective/fn_ExternalMap |
| POL (SOL 領域内) | k64 | 領域内の MULT・BEND・CAVI のスピン (Emittance[]・追跡とも)。スピン追跡は損失の処理なし | 07_spin_solenoid_special/flag_POL_RADPOL |
| POL (SOL 領域内) | v2 | **Emittance[]・EMIT の POL の量 (SpinTune・Polarization など) は領域の中だけ K64 の式 (B81 の bsi) のまま**: 軸上のソレノイドでは追跡と 2e-5 で合うが、ずれて傾いたソレノイドでは SpinTune で 2.3e-4、n0 の向きで 4e-3 rad 違う (`spin.tradke_spin` の `xpzb` と 6 次元の微分 dxpzb を直す必要があり未着手)。領域内の APERT のスピン追跡 (NotImplementedError)、領域内の MULT・CAVI の POL (Emittance[]) は未確認。track6 の経路 (POL + RAD) での領域内の TRPT の加速 (RF は入るが基準運動量は変えない) は未対応。MAP のスピンは未対応 | |
| `Polarization` | v2 | 軌道の振幅について 1 次 (SLIM 相当)。追跡との一致は 3 例で 0.5〜1.6 倍 (16 粒子、統計 ±25%、確かめたのは小さな 1 つのリング)。高次のシンクロトロン側帯共鳴、運動論的偏極 (b·γ∂n/∂γ の項) なし。エネルギー幅の大きいリング (aγσδ/νs が 1 程度以上) では側帯の寄与が増えるはず。K64 の `spdepol` (2 次の摂動) は大きさが追跡と 3〜5 桁合わないので置き換え (直すには導出が要る)。`Polarization2/4/6` は実装していない。POL のときは四極・MULT を細かく分割するので Emittance[] が遅くなる | |
| TCAVI | 両 | `RANKICK`・`RANPHASE`。`Emittance[]` は L=0 のみ gs と一致 (L≠0 は gs が壊れる: B36。pysad は有限値を出す) | |
| UND | k64 | 磁場なしのドリフトのみ (K64 が磁場を設定できない) | |
| UND | v2 | 端の磁極の整形なし (軌道は x' = −(F0/kz) sin(kz s) で中心が −F0/kz² ずれる)、`BX`・`BY` という名前は K64 にも無い、PHI≠0 の入口はベクトルポテンシャルの飛びが蹴りになる (K64 の式のまま)。`SLICE` は `POLE` の 20 倍以上が条件。挿入光源のあるリングでは未確認 | example/v2_C/C1_undulator.sad |
| PHSROT | v2 | スピンを動かさない (リング全体を PHSROT 1 個で表すとスピンチューンは 0 になり偏極の計算に使えない)。拡散行列に負の固有値があると ValueError。分散つきの Twiss を CALC で読むと BX・EX は入力と 1e-3 ほど違う (H 行列の定義の差)。(z, δ) の回転は向きが逆 (縦のチューンは −PZ と出る) | example/v2_C/C2_phsrot.sad |
| `APERT` | k64 | SOL 領域内は NotImplementedError | 05_tracking/elem_APERT_LOSSAMPL |
| `APERT` | 両 | 再入射 (`JDPX`・`JDPY`) は粒子 12 個以上で gs と乱数列が違う。`LOSSAMPL`・RAD・FLUC つきの TRPT、MAIN の TRACK の損失報告、DAPERT での APERT は gs と未照合 | 05_tracking/elem_APERT_LOSSAMPL |
| LOSSMAP | 両 | MAP 要素 (ExternalMap) が落とした粒子の要素番号は 0。空間電荷が落とした粒子の要素番号も。POL なしで高速追跡に乗らない要素の並びは `LOSSMAP: … 記録できない` の SADError。落ちた要素は判定 (tapert) のときの l_track なので 15 要素ごとの判定点か APERT の番号 | 05_tracking/flag_LOSSMAP |
| RADTAPER | 両 | CALC (4 次元) のテーパー (上の §1)、SOL 領域内の要素、RFSW OFF (見つからないときは `???-RADTAPER-closed orbit with taper not found` を出してテーパーなしで続ける)。追跡を Emittance[] の閉軌道から始めても 1 周で x が 4e-7 m、z が 2e-6 m 動く (K64 の追跡と tturne の放射の入れ方の違い。k64 の追跡は gs と 1e-12 で合う。直していない)。gs と比べるときは gs 側で `CODCONV = 1e-4;` が要る | 04_emittance_radiation/flag_RADTAPER |
| RF のある MULT (`tmultiacc`) | 両 | `F1`・`FRINGE`。`AUTOPHI` は k64 で NotImplementedError、v2 で実装。領域内の CAVI・MULT: `V1 V20 V11 V02`・ANGLE≠0 の MULT・領域の外と bz = 0 の領域の RF つき MULT の Emittance[] は未対応 (v2 の C2 の記録に「bz = 0 の領域の中の CAVI・MULT は領域の外と同じ写像」という記述もあり、**記録内の表現に食い違いがある: 要確認**) | |
| SOL 領域 (v2) | v2 | 使えない: SEXT・COORD・BEAMBEAM など (`SolRegionError`)、領域内の BEND の `K1`・`DX`・`DY`、追跡を領域の途中から始める・途中で終える (`ソレノイドの終端がない` のエラー)、領域より上流の幾何に軸まわりの回転 (GCHI3≠0) があるラティスで追跡を先頭以外から始めた場合 (CHI3 が入らない)。未確認: ソレノイドの軸が設計軌道と 10 mrad 以上ずれる場合、領域内の SOL (BOUND=0) にずれ・回転を書いた場合、GEO のない領域に手で書いた CHI の符号の約束、bz ≥ 3 1/m を超える場、領域内の空洞の wake の横の蹴り (gs と未比較)、加速 + RAD の組み合わせ、領域内の RF つき MULT の RAD。放射は bz·L/ndiv ≤ 0.02 に刻んで N5 (B80) の誤差を避けている | |
| SOL (両モード) | 両 | BOUND=0 で BZ≠0 の SOL (領域の外のソレノイド) の Twiss は NotImplementedError。QUAD の FRMD の向きの入れ替えは未実装。行列は数値微分で解析行列の gs と 1e-10 程度の差。RADCOD の閉軌道探索中の分割の固定は未確認 | 07_spin_solenoid_special/elem_SOL |
| 空間電荷 (v2) | v2 | **使えない**: `ANGLE≠0`・`VOLT≠0` の MULT と SOL 領域内の MULT の SPAC (メッセージを出して空間電荷なしで追う)、QUAD・BEND・SEXT の中の SPAC (偏向磁石のあるリングでは SPAC は全周には入らない。全周なら PSPAC + SPCH か WSPAC)、複数バンチ (`NBUNCH`)、壁の抵抗、円でないダクト、MAIN の `TRACK`・`DynamicApertureSurvey`、スピン (`POL`) の追跡。`NORFSW` はエラー (`RADIUS<0` の絞りだけなら走る)、粒子の z が全部同じでもエラー。**未確認**: ウェーク・放射との組み合わせ (走るが独立な根拠なし)、放射 + 空間電荷の MULT、SOL 領域内の四極、バンチの静止系で相対論的な運動 (γ·x' ≳ 1)。**`TrackParticles` が `MASS`・`CHARGE` を追跡に渡さないので電子以外は SAD スクリプトから使えない** (式は 100 MeV の陽子で 0.5% を確認。要判断)。精度: 4096 粒子の既定のメッシュで中心の勾配 (チューンシフト) が 1〜3% 小さめ。条件: 電子 5〜50 MeV、1e8〜1e12 個/バンチ、ビーム半径 0.1〜7 mm、ダクト 3 mm〜3 cm。gs との比較は SPAC の包絡・PSPAC・WSPAC の 4 例だけ。PSPAC の直線の例で gs は x と y が 16% 違う (原因は調べていない) | 05_tracking/flag_SPAC |
| WSPAC | v2 | ビームの大きさは追跡の間じゅう固定 (強-弱。エミッタンスの増大・コヒーレントな振動は扱えない)。垂直のエミッタンスが 0 の平らなリングは Emittance の経路で使えない (`WSPACBEAM` を与える)。`WSPACBEAM` は σz 一定・x-z の相関なし・MINCOUP なし。TRPT の輸送路ではエミッタンスの断熱減衰を入れていない。縦の力は既定で「軸上で 0」(K64 と同じ)。SOL 領域内・`SELFCOD` + 閉軌道のあるリング・放射つきの追跡は走るが値の確認なし。gs で `EMIT` を `NOWSPAC` で行ってから `WSPAC` で追跡すると NaN (K64 は「WSPAC つきの EMIT」が前提) | |
| MEASURE・TRCOD (v2) | v2 | gs との比較が無い。`TRCOD FIXCOD`・`MVCOD`・`STNCOD` の基準は REFERENCE と同じ参照光学に置いた (K64 の `twiss(*,ndim,*)` は利用者が読む手段が無い。`measure._ref_rows` の 1 か所で替えられる)。`ALIGN` の対象に DECA DODECA MULT を足した (K64 は BEND QUAD SEXT OCT だけ。`measure._ALIGN_V2` を `_ALIGN_K64` にすれば K64 のとおり)。SOL・MONI は動かさない。区間に BEND を含む ALIGN は平行移動との一致が 4e-10 m まで。GETMEASURE の対象の関数 (DX・DY) は K64 の関数番号 32・34 (今の番号では PEPY・PZPX) から確定できず式から判断 (B103)。**`Twiss["DETR"]` が CALC の後で常に 0 になる pysad の不具合 (両モード、未修正。docs/47 §(e)) のため、モードが入れ替わる強い結合のリングの MEASURE は未確認** (輸送路で先頭 MARK に DETR を書く場合は正しい)。TWAKE / LWAKE ON は z の分布を作る (K64 は全粒子 z = 0。一様・階段状の wake の重心と広がりだけ確認)。MEASURE の TRPT は追跡の間だけ ON にして戻す (K64 は戻さない。要判断で決定済み)。MAIN の `TRACK TURNS=-1` は粒子種の語 (PROTON など) を見ない。TRGAUSS OFF の分布と表示の桁 (autofg) は gs の実出力と比べていない。TRCOD の途中の統計 (K64 は標準エラーに出す) は出さない。**QUIET の重心合わせ (B87) と TWAKE/LWAKE の z = 0 (B88) は K64 の誤りと判断して直した推定** | example/v2_C/C2m_putget_measure.sad |
| RADLIGHT (v2) | v2 | SOL 領域を含む区間 (SADError)、COORD による幾何の変化、曲げのある MULT の中の点 (出口だけ)。DX・DY・F1・FRINGE のある BEND は中の点だけ 1e-6 m 程度ずれる。高エネルギーでは点数が γ に比例して増える (3 GeV・0.26 rad の BEND で 2.5 万点)。複数周の軌道で区間の始めの点が重なり `RadiationField` が失敗する限界 (k64 モードでは確認済み) が v2 でも残るかは未確認。gs との比較は無い | 07_spin_solenoid_special/flag_RADLIGHT_fn_RadiationField_fn_RadiationSpectrum |
| DELCOR (v2) | v2 | corlength≠0 の規格化 (delx を DX の rms と読む) は作者の文書による裏づけが無い推定 (B99。K64 の値が必要なら k64 モードか delx に 1/√2 を掛ける)。CAVI の符号 (+delx) は決定済みだが K64 が SOL と同じ形の式を写した意図があるなら逆。rms は多数回の平均で、GCUT を掛ければその分小さい。相関長が格子間隔 (周長/nc) より短いと線形補間で rms が下がる | |
| RADINT・Emittance[]・EMIT (v2) | v2 | RADINT は MULT・SOL・BEND の K0 が積分に入らない (MULT・SOL は表示が出る)、結合のあるリングの V 面 (I4V) は未確認 (B107)、空洞は VOLT を単純に足す (PHI を見ない。ベクトル和にするかは仕様の判断)。**ANGLE=0 の BEND (補正磁石) と K0 つき MULT の端の角は K64 のまま** (B108。Jx の端の項が 2 倍になる作りと見られるが独立な値が無い)。既定の分割 (EMITDIV=1) では Emittance[] の Jx が解析値より 1.5e-3 大きい (`EMITDIV=20;` で 1e-4。Jx+Jy+Je は 3.9990 で 4 との差 2.5e-4)。CALC (tbende) と追跡 (tbend) の、ずれた BEND (矩形磁石を DX=2 mm ずらすと CALC のチューンが 2.9e-5 変わる。追跡が正しい。tbende の 2 次の幾何の項で場所は未特定)。K1=0 の BEND の追跡のスピン (tbend + tradk) が Emittance[]・tbendi (v2) と 2.2e-5 rad 違う (放射の分割が粗い K64 の精度)。tradke の dbeam の b = bzh·0.5 (B106)、tbgrote の trans(5,6) (B109)、tthine の bsir0 (B110) は直していない (要判断) | |
| 追跡 (RAD) の値 | v2 | FLUC の乱数列は K64 と同じだが、v2 では端の角のある BEND の放射の大きさが変わるので値は K64 と違う (k64 モードは gs と一致) | |
| BEAMBEAM の傾き (v2) | v2 | 傾き回転 (B34) の修正は共変性の確認のみ (t = 0.3 の 1 例。gs は座標が潰れるので比較できない)。実機の衝突点の規模では未確認。強ビームに x-y 結合があるときだけ結果が変わる | |
| `tsteer` | 両 | `F1/FB1/FB2` つきは NotImplementedError | |
| `FindRoot` (v2) | v2 | 変数の範囲 `{var, ini, {min, max}}` は複素では無視。Abs・Conjugate・Re・Im を含む式は収束を保証しない。実数の初期値で複素の値になる式は `General::wrongval` | |
| 特殊関数 (v2) | 両 | PolyLog の負の非整数 s が大きい (s < -10 程度) ときの精度は未確認。HurwitzZeta の q が非正の整数 (`HurwitzZeta[2,-1]` は gs が 2.6449、pysad は NaN) | |
| `ElementValues` (`Couple.n`) | 両 | `VISIT` / `BYE` の退避・復元、フィット変数との連動 (`CouplingMatrix`・`CoupledKeys`) | 13_packages/fn_Couple |
| `NOKEEPEXP` | 両 | `tfvcomp` は常に偽 (pysad の `FREE` はファミリー単位で成分ごとの変数が無い。成分ごとの FREE を実装するときに `FFS._tfvcomp` を直す) | |
| `RENUMBER` | 両 | 実装済み (INDEX 属性を含む)。ただし RENUMBER のあと `LINE["POSITION","QF.1"]` が `LINE["NAME"]` の表示と食い違う gs の癖は写さない (pysad は `LINE["NAME"]` の表示と一致する位置を返す) | |
| SADScript の細部 | 両 | `LINE["POSITION",{{1,2},{3,4}}]` (入れ子のリスト。gs は入れ子のまま返し pysad は平らにする)、`LINE["POSITION",0.5]` (gs は 244、pysad は 1)。`f[x_]:=Hold[x]` の仮引数 x や `Hold[1+2]` の算術の評価 (gs は値を埋め込む / 評価する)。`LINE["LENGTH",loc]` (引数 1 個に限る K64 の仕様で、位置の取り方の案内だけ足した: B24)。`BeginPackage` の戻り値 (`{$ContextPath,"`"}`) の gs との一致は v2 の記録では確かめていない (要確認) | 13_packages/fn_OptimizeOptics ほか |
| 文脈 (`BeginPackage`・`Begin`) | 両 | 利用者のコードでは実装済み。差 (意図的): pysad の ``System` `` が K64 より広い (`AccuracyGoal`・`Plot` などのオプション名を含む)、K64 が `$ContextPath` に無い文脈の記号に文脈名を付けて表示するのに対し pysad は文脈つきの名前だけ省く、``Symbol["Stat`x"]`` は K64 が値を返し pysad は記号を返す。SAD 自身の Packages の読み込みは従来どおり文脈を無視 | 08_scripting/fn_Get_Needs_BeginPackage |
| 台本の読み方 | 両 | `Ffs;` (大文字小文字まじり)・`ffs use ring;` は gs では読めないが pysad は読む。`FFS` (`;` なし) の次の行が USE 以外の文のとき gs はエラーにして行の残りを捨てるが pysad は FFS 部として実行する。MAIN の `A = 3 B = 4;` は gs では A = 3 (エラーつき) だが pysad は代入しない。要素定義の中に直接書いた MAIN の関数 (`L=2*SIN(A)`) は従来の読み取りのまま (未確認)。`LINE[{"EMITY",…},1]` は gs が MAIN の値、pysad が MARK に書いた値 (0) | |
| 描画 | 両 | `CanvasSymbol`、`Canvas3D*`、`UseBG->True` の文字背景。`OpticsPlot` の text 項目数 (gs 69、pysad 74。Tk のフォントによる `RemoveOverlap` の差)。`BeamPlot` (gs は `SIGXX=0` で楕円を描かない) | 11_graphics/fn_OpticsPlot、fn_BeamPlot |
| GUI | 両 | 最初の `Window[]` の親 (gs は Toplevel、pysad はルート窓)。起動時の窓の再現が要る | 12_gui_tkinter/README |
| `Packages/CSR.n` | 両 | 全体は読み込まない (必要な部分だけ) | 06_collective/fn_CSRInit_CSRTrack |

## 3. 未確認 (コードはあるが gs と数値で比べていない)

| 項目 | モード | 状況 | 例 |
|---|---|---|---|
| v2 モードの機能全般 | v2 | 値の根拠は解析式・独立な数値計算で、gs との比較は無い (上の §1a・§2)。特に MEASURE・TRCOD (K64 のソースの読みのみ)、SOL 領域内の BEND・MULT・CAVI・ずれ・F1 (ローレンツ力の積分)、UND・PHSROT、`Polarization` の有限値 (FLUC つきのスピン追跡との比較 3 例)、領域内の AUTOPHI (gs と一致するのは TRPT の場合だけ)、RADLIGHT の 3 次元の軌道 | example/v2_B、example/v2_C |
| `FITP>1` の最終解、`tffamsetup` (fam 点) | 両 | gs が `GO` で再起動することがある。自己検証のみ | 03_matching/FFS_FITPOINTS |
| `CSROYMatrix` | 両 | gs が wrongtype を返す (B60)。比較なし。v2 の `CSRSetupOY` の zj・dpj の修正は蛙跳びの関係で確認 | |
| `tbrote` の `srot` 更新、`tbendebody` の `bsi2` | k64 | ソースの読みだけで、gs の個別確認なし。v2 は行列積・入口と出口の縦磁場に直した (Emittance[] と追跡で 1e-5 以下、Thomas-BMT との比較) | |
| `RadiationField` | 両 | gs が SIGSEGV。速度場との比較のみ。v2 の軌道はピークの場と Liénard-Wiechert の式で 0.25% | 07_spin_solenoid_special/flag_RADLIGHT_fn_RadiationField_fn_RadiationSpectrum |
| `SynchroBetaEmittance` | 両 | 平面リングのみ照合。結合のあるラティスは gs が落ちる | 04_emittance_radiation/fn_SynchroBetaEmittance |
| `measure.py` の `tlum`・`tluma` | 両 | gs で MEASURE が動かない。解析値のみ。K64 の MEASURE・TRACK からは呼ばれず、`ttstat_lines(…, lum=True)` で出せる | |
| `tfloor.f` の `cfunc` (`xsincos` など) | 両 | 個別実装 | |
| `Optics.n` の `PeriodicR`、`CorrectOrbit` の Micado・`CorrectLocal` | 両 | 未確認 | 13_packages/fn_Optics_TransferMatrix、fn_CorrectOrbit_Monitor_Steer |
| `BeamPlot`・`ColumnPlot`・`LaTeX`・`PointSymbol` の見た目 | 両 | 動くことだけ確認。LaTeX は `$HaveLaTeX` が空の状態での一致 | 11_graphics |
| Tkinter の見た目・操作 | 両 | 色・大きさ・フォントは SAD と比べていない | 12_gui_tkinter |
| `dapert.n`、`element.n`、`emit.n`、`beamline.n`、`functions.n`、`Misc.n` | 両 | ファイル単位で gs と比べた記録がない | |
| 乱数列 (FLUC つきの追跡) | 両 | δ の rms が gs 4.21e-4、pysad 4.52e-4 (乱数列が違う) | 04_emittance_radiation/flag_RFSW_RAD_RADCOD_FLUC |
| 乱数の消費の NaN の扱い (B17) | 両 | 決定済み: 現状のまま。`asinz(NaN)`・`gfloor`・`sqrt1` の下限は、発散した粒子が落ちるまでの座標の埋め方にだけ効く。「NaN になった粒子をその場で落とす」に変えると落ちる要素・周が K64 と変わる (LOSSMAP の値)。マッチングの fuzz (tffsmatch) は未確認 | |
| k64 モードの最初の呼び出しのスピン追跡の損失 | k64 | K64 の B33 の癖 (落とす印 0 の粒子が後で入れ替えで戻る) は track6 の経路では近似 (gs は要素 60、pysad は 15)。高速追跡 (track6n) は全桁一致 | |

## 4. 既知の差 (未解決。原因は未調査)

| 差 | モード | 状況 |
|---|---|---|
| BEND の `F1=0.1 FRINGE=3` + RAD の追跡 | 両 | 相対 0.7 の差 |
| L≠0 の CAVI + `DX`・`DY` | 両 | 1e-3 の差 |
| SOL 内の四極 + RAD | k64 | gs が RAD ON・BZ≠0 の領域内で四極の力を入れていない疑いが強い (gs 側の欠陥。B68)。v2 の領域内の四極は物理 (Thomas-BMT・運動方程式の積分) に合わせた |
| EMIT のあとに `USE` する | 両 | gs は EMIT の効果が追跡に残り、pysad は捨てる。v2 でも未解決 (原因未調査) |
| `pm1.sad` の NaN (K3 が強くずれの大きい MULT リングの CALC の閉軌道が gs と別の解)、複数セグメントの `RadiationField` の失敗 | 両 | 未調査 |
| `LOSSAMPL`・RAD・FLUC つきの TRPT、MAIN の TRACK の損失報告、DAPERT での APERT | 両 | gs と比べていない |
| `emit6d._tquade` (Python 版) の解析ヤコビアンが x のずれた軌道で δ の微分を落とす | 両 | `emit_extra`・`calc6d`・`touschek` で 1e-5 ずれる (`tturne` の経路は 1e-12)。直したかは今回コードで確かめていない (未確認) |
| PSPAC の直線の例で x と y が gs と 16% 違う | v2 | 5 MeV・1e11・半径 1 mm の一様なビーム、`D(0.5) SPCH D(0.5) SPCH` で gs は rms x 1.461 mm・y 1.259 mm、v2 は 1.228 / 1.226 mm、薄いキックの解析式は 1.248 mm。原因は調べていない |

## 5. gs 側の問題で比べられないもの

| 項目 | 状況 |
|---|---|
| `For` | gs は 4 番目の引数 (本体) を実行しない (B71)。k64 モードの pysad は gs と同じ (本体を実行しない。以前の pysad は実行していた)、v2 モードは本体を実行する |
| Bessel・HermiteH・Ei などの一部 | gs が誤った値を返す (`BesselJ`(100) など。B43〜B47、B53、B54)。pysad は scipy の値 (真値)。`InverseErf` だけ k64 モードが K64 の値、v2 モードが真値 (B47) |
| gs が落ちる | `RadiationField` の SIGSEGV (B49)、`TRCOD` の SIGSEGV と MEASURE の無出力・MAIN の TRACK (B42)、`TCAVI` の L≠0 (B36)、`SynchroBeta` の結合ラティス、PSPAC の既定 (`PSPACNZ=1`、B93)、`SPAC` で粒子数が少ない (2〜33 個などで SIGSEGV・NaN。**原因は粒子数が少ないことで RADIUS の大きさではない**: spkick.f の spmesh の添字。B98)、BEAMBEAM が SOL 領域内にある配置 (B35)、`RenumberElement` (B51) |
| gs が誤った値を返す | SPAC の力が 1/γ 倍に弱い (B89。解析式の K/γ に近い)、領域内の BEND の偏向をほぼ無視する (B69)、SOL (BZ≠0) 領域内の四極 + RAD (B68) |
| `CSROYMatrix` | gs が wrongtype を返す (B60) |
| `draw` コマンド | gs は表示環境がないと出力せず、`FFS["draw ...",12]` は内部エラー。pysad は無視する |

## 6. 実装・解消して一覧から外したもの (現状)

記録 (`history/v2_changes_*.md`) と試験で確かめた。「両」はモードによらない。

| 項目 | モード | 現状 |
|---|---|---|
| `LOSSMAP` | 両 | `TrackParticles` の戻りに 2 行 (落ちた要素・落ちた周) を足す。gs の例と一致 |
| `RADTAPER` | 両 | Emittance[]・EMIT・追跡・CALC6D に入る。gs (CODCONV=1e-4) と x 6 桁、U0 10 桁、チューン全桁 |
| MULT の `SK0 SK1 SK2 …` (スキュー成分) | 両 | 以前は黙って無視していた (pysad の不具合)。gs と追跡 2e-13、Twiss 1e-9 |
| `ResponseMatrix` の TRPT、`CalculateOptics` の TRPT の加速、加速区間の分散 | 両 | gs と全桁 |
| `NOKEEPEXP`、`RENUMBER` の INDEX 属性、`LINE["INDEX",…]` | 両 | gs の実測と一致 |
| `BeginPackage`・`Begin` の文脈 | 両 | 利用者のコードで実装 (名前空間を作る)。以前は無視 (全部大域) |
| `Shared`・`SharedSize`・`OpenShared` | 両 | fork した子と共有 (無名の mmap)。`SharedSize[{1,2,3}]` = 32 (gs と同じ) |
| 複素引数の Bessel・超幾何・Airy | 両 | 実数のパラメータで実装 (1e-12〜1e-13) |
| PolyLog の非整数次の大きい \|z\| (z < -12・z > 60 の精度、\|z\| ≥ 150 の NaN) | 両 | `_polylog_far` で相対 1e-11 (s ≤ 7.9、\|z\| ≤ 1e12) |
| `LegendreP` の奇数 m の符号 | 両 | 以前は符号が逆 (pysad の不具合)。+.2749545416973507 |
| `EMITDIV`・`EMITDIVB`・`EMITDIVQ` | 両 | 以前は無視していた。Emittance[] に効く (MAIN の `EMITDIV = 20;` も) |
| MAIN の `NAME = 式;` の代入、`FFS;` の書き方 | 両 | gs の動作に合わせた (32 個の式を gs と 2e-14) |
| `For` の本体 | k64 / v2 | 上の §5 |
| RADINT と EMIT の Jx の 2.4% の差 (0.9415 と 0.9651) | v2 | 原因は K64 の Emittance[]・追跡が偏向磁石の端の角で放射の行路を変えないこと (B75)。**RADINT が正しく K64 の EMIT が誤り**。v2 は Emittance[] を直し、E1=E2=0.5 で解析値 1.005751 に対し 1.00586 (`EMITDIV=20` で 1e-4)。k64 モードは 0.9651 のまま |
| `BUCKET HIGHT` (`RADINT`) と `BucketHeight` (`Emittance[]`) の 0.010045 と 0.010035 | v2 | 原因は K64 の MomentumCompaction の誤り (RFSW OFF の入口の分散の項の欠落、RFSW ON の μ²/12。B77)。v2 は 0.010047 と 0.010045 (残りの 2e-4 は空洞と放射の分)。k64 モードは K64 のまま |
| RADINT の量子寿命・線形クロマティシティ・BEND ごとの半径・六極・単精度の定数・CAVI | v2 | 直した (B06・B07・B12・B13・B72・B73)。k64 モードは K64 のまま |
| DELCOR の MARK・CAVI の書き込み先、波数 −1 の成分 (B40・B62) | v2 | 直した。k64 モードは K64 のまま (CAVI は HARM・PHI に入る) |
| `PolynomialFit` の平均 0 のデータ (B18)、`InverseErf` (B47) | v2 | 直した。k64 モードは K64 の値 |
| MOMENTUM を設定しない計算の警告 (B67) | v2 | 1 度だけ警告 (値は変えない) |

## 7. 優先度 ([16_fortran_status.md](16_fortran_status.md) §9.4 から)

| 優先 | 項目 | 理由 |
|---|---|---|
| 中 | 誤差入り 6D emit と DA の反復 (性能) | 反復が 600 s を超える |
| 中 | 光学・Emittance への空間電荷 (WSPAC の自己無撞着なビーム行列) | v2 の WSPAC は追跡だけ。入れる判断が要る (docs/47 (e)) |
| 中 | Emittance[]・EMIT の POL の量の SOL 領域内 (ずれて傾いたソレノイド) | v2 で追跡のスピンは通るが POL の量は K64 の式のまま |
| 中 | §4 の未解決の差 | 原因の切り分け |
| 中 | v2 の要判断 (docs/47 (e)) | ANGLE=0 の BEND・K0 つき MULT の端の角 (B108)、RADINT の空洞の位相、UND の光学、`TrackParticles` が `MASS`・`CHARGE` を追跡に渡さない (電子以外が追えない) ほか |
| 低 | `tdgeo` の四極の微分、`wfres` 系、`tphrec` | 必要になったら |
| 不要 | MON/MCOD/TCOD の機能、Obsolete・非リンクのファイル、`tdrwdt.f` の粒子追跡 | K64 で無効 |

## 8. 注意

- 条件: 例の数値は、サンプルのラティス (`example/00_lattice/FODO_ring.sad`: 3 GeV、周長 105.6 m) と、PF・ASTRID2 のラティス (`example/pf`、`example/astrid`) で得たものです。v2 の機能の数値は、各記録の条件 (上の冒頭) の範囲です。
- 表示や数字の桁の微妙な差は、ここには入れていません。完走することと、意味のある値が合うことを基準にしています。
- 一覧は変わります。項目を直したら、この文書と [16_fortran_status.md](16_fortran_status.md) §9 を一緒に更新してください。
