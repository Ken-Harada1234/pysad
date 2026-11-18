# 46. オリジナル SAD (K64) に対する pysad の未対応・未確認の一覧

最終更新: 2026-11-18 (v2.1.0 の後の現状のスナップショット)

この文書は、オリジナルの SAD (K64 版、`src/*.f`・`Packages/*.n`・SADHelp) に対して、pysad が**まだ対応していないもの**と**対応したが gs で確かめていないもの**を、1 か所にまとめたものです。

- **互換モード (v2.0.0 から)**: pysad は 2 つのモードを持つ。`k64` は K64 (gs) と同じ動作 (v1 と同じ。gs との数値比較はこのモード)、`v2` (既定) は K64 のバグを直し機能を足した動作。環境変数 `PYSAD_COMPAT=k64|v2` か `pysad.compat.set_mode()` で切り替える ([47_v2_changes.md](47_v2_changes.md))。**各項目の「モード」の欄が `k64` ならそのモードだけの話、`v2` なら v2 モードだけの話、`両` ならモードによらない。**
- v2 の値の根拠は、解析式・独立な数値計算 (ローレンツ力・Thomas-BMT の数値積分など)・保存則・別のコードとの一致で、gs ではない。**v2 で実装・修正した機能の多くは、gs が動かない・K64 が誤っているので gs と比較できない。** 「未確認」にはその旨を書いた。使った条件は電子、0.2〜4.2 GeV のリングと輸送路 (空間電荷は 5〜50 MeV、RADLIGHT は 10.2 MeV)、ソレノイド bz = 0.3〜3 1/m (一部 20 1/m まで) で、**実機の強い結合・挿入光源のずれのあるリングでは確かめていない**。大きなリング A (3 GeV、周長 603 m、`pysad/tests/kek-ls_v03_12.lat`)・PF リング・ASTRID2 (MPW 入りを含む) では多くの項目を gs と比べた。大きなリング B (8 GeV、周長 1436 m) は数値だけ。
- 出所: [16_fortran_status.md](16_fortran_status.md) §9、[history/work_unimpl.md](history/work_unimpl.md)、[21_sad_bug_report.md](21_sad_bug_report.md)、v2 の作業の記録 `history/v2_changes_*.md` (v2.0.0) と `history/v2next_*.md` (v2.1.0: A1・A2・B1・B2・B3・C・D1・D2・E)、[47_v2_changes.md](47_v2_changes.md) の v2.1.0 の (d) 要判断・未解決、`example/` の各例の「gs と pysad の比較」・「未確認・限界」節、[34_sad_packages.md](34_sad_packages.md) ほか各解説文書。
- 基準は gs (SAD K64) の実動作です。gs が落ちる・誤った値を返すものは「gs 側の問題」(§5) に分けています。K64 のバグの番号は [21_sad_bug_report.md](21_sad_bug_report.md) の B01〜B134 (B72〜B110 は v2.0.0、B111〜B134 は v2.1.0 の作業で見つけたもの)。
- `例` の欄は、その項目を示す `example/` の例です (パスは `example/` からの相対)。v2 の機能の例は `example/v2_B/`・`example/v2_C/`・`example/v2_next/`。

凡例: 「未翻訳」は pysad に物理が入っていない。「一部未対応」は一部の入力・経路だけ未対応。「未確認」はコードはあるが gs と数値で比べていない。「モード」は上のとおり。

## 1. 未翻訳

| 項目 | モード | 理由・状況 | 例 |
|---|---|---|---|
| `wjfit`、`wfres`、`wfbin`、`wiord`、`wjdraw` (MAIN の TRACK の Fourier 解析、FOURIER フラグ) | 両 | FFS から使えない。gs も SIGSEGV (`tracka` → `tinip`) で出力が取れず、直訳しても検証できない。pysad は FOURIER を無視し `wtunem` の表だけ出す。チューンの抽出は NAFF (`pysad/naff.py`) を足した | v2_next/B3_naff_tunes.py |
| 光子 (`tphrec`) のうち SOL 領域内の要素・COORD・TCAVI・MAP のあるラティス | 両 | BEND・QUAD・MULT・SEXT・OCT・磁場つきドリフト・BEAMBEAM + BEND は記録する (v2.1.0)。SOL 領域内は記録しない (乱数だけ K64 と同じに引く。`_ph_set` を領域内の要素に足し、領域の座標系の `geo_frac` の使い方を確かめる必要)。COORD・TCAVI・MAP のあるラティスは従来の高速追跡で BEND があっても光子を記録しない | v2_next/B3_photons_bend.sad |
| TrackParticles の F1 (`tsteer` (ANGLE=0 の BEND)・`tmulti` (MULT) の `F1/FB1/FB2`) | 両 | NotImplementedError。Emittance[]・EMIT・CALC は実装済み (§6) | |
| CALC (4 次元) の RADTAPER | 両 | **実装しない**: K64 の 4D の CALC は RADTAPER を見ない (gs で確認)。CALC6D は RADTAPER を入れる | 04_emittance_radiation/flag_RADTAPER |
| MON・MCOD・TCOD の機能 | 両 | K64 にコマンドが無く、古い SAD の MON の書式・意味を確かめられる資料が無い。v2 は「K64 に無いコマンド」と警告して文全体を読み飛ばす (定義はしない)。軌道補正は CorrectOrbit パッケージ ([34_sad_packages.md](34_sad_packages.md)) | |
| `DRAW` の中の粒子追跡 (`tdrwdt.f`) | 両 | K64 でも動かない死んだコード (`trsize` が常に偽。B105) | |
| `Out[n]:=` の gs の表示 | 両 | gs は未定義の記号 1 つの文で前の `Out[n]` を繰り返すなど値と関係のない行を出す。K64 の表示の癖として写さない | |
| 複素の特殊関数で gs が誤る点の再現 (B118〜B120) | k64 | k64 モードも真値を返す (gs の誤った値を再現しない)。再現には K64 の `gamma.f` (複素の不完全ガンマ `cgamma2`・`czconvb`・`cpolylogl`) の直訳が要る。要判断 | v2_next/B2_special_functions_complex.sad |
| `DynamicApertureSurvey` の `DAPWIDTH` の窓 (既定 7) | 両 | pysad は既定 51 で全点を追跡する。gs は x の振幅を上から 7 個ずつ追跡し、続けて生き残れば残りは追跡せず `*` で表示する。`trackd.f` の追跡順の直訳が要る (`Maximum number of particles` も gs 42、pysad 306 と違う) | |
| `TRCOD` の FIFO での待ち合わせ (`PUTMEASURE`・`GETMEASURE`)・`ORBITJITTER` のパイプ | v2 | 普通のファイルで同じ形式を読み書きする (待たない)。Linux・WSL の FIFO は未確認 (試験は普通のファイルだけ) | |
| K64 の `qwspac` (4 次元の CALC の空間電荷) の直訳 | v2 | v2 は CALC6D と同じ 6D の経路で入れる。gs が NaN を返す機能 (B133) なので直訳の検証ができない | v2_next/B1_wspac_emittance.sad |

## 1a. モードで状態が違うもの (k64 モードでは未対応・v2 モードで実装済み、または結果が違う)

v2 モードの実装は、特記がなければ**gs と比べていない** (gs が動かない・K64 が誤っているため)。根拠は解析式・数値積分。v2.1.0 で k64 と v2 の結果が変わった項目の一覧は [47_v2_changes.md](47_v2_changes.md) の v2.1.0 の (a) (B111・B112・B114〜B117・B121〜B123・B125〜B127・B129〜B133 ほか)。

| 項目 | k64 モード | v2 モード | 例 |
|---|---|---|---|
| MEASURE・TRCOD (FIXCOD・MVCOD・STNCOD・ALIGN)・LOGMEASURE・ORBITJITTER・PUTMEASURE・GETMEASURE、MAIN の `TRACK TURNS=-1` (`trackb.f`・`tracke.f`・`trcoda.f`・`ttinit.f`・`ttstat.f`) | 未対応 (警告して無視)。gs で TRCOD が SIGSEGV、MEASURE は出力なし (B42) | 実装 (K64 のソースの読み。**gs との比較は 1 件も無い**。解析式と pysad の独立な経路で確認: docs/03、C の記録 `history/v2next_C.md` §3.1) | example/v2_C/C2m_trcod_reference.sad、v2_next/C3_trcod_orbit.sad |
| SOL 領域内の BEND・MULT・CAVI・APERT (`tsole1`・`tmulte1`)、SOL の `DX DY DZ CHI1..3`・`F1`、領域の幾何 | 未対応 (NotImplementedError)。BEND は K64 の tsol.f どおりに実装したが gs が偏向をほぼ無視して合わず取り下げ (B69) | 実装 (ローレンツ力の数値積分との一致。docs/44 §7a) | example/v2_C/C3_sol_region.sad、C2a_sol_tilt.sad |
| 領域内の CAVI・RF つき MULT の TRPT の加速・wake・`AUTOPHI` | 未対応 | 実装 (加速は運動方程式の積分と一致。AUTOPHI の意味は K64 のソースの読みで、gs と一致するのは TRPT の場合だけ) | example/v2_C/C2b_sol_rf.sad |
| 空間電荷 (SPAC・PSPAC・SELFCOD) の追跡 | 追跡に入らない (K64 の SPAC は力が 1/γ 倍に弱い: B89 など。ビット単位で再現する翻訳はしていない)。WSPAC は Emittance[]・EMIT に入る (gs と全桁) | 実装 (SPAC は解析式・Laslett の式、PSPAC は gs と 3〜4 桁、WSPAC は gs と 1%。粗いメッシュの PSPAC は警告。WSPAC の CALC・CALC6D にも入る) | example/v2_C/C2s_mult_radius.sad、C2s_wspac.sad、05_tracking/flag_SPAC |
| `PHSROT` | 追跡は NotImplementedError、光学は恒等 (K64 で動作しない: B38) | 実装 (意図された写像。固有値・シンプレクティック性・平衡エミッタンスで確認) | example/v2_C/C2_phsrot.sad |
| UND の磁場 (`FBX`・`FBY`) | 磁場なしのドリフトのみ (K64 が設定できない: B37)。光学は恒等 | 実装 (解析式で確認)。磁場がなくても光学にドリフトとして入る (要判断: K64 の「光学では恒等」を仕様と見るなら磁場つきのときだけにする選択もある) | example/v2_C/C1_undulator.sad |
| `Polarization` (平衡偏極度) の有限値 | K64 が NaN (B55)・値が大きすぎる (B125)。逆向きの BEND (挿入光源) で緩和率が相殺する (B116) | 有限値 (平面リング 92.376%。誤差のあるリングは側帯の項つきで FLUC の追跡の 0.43 倍に過小評価: §2。逆向きの BEND は符号に依らない率)。docs/43 | 07_spin_solenoid_special/flag_POL_RADPOL、v2_next/B3_pol_sideband.sad、A2_pol_mpw.sad |
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
| MULT の端の `FB1/FB2` の `tblfre` (`multae.py`)・AE1・AE2 | 両 | Emittance[]・CALC の MULT の F1・F2・FRINGE は実装 (§6)。`FB1/FB2` の `tblfre` は NotImplementedError、AE1・AE2 (端の角: `tmulke` の apsi) は pysad が MULT では読まず黙って無視する。そのため `tbende` から tmulte に回る K1≠0・K0≠0 の ANGLE=0 の BEND の AE も無視する。逆向きの要素 (`direlc < 0`) の F1 の入れ替えは未対応 | v2_next/E_mult_f1.sad |
| MAP | 両 | `GEO`・`CompiledMap`。MAIN の TRACK・DA・Touschek は ExternalMap を呼ばない。MAP のスピンは未対応 | 06_collective/fn_ExternalMap |
| POL (SOL 領域内) | k64 | 領域内の MULT・BEND・CAVI のスピン (Emittance[]・追跡とも)。スピン追跡は損失の処理なし | 07_spin_solenoid_special/flag_POL_RADPOL |
| POL (SOL 領域内) | v2 | 領域内の APERT のスピン追跡 (NotImplementedError)。領域内の CAVI の POL (Emittance[]) は VOLT=0 だけ確認 (加速する CAVI は一周でエネルギーが元に戻るので解析式は縦の運動を含む)。bz ≥ 3 1/m の QUAD・MULT のスピンは放射の分割 (EPS 由来) で 1e-4 台。領域内の QUAD の端の ψ は追跡と同じ模型 (hard edge の Bs なし。軌道が数 mm の範囲で確認、cm のずれは未確認)。領域内の四極の DX・DY のある POL の追跡 (track6) は NotImplementedError。領域内の光子は記録しない。大きなリング A への SOL の確認は小リングのみ。領域の外の B∥ (垂直な軌道の傾きの中の偏向磁石の B_y y') は K64 も追跡も入れていない | |
| POL + TRPT + 加速 (track6) | v2 | SOL 領域の外の VOLT≠0 の CAVI だけ。SOL 領域内の加速・RF つき MULT・TCAVI の加速・リングの何周もの TRPT (K64 の TRPT は輸送路用) は SADError。スピン追跡の損失の処理 | v2_next/D1_spin_accel.py |
| `Polarization` | v2 | 軌道の振幅について 1 次の SLIM に、スピンの歳差の位相の変調 (シンクロトロン側帯) を厳密に取った項を足した (B125)。**脱偏極率を FLUC つきの追跡より小さく見積もる**: 4.2 GeV・12 MV の 12 セル FODO で追跡の 0.43 倍 (480 粒子、3.7σ)、線形の軌道 + 厳密なスピン回転の MC とは 10% 以内。残りの 2 倍は軌道の非線形 (全座標の蹴りへの応答が線形モデルの 2.07 ± 0.24 倍。δ 1.9 倍、水平 1.3 倍、垂直 1.0 倍。振幅 2σ 以上で急に増える。どの非線形かは未分離)。低エネルギー (ASTRID2) は応答の比 1.0。目安: 変調指数 aγσδ/νs ≲ 0.3 で線形の式と 1〜15% の差、≳ 0.5 では追跡 (FLUC) で確かめる。限界の表は docs/43。結合項の 2 次の項は入れていない (MC に対し −10%)。運動論的偏極なし。不安定なリングでは線形の式に戻す (大きなリング A に鉛直の四極の誤差 rms 0.2・0.4 mm で光学が不安定)。n が傾く・非平面のリングで緩和率 Σ r (b·n)² が最適かは未確認。`Polarization2/4/6` は実装していない。POL のときは四極・MULT を細かく分割するので Emittance[] が遅い | v2_next/B3_pol_sideband.sad |
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
| 空間電荷 (v2) | v2 | **使えない**: `ANGLE≠0`・`VOLT≠0` の MULT と SOL 領域内の MULT の SPAC (メッセージを出して空間電荷なしで追う)、QUAD・BEND・SEXT の中の SPAC (偏向磁石のあるリングでは SPAC は全周には入らない。全周なら PSPAC + SPCH か WSPAC)、複数バンチ (`NBUNCH`)、壁の抵抗、円でないダクト、MAIN の `TRACK`・`DynamicApertureSurvey`、スピン (`POL`) の追跡。`NORFSW` はエラー (`RADIUS<0` の絞りだけなら走る)、粒子の z が全部同じでもエラー。**未確認**: ウェーク・放射との組み合わせ (走るが独立な根拠なし)、放射 + 空間電荷の MULT、SOL 領域内の四極、バンチの静止系で相対論的な運動 (γ·x' ≳ 1)。陽子など電子以外は `MASS`・`CHARGE` が追跡に渡るようになった (v2.1.0。空間電荷の式は 100 MeV の陽子で 0.5% 確認、SAD スクリプトからの使用は未確認)。精度: 4096 粒子の既定のメッシュで中心の勾配 (チューンシフト) が 1〜3% 小さめ。**PSPAC の粗いメッシュ** (横の rms < 1.5 セル、縦の rms < 1.5 セル、粒子の 1% 以上がメッシュの外) は警告が出て、定量値は gs と違う (v2 は 1 セル以下で gs より解析式から遠い。K64 の形状関数 CIC・双 3 次・「ほぼ最近接」を入れていない。メッシュが十分細かければ解析式と 1〜2%、gs と 1%)。条件: 電子 5〜50 MeV、1e8〜1e12 個/バンチ、ビーム半径 0.1〜7 mm、ダクト 3 mm〜3 cm。gs との比較は SPAC の包絡・PSPAC・WSPAC の 4 例。gs は PSPAC の状態を C の static 変数に持つので、1 セッションで TrackParticles を繰り返すと 2 回目から PSPAC が働かない (B96。gs と直接比べるときは 1 セッション 1 回) | 05_tracking/flag_SPAC、v2_next/A1_pspac_resolution.sad |
| WSPAC | v2 | **追跡の WSPAC** (`LINE["SIG…"]`) のビーム行列は空間電荷なし (ゼロ電流) の平衡値・MINCOUP なしのまま。ビームの大きさは追跡の間じゅう固定 (強-弱。エミッタンスの増大・コヒーレントな振動は扱えない)。垂直のエミッタンスが 0 の平らなリングは Emittance の経路で使えない (`WSPACBEAM` を与える)。`WSPACBEAM` は σz 一定・x-z の相関なし・MINCOUP なし。TRPT の輸送路ではエミッタンスの断熱減衰を入れていない。縦の力は既定で「軸上で 0」(K64 と同じ)。**Emittance[]・EMIT・CALC6D・CALC (RFSW ON) の WSPAC は実装済み (§6)**: 強い空間電荷 (0.2 GeV で PBUNCH ≳ 1e9 は反復が収束しない・不安定)、RF なし (σz が決まらない) の CALC、TouschekLifetime との併用、要素ごとのレンズの位置 (K64 の tintrb は磁石の中で厳密でない近似)、実機の大きなリングは未確認。SOL 領域内・`SELFCOD` + 閉軌道のあるリング・放射つきの追跡は走るが値の確認なし。gs で `EMIT` を `NOWSPAC` で行ってから `WSPAC` で追跡すると NaN (K64 は「WSPAC つきの EMIT」が前提) | v2_next/B1_wspac_emittance.sad |
| MEASURE・TRCOD (v2) | v2 | gs との比較が無い (C の独立な基準 (閉軌道の式、固有エミッタンス、ビームサイズ) で誤りは見つからなかった)。`TRCOD FIXCOD`・`MVCOD`・`STNCOD` の基準は REFERENCE と同じ参照光学に置いた (`measure._ref_rows` の 1 か所で替えられる)。`ALIGN` の対象に DECA DODECA MULT を足した (K64 は BEND QUAD SEXT OCT だけ。`measure._ALIGN_V2` を `_ALIGN_K64` にすれば K64 のとおり)。SOL・MONI は動かさない。区間に BEND を含む ALIGN は平行移動との一致が 4e-10 m まで。GETMEASURE の対象の関数 (DX・DY) は K64 の関数番号 32・34 (今の番号では PEPY・PZPX) から確定できず式から判断 (B103)。**`Twiss["DETR"]` の不具合は直した (v2.1.0) が、DETR ≥ 1 (モードが入れ替わる) の INS と MEASURE は未確認** (輸送路で先頭 MARK に DETR を書く場合は正しい)。TWAKE / LWAKE ON は z の分布を作る (K64 は全粒子 z = 0。一様・階段状の wake の重心と広がりだけ確認)。MEASURE の TRPT は追跡の間だけ ON にして戻す (K64 は戻さない。要判断で決定済み)。MAIN の `TRACK TURNS=-1` は粒子種の語 (PROTON など) を見ない。TRGAUSS OFF の分布と表示の桁 (autofg) は gs の実出力と比べていない。TRCOD の途中の統計 (K64 は標準エラーに出す) は出さない。周期解を matched 分布として 1 周追う試験は δ が大きいと崩れる (色収差。DP = 1e-3 で 40%)。QUIET の重心合わせ (B87) と TWAKE/LWAKE の z = 0 (B88) は K64 の誤りと判断して直した推定 | example/v2_C/C2m_putget_measure.sad、v2_next/C3_trcod_orbit.sad |
| RADLIGHT (v2) | v2 | SOL 領域を含む区間 (SADError)、COORD による幾何の変化、曲げのある MULT の中の点 (出口だけ)。DX・DY・F1・FRINGE のある BEND は中の点だけ 1e-6 m 程度ずれる。高エネルギーでは点数が γ に比例して増える (3 GeV・0.26 rad の BEND で 2.5 万点)。`RadiationField` の振幅・位相・スペクトル・全エネルギーは Liénard-Wiechert の厳密な数値と合う (1 GeV の BEND、点間隔 1/γ の角度あたり 16 点以上で振幅 3e-3、64 点で 2e-4。観測点がスペクトルの波長で Fraunhofer の距離に無いと光源の大きさ ρ/γ の分だけ物理的に違う)。RAD が ON のため軌道が放射で変わる点に注意 (同じ軌道を別の計算に渡すとき)。gs との比較は無い (gs は SIGSEGV: B49) | 07_spin_solenoid_special/flag_RADLIGHT_fn_RadiationField_fn_RadiationSpectrum、v2_next/C3_radfield_bend.sad |
| DELCOR (v2) | v2 | corlength≠0 の規格化 (delx を DX の rms と読む) の根拠は確定 (K64 は複素数の ⟨\|c\|²⟩ = 1。corlength = 0 の枝は成分ごとで約束が揃っていない) が、作者の意図かは推定 (B99。K64 の値が必要なら k64 モードか delx に 1/√2 を掛ける)。CAVI の符号 (B40) は追跡に効かず観測できない (+delx のまま)。rms は多数回の平均で、GCUT を掛ければその分小さい。相関長が格子間隔 (周長/nc) より短いと線形補間で rms が下がる | |
| RADINT・Emittance[]・EMIT (v2) | v2 | RADINT は MULT・SOL が積分に入らない (表示が出る)、ANGLE≠0 の BEND の K0 は K64 のまま、ANGLE=0 の K0 の極は入るが εx は Emittance[] と 3.5%・Jx は 0.7% 違う (極の中の分散の伝搬を扇形の行列で近似し、K0 の場の端は無視。B115)。SEXT のフィードダウン・BEND の DX と K1 の勾配の軌道のずれは入れていない (B126)。軌道のずれで経路長が変わるとき RADINT (δ = 0 の約束) と Emittance[] の J が違う。結合リングで νy ≈ νs のときの J は RADINT と追跡・Emittance[] が 0.015 違う (モードが混ざり定義できない条件)。K1≠0 の BEND の傾き θ ≠ 0 の CKL1Y の数値の確認が無い。空洞は VOLT を単純に足す (PHI を見ない。ベクトル和にするかは仕様の判断)。既定の分割 (EMITDIV=1) では Emittance[] の εx が低エミッタンスのリング (大きなリング A) で 2.1% 大きく Jx が 1.5e-3 大きい (`EMITDIV=5` で εx 8e-5、`EMITDIV=20` で 6e-6・Jx 1e-4)。Jx・Je・τ・σδ は EMITDIV を上げても U0/E のオーダーの差が残る (RADINT は U0/E の 1 次の理論)。**v2 の Emittance[] と RADINT が矩形磁石で 1.4e-3 違う** (大きなリング B: v2 6.9651e-9、RADINT 6.9554e-9、k64 は 1e-4 で整合。v2 の修正の側の問題かもしれない: 未調査)。**v2 の Emittance[] (RADCOD) は gs と 1e-4〜1e-3 ずれる** (B111。k64 は全桁一致)。K1≠0 の BEND の CALC の解析行列は y ≠ 0 の軌道で自身の写像の数値微分と y = 2 mm で 1.5e-6 違う。K1=0 の BEND の追跡のスピン (tbend + tradk) が Emittance[]・tbendi (v2) と 2.2e-5 rad 違う (放射の分割が粗い K64 の精度)。RADCOD の 1 周の解析 T と写像の数値微分が 2e-3 違う (K64 の tradke の微分が近似、k64・v2 とも)。tthine の bsir0 (B110) は直していない (影響なし)。K1≠0 の BEND の四極が半径方向の面で始まる模型は K64 の仕様 (直さない)。ANGLE=0 の BEND の AE1・AE2 の残る 3〜4% (傾き 5 mrad の極。30 mrad では 2〜3%) | |
| `SynchroBetaEmittance` | 両 | 大きなリング A (3 GeV、603 m、νs = 0.004) で独立な値と 40 倍違い、モード数で収束しない (B128。原因は未特定)。v2 は収束していないとき警告 (AzimuthalModes m と 1.5 倍の m が 20% 以上違うとき)。PF・ASTRID2・大きなリング B は Emittance[] と 0.5〜1%。結合のあるリング (24 セル FODO + skew 四極) では gs も落ちず k64 は 7 桁一致、Emittance[] と独立な 6x6 の平衡エミッタンスとも一致。`AzimuthalModes->1` は pysad が IndexError (gs は 8.6e-7 を返すが意味がない) | 04_emittance_radiation/fn_SynchroBetaEmittance、v2_next/C3_synchrobeta_coupled.sad |
| 追跡 (RAD) の値 | v2 | FLUC の乱数列は K64 と同じ (`NPARA=n;` で gs の子プロセスの分け方も再現)。v2 では放射の式 (B111・B112)・端の角のある BEND (B75) が変わるので値は K64 と違う (k64 モードは gs と一致) | |
| BEAMBEAM の傾き (v2) | v2 | 傾き回転 (B34) の修正は共変性の確認のみ (t = 0.3 の 1 例。gs は座標が潰れるので比較できない)。実機の衝突点の規模では未確認。強ビームに x-y 結合があるときだけ結果が変わる | |
| `tsteer` (ANGLE=0 の BEND) | 両 | `F1/FB1/FB2` つきの TrackParticles は NotImplementedError (§1)。ステアリングに `ROTATE` と `DX`/`DY` を同時に付けると gs と CALC のチューンが 1e-4 違う (F1 に関係なく既存。K64 はこのとき xy を回す行列を返す癖があり、pysad も似た値を返すが一致しない。未調査)。追跡 (`tsteer`) の端の放射の AE1・AE2 は未確認 | |
| `FindRoot` (v2) | v2 | 変数の範囲 `{var, ini, {min, max}}` は複素では無視。Abs・Conjugate・Re・Im を含む式は収束を保証しない。実数の初期値で複素の値になる式は `General::wrongval` | |
| 特殊関数 | 両 | 複素の PolyLog・HypergeometricU・Zeta・HurwitzZeta は実装済み (§6)。使えない: 複素の次数の PolyLog、Re a ≤ 0 の複素の a・b の HypergeometricU (SADError)、z = 0 の U、\|s\| > 500 (Zeta・HurwitzZeta)、\|z\| > 1e8 (PolyLog)・\|z\| > 30 (HypergeometricU)。PolyLog の負の非整数 s が大きい (s < −10 程度) ときの精度は未確認。HurwitzZeta の q が非正の整数 (`HurwitzZeta[2,-1]` は gs が 2.6449、pysad は NaN) | |
| `ElementValues` (`Couple.n`)・`FREE` の成分 | 両 | フィット変数との連動 (`CouplingMatrix`・`CoupledKeys`)、COUPLE で結合した要素の成分、`USE NOEXPAND` の展開済みの成分の値の保存は未確認。**先頭でない成分の既定のキーの変数は作らない** (gs は値を持たない変数として数える: 反復の行の変数の数が 1 つ違い、その変数だけで GO すると pysad は `?No variable.`。`VISIT`・`BYE` の退避・復元は実装済み。`BYE ALL` は gs が SIGSEGV で比べられず、VISIT の未定義の LINE のメッセージ (gs は `Missing beamline in VISIT.`) は違う) | 13_packages/fn_Couple、v2_next/B2_free_component.sad |
| `RENUMBER` | 両 | 実装済み (INDEX 属性を含む)。ただし RENUMBER のあと `LINE["POSITION","QF.1"]` が `LINE["NAME"]` の表示と食い違う gs の癖は写さない (pysad は `LINE["NAME"]` の表示と一致する位置を返す) | |
| SADScript の細部 | 両 | `LINE["POSITION",{{1,2},{3,4}}]` (入れ子のリスト。gs は入れ子のまま返し pysad は平らにする)、`LINE["POSITION",0.5]` (gs は 244、pysad は 1)。`f[x_]:=Hold[x]` の仮引数 x や `Hold[1+2]` の算術の評価 (gs は値を埋め込む / 評価する)。`LINE["LENGTH",loc]` (引数 1 個に限る K64 の仕様で、位置の取り方の案内だけ足した: B24)。`BeginPackage` の戻り値 (`{$ContextPath,"`"}`) の gs との一致は確かめていない (要確認)。除算 `a/b` は gs が `a*(1/b)` で最後の桁が 1 ulp 違うことがある (`Log[100]/Log[10]` は gs 1.9999999999999998、pysad 2。評価器全体は変えていない: 全試験の基準値に影響)。`Through[{f,g}[1,2]]` はどちらもエラーだが文言が違う。`PROFILE` の代入は gs も `General::invset`。MARK の Twiss のキー (AX など) は CALC のあとの SetElement で gs は MAIN の 0、pysad は CALC が書き戻した値。`Element["K1","QF"]=0.9` の読み戻しが 1 ulp (0.8999999999999999) | 13_packages/fn_OptimizeOptics ほか |
| 文脈 (`BeginPackage`・`Begin`・`$Context`・`WithContext`) | 両 | 利用者のコードでは実装済み (代入・`WithContext`・入れ子は gs と全行一致)。差: pysad の ``System` `` が K64 より広い (`AccuracyGoal`・`Plot` などのオプション名を含む)。文脈つきの `Names` の規則が決められない (gs は状況で `{}` と `{"A`"}` を返す)。init.n が大域に作る名前は 51 個しか調べていない (gs の `Names` が列挙しない)。`Off[f::tag]` による抑止 (pysad の `Off` は何もしない)。`Message` の `Tkinter::undef` は出さない (pysad の Tcl 橋が gs と違う経路で通るため)。SAD 自身の Packages の読み込みは従来どおり文脈を無視 | 08_scripting/fn_Get_Needs_BeginPackage、v2_next/B2_context.sad |
| 台本の読み方 | 両 | `Ffs;`・`FFs;` は gs と同じ構文エラー、LINE 名は書かれたとおり (`USE ring` は ring が無いとエラー。gs の「USE の失敗のあとの仮の 1 要素の LINE」は作らない)。`FFS` (`;` なし) の次の行が USE 以外の文のとき gs はエラーにして行の残りを捨てるが pysad は FFS 部として実行する。MAIN の `A = 3 B = 4;` は gs では A = 3 (エラーつき) だが pysad は代入しない。MAIN の要素定義の式は gs の文法で読む (v2.1.0)。gs が構文エラーにするが pysad が受ける形: `45DEG`・`SIN[0.1]`・未定義の名前 (`L=A`)。gs が受けるが pysad が構文エラーにする形: 末尾の `1+`。MAIN の変数を `E`・`Pi`・`Degree` で代入しても要素定義の式の `E` は SADScript の定数。`LINE[{"EMITY",…},1]` は gs が MAIN の値、pysad が MARK に書いた値 (0) | v2_next/D2_main_expr.sad |
| `Parallelize` | 両 | **Windows は逐次** (`multiprocessing` の spawn は Interp・FFS・ラティスを pickle できない)。POSIX は fork で並列。numba と fork の併用は実際には試せていない (WSL に numba が無く Windows に fork が無い)。(i) numba を `parallel=True` にする、(ii) 初回のコンパイルを子で行ってキャッシュを複数の子が同時に書く、(iii) OpenMP 版の BLAS の numpy、のどれかなら未確認。親でカーネルを一度走らせてから呼ぶ | v2_next/B2_parallelize.sad |
| `Packages/*.n` の関数 | 両 | `functions.n`・`Misc.n`・`element.n`・`beamline.n`・`emit.n`・`dapert.n` は小さな台本を gs と行ごとに突き合わせた (v2.1.0)。残った差: `DAPWIDTH` の窓 (§1)、`DynamicApertureSurvey` の `Axes->"XZ"` で z が 3 点のときの gs のエラー、`SymplecticInverse` の `-0` の表示、`Names`、不正な Axes のメッセージ | v2_next/B2_packages_compare.sad |
| `CorrectOrbit` の Micado・`CorrectLocal`・`PeriodicR` | 両 | gs と比較済み (PeriodicR は 12 成分が 3e-16、CorrectLocal の補正量が 1e-14、Micado は k64 が gs と 1e-6、v2 は B129 の修正)。未確認: 水平・2 平面同時の Micado、`Bump`・`SVDmode->1` つきの Micado | v2_next/D2_micado.sad |
| 描画 | 両 | `CanvasSymbol`、`Canvas3D*`、`UseBG->True` の文字背景。`OpticsPlot` の text 項目数 (gs 69、pysad 74。Tk のフォントによる `RemoveOverlap` の差)。`BeamPlot` (gs は `SIGXX=0` で楕円を描かない) | 11_graphics/fn_OpticsPlot、fn_BeamPlot |
| GUI | 両 | 最初の `Window[]` の親 (gs は Toplevel、pysad はルート窓)。起動時の窓の再現が要る。Tk 9.0 (Python 3.14) では `cget` が長さを pt で返す (`-padx` 4p など) ので、SAD の `.n` が `cget` の値を数として使うと壊れる (未確認。見つかったら `winfo fpixels` で直す) | 12_gui_tkinter/README |
| `Packages/CSR.n` | 両 | 全体は読み込まない (必要な部分だけ) | 06_collective/fn_CSRInit_CSRTrack |

## 3. 未確認 (コードはあるが gs と数値で比べていない)

| 項目 | モード | 状況 | 例 |
|---|---|---|---|
| v2 モードの機能全般 | v2 | 値の根拠は解析式・独立な数値計算で、gs との比較は無い (上の §1a・§2)。特に MEASURE・TRCOD (K64 のソースの読みのみ)、SOL 領域内の BEND・MULT・CAVI・ずれ・F1 (ローレンツ力の積分)、UND・PHSROT、`Polarization` の有限値 (線形の軌道 + 厳密なスピン回転の MC と FLUC つきのスピン追跡)、領域内の AUTOPHI (gs と一致するのは TRPT の場合だけ)、RADLIGHT の 3 次元の軌道、v2.1.0 の修正の値 (B111・B114〜B117・B125〜B127・B130 ほか: 独立な式・厳密な幾何で確認し、gs とは k64 モードの再現だけ比べた) | example/v2_B、example/v2_C、example/v2_next |
| 大きなリング B (8 GeV、周長 1436 m) | 両 | 数値だけ確認 (CALC・RADINT・Emittance[]・COD は k64 が gs と 1e-10、v2 は εx +0.13%)。U0/E が大きいリングの RADTAPER・tradke の符号 (B111)、陽子の大きなリング、POL (四極 60 台に 0.1 mm の誤差を入れると光学が不安定で NaN) は未確認 | |
| `FITP>1` の最終解、`tffamsetup` (fam 点) | 両 | FITP>1: pysad の最終解は独立な最小化と 5e-8 で一致、gs は最小値の手前で止まる (残差が相対 5e-5 大きい。`Overflow & closed orbit not found` の試行点での打ち切り。今の gs は再起動しなかった)。fam 点: 初期条件の Courant-Snyder の不変量と、八極の薄いキックのヤコビアンによるチューン (1e-4) は確認。有限振幅の周期解の条件そのもの、振幅依存チューンシフトを目標にした FIT の収束、結合・六極のあるラティスの fam 点は未確認 (追跡の NAFF か六極の解析的な振幅依存チューンシフトが要る) | 03_matching/FFS_FITPOINTS |
| `CSROYMatrix`・`OYStability`・`OYScan` | 両 | `CSROYMatrix` の行列そのものは、手で作った入力 (nr=64、nj=3、nphi=2) で gs と 4e-20 で一致。gs で `CSRSetupOY` の出力から作ると wrongtype (原因は特定: `tfcsroysetup` が zj の実数リストを head から書く。B60)。`OYStability`・`OYScan` は gs が完走せず比較なし。pysad は両モードで `CSRSetupOY` の出力を受ける (k64 で gs の wrongtype を再現するかは要判断)。nj・nphi が大きい手作りの入力は未確認 | 06_collective/fn_CSRInit_CSRTrack、v2_next/D2_csroymatrix.sad |
| `tbrote` の `srot` 更新、`tbendebody` の `bsi2` | k64 | ソースの読みだけで、gs の個別確認なし。v2 は行列積・入口と出口の縦磁場に直した (Emittance[] と追跡で 1e-5 以下、Thomas-BMT との比較) | |
| `measure.py` の `tlum`・`tluma` | 両 | gs で MEASURE が動かない。解析値のみ。K64 の MEASURE・TRACK からは呼ばれず、`ttstat_lines(…, lum=True)` で出せる | |
| `tfloor.f` の `cfunc` (`xsincos` など) | 両 | 個別実装 | |
| `BeamPlot`・`ColumnPlot`・`LaTeX`・`PointSymbol` の見た目 | 両 | 動くことだけ確認。LaTeX は `$HaveLaTeX` が空の状態での一致 | 11_graphics |
| Tkinter の見た目・操作 | 両 | 色・大きさ・フォントは SAD と比べていない。Python 3.14 (Tk 9.0.4) では、ウィジェットの木・幾何・フォント・色を 3.12 (Tk 8.6.13) と全ウィジェットで比べ、標準出力は Tk の版の行以外同じ。Checkbutton・Radiobutton の幅が 2〜3 px 大きい。実際の画面の撮影とマウス操作は未確認 | 12_gui_tkinter |
| Python 3.14 (Tcl/Tk 9.0) | 両 | v2.1.0 で `test_tclbridge.py`・`test_plot.py`・`test_gui_lang.py` の 25 件と、`test_v2_D2_*`・`test_spin.py`・`test_corr.py`・`test_csr.py`・`test_fam_points.py` の 90 件が通る (1 skip)。全試験は 3.14 では走らせていない。3.10・3.13 は未確認 (3.10 は同じ `nan+nanj` の分岐なので通る見込み) | |
| 乱数の消費の NaN の扱い (B17) | 両 | 決定済み: 現状のまま。`asinz(NaN)`・`gfloor`・`sqrt1` の下限は、発散した粒子が落ちるまでの座標の埋め方にだけ効く。「NaN になった粒子をその場で落とす」に変えると落ちる要素・周が K64 と変わる (LOSSMAP の値)。マッチングの fuzz (tffsmatch) は未確認 | |
| k64 モードの最初の呼び出しのスピン追跡の損失 | k64 | K64 の B33 の癖 (落とす印 0 の粒子が後で入れ替えで戻る) は track6 の経路では近似 (gs は要素 60、pysad は 15)。高速追跡 (track6n) は全桁一致 | |
| MASS・CHARGE を与えた粒子の Touschek・Intrabeam・RADINT | 両 | `rclassic = charge²·elradi` は電子の値のまま (B112)。陽子での確認なし。放射の β の因子も直していない (β < 0.99 の粒子の放射は U0 の 1/β² ほど大きく出る) | v2_next/A1_proton_ring.sad |
| `Twiss["DETR"]` ≥ 1 | 両 | モードが入れ替わる INS・MEASURE は未確認 (DETR < 1 のリングは R1 = 0.3 まで 0.24〜0.43 で入れ替わらなかった。入口を DETR > 1 の生の値で与える輸送路は gs と pysad で入口の規格化が違い比較できなかった) | v2_next/A1_n12_ins_detr.sad |

## 4. 既知の差 (未解決・または意図した差)

| 差 | モード | 状況 |
|---|---|---|
| SOL 内の四極 + RAD | k64 | gs が RAD ON・BZ≠0 の領域内で四極の力を入れていない疑いが強い (gs 側の欠陥。B68)。v2 の領域内の四極は物理 (Thomas-BMT・運動方程式の積分) に合わせた。現在の gs と bench_level9 の保存値が違う |
| `pm1.sad` の NaN (K3 が強くずれの大きい MULT リングの CALC の閉軌道が gs と別の解) | 両 | 原因未調査 |
| `LOSSAMPL`・RAD・FLUC つきの TRPT、MAIN の TRACK の損失報告、DAPERT での APERT | 両 | gs と比べていない |
| v2 の Emittance[] (RADCOD) と gs | v2 | 1e-4〜1e-3 ずれる。**意図した差** (tradke の符号の修正。B111。k64 は全桁一致)。大きなリング A (RADCOD) の εx +1.9e-4、ダンピング率 +2e-4、U0 +1e-7 |
| v2 の Emittance[] と RADINT (矩形磁石) | v2 | 1.4e-3 (大きなリング B)。v2 の修正の側の問題かもしれない (未調査) |
| PSPAC の gs との x と y の 16% の非対称 | v2 | 再現しない (v2.1.0。同じ記述で gs は x 1.2338・y 1.2356 mm で対称、v2 は 1.2270・1.2282 mm、解析式 1.2481 mm)。元の例の乱数のビームか PSPAC の状態の持ち越し (B96) による可能性が高い。粗いメッシュの差の原因はメッシュの分解能で、警告を足した |
| ASTRID2 の MPW・K0 の極の RADINT | v2 | εx が Emittance[] と 3.5%、Jx が 0.7% 違う (B115) |
| 追跡の FLUC の脱偏極 (v2 の Polarization の過小評価) | v2 | 4.2 GeV・12 MV で追跡の 0.43 倍 (§2) |

## 5. gs 側の問題で比べられないもの

| 項目 | 状況 |
|---|---|
| `For` | gs は 4 番目の引数 (本体) を実行しない (B71)。k64 モードの pysad は gs と同じ (本体を実行しない。以前の pysad は実行していた)、v2 モードは本体を実行する |
| Bessel・HermiteH・Ei などの一部 | gs が誤った値を返す (`BesselJ`(100) など。B43〜B47、B53、B54)。pysad は scipy の値 (真値)。`InverseErf` だけ k64 モードが K64 の値、v2 モードが真値 (B47) |
| 複素の `PolyLog`・`HypergeometricU`・`HurwitzZeta` | gs が誤る点がある (B118〜B120: 非整数次の PolyLog の Re z < 0・Im z ≤ 0・\|z\| > 1、HypergeometricU の b が整数、HurwitzZeta の複素の s で a の畳み方によって誤る)。pysad は両モードで真値 |
| gs が落ちる | `RadiationField` の SIGSEGV (B49。ほかの引数の渡し方でも落ちる)、`TRCOD` の SIGSEGV と MEASURE の無出力・MAIN の TRACK (B42)・`ON FOURIER` の TRACK の SIGSEGV (`tinip`)、`TCAVI` の L≠0 (B36)、PSPAC の既定 (`PSPACNZ=1`、B93)、`SPAC` で粒子数が少ない (2〜33 個などで SIGSEGV・NaN。**原因は粒子数が少ないことで RADIUS の大きさではない**: spkick.f の spmesh の添字。B98)、BEAMBEAM が SOL 領域内にある配置 (B35)、`RenumberElement` (B51)、`BYE ALL` (SIGSEGV)。**`SynchroBeta` の結合ラティスで gs が落ちるという以前の記述は再現しなかった** (24 セル FODO + skew 四極で gs は落ちず k64 と 7 桁一致。SOL・強い結合で落ちるかは未確認) |
| gs が誤った値を返す | SPAC の力が 1/γ 倍に弱い (B89。解析式の K/γ に近い)、領域内の BEND の偏向をほぼ無視する (B69)、SOL (BZ≠0) 領域内の四極 + RAD (B68)、WSPAC のチューンシフトが 3〜4% 強い (B127)、WSPAC ON の CALC・CALC6D が NaN (B133)、`Polarization` (B55・B125)、光子の偏光 ξ1 が ±1 を超える (B123)、`Union` の `SameTest` (B121)、`CSRSetupOY` の出力で `CSROYMatrix` が wrongtype (B60) |
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

## 7. 優先度 ([16_fortran_status.md](16_fortran_status.md) §9.4 から)

| 優先 | 項目 | 理由 |
|---|---|---|
| 中 | 誤差入り 6D emit と DA の反復 (性能) | 反復が 600 s を超える |
| 中 | `Polarization` の追跡との差 (脱偏極率の過小評価) | 非線形の軌道の効果 (2σ 以上の振幅の応答) を式に入れる。多数の試料 (数百) で裾の応答を測る必要 |
| 中 | v2 の Emittance[] と RADINT の矩形磁石の 1.4e-3 の差、v2 の Emittance[] (RADCOD) と gs の 1e-4〜1e-3 の差 | v2 の修正の側の問題か確かめる (§4) |
| 中 | `SynchroBetaEmittance` の大きなリングの不一致 (B128) | K64 の式の導出 (temits) か長時間の追跡との比較 |
| 中 | TrackParticles の F1 (`tsteer`・`tmulti`)、MULT の AE1・AE2・`FB1/FB2`、SOL 領域内の光子・領域内の四極のずれの POL (track6) | K64 のソースの直訳 |
| 中 | v2 の要判断 (docs/47 v2.1.0 の (d)) | 複素の特殊関数の k64 の値、tradke の符号・RADTAPER の重み、RADINT の空洞の位相、UND の光学、`DAPWIDTH` の窓 ほか |
| 低 | `wfres` 系 (MAIN の FOURIER) | 必要になったら。NAFF で代用 |
| 不要 | MON/MCOD/TCOD の機能、Obsolete・非リンクのファイル、`tdrwdt.f` の粒子追跡 | K64 で無効 |

## 8. 注意

- 条件: 例の数値は、サンプルのラティス (`example/00_lattice/FODO_ring.sad`: 3 GeV、周長 105.6 m) と、PF・ASTRID2 のラティス (`example/pf`、`example/astrid`) で得たものです。v2 の機能の数値は、各記録の条件 (上の冒頭) の範囲です。
- 表示や数字の桁の微妙な差は、ここには入れていません。完走することと、意味のある値が合うことを基準にしています。
- 一覧は変わります。項目を直したら、この文書と [16_fortran_status.md](16_fortran_status.md) §9 を一緒に更新してください。
