# 08_tracking.md から取り除いた履歴の記述

08_tracking.md (現在の仕様の解説) から、版・経過に結びついた記述を取り除いた。下は取り除く前の元の行を、節の見出しつきでそのまま並べたもの (行番号は整理前のもの)。

## # 08. TrackParticles — SAD の tturn 系粒子追跡と pysad の track6.py / track6n.py

整理前の 3 行目:

````text
最終更新: 2026-11-18 (v2.1.0 の後) (現状のスナップショット)。SAD ソースの行番号は 2026-10 時点。
````

## ### 1.2 NPARA > 1 での fork (tftrack.f:159-196)

整理前の 45 行目:

````text
- 子プロセスごとに `kseed = kseed + 2` として `tfaddseed(kseed)` を呼び、`SeedRandom[現在の種 + kseed]` で乱数の種をずらす (tftrack.f:629-)。したがって FLUC の乱数列を SAD と pysad で一致させるには **`NPARA = 1` をスクリプトに入れる**か、**`NPARA=n;` で同じ分け方を再現する** (v2.1.0。pysad の既定は NPARA=1)。子 k (k = 1, 2, …、条件 `ipn + np1 < np` かつ `npr < npara`) は乱数の種 = 親の種 + 2k で粒子 `[ipn, ipn + np1)` を追い、親は残りを元の種で追う (有効な並列数 = min(NPARA, np/9 + 1, ne·np/3000 + 1))。pysad の分け方は RAD + FLUC の追跡だけ (POL・RADLIGHT・PHOTONS・WAKE は分けない)。FODO リング (np = 200、1・10・100 周) の rms dp と最終の乱数の種が NPARA=1 と 4 の両方で gs と全桁一致、np = 16 (有効 npara = 2) の粒子ごとの dp も 12 桁一致、PF リング (760 要素、6 極つき、np = 40) も一致 (k64 モード。v2 は放射の式の修正で rms が 6e-6〜1e-3 違うが、乱数の引き方・最終の種は同じ。例 [`example/v2_next/C3_fluc_nparallel.sad`](../example/v2_next/C3_fluc_nparallel.sad))。
````

## ### 9.7 v2 モードで変わった追跡 (k64 モードは K64 のまま)

整理前の 367 行目:

````text
- **`MASS`・`CHARGE`** (v2.1.0。pysad の不具合の修正で、両モード): `TrackParticles`・Emittance[]・CALC6D に渡る (以前は電子の質量・電荷 1 のまま)。FFS の中で与える (`CHARGE` は FFS に入るたびに 1 に戻る: gs と同じ)。運動学 (ドリフトの z = L (β1/β0 − 1)、1e-12)、ソレノイドのラーモア回転 (電荷 1・−1・2、5e-12 m)、Qs (30 GeV/c の陽子、1e-4) が解析式と合い、k64 は gs と全桁一致 (陽子の小リング、大きなリング A の 3 GeV/c の陽子)。放射の古典半径は v2 だけ質量に比例 (B112)。使えない・未確認: 放射の β の因子 (β < 0.99 の粒子の放射は U0 の 1/β² 倍)、TouschekLifetime・Intrabeam・RADINT の古典半径 (電子の値のまま)、陽子での確認は 3 GeV/c・30 GeV/c のリングだけ。例 [`example/v2_next/A1_proton_ring.sad`](../example/v2_next/A1_proton_ring.sad)。
````

整理前の 368 行目:

````text
- **光子 (PhotonList。RAD・FLUC・PHOTONS ON)** (v2.1.0): `tdusrnpl`・`radangle1`・`tphrec`・`tphotonlist` の直訳。行 = [E (eV)、gx gy gz、ux uy uz、xi1 xi2 xi3、粒子番号、要素番号]。対応: BEND (K1=0: `tbend` の刻みごと、K1≠0: `tbendi` の刻みごと + 最後。ROTATE・DX・DY・DROTATE)、QUAD、MULT (K0・K1・K2・DX・DY・ROTATE)、SEXT・OCT (DX・ROTATE)、磁場つきドリフト、BEAMBEAM + BEND (TRPT つきでも)。k64 は gs と全列一致 (光子 358 個: BEND と QUAD、325 個: MULT・QUAD・BEND、71 個: BEAMBEAM + BEND + TRPT と追跡後の乱数の状態。E 列は 1.4e-8 違う: K64 が γ でなく βγ を掛けるため。B132)。**v2 は ξ1 の誤りを直す** (B123: |ξ1| ≤ 1) し、E 列に γ を掛ける。解析式 (v2): 光子数 (5/(2√3)) α γ θ が統計の 4σ 以内 (24 セル FODO・3 GeV・100 粒子で −0.10%)、⟨u⟩/u_c = 8/(15√3) = 0.3079 (0.3059 ± 0.0023)、⟨u²⟩/u_c² = 11/27、スペクトルの KS 検定、位置が円弧の上 (2 mm 以内)、SEXT・OCT・MULT の光子数。BEND・MULT・SEXT・OCT・QUAD・BEAMBEAM があるラティスは Python の追跡 (`track6`) になり遅い (24 セルで 1 粒子 1 周 1 ms 程度)。使えない・未確認: **SOL 領域内の要素・ドリフト (磁場なし) の光子は記録しない** (乱数だけ K64 と同じに引く)、COORD・TCAVI・MAP があるラティス (従来の高速追跡で BEND があっても光子なし)、POL の追跡との併用、QUAD の光子の位置は gs と未照合 (BEND は全桁)、強い6 極・8 極で不安定になる条件。track6 の `_ndivrad` が TRPT の 10 倍を見ていなかったのを直した (POL・RADPOL・光子の追跡の TRPT の結果が gs に近づく)。例 [`example/v2_next/B3_photons_bend.sad`](../example/v2_next/B3_photons_bend.sad)・`D1_photons_mult.sad`・`D1_photons_bb_bend.sad`。
````

整理前の 369 行目:

````text
- **BEND の F1・FRINGE の RAD の追跡** (v2.1.0。pysad の不具合の修正。B113): K64 の tfgeo の inext・iprev を再現した (`track6n.bend_links`)。先頭が長さ 0 の MARK の輸送路の最後の BEND は自分自身に繋がり、F1 の放射の断片が無い。以前の pysad は断片を常に入れていた。gs と 1e-6 (8 桁表示): 1 個の BEND (INS) の 1 周 (δ = −2.17142054679e-5)、隣り合う BEND の 12 セルのリング、K1≠0 の BEND の対 (`tbendi`)、POL (Python の追跡)。大きなリング A・PF・ASTRID2 は F1 なしで一致。
````

整理前の 370 行目:

````text
- **CAVI の DX・DY・ROTATE と `V1 V20 V11 V02`** (v2.1.0。pysad の不具合の修正): 追跡の CAVI に K64 の TENT・TEXIT を入れ、パーサが V1・V20・V11・V02 を捨てていたのを直した。L≠0 の CAVI + DX・DY の差 1e-3 が解消 (gs と全桁。大きなリング A で空洞に DX = 2e-4・DY = −3e-4、QDA2 を 0.03 rad 回して RAD の 5 周追跡が差 0)。使えない: WakeFunction のある CAVI・SOL 領域内の CAVI の V1 V20 V11 V02 (NotImplementedError)。
````

整理前の 372 行目:

````text
- **POL + TRPT + 加速** (v2.1.0。v2 だけ): SOL 領域の外の VOLT≠0 の CAVI (DX・DY・ROTATE も) を track6 で通す。断熱加速の解析式と 5e-7 (1 GeV、扇形の偏向磁石の間に 0・0.2・0.5 GeV)。k64 は従来どおり明示のエラー。使えない: SOL 領域内の加速、RF つき MULT・TCAVI の加速、リングの何周もの TRPT (K64 の TRPT は輸送路用)、スピン追跡の損失の処理。
````

整理前の 373 行目:

````text
- **PSPAC のメッシュの警告** (v2.1.0。v2 だけ): 横の rms が `PSPACDX`・`PSPACDY` の 1.5 セル未満 (キックが 10% 以上小さい)、縦の rms が `PSPACDZ` の 1.5 セル未満 (`PSPACNZ ≥ 4` のとき)、粒子の 1% 以上がメッシュの外のとき `???-PSPAC-Warning-…` を出す (結果は変えない)。メッシュが十分細かいとき (a ≥ 4 DX、一様なバンチの lb ≥ 10 DZ または `PSPACNZ<4`) は解析式と 1〜2%、gs と 1%。粗いメッシュは gs と違う (K64 は電荷を横 CIC・縦「ほぼ最近接」で置き、v2 は 2 次の B スプライン)。gs の x と y の 16% の非対称は再現しなかった。gs は PSPAC の状態を C の static 変数に持つので、1 セッションで TrackParticles を繰り返すと 2 回目から PSPAC が働かない (gs と比べるときは 1 セッション 1 回)。例 [`example/v2_next/A1_pspac_resolution.sad`](../example/v2_next/A1_pspac_resolution.sad)。
````

整理前の 386 行目:

````text
設定変数: `PBUNCH` (既定 1e10)、`SPACSTEP`、`SPACNX/NY/NZ`、`PSPACNX/NY/DX/DY/NZ/DZ/NTURNCALC`、`WSPACBEAM`、`WSPACRADIUS` (詳細は [`docs/history/v2_changes_C_spch.md`](history/v2_changes_C_spch.md)・`v2_changes_C2_spch.md`)。フラグは FFS の最上位に書く (`Do[...]` の中では切り替わらない)。粒子には z の広がりが要る。`RFSW` を ON にする (`NORFSW` は区間ごとに z を 0 にするので v2 は空間電荷をエラーにする。`RADIUS<0` の絞りだけなら `NORFSW` でも走る)。例 `example/05_tracking/flag_SPAC` (作り直した。以前の「gs でも SPAC の効果が出ない」は `NOSPAC;`・`SPAC;` を `Do[...]` の中に書いてフラグが切り替わっていなかったため)、[`example/v2_C/C2s_mult_radius.sad`](../example/v2_C/C2s_mult_radius.sad)・`C2s_wspac.sad`。
````

整理前の 389 行目:

````text
- 使えない・未確認: **追跡の WSPAC のビーム行列は空間電荷なし・MINCOUP なしの平衡値** (Emittance[]・EMIT・CALC6D・CALC (RFSW ON) には空間電荷が入る: docs/07。ビームサイズ `LINE["SIG…"]` は従来どおり)、ANGLE≠0・VOLT≠0 の MULT と SOL 領域内の MULT の SPAC (メッセージを出して空間電荷なしで追う)、QUAD・BEND・SEXT の中の SPAC (偏向磁石のあるリングでは SPAC は全周には入らない。全周なら PSPAC + SPCH か WSPAC)、複数バンチ (`NBUNCH`)、壁の抵抗、円でないダクト、MAIN の TRACK・`DynamicApertureSurvey`・スピン (`POL`) の追跡、ウェーク・放射との組み合わせ (走るが独立な根拠なし)、バンチの静止系で相対論的な運動 (γ·x' ≳ 1)。電子以外の粒子: 式は 100 MeV の陽子で 0.5% 確認した。`TrackParticles` は `MASS`・`CHARGE` を追跡に渡す (v2.1.0。SAD スクリプトからの空間電荷の使用は未確認)。マクロ粒子は 1000 個以上を勧める (チューンシフトを数 % で見るなら 4000 個以上)。
````

