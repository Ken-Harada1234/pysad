# 21_sad_bug_report から取り除いた履歴の記述

`docs/21_sad_bug_report.md` を現在のスナップショットだけの解説にするため取り除いた文章。元の文脈 (節の見出し) つきで、文章はそのまま。

## ## 0. 一覧 (重要度順)

取り除いた部分: 「**v2.2.0 の更新**」 (置き換え後: 「**追記**」)

元の行 (整理前の行 60):

> | B42 | `TRCOD` が gs (この環境の K64 ビルド) で SIGSEGV。`MEASURE` は何も出力しない | trcoda.f / trackb.f | TRCOD・MEASURE の参照が取れない | 低 | 確定 (gs で再現) | v2 で MEASURE・TRCOD・LOGMEASURE・MAIN の `TRACK TURNS=-1` を実装 (ソースの読み + 解析式・pysad の独立な経路で確認。gs との比較は無い)。k64 モードは未対応 (無視) **v2.2.0 の更新**: gs の 10-19 ビルドは MEASURE を出力する (先頭の要素 = 入口の MARK の位置だけ。それ以外は無出力、TRCOD は SIGSEGV)。k64 は先頭の位置の MEASURE・TRCOD・LOGMEASURE を gs と同じ統計 (4 桁) で実装した (B151・B152・B154)。先頭以外は無出力・警告。v2 の追跡を通した値は gs と比べられない |

## ## 0. 一覧 (重要度順)

取り除いた部分: 「**v2.2.0 の更新**」 (置き換え後: 「**追記**」)

元の行 (整理前の行 129):

> | B108 | ANGLE=0 の BEND (tsteee・tsteer、K0 だけの補正磁石) と K0 つき MULT の端の角は、B75 と同じ「端の蹴りを偏向角に入れるが行路は変えない」作りのはず | tsteee・tsteer 等 | 端の角つきの補正磁石 (AE1・AE2 つきの挿入光源の極) の放射 | 中 | 確定 (硬い端の一様磁場の極の厳密な幾何 (円弧と直線) とロビンソンの定理 Je = 2 + dln(Σ B_i² l_i(δ))/dδ。24 セルの FODO に極 3 個の組 (K0・−2K0・K0) を入れ、AE = (0.3, 0.3) で ΔJx = 厳密 +0.01440、v2 +0.01487 (1.032)、K64 +0.02992 (2.08 倍)。(−0.2, 0.4) で 1.037 と 2.16 倍。K64 は端の項が約 2 倍 = B75 と同じ誤り) | v2 で修正 (`_tsteee` の `tradke` に tin = tan(AE1)・tout = tan(AE2) を渡す。ANGLE≠0 の `tbende` と同じ)。k64 モードは K64 のまま。残る 3〜4% は EMITDIV = 20 の有限の分割と 2 次の項 (傾き 30 mrad の極では 2〜3%)。K1≠0・K0≠0 の ANGLE=0 BEND と MULT の端の角は pysad が AE を読まず未確認 **v2.2.0 の更新**: MULT・K1≠0・K0≠0 の ANGLE=0 BEND の端の角 AE1・AE2 を実装し (CALC・Emittance[]・追跡。k64 は gs と全桁一致)、放射の行路も v2 で修正した (B143)。追跡 (`tsteer`) の端の放射も v2 で修正 (ΔJx は厳密な幾何の 1.03〜1.04 倍)。傾き 30 mrad の極の 0.02〜0.03 は K0 に依らない 2.7% の系統差 (AE に比例、極の長さ・間隔に依る。原因未特定) |

## ## 0. 一覧 (重要度順)

取り除いた部分: 「**v2.2.0 の更新**」 (置き換え後: 「**追記**」)

元の行 (整理前の行 132):

> | B111 | 放射の反作用の角度の項の符号が tradke (Emittance[]・閉軌道: `ddpx = +0.5·dpx·dg`) と tradkn (追跡: `−0.5·dpx·dg`、dg < 0) で逆。エネルギーを失う間は同じ磁場が運動量の小さい粒子をより強く曲げる (BEND の中の平均 δ = dg/2 → 偏向角 +\|dg\|/2 倍) ので追跡の符号が物理的に正しい | photons.f (tradke:1417、tradkn:1247) | RADCOD ON の閉軌道が追跡で閉じない (3 GeV・周長 105.6 m の FODO で 1 周 x −4.3e-7 m、z 2.4e-6 m。RADTAPER ありで 2e-6・1e-5 m)。閉軌道・U0・チューン・ダンピング・エミッタンスが相対 1e-5〜1e-3 動く (大きなリング A (3 GeV、603 m): εx +1.9e-4、ダンピング率 +2e-4) | 中 | 確定 (gs で再現。符号は解析的な導出と、反転後に閉軌道が 1 周で 1e-9 m で閉じることによる。作者の意図の裏づけは無い) | v2 で修正 (tradke の ddpx・ddpy と行列 D・B の微分の符号を反転。RADTAPER の追跡の BEND の重みも tturne と同じ (3δ_l + δ_l+1)/4 にした: これは「閉じる」ことを根拠にした選択で、どちらの重みが物理的かは未判断)。k64 モードは再現。v2 の Emittance[] (RADCOD) は gs と相対 1e-4〜1e-3 ずれる **v2.2.0 の更新**: 符号は確定 (放射反作用つきのローレンツ力の数値積分と v2 の tradke が 6 桁一致。k64 は 1 台の偏向磁石・既定の EMITDIV = 1 で Δpx が 9% 違う)。RADTAPER の BEND の係数は (δ_l + δ_l+1)/2 が物理的なので、v2 は Emittance[] も追跡もこれにした (3/4 は K64 の粗い離散化に合わせた係数)。v2 の Emittance[] が gs とずれるのは K64 の誤差。微分の誤りは B138 |

## ## 0. 一覧 (重要度順)

取り除いた部分: 「**v2.2.0 の更新**」 (置き換え後: 「**追記**」)

元の行 (整理前の行 133):

> | B112 | 放射の古典半径が粒子の質量によらない (`re0 = elradi`、電子の値。`rclassic = charge²·re0`)。追跡の `tradkn` は光子数 (5/(2√3)) α γ θ と臨界エネルギーから δ を下げる超相対論的な式で、β の因子が無い (U0_K64 = U0/β³) | tfsetparam.f、photons.f (tradkn) | MASS≠電子の粒子の放射: 30 GeV/c の陽子 (γ = 32) の U0 が古典式 6.88e-4 eV に対し K64 は 1.2651 eV (mp/me = 1836 倍)。β の因子は γ = 32 で +0.15% (β < 0.99 では 1/β² ほど大きい) | 中 | 確定 (K64 の値は gs と全桁一致。古典式 U0 = (2/3) r_p mc² β³γ⁴ Σθ²/L と v2 の追跡の 1 周の δ の変化が 1e-6 で一致) | v2 で修正 (`re0 = elradi·elmass/amass`)。β の因子は直していない (γ ≫ 1 の粒子だけ意味を持つ)。k64 モードは再現。TouschekLifetime・Intrabeam・RADINT の `rclassic` は K64 のまま (陽子での確認なし) **v2.2.0 の更新**: β の因子は直していない (要判断)。k64 は gs と全桁 (1.5・3・30 GeV/c の陽子の U0)。v2 の U0 は古典式のちょうど 1/β³ 倍 (1.5・3・30 GeV/c で 64%・15%・0.15%。γ ≫ 1 では 3/(2γ²)) |

## ## 0. 一覧 (重要度順)

取り除いた部分: 「**v2.2.0 の更新**」 (置き換え後: 「**追記**」)

元の行 (整理前の行 137):

> | B116 | `Polarization`・`PolarizationTime` は、逆向きの磁場 (挿入光源の交互の極) があっても符号つきの率 sdamp = n·gintd = Σ r_i (b_i·n) を緩和率に使う。逆向きの極が相殺する | temit の偏極の量 | ASTRID2 + MPW で P = 92.376% (一定。独立な式は 70.16%)、τ = 28842 s (独立な式 21886 s。実際より 31% 長い)。K64 の RADPOL 追跡の緩和は Σ r (b·n)² の側 | 中 | 確定 (独立な式 P = (8/5√3)⟨ρ⁻³⟩/⟨\|ρ\|⁻³⟩、1/τ = (5√3/8) r_e ƛ_C c γ⁵⟨\|ρ\|⁻³⟩ と、K64 の追跡の緩和率の比 (MPW あり/なし 1.311) で確認) | v2 で修正 (緩和率に符号に依らない Σ r_i (b_i·n)² (`sdamp_rel`) を使う。P = 0.70190、τ = 21915 s (3.9e-4、1.3e-3 の差))。k64 モードは再現 (P は NaN: B55)。n が傾く・非平面のリングで Σ r (b·n)² が最適かは未確認 **v2.2.0 の更新**: n が傾く・非平面のリングは厳密式 Σ r (1 − 2/9 (n·v)²) を v2 に入れた (`spin.ST_EXACT`。垂直の偏向・傾いた n・水平磁場の MPW で独立な Derbenev-Kondratenko の数値積分と τ が 1e-5〜2e-3)。分子の運動論的偏極 (`spin.KINETIC`)、脱偏極率の放射の刻みごとの積分 (`spin.STEP_RATE`) も入れ、P は 0.1〜0.6% (脱偏極率が分母の 15% 以上では 2〜3%) で合う。A2 の Σ r (b·n)² は n が b に垂直で 0 になる (水平磁場の MPW で τ が 32% 短い) |

## ## 0. 一覧 (重要度順)

取り除いた部分: 「**v2.2.0 の更新**」 (置き換え後: 「**追記**」)

元の行 (整理前の行 139):

> | B118 | `HurwitzZeta[s,a]` の s が複素で、a を 1 以下に畳んだ値が 0.5 以上 1 未満 (a = 0.5・1.5・3.5 など。0.25・1・2・10 は正しい) のとき誤る。`Zeta[s]` の臨界線上の零点は 1e-16 の雑音 | gamma.f | 複素の s の HurwitzZeta: 相対 1e-6〜1e9 の誤り (例 `HurwitzZeta[2+I,.5]` = 3.9549+2.1791 I、真値 3.5075+2.0314 I)。Zeta 40 点・HurwitzZeta の他の点は 1e-13 以内で正しい | 低 | 確定 (別に書いた Euler–Maclaurin と 1e-10、Zeta は Borwein のアルゴリズム・関数等式と 1e-11) | v2 は真値を返す (`specfn_cx.py`)。**k64 モードも真値を返す** (gs の誤った値は再現しない: 再現には複素の不完全ガンマなどの直訳が要る。要判断)。gs の値が必要なら gs を使う **v2.2.0 の更新**: 複素の特殊関数の gs の誤りの全体は B158 (732 式の走査。追加した 15 関数で 641 式が gs と 1e-9) |

## ## 0. 一覧 (重要度順)

取り除いた部分: 「**v2.2.0 の更新 (原因を特定して v2 で修正)**」 (置き換え後: 「**追記 (原因を特定して v2 で修正)**」)

元の行 (整理前の行 149):

> | B128 | `SynchroBetaEmittance` が大きなリング A (3 GeV、周長 603 m、νs = 0.0041、RF 1.4 MV) で 6.18e-9 を返し、Emittance[] の 1.55e-10 と 40 倍違う。AzimuthalModes 3〜8 で変わらず、12・20・30・32 以上で 3.19e-9・1.46e-9・8.50e-10・7.87e-10 と動き続ける。νs の入力・エネルギー・6 極の強さにも依らない | SynchroBeta (temits) | 大きなリングで使えない (独立な 6x6 の平衡 1.5173e-10 と、RAD + FLUC の 6 次元追跡 (60 粒子、4000 周) の 1.6e-10 と合わない)。PF・ASTRID2・SP8 (大きなリング B、8 GeV) は Emittance[] と 0.5〜1% で合う | 中 | 確定 (v2.2.0 で原因を特定: 右端の欄の更新。v2.1.0 までは「疑い」で、gs も pysad k64 も 6.18e-9 だった) | v2 は値を変えず、AzimuthalModes m と 1.5 倍の m の emitx が 20% 以上違えば警告を出す。k64 モードは gs と同じ (警告なし)。使うときは AzimuthalModes を変えて収束を確認し Emittance[] と比べる。`AzimuthalModes->1` は pysad が IndexError (gs は 8.6e-7 を返すが意味がない) **v2.2.0 の更新 (原因を特定して v2 で修正)**: K64 の `tefsetup` が 1 周の写像 X' = T X + h1 の源泉 h1 に閉軌道 c(δ) (1 次の分散 η δ を含む) をそのまま使う (固定点が c にならず、2 次元の回転では平均が c の半分になる)、`tesumb` が x_b の 2 次のモーメントに η η σδ² だけを足し相関 η σδ² (`beam(16:19)`) を入れない。リングの始点の分散 η ≠ 0 で εx が始点に依り数倍〜40 倍過大 (24 セル FODO の始点を QF・BEND・BEND の途中にして 6.54・6.54・4.74 倍、gs でも同じ。η が小さい PF・ASTRID2 は 1%)。確度は確定。v2 は源泉から 1 次の分散を引き、軌道の δ 依存の効果を解析的な 2 項に分けて入れる (`tesumb_v2`): Emittance[] の 1.0003〜1.0005 (FODO)・1.0084 (大きなリング A、AzimuthalModes 3。32 で 1.0011)・1.0045 (SP8 (大きなリング B))・1.005 (PF)・1.006 (ASTRID2)。k64 は K64 のまま (gs と 2e-7)。警告 (上記) は v2 の大きなリング A では出なくなった |

## ## 0. 一覧 (重要度順)

取り除いた部分: 「**v2.2.0**:」 (置き換え後: 「**追記**:」)

元の行 (整理前の行 150):

> | B129 | `CorrectOrbit` の `Solver->{"Micado",n}` (`Condition` なし) が、選んだ nx 本の補正量でなく、全補正量の最小二乗解の先頭 nx 個を使う。ハウスホルダー変換のあとの `call tsvd(a,b,x,eps,svd)` で、新しい tsvd は配列の大きさを `size()` で取る (コメントアウトされた旧 `tsolva(a,b,x,nx,nx,mdim,eps)` は nx を渡していた) ので全 n 列を解き、x(nx+1:n) = 0 にする | pmicad.f | 24 セル FODO (BPM 48 台、垂直ステアラ 24 台) の補正後の垂直 rms が補正前 9.41e-5 m に対し Micado n=1 で 1.46e-4 (悪化)、n=3 8.2e-5、n=6 1.19e-4、n=8 1.75e-4 (本数を増やしても単調に下がらない) | 中 | 確定 (gs で再現。選んだ nx 列の最小二乗 `lstsq` と v2 が 1e-9、選び方は独立な貪欲法と同じ) | v2 で修正 (上三角の nx x nx だけを解く。n=1 4.74e-5 … 24 本 2.06e-5 と単調に下がる)。k64 モードは gs と 1e-6 以内で再現。`Condition` つき (CorrectLocal) の Micado は選んだ列だけで `psvdcond` を解く元の分岐で K64 でも正しい。水平・2 平面同時・`Bump`・`SVDmode->1` つきは未確認 **v2.2.0**: 水平・2 平面同時の Micado も gs の k64 と 5e-6 で一致し、v2 は単調に下がる。`SVDmode->1` は B157 |

## ## 0. 一覧 (重要度順)

取り除いた部分: 「(以前から。モードによらない)」 (置き換え後: 「(モードによらない)」)

元の行 (整理前の行 61):

> | B43 | `Hypergeometric2F1[.5,.5,1,.99]` が 2.348414 (真値 2.352716)。`BesselJ[.5,30.]` が -.143924709693 (真値 -.143929653370)。級数・漸近展開の打ち切りが粗い | gamma.f / hg.f (hgrr, cbesj = confhg0) | 引数が 1 に近い 2F1、x が大きい Bessel で相対 1e-3〜3e-5 の誤差 | 低 | 確定 (gs で再現) | pysad は scipy の値 (真値) を返す (以前から。モードによらない) |

## ## 0. 一覧 (重要度順)

取り除いた部分: 「(以前から。モードによらない)」 (置き換え後: 「(モードによらない)」)

元の行 (整理前の行 66):

> | B48 | `FixedPointList` は gs で未評価のまま返る (autoload の表に `FixedPoinstList` と綴り違い (Misc.n の SetAttributes と init.n)) | Packages/init.n:597、Misc.n:23 | FixedPointList が使えない | 低 | 確定 (gs で再現) | pysad は Misc.n の定義どおり動かす (以前から。モードによらない) |

## ## 0. 一覧 (重要度順)

取り除いた部分: 「k64 は従来のまま」 (置き換え後: 「k64 は K64 のまま」)

元の行 (整理前の行 195):

> | B186 | `CoshIntegral[x]` (x < 0) は実部だけ (真値 Chi(|x|) + iπ)。`SinhIntegral`・`CoshIntegral` に複素数を渡すと虚部を捨てた実数 (pysad) | D.n・pysad | x < 0、複素引数 | 低 | 確定 (mpmath) | v2 は複素数に対応 (Shi(z) = −i Si(iz)、Chi(z) = Ci(iz) − iπ/2、Re z < 0 かつ Im z > 0 は + 3πi/2)。k64 は従来のまま。**gs 照合: gs は正しい** (`CoshIntegral[-2.5]` = 3.5244+πi、複素の Shi・Chi も 1e-15 で一致) **= pysad の誤り。k64 も直した** |

## ## 0. 一覧 (重要度順)

取り除いた部分: 「k64 は従来のまま」 (置き換え後: 「k64 は K64 のまま」)

元の行 (整理前の行 196):

> | B187 | `LogIntegral[z]` (z < 0 の実数) が NaN、複素数はエラー (pysad)。真値 li(z) = Ei(Log z) | D.n・pysad | z < 0・複素の LogIntegral | 低 | 確定 (mpmath) | v2 は Ei(Log z) (複素)。k64 は従来のまま。**gs 照合: gs は正しい** (`LogIntegral[-2.]`・複素) **= pysad の誤り。k64 も直した** |

## ## 0. 一覧 (重要度順)

取り除いた部分: 「k64 は従来のまま」 (置き換え後: 「k64 は K64 のまま」)

元の行 (整理前の行 200):

> | B191 | `HypergeometricU[a,b,z]` (実数) の scipy の値は z = 10 付近で相対 7e-11 の誤差 (`HypergeometricU[.5,.5,10.]`) | scipy | z ≳ 10 | 低 | 確定 (mpmath) | v2 は a > 0 を積分表示 (`specfn_cx.hyperu_c`、1e-13)。k64 は従来のまま。**gs 照合: gs は z = 10 で正しい (`HypergeometricU[.5,.5,10.]` が全桁)、z = 30 は gs が誤り** (`HypergeometricU[2.,2.5,30.]` = .001111、真値 .0010766)。pysad の scipy の 7e-11 は小さい差。v2 のみ |

## ## 0. 一覧 (重要度順)

取り除いた部分: 「k64 は従来のまま」 (置き換え後: 「k64 は K64 のまま」)

元の行 (整理前の行 201):

> | B192 | `HermiteH[n, z]` の複素の z がエラー (pysad)。`LegendreP[ν, m, x]` の非整数 ν・正の整数 m (−1 < x < 1) が NaN (正則化した 2F1 の Γ(1−m) の極。例 `LegendreP[2.5,1,.3]` = −0.1292) | D.n・pysad | 複素の Hermite、非整数 ν の Legendre | 低 | 確定 (mpmath、P_ν^m = (−1)^m (1−x²)^{m/2} d^m P_ν/dx^m の数値微分) | v2 は漸化式・(ν)_m の積を使う公式。k64 は従来のまま。**gs 照合: gs が誤り** (`HermiteH[3,1+I]` が NaN、`LegendreP[2.5,1,.3]`・`LegendreP[3.5,2,.3]` が 0)。k64 は NaN のまま |

## ## 0. 一覧 (重要度順)

取り除いた部分: 「k64 は従来のまま」 (置き換え後: 「k64 は K64 のまま」)

元の行 (整理前の行 204):

> | B195 | `LerchPhi[z,s,a]` の z < −1 が NaN、0.9 ≤ |z| < 1 は 20000 項の級数で精度が落ちる | pysad | |z| ≥ 0.9 の実数 z | 低 | 確定 (mpmath、関数等式 Φ = a^{−s} + zΦ(a+1)) | v2 は積分表示 (z < 1、s > 0、a > 0)。k64 は従来のまま。**gs 照合: gs は正しい** (`LerchPhi[-3,2,1]`・`LerchPhi[.95,2,1]` が 1e-16) **= pysad の誤り。k64 も直した** |

## ## 0. 一覧 (重要度順)

取り除いた部分: 「k64 は従来のまま」 (置き換え後: 「k64 は K64 のまま」)

元の行 (整理前の行 205):

> | B196 | `PolyLog[2,z]` が小さい z で spence(1−z) の桁落ちで相対 8e-8 (z = 1e-10)、`PolyLog[−n,z]` の 0.8 < |z| < 1 は相対 2e-11・|z| > 1 は NaN (`PolyLog[-2.,2.]` = −6)、`HarmonicNumber[x]` が小さい x で ψ(x+1)+γ の桁落ち (x = 1e-10 で相対 7e-8)、`Binomial[n, k]` の負の整数 k が NaN (真値 0) | pysad・scipy | 小さい z・x、負の整数次 | 低 | 確定 (mpmath、Li_{−n} の有理式) | v2 は級数・有理式 Li_{−n}(z) = N_n(z)/(1−z)^{n+1}・ζ(k) の級数・0。k64 は従来のまま。**gs 照合: `PolyLog[-3,-.9]`・`PolyLog[-2,2.]` と `Binomial[-5,-2]` = 0 は gs が正しい = pysad の誤り (k64 も直した)。gs が誤り: `PolyLog[2,1e-10]` が相対 2e-4、`HarmonicNumber[1e-10]` が相対 4e-6、`Binomial[3.5,-1.]` = .2222 (真値 0)**。この 3 つは v2 のみ |

## # 21. SAD K64 のバグ報告 (pysad 移植中に見つかったもの)

取り除いた部分: 「最終更新: 2026-11-22 (現状のスナップショット。B01〜B163・B220〜B227・B230〜B242・B290〜B292 (B290〜B292 は v2.4.1 の ANGLE=0 の K1≠0・K0≠0 の BEND・6 極のフリンジ・MULT の端の縦収束。B163 は v2.4.0 の放射・エミッタンスの作業、B220〜B227 は v2.4.0 D2 の端の幾何、B230〜B242 は v2.4.0 D3 のスピン・放射・SOL)。B72〜B110 は v2.0.0 の作業、B111〜B134 は v2.1.0 の作業、B135〜B159 は v2.2.0 の作業、B160・B161・B162 は v2.3.0 の作業 (T1・U1・T2) で見つけたもの (B160 は U1 で「誤りではなく K1 の暗黙の定義」に再分類、B161 は U2 で v2 で直した、B162 は Touschek の損失率の式の適用範囲で v2 は Piwinski 型の式を既定にした)。B40・B42・B60・B81・B83・B99・B106〜B109 は v2.1.0 で、B42・B108・B111・B112・B116・B118〜B120・B128・B129 は v2.2.0 で更新)」

元の行:

> 最終更新: 2026-11-22 (現状のスナップショット。B01〜B163・B220〜B227・B230〜B242・B290〜B292 (B290〜B292 は v2.4.1 の ANGLE=0 の K1≠0・K0≠0 の BEND・6 極のフリンジ・MULT の端の縦収束。B163 は v2.4.0 の放射・エミッタンスの作業、B220〜B227 は v2.4.0 D2 の端の幾何、B230〜B242 は v2.4.0 D3 のスピン・放射・SOL)。B72〜B110 は v2.0.0 の作業、B111〜B134 は v2.1.0 の作業、B135〜B159 は v2.2.0 の作業、B160・B161・B162 は v2.3.0 の作業 (T1・U1・T2) で見つけたもの (B160 は U1 で「誤りではなく K1 の暗黙の定義」に再分類、B161 は U2 で v2 で直した、B162 は Touschek の損失率の式の適用範囲で v2 は Piwinski 型の式を既定にした)。B40・B42・B60・B81・B83・B99・B106〜B109 は v2.1.0 で、B42・B108・B111・B112・B116・B118〜B120・B128・B129 は v2.2.0 で更新)

## # 21. SAD K64 のバグ報告 (pysad 移植中に見つかったもの)

取り除いた部分: 「> **v2.0.0 から**: pysad は互換モードを持つ。」

元の行:

> > **v2.0.0 から**: pysad は互換モードを持つ。既定の `v2` は K64 のバグを直した動作、`PYSAD_COMPAT=k64` (または `pysad.compat.set_mode('k64')`) は v1 と同じ K64 (gs) の動作。

## # 21. SAD K64 のバグ報告 (pysad 移植中に見つかったもの)

取り除いた部分: 「は v1 と同じ K64 (gs) の動作。」

元の行:

> > **互換モード**: pysad は互換モードを持つ。既定の `v2` は K64 のバグを直した動作、`PYSAD_COMPAT=k64` (または `pysad.compat.set_mode('k64')`) は v1 と同じ K64 (gs) の動作。

## # 21. SAD K64 のバグ報告 (pysad 移植中に見つかったもの)

取り除いた部分: 「> 一覧は [`47_v2_changes.md`](47_v2_changes.md)。B72〜B110 は v2.0.0 の作業の記録 (`history/v2_changes_*.md`)、B111〜B134 は v2.1.0、B135〜B159 は v2.2.0 の作業の記録 (`history/v2next_*.md`) の「新しく見つけたもの」を採番したもの。B160・B161 は v2.3.0 の T1・U1 (`history/v23_T1.md` §6、`history/v23_U1.md`)。」

元の行:

> > 一覧は [`47_v2_changes.md`](47_v2_changes.md)。B72〜B110 は v2.0.0 の作業の記録 (`history/v2_changes_*.md`)、B111〜B134 は v2.1.0、B135〜B159 は v2.2.0 の作業の記録 (`history/v2next_*.md`) の「新しく見つけたもの」を採番したもの。B160・B161 は v2.3.0 の T1・U1 (`history/v23_T1.md` §6、`history/v23_U1.md`)。

## ## 0. 一覧 (重要度順)

取り除いた部分: 「**v2.4.1 で影響の大きさを測った**」

元の行:

> | B58 | `tbendebody`・`tbendebody0` の `bsir0` の更新で `bsi2` の項は enarad が偽の最後のスライスにしか付かず、常に捨てられる (`bsi2=1` は n==n2 のときだけで、そのとき `enarad .and. n /= n2` は偽)。意図は偏向磁石の出口の bsi の補正とみられる | tbende.f:215・313、535-537 | 縦のオフセットがある偏向磁石の POL で出口側の bsi が入らない。**v2.4.1 で影響の大きさを測った**: 12 セル FODO (4.2 GeV、縦ずれ ±3 mm の四極、閉軌道 y 6 mm) の gs で Emittance[] の SpinTune が gs の追跡 (閉軌道上の 1 粒子の毎周のスピンの回転角) と 3.8e-4 違う (縦ずれなし 3e-5、1 mm で 1.4e-5。`history/v241_pol_fluc.md`) | 中 | ソースの読み + **gs の追跡との比較**: 出口の項だけを k64 の Emittance[] に足すと差が 3.9e-4 → 2.7e-5 に減る (pysad の emit6d を一時的に書き換えた試験。恒久の変更ではない) | v2 で修正 (入口で +(y/ρ + k x y)、出口で −。スピン軌道結合にも入れる。追跡の tbendi も同じ)。k64 モードは K64 と同じ (出口の項を捨てる) |

## ## 0. 一覧 (重要度順)

取り除いた部分: 「(手で作った正しい形の実数リストなら gs でも動く: v2.0.0 の記述「手で作った実数リストでも」は誤り)」

元の行:

> | B60 | `CSROYMatrix[{wi,J,omega,zj,dpj,rhoj},{dz,inten,nus,sigz,sige}]` が gs (この K64 ビルド) で、`CSRSetupOY` の出力を渡すと General::wrongtype を返す (手で作った正しい形の実数リストなら gs でも動く: v2.0.0 の記述「手で作った実数リストでも」は誤り)。原因は D2 の調査で特定: `tfcsroysetup` が 4 番目 (zj) の実数リストへ `c_f_pointer(c_loc(rlist(kazjn)), zjn, [mphi,nj])` で書くとき、先頭の語 `rlist(kazjn)` (リストの頭 = head) から書き始める (要素は `rlist(kazjn+1)` から書くべきで、ほかの 4 つの配列は `rlist(k+1:…)` で正しい)。head が zjn(1,1) で上書きされ、`tfcsroymat` の `tfreallistq(las%dbody(4), lazj)` が偽になる (B59 の `dpjn` が `zjn` と同じ記憶域も同じ箇所)。登録名は `CSRSetupOY`・`CSROYMatrix` (tfinitn.f:1044-1045。ソース内の関数名 tfcsroysetup・tfcsroymat と綴りが違う) | csroy.f:tfcsroysetup・tfcsroymat | `CSRSetupOY` の出力から `CSROYMatrix` を作れない。`OYStability`・`OYScan` は gs で完走しない | 低 | 確定 (gs で再現し、`CSRSetupOY` の出力の各要素の `Head` を調べて原因を特定: zj の Head が実数 1.2952814481588496e-5) | pysad の `CSROYMatrix` は両モードで動く (`CSRSetupOY` の出力も受ける: gs の wrongtype は再現しない。直していない: 要判断)。手で作った入力 (nr=64、nj=3、nphi=2) の 6x6 の行列 (36 要素) は gs と相対 4e-20 で一致 |

## ## 0. 一覧 (重要度順)

取り除いた部分: 「k64 モードは従来どおり TypeError」

元の行:

> | B61 | `FindRoot` に複素の初期値を与えると壊れた結果になる (`FindRoot[z^2+1==0,{z,1.+1.*I}]` が `{z->Restrict[(.25+I),-Infinity,Infinity],Residual->0}`、`Exp[z]==2` も同様) | tffindroot.f:tfsetupvars (352-435。初期値は実数だけ受け付ける作りで、:414-418 が複素に `General::wrongval` を出す経路。`Restrict[…]` が返る理由は未調査) | 複素の根が求まらない (実数の初期値なら正常: `FindRoot[x^2-2==0,{x,1.}]` は gs と pysad が 1e-14 で一致) | 低 | 確定 (gs で再現、2026-11-17。原因の特定はソースの読みのみ) | v2 は複素の初期値を複素の Newton 法で解く。k64 モードは従来どおり TypeError (gs の壊れた値は再現しない。実数の初期値は両モードとも gs と一致) |

## ## 0. 一覧 (重要度順)

取り除いた部分: 「。以前の pysad は S 座標だけに入れていた)」

元の行:

> | B63 | SOL (BOUND≠0) の長さが周長に入らない (omega0・MomentumCompaction・空洞の位相・偏極時間に効く)。SOL に GEO=1 がないと S 座標が途切れて周長がさらに狂う | 周長の計算 (場所は未調査) | 周回周波数・運動量圧縮・偏極時間が SOL の長さの分だけずれる。GEO なしでは偏極時間が約半分、MomentumCompaction が約 100 倍 | 中 | 確定 (gs で確認: PLANAR + SOL で、SOL を除いた周長 45.6 のときと偏極時間が 15 桁一致)。意図か誤りかは未判定 | バグではないと判断 (SOL に L のキーワードは無い: SADHelp)。pysad は SOL の L を軌道・光学・S 座標・周長に入れない (gs と同じ。モードによらない: `test_sol_has_no_length`。以前の pysad は S 座標だけに入れていた)。GEO なしの S 座標の途切れは再現しない |

## ## 0. 一覧 (重要度順)

取り除いた部分: 「(Break・Continue も効く)。3 引数は本体なしのループ。k64 モードは gs と同じ (本体を実行せず、3 引数は未評価。以前の pysad は本体を実行していた)」

元の行:

> | B71 | 組み込みの `For[ini,cond,incr,body]` は 4 番目の引数 (本体) を実行せず、3 番目 (incr) だけを繰り返す (`For[i=0,i<3,i++,Print[i]]` は何も印字せず i=3)。`Packages/functions.n` の定義 `For[`ini_,`cond_,`incr_,`body_]:=(ini;While[cond,body;incr];)` は本体を実行するはずで、同じ式を別名 (For2) で定義すると gs も実行する。3 引数の `For` は未評価のまま | `Packages/functions.n:34` (読み込み時の名前解決と思われるが原因は未特定) | `For` の本体が動かない (`Do`・`While` は正常) | 低 | 実測 (gs で確認。`local/wake/` の f1〜f5 相当) | v2 は定義どおり本体を実行する (Break・Continue も効く)。3 引数は本体なしのループ。k64 モードは gs と同じ (本体を実行せず、3 引数は未評価。以前の pysad は本体を実行していた) |

## ## 0. 一覧 (重要度順)

取り除いた部分: 「(CAVI だけ。行をまたいで書ける)。k64 モードは従来どおり |」

元の行:

> | B66 | 旧式の `WAKE` コマンドは `tfwake.f` を呼ばず `WAKECOMMAND[]` (init.n の AutoLoad) を評価する。`Packages/WAKECOMMAND.n` は配布に無く、fileopen のエラーを出して行の残りを捨てる | tffsa.f:976 (docs の旧記述は 972)、init.n:696 | WAKE コマンドが使えない (tfwake.f は死んだコード) | 低 | 確定 (gs で再現) | v2 は旧式の `WAKE L|T 要素 s1 w1 … ENDWAKE` を WakeFunction の表として登録 (CAVI だけ。行をまたいで書ける)。k64 モードは従来どおり |

## ## 0. 一覧 (重要度順)

取り除いた部分: 「k64 モードは PHSROT が従来どおり動かない」

元の行:

> | B79 | phsrot.f の拡散の乱数は v = U·√λ (ベクトル) を成分ごとに正規乱数へ掛けてから N を掛ける。U が対角でないと共分散が N B Nᵀ にならない | phsrot.f | PHSROT の拡散 (RAD かつ FLUC) の共分散。K64 では PHSROT が動かない (B38) ので現状は効かない | 低 | ソースの読み (gs では PHSROT が動かず確認できない) | v2 で修正 (N U (√λ ⊙ ξ))。k64 モードは PHSROT が従来どおり動かない |

## ## 0. 一覧 (重要度順)

取り除いた部分: 「B130 で v2.1.0 が修正」

元の行:

> | B81 | スピン追跡 sprot の縦磁場: `bz = bz0/pzm + bsi` で bsi = bz·al を p0 で規格化したまま z 軸の向きに足す。(1) δ≠0 で (1+a)B∥ の歳差が 1/(1+δ) 倍にならない (δ = −5e-3、0.12 rad で 6e-4 rad)、(2) 軌道が傾くと bsi·ẑ の運動量に垂直な成分が B⊥ に二重に入る、(3) 双極 + ソレノイド (tdrift) の `bsi0 − ak0x·y − ak0y·x` が T-BMT と 3.6e-3 rad 合わない | sprot (SOL 領域・UND の中のスピン追跡) | SOL 領域・UND の中のスピンの歳差 (POL + RAD の追跡) | 中 | ソースの読み + Thomas-BMT の数値積分 (v2 は 1e-7〜1e-5、K64 の式は 6e-4〜4e-3 ずれる)。gs での照合は記録に無い | v2 で修正 (SOL 領域と UND では運動量の向きの磁場の積分 ∫B∥dl/(Bρ) を運動量の向きに置く。k64 モードは再現)。Emittance[]・EMIT の POL の量 (同じ縦磁場の足し方) は B130 で v2.1.0 が修正 (k64 モードは K64 の式のまま) |

## ## 0. 一覧 (重要度順)

取り除いた部分: 「v2.1.0 の確認 (A2):」

元の行:

> | B83 | スピン追跡 sprot の歳差の角が (1+aγ)·2 sin(θ/2) (外積 sin θ を中点の pz ≈ cos(θ/2) で割る)。刻みが粗いと θ²/24 小さい (0.26 rad の BEND を 3 刻みで 3e-4) | sprot (すべての要素の追跡のスピン) | 平面リングの 1 周のスピンの回転角が 2π aγ から 4e-3 rad ずれる | 低 | ソースの読み + 解析式。v2.1.0 の確認 (A2): 24 セル FODO (3 GeV、RF 4 MV) の 1 周のスピン角が k64 −1.1940907 rad、v2 −1.1935299 rad (差 5.6e-4 rad/周)。偏向磁石に `EPS=1E-5` を付けると gs −1.1935395510、v2 −1.1935395636 (1e-8 で一致)。k64 の誤差は刻みを細かくすると消え、v2 は 57 倍速く収束する。example/07 の 500 周後の保存値 0.8794 rad (gs・k64) は v2 の 1.1516 rad と違うが、**この K64 の刻みの誤差が 500 周に積算した値で物理の値ではない** (EPS=1E-5 で gs 0.3969、v2 0.3970) | v2 で修正 (厳密な θ。k64 モードは再現) |

## ## 0. 一覧 (重要度順)

取り除いた部分: 「(従来は NotImplementedError。k64 モードは従来どおり)」

元の行:

> | B85 | `phis=merge(0.d0,merge(phic,w*trf0,autophi),trpt)` は TRPT でないとき phis = phic にするので dphis = 0 になり、足したばかりの DZ·w も PHI も効かない。光学 (tcave.f: phii = phic) とも合わない。tcave.f は autophi のとき dvcacc・ddvcacc を足さないので、空洞がすべて AUTOPHI だと実効の RF の傾きが 0 | tcav.f、tmulti.f、tcave.f | AUTOPHI のリングの追跡 (位相の差が消える)・Emittance[] の RF の傾き。SADHelp に AUTOPHI の説明は無い (意図は未確認) | 低 | ソースの読み。TRPT の場合は gs と全桁一致 (v2 の定義 = K64 の TRPT の動作) | v2 は TRPT のときの形 (phis = 0) にそろえ、追跡・CALC (TRPT の加速)・Emittance[] に AUTOPHI を実装 (従来は NotImplementedError。k64 モードは従来どおり) |

## ## 0. 一覧 (重要度順)

取り除いた部分: 「k64 モードは SOL のずれ・回転が従来どおり未対応」

元の行:

> | B86 | tsol.f の trots の入口は z に DZ だけを足し、出口は DZ·(1 − dv) を足す | tsol.f (trots) | 入口の DZ≠0 (GEO が出口側の SOL にあるとき) の z の差 dv·DZ (3 GeV で 1e-13 m 以下) | 低 | ソースの読み (参考。gs では未確認) | v2 は入口・出口とも (1 − dv)。k64 モードは SOL のずれ・回転が従来どおり未対応 |

## ## 0. 一覧 (重要度順)

取り除いた部分: 「(v2.1.0 の C: K64 は」

元の行:

> | B99 | tvcorr.f の corlength≠0 の枝は `s = 1/sqrt(4Σw²+2)` で、複素数 c の絶対値の 2 乗平均を 1 に規格化する。tcorr.f は DX = Re(c)·delx、DY = Im(c)·dely と使うので DX の rms は delx/√2。同じサブルーチンの corlength=0 の枝は成分ごとの規格化 (1/2 なし) で整合しない | tvcorr.f | corlength≠0 の DELCOR のずれの rms が delx の 1/√2 になる | 低 | 規格化の根拠は確定 (v2.1.0 の C: K64 は複素数の ⟨\|c\|²⟩ = 1 にする。corlength = 0 の枝は実部・虚部それぞれ分散 al/n で、2 つの枝で実部・虚部の分散の約束が揃っていない。試験: k64 は ⟨\|c\|²⟩ = 1.0、v2 は 2.0)。「DX の rms = delx」が作者の意図かは推定のまま (作者の文書による裏づけは無い) | v2 で修正 (s = 1/sqrt(2Σw²+1)。DX の rms = delx)。k64 モードは K64 のまま (delx に 1/√2 を掛けても同じ) |

## ## 0. 一覧 (重要度順)

取り除いた部分: 「、gate なし。v2.0.0 は断片を常に入れていた pysad の不具合で、K64 の動作に合わせた)」

元の行:

> | B113 | tfgeo の `inext`・`iprev` (「同じ型 (BEND) の次の要素。間は長さ 0 の要素だけ」を繋ぐ。リングの終わりは先頭に回って探す) により、先頭が長さ 0 の MARK の輸送路の最後の BEND は自分自身に繋がる (iprev = inext = 自分)。tbend・tbendi は `iprev == 0 .and. fb1 /= 0` で入口側、`inext == 0 .and. fb2 /= 0` で出口側に F1 の放射の断片を入れるので、この BEND は F1 の断片の放射が無い。隣り合う BEND の間も同様 | tfgeo.f、tbend.f、tbendi.f | F1・FRINGE のある BEND の RAD の追跡 (1 個の BEND の輸送路で 3%、隣り合う BEND の 12 セルのリングで 3.6%)。意図かどうかは未確認 | 低 | 確定 (gs で再現: 1 個の BEND (INS) の 1 周 δ = −2.17142054679e-5)。意図は疑い | pysad は K64 の動作を再現 (`track6n.bend_links`、両モード、gate なし。v2.0.0 は断片を常に入れていた pysad の不具合で、K64 の動作に合わせた) |

## ## 0. 一覧 (重要度順)

取り除いた部分: 「。v2.4.0: 極の中を」

元の行:

> | B115 | RADINT は ANGLE=0・K0≠0 の BEND (挿入光源の極: 基準軌道は曲がらず磁場だけ) を積分に入れない | INTGRL.f (RADINT) | ASTRID2 + MPW (0.58 GeV、周長 45.7 m) で U0 6.17 keV (Emittance[] 7.80 keV)、εx 1.268e-8 (Emittance[] 1.056e-8) | 中 | 確定 (gs で再現。I2 は Σ(ANGLE+K0)²/L と 1e-6) | v2 で修正 (K0 を曲げ角として扱う。端の角なし、TOTAL BENDING ANGLE には足さない)。U0・I2・I3・σδ は使える (U0 1.1e-3、σδ 2.2e-3)。v2.4.0: 極の中を直線の基準軌道の中の場として積分し (扇形の h² の収束・I4 の h³ の項を入れない)、CALC の分散が極の分散の発生 (η' += K0 δ) を含まない分を 1 周の周期解で補正した。I1〜I5 は独立な積分と 5e-5、εx・Jx は RF なしの Emittance[] と 1e-4・2e-4 で一致 (従来の 3.5%・0.7% の原因: 極の扇形近似 + CALC の分散。RF あり Emittance[] との +0.7% (Jx) は RF の位置の分散の物理的な効果で、独立な 4 次元の写像が再現)。k64 モードは再現 |

## ## 0. 一覧 (重要度順)

取り除いた部分: 「(従来の 3.5%・0.7% の原因: 極の扇形近似 + CALC の分散。」

元の行:

> | B115 | RADINT は ANGLE=0・K0≠0 の BEND (挿入光源の極: 基準軌道は曲がらず磁場だけ) を積分に入れない | INTGRL.f (RADINT) | ASTRID2 + MPW (0.58 GeV、周長 45.7 m) で U0 6.17 keV (Emittance[] 7.80 keV)、εx 1.268e-8 (Emittance[] 1.056e-8) | 中 | 確定 (gs で再現。I2 は Σ(ANGLE+K0)²/L と 1e-6) | v2 で修正 (K0 を曲げ角として扱う。端の角なし、TOTAL BENDING ANGLE には足さない)。U0・I2・I3・σδ は使える (U0 1.1e-3、σδ 2.2e-3)。極の中を直線の基準軌道の中の場として積分し (扇形の h² の収束・I4 の h³ の項を入れない)、CALC の分散が極の分散の発生 (η' += K0 δ) を含まない分を 1 周の周期解で補正した。I1〜I5 は独立な積分と 5e-5、εx・Jx は RF なしの Emittance[] と 1e-4・2e-4 で一致 (従来の 3.5%・0.7% の原因: 極の扇形近似 + CALC の分散。RF あり Emittance[] との +0.7% (Jx) は RF の位置の分散の物理的な効果で、独立な 4 次元の写像が再現)。k64 モードは再現 |

## ## 0. 一覧 (重要度順)

取り除いた部分: 「実軸上の z > 1 は従来どおり切断線の下側」

元の行:

> | B119 | 非整数次の `PolyLog[s,z]` が複素の z で、Re z < 0・Im z ≤ 0・\|z\| > 1 のとき誤る。整数次と第 3 象限以外は正しい | gamma.f (`cpolylogl`) | 例 `PolyLog[1.5,-3.28-3.54 I]` = 3.47+9.43 I (真値 −2.09−1.06 I) | 低 | 確定 (冪級数・spence・積分表示・整数次の反転公式と 2e-11 で比較) | v2 は真値 (冪級数・μ = ln z の展開・倍角の公式の再帰)。k64 モードも真値 (B118 と同じ扱い。要判断)。負の実数の z は上側 (arg z = +π)、実軸上の z > 1 は従来どおり切断線の下側 |

## ## 0. 一覧 (重要度順)

取り除いた部分: 「(v2.2.0 で原因を特定: 右端の欄の更新。v2.1.0 までは「疑い」で、gs も pysad k64 も 6.18e-9 だった)」

元の行:

> | B128 | `SynchroBetaEmittance` が大きなリング A (3 GeV、周長 603 m、νs = 0.0041、RF 1.4 MV) で 6.18e-9 を返し、Emittance[] の 1.55e-10 と 40 倍違う。AzimuthalModes 3〜8 で変わらず、12・20・30・32 以上で 3.19e-9・1.46e-9・8.50e-10・7.87e-10 と動き続ける。νs の入力・エネルギー・6 極の強さにも依らない | SynchroBeta (temits) | 大きなリングで使えない (独立な 6x6 の平衡 1.5173e-10 と、RAD + FLUC の 6 次元追跡 (60 粒子、4000 周) の 1.6e-10 と合わない)。PF・ASTRID2・SP8 (大きなリング B、8 GeV) は Emittance[] と 0.5〜1% で合う | 中 | 確定 (v2.2.0 で原因を特定: 右端の欄の更新。v2.1.0 までは「疑い」で、gs も pysad k64 も 6.18e-9 だった) | v2 は値を変えず、AzimuthalModes m と 1.5 倍の m の emitx が 20% 以上違えば警告を出す。k64 モードは gs と同じ (警告なし)。使うときは AzimuthalModes を変えて収束を確認し Emittance[] と比べる。`AzimuthalModes->1` は pysad が IndexError (gs は 8.6e-7 を返すが意味がない) **追記 (原因を特定して v2 で修正)**: K64 の `tefsetup` が 1 周の写像 X' = T X + h1 の源泉 h1 に閉軌道 c(δ) (1 次の分散 η δ を含む) をそのまま使う (固定点が c にならず、2 次元の回転では平均が c の半分になる)、`tesumb` が x_b の 2 次のモーメントに η η σδ² だけを足し相関 η σδ² (`beam(16:19)`) を入れない。リングの始点の分散 η ≠ 0 で εx が始点に依り数倍〜40 倍過大 (24 セル FODO の始点を QF・BEND・BEND の途中にして 6.54・6.54・4.74 倍、gs でも同じ。η が小さい PF・ASTRID2 は 1%)。確度は確定。v2 は源泉から 1 次の分散を引き、軌道の δ 依存の効果を解析的な 2 項に分けて入れる (`tesumb_v2`): Emittance[] の 1.0003〜1.0005 (FODO)・1.0084 (大きなリング A、AzimuthalModes 3。32 で 1.0011)・1.0045 (SP8 (大きなリング B))・1.005 (PF)・1.006 (ASTRID2)。k64 は K64 のまま (gs と 2e-7)。警告 (上記) は v2 の大きなリング A では出なくなった |

## ## 0. 一覧 (重要度順)

取り除いた部分: 「k64 モードは従来どおり空間電荷を入れない」

元の行:

> | B133 | WSPAC ON の CALC・CALC6D は入口以外の Twiss が NaN になる (0.2 GeV のリングで PBUNCH = 1e4 でも `Horizontal/Vertical Unstable` と NaN) | qtwiss.f (qwspac)・qtwissfrac | WSPAC の光学の CALC が使えない (Emittance[]・EMIT は使える: B127) | 中 | 確定 (gs で再現。原因は未調査) | v2 は CALC6D と RFSW ON の CALC (4D) に空間電荷を入れる (6D の経路。Emittance[] のチューンと 2e-6、一次のラスレットと 0.15%)。K64 の `qwspac` の直訳ではないので gs との比較はできない。k64 モードは従来どおり空間電荷を入れない (gs が NaN) |

## ## 0. 一覧 (重要度順)

取り除いた部分: 「~~残る 2.7% (AE に比例、極の長さ・間隔に依る) は原因未特定~~ → **[v2.4.0 D2 の訂正]** K64 の誤りではなく、基準側の四極のフリンジ未導入と U0 の定義 (閉軌道のエネルギーの沈み込み) の差だった (B222)」

元の行:

> | B143 | `tmulte`・`tmulti` (MULT) の端の角 AE1・AE2 の蹴りは、放射の場の中の行路を変えない (B75・B108 と同じ型) | tmulte.f・tmulti.f | 端の角つきの MULT (勾配つきの極を含む) の放射: 端の項が約 2 倍 (厳密な幾何との ΔJx の比 K64 : v2 = 2.0 : 1.03) | 中 | 確定 (硬い端の磁石の厳密な幾何 + ロビンソンの定理。表は `history/v2next_P3.md` §3.3) | v2 で修正 (追跡 `tsteer`・`tmulti` の `zr0` に x tan(e)、Emittance[] の `_tmulte` の `tradke` に tin・tout)。厳密な幾何と 3〜4%。~~残る 2.7% (AE に比例、極の長さ・間隔に依る) は原因未特定~~ → **[v2.4.0 D2 の訂正]** K64 の誤りではなく、基準側の四極のフリンジ未導入と U0 の定義 (閉軌道のエネルギーの沈み込み) の差だった (B222)。ANGLE=0 の BEND の極は D2 で厳密解と減衰率 1.2e-5 以内 (B220・B221)。k64 モードは K64 のまま。B108 の MULT への拡張 |

## ## 0. 一覧 (重要度順)

取り除いた部分: 「**v2 のオプション (v2.3.0 U2):」

元の行:

> | B160 | **誤りではなく K1 の暗黙の定義 (U1 で再分類)**: K1≠0 の BEND の本体は、中央面の磁場 B_y = B0 + (K1/h) ln(r/ρ0)、B_r = K1 y ρ0/r (= (1 + hx)∂b/∂x が一様、B''(ρ0) = −h B'(ρ0)、x y² の項なし) の中の運動を厳密に実装している。elegant (CSBEND K2 = 0) など B' が x に一様な定義とは K1 が同じ記号でも別の磁石で、本体に K2 = −h K1 の6 極成分を足したのと同じ。SADHelp の BEND の K1 は B'L/Bρ とだけ書き B'' を指定しない | BEND の本体の写像 (`tbend.f` tbendiinit・tbendicorr・tbendibody、`tbendi.f`、CALC の `tbende.f` tbendecorr) | オフモメンタムのチューン (ξ・∂²ν/∂δ²)、LMA・DA・クロマティック ADTS を B' 一様の定義と比べるとき。δ = 0 の光学・エミッタンスには効かない。n = 0.4 の扇形磁石 16 台 (3 GeV) の dν/dδ は gs・pysad・独立な積分 'log' が (−0.17213、+0.21082)、B' 一様 (独立な積分・解析式・elegant) は (−0.60246、+0.73786)。ASTRID2 の自然な ξy は SAD −11.8、elegant −8.2 | 低 (仕様・癖。SAD → 他コードの変換は K2 = −h K1 を足す。docs/11 A59) | 実装の同定は確定 (ソース、導出、独立なローレンツ力の積分と 16 台の非線形の写像が 1e-10、MULT (ANGLE≠0) の別経路が 1e-12 で同じ)。「定義であって誤りではない」は中〜高 (BEND の K1 の B'' は文書にない)。T1 の「疑い〜確定の間」を更新 | 直さない (k64・v2 とも既定は K64 と同じ 'log')。**v2 のオプション (v2.3.0 U2): `compat.set_bend_k1_def('lin')` (環境変数 `PYSAD_BEND_K1_DEF=lin`) で B' 一様の定義 ('log' + 6 極 K2 = +hK1) にできる (k64 では無視)。n = 0.4 の 16 台の dν/dδ が (−0.60246、+0.73786) になり、独立な積分・解析式・elegant (K2 = 0) と一致。ASTRID2 の自然な ξ も elegant (K2 = 0、EDGE_ORDER = 2) と 4 桁で一致。`history/v23_U2.md`** CALC の tbendecorr の dk の不整合は B161 |

## ## 0. 一覧 (重要度順)

取り除いた部分: 「**v2 で直した (v2.3.0 U2、`compat」

元の行:

> | B161 | CALC・Emittance[] の BEND (K1≠0) の解析行列 `tbendecorr` の dk = −K1 xr (1 − xr (1 − 2 xr)) (xr = x/ρ0) が、同じ関数の kick の導関数 −K1 xr (1 − xr/2 + xr²/3) と xr の 2 次の項が違う (xr = 0.001 で −0.05%、0.1 で −3.5%、0.3 でたまたま一致)。写像 (kick) と行列の不整合 | `tbende.f` tbendecorr の 86 行 (kick)・87 行 (dk) | CALC・Emittance[] の転送行列だけ (追跡は kick を使うので影響を受けない)。スライス 1 cm・K1 = −0.06 /m²・xr = 1e-3 で t[1,0] に 3e-10、リングのチューンで 1e-8 程度以下と見積もった (全体の量は見ていない) | 低 | 確定 (U2: CALC の行列が kick の数値微分と合わないことを `bendecorr` で確認、K1 の強い扇形磁石のリングで CALC のチューンが独立な積分と合わない: δ = 1% で 3e-5、2% で 1e-4、4% で 5e-4。K64 の作者の意図は不明。実際のリングへの影響は小さい: ASTRID2・大きなリング A で ξ に 4e-6、∂²ν/∂δ² に 0.01、ν(δ = 2%) に 5e-7、Emittance[] の EMITX・SIGE に 1e-9 (相対) 以内) | **v2 で直した (v2.3.0 U2、`compat.bend_dk_fix()`、`perfkern_bend.py` `bendecorr` の dk を kick の導関数 −K1 xr (1 − xr/2 + xr²/3) に)。k64 は K64 のまま。** `history/v23_U1.md` §8、`history/v23_U2.md` |

## ## 0. 一覧 (重要度順)

取り除いた部分: 「6D 固定点を Newton 法で解き (従来の動作)、」

元の行:

> | B163 | **K64 の `tcod` (RFSW ON・RADCOD ON の 6D 閉軌道) は、最初の Newton 法で運動量の誤差が大きい (`abs((dz−dtrf0·v0)·dvcacc/pgev) + abs(ddp) > dpthre`) と横方向の修正 `dcod(1:4)` を 0 にし、r が下がらないとその `dcod` を使い回して `fact` を半分にし続け、`fact < fmin` で `fndcod = .false.` のまま終わる。`temit` は `???-Emittance[]-closed orbit not found.` を出して、収束していない cod から Emittance[] の計算を続ける** (Stable = 0)。強い 6 極のあるリングで起きる (DDBA_15、ノーマルセルのリング: 6 極の強さを 0.5 倍にすると収束する)。結果は閉軌道が固定点でない状態で出るので、εx・σE・αc・U0 が数 1e-3 ずれる | tcod.f (`dcod(1:4)=0` と `r >= r0` の分岐の組み合わせ)、temit.f:148-152 | DDBA_15 (5 GeV、731 m、6 極 600 台) の εx が gs 9.6006e-10、収束した閉軌道の値 9.632e-10 (+0.33%)。設計者のメモの値 9.6310e-10 は収束した値と一致 | 中 | 確定 (gs の出力に "closed orbit not found"、pysad の `tcod` のトレースで同じ停止を再現、k64 の再現で gs と全桁一致) | k64 は gs と同じに収束しない軌道から続行し、メッセージを出す (`Stable` = 0)。v2 は放射込みの一周写像の 6D 固定点を Newton 法で解き (従来の動作)、メッセージは出ない。試験 `test_k64_tcod_notfound.py` |

## ## 0. 一覧 (重要度順)

取り除いた部分: 「は v2.4.0 D3 で対応 (B230・B231・B234)」

元の行:

> | B165 | **磁石の本体の運動量の向きの磁場 (1+a)B∥ ds が入らない (B136・B137 の原因の特定)**: K64 の追跡 (tradk の sprot)・Emittance[] (tradke) は、運動量の向きの変化 (B⊥) と端の縦磁場 ψ (BEND・QUAD の入口 +ψ、出口 −ψ。ψ = y/ρ + k x y) だけでスピンを歳差させ、本体の中の B∥ = B·p̂ の積分を入れない。静磁場の ∫B∥ ds = ∫B·dr = ψ_出口 − ψ_入口 (スカラーポテンシャルの差) なので、本体 + 両端の B∥ の合計は 1 次で 0 になるはずが、K64 は端の ψ だけが残る (v2 で B58 により出口の ψ も入れたので、k64 の「出口の ψ の欠落」が偶然 B∥ の欠落を隠していた入口 y = 0 の場合を除き、k64・v2 とも合わない)。B136 の「厳密な値 (一様磁場が出口の面で急に切れる非 Maxwell 場) との差 (2+a)θ py」のうち θ py の部分は、この端の ψ (非 Maxwell の基準には無い) | tradke・tradk・tbend の本体 | py = 5e-3 の 15° 扇形 BEND (0.5 GeV) で R[0,1] が追跡 −4.08e-3・Emittance[] −4.07e-3、Maxwell 場つきの硬い端の Thomas-BMT の積分 −2.77e-3。dR/dpy の 3 要素 (R01・R10・R12・R21) の誤差 0.26 (R01・R10 の 47%) | 中 | 確定 (独立な基準: スカラーポテンシャル Φ = −b y G から作った Maxwell 場の中の Thomas-BMT の数値積分 `tests/_spin_maxwell_indep.py`。フリンジ幅 w = 4e-3〜5e-4 m で 1e-6 以内に収束。軌道 (y・py) は追跡と 1e-9 で一致) | **v2 で修正 (BEND の K1=0・K1≠0 の本体)**: 各放射・スピンの刻みに (ψ_終 − ψ_始)/(1+δ) を運動量の向きの磁場として足す (追跡 `p.bpadd`・Emittance[] `st.bpot`)。単一 BEND で dR/dpy の誤差: 追跡 R01 −8e-4 (K64 −0.26)、R12・R21 1e-2、Emittance[] (EMITDIV = 256) は全要素 4e-4。残る追跡の差は、端の ψ を刻みの中点の運動量の向きで分解する K64 の近似。k64 は K64 のまま。**QUAD・MULT・SOL 内の四極の本体の B∥ は v2.4.0 D3 で対応 (B230・B231・B234)。SEXT・ステアリングは端の ψ も本体も入れないので和が 0 で元から正しい。追跡の R12・R21 の 1e-2 の差は、端の ψ の分解でなく K64 の刻みの分け方の 1 次誤差で、刻みの偏向角を 5e-3 rad 以下にして 6e-4 (B235)**。切り替え: `compat.set_spin_bpar(False)`・`PYSAD_SPIN_BPAR=0`。試験 `test_v24_D_spin_bpar.py` |

## ## 0. 一覧 (重要度順)

取り除いた部分: 「(v2.4.0 F2。`->False` で」

元の行:

> | B200 | **Touschek の損失率の基準 δ_old = 0 (K64 の `TouschekTable`・`touckf`、Piwinski 型、elegant `touschekLifetime` も同じ) は、バンチ内の δ の分布を無視している**。散乱後の 2 粒子の δ は δ_c ± Δ (δ_c = 対の重心の δ ~ N(0, σδ²/2)、Δ は相対運動量の散乱だけで決まる) で、損失の条件 (δ' > 受容) は a ∓ δ_c を基準に判定すべきだが、K64 は元の δ = 0 の粒子を散乱させる kick の判定 (δ_new − 元の自分の δ) のため、a が σδ の 10 倍程度でも損失率が 1〜5% 小さく (寿命が長く) 出る。相対運動量の縦成分 q_s の効果も同じ基準で落ちる。a/σδ = 9 で ASTRID2 は損失率 +4.2〜4.5% (δ_c + q_s)、14 で +1.9%、23 で +0.7%、KEK-LS は 8 で +4.5%、13 で +1.8% | `tintrb.f` touckf・`TouschekTable` (δ_old = 0 の kick)、emit.n `TouschekF` | Touschek 寿命の絶対値 (受容が σδ の 5〜15 倍のとき 1〜5%)。受容が σδ の 3 倍以下では打ち切りの問題で意味がない | 低〜中 | 確定 (実験室系の 2 粒子の対のモンテカルロ (厳密な相対論の運動学) で δ_c の効果が理論と 2% 以内、q_s の効果を含む全体が 3 次元の積分と 0.3% 以内。`v24_A.md` A3・A10、`example/v2_4/A3_mc_bunch.py`) | k64 は K64 のまま (δ_c = 0)。**v2 は `TouschekLifetimeLMA[]`・`TouschekLifetime[INF, INF, nz]` (Piwinski 型) の既定を、バンチ内の δ の分布 (BunchDelta) と相対運動量の縦成分 (BunchQs) つきにした** (v2.4.0 F2。`->False` で δ_c = 0 の基準に戻せる)。elegant との比較 (T1・T2) は「δ_c = 0 の基準との差」になる (`docs/48` §4.4)。試験 `test_v24_F2_defaults.py` |

## ## 0. 一覧 (重要度順)

取り除いた部分: 「(v2.4.0 F2。明示すれば」

元の行:

> | B201 | **`DynamicApertureSurvey` の `ExpandElementValues` の既定 True は、補正磁石の成分ごとの値 (K0・位置ずれ DX・DY) をファミリーの代表の値に揃えて消す**。補正した状態のリング (軌道補正で成分ごとに K0 を変えた状態) の DA を、補正のない状態の値で測ってしまう (help の注意 'Machine errors may be reset' のとおり。既定のまま使うと誤差・補正の研究で DA が実際と違う): DDBA v07_14 の例で補正後の DA は 29〜36σ だが、既定では補正が消えて 0〜38σ の偶然の値になる (gs・pysad とも)。小さな FODO (QF.3 の DX = 5 mm) でも既定 True は DA が 15 → 51 (ずれが消える) | Packages/dapert.n の既定 `ExpandElementValues->True` | 成分ごとに値を変えた状態 (誤差・補正) のリングの DynamicApertureSurvey | 中 | 確定 (gs で再現: `ExpandElementValues->False` で DA は 34σ・29σ、既定 True で 0〜1σ。`v24_A.md` A5) | k64 は既定 True のまま (gs と同じ)。**v2 は既定を False にした** (v2.4.0 F2。明示すれば両モードともその値。`->True` で従来どおり揃える)。試験 `test_v24_F2_defaults.py`・`test_v24_A5_dapert_options.py` |

## ## 0. 一覧 (重要度順)

取り除いた部分: 「v2.4.1 の B292 で θ 依存と ψ に直し」

元の行:

> | B224 | `tradke` の D・B の更新 (photons.f の `trans`・`transr` から区間の写像 tr2 = T·transr⁻¹ を作り、`tinv6(r) = −J rᵀ J` (シンプレクティックの逆行列の式) を使う) は、区間の始点・終点の行列がシンプレクティックであることを前提にしている。K64 の写像は端の面が薄い蹴りで近似的にシンプレクティックなので問題が出ないが、厳密な硬い端 (B223) の MULT の端の写像は K1 ≠ 0 で運動学的運動量の写像として厳密である代わりにシンプレクティックでない (K1 = 0.05・AE = 0.3 で 5e-4、K1 = −0.01・K0 = 0.005 で 3e-2) | photons.f の tradke (pysad: `perfkern.tradke_db`) | 逆行列の式が効かず、δ の行 (dδ/dy) が数値微分の 2〜3 倍ずれ、K1 < 0 の MULT の AE のある Emittance[] の ΔJx が −0.011 になっていた。MULT の厳密な端には縦方向の収束 (Maxwell) も無かった (K64 の薄い蹴りの縦成分 Δpy = fe·Im(cx) が抜けていた。ANGLE=0 の BEND の端の鎖と y・py が 8e-4 違った) | 中 (K1 ≠ 0 の MULT の極) | 確定 (数値微分との一致、ANGLE=0 の BEND の端の鎖との比較) | v2 で修正: `tradke_db` は v2 で厳密な逆行列 (`np.linalg.inv`) を使い、MULT の端に薄い蹴りの縦成分を足した (BEND の鎖と 2e-5 以内。v2.4.1 の B292 で θ 依存と ψ に直し、独立な Maxwell 場の基準と 2e-8 以内)。δ の行が数値微分と 1e-17 で一致し、K1 < 0 の MULT の ΔJx は 1e-5 以下。k64 は K64 のまま |

## ## 0. 一覧 (重要度順)

取り除いた部分: 「**v2 で直した (v2.4.1)**」

元の行:

> | B226 | **K64 の四極の本体の積分器 (tsolqu・tsolque) の刻みの限界**: [非線形ドリフトの半分 (tqente)] + [線形の本体 (厳密)] + [半分] の Strang 分割の 2 次の誤り。刻み数 1 + int(\|L·K1L\|/0.1) (CALC)・/0.2 (追跡) は軌道ずれの大きさに依らず、誤りは (aln²·k)·p² (p = 軌道の傾き) に比例する | tsolque.f・tsolqu.f | 軌道が x = 3 mm・px = 2 mrad ずれた四極 (L = 0.3 m、K1 = −0.6) の 1 台のヤコビアンが厳密な写像と、CALC 1.8e-7・追跡 7e-7 違う (写像は 1e-9)。追跡と CALC のチューンが K0² で増える差 (K0 = 0.006 で 2.3e-6、0.012 で 9e-6) の原因 | 低 (軌道ずれ 100 µm 以下で 1e-9 以下) | 確定 (四極 1 台の厳密な写像 (RK4) との比較) | **v2 で直した (v2.4.1)**: 原因は Strang 分割の 2 次の誤り (H = H_lin + H_nl の D(h/2) B(h) D(h/2)。D = 運動学的な非線形ドリフト、B = 厳密な線形の本体) と、刻みの数が軌道のずれに依らないこと。v2 は **4 次の Yoshida の合成** S(w1 h) S(w0 h) S(w1 h) (w1 = 1/(2−2^(1/3))、w0 = 1−2 w1 < 0。D・B は厳密な流れなので 4 次のシンプレクティックな写像) にした。刻みの基準 ε は追跡が 0.2 → 0.1、CALC が 0.1 → 0.05 (B243)。11 か所に重複していた刻みの構造は `pysad/quadbody.py` に集約 (`group()`・`ndiv()`・`make_track()`。追跡 python 版・numba 版・MULT・放射つき・CALC python 版・numba 版 (`_quad_full`・`k_quad`)・解析ヤコビアン・MULT の `tsolque_pieces`・emit6d の放射の刻み)。**k64 は従来とビット単位で同じ** (追跡 python・numba、CALC・Emittance[] の numba・python、MULT、放射あり・なしの 146 量を従来の版と repr 比較して不一致 0)。独立な基準 (運動方程式と変分方程式の DOP853 積分): 3 mm・2 mrad ずれた四極 (L = 0.3、K1L = −0.6) のヤコビアンの誤り CALC 1.8e-7 → 5e-10、追跡 7.1e-7 → 7.9e-9。1 要素の写像の誤り (10 mm・5 mrad): 3.4e-8 → 1.3e-10 (L = 0.3)、4.9e-9 → 4e-11 以下 (L = 0.5、K1L = 1.5)、2.2e-9 → 1e-11 以下 (L = 1.0、K1L = 1.5。刻み 1 つあたりの評価は 3 倍、刻みの数は 1〜2 倍)。刻み 2 倍で誤りは k64 が 1/4 (2 次)、v2 が 1/16 (4 次)。**残り**: 追跡と CALC のチューンの差 (ステアリング K0 = ±0.006 の FODO。k64 は +2.3e-6、0.012 で +9.2e-6) は v2 で −6.3e-8 (0.012 で −2.5e-7) に減る (36 分の 1 だが 0 ではなく符号が変わる。CALC と追跡の刻みの基準の残りの違いと、1 刻みあたりの誤りが原因と推定。確認していない)。k·aln² > 1 の強い 1 刻み (L = 1 m・K1L = 1.5 で EPS を大きくして 1 刻みにした場合など) は Yoshida の方が Strang より悪い (6.8e-7 対 4.5e-7) ので、EPS を大きくしてはならない。v2 の追跡は四極 1 台の評価が numba で 2〜2.5 倍 (1 台 0.05〜0.1 µs → 0.1〜0.22 µs。刻みあたり本体 3 回・D 4 回で刻みの数は 1〜2 倍)。リング全体の追跡時間は測っていない。試験 `test_v241_quad_integrator.py`。SOL 内の四極 (bz ≠ 0 の `tsolqu`・`solemit`) と MULT ANGLE ≠ 0・偏向磁石の本体 (tbendebody) は別の積分器で、この変更の対象外 |

## ## 0. 一覧 (重要度順)

取り除いた部分: 「以前は K1 < 0 の MULT の AE 依存の差と見えていた」

元の行:

> | B227 | 追跡の放射の刻み (`ndivrad`、要素の EPS の既定) が粗く、放射つきの追跡の一周写像から求めた Jx (減衰の分配) が刻みの誤りで 0.27% 小さい | tradk を呼ぶ tsteer・tbend・tmulti ほか (K64 の追跡の刻み) | 24 セルの FODO + 3 極 (3 GeV): 追跡の Jx は既定の刻みで 0.9272、EPS = 0.1 で 0.9292、0.01 で 0.9294。Emittance[] は EMITDIV = 1〜80 で 0.9296 に収束 (BEND の極でも MULT (K1 = 0・±0.01) でも同じ差)。以前は K1 < 0 の MULT の AE 依存の差と見えていた | 低 | 確定 (EPS の走査) | 直していない (v2 でも K64 の刻み。細かくすると放射つき追跡が遅くなる。EPS を台本で与えれば近づく) |

## ## 0. 一覧 (重要度順)

取り除いた部分: 「(F1・FB・AE が無いとき。以前は常に Python の経路)」

元の行:

> | B290 | **ANGLE=0・K1≠0・K0≠0 の BEND (B02) の原因 = tbend の tmulti 呼び出しと tbende の tmulte 呼び出しの引数の誤り**。(1) 追跡 `tbend` (tbend.f:584-593) は `tmulti` に `radius = 0.d0`・`rtaper = 0.d0` を渡す (tturn.f の MULT は rtaper = 1、これは K・端のフリンジに掛ける係数。`akr(0) = (akr0(0)+…)·rtaper`・`akr(1:nmmax) = akr0(1:nmmax)·rtaper`) ので K0・K1 の蹴りがすべて 0 になる。さらに `nmmax = merge(1, 0, ak == 0.d0)` は判定が逆 (ak ≠ 0 で 0)。(2) CALC・Emittance[] の `tbende` (tbende.f:403) は `tmulte` に引数を 2 個 (alg・phig) 渡し忘れる (tmulte は `…theta, dtheta, alg, phig, eps0, enarad, fringe, f1in…` の 27 個を取るが、呼び出しは 25 個。以降の実引数が 2 つずつずれて logical に実数が入る) | tbend.f:584・tbende.f:403 | **gs の追跡はこの BEND を長さ L のただのドリフトにする** (K0・K1 とも無視。x = 1.02e-3・px = 1e-4 の粒子が L = 0.2 m で変化なし、同じ K0・K1・L の MULT は px が −8e-5 など変わる)。**CALC (Twiss・Emittance[]) は SIGSEGV で落ちる** (USE・CELL は通る。B02 の原因。F1 つきも同じ) | 中 (補正用の複合機能磁石に使うと CALC が落ちる) | 確定 (実行: ANGLE=0・K1 = 0.05・K0 = 0.0008・L = 0.2 の BEND の 1 要素の追跡、同じ MULT との比較。ソースの引数の数) | **直す (v2・k64 とも)**: gs の誤った値 (ドリフト・SIGSEGV) は再現せず、tbend が tmulti に渡すはずだった値 (akm(0) = K0、akm(1) = K1·L、rtaper = 1) で追跡・CALC・Emittance[] を計算する = **同じ K0・K1・L の MULT と同じ** (gs の MULT と追跡・CALC・Emittance[] が全桁一致。F1 は MULT の FB1・FB2 = F1 + FB)。v2 は端の角 AE も追跡に使う (B142)。numba の追跡 (K_MULT) に対応 (F1・FB・AE が無いとき。以前は常に Python の経路)。独立な基準: 一様な K0・K1 の運動方程式と 4e-9 (k64)・6e-10 (v2)。docs/46 に gs の動作を書く。試験 `test_v241_A_bend_k1_k0.py` |

## ## 0. 一覧 (重要度順)

取り除いた部分: 「k64 は従来どおり入れない。試験は既存」

元の行:

> | B320 | **WSPAC ON の CALC・CALC6D の NaN と 4 次元 CALC の `Horizontal/Vertical Unstable` の原因 (B133 の追究)**: (1) K64 の `twspace` (tintrb.f) は `sy = sqrt(sx − a/2)` (σ² の主軸の小さい方) を、垂直のエミッタンスが 0 の (結合のない) リングで丸めの負の値 (εy ≈ −9e-31) から計算して NaN にする。kick の f = 2 N r0 L/(p1² h1 (sx+sy) √(2π) σz) が NaN。結合を入れて εy ≠ 0 にすると gs の CALC6D は有限の値を返す (0.2 GeV・PBUNCH = 1e8・四極 1 本を 45° 回転して εy = 2.2e-11: NX 2.8119763、NY 1.9003734、NZ −0.0176484)。(2) 結合があっても **4 次元の CALC (`qwspac`) は PBUNCH = 1e4 の弱い空間電荷でも `Horizontal/Vertical Unstable` と報告して Twiss を更新しない** (gs で確認: PBUNCH = 1e4・1e6・1e8 とも)。`qwspac` の直訳は、基準が存在しないので作らない (原因は qtwiss.f の呼び出し `rlist(ifsize+(l1−1)*21)` の大きさの表か、qwspac の行列の掛け方 (`trans` に右から 1 回 K を掛け、coup = F では (3,2)・(4,2)・(1,4)・(2,4) を 0 にする) と推定、確かめていない) | tintrb.f (twspace・qwspac)・qtwiss.f | WSPAC の光学の CALC (4 次元) が使えない。CALC6D は結合のあるリングだけ | 中 | 確定 (gs で再現。(1) は原因まで特定、(2) は原因未特定) | v2 は CALC6D と RFSW ON の CALC (4D) に空間電荷を 6D の経路で入れる (`spch.WSpace`・`calc6d`)。εy = 0 では `SpaceChargeError` で理由を出して止める (NaN を返さない)。k64 は従来どおり入れない。試験は既存の WSPAC の試験 |

## ## 0. 一覧 (重要度順)

取り除いた部分: 「(v2.4.0 で確認)。gs の追跡は」

元の行:

> | B321 | **領域の外の SOL (BOUND=0、BZ≠0) を、追跡 (tsol) は領域の始点として扱い、幾何 (tsgeo)・光学は扱わない** (K64 の内部の不整合): `tturn.f` の icSOL は `tsol` を呼び、`tsol` は始点の SOL の BOUND を見ずに次の BOUND≠0 の SOL (無ければ `???-TRACK-?Missing end of solenoid` を出して格子の終わり) までを領域として追う。一方 `tsgeo.f` は BOUND≠0 の SOL だけを始点とみなして ` Missing BOUND of Solenoid` を出し、領域の幾何を作らない。CALC は SIGSEGV (v2.4.0 で確認)。gs の追跡は、始点を BOUND=0 にしても BOUND=1 にしても 18 個の座標が全桁一致 (1.7 GeV、BZ = 0.8 T、QUAD ×4) | tsol.f・tsgeo.f・tturn.f | BOUND=0 の SOL が領域の外にあるラティスで、追跡と幾何・光学が別の領域を見る | 低 | 確定 (gs で再現) | **追跡は gs と同じ (BOUND=0 でも次の BOUND≠0 の SOL までを領域として追う。両モード)**。Twiss・CALC・幾何は従来どおり NotImplementedError (基準が無い)。試験 `test_v2_F4_sol_bound0.py` (gs と 10 桁、k64・v2) |

## ## 0. 一覧 (重要度順)

取り除いた部分: 「Twiss・CALC・幾何は従来どおり NotImplementedError」

元の行:

> | B321 | **領域の外の SOL (BOUND=0、BZ≠0) を、追跡 (tsol) は領域の始点として扱い、幾何 (tsgeo)・光学は扱わない** (K64 の内部の不整合): `tturn.f` の icSOL は `tsol` を呼び、`tsol` は始点の SOL の BOUND を見ずに次の BOUND≠0 の SOL (無ければ `???-TRACK-?Missing end of solenoid` を出して格子の終わり) までを領域として追う。一方 `tsgeo.f` は BOUND≠0 の SOL だけを始点とみなして ` Missing BOUND of Solenoid` を出し、領域の幾何を作らない。CALC は SIGSEGV。gs の追跡は、始点を BOUND=0 にしても BOUND=1 にしても 18 個の座標が全桁一致 (1.7 GeV、BZ = 0.8 T、QUAD ×4) | tsol.f・tsgeo.f・tturn.f | BOUND=0 の SOL が領域の外にあるラティスで、追跡と幾何・光学が別の領域を見る | 低 | 確定 (gs で再現) | **追跡は gs と同じ (BOUND=0 でも次の BOUND≠0 の SOL までを領域として追う。両モード)**。Twiss・CALC・幾何は従来どおり NotImplementedError (基準が無い)。試験 `test_v2_F4_sol_bound0.py` (gs と 10 桁、k64・v2) |

## ## 0. 一覧 (重要度順)

取り除いた部分: 「**v2.4.1 で修正** (`_tthine」

元の行:

> | B270 | **pysad 自身の誤り (K64 のバグではない): 台本の `Emittance[]` (`emit_extra._turn_record`) が、L ≠ 0 の 6 極・8 極に v2 の B240 (追跡と同じ 1 刻みの tradke) を渡さない** (`_tthine(st, el)` が enarad=False)。`emit6d.emittance6d` を直接呼ぶ経路は B240 が入るので、同じリングで台本の `SpinTune`・`PolarizationVector`・`Polarization`・`PolarizationTime` と `emittance6d` の値が違った | pysad/emit_extra.py | v2 で POL ON の台本の Emittance[] (6 極つきリング、軌道が 6 極の中心からずれるとき): ASTRID2 (縦ずれ rms 2 mm の四極、閉軌道 y 12 mm) の SpinTune が 0.31128 (正しくは 0.31491、差 3.6e-3)、PF (rms 0.4 mm) で 1e-3、小リング (6 極 K2L = 3、y 2.5 mm) で 3.2e-4。k64 は K64 と同じ (6 極の放射なし) で変わらない | 中 | 確定 (台本の `SpinTune` が `emittance6d` の nup と 1e-9 以内。追跡 (v2) のスピンチューンとの差は FODO 1.2e-5・ASTRID2 1.6e-5・PF 4e-5 に縮んだ) | **v2.4.1 で修正** (`_tthine(st, el, enarad=True)`。試験 `test_v241_pol_fluc.py::test_b270_*`) |

## ### B02 ANGLE=0 かつ K1≠0 かつ K0≠0 の BEND で FFS 初期化が SIGSEGV

取り除いた部分: 「**特定 (v2.4.1、B290)**」

元の行:

> - 原因: **特定 (v2.4.1、B290)**: 落ちるのは初期化でなく CALC (USE・CELL・追跡は通る)。`tbende` の `tmulte` 呼び出しの引数が 2 個足りない。追跡 `tbend` は `tmulti` の rtaper に 0 を渡すので、この BEND はただのドリフトになる。

## ### B02 ANGLE=0 かつ K1≠0 かつ K0≠0 の BEND で FFS 初期化が SIGSEGV

取り除いた部分: 「- pysad での扱い: v2.4.1 から追跡・CALC」

元の行:

> - pysad での扱い: v2.4.1 から追跡・CALC・Emittance[] とも同じ K0・K1・L の MULT と同じ値で計算する (numba の追跡も)。gs は CALC が落ちる・追跡がドリフトなので、gs との照合は MULT で行う (全桁一致)。B290 を見よ。

## ### B63-B71 後から見つかった癖・疑い

取り除いた部分: 「v2.0.0 のときの 3 件のうち、BEND の F1=0.1 FRINGE=3 の RAD ON 追跡の差 (B113) と L≠0 の CAVI の DX・DY の差 (CAVI の TENT・TEXIT が pysad に無かった) は v2.1.0 で解決した (gs と全桁一致」

元の行:

> - B69: どちらが正しいか決められない。pysad 側の誤りの可能性も残る。pysad と gs が合わない未解決の差のうち、B68・B69 以外は、SAD の誤りか pysad の誤りかも未判定なのでこの表に載せていない (docs/history/work_unimpl.md)。v2.0.0 のときの 3 件のうち、BEND の F1=0.1 FRINGE=3 の RAD ON 追跡の差 (B113) と L≠0 の CAVI の DX・DY の差 (CAVI の TENT・TEXIT が pysad に無かった) は v2.1.0 で解決した (gs と全桁一致、`history/v2next_A1.md`)。K3 が強い MULT リングの CALC (qcod) の閉軌道は未解決のまま。

## ### B72-B162 v2.0.0・v2.1.0・v2.2.0・v2.3.0 の作業で新しく見つけたもの

取り除いた部分: 「### B72-B162 v2.0.0・v2.1.0・v2.2.0・v2.3.0 の作業で新しく見つけたもの」

元の行:

> ### B72-B162 v2.0.0・v2.1.0・v2.2.0・v2.3.0 の作業で新しく見つけたもの

## ### B72-B162 作業で新しく見つけたもの

取り除いた部分: 「(B111〜B134 は v2.1.0。記録の項目名との対応は」

元の行:

> B72〜B162 は表 (§0) に症状・場所・影響・重要度・確度・pysad の扱いを書いた。詳細 (根拠の数値、試験、使った条件) は次の記録にある (B111〜B134 は v2.1.0。記録の項目名との対応は `local/wake/v21_bnum.md`)。

## ### B72-B162 作業で新しく見つけたもの

取り除いた部分: 「- B111〜B134 (v2.1.0):」

元の行:

> - B111〜B134 (v2.1.0): `history/v2next_*.md` の各担当の「K64 の誤り・新しく見つけたもの」。B111・B112・B113 は A1 (項目 7・2・4)、B114・B115・B116 は A2、B117 は B1 (項目 2)、B118〜B122 は B2、B123・B124 は B3、B125 は B3 と D1、B126 は C (B107 の 2 つ目)、B127・B128・B132 は D1、B129 は D2、B130・B131 は E、B133・B134 は B1。B40・B60・B81・B83・B99・B106〜B109 の更新は C・D2・E・A2。

## ### B72-B162 作業で新しく見つけたもの

取り除いた部分: 「- B135〜B159 (v2.2.0):」

元の行:

> - B135〜B159 (v2.2.0): `history/v2next_P1.md`〜`P5.md` の各担当の「K64 の誤りの候補」。B135・B136・B137 は P1、B138・B139 は P2、B140〜B145 は P3、B146〜B150 は P4、B151〜B159 は P5。B42・B108・B111・B112・B116・B128・B129 の更新は P5・P3・P2・P1。B128 は P2 が原因を特定した (SB の式の誤り) ので新しい番号は付けていない。P4 の「DAPWIDTH の窓を pysad が追っていなかった」は pysad 自身の不具合で、B 番号を付けていない (docs/47 (c))。

## ### B72-B162 作業で新しく見つけたもの

取り除いた部分: 「- B160 (v2.3.0):」

元の行:

> - B160 (v2.3.0): `history/v23_T1.md` §3.2・§6 (K1≠0 の BEND のオフモメンタムの光学。gs・解析式・elegant との比較、[`example/v2_3/T1_gs_bend_k1.py`](../example/v2_3/T1_gs_bend_k1.py))。U1 で K1 の暗黙の定義 (誤りではない) に再分類: `history/v23_U1.md`。

## ### B72-B162 作業で新しく見つけたもの

取り除いた部分: 「- B161 (v2.3.0):」

元の行:

> - B161 (v2.3.0): `history/v23_U1.md` §8 (CALC の `tbendecorr` の dk と kick の導関数の xr² の項の不整合。疑い)。

## ### B72-B162 作業で新しく見つけたもの

取り除いた部分: 「- B162 (v2.3.0):」

元の行:

> - B162 (v2.3.0): `history/v23_T2.md` (Touschek の損失率の式 `touckf` の適用範囲。独立な式 `touschek_piwinski` との比較、elegant との比較、`example/v2_3/T2_*`)。

## ### B72-B162 作業で新しく見つけたもの

取り除いた部分: 「- B163 (v2.4.0):」

元の行:

> - B163 (v2.4.0): `history/v24_RAD.md` (K64 の tcod が収束せず、gs が収束していない閉軌道から Emittance[] を計算する。DDBA_15 の εx の gs との 0.33% の差の原因)。

## ### B72-B162 作業で新しく見つけたもの

取り除いた部分: 「(v2.4.0 担当 D3、B165」

元の行:

> - B230〜B242 (v2.4.0 担当 D3、B165・B167 の続きと B145・B146): `history/v24_D3.md`。B230・B231 は QUAD・MULT の本体の B∥ (SEXT・ステアリングは元から正しい)、B232 は RF の同期位相 (エネルギー損失・電荷)、B233 は PhotonList の β、B234 は SOL 内の四極の B∥、B235 は追跡の BEND のスピンの刻みの 1 次誤差、B236 は B145 の再分類 (誤りではない)、B237 は `Bend.transfer_matrix` (pysad 自身の誤り)、B238 は SOL 内の加速 + POL、B239 は短い磁石・硬い端の放射の限界 (直さない)、B240 は Emittance[] の 6 極の放射・スピン、B241 は追跡の MULT の ROTATE のスピン、B242 は端の ψ の (1+δ) 規格化。

## ### B72-B162 作業で新しく見つけたもの

取り除いた部分: 「(v2.4.0 担当 E)」

元の行:

> - B180〜B199 (v2.4.0 担当 E): 組込の数学関数・特殊関数・リスト関数の監査 (mpmath との 4000 点超の比較と、gs との 600 式の照合)。記録は `history/v24_E.md`、試験は `pysad/tests/test_v24_E_audit.py` (参照値 `specfn_audit_ref.json`)。**複素の次数の `PolyLog`・`ExpIntegralE` (従来は SADError・未評価。gs は対応していて `PolyLog[1.5+I,.3]` = .3270−.0262 I) は k64・v2 とも実装 (|Log z| ≤ 5、`specfn_cx2.polylog_cs`・`expn_cn`)。** 「確度」は mpmath (collective312 の venv、`mp.dps = 40`) の主値と一致するかで判定。v2 の主値は Mathematica と同じ切断線の取り方 (ArcSin[x>1] は下側など)。

## ### B72-B162 作業で新しく見つけたもの

取り除いた部分: 「(従来は SADError・未評価。gs は対応していて」

元の行:

> - B180〜B199 (担当 E): 組込の数学関数・特殊関数・リスト関数の監査 (mpmath との 4000 点超の比較と、gs との 600 式の照合)。記録は `history/v24_E.md`、試験は `pysad/tests/test_v24_E_audit.py` (参照値 `specfn_audit_ref.json`)。**複素の次数の `PolyLog`・`ExpIntegralE` (従来は SADError・未評価。gs は対応していて `PolyLog[1.5+I,.3]` = .3270−.0262 I) は k64・v2 とも実装 (|Log z| ≤ 5、`specfn_cx2.polylog_cs`・`expn_cn`)。** 「確度」は mpmath (collective312 の venv、`mp.dps = 40`) の主値と一致するかで判定。v2 の主値は Mathematica と同じ切断線の取り方 (ArcSin[x>1] は下側など)。

## ## 6. 本レポートに含めなかったもの (参考)

取り除いた部分: 「v2.0.0 の作業で見つかった pysad 自身の不具合」

元の行:

>   v2.0.0 の作業で見つかった pysad 自身の不具合 (MULT のスキュー成分、`CalculateOptics` の TRPT、加速区間の分散、`LegendreP` の符号、EMITDIV ほか) は [`47_v2_changes.md`](47_v2_changes.md) の (d)。
