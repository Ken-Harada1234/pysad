# 20_translation_candidates から取り除いた履歴の記述

`docs/20_translation_candidates.md` を現在のスナップショットだけの解説にするため取り除いた文章。元の文脈 (節の見出し) つきで、文章はそのまま。

## # 20. pysad 翻訳候補と状況 (SAD の機能別)

取り除いた部分: 「最終更新: 2026-11-18 (v2.1.0 の後)」

元の行:

> 最終更新: 2026-11-18 (v2.1.0 の後)

## ## 1. FFS のコマンド (`tffsa.f` ほか)

取り除いた部分: 「v2 モードは MARK・CAVI も DX・DY に書き、corlength≠0 の規格化を直した。」

元の行:

> | DELCOR | `tcorr.f`、`tvcorr.f`、`tffsa.f:982` | 済み (k64 モードは gs と一致: CAVI は HARM・PHI に書く K64 のとおり、OCT も gs と一致。v2 モードは MARK・CAVI も DX・DY に書き、corlength≠0 の規格化を直した。docs/31) |

## ### 1.1 ON/OFF のフラグ

取り除いた部分: 「光学・Emittance への空間電荷は v2.1.0 で実装 (Emittance[]」

元の行:

> | WSPAC / SPAC / SELFCOD / PSPAC | k64 モード: 受理のみ (空間電荷は追跡に入らない)。**v2 モード: `TrackParticles` に入る** (SPAC = DRIFT・MULT の RADIUS、PSPAC = SPCH 要素、WSPAC = ビーム行列のガウス分布の力、SELFCOD)。光学・Emittance への空間電荷は v2.1.0 で実装 (Emittance[]・EMIT は両モード、CALC・CALC6D は v2 だけ。docs/07)。docs/08 §9.7。tspac.f の bbkick1 と tintrb.f の twspu は `GaussianCoulomb(Fitted)` として翻訳済み (docs/32) |

## ## 2. SADScript の関数

取り除いた部分: 「複素 FindRoot は v2 モードで実装 (gs が壊れた値: B61。k64 モードは TypeError)。複素の PolyLog・HypergeometricU・Zeta・HurwitzZeta は v2.1.0 で両モードで実装」

元の行:

> | PolynomialFit / Fit / FitEmit / FitGaussian / Spline / NIntegrate / Integrate / FindRoot / D | P (Fit.n、D.n) | | 済み (`fitroot.py`。gs と 5e-10〜1e-14。ベクトル独立変数は対応 (gs が誤値: B52)、複素 FindRoot は v2 モードで実装 (gs が壊れた値: B61。k64 モードは TypeError)。複素の PolyLog・HypergeometricU・Zeta・HurwitzZeta は v2.1.0 で両モードで実装 (`specfn_cx.py`。gs が誤る点は B118〜B120。docs/33)。NIntegrate は複素数の被積分関数も可。D・Integrate は D.n を動かす。docs/33・34) |

## ## 2. SADScript の関数

取り除いた部分: 「一部対応 (v2.1.0: BEND・QUAD」

元の行:

> | PHOTONS (`PhotonList`) | | `tftrack.f:144,315`、`photons.f` の kphlist | 一部対応 (v2.1.0: BEND・QUAD・MULT・SEXT・OCT・磁場つきドリフト・BEAMBEAM + BEND の光子に対応 (k64 は gs と全列一致)。SOL 領域内・COORD・TCAVI・MAP のあるラティスの光子は未対応。RADLIGHT の軌道と RadiationField・RadiationSpectrum は docs/45) |

## ## 3. 要素・物理モジュール

取り除いた部分: 「v2.1.0 で SOL 領域内の四極の DX・DY の解析微分 (`match6d.sol_dcod`) を足した。」

元の行:

> | 幾何 (LINE["GX"/"GY"/"GZ"/"CHI1".."CHI3"]、fit の GX) | `tfgeo.f` (411 行)、`tdgeo.f` | 済み (geo.py、4.8e-13。fit の幾何の微分は、DRIFT の L と BEND の ANGLE・L だけ解析微分 (tdgeo)、ほかは中心差分。v2.1.0 で SOL 領域内の四極の DX・DY の解析微分 (`match6d.sol_dcod`) を足した。CHI3 の微分は K64 でも領域の終端のリセットと合わず fit は収束しない: B124) |

## ## 3. 要素・物理モジュール

取り除いた部分: 「wake の重み (B84) は v2 で直した |」

元の行:

> | WAKE | `tffswake.f`、`txwake.f` | 追跡と光学 (WAKEOPT) が gs と一致 (docs/23・38)。旧式の WAKE コマンド (tfwake.f) は K64 が WAKECOMMAND.n (配布に無い) を呼ぶので k64 モードは同じ経路。v2 モードは `WAKE … ENDWAKE` の表を読む (B66)。wake の重み (B84) は v2 で直した |

## ## 5. 見送り・未対応とその理由

取り除いた部分: 「([`pysad/naff.py`](../pysad/naff.py)、v2.1.0) で代用」

元の行:

> | `wjfit`・`wfres`・`wfbin`・`wiord`・`wjdraw` | MAIN の TRACK の Fourier 解析で、FFS から使えない。gs も SIGSEGV で参照が無い。チューンの抽出は NAFF ([`pysad/naff.py`](../pysad/naff.py)、v2.1.0) で代用 |

## ## 5. 見送り・未対応とその理由

取り除いた部分: 「光学・Emittance への空間電荷は v2.1.0 で実装 (追跡の WSPAC」

元の行:

> | `tspac.f` の追跡側 (空間電荷) | **v2 モードで実装済み** (SPAC・PSPAC・WSPAC。docs/08 §9.7)。k64 モードは入れない (K64 の SPAC は B89・B90・B92 の誤りを持つ)。光学・Emittance への空間電荷は v2.1.0 で実装 (追跡の WSPAC のビーム行列は空間電荷なしのまま) |

## ## 5. 見送り・未対応とその理由

取り除いた部分: 「**実装済み** (v2.1.0。両モード。」

元の行:

> | 複素引数の PolyLog・HypergeometricU・Zeta | **実装済み** (v2.1.0。両モード。`specfn_cx.py`。独立な参照と 1e-10〜1e-11。gs の誤りの点でも真値) |

## ## 5. 見送り・未対応とその理由

取り除いた部分: 「`Parallelize` は v2.1.0 で POSIX は fork で並列 (Windows は逐次)。`$Context` の代入・`WithContext`・VISIT/BYE も v2.1.0 で実装 |」

元の行:

> | `System`・`Fork`・`Shared`・`BeginPackage` の文脈・`NOKEEPEXP`・`Out[n]:=` | `System`・`Fork` (POSIX)・文脈は **v2 モードで実装済み** (docs/04)、`Shared`・`NOKEEPEXP`・RENUMBER の INDEX は両モード。`Out[n]:=` の gs の表示 (前の値の繰り返し) は K64 の癖として写さない。`Parallelize` は v2.1.0 で POSIX は fork で並列 (Windows は逐次)。`$Context` の代入・`WithContext`・VISIT/BYE も v2.1.0 で実装 |

## ## 6. 未解決の差 (pysad と gs が合わない。原因は確定していない)

取り除いた部分: 「(BEND の `F1=0.1 FRINGE=3` を RAD ON で追跡 (相対 0.7) と L≠0 の CAVI の DX・DY (相対 1e-3) は v2.1.0 で解決: B113・CAVI の TENT・TEXIT)、SOL」

元の行:

> (BEND の `F1=0.1 FRINGE=3` を RAD ON で追跡 (相対 0.7) と L≠0 の CAVI の DX・DY (相対 1e-3) は v2.1.0 で解決: B113・CAVI の TENT・TEXIT)、SOL (BZ≠0) 内の四極 + RAD の追跡で現在の gs が保存済みの参照 (10-14) と別の値を返す (B68)、

## ## 6. 未解決の差 (pysad と gs が合わない。原因は確定していない)

取り除いた部分: 「は v2 の作業で原因が分かり解消した (K64 の Emittance[] の誤りと MomentumCompaction の誤り: B75・B77。v2 モードで直した。docs/07)。「EMIT のあとに `USE` すると gs は EMIT の効果が追跡に残る」差は v2.1.0 で解決した (K64 の仕様: `trf0` を USE で持ち越す)。」

元の行:

> なお、「RADINT と EMIT の Jx の 2.4% の差」と「`BUCKET HIGHT` と `BucketHeight` の差」は v2 の作業で原因が分かり解消した (K64 の Emittance[] の誤りと MomentumCompaction の誤り: B75・B77。v2 モードで直した。docs/07)。「EMIT のあとに `USE` すると gs は EMIT の効果が追跡に残る」差は v2.1.0 で解決した (K64 の仕様: `trf0` を USE で持ち越す)。
