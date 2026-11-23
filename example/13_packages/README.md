# 13_packages: SAD の Packages/*.n

`correction.n`・`bump.n`・`optimize.n`・`D.n`・`rational.n`・`tableform.n`・`Couple.n`・`LDSolve.n`・`Optics.n` を使う例。軌道補正とバンプの例は、誤差入りのラティス `00_lattice/errors_ring_pkg.sad` (FODO_ring と同じ光学に、BPM・ステアラ・四極の据え付け誤差を付けたもの。`errors_ring.sad` には BPM・ステアラが無いので、別名で作った) を使う。

一致の記号: ○ 一致 / △ 差あり (md に理由) / ×gs 失敗 / ×pysad 未対応・失敗

| 例 (NAME.sad・.gs.txt・.pysad.txt・.md) | 内容 | gs と pysad |
|---|---|---|
| [fn_CorrectOrbit_Monitor_Steer](fn_CorrectOrbit_Monitor_Steer.md) | Monitor・Steer・CorrectOrbit (垂直・水平) | ○ |
| [fn_MakeBump](fn_MakeBump.md) | Bump・MakeBump (垂直の閉バンプ) | ○ (目標の 1 mm に一致、12 本) |
| [fn_OptimizeOptics](fn_OptimizeOptics.md) | OptimizeOptics (チューンの調整) | ○ (Output->0 で反復を出さない) |
| [fn_D_Integrate](fn_D_Integrate.md) | D (記号微分)・Integrate | ○ (Integrate の記号積分は両方とも評価しない) |
| [fn_Rational](fn_Rational.md) | Rational・ContinuedFraction・SmallPrime | ○ (和・積を評価する) |
| [fn_TableForm](fn_TableForm.md) | TableForm | ○ |
| [fn_Couple](fn_Couple.md) | Couple.n (ElementValues) | ○ (ElementValues に対応) |
| [fn_LDSolve](fn_LDSolve.md) | LDSolve | ○ |
| [fn_Optics_TransferMatrix](fn_Optics_TransferMatrix.md) | TransferMatrix・BeamMatrix・OrbitGeo (Optics.n) | ○ |
| fn_CorrectOrbit_Micado | CorrectOrbit の Micado (`Solver->{"Micado",n}`)。12 セルの FODO (3 GeV、周長 52.8 m)、BPM 24 台、水平・垂直ステアラ各 12 台、四極の据え付け誤差 DX・DY を式 (SIN, PI) で与える | k64 モードで gs と 3e-8 以内 (どちらも n に対して単調に下がらない)。v2 は n が増えるほど単調に下がる (K64 の `pmicad` の誤り) |
| fn_CorrectOrbit_Micado_variants | Micado の変種 (水平・2 平面同時・`SVDmode->1`・`Bump->`)。ラティスは fn_CorrectOrbit_Micado と同じ | k64 モードで gs と 5e-6 (B129 の K64 の誤りまで再現)。v2 は n で単調に下がる。試験 `test_v2_P5_micado_fitp_oy.py` |
| fn_CorrectOrbit_Bump | MakeBump と CorrectOrbit の `Bump->`。ラティスは fn_CorrectOrbit_Micado と同じ | MakeBump は gs と 1e-12。`Bump->` は K64 の誤りで gs も pysad もエラー |
| fn_CorrectLocal | CorrectLocal (局所バンプで領域 {BA3, BA8} の外の軌道は変わらない)。12 セルの FODO | 24 BPM の DY が gs と 1e-14 以内、補正量 12 個も同じ (v2・k64 とも) |
| fn_PeriodicR | Optics.n の PeriodicR (4×4 の一周行列から結合行列 R と x・y の面の行列 A・B) | 6 行列 × 12 成分が gs と 3e-16 以内 (NaN になる 2 行列も同じ) |
| fn_CSROYMatrix | CSROYMatrix を、手で作った正しい形の実数リスト (36 要素) で | gs と相対 1e-12〜4e-20 で一致。`CSRSetupOY` の出力の 4 番目の Head が gs では数、pysad では List (B60 の原因) |
| fn_OYStability_pipeline | OYStability の後半 (CSROYMatrix → Eigensystem → 成長率)。CSRSetupOY の 6 つのリスト (nr=256、nj=6、nphi=3、dz=2e-5) を書き込む | gs と 1e-12 (安定 −damp)。通しの OYStability・OYScan は gs が動かず未確認 |
| fn_OYStability_pipeline_unstable | 上と同じ。不安定な強度 (aw=1e-3) の場合 | gs と 1e-12 (成長率 0.0272796226) |
| fn_packages_compare | Packages/functions.n・Misc.n・element.n・beamline.n の主な関数を gs と比べた小さな台本 | gs の出力は `.gs.txt`、pysad の出力は `.pysad.txt` |

`SetCoupledElements`・`CorrectLocal`・`PeriodicR` などは試していない。
