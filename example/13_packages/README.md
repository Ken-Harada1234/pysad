# 13_packages: SAD の Packages/*.n

`correction.n`・`bump.n`・`optimize.n`・`D.n`・`rational.n`・`tableform.n`・`Couple.n`・`LDSolve.n`・`Optics.n` を使う例。軌道補正・バンプ・最適化の例は、誤差入りのラティス `00_lattice/errors_ring_pkg.sad` (FODO_ring と同じ光学に、BPM・ステアラ・四極の据え付け誤差を付けたもの。`errors_ring.sad` には BPM・ステアラが無いので別名で作った) を使う。

一致の記号: ○ 一致 / △ 差あり (md に理由) / ×gs 失敗 / ○ (2026-11-17 に ElementValues を対応)・失敗

| 例 (NAME.sad・.gs.txt・.pysad.txt・.md) | 内容 | gs と pysad |
|---|---|---|
| [fn_CorrectOrbit_Monitor_Steer](fn_CorrectOrbit_Monitor_Steer.md) | Monitor・Steer・CorrectOrbit (垂直・水平) | ○ |
| [fn_MakeBump](fn_MakeBump.md) | Bump・MakeBump (垂直の閉バンプ) | ○ (2026-11-17 に修正。目標の 1 mm に一致、12 本) |
| [fn_OptimizeOptics](fn_OptimizeOptics.md) | OptimizeOptics (チューンの調整) | ○ (2026-11-17 に Output->0 で反復を出さないよう修正) |
| [fn_D_Integrate](fn_D_Integrate.md) | D (記号微分)・Integrate | ○ (Integrate の記号積分は両方とも評価しない) |
| [fn_Rational](fn_Rational.md) | Rational・ContinuedFraction・SmallPrime | ○ (2026-11-17 に和・積を評価するよう修正) |
| [fn_TableForm](fn_TableForm.md) | TableForm | ○ |
| [fn_Couple](fn_Couple.md) | Couple.n (ElementValues) | ×pysad 未対応 |
| [fn_LDSolve](fn_LDSolve.md) | LDSolve | ○ |
| [fn_Optics_TransferMatrix](fn_Optics_TransferMatrix.md) | TransferMatrix・BeamMatrix・OrbitGeo (Optics.n) | ○ |

`SetCoupledElements`・`CorrectLocal`・`PeriodicR` などは試していない。
