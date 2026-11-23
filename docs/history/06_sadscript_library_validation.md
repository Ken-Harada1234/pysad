# 06_sadscript_library の検証の記録 (gs との突き合わせ。旧 33・34・32 から移したもの)

## 旧 33 §1
「K64 の実装 (gamma.f、hg.f、Packages/D.n) の値と、gs で 458 式を突き合わせて一致 (`tests/specfn_ref.json`、`tests/test_specfn.py`)。」(本文には 1 行で残した)

## 旧 33 §2
- 「gs で `Fit` の導関数が解析的な場合との差は線形な式で相対 1e-9、FindRoot の収束値で 1e-14。FitGaussian の出力は gs と 5e-10 で一致。」
- 「NIntegrate … 収束した値は gs と 1e-14 で合う。… (gs と虚部 1e-15 で一致。実部は 1e-16 の雑音)。」
- 「gs と全桁一致を確認 (2026-11-16 夜): `Spline[{{x,y},…}]` (自然スプライン)、`Sort[list, Greater]`、`Scan`。」
- 「`FitGaussian` … gs と同じ `List` 5・キャンバス 74 項目」「`PolyLog` … (gs と相対 1e-14)」
- 複素の特殊関数: 「gs の値を調べた点は乱数 230 点 (Zeta 40、HurwitzZeta 80、PolyLog 60、HypergeometricU 50) と 26 点の格子。gs が正しい 18 点は gs と相対 2e-13。」

## 旧 33 §3
「Amoeba の移植。gs と 1e-13 (Rosenbrock、500 反復) で一致。」「`OptimizeOptics` は … gs と一致する (docs/34)。」

## 旧 33 §5
「(ARRAY、Len、StringInsert、Diagonal、Expand は gs と一致)」「`SingularValues` `MatrixBalance` (gs と一致)」

## 旧 32 GaussianCoulomb
- 「15 点 (軸上・x 軸上・y 軸上・負の座標・σx<σy・σx≈σy (1e-5)・大きい引数 (漸近展開)・小さい引数) の実部・虚部が gs と 12 桁で一致 (表示の 12 桁が全部同じ)。」
- 「gs と相対 1e-9 以内で一致 (試験 `test_gaussian_coulomb_fitted_matches_gs`)。… 0 とみなして gs と一致した。」

## 旧 34 (SAD の Packages)
- 「`MakeBump` は gs と同じバンプ (ybump=.0010000141、使うステアラ 12 本)」「QD の K1 は -0.56、NX=6.788257 (gs と一致)」「`Rational` … gs と同じ」
- gs との比較 (`tests/test_corr.py`: 12 セル、24 BPM、24 ステアラ、DX/DY の誤差つき QUAD): 「応答行列 `ResponseMatrixG` の行、CorrectOrbit の補正量 (Y、X、Micado 6 本) は gs と全桁一致 (補正量 1e-7 相対)。OptimizeOptics の結果 (K1、チューン) も一致。」「ResponseMatrix の TRPT … gs と全桁。」
- LDSolve・Optics.n (2026-11-16 夜に gs と比較): 「Optics.n の `TransferMatrix`・`BeamMatrix`・`NormalCoordinates`・`OrbitGeo`: gs と表示桁で一致 (gs は 1e-16 の雑音の要素を出し、pysad は 0)。LDSolve の終点 (cos 1, −sin 1) と 3 通りの呼び出しが gs と一致 (試験 `test_evalfixes.py` の 2 件)。`PeriodicR` は gs と比べた。」
- 「小さな台本の gs との突き合わせ (小さな台本を gs と行ごとに突き合わせた。k64 モード)。数値は 1e-12 を超える差だけを「差」とした。」の表の「比べた関数」「一致の状態」の列:

| ファイル | 比べた関数 | 一致の状態 |
|---|---|---|
| `functions.n` | `Log2$`、`InverseErfc`、`Names`、`RotateRight`・`RotateLeft`、`Union` (`SameTest`)、`ListQ`、`MatrixQ`、`Diagonal`、`Count`、`MemberQ`・`FreeQ` (レベルつき)、`Nest`・`NestList`・`Fold`・`FoldList`、`Through`、`FractionalPart`、`Positive` 系、`IntegerQ`、`Subtract`・`Divide`、`Xor`、`Append`・`Prepend`、`HeldPart`、`Evaluate`、`StringQ`・`Characters`・`StringDrop`・`StringInsert`、`Timing`、`For` | `Log2$`・`Union` (`SameTest`)・`Names` 以外は一致 |
| `Misc.n` | `NearlySameQ`、`FixedPoint`、`FixedPointList`、`FunctionDefinition` | `FunctionDefinition` 以外は一致 |
| `element.n` | `Element` の読み出し (LENGTH・NAME・POSITION・COMPONENT・DEFAULT・VALUE・TYPE・K1・L・ANGLE・GEO・DK・EXPAND)、代入 (`=`・`+=`・`-=`・`*=`・`/=`・リスト形・`VALUE`)、`PROFILE`、`SetElement` | 読み出しは一致 |
| `beamline.n` | `BeamLine` の演算、`ExtractBeamLine`、`WriteBeamLine`・`PrintBeamLine` (FFS 形式・MAIN 形式・`Name` のリスト)、`Fill$String`、`Remove$Paren` | `BeamLine`・`ExtractBeamLine`・`WriteBeamLine` は全行一致 |
| `emit.n` | `Emittance[]` の戻りの規則の名前・並び・値、`SymplecticJ`・`SymplecticInverse`、`TouschekF`、`TouschekLifetime` の未定義のメッセージ、`Beamstrahlung*` 4 つ | 値は 1e-12 以内 (3 GeV、周長 41.8 m 余りの FODO 8 セル。Emittances の y は gs が 2e-38、pysad は 0 の雑音) |
| `dapert.n` | `DynamicApertureSurvey` (`Output`・`Axes`・`ExpandElementValues`、不正な `Axes`)、`AccelerateParticles`、`TrackParticles` | `AccelerateParticles` は 1e-14、`TrackParticles` も一致、`DynamicApertureSurvey` の Score は一致 |
| `Optics.n` の `PeriodicR` | 4x4 の一周行列 6 通り (結合あり・なし・半整数に近い・NaN になる 2 通り) | R・A・B の 12 成分が 3e-16 以内 |
| `correction.n` の `CorrectLocal` | 12 セル FODO (3 GeV、周長 52.8 m、BPM 24 台)、領域 {BA3, BA8}・垂直 | 補正量 12 個と 24 BPM の DY が 1e-14 以内。領域の外の軌道は 1e-11 以内で変わらない |
| `correction.n` の `CorrectOrbit` の Micado | 24 セル FODO (BPM 48 台、垂直ステアラ 24 台)、`Solver->{"Micado",n}` | k64 は gs と 1e-6 以内で一致 |

- 「`dapert.n` ほか gs との主な差の一覧」の表 (本文 §5.4 の「残った差」と重複するので本文では省いた):

| 項目 | gs | pysad |
|---|---|---|
| `DynamicApertureSurvey` の `TURNS` の既定 | −1 | −1 |
| `DAPWIDTH` の既定と窓 | 7。窓の追跡順 | 51。全点を追跡 (§未対応・限界) |
| `Emittance->False` | DampingRate 以降の規則を返さない | 同じ |
| `Union[…, SameTest->f]` | 先頭の要素だけ (K64 の誤り) | k64 は同じ、v2 は真値 |
| `Names["組み込みの名前"]` | `{}` (大域に作った記号だけ) | k64 は同じ、v2 は組み込みも含める |
| `WriteBeamLine[…, Name->"L2"]` | エラー (K64 の誤り) | k64 は同じ、v2 は受ける |

## 旧 34 の未対応・限界
「`SingularValues`・`MatrixBalance` … 特異ベクトルの符号は各 u の先頭成分が負になるよう揃えた (gs と一致を 3 行列で確認)。」

## 旧 32 BBBrem1 (14_beambeam 側の担当へ)
「検証 (`test_bbbrem.py`): `SeedRandom[3]` のあと 6 事象、sqrt(s)=10,20,30,40 GeV の 4 事象で、乱数の列も含めて相対 1e-15 (最後の桁の丸め) で一致。」
