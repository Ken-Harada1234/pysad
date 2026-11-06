# SAD の解析微分 (tffsqu / qdcell) の pysad への移植 — 作業ログと解説

作成: 2026-10-15。対象: `pysad/match6d.py` (新規)。既存ファイルは変更していない。

## 0. 結果の要約

- SAD の既定経路 (`FFS$NumericalDerivative=False`) の解析微分 `tffsqu → qdcell → qdtwis → qdquad/qdbend`,
  `qdmdia`, `qgettru/tftmatu/qtentu/qchg` を `pysad/match6d.py` に直訳した。
- match6.py の数値微分 (`qu = (ddf1−ddf2)/2/dvkc/wvar`) と全成分を照合: Level 13B (RING, Q* 7 変数 × 6 条件) で
  最大相対差 2.0e-8、ESRF ARCA (Q* + BPI6/BPI7 K1 の 9 変数 × 14 条件、EX/EPX/EY/EPY/NX/NY/AX/AY/BX/BY) で 4.1e-7。
  (対称点でゼロになる成分と、閉軌道探索の収束誤差が数値微分に現れる EPX 行は雑音として除外。§5)
- `AnalyticMatcher` (match6.Matcher の go() の微分ブロックだけを置き換えたサブクラス) で
  - Level 13B: SAD と同じ 2 反復、Matched 3.3724E-10、K1 は 10 桁一致。
  - esrf/matching (matching1): SAD の反復列 `3 2.9742E+07 (NEWTON) 0.4561 / 4 7.2238E+06 / 5 5.1917E+05 / 6 7.340 /
    7 0.1633 / Matched 8.7338E-09` を**そのまま再現**し、9 変数の K1 が 12〜13 桁一致 (§6)。
- 途中で見つけた SAD 側の事実 (pysad 全体に関わる):
  1. SAD の `sinc(x)` は `x·cos x − sin x` (`sinhc` は `x·cosh x − sinh x`)。sin x / x ではない (src/tfloor.f)。
  2. この ESRF スクリプトでは `cell calc emit` の Emittance[] が閉軌道を見つけられず **EMITX=EMITY=NaN** になる。
     tffscalc の重み計算 (`em=max(emminv,|emx|+|emy|)` …) は gfortran の MAX/MIN が NaN を無視する性質により
     `emxx = emyy = MINCOUP/(1+MINCOUP)·emminv = 0.01/1.01·1e-15 ≈ 9.9e-18` になり、EX/EPX の重みが ~1e6〜1e7 になる。
     SAD の `calc` 直後の MatchingResidual 4.1205E+07 はこの重みで説明がつく (pysad で emx=emy=9.9009901e-18 とすると
     4.1205E+07 を再現)。pysad の ffscmd は EMIT 失敗後も EMITX を MARK の値 (1.8e-9) のままにするので、ここが
     esrf/matching の「別解に着く」真の原因だった (微分の種類ではない)。
  3. SAD の数値微分経路 (`FFS$NumericalDerivative=True`) は同じ問題で 47 反復 (DESCEND 混じり) になる。
     SAD の qcell は RAD/RADCOD ON で閉軌道が収束しきらず数値微分が荒れるため。pysad の数値微分 (match6) は
     解析微分と同じ反復列を出す (pysad の光学は決定的なので)。

## 1. SAD 側の呼び出し階層と qu の定義

```
tffsmatch (src/tffsmatch.f 395-403)
  call tffsqu(nqcol,nqcol1,nvar,nqumax,ifquw,ifqu,free,irtc)      ← 解析微分 (条件 1..nqcol1)
  rlist(kqu:kqu+nqcol1-1) = rlist(...)*wiq(1:nqcol1)/wvar(kc)       ← 重みと変数スケール
  (nqcol > nqcol1 なら FitFunction の行を v+eps1/wvar の片側差分で: (df−df1)/eps1)   420-460 行
tffsqu (src/qdcell.f 1253-)
  for k in 1..nlat-1 (free なファミリーの要素):
    for ii in 変数: coup = couple(k)*wk (通常 1)
      for kq in 条件: lp = lfp(1,kq) (= kpe), lfp(2,kq) (= kpb, 2 点条件のみ)
        qdcell(dtwiss, k, lp, idp, iv, ...)      ← 要素 k・変数 iv による位置 lp の Twiss 微分
        rlist(kqu) += coup*dtwiss(kf)  (lfp(2) 側は −coup)
qdcell (src/qdcell.f 577-)
  qdtwis(...)                                   ← 輸送 (k 入口固定) の微分
  cell なら周期解の補正 (後述)
  dgammag                                       ← GMX/GMY, BMAG
qdtwis (src/qdtwis.f 1-)
  要素型と iv で分岐: qddrif / qdbend (ANGLE, K1) / qdquad (K1, ROT) / qdthin / qdtrans (数値) / qdini (MARK 先頭)
  tftmatu(k+1 → l), qgettru(k → l, disp, fold) で l まで運び β, α, ψ, η, R の微分に変換
```

qu の定義 (SAD も match6 も同じ): `qu[条件, 変数] = ∂f/∂v · wiq[条件] / wvar[変数]`。
match6 の数値微分は `df = 目標 − f` を使うので `(ddf1−ddf2)/2/dvkc = (df(v−dv) − df(v+dv))/2dv = +∂f/∂v` と符号が合う。
f は条件関数の「比較に使う形」: BX/BY は log β (→ dtwiss(bx) = dβ/β)、AX/AY は atan α (→ dtwiss(ax) = dα/(1+α²))、
NX/NY は ψ [rad]、EX.. は η そのもの。tdfun1 (qdcell.f) の形と対応している。

## 2. 数式 (qdtwis + qdcell)

記号: 位置 l の Twiss (β_l, α_l, ψ_l)、一周の位相 μ = ψ(nlat) − ψ(1)、要素 k の微分行列 dM_k = ∂M_k/∂v。
SAD の 4×5 行列は `trans(1:4,1:4)` が (x, px, y, py) の転送行列、`trans(:,5)` が分散列。

### 2.1 要素の微分行列 (qdquad, qdbend)

`dtrans = dM_k · S_k`。S_k は k 入口の R 変換 (qtentu: 正規モード → 物理、R=0 なら単位行列) で、5 列目は η_k
(qtentu の disp=.true.) なので `dtrans(:,5) = dM_k·η_k + (閉軌道に依存する項 a16, a26, …)`。
QUAD K1 (qdquad, iv=2): k = K1/(L·pr), φ = √|k|·L, dφ/dK1 = φ/(2K1) として
```
a11 = −sinφ·dφ,  a12 = (cosφ − sinφ/φ)·dφ/√k,  a21 = −(cosφ + sinφ/φ)·dφ·√k   (集束面; 発散面は sinh/cosh)
```
これは M = [[cosφ, sinφ/√k],[−√k sinφ, cosφ]] の K1 微分そのもの。
BEND K1 (qdbend, iv=8): 複合関数セクター (aind = K1·ρ0/φ) の x/y 行列の微分。`a12 = ρx/φx·sinc(φx)·dφx` の
`sinc` は SAD の定義 `sinc(x) = x·cos x − sin x` なので a12 = ρx(cos φx − sin φx/φx)·dφx = ∂(ρx sin φx)/∂K1。
(最初 sin x / x と誤訳して BPI6/BPI7 列が 10〜300% ずれた。§6 の反復列再現が崩れたことで発見。)
端部角 E1/E2 は `tan(ψ1·φ0)` (ψ1 は ANGLE に対する比) で ak1 = −tan/ρ の薄レンズとして前後に掛かる。

### 2.2 輸送の微分 (qdtwis)

```
trans  = tftmatu(k+1 → l)              物理座標の転送行列 (Twiss から qgettru で作り、入口側 R を右から掛ける)
trans1 = trans(:,1:4) · dtrans         = d(M_{k→l}) (k の入口 Twiss を固定したときの微分)
trans  = qgettru(k → l, disp=.true., fold=.true.)   正規モード行列 X (k≥l なら一周をまたぐ: dpsi += μ)
x11.. = trans/√(1−det R_l),  dx11.. = (trans1 − x·d√)/√
```
β_l = x11²β_k − 2x11x12α_k + x12²γ_k 等から
```
dβ_l/β_l = 2(−(dx11 x12 + x11 dx12)α_k + dx11 x11 β_k + dx12 x12 γ_k)/β_l
dα_l     = (dx11x22 + x11dx22 + dx12x21 + x12dx21)α_k − (dx11x21 + x11dx21)β_k − (dx12x22 + x12dx22)γ_k
dψ_{k→l} = (dx12(β_k x11 − α_k x12) − x12(β_k dx11 − α_k dx12)) / (x12² + (β_k x11 − α_k x12)²)
dη_l     = trans1(1,5) (+ R の変化の項)、 dR_l = (trans1 の 3,4 行と x の組合せ: qdtwis.f 148-200)
```
dtwiss(ax) は最後に (1+α_l²) で割って atan の微分にする。l = nlat なら TRX/TRY (trace) の微分も作る。

### 2.3 周期解の補正 (qdcell, CELL のとき)

位置 l の一周行列 M_l = S_l⁻¹·diag(X, Y)·S_l (X = [[cosμ+α sinμ, β sinμ],[−γ sinμ, cosμ−α sinμ]])。
要素 k だけが変わるときの一周行列の微分は dM_l = d(M_{k→l})·M_{l→k} で、これは §2.2 の輸送微分
(dβ/β, dα, dψ_{k→l}) を使って
```
dx11 = (dβ/2β)x11 + (−sinμ + α cosμ)dψ,     dx12 = (dβ/2β)x12 + β cosμ dψ
dx21 = −x21(dβ/2β) + (−(α sinμ + cosμ)dα − (1+α²)cosμ dψ)/β,   dx22 = −(dβ/2β)x22 − sinμ dα − (sinμ + α cosμ)dψ
dM_l = dS⁻¹·D·S + S⁻¹·dD·S     (dS⁻¹ は輸送微分の dR から)
```
と書ける (B_l R(μ) B_l⁻¹ の形で検算済み: dB_l R B_l⁻¹ + B_l R' B_l⁻¹ dψ)。
次に qdmdia で一周行列の対角化行列 R の微分 dr1..dr4 (正常分岐: c² = ½ + ½|TrP−TrQ|/√((TrA−TrB)²), r = ±(M の非対角ブロック)/(c√…))、
それから `A = S·M_l·S⁻¹` の 2×2 ブロックの微分 da, db を
`dmd = dS·(M S⁻¹) + S·(dM S⁻¹) + S·(M dS⁻¹)` で作り、
```
tr A = 2cosμ      → dμ    = −(da11+da22)/(2 sinμ)
a11−a22 = 2α sinμ → dα_l  = (½(da11−da22) − α cosμ dμ)/sinμ
a12 = β sinμ      → dβ/β  = (da12 − β cosμ dμ)/(β sinμ)
```
位相 ψ_l (先頭からの位相) の微分:
- k < l: 要素は 1→l の間にある。l→終端の転送行列は固定で、その位相進み ψ_{l→end} = μ − ψ_l は l の Twiss の変化で
  `dψ_{l→end} = −sin ψ'[(cos ψ' + α sin ψ')dβ/β − sin ψ' dα]` (ψ' = μ − ψ_l) だけ変わるので
  `dψ_l = dμ + sin ψ'[(cos ψ' + α sin ψ')dβ/β − sin ψ' dα]`。
- k ≥ l: 1→l の転送行列が固定。`dψ_l = −sin ψ_l[(cos ψ_l − α sin ψ_l)dβ/β + sin ψ_l dα]`。
閉軌道: 正規モードで `dx = (I−X)⁻¹(dM·x + δx)` (detimx = 4 sin²(μ/2))。分散: `dη = (I−X)⁻¹·d(分散列)` に
R の変化 `ds·η_phys` を加える。GMX/GMY は dgammag。

### 2.4 tffsqu の組み立て

変数 (ファミリー, キーワード) に属する要素インスタンス k ごとに、条件の位置 lp について qdcell を 1 回計算し
(同じ lp の条件には使い回す: SAD の col フラグ)、`qu += coup·dtwiss(kf)`。2 点条件 (kp≠kp1) は
`dtwiss(kpe) − dtwiss(kpb)` (trdfun の ve − vb、log の場合も dβ/β の差で同じ)。LENG は要素長 L の変数だけ ±1。
FitFunction の行 (iqcol<0) は tffsmatch 側で片側差分。

## 3. 対応表

| SAD (src) | 行 | match6d.py | 備考 |
|---|---|---|---|
| qdcell.f `tffsqu` | 1253- | `tffsqu`, `analytic_qu` | coup=1 (COUPLE, 境界分数 wk は未対応)、CouplingMatrix 無し |
| qdcell.f `qdcell` | 577- | `qdcell` | tm/dtm/ts/dts/tds/da/db は行列積で書いた (SAD の展開式と代数的に同一、s の 0 成分の省略) |
| qdcell.f `dgammag`, `dbmag2` | 1227- | `dgammag`, `_finish_gamma` | BMAG の参照 Twiss (twiss(ll,−1,·)) は無いので BMAG は未対応 |
| qdtwis.f `qdtwis` | 1- | `qdtwis` | qddtwiss (閉軌道の数値伝搬 ctrans) は CELL では DZ/DDP 条件以外で呼ばれないので未移植 |
| qdtwis.f `qdini` (MARK 先頭の変数) | 311- | — | 未移植 (dtwiss=0) |
| qdtwis.f `qdtrans` (要素写像の数値微分) | 562- | — | 未移植。未対応の変数型は `AnalyticMatcher._numeric_columns` (match6 と同じ中心差分) |
| qdquad.f `qdquad`, `qdrotate` | 1-, 107- | `qdquad`, `qdrotate` | L=0 の QUAD (qdthin) は未移植 |
| qdbend.f `qdbend` | 1- | `qdbend` | ANGLE (iv=2), K1 (iv=8)。sinc/sinhc は SAD 定義 |
| qdmdia.f `qdmdia` | 1- | `qdmdia` | 縮退分岐も直訳。label 1300 の dr4 未設定 (SAD のバグ) は dr4=0 にした |
| qgettr.f `qgettru` | 21- | `qgettru` | gammab は一定 (gr=1) |
| tftmat.f `tftmatu` | 22- | `tftmatu` | |
| qtent.f `qtentu` | 125- | `qtentu` | |
| qchg.f `qchg` | 1- | `qchg` | |
| toplvl.f `tmultr45` | 1516- | `tmultr45` | |
| tfloor.f `sinc`, `sinhc` | 369-, 383- | `_sinc`, `_sinhc` | x cos x − sin x |
| tffs.f `tffitcode` | 104- | 定数 `AX..LENG` | utwiss の添字 |
| tffsmatch.f 320-470 (微分ブロック) | | `AnalyticMatcher._derivative`, `go` | go() は match6.Matcher.go の写し |

pysad 側のデータの取り方 (`DerivContext`, `build_utwiss`):
- `lattice.twiss_table()` の行 `rows[l-1]` = SAD の位置 l (要素 l の入口、rows[nlat-1] が終端 $$$)。
  utwiss(l, mfit) に詰め替える (`build_utwiss`)。
- 変数 (family, key) → SAD の iv: QUAD K1→2, ROTATE/ROTATION→4, BEND ANGLE→2, K1→8 (`_IV_TABLE`)。
- 要素パラメータの単位: K1 は SAD の積分値 (`match6._sad_get` で el.K1·L)、BEND の E1/E2 は pysad が rad で持つので
  ANGLE で割って SAD の比に戻す (parser.py の E1 変換が推測ベースなのはそのまま)、ROT は `ROTATION` 属性、
  K0 (BEND の追加曲げ) は `K0` 属性、phib0 = ANGLE + K0。
- ANGLE=0 の BEND (VD/HD) は parser が Drift にする (`_sad_type='BEND'`) → qdbend は al=0 で 0。

## 4. match6d.py の使い方

```python
from pysad.match6d import AnalyticMatcher, analytic_qu, compare_derivatives
m = AnalyticMatcher(lat)        # match6.Matcher と同じ API。go() が解析微分
qu, numeric_cols = analytic_qu(m, ids, rows, wvar)   # go() の数値微分ブロックの代替 (ids = dfs のキー列)
qa, qn, ids, worst = compare_derivatives(m)          # 数値微分との照合 (雑音成分は除外して最大相対差を返す)
```
`AnalyticMatcher.numerical_derivative = True` にすると match6 と同じ数値微分 (FFS$NumericalDerivative=True 相当)。
ffscmd で使うにはモンキーパッチ `pysad.ffscmd.Matcher = AnalyticMatcher` (run_matching_analytic.py 参照)。

## 5. 照合結果 (bench_match6d_check.py)

`PYTHONIOENCODING=utf-8 python -u bench_match6d_check.py`

| ケース | 条件×変数 | 最大相対差 | 雑音として除外 |
|---|---|---|---|
| Level 13B (RING; NX, NY, MKC BX/AX/BY/AY; free Q*) | 6×7 | 2.0e-8 | 14/42 (MKC の AX/AY: 対称点で 0。解析 ~1e-12、数値 ~1e-9) |
| ESRF ARCA (MKS BY, $$$ NX/NY, MKC AX/AY/EPX/BY, MKSFA1 AY/BY/EX/BX, MKSFA1.2 NX, MKSDC1 EY/EPY; free Q*, BPI6 K1, BPI7 K1) | 14×9 | 4.1e-7 | 45/126 (MKC AX/AY/EPX: 対称点で 0。EY/EPY: 全部 0) |

雑音判定: 数値微分を刻み dvkc と 2dvkc で計算し 1% 以上違う成分。EPX@MKC は重み 246 で数値側が 10〜40 (非重み 0.1)
出るが、これは閉軌道探索 (COD 反復) の収束誤差が差分で増幅されたもので、解析側は ~1e-3 (非重み 1e-5)。
相対差 1e-6 以下の判定は、対称でない成分すべてで満たされた (13B の最大 2e-8 は数値微分の打ち切り誤差の桁)。
qdquad の K1 微分行列と pysad Quadrupole.transfer_matrix の差分、qdbend の K1 微分行列と Bend.transfer_matrix の差分も
全成分一致 (sinc 修正後)。

## 6. esrf/matching の再現実験 (run_matching_analytic.py)

SAD 参照: <work-esrf>/matching1.sad.out (既定 = 解析微分)。追加で確認用に <work-esrf>/matching1_dbg.sad
(`cell calc emit` 後に EMITX/EMITY を Print、GO 前に `calc` と MatchingResidual を Print) と
matching1_dbgnum.sad (`FFS$NumericalDerivative=True`) を WSL の SAD で実行 (出力 *.sad.out、run_dbg.sh)。

SAD の出力: `EMIT0 {NaN,NaN,NaN}`、`RES0 41205475.56577`、反復 3 2.9742E+07 (0.4561) … 7 0.1633 → Matched 8.7338E-09。
数値微分版: 同じ 3〜7 の後 18 DESCEND 3.8496E-02 … 47 DESCEND 1.9948E-08 → Matched 9.2554E-09 (47 反復)。

pysad: `PYSAD_EMITX=9.9009901e-18 python run_matching_analytic.py` (GO 直前に EMITX=EMITY=9.9009901e-18 を差し込み、
SAD の NaN 重みと同じ emxx を作る):
```
       3        0   2.9742E+07  (NEWTON)   0.4561          9
       4        0   7.2238E+06  (NEWTON)    1.000          9
       5        0   5.1917E+05  (NEWTON)    1.000          9
       6        0    7.340      (NEWTON)    1.000          9
       7        0   0.1633      (NEWTON)    1.000          9
 Matched. ( 8.7338E-09)
```
| 変数 | SAD | pysad 解析微分 | pysad 数値微分 (PYSAD_NUMDERIV=1、同じ重み) |
|---|---|---|---|
| QFA1H K1 | .40014410823165847 | 0.4001441082316458 | 0.4001441095561435 |
| QDA2 K1 | −.5641039881930263 | −0.5641039881929922 | −0.5641039898980917 |
| QDA3 K1 | −.43699540131878234 | −0.4369954013188038 | −0.4369954013144212 |
| QFA4 K1 | .38110158014084466 | 0.3811015801408816 | 0.3811015806377782 |
| QDA5 K1 | −.4011725279548499 | −0.4011725279548886 | −0.4011725285008895 |
| QFA6H K1 | .8857122301470066 | 0.8857122301470607 | 0.8857122302530256 |
| QFA8H K1 | 1.141437202315654 | 1.141437202315554 | 1.141437199990928 |
| BPI6 K1 | −.949639998526657 | −0.9496399985266456 | −0.9496399978870261 |
| BPI7 K1 | −.6613806439882316 | −0.6613806439881901 | −0.6613806432900691 |

解析微分版は 12〜13 桁一致、pysad の数値微分版も同じ反復列 (Matched 8.7339E-09) で 9 桁一致。
つまり esrf/matching で pysad が別解に着いていた原因は微分ではなく **EMITX=NaN 時の重み** だった
(重みが 1e5 倍違うと最小二乗の解が違う)。EMITX を MARK の値 1.8e-9 のまま (pysad の現状) にすると、
解析微分でも 3 0.2449 (0.3898) … 12 3.4647E-08 → Matched 8.0925E-09 の別解 (QFA1H 0.38619…) になる。

## 7. match6.py / ffscmd.py への統合案 (既存ファイルは今回未変更)

1. **match6.Matcher.go の微分ブロックをメソッドに切り出す**: `qu = self._derivative(ids, dfs, rows, wvar)` とし、
   既定実装を `match6d.analytic_qu` + FitFunction 行の片側差分 + 未対応列の中心差分 (AnalyticMatcher._derivative の中身)
   にする。`numerical_derivative` 属性 (FFS$NumericalDerivative) で従来の中心差分に切り替え。
   こうすれば AnalyticMatcher の go() の写しは不要になり、ffscmd は変更なしで解析微分になる。
   (`self.last_qu` に (ids, qu, wvar) を残しておくとデバッグに便利。)
2. **ffscmd: EMIT 失敗時に EMITX/EMITY/EMITZ=NaN** にする (SAD temit.f 503 行: 閉軌道が見つからなくても
   rsetgl1('EMITX', eemx) まで進み NaN が入る)。pysad は現状 MARK の値を保持しているので SAD と重みが違う。
3. **match6._tffscalc の重み計算を gfortran の MAX/MIN 意味論で**: `max(a, b)` は「b > a なら b、さもなくば a
   (NaN は無視)」。em = max(emminv=1e-15, |emx|+|emy|)、coum = min(1, max(MINCOUP=0.01, 1e-3, 1/(…)))、
   coum/(1+coum)、emxx = max(emx, coum·em)。match6 は現状 `em = max(1e-30, …)` と書いており SAD の emminv=1e-15、
   coumin=MINCOUP(0.01, tfsetparam で rgetgl1) と違う (emx が正常なら影響しないが NaN/0 のとき違う)。
   SAD 本来の式: tffscalc.f 442-452。
4. `FFS$NumericalDerivative` を ffscmd のグローバルとして受け、Matcher に渡す。
5. 残り: qdtrans (MULT/SOL/L 等の数値微分列: いまは match6 と同じ中心差分で代用。SAD は要素 1 個の写像だけを
   数値微分するので結果はほぼ同じ)、qdthin (L=0 QUAD)、qdini (先頭 MARK の変数)、qddtwiss (DZ/DDP 条件と INS)、
   BMAG、CouplingMatrix、COUPLE (coup≠1)、境界分数 wk、gammab (加速) — いずれも esrf スクリプトには出ない。

## 8. ファイル

- `pysad/match6d.py` — 本体 (日本語コメント)。
- `bench_match6d_check.py` — 13B と ARCA の照合 + 13B の AnalyticMatcher 実行。
- `run_matching_analytic.py` — esrf/matching を AnalyticMatcher で実行 (環境変数 PYSAD_EMITX, PYSAD_NUMDERIV)。
  出力 <work-esrf>/matching1.pysad_analytic.out, matching1.pysad_numeric.out。
- `<work-esrf>/matching1_dbg.sad`, `matching1_dbgnum.sad`, `*.sad.out`, `run_dbg.sh` — SAD の確認実行
  (`wsl -d Ubuntu -e bash "<repo>/<work-esrf>/run_dbg.sh"`)。実行で <work-esrf>/datas.txt に追記される (スクリプトの OpenAppend)。
