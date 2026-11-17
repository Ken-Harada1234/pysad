# 40. CSR (コヒーレント放射) の解説

最終更新: 2026-11-17

偏向磁石の中でバンチが出す CSR (coherent synchrotron radiation) のウェイクを、SAD K64 で計算・追跡する道具の解説。
対象は 2 層ある。K64 本体の組み込み関数 (Fortran: `csrtrack.f`、`csroy.f`、`csrtest.f`) と、`Packages/CSR.n` の 2 つのクラス `CSR` (ダクトの断面からインピーダンスを求める) と `CSRStability` (ウェイクを使う追跡と安定性)。
pysad は組み込み関数を `pysad/csr.py` に直訳し、`CSR.n` は書き直さず SAD のソースのまま pysad のインタプリタで動かす (`pysad/pkgload.py`)。
関数ごとの照合の一覧は `docs/36_spin_csr.md` の CSR 節。ここは初めて読む人向けに、全体の流れ・使い方・物理・癖・限界をまとめる。

## 1. 何ができるか

| 層 | 名前 | 内容 |
|----|----|----|
| 組み込み | `CSRInit` | ダクト断面の格子上で、バンチの電荷分布が作る静電場 (2 成分) を解く |
| 組み込み | `CSRMatrix` | ダクト断面の格子上の差分行列 (波動方程式の固有値問題の行列)。mode で種類を選ぶ |
| 組み込み | `CSRConvert` | 断面の格子 (列ごとの表) と 1 次元ベクトルの相互変換 (K64 内部名 tfcsrvec) |
| 組み込み | `CSRTrack` | ウェイクの kick を入れたマクロ粒子の縦方向追跡 (z, dp) |
| 組み込み | `CSRHaissin` | Haissinski 方程式 (ウェイクのある定常分布) を Newton 法で解く |
| 組み込み | `CSRSetupOY`、`CSROYMatrix` | Oide-Yokoya の線形安定性の行列を作る (作用角の格子) |
| `CSR.n` の `CSR` | `MakeZL`、`Ez`、`ZL`、`WL`、`ZThre`、`SkinDepth` ほか | ダクト断面と曲率半径から CSR インピーダンス ZL(k) を作る |
| `CSR.n` の `CSRStability` | `InitializeParticles`、`Track`、`ScanTrack`、`OYStability` ほか | ZL を使ったバンチ長・エネルギー幅の時間発展と閾値の探索 |

全体の流れ: ダクト断面 → `CSRMatrix`・`CSRInit`・固有値分解で ZL(k) → (ファイルに保存) → `CSRStability` が ZL を格子 (FFT 形式) に直し、`CSRHaissin` で定常分布、`CSRTrack` で追跡、`CSRSetupOY`・`CSROYMatrix` で安定性。

## 2. 使い方 (動く SAD スクリプト)

以下は `pysad/tests/test_csr.py` などから取った。2.1〜2.3 は gs (K64) と照合済み。

### 2.1 CSRTrack: ウェイク付きの追跡

```
FFS;
$FORM="18.14";
SeedRandom[12];
nr=256; dz=2e-5; np=2000;
zl=Table[0.02*Cos[0.3 i]/(1+i/10.), {i,1,nr}];      (* 実 FFT 形式の ZL。ここは試験用の人工値 *)
z=3e-4*GaussRandom[np]; dp=1e-3*GaussRandom[np];
r=CSRTrack[z,dp,zl,{dz,4,0.5,0.1,-0.5,0,0,1e-3,1e-3},5];
Print["R ",r];
Exit;
```

`CSRTrack[z, dp, zl, {dz, ndiv, inten, r56, r65, frf, phirf, damp, sigmae}, nturn]`:
- `z`, `dp` は粒子ごとの縦位置とエネルギー偏差の表。**その場で書き換わる** (K64 はリスト本体を更新する)。
- `zl` は 2 のべき乗の長さ `nr` の実 FFT 形式のインピーダンス (k = 0 から Nyquist まで、実部・虚部の組。Nyquist は 2 番目の要素に入る)。
- `dz` はビン幅 (m)、`ndiv` は 1 周を何分割してウェイクを当てるか、`inten` は粒子数に比例する強度 (内部で `inten/粒子数` を 1 粒子の重みにする)。
- `r56`、`r65` は 1 周の縦方向の行列要素 (m、1/m)。`frf`、`phirf` は RF の周波数 (Hz) と位相 (rad)、0 なら RF なしで `r65*z` の線形復元力。`damp` は 1 周の放射減衰、`sigmae` は平衡エネルギー幅。
- 戻りは各周回の {σz, σδ} を順に並べた長さ `2·nturn` のリスト。σz が `nr·dz/4` を超えるか σδ が上限を超えると発散とみなし、その周以降は 0 を返す (乱数は ffs と共通の `SadRandom`)。

### 2.2 CSRHaissin: 定常分布

```
rho=CSRHaissin[zl,{dz,1e-8,-0.3,3e-4}];     (* {dz, inten, r65, sigz} *)
Print[Length[rho]," ",Max[rho]];
```

ウェイクのある定常な縦分布 (長さ `nr` の表、前半だけが非零) を返す。解けなければ `General::toomany` で停止する。

### 2.3 Oide-Yokoya の行列

```
s=CSRSetupOY[rho,zl,{dz,6,3,1e-8,0.05,-0.3,3}];       (* {dz, nj, nphi, inten, r56, r65, zspan} *)
m=CSROYMatrix[s,{dz,1e-8,0.01,3e-4,1e-3}];             (* {dz, inten, nus, sigz, sige} *)
Print[Dimensions[m]];                                  (* {18,18} = nj*nphi の正方行列 *)
```

`CSRSetupOY` は作用角 J の `nj` 点の格子、方位角のモード `nphi` で線形化する。戻りは {wi, J, ω, zj, dpj, rhoj}。**`CSROYMatrix` は K64 (gs) では動かない** (4 節の B60)。pysad は動く (行列を返すが gs との照合はできない)。

### 2.4 ダクト断面の行列と電場

```
a=CSRMatrix[{3,4,4,3},{1e-3,2e-3,1.1,1.3,0.9,0.7,0.1},3e3,5.0,2,22];
v=CSRConvert[Range[14],{3,4,4,3}];                      (* ベクトル → 列ごとの表 *)
r=CSRInit[{{0.5e-3,0.5e-3},4,1e-3,{4,6,8,8,8,8,6,4}}];   (* {{dx,dy}, 中心列, rho, pipe} *)
```

- `pipe` は x 方向の各列の y 方向の格子点数の表 (中心線から上半分の断面)。`CSRConvert` は点を列ごとに上から数えて並べる。
- `CSRMatrix[pipe, {dx, dy, f0, f1, f2, f3, tilt}, omega, rho, center, mode]`: `omega` は ω/c = k (1/m)、`rho` は曲率半径 (0 は直線部)、`center` は軌道の列。`f0..f3` は境界の補正係数 (fudge)。mode は 22・33・23・32 (縦・横の 4 つのブロック)、110・220・330 (抵抗壁の減衰用。第 1 列だけのベクトル)、1000 (らせん断面)。
- mode 1000 (`csrhelmat`) は境界の点列を中点で 2 倍に細かくして境界要素法で各セルの電場への応答を作る (2nm × 2nm の行列)。K64 の添字のずれまで写してある。
- `CSRInit` の第 1 要素に実数 1 つを渡すと、K64 は dx を読まず 0 のまま (バグ) で常に `wrongtype`。必ず `{dx, dy}` で渡す。

### 2.5 CSR.n: ダクトから ZL を作る

`CSR`、`CSRStability` は `init.n` の AutoLoad の表にあり、pysad も最初に使われたときだけ `Packages/CSR.n` を読む (`PYSAD_SAD_PACKAGES` か、pysad の隣の `SAD-master/Packages` が必要)。

```
FFS;
$FORM="18.14";
PageWidth=1999;
c=CSR[NParallel->1,Nk->4,BunchLength->3e-4,BendingRadius->10,PipeWidth->0.04,Center->0.02,MeshRatio->{2,2}];
Print["A ",c@SkinDepth[1e3]," ",c@EffectiveSkinDepth[2e3]," ",c@EffectiveCond[2e3]];
Print["C ",c@StandardKList[]];
Print["D ",c@VariableMesh[3000.,10.]];
Print["E ",c@ZThre[1e-3,0.01,1e-3]];
z=c@MakeZL[{{10.,1.}},2000.,6000.,3];                  (* 約 7 分 (pysad) *)
Print["Z ",z];
Exit;
```

(`MakeZL` は `test_csr_makezl_matches_gs` で gs と照合した。)

主な使い方:
- `CSR[オプション...]`: ダクトを作る。主なオプションは `BunchLength` (σz)、`BendingRadius`、`PipeWidth`、`PipeHeight`、`PipeFun` (x の関数としての半高さ。既定は定数 0.05)、`Center` (軌道の x 位置、既定は幅の半分)、`MeshRatio`、`Nk` (波数の点数)、`NParallel` (並列数。`Parallelize` を使う)、`Resistive` (抵抗壁)、`Material`・`Cond`、`CacheFile` (途中結果の保存)、`Adaptive` (波数の適応刻み)。
- `MakeZL[l, k1, k2, n]`: `l = {{ρ1, s1}, {ρ2, s2}, ...}` は 偏向の曲率半径と長さの並び (ρ=0 は直線部) で、波数 k1〜k2 の `n` 点の ZL を作る。戻りは `{k, ZLr, ZLi}` の表 (単位は ZL の内部単位。Ω に直すには `ZUnit` を掛ける)。`MakeZL[l]` (引数なし) は `StandardKList[]` (k = (i/Nk)^(4/3) · 3.5/σz) を使う。
- `ZL[omega, l]`: 1 つの波数の ZL。`Ez[omega, rho, s]`: 定常な曲率半径 ρ の偏向の長さ s での縦電場 (`{k, Re, Im}`)。引数を 1 つ減らした `Ez[z, sz]` は `MakeEz` が作った表から z での縦電場 (Gauss バンチ) を数値積分で返す。
- `WL[z, sz]`: ZL から時間領域のウェイク関数 (V/pC)。`WLPlot`・`ZLPlot`・`PipePlot`・`EzPlot`・`Er0Plot`・`WLPrint`・`ZLPrint` は描画・表示。
- `ZThre[acoulomb, r56, sigdp]`: 簡易な不安定閾値の式 (`sqrt(π/2)·(p/m)·Z0/r_e/(N)·|R56|·σδ²·c`)。
- `WriteZL[file, ...]`: ZL を `{"BunchLength"->σz, "ZL"->{ZLr, ZLi}}` の形でファイルに保存。`ReadWake[file]`: 外部のウェイク表 (z, W) から ZL を作る。`MakeResZL[a]`: 抵抗壁の解析式の ZL。

### 2.6 CSRStability: 追跡と安定性 (書き方の骨格。通しの照合は未)

```
(* 未実行の骨格: ZL ファイルは MakeZL → WriteZL で作ったもの *)
s=CSRStability[BunchLength->3e-4, MomentumSpread->1e-3, DampingRate->1e-4,
   SynchrotronTune->-0.01, MeshSize->3e-5, Meshes->512, Division->4];
s@SetupZL[{{"zl.dat",1}}];                       (* ZL の読み込み (複数ファイルの和) *)
{z,dp}=s@InitializeParticles[20000,aw];          (* Haissinski の分布から z、dp は Gauss *)
s@Track[aw,z,dp,2000,100];                       (* 2000 周、100 周ごとに描画 *)
s@ScanTrack[awlist,np,{t1,t2}];                  (* 強度のスキャン *)
r=s@OYStability[aw,nj,nphi];                     (* {安定か, 成長率, 固有値, 固有ベクトル} *)
```

- `Constructor` は `r56 = 2 sin(μ)/(1+cos μ) · σz/σe`、`r65 = −4 sin²(μ/2)/r56` (μ = 2π·νs。νs は α>0 なら負で与える) とし、ZL の格子の幅 `zspan = dz/σz · nr/2` を決める。
- `GetZL` は ZL を k = 2π/(dz·nr) の刻みの格子に直し、ビン幅の平均化 `exp(−(k dz)²/12)` を掛け、バンチ長補正 `min(exp((kσ)²/2), enh)` を掛けた実 FFT 形式の表にする。`SetupZL` はそれを `SetZl` で `Zl` に入れる。
- `Haissinski->True` (既定) だと `InitializeParticles` は `CSRHaissin` の解から z を作る。解けなければ Gauss にフォールバックする。
- `Track[aw, z, dp, nt, ntplot]` は `CSRTrack` を `ntplot` 周ずつ呼び σz/σz0、σδ/σδ0 を描く。`ScanTrack` は強度ごとに `CSRTrack` を走らせ、後半の平均と揺らぎを描く。`ScanBL` は バンチ長と νs のスキャン。`TrackAnimate` は `ListPlot` のあと `System["import -silent -window … -format -jpeg …"]` (ImageMagick) と `System["mkdir …"]` を呼ぶので、`System` を実行しない Windows の pysad では画像を作れない (K64 の `CSR.n` のソースで確認、2026-11-16)。
- `OYStability[aw, nj, nphi, Range->4, BLFactor->1]`: `CSRHaissin` → `CSRSetupOY` → `CSROYMatrix` → `Eigensystem`。成長率は `max(|Im √ev| − damp·clip(|Re √ev|/(2π νs), 1, ∞))`、負なら安定。`GuessOYStability`・`GuessOYStability1`・`OYScan` は BLFactor を動かして閾値を探す (割線法/FindRoot)。

## 3. 物理と K64 の式

### 3.1 ダクト断面のモデル (CSR.n と csrtest.f)

ダクトは矩形ではなく「x の各列の高さ」`pipe[i]` (`PipeFun[x]/dy` を丸めた格子点数) で表す断面で、軌道は `Center` の列にある。`dx = PipeWidth/mx`、`dy = min(dx·MeshRatio1/MeshRatio2, MaxMeshSize)`。メッシュ数 `mx` は `VariableMesh[ω, ρ]` が `(|ρ|/ω²)^(1/3)` に比例して決め、`mxmin=20`〜`mxmax=3201` に制限する (ω が大きいほど細かい)。

波数 ω (=k) ごとに、縦 (22, 33) と横 (23, 32) の場の成分を 1 つの行列 `[[a22, a23],[a32, a33]]` にまとめ (`CSRDataCache`)、`RealEig` で固有値 ev・固有ベクトル u を求める。`ρ = 0` (直線) の行列は k=`omega0` (10000) で 1 度だけ作り、`ev·(omega0/ω)` で他の k に使い回す (`CSRStraight`)。
`CSRMatrix` の係数 (ac, ce, cg) は `ωρ` を含む (`0.5ω/((ωρ)²−1)` など)。曲がったダクトでは `1/(ωρ)² → 1` のところに特異点がある。
境界付近の格子には補正係数 `fudge0〜3` が掛かる。`CSR.n` の既定 (`fudge1=-1.83985`、`fudge3=1.30899`) は円形と角形ダクトの固有値に合わせた当てはめ値であり、物理の導出値ではない。

### 3.2 固有値・共鳴減衰・インピーダンス

- 固有値が複素になる (`ComplexEig`/`RealEig`) と、虚部を捨てて逆行列の許容誤差を `svdcmp` (2e-2) に緩める。固有値が近すぎる (`epsev=1e-4`) ときは `sttaper` のテーパーを加えて最大 10 回やり直す。
- 抵抗壁 (`Resistive->True`) は `ResDamp` が、壁の格子点を拾う mode 220・330 のベクトルの固有モードへの射影から減衰率 dev を作り、`EffectiveSkinDepth` (TiN 被膜の 2 層の表皮厚) を掛ける。ZL では `exp(−dev·s)` の減衰として現れる (`ZL` の `a=dev*s`)。これが「共鳴の減衰」。
- 初期の静電場は `CSRInit` が解き、`CSRConvert` でベクトルにして固有ベクトルに射影する。
- 区間 (ρ, s) の ZL は、固有モードごとに `cos(ev·s)`、`sin(ev·s)` の積分 (`coswi=sin w/ev`、`sinwi=2 sin²(w/2)/ev`、`w=ev·s`) で場を進め、偏向を渡るたびに初期条件を引き継ぐ (`ZL` の `Do` ループ)。最終の縦インピーダンスは軌道の列 (`cent0`、`cent1`) の場の差分を `Wnorm=1e-12/(2π ε0)` で規格化したもの。
- 全体の波数の表を FFT の格子に載せるときは、負の k に対称に延ばし (実部は偶、虚部は奇)、自然スプラインで補間する。

### 3.3 ウェイクの追跡 (csrtrack.f)

1 周を `ndiv` 回に分け、各回で (1) `z += (r56/ndiv)·dp`、(2) 粒子を `nr` のビンに線形配分して密度 rho を作る、(3) rho の実 FFT に ZL を掛けて戻しウェイクを得る (`csrwake`)、(4) 3 次スプラインで各粒子の位置のウェイクを補間し `dp += (inten/np/ndiv)·W`、(5) 分割の真ん中で `r65·z` (RF ありなら `vr·sin(ωz+φ)−dpr`) の復元力、(6) `ad=1−2·damp/ndiv` の減衰と `bd` の量子励起の乱数、を行う。周ごとに重心を引いて σz、σδ を記録する。RF ありは粒子の分離 (|z| が RF バケットの半周期を超える) や `dp<−DPMAX` の損失を扱う。

### 3.4 Haissinski と Oide-Yokoya

- Haissinski: ρ(z) = exp(u(z)) と置いた非線形方程式 `u = u0 − aw·(W*ρ)` (u0 は線形復元力の Gauss) を、規格化条件を付けた Newton 法 (しきい 1e-5、最大 100 回。ヤコビアンは `tsvd` の最小ノルム解) で解く。`aw = inten/(r65·σz²)`。
- Oide-Yokoya: 定常分布のまわりの線形摂動を、作用 J・方位角 φ の格子で展開した 1 周写像に直す。`CSRSetupOY` は各 J の軌道を `dt = min(1, 0.01/√|r56·r65|)` の刻みで周回して ω(J)、位相空間の座標を作り、`CSROYMatrix` が (nj·nphi)² の行列を返す。その固有値 ev の平方根の虚部が成長率で、`OYStability` は減衰 `damp` を差し引いて安定を判定する。

## 4. K64 の癖

(`docs/21_sad_bug_report.md` の B59・B60。いずれも pysad は K64 どおりに再現する。)

- **B59**: `CSRSetupOY` は `dpjn` を `zjn` と同じ記憶域に割り当てる。返す 4 番目の表 (zj) は z ではなく dp の値、5 番目 (dpj) は 0 のまま。続く `CSROYMatrix` は zj を位置として使うので、行列は物理的には狂う。pysad は同じ上書きを写した。
- **B60**: `CSROYMatrix` が、`CSRSetupOY` の出力でも手で作った実数リストでも `General::wrongtype` を返す (この K64 ビルドで再現)。登録名は `CSRSetupOY`・`CSROYMatrix` だが、ソース内の関数名は `tfcsroysetup`・`tfcsroymat` で綴りが違う。原因は未調査。pysad はソースどおりに翻訳して動くが、gs と照合できない。したがって `OYStability` と `OYScan` は gs では完走しない可能性が高い (gs で確かめていない)。
- `CSRInit` の実数 1 つの dx は常に `wrongtype` (dx を読まない)。
- `CSRMatrix` mode 1000 は、K64 が `phix(0:nphimax)` を `csrphimesh` の `phix(*)` に渡すため、点の添字が 1 つずれる (pysad も同様)。また `{4,4,4,4}` などの断面は gs 自体が SIGABRT/SIGSEGV で落ちる (参照なし)。
- `CSR.n` の `FreeMem[]` は `top` の出力 (`OSTYPE` で Switch) を読むので、非対話の WSL では `OSTYPE` が空で `OpenRead` が失敗する。gs の参照を取るときは WSL に `top` を入れるか、`Mem: ... free ...` を出す代用を PATH に置き、`OSTYPE=linux` を設定した。pysad は Windows などの `top` のない環境では代用の行を内部で作る。
- `CSR.n` は GUI・描画 (`ListPlot`、`Plot`、`HistoPlot`、`WidgetInfo`) を呼ぶ。`pkgload.py` の許可リストには `CSR.n` と `ListPlot.n`・`Plots.n` が入っている (描画の可否は docs/37)。

## 5. gs との一致の度合い

すべて SAD K64 (gs、2026-10-19 ビルド) の出力との比較 (試験は `pysad/tests/test_csr.py`)。

| 項目 | 一致 |
|----|----|
| `CSRTrack` (RF あり・なし、発散後 0、z・dp のその場の更新) | 相対 1e-11 |
| `CSRHaissin` | 相対 1e-11 |
| `CSRSetupOY` の J・ω・rhoj | 1e-9〜1e-7 (周回の離散計算のため) |
| `CSRMatrix` mode 22・33・23・32・110・220・330 | 全要素 |
| `CSRConvert` | 一致 |
| `CSRInit` (pipe 8 列と 5 列) | 相対 1e-10 |
| `CSRMatrix` mode 1000 ({2,2,2} と {5,5,5,5,5,5}) | 相対 1e-12〜1e-13 |
| `CSR.n`: `SkinDepth`、`EffectiveSkinDepth`、`EffectiveCond`、`PathDifference`、`StandardKList`、`VariableMesh`、`ZThre`、`Tweak` | gs と一致 |
| `CSRData` (8 列のメッシュ) の固有値 | 相対 1e-15 (固有ベクトルは縮退で並びが違うが物理量は同じ) |
| `Ez` | 相対 5e-10 |
| `MakeZL` (`{{10.,1.}}`、2000〜6000 /m の 3 点、PipeWidth 0.04、BunchLength 3e-4) | 相対 1e-10 以内 |

この文書の作成時に `py -3.12 -m pytest pysad/tests/test_csr.py` を `PYSAD_SAD_PACKAGES` 付きで実行して 6 通過・1 skip (skip は `PYSAD_SLOW=1` が要る MakeZL) を確かめた。

試験で使った量は小さい値の人工的な ZL・断面で、実機のパラメータ (エネルギー、実際のダクト、周長) での物理的な妥当性の検証ではない。ZL の `MakeZL` の照合は 3 周波数・1 区間だけである。

## 6. 遅い箇所

- `MakeZL` (上の 3 点) は pysad で約 7 分 (試験は `PYSAD_SLOW=1`)。メッシュ数が波数とともに増え (`mxmax=3201`)、各 k で `CSRMatrix` と固有値分解 (密行列) が走る。k の点数 `Nk` の既定は 32 (+ `StandardKList`)、実用の ZL の作成は桁違いに長くなる (未計測)。
- `CSR[...]` の `Wait[]` は `FreeMem[]` (K64 では `top` を呼ぶ) と `Sleep` で空きメモリと並列数を見て待つ (`NParallel` 本)。gs はこの環境に `top` が無く `FreeMem` が OpenRead エラーで失敗する。pysad の `FreeMem` は Python (psutil) で代替して 1.0 を返す。pysad は fork しないので、並列を多くして速くなるかは確かめていない (未確認)。
- `CSRTrack` と `CSRHaissin` の中の `tsvd` (ヤコビアン) は numpy で速い。`CSRSetupOY` は nj・nphi に比例して重い。

## 7. 限界・未確認

- `ScanTrack`・`OYStability`・`OYScan`・`GuessOYStability*`・`ScanBL`・`TrackAnimate` (画像の生成は pysad では不可、上記)・`SetupZL`・`GetZL`・`CSRStability` のコンストラクタ・`ReadWake`・`MakeEz` の長い流れは、gs と照合していない (pysad での完走の試験もない)。`OYStability` は B60 により gs では参照が取れない。
- 描画 (`WLPlot` `ZLPlot` `PipePlot` `EzPlot` `Er0Plot`) の出力は確認していない。
- `CSROYMatrix` の値は gs と未照合 (ソースどおりの翻訳)。
- 抵抗壁 (`Resistive->True`、`ResDamp`) と `Adaptive->True` の経路は gs と照合していない。
- 中心がダクトの中央でない (`Center≠PipeWidth/2`) 場合、`PipeFun` が定数でない断面、`MeshRatio` が 1 でない場合の一致は、`MakeZL` の 1 例 (`Center=0.02`、`PipeWidth=0.04`、`MeshRatio->{2,2}`、一定の高さ) 以外は未確認。
- `CSRMatrix` mode 1000 は小さな断面 2 例のみ。{4,4,4,4} など gs が落ちる断面は参照がない。
- 物理モデルとしての CSR の公式 (ウェイク・Haissinski・OY の式) の独立な導出や、実験との比較はしていない。K64 の実装の直訳であり、`fudge` 係数は当てはめ値。

## 8. 関連文書

- `docs/36_spin_csr.md`: 関数ごとの照合の記録 (CSR 節)。
- `docs/21_sad_bug_report.md`: B59、B60 の K64 の癖。
- `docs/23_wake.md`: 進行波・ウェイク一般 (WAKE)。
- `docs/37_plot_canvas.md`: 描画 (`ListPlot`、`Plot`、Canvas)。
- `docs/22_gui_tclbridge.md`: GUI。
- `docs/34_sad_packages.md`: `Packages/*.n` の自動読み込みの仕組み。
- `docs/history/work_unimpl.md`: 未解決の差と最新の記録。
- コード: `pysad/csr.py`、`pysad/pkgload.py`。K64 のソース: `Packages/CSR.n`、`src/csrtrack.f`、`src/csroy.f`、`src/csrtest.f`。
- 試験: `pysad/tests/test_csr.py`。
