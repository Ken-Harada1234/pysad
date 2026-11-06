# 09. DynamicApertureSurvey — dapert.n → tfdapert → trackd と pysad の dapert.py

最終更新: 2026-11-06 (現状のスナップショット)。

> 出典 (すべてソースで確認): `Packages/dapert.n`, `src/tfdapert.f`, `src/trackd.f` (trackd, tinip1, tpdamp, rad62a),
> `src/tday.f` (rad62), `src/tapert.f`, `src/tturn.f`, `src/tconv.f`, `src/initb1.f`;
> pysad 側は `pysad/dapert.py`, `pysad/ffscmd.py` (fn_dasurvey), `pysad/track6n.py` (turn, k_apert)。

## 1. 呼び出し経路

```
DynamicApertureSurvey[range, nt, opt...]        Packages/dapert.n  (DALib クラスのメソッド)
  └─ DynamicAperture[range, lfno, ref, {phix,phiy,phiz}, damp, dampenough, {iv1,iv2}]
       └─ tfdapert (tfdapert.f)      引数検査、NP の一時変更、codin の保存
            └─ tfdapert1              tfsetparam, tclrpara, dapert=.true., 終了後 tftclupdate(7)・FPE 検査
                 └─ trackd (trackd.f) 粒子の配置・追跡・損失周の記録・地図の表示・ResultOfDAPERT の構築
                      ├─ tinip1       規格化振幅 → 物理座標
                      ├─ tturn        TrackParticles と同じ 1 周ルーチン (損失判定 tapert を含む)
                      └─ tpdamp       Damping / DampEnough のとき各周の後に呼ぶ
```

pysad では `ffscmd.fn_dasurvey` が dapert.n + tfdapert の役を、`dapert.dynamic_aperture` が trackd の役を担い、
1 周は `track6n.turn(..., apert=...)` で回す。

## 2. SADScript 側 (Packages/dapert.n)

### 2.1 既定オプション

| オプション | 既定 | 意味 |
|---|---|---|
| Output | 0 | 出力ユニット。0 は実質 /dev/null、6 で標準出力、−1 は tfdapert で icslfno() に置換 |
| ExpandElementValues | True | True なら先に `Element["EXPAND"]` |
| ReferenceOrbit | Default | 長さ 6 の実数リストを与えると tfdapert が `codin` を一時的に置き換える |
| PhaseX, PhaseY | 0 | 走査軸の初期位相 (§4.3) |
| PhaseZ | −π/2 | 同上。既定では z=0, δ=a になる |
| Axes | "ZX" | Mode 表 (下) |
| Damping | False | 0 以外なら tpdamp で各周に放射減衰を掛ける |
| DampEnough | 0 | 0 以外なら「十分減衰した」粒子を早期に生存扱い (tpdamp) |

Mode 表: `"ZX"→{3,1}`, `"ZY"→{3,2}`, `"XY"→{1,2}`, `"XZ"→{1,3}`, `"YX"→{2,1}`, `"YZ"→{2,3}` (小文字は大文字化、他は `DALib::wrongmode`)。
1 文字目が **ivar1 = 値リストで与える軸**、2 文字目が **ivar2 = 51 点で走査する軸**、残りが ivar3。

### 2.2 本体の処理

1. `TURNS0=TURNS; TURNS=nt; NP0=NP` と退避してから一時設定。
2. `n1 = Min[201, Length[range[[ivar1]]]]`、`NP = Min[4000, n1*Min[DAPWIDTH,51]*2]` (DAPWIDTH の既定は initb1.f で 7)。
3. `Write[lfno, " Turns =", nt, " Maximum number of particles =", NP]` を Output ユニットに出す。
4. `score = DynamicAperture[...]` を呼び、NP, TURNS を元に戻す。
5. 戻り値は `{score, ReplacePart[ResultOfDAPERT, range, 1]}` — ResultOfDAPERT の第 1 要素 (min/max の組) を
   ユーザーの range そのものに置き換えたものを返す (グローバルの `ResultOfDAPERT` 自体は trackd が作ったまま)。

## 3. tfdapert (tfdapert.f)

| 項目 | 内容 |
|---|---|
| 引数 | 7 個固定。第 7 引数 `{v1,v2}` から `iv(1), iv(2)` を取り、`iv(3) = 6 − iv(1) − iv(2)` |
| range | 長さ 3 のリスト。`range[[iv(2)]]`, `range[[iv(3)]]` は `{min,max}` (→ `range(1:2,1)`, `range(1:2,2)`)、`range[[iv(1)]]` は実数リスト r1 で `n1 = min(200, 長さ)` |
| Output | 実数。−1 なら `icslfno()` |
| ReferenceOrbit | 長さ 6 の実数リストなら `codsave=codin` の上で `codin` を差し替え、終了後に復元 |
| NP | `np0 = max(2000, n1, n1*min(51, nint(DAPWIDTH)))` を設定し `NP` にも書き込む。終了後に復元 |

tfdapert1 は `tfsetparam`, `tclrpara` の後 `dapert=.true.` にして trackd を呼ぶ。このフラグで tapert の損失報告
(`tapert_report_dropped`) が抑止される。終了後 `tftclupdate(7)`、FPE があれば `DynamicApertureSurvey-FPE` を表示。
dapert.n の 201 と tfdapert の 200、trackd の `n1p0=256` と上限が三重になっている点に注意 (実効は 200)。

## 4. trackd のアルゴリズム (trackd.f)

### 4.1 定数と格子

| 名前 | 値 | 意味 |
|---|---|---|
| n2p | 51 | ivar2/ivar3 を走査する点数 (地図の 1 行 = 51 文字) |
| n1p0 | 256 | ivar1 の値の上限 |
| maxturn | 2^29 | ntloss の「未起動」印 |
| nw | 16 | 全スロットが埋まっていないとき、この周数ごとに配置ループへ戻る |
| minnp | 16 | 並列数の上限 `np0/minnp` に使う |

- `a2step = (a2max−a2min)/50`, `a3step = (a3max−a3min)/50`。j 番目 (1..51) の振幅は `a2min + a2step*(j−1)`。
  ivar2 と ivar3 は **同じ j** で同時に動く (ESRF の `{1,400},{1,1}` では y は常に 1)。
- `a1i(1:n1p) = r1` (n1p = min(256, n1))、a1min/a1max はその min/max。`n1==0` の経路 (`tftrad` 用、r1(3) を刻みにする) は DynamicApertureSurvey では使われない。
- `emx = sqrt(|EMITX|+|EMITY|)`、`emz = sqrt(|EMITZ|)` (RFSW) または `|SIGE|` (RFSW でないとき)。
- `muls` (1 文字あたりの周数): muld=1, mule=1 から `muld*mule*62 < nturn` の間 1→2→5→(1, mule×10) と進める。
  nturn=20→1, 100→2, 1000→20, 5000→100。

### 4.2 粒子スロットの管理

- `npara = min(NPARA, 256, np0/16)`。npara>1 なら `itffork` で子プロセスを作り、子は `tfaddseed` で乱数の種をずらす。
  `npmax = max(1, min(n1p*ncons/npara, np0))` (ncons = max(1, nint(DAPWIDTH)))。ntloss は共有メモリ (`ktfallocshared`)。
- 配列: `kptbl(k,1)` = 粒子 ID k → 配列位置、`kptbl(ip,2)` = 配列位置 → ID。`kzx(1:2,k)` = 粒子 k の (i,j)、0 なら空き。
  `mturn(k)` = 粒子 k が完了した周数。`nxm(i)` = 行 i で「損失が見つかった最小の j」(初期 52)。
- **LOOP_K (配置ループ)**: 空きスロット k ごとに、`i=1..n1p`、`j=nxm(i)−1 → 1` の順で未起動 (`ntloss==maxturn`) の点を探す。
  `ntloss(i,j) < nturn` (損失済み) を見つけると `nxm(i)=j` に窓を下げ、`j < nxm(i)−ncons` で行を打ち切る。
  つまり **DAPWIDTH は「損失点から下へ何点まで同時に追跡するか」の窓幅**で、DAPWIDTH=51 なら行全体を一度に起動する。
  起動時に `ntloss(i,j)=nturn` (暫定「生存」)、`kzx=(i,j)`、`mturn=0` を入れ、座標を設定して tinip1 を呼ぶ。
  起動する点がなくなると空きスロットを末尾の粒子と入れ替えて `np1` を減らす (詰め替え)。
- **追跡ループ**: `np>0` の間 `tturn(np,1,nlat,...)` を呼び、(Damping/DampEnough なら tpdamp)、`n` を進め、
  起動中の粒子ごとに「配列位置 ≤ np なら生存 → mturn++ (nturn に達したら ntloss=nturn で解放)」
  「np より後ろなら損失 → `ntloss = mturn` (= 損失前に完了した周数、1 周目損失なら 0) で解放」。
  解放があれば配置ループへ戻る。なければ、スロットに空きがあるときだけ nw=16 周ごとに配置ループへ戻る。
- 子プロセスは追跡後 `exit_without_hooks`、親は `waitpid_nohang` で待つ (600 回 (約 600 s) 応答がなければ `tkill`)。

起動順が i=1..n1p、j=51→1 (大振幅から) であることは、RAD+FLUC の乱数列を pysad で再現するために重要 (§6)。

### 4.3 初期座標 (tinip1 と位相)

trackd が設定する規格化座標 (振幅 a は `Sqrt[2J/(EMITX+EMITY)]` 単位):

| 軸 | ivar1 のとき | ivar2 / ivar3 のとき |
|---|---|---|
| X | x=a, px=0 | x = a cosφx, px = −a sinφx |
| Y | y=a, py=0 | y = a cosφy, py = −a sinφy |
| Z | z=0, δ=a | z = a cosφz, δ = −a sinφz (PhaseZ=−π/2 → z=0, δ=a) |

tinip1: `xa = (x,px,y,py,z,δ) × (emx,emx,emx,emx,emz,emz)` → `tmap(xa,−1)` (= temw の `r` を掛ける: 正規モード→物理)
→ `+ codin` (emit の閉軌道) → `tconv(−1)` (px,py を (1+δ) で割る追跡内部表現、dv を δ から計算し dvfs を足す)。
DampEnough 用の閾値 `aenox/y/z = (振幅²)·em²·dampenough` もここで決め、ivar3 の軸は `big` (判定しない)。

### 4.4 損失判定 (tturn / tapert)

- tturn は最初の要素、その後 15 要素ごと (`la1=15`)、周末に `tapert` を呼ぶ (ax=ay=0、`xl=yl=−LOSSAMPL`, `xh=yh=LOSSAMPL`)。
- tapert 冒頭の `limittrack`: x,y を ±3·LOSSAMPL、px,py を ±0.9999、z を ±1e10 にクランプ (NaN も境界値に)。
- 損失条件は `.not.(−alost<x<alost .and. −alost<y<alost .and. |z|≤LOSSDZ)` — NaN を落とすため否定形で書かれている。
- 落ちた粒子は末尾の生きている粒子と入れ替え (`kptbl(*,1:2)` の写像も更新)、np を減らす。
- RFSW でなければ tturn の冒頭で毎周 `z(1:np)=0`。

### 4.5 tpdamp (Damping / DampEnough)

各周の後、`tconvm` で正準化 → `−codin` → `tmap(+1)` (物理→正規モード `ri`) し、3 面の振幅² がすべて aeno* 未満なら
`mturn=nturn−1` (次周で生存確定)。Damping なら `dampx=exp(−T0/τx)` 等 (T0=2π/omega0, τ は taurdx/y/z) を掛けて戻す。

### 4.6 集計・表示・ResultOfDAPERT

- 行 i の `nsc` = 最初に `ntloss(i,k) < nturn` となる k−1 (なければ 51)。`Score = Σ nsc`。
- 表示 (lfno): 見出し `Range ... / Display: muls turns/character / "     NZ     0----|----1----..."`、
  各行 `(1x,f8.2,i3,1x,51a1)` = a1i, nsc, `rad62a(ntloss/muls, j)`。rad62 は `0-9 A a B b ... Z z` の 62 文字 (範囲外 `*`)、
  rad62a は値 0 のとき j が 11,21,31,41,51 なら `'.'`、それ以外 `' '`。
- `ResultOfDAPERT = {{{a2min,a2max},{a3min,a3max},{a1min,a1max}}, {{a1i, nsc, {ntloss(1..51)}}, ...}}`。

## 5. pysad の対応 (dapert.py / ffscmd.fn_dasurvey)

| SAD | pysad | 備考 |
|---|---|---|
| dapert.n の既定・Mode 表・TURNS/NP の一時設定・"Turns = …" 行 | `fn_dasurvey` | "Turns =" 行と地図は Output->6 のときだけ出す (0 は捨てる) |
| tfdapert の iv の並べ替え、range の解釈 | `fn_dasurvey` → `range3=[[a2min,a2max],[a3min,a3max]]`, `r1` | n1 の上限は 201 (SAD 実効 200) |
| emx, emz | `fn_dasurvey` | EMITX/EMITY/EMITZ/SIGE は直前の emit (無ければ `run_emit(quiet=True)`) |
| trackd の起動順・座標 | `dynamic_aperture` | i=1..n1p、j=51→1 の順に全 n1p×51 粒子を一度に並べる (= npmax ≥ n1p·51 の SAD) |
| tinip1 | `rmat @ (amp*em) + codin` → `track6n.make_state` | make_state が tconv(−1) に相当 |
| tturn + tapert | `track6n.turn(..., apert=(kid, lost_turn, it, LOSSAMPL, LOSSDZ))` | `k_apert` = limittrack + tapert (ax=ay=0)。落ちた粒子は末尾と入れ替え、`kid` (配列位置→ID) も入れ替える |
| 損失周 | `lost_turn[ID] = it` (完了周数) → `ntloss` | 未損失は nturn |
| RFSW でない z=0 | `S[:n,4]=0` を周末に | 次周冒頭に相当 |
| 乱数 | `track6.RNG.seed` を配列で渡し、終了後に戻す | SeedRandom[] と同じ列 |
| rad62 / rad62a / muls / 書式 | 同名関数、f-string で `f8.2 i3` | 153 文字の地図が一致 |
| ResultOfDAPERT | `g['ResultOfDAPERT']`、戻りは `[score, [range, rows]]` | ReplacePart と同じ |
| 進捗 | `progress` 引数 / 環境変数 `PYSAD_PROGRESS` / グローバル `PYSAD$PROGRESS` | N 周ごとに生存数と経過秒を表示 |
| 落ちる粒子の乱数 | `track6n` の `sqrt1`/`pxy2dpz` (max(1e-200,·))、`asinz(NaN)=-π/2`、`norm2_3`、`gfloor` (NaN で gmin)、`rot1` (th2=0 は回さない)、`Prog.la_count` (DRIFT/MARK/MONI は la を減らさない) | 発散した粒子が tapert で除かれるまでの乱数の引き方 (θ=π/2 なら約 100 個、NaN なら 1 個) を SAD と揃える。SAD デバッグ版との再同期テスト (`<work-scratch>/r2_resync.py`) で 30 周分の seed が一致 (`history/work_loss_rng.md`, `history/work_loss_rng2.md`、docs/08 §7) |

**未対応**: `Damping->True` (SADError)、`DampEnough` (無視)、`ReferenceOrbit` (読まない)、`ExpandElementValues`、
DAPWIDTH<51 の窓 (常に行全体を起動するので、損失が多いときの起動順と NP 上限が SAD と異なる)、NPARA>1 の fork、
nw=16 周ごとの詰め替え (全粒子を同時に起動するので不要)、Output が 6 以外のユニット、tftrad の n1=0 経路。

## 6. 照合結果 (bench_level16*)

| ベンチ | 条件 | 結果 |
|---|---|---|
| 16a | ESRF 20 セル、NORAD/NOFLUC、`{{1,400},{1,1},{-0.5,0,0.5}}`、100 周、DAPWIDTH=51 | 地図 3 行×51 文字と Score 63、nsc=21×3 が一致。ntloss も track6n 版で一致 (z=0.5 行 j=35 の 1 点は損失境界で SAD 14、muls=2 なので文字は同じ `7`)。SAD 0.52 s (CPU)、pysad は numba 版 track6n で数秒 (初回は JIT 込み。純 Python の track6 では 6.5 分) |
| 16g / 16g_np1 | RAD+FLUC、`{{1,20},{1,1},{0}}`、20 周、SeedRandom[17] | 損失なし。NPARA=1 の SAD と最終 seed (221051865) まで一致。NPARA 既定 (4) の SAD は fork で乱数列が分かれ seed 3308107025 |
| 16f / 16f_np1 | RAD+FLUC、16a と同じ range、100 周 | Score 63・nsc は一致するが、損失粒子の周数 (地図の中央部) が一部違い、最終 seed も違う (SAD 3824086097、pysad 4221456353)。理由は下記 (カオス) |
| 16l_j22/25/27 | 単粒子を TrackParticles で 1,2,4,10 周 | 12 件すべて最終 seed 一致 (`history/work_loss_rng.md`) |
| 20 small (`bench_level20_dasv_small.sad`) | esrf リング、RAD+FLUC、`{{0,400},{1,1},{-0.5,0,0.5}}`、100 周、SeedRandom[17] | Score 63・nsc=21×3 一致。1・2 周目の損失周は SAD と完全一致、3 周目以降は一部違い、最終 seed は不一致 (SAD 966470477、pysad 620061457)。再同期テスト (SAD デバッグ版との要素ごとの突き合わせ) は 30 周・93 030 記録で seed 全一致 |

- 損失周・最終 seed が一致しない理由 (`history/work_loss_rng2.md`): 落ちる粒子は発散している間も光子の乱数を引き続け (θ=π/2, anp≈94 なら約 100 個、NaN なら 1 個)、
  発散の途中の値は丸め誤差の 1 ulp で別の分岐 (有限か inf/NaN か) に入る。SAD 自身も 1 周目の初期座標に相対 1e-16 の摂動を入れるだけで最終 seed が変わる (12 回中 2 回)ので、
  pysad (numba) と SAD (gfortran -Ofast + glibc) の演算順序・libm の差 (1 要素あたり 1〜数 ulp) が 600 要素で 1e-11 に増えると、落ちる粒子 1 個の分かれ目が逆になって以後の seed が全部変わる。
  NaN の扱いの論理 (`rot1`、`gfloor`、`xsincos`) は再同期テストで SAD と揃えてあるので、残りは丸めの運。判定は再同期テストの seed 一致で行う。
- esrf の `DAsv[]` (modules_general.n): DAPWIDTH=51、z を `Range[-4.5dp, 4.6dp, dp/2]` の 19 点、5000 周、x と y の 2 回。
  1 回分は 19×51 = 969 粒子 × 5000 周 × 要素数 (ESRF 20 セルで数千) ≈ 1.7e10 粒子・要素。176 ns/粒子・要素なら **約 50 分** (損失で粒子が減る分は短くなる)。

## 7. 使い方

```
cd "<repo>"
python run_sad.py pysad/tests/kek-ls_v03_12.lat bench_level16a.sad        # ラティス + スクリプト (SAD の gs と同じ連結)
set PYSAD_PROGRESS=100 && python run_sad.py lat.lat da.sad   # 100 周ごとに "[dapert] turn n/N alive k/M (s)" を表示
```

SAD 側は WSL で `<SAD>/bin/gs bench_level16a.sad`。RAD+FLUC の参照を取るときは
スクリプトに `NPARA = 1;` を入れる (既定 4 では fork で乱数列が分かれ、出力も重複する)。
スクリプト例 (bench_level16a.sad の FFS 部):

```
FFS; USE RING; CELL; CALC; RFSW; NORAD; NOFLUC; emit;
DAPWIDTH = 51;
res = DynamicApertureSurvey[{{1,400},{1,1},{-0.5,0,0.5}}, 100, Output->6];
Print["SCORE ", res[1]];  Print["RES ", res[2]];  Print["RDAP ", ResultOfDAPERT];
```
