# 作業ログ: ROTATE (ROT/TILT) 引数への全面対応 (2026-10-16)

## 1. 結論

SAD の `ROTATE` (= `ROT` = `TILT`) を BEND・QUAD・SEXT/OCT/DECA・MULT のすべてで、CALC (4D)・CALC6D・Emittance[]・TrackParticles・GEO の全経路で SAD K64 と一致させた
(参照ベンチ bench_level25a〜f, 25h と 1 要素の写像テスト)。

| 要素 | 状況 |
|---|---|
| (1) ANGLE=0 の BEND に ROTATE=−90 DEG (垂直ステアリング VD/PV/VS) と K0 | parser が Drift にしていたのをやめ、`Bend` (ANGLE=0) として保持。CALC・emit・追跡とも一致。pf2 の v9_25 (VD L=.15, PV/VS L=0) でも一致 |
| (2) skew 四極 (QUAD の ROTATE + K1) | 元から一致 (tsolrote の経路)。**L=0 の QUAD (pf2 の SQ: `L=0 K1 ROTATE=45 deg`) は K1 が無視されていた** → tquade の al=0 → tthine(nord=4) の薄い四極として CALC・emit・追跡に追加 |
| (3) ROTATE 付きの BEND (ANGLE≠0, 5 DEG / 90 DEG / K1 付き / K0) | 一致。回転とずれ・DROTATE の組み合わせ順を tchge の直訳に直した (ずらす→回す、行列の作り方まで) |
| (4) 六極・八極・十極・MULT の ROTATE | 一致。ずれ DX/DY との組み合わせで SAD の癖 (後述 §3) があり、それも再現 |

付随して直したもの: FFS の `LINE["ROTATE",BEND]` の読み書き (Bend.ROTATE)、`DTHETA` の格納先 (BEND/MULT は DROT)、FFS で ROTATE/ANGLE を変えたときの GEO、
`emit_extra`/`touschek` に複製されていた要素ループ、numba の ANGLE=0 BEND (tsteer) カーネル、L=0 の QUAD の追跡 (`NotImplementedError('長さ 0 の四極')` の解消)。

## 2. SAD の仕様 (src で確認)

- **ROT の意味**: `ROTATE`=`ROT`=`TILT` は 1 つの値 (`ky_ROT_*`)。BEND だけ別に `DROTATE` (`ky_DROT_BEND`)。MULT は `ROT` + `DROT` + `akang` を `theta2` にまとめて tsolrote の chi3 に渡す (tffs.f:1398)。
- **回転の向きと順序** (tchge.f): 入口は「`theta` で (x,px),(y,py) を回す (x'=c x − s y, y'=s x + c y) → 回した座標系のずれ (c·dx−s·dy, s·dx+c·dy) を引く → tbrote」。
  つまり **ずらす→回す** (ずれは実験室系で与える)。出口はその逆。四極の tsolrote も同じ順 (x −= dx; y −= dy; 回す)。BEND は ANGLE≠0 なら `th=theta` (DROTATE は tbrote で BEND の中点まわり)、phig=0 (ANGLE=0) なら `th=theta+dtheta`。
- **tbende.f の分岐ごとの theta の渡し方** (`lattice._bend_chg_args`):

  | 条件 | 呼ばれるもの | tchge の theta |
  |---|---|---|
  | ANGLE≠0 | tchge(dx,dy,dz,theta,dtheta,dchi2, alg=L/2, phig=ANGLE/2) | theta (DROTATE は tbrote) |
  | ANGLE=0, K1=0 | tsteee (CALC/emit) / tsteer (追跡) | CALC: **theta だけ (DROTATE は渡らない)**、追跡: **theta+dtheta** |
  | ANGLE=0, K1≠0, K0=0 | tquade / tquad | CALC: theta、追跡: theta+dtheta |
  | ANGLE=0, K1≠0, K0≠0 | tmulte / tmulti | theta+dtheta (+akang) |

  tsteee/tthine は **K0=0 (ak=0) だとドリフトを返して tchge を呼ばない** (ずれ・回転とも効かない)。L=0 の ANGLE=0 BEND は tthine(nord=2, ak=K0) で `px −= K0` (CALC)。
- **DTHETA の格納先** (terror.f): DROT のある型 (BEND, MULT) は DROT、なければ ROT (QUAD, SEXT...)。
- **tquade の al=0**: `tthin(nord=4, al=0, ak=K1 (積分値), dx, dy, theta, krad=.false., fringe=.false.)`。K1<0 の −90° 回転 (akang) は入らない。

## 3. SAD の癖: tchge の入口で回転が行列に入らない (再現済み)

tchge.f の入口で `th≠0` のとき回転行列 `trans1` を作って `transa=>trans1` とし、続けて `dx_rot≠0 (または dz≠0, dtheta≠0)` なら **同じ `transa` を tbrote に渡す**。
tbrote → tbgrote (および tdrife1) は `trans` を「掛ける」のではなく「上書き」する (tinitr してから代入) ので、**入口では COD は回るが行列には回転が入らない**。
出口 (`th≠0`) は `tmultr5(transa,trans1)` で掛けるので回転が入る。結果、一周行列は R^T·M (非対称) になり、**th≠0 かつずれ≠0 の要素 (BEND, 六極などの薄い多極) は結合 (R1..R4)・チューンが物理とずれる**。
R1 が 1e-4 → 0.26 になるほどで、ESRF のリングに rotate 30° + DX=1e-3 の六極 1 台を入れるとチューンが 0.01 動く。

- 再現: `lattice._tchge_jac` (入口で th≠0 かつ tbrote を呼ぶとき、tbrote だけの行列 = 単位行列 ((5,6) だけ dpzi(pr+1)/pzi) か、phig≠0 なら tdrife1 (長さ dz1=dxr·sin φ の BEND 中点ドリフト) を含む合成の数値微分)。
  th=0 のときは普通に (回転なしで) tbrote の行列。1 要素の写像テスト (bench_level25g) で BEND (K1=0/K1≠0)・六極とも 1e-11 で一致 (この癖を入れる前は 1e-7)。
- 適用先: BEND (`_bend_element_cod_matrix`、emit6d/calc6d の `_tbende`)、薄い多極 (`lattice._thin_cod_matrix`、emit6d `_tthine`/`_tthine_gen`)。四極・MULT の tsolrote には癖はない。
- 追跡 (TENT/TEXIT, tbshift) は座標だけなので影響なし。
- **帰結**: この癖のある要素 (th≠0 かつずれ≠0) が入ったリングでは SAD の閉軌道 (qcod は itmax=10、tcod は緩い収束) が収束しないことがあり、**Emittance[] の δ・チューンは SAD と一致しない** (CALC の Twiss と追跡は一致する)。
  例: 六極 1 台 (ROTATE=30 DEG DX=1e-3 DY=−5e-4) だけで SAD の δ=−1.2e-9, z=−1.8e-7 (z が閉じない) に対し pysad は z を閉じる δ=−1.3e-6 (通常のリングでは SAD も z を閉じて一致する)。
  K0 のキックで軌道が 1e-4 に達するリングでは CALC の COD も itmax=10 の打ち切りで 2e-5 ずれることがある。この組み合わせの Emittance[] は参考値 (bench_level25f)。

## 4. SAD が落ちる入力

- **`BEND` に ANGLE=0 かつ K1≠0 かつ K0≠0 を同時に与える** (L>0 でも ROTATE の有無にも依らない) と、SAD K64 (WSL ビルド) は FFS の初期化 (`USE RING` 後の最初の CALC 前) で **SIGSEGV (exit 139)** する。K1 だけ (K0=0)、K0 だけ (K1=0) は動く。
  → 参照値が取れないので、ANGLE=0, K1≠0, K0≠0 の BEND の SAD 照合は不可能 (work_bend_k0.md §7 の項目はこれで「SAD 側が不可」と確定)。pysad は tmulte (CALC/emit) で計算し、追跡は NotImplementedError のまま。

## 5. 変更点 (ファイル別)

- `parser.py`: `ANGLE=0 かつ |ROTATE|>0.1` の BEND を Drift にしていた分岐を削除。`Bend(ANGLE=0, ROTATE=…)` として保持 (`_sad_type='BEND'` の仕組みは他の用途のために残す)。
- `elements.py`: `Bend.constant_kick` (ROTATE を垂直キックとして足していた基底クラスの `d[3] += ROTATE` を継がない。旧 ffs.py の線形 COD 用)。
- `lattice.py`:
  - `_bend_theta` (ROT = ROTATE+ROTATION)、`_bend_chg_args` (§2 の表)、`_bend_noerr` (誤差なしのコピー)、`_tchge_cod` / `_tchge_jac` (tchge の直訳)、
    `_bend_core_cod_matrix` / `_bend_element_cod_matrix` (誤差なしの中身 + 外側の tchge)、`_thin_cod_matrix` (薄い多極 + tchge)。
  - 削除: `_roll_matrix` (ROTATE を外側の R6 サンドイッチにしていた。DX/DY との順序が違った)、`_bend_err_wrap`。
  - 薄い四極 (L=0 の QUAD): `_one_turn_cod_matrix` の Quad 分岐と `_has_offset` (`transfer_matrix` は L=0 で恒等なので近道しない)。
- `emit6d.py`: `_tbende` を「tchge (入口) → `_tbende_core` → tchge (出口)」に。`_tthine` / `_tthine_gen` (薄い多極と薄い四極)、`_tquade` の L=0。
  `calc6d.py`, `emit_extra.py`, `touschek.py` は要素ループの複製を `_tthine` 呼び出しに置き換え (複製が古い `_thin_offset_map` を使っていたため、癖を入れても Emittance[] に届かなかった)。
- `track6.py`: BEND の theta = ROTATE+ROTATION (+ tsteer は DROTATE)、L=0 の四極 (`_tquad_ens`: tthin(nord=4))。
- `track6n.py` (numba): `K_STEER` カーネル (tsteer: ANGLE=0 BEND の L>0, K0≠0 で、ずれ・ロール・放射込み。以前は track6 へフォールバック)、`K_STEER0` にずれ・ロール、L=0 の四極 (K_THIN の kord=1)。BEND の theta は ROTATE を含む。
  `track6n_ext.py` の step_kind 呼び出しを新しい引数に合わせた。
- `ffscmd.py`: `LINE["ROTATE",BEND]` は `Bend.ROTATE` (他の型は ROTATION)、`DROTATE` は DROTATE 属性、`DTHETA` (BEND/MULT → DROTATE、他 → ROTATION で、設計値 + 誤差)、`_dump`、
  FFS で ROTATE/ANGLE を変えたとき `el._raw` も更新 (geo.py は `_raw` を先に見る。以前は CALC の GEO が古い値のまま)。
- `match6d.py`: BEND の qdbend の theta に ROTATE を足した (解析微分。VD は Bend になったので回転つきで微分される)。
- `bench_level12_pysad.py`: 読み込みの ROTATE 照合に Bend.ROTATE を足した (VD が Drift ではなくなったため)。
- 新規テスト `pysad/tests/test_rotate.py` (16 件: parser、tchge の行列と癖、垂直ステアリングの向き、25a〜d の CALC を SAD の .result と照合)。

## 6. 参照ベンチ (SAD K64 を WSL で実行して .result を保存、pysad は `python run_sad.py bench_level25X.sad` → `_pysad.result`)

ラティスは esrf (bench_level23f の ARCA に要素を挿入) と pf2 v9_25。比較は `tmp25/cmp3.py` / `tmp25/summary.py` (Print の数値行を突き合わせ、絶対差)。生成は `tmp25/make25*.py`。

| ベンチ | 内容 | 結果 |
|---|---|---|
| 25a | `VL=(L=.15 ROTATE=-90 DEG)` と `PV=(L=0 ROTATE=-90 DEG)` を各セルに挿入。`LINE["K0","VL.1"]=1e-4`、`PV.1`=−2e-4。CALC の DY/DPY/EY/NY ほか、LINE["GX"/"GY"]、Emittance[]、TrackParticles NORAD/RAD 1,20 周 | CALC 9e-10 (3 段)、GEO 4e-13、emit εy 2.7e-14 (相対 1e-8)、TUN 1.6e-8、追跡 7e-16 / RAD 3e-13 |
| 25b | skew 四極 `QS1` (ROTATE=1 DEG)、`QS2` (K1<0, −2 DEG)、`QS3` (FFS で ROTATE=.7 deg)、ほか QDA3 の ROTATE を FFS で変更。R1..R4, DETR, 全 Twiss 列 | CALC 2e-12、εy 相対 2e-12、TUN 5e-14、追跡 1e-14 |
| 25c | ROTATE 付き BEND: 5 DEG の対 (BR1 は E1/E2=.5、BR2 は ANGLE 逆符号)、90 DEG (垂直偏向、E1/E2) の対、K1 付き (3 DEG / −70 DEG)。T0、GEO (GX,GY,GZ,GCHI1..3)、emit、FFS で ROTATE と K0 を変更 (T1) | CALC 6e-10、GEO 8.5e-14、emit 相対 8e-10、TUN 4e-9、追跡 1e-15 / RAD 4e-13 |
| 25d | SEXT (L>0 ROTATE=10/30 DEG、L=0 −20 DEG)、OCT (L>0, L=0)、DECA、MULT (K1/K2/K3 + ROTATE、L=0 の K0 + DX/DY + ROTATE)。FFS で ROTATE と DX/DY を変更 | CALC 7e-9 (BX≈10)、emit 相対 3e-12、TUN 1.5e-11、追跡 3e-16 |
| 25e | 組み合わせ (癖が出ない組): ROTATE+K0 の BEND、DX/DY+K0 の BEND、K1 付き BEND の ROTATE=−40 DEG +K0、DX/DY+DROTATE、DROTATE+E1/E2、垂直ステアリング (L>0, L=0) + K0、水平ステアリング + DX/DY、QUAD ROTATE+DX/DY、MULT K0+ROTATE+DX/DY。FFS で変更した T1 も | CALC 3.6e-9、GEO 5e-13、追跡 1e-16、emit 相対 1.6e-6 / TUN 4.5e-8 (δ の収束の差: 25c/23g と同じ程度) |
| 25f | 癖が出る組 (th≠0 かつずれ≠0): BEND (ROTATE=20 DEG, DX/DY, K0) の対、垂直ステアリング (L>0, L=0) に DX/DY + K0、六極 (ROTATE=30 DEG, DX=1e-3, DY=−5e-4)。CALC の R1..R4 ほか全列 | CALC 9.6e-10、追跡 1e-16。Emittance[] は §3 のとおり一致しない (参考値) |
| 25h | pf2 v9_25 (USE RING): VD (L=.15)、PV (L=0)、VS (L=0) の K0 (ファミリー展開) を 1 つずつ + 3 つ同時、SQ (L=0 skew 四極) の K1、追跡 | CALC 6.4e-10 (V, P, S, M)、SQ 3.5e-12、追跡 1e-16、emit 相対 1.5e-6 / TUN 3e-8 |
| 25g | 1 要素の INS (`tmp25/maptest.py`): 入口に大きな COD (DX=1e-3, DPX=5e-4, DY=−5e-4, DPY=2e-4) を与え、出口の COD・Twiss・R1..R4 を比較。13 種 (BEND/K1 付き BEND の ROTATE/K0/DX/DY/DROTATE、垂直ステアリング、QUAD、SEXT、MULT) | すべて ≤ 5e-11 (出力 10 桁の表示桁以内を含む)。癖を入れる前は BEND (ずれ+ROTATE) が 1e-7 |

<work-pf2> の垂直バンプ系 (`grep -l "VS\|PV\|VD" <work-pf2>/*_r`) は、実際に垂直補正コイルへ K0 を入れて動かす SAD 出力が揃っているスクリプトがなかった (correct_cod_9_16 は `dtny` に K0 を入れるが、以前の記録どおり別の理由 (CalculateOptics, 長さ 0 の四極の追跡) で pysad が止まっていた。今回 L=0 の四極の追跡は入れたが、CalculateOptics の未実装は残る)。
代わりに pf2 の lattice (v9_25) で 25h を作り、VD/PV/VS の K0 と SQ の K1 を直接照合した。

## 7. 回帰

- 作業前に 13a/13b/14a/16a/16b/18a/19a/22a/23a の pysad 出力を `base_pre/` に保存し、作業後 (`base_post/`) と diff: 13b/14a/16a/19a/22a/23a は完全一致 (経過時間・警告の行番号を除く)、
  16b は 1e-11、18a は 1e-10 → 1e-13 (ESRF の VD が Drift から Bend になり、理想のリングの Twiss が SAD に近づいた)。
- `QUICK=1 PYTHONIOENCODING=utf-8 python run_regression.py`: 28 ベンチ **FAIL 0** (3c の EMITY も PASS)。QUICK で飛ばす 7, 12t, 13a, 16a も実行 → 7: PASS (5e-11)、12t PASS、13a PASS* (13A-2 の既知)、16a PASS。
  (最初の実行で [12] の「ラティス読み込み ROTATE」が FAIL: VD が Bend になり ROTATE の照合式 (`ROTATION`、Drift は除外) が合わなくなったため。照合式を Bend.ROTATE に直した。)
- `python -m pytest -q pysad/tests --ignore=pysad/tests/test_bend_feeddown.py --ignore=pysad/tests/test_coord.py --ignore=pysad/tests/test_sext_radiation.py`: 85 passed → test_rotate.py を足して **101 passed**。

## 8. 未解決・既知

- th≠0 かつずれ≠0 の要素を含むリングの Emittance[] / COD (§3)。CALC の Twiss と追跡は一致。
- ANGLE=0, K1≠0, K0≠0 の BEND: SAD が落ちるので照合不可 (§4)。追跡は NotImplementedError。ANGLE=0, K1≠0, K0=0 の BEND の追跡 (tquad の theta+dtheta) も NotImplementedError のまま (CALC/emit は tquade で実装済み、未照合)。
- CAVI/TCAV/MONI の ROTATE: SAD は ROT を持つ (`ky_ROT_CAVI` …) が pysad は無視 (CALC/emit/追跡)。ソレノイドの ROTATE (trots) も同様。
- ANGLE=0 の BEND の F1/FB1/FB2 フリンジ (tsteee/tsteer) は NotImplementedError のまま。
- Emittance[] の δ は SAD の tcod の収束 (z の 1e-7 打ち切り) の分だけ揺れる (25e/25h の TUN 3e-8〜4.5e-8、εz 1.5e-6 相対)。回帰の許容 1e-6 をわずかに超えるので、25e/25h は「δ の収束差」として記録。

## 9. 作業用ファイル (`tmp25/`)

`gen25.py` (23f のラティスに要素を挿入する生成関数)、`make25a〜f.py`、`bis.py` / `bis*.py` (SAD と pysad を両方走らせて要素単位で二分する)、`maptest.py` (§6 の 25g)、
`dbg_*.py` (一周行列・COD の残差)、`cmp.py` / `cmp3.py` / `summary.py` (比較)、`runsad.sh` (SAD の実行: PowerShell から `wsl -d Ubuntu -- bash -lc` でクォートが外れるので `;` でつなぐ)。
