# 10. ベンチマーク一覧 — SAD (K64) と pysad の照合

最終更新: 2026-11-06 (現状のスナップショット)。

ベンチは `<repo>/` の `bench_level*.sad` (SAD スクリプト) / `.result` (SAD の出力) / `*_pysad.py` (pysad の比較スクリプト) / `*_run.result` (pysad の出力)。
回帰テストの実行スクリプトは `run_regression.py` (§4)。レベルの階層 (Mainlevel / FFS の構造) は `<work-olddocs>/sad_level_structure.md`。

場所: すべて `<repo>\` (以下「`<repo>`」)。pysad 本体は `<repo>\pysad\`。

---

## 1. 方針

- SAD (K64 fork, `<SAD>/bin/gs`, Ver.1.4.4.2k64) を WSL Ubuntu で実行して参照値 (`bench_levelN.result`) を取り、
  pysad を **SAD の Fortran の直訳** として作り、同じ入力で同じ数が出ることを確かめる (再導出はしない)。
- 一致の目安: 追跡 (TrackParticles) は座標の絶対差 ≤ 1e-9 (SAD の印字精度 1e-10〜5e-11)、オンモメンタムの CALC は 1e-9、
  Emittance[] は相対 1e-6、DA は地図 (文字列) と Score の一致。閉軌道・結合・オフモメンタムを経由する CALC は
  SAD の qcod の緩い収束判定 (conv = min(1e-10, …)) と pysad の数値微分の丸めが入るので相対 1e-7 を目安にする (§3.1)。
- レベルは「SAD の機能」の階層で、下のレベルから順に一致させる。

---

## 2. 実行方法

### 2.1 SAD (参照値の生成)

```bash
wsl -d Ubuntu -- bash -lc 'cd "<repo>" && <SAD>/bin/gs bench_levelN.sad > bench_levelN.result'
```

- `bench_levelN.sad` は **ラティス (Mainlevel の素子定義・LINE) と FFS 部 (`FFS;` 以降の SADScript) を 1 ファイルに連結したもの**。
  実ラティス (Level 12 以降の ESRF v03_12) では `pysad/tests/kek-ls_v03_12.lat` の内容を先頭に貼り、`FFS;` の後にスクリプトを書く。
  SAD の `gs` は複数ファイルを連結して渡すのと同じ (`gs lattice.lat script.sad` でもよい) が、リダイレクトで参照を取るときは
  `gs lat script` の 2 引数が標準入力待ちで止まるので 1 ファイルに連結する。
- **`NPARA=1` を FFS 部の先頭に入れる。** SAD の `NPARA` の既定値は `Packages/init.n` で 4 で、`TrackParticles` は粒子 9 個以上
  (npparamin=9) で fork して乱数列が分かれ、また fork した子プロセスが出力を重複して書く。参照値を取るときは NPARA=1 にする
  (`bench_level16g_np1.sad` など `_np1` が付いたものは NPARA=1 版)。
- 出力桁数は `$FORM="15.10";` (10 桁)。`PageWidth = 1999;` で折り返しを抑える。Print の出力は `Print["TAG ", {...}]` の形にして、
  結果ファイルの行頭 `TAG {...}` を pysad の出力と突き合わせる。Level 24 系の参照は WSL の gs を pty 経由 (`script -qec`) で取る。
- FFS で使う LINE の先頭要素は MARK でなければならない (`tffs.f:2989`)。`LINE RING=(10*(CELL CELL2))` のような入れ子の繰り返しは構文エラー
  (`LINE CELLP=(CELL CELL2); LINE RING=(10*CELLP);` にする)。Do[] などの中で FFS コマンドを使うときは `FFS["cell calc;"]` で包む
  (出力ユニットを省略すると表示なし)。
- SAD (WSL 版) は FITP>1 の GO (bench_level15a, esrf/chrom) で `*** Run time Environment` からスクリプトを再実行してしまい
  参照値が取れない (NPARA=1 でも同じ)。
- SAD のクラッシュ (SetElement 周辺の SIGSEGV) に注意: ファイルへのリダイレクトで `SetElement[QF,QUAD,K1->0.3]` が落ちる。
  SAD のバグ・未初期化のためベンチから除いてあるもの: `Delete[l,0]` (SIGABRT)、`Level[..,{-1}]`、`DeleteCases` の第 4 引数、
  `Insert[l,x,0]`、stabz/idim=3 の stab・trace。
- 誤差は FFS コマンド `DELX/DELY/DK/DDK/DTHETA σ 要素` (QD.1 で 1 台、QD* で全台に別々の乱数、SeedRandom[n] で種固定、DUMP で一覧)。
  `LINE["DX","QD.1"]=v` → CALC では **全 QD に展開** される (QD.1 はファミリーの代表)。

### 2.2 pysad

Windows (Git Bash / PowerShell)、`<repo>` で:

```bash
PYTHONIOENCODING=utf-8 python -u bench_levelN_pysad.py            # Level 1〜12 の比較スクリプト (SAD の値を内蔵/または .result を読む)
PYTHONIOENCODING=utf-8 python -u bench_level12_track_fast.py      # 追跡の numba 版 (track6n)
PYTHONIOENCODING=utf-8 python -u run_sad.py bench_level13a.sad    # SAD スクリプトをそのまま実行 (Level 13 以降)
PYTHONIOENCODING=utf-8 python -u run_sad.py pysad/tests/kek-ls_v03_12.lat <work-esrf>/chrom   # ラティスとスクリプトを連結 (SAD の gs と同じ)
PYTHONIOENCODING=utf-8 python -u run_regression.py                # 回帰テスト (全部)。QUICK=1 で遅いものを省略
```

- `bench_levelN_pysad.py` は pysad の API (`Lattice`, `twiss_table`, `emit6d.emittance6d`, `track6`, `match6.Matcher` …) を直接呼び、
  SAD の値 (スクリプトに埋め込み、または `bench_levelN.result` から読む) との差を印字する。出力は `bench_levelN_pysad.result` に保存してある。
- `run_sad.py lattice script` は SADScript インタプリタ (`pysad/sadscript.py`) + FFS コマンド層 (`pysad/ffscmd.py`) で
  SAD スクリプトをそのまま実行する。MAIN 部 (素子定義) は `parser.py`、`FFS;` 以降は上記。出力は `bench_levelN_run.result`。
- 遅いもの: Level 7 (K1≠0 BEND の純 Python 追跡、約 2 分)、`bench_level12_track_pysad.py` (純 Python、約 1 分)、
  `run_sad.py bench_level13a.sad` (13A-2 の解なしの反復が 3〜4 分)、`bench_level16a.sad` (DA、numba、数十秒〜)。
  `run_regression.py` は `QUICK=1` でこれらを飛ばす。
- 比較の補助スクリプト: `cmp_tags.py` (タグ行・ブロックの比較、Level 24)、`cmp_lines.py` (Print のラベル付き数値行を順に比較)、
  `cmp_emit_report.py` (EMIT 報告の行単位比較)、`cmp_radint.py` (RADINT 出力ブロック)、`cmp_blocks.py --noise 1e-9` (DISP などのブロック比較。
  SAD の 1e-17 などの雑音と pysad の 0 は同じ扱い)。pf2 の matching の反復列比較は `<work-pf2>/cmpfit.py` (NEWTON/Matched/Residual の行だけを抜き出して比べる) と `*_r.sad.out`。

---

## 3. 全レベルの一覧

「現在の結果」は `run_regression.py` の実行値と各ベンチの `.result` / `_pysad.result` による。
SAD 側の表示桁 (10 桁、追跡は 5e-11) が比較の下限。後ろの 2 列は、そのベンチが照合している pysad の実装と SAD の仕様 (癖を含む)、および未照合の項目。

| Level | 内容 (何を照合) | スクリプト (SAD / pysad) | 現在の結果 (一致度) | 関係する pysad の実装・SAD の仕様 | 未照合の項目 |
|---|---|---|---|---|---|
| 1 | 20 セル FODO (2.5 GeV, 140 m) の線形光学 (入口・QD.1 の β, α, η, Q)、COD (QD DX=1e-4)、RADINT (I1〜I5, εx, U0)、Emittance[] (εx, σE, α, U0) | bench_level1.sad / bench_level1_pysad.py | Twiss 表示 6 桁で一致、COD x −110.8245 µm (一致)、I1〜I5 7 桁一致、Emittance[] 相対 ≤3e-9 | periodic_optics は周期分散を含む。INTGRL.f (SYMPS, BEDGE/TEDGE) の直訳。`LINE["DX","QD.1"]` は CALC で全 QD に展開 (SAD の既定): `set_element` (ファミリー) と `set_error` (1 台) を区別。SAD の MomentumCompaction は Emittance[] (6D, −R56/C) の値で RADINT の I1/C とは別 | — |
| 2 | オフモメンタム Twiss (DP0 = 0, ±0.01, ±0.02)、K1 付き BEND (RINGK) の線形光学・RADINT・Emittance[] | bench_level2.sad / bench_level2_pysad.py | Twiss 表示 6 桁で一致 (β 相対 ≤3e-8)、RADINT 一致、Emittance[] 相対 ≤5e-9 | 分散は qcell.f の δ 列 (a13, a23) から周期解を出す式。四極のハードエッジフリンジ tqfrie、tsolque の分割 ndiv=1+int(\|L·K1L\|/0.1)、qcod の直訳 (収束判定・ダンピング・初期値は直前の閉軌道 twmov)。SAD の off-momentum β は初期値で相対 1e-6 動くので SAD と同じ順で呼ぶ | 不安定な dp (RINGK dp=−0.01) で SAD は 1/cosμ の値を返すが pysad は例外 (qcell の intres/halfres 未移植) |
| 3A | 六極 SF/SD と色収差: DP0 スキャン 21 点のチューン、3 次当てはめ (補正前・SAD の補正後 K2) | bench_level3a.sad / bench_level3a_pysad.py | チューン 21 点 ≤1.5e-10、ξ (c1) 差 ≤1.5e-6 | 薄レンズ多極の COD 写像 tthine (L/6–2L/3–L/6)。SAD 組み込みの CHRO/RADINT の色収差は BEND 成分を含まないので使わない。SAD の PolynomialFit (Packages/Fit.n) は最高次以外の係数が壊れるので、ここでは自作 MyPolyFit で当てはめる | 色収差補正 (FIT) の K2 は 5e-4 の差 (評価関数の違い、ユーザー判断で許容: 同じ K2 を入れれば結果は一致) |
| 3B | 六極を通る COD (DELX QD.1, DELY SF.1 を 1 台だけ)、COD 下の β/α/η/Q/EY/R1–R4、色収差 | bench_level3b.sad / bench_level3b_pysad.py | 全項目が許容内 (閉軌道・Twiss・Q・R1–R4・ξ、最大 4e-13) | tchge.f (ずれ・回転)、結合時の分散 (R1–R4 と c=√(1−det R) で対角化)。ドリフトなしで続く要素は上流と一緒に動くので QF–SF 間に DQS=0.1 m を置く | 結合時の途中要素の位相 |
| 3C | x-y 結合 (DTHETA: QF.3 ROTATE, B.5 DROTATE)、閉軌道・光学・R1–R4・εx, εy, σE, U0 | bench_level3c.sad / bench_level3c_pysad.py | 閉軌道 相対 ≤1e-8、β/α/Q ≤1e-8、EX/EPX 3.2e-11 / 6e-12、Emittance[] 相対 ≤7.8e-7 | 回転の向きは x′ = cosθ x − sinθ y に統一、四極のロールは「ずらす→回す」、BEND のロールは中点まわり (LGEO=L/2, tfgeo.f)。6D 行列の z 行はシンプレクティック条件で組む (数値微分のままだと 6.6e-5 非シンプレクティックで εy が 2% 狂う)。QD K1 は −1.7×0.3 (Qy=8.0095 だと整数共鳴に近く σE が NaN) | ずれ・ロール付き BEND の行列は数値微分のまま (`history/work_ins_offset.md`) |
| 3D | 強さの誤差 (DK QF.2, DK QD* の乱数, DDK B.3)、閉軌道・光学・εx, σE, U0 | bench_level3d.sad / bench_level3d_pysad.py | 閉軌道 相対 ≤4e-8、Twiss 相対 ≤2e-8、Emittance[] ≤3.3e-8 | DK は既定キーワード (QUAD は K1) への相対誤差、DDK は K0 への絶対誤差 (terror.f)。BEND の K0 (phib = ANGLE + K0, tbende.f) | — |
| 4 | RF 空洞 (RFSW ON): (1) RADCOD OFF、(2) RADCOD ON の Emittance[] (εx, εz, σE, σz, Qx/y/z, 減衰率, U0, Vc, trf0, α, バケット高さ, 放射込み閉軌道) | bench_level4.sad / bench_level4_pysad.py | (1) 相対 ≤2.1e-9、(2) 相対 ≤1.1e-7 (α_c 8e-10)、閉軌道 x の差 1.5e-10 m (SAD tcod の収束誤差内) | tcave.f (長さ 0 の空洞)、temit.f tecalc の synchm (縦モード)、photons.f tradke の radcod 分岐。RF なしでは z が閉じる δ を選ぶ。SAD は z のずれを trf0 側へ移すので z と trf0 の分け方が違う (z/β+trf0 は一致) | — |
| 5 | 輸送ライン (INS) とフリンジ: A) F1=0.1, FRINGE=3、B) F1=0、C) DISFRIN=1。終端と要素 6 の Twiss・軌道・位相進み | bench_level5.sad / bench_level5_pysad.py | 軌道 ≤1.3e-11、Twiss 絶対 ≤2.6e-7 (β 〜10 m)、NX/NY ≤2.7e-9 | ins_twiss (入口の Twiss・軌道から非線形写像で伝搬)。tqlfre (線形フリンジ)、tbendef1 (BEND の F1)、tsetfringepe。入口の Twiss・軌道は先頭 MARK のパラメータで与える。Do[] の反復子に文字列リストは使えない、USE の引数に変数は使えない | K1≠0 BEND の F1、逆向き要素の F1/F2 |
| 6 | 粒子追跡 TrackParticles (5 粒子 × 1/10/100 周): RAD OFF、RAD ON (FLUC OFF) | bench_level6.sad / bench_level6_pysad.py | ≤5e-11 (SAD の印字精度) | 追跡は CALC/emit と別ルーチン → `pysad/track6.py` に直訳 (tdrift, tquad/tsolqu, tthin+ttfrin, tbend/tbendr, tcav, tradk)。内部座標 px/(1+δ) と粒子ごとの dv | FLUC → Level 10、K1≠0 BEND・誤差 → Level 7、MULT/SOL → 8〜11 |
| 7 | 追跡: K1=−0.05 の複合機能 BEND のリング (7-1)、誤差 DELX QD.1 / DELY SF.1 / DTHETA QF.3, B.5 / DELX B.2 (7-2) | bench_level7.sad / bench_level7_pysad.py (約 2 分) | ≤1.7e-10 (z のみ、trf0 の表示桁の差)、他 ≤5e-11 | tbendi.f (module bendib: rhoe は入口で 1 回だけ決める)、tsolrot、TENT/TEXIT.inc、tbshift/tbrot。DUMP は 5 桁なので LINE["DX",…] で 10 桁を取る。誤差の乱数は既定の種 17 の GaussRandom の列そのもの。emit6d: 誤差で周長が変わると δ≠0 の閉軌道 | — |
| 8 | 追跡: 四極の F1=0.05/FRINGE=3 (8-1)、長さ 0.5 m の空洞 (8-2, QF K1=1.76×0.3)、六極を MULT (K2) に (8-3) | bench_level8.sad / bench_level8_pysad.py | ≤5e-11 | tsetfringep の係数 (f1=−\|K1/L·F1²\|/24, f2=\|K1/L\|·F2)、放射ありでは tradk に f1in/f1out をそのまま渡す (SAD のまま)。tcav の L≠0 (vnominal=0 で分割 1)、tmulti (eps00m, ampmaxm)。emit6d も tcave L≠0・tmulte・tqlfre を持ち trf0 を自前で出す | MULT の ANGLE≠0 (tmulta)、SK |
| 9 | 追跡: ソレノイド領域 SOL BZ=2 T (BOUND=1, GEO=1) に D1, QD を含む。RAD OFF / ON | bench_level9.sad / bench_level9_pysad.py | ≤1.5e-10 (z, trf0 は SAD 表示値)、他 ≤5e-11 | tsfrin、tdrift_solenoid、tsolqu/tsolqur の bz 分岐、tfbzs (bzs = charge·BZ/(p0/c))。**SAD の癖: tsolqur (RAD ON, bz≠0) の最後の半ステップが最後に tzsetparams を呼んだ粒子の bzp を全粒子に使う** → 粒子の並びで結果が 1e-6 変わる。pysad は要素外側・粒子内側のループで再現 | SOL の DX/DY/DZ/CHI (trots)、F1、BOUND=0、SOL 内の BEND/MULT/CAVI、emit6d の SOL |
| 10 | 追跡: RAD + FLUC (SeedRandom[17])、1/10/100 周、追跡後の乱数の状態 | bench_level10.sad / bench_level10_pysad.py | ≤4.9e-11、乱数の最終状態 3 回とも一致 | SAD の乱数 (random_plugin_sad.c: 32 bit 合同法 seed←48828125·seed、Box–Muller、GCUT) を `pysad/sadrandom.py` に直訳 (15 桁一致)。FLUC の放射 (tradkfn/tradkf1, tdusrn/drndsr)。乱数を引く順: 要素外側・粒子内側、放射ありの BEND/四極/MULT はスライスごとに全粒子の本体→全粒子の tradk | — |
| 11 | 追跡: MULT の K0/K1/K2 とずれ・回転 (ANGLE=0)、閉軌道が δ=1.4e-3 にずれる場合 | bench_level11.sad / bench_level11_pysad.py | ≤5.9e-11、emit の trf0 2.6e-10 | akang (K1<0 で 90° 回転)、tbfrme、tsoldz、tsolque の dx0/dy0、ttfrin nord=2・ttfrins。**SAD の癖: 放射なしの MULT (tsolqum) の z の式は dv の項を括弧の中に入れており tsolqu と一致しない** (δ=1e-3 で 1 周 5e-10) → そのまま再現 (`_tsolqu(qum=True)`) | ANGLE≠0 の MULT、SK、MULT の F1/RF |
| 12 | 実ラティス ESRF 型 3GLS (pysad/tests/kek-ls_v03_12.lat: 3 GeV, 603.12 m, 3563 要素): 読み込み (名前・型・L・K1・ANGLE・K2・K3・E1・ROTATE)、CELL CALC ($$$ と MKC の Twiss, NX/NY)、emit (RFSW)、追跡 5 粒子 × 1/20 周 (NORAD / RAD / RAD+FLUC) | bench_level12.sad, bench_level12_track.sad / bench_level12_pysad.py, bench_level12_track_pysad.py (純 Python, 数分), bench_level12_track_fast.py (numba) | 読み込み 全要素一致 (≤5e-11)、CELL CALC ≤4e-11、emit は SAD の表示桁 (εx・α は 5 桁 → 相対 ≤1.8e-5)、追跡 ≤4.7e-11、乱数一致、track6 と track6n の差 4e-18 | Cavity.PHI の既定値は 0。parser: SEXT K2 は per-length、OCT/DECA。CALC では空洞は RFSW に依らずドリフト (4 次元 CALC)。tthin の kord≥3、tsteer。tbendi の最後の tradk は全粒子のスライス後にまとめて呼ぶ (FLUC の乱数の割り当て) | FIT/FREE/GO → 13、SADScript → 14、DA → 16、Touschek、DISP/GEO |
| 13A | マッチング FIT/FREE/GO (小リング): 13A-1 NX/NY を QF, QD で、13A-2 + ORG の BX/BY を QF, QD, QF2 で (解なし) | bench_level13a.sad / bench_level13a_pysad.py (API)、`run_sad.py bench_level13a.sad` | 反復行・Matched 行が同じ、13A-1 は残差 5.0759e-16 = SAD (解析微分; 数値微分では 5.0763e-16 = `FFS$NumericalDerivative=True` の SAD と同じ)。13A-2 は SAD も収束せず Residual 2.5977e-3 で打ち切り (pysad 同じ)。**解析微分 (既定) では 13A-2 の最終点が K1 で 3e-10 (相対)、BX@ORG で 4.7e-8 (絶対) SAD と違う** (残差の 1e-15 の大小で決まる隣の点; 数値微分にすると SAD と 10 桁一致。run_regression では KNOWN。history/work_integration.md §2) | `pysad/match6.py` = tffsmatch/tffscalc/qdcell (tdfun…)/twfit/tweigh/tfsolv/tgrad/fmincube の直訳。微分は既定で解析微分 (match6d.py。`FFS$NumericalDerivative=True` で数値微分) | FITP>1 の照合 (SAD が再起動)、INS のマッチング、R1–R4/DX の tdfun、CHI/G* |
| 13B | マッチング (ESRF v03_12): NX 47.6 NY 17.6, MKC の BX/AX/BY/AY, free Q* (7 ファミリー, 28 台) | bench_level13b.sad / bench_level13b_pysad.py、`run_sad.py bench_level13b.sad` | 2 反復 Matched、残差 3.37240e-10 = SAD (解析微分; 数値微分では 3.37243e-10)、K1 は SAD の表示 10 桁と一致、Twiss ≤4e-11 | `LINE RING=(10*(CELL CELL2))` は構文エラーなので CELLP を介す | 同上 |
| 14A | SADScript インタプリタ + FFS 層: v03_12 + esrf/calc_chromaticity.sad をそのまま実行 (Get, OpenWrite, Table, Do, DP0 スキャン + FFS["cell calc;"], PolynomialFit) | bench_level14a.sad (+ _dbg, _dbg2: 4 次元 CALC の空洞の扱いの確認) / `run_sad.py bench_level14a.sad` | 21 点のチューン (DP0 スキャン) SAD と 1e-10、PolynomialFit は SAD の式どおり (係数 相対 ≤5e-6)、EMIT の U0/α/εx/σE 表示桁で一致 | `pysad/sadscript.py` (Pratt 構文解析、autofg.f の書式、Fit.n の PolynomialFit を評価: 平均 0 のとき xl=1 になる癖も再現)、`pysad/ffscmd.py`。tffsmatch の fuzz 論理 (残差悪化時の NaN を gfortran の min/max が無視) を再現。この SAD では MON/MCOD/TCOD は何もしない | **EMIT (RADCOD ON) の dz が SAD 20.618669 mm に対し pysad 20.619466 mm (8e-7 m、z と合わせても 3e-7 m)** (run_regression では KNOWN)。esrf/matching は数値微分のため SAD と別の解 (SAD は途中で残差 3e7 を経由) |
| 15A | オフモメンタム条件 (FITP, "NX nx 2") による色収差補正 (ESRF) | bench_level15a.sad / `run_sad.py bench_level15a.sad` (→ bench_level15a_run.result) | **参照値なし**: SAD (WSL) が FITP>1 の GO で再起動する | tffscalc の idp ループ・OffMomentumWeight・INTRES 等の安定判定を直訳 (未照合) | 全部 |
| 16A | DynamicApertureSurvey[{{1,400},{1,1},{-0.5,0,0.5}}, 100, Output->6] (NORAD, RFSW, Axes "ZX") の地図・Score・戻り値 | bench_level16a.sad / `run_sad.py bench_level16a.sad` | 地図 3 行 153 文字すべて一致、Score 63 = 63。損失周数 (RES/RDAP) は損失境界の 1 粒子 (SAD 14 / pysad 15、muls=2 なので地図の文字は同じ) を run_regression で KNOWN 扱い (track6n 版では一致、docs/09 §6)。EMITS 相対 ≤1e-9 | `pysad/dapert.py` (dapert.n → tfdapert.f → trackd.f): tinip1 の初期値、tapert の損失判定 (15 要素ごと, LOSSAMPL/LOSSDZ, NaN)、rad62 の文字表示、粒子のスロット順 (i=1..n1p, j=51→1)。`pysad/track6n.py` (numba, 176 ns/粒子・要素) で数秒 (純 Python の track6 では 6.5 分) | Damping->True、DAPWIDTH の窓 |
| 16B〜E | Emittance[] の全桁比較 (16b)、一周行列と拡張 Twiss (16c)、Emittance[Matrix->True] の要素ごとの転送行列 (16d: 先頭 180 要素, 16e: 全要素) | bench_level16b〜e.sad / (調査用。比較は対話的に実施) | εx, εz, σE, σz, νx, νy, νs, α_c, バケット高さ, trf0, U0 が相対 ≤1.3e-11 | emit 経路の薄い多極 (tthine) の中のドリフトは dv 項 (tdrife0 の z −= L·dv) を含む。数値微分の刻みは 1e-7 (K1=0 BEND の (1,6) 成分の分岐境界に 1e-6 が当たるため; §3.1) | — |
| 16F〜L | RAD+FLUC の DA と乱数の順 (16f: 16a を RAD FLUC で、16g/h: 短い range、16i/j: DA の初期値 51 粒子を TrackParticles で、16k_*: 粒子数・周回数のスキャン、16l_j*: 単粒子で損失周の切り分け、_np1: NPARA=1 版) | bench_level16f〜l*.sad / bench_level16f/g/h_run.result | **損失のない場合は最終 seed まで一致** (16g_np1, 16k_n12/32_np1)。16l の 12 件 (単粒子、損失あり) は最終 seed 一致。16f (損失あり) と 20 small は地図 (Score 63) が SAD と一致するが最終 seed は不一致 (丸め誤差のカオス: SAD 自身も 1e-16 の摂動で変わる)。再同期テスト (`<work-scratch>/r2_resync.py`) では 30 周 seed 全一致 | NPARA の既定値 4 では粒子 9 個以上で fork し乱数列が分かれる (NPARA=1 にする)。落ちる粒子の乱数消費 (NaN になる途中の tdusrn が 10000 回ループ、gfortran の MAX/MIN の NaN 挙動、`rot1`/`gfloor`/`xsincos`) は SAD に合わせてある (history/work_loss_rng.md, history/work_loss_rng2.md、docs/08 §7) | 自由に走らせた最終 seed の一致 (カオスのため基準にしない) |
| 22 | INS の初期条件・軌道オフセット (kek-ls_v03_12 のリング): 22a CELL → INS で DX/DY/DPX を MARK に与える (RAD OFF)、22b CELL なしの INS (AXI/BXI、EX.., R1–R4、軌道)、22c RAD ON で CELL → INS に切り替えるだけ (DP0=0.003 も)、22d RAD ON の小さなオフセット、22e DP0=0.004・四極ずれ (閉軌道≠0)・DZ (`Print` の 11 行の最大相対差、SAD の `$FORM` 12 桁)、22f 全要素 k=1..60 の Twiss、22g 単独要素 (BEND/QUAD/SEXT/DRIFT) | bench_level22a〜g.sad / `run_sad.py bench_level22a.sad` ほか (`_pysad.result`) | 22a 2.4e-10、22b 8.7e-11、22c 1.5e-10、22d 2.9e-11、22e 1.4e-11 (許容 1e-10)、22f 3.4e-10、22g 2e-12。esrf `amp_dep_ts` (`<work-esrf>/amp_dep_ts_q.sad`) の 1 点目 NX 47.5800386898 / SAD 47.5800386899 | **CELL 後の twmov (周期解を先頭 MARK に書き戻す)、qcod の初期値は MARK の DX..DDP、INS の軌道は MARK の DX..DZ、`DP0`=`LINE["DDP",1]`**。SAD の BEND 行列は解析式で、写像の厳密な微分とずれる (px, py≠0 で 2e-8) | ずれ・ロール付き BEND の解析行列 (INS + CALC6D は 30A、ずれ・ロール付き四極・六極は 26B で照合) |
| 24a | SADScript のリスト・文字列・行列関数 (T/U/V/W、約 600 の Print) | bench_level24a.sad / `run_sad.py` + `cmp_tags.py` | 文字列一致 ($FORM="15.10"): 559 タグ一致 (SAD に無い pysad の拡張 10 タグを除く) | `pysad/sadscript.py` の組み込み関数 | — |
| 24b / 24b2 | CHROMATICITY (ESRF、960 要素) / CHROMATICITY (MULT の K1/K2、DP0≠0) | bench_level24b.sad, bench_level24b2.sad | ブロックの全行一致 | FFS コマンド CHROMATICITY | — |
| 24c | REPEAT/UNTIL、IF/ELSEIF/ELSE/ENDIF | bench_level24c.sad | 19 タグ一致 | FFS の制御構文 | — |
| 24d / 24d2 | CalculateOptics (4 次元) / CalculateOptics (ESRF、6 次元) | bench_level24d.sad, bench_level24d2.sad | 数値 (絶対 1e-8 + 相対 1e-7: R1.. の 1e-16 級の雑音を許す): 29 タグ / 6 タグ一致 | CalculateOptics | — |
| 24f | Element/SetElement/FLAG/BeamLine/WriteBeamLine/USE BeamLine | bench_level24f.sad | 164 タグ (WriteBeamLine は `L001 begin…end` のブロック)、`SetElement["ORG"]` の MARK の値 (G004) を除き一致 | Element 系の関数 | `SetElement["ORG"]` の MARK の値 (G004) |
| 26A | EMIT 報告 (temitprint) の全行: ESRF、RFSW × RADCOD の 4 通り × 37 行 | bench_level26a.sad / `cmp_emit_report.py` | 不一致行 0 (0 同然の雑音 \|値\|<1e-11 は許す) | temitprint の直訳 (docs/history/work_emit_report.md) | — |
| 26B | 四極 FRINGE=1,2,3 (tqlfre)、ずれ・ロール付き四極、六極のずれ・ロール: CALC (A/B/B2) と Emittance[] (C) | bench_level26b.sad | 1.3e-9 (CALC は 5e-13; 許容 3e-9) | 四極の解析行列 (docs/history/work_emit_report.md) | — |
| 26C | EMIT 報告の Orbit dilation (FSHIFT=1e-6, RFSW+RADCOD) | bench_level26c.sad | 不一致行 0 | 同上 | — |
| 26D | RFSW OFF + RADCOD ON の `Emittance[]` と CALC6D (ESRF) | bench_level26d.sad | 2.2e-11 (許容 1e-9) | 同上 | — |
| 27A | MAIN の `TRACK` (初期分布 r·振幅+codin、smear、tune (wtunem)): RFSW ON/RAD OFF、RAD ON/FLUC OFF、RAD ON/FLUC ON、6 粒子 250 周。SAD K64 は MAIN の TRACK が動かない (SIGSEGV) ので、参照は同じ手順を SADScript で書いた版 | `bench_level27a.sad` (SAD)、`bench_level27a_main.sad` (pysad)、`bench_level27a_pysad.py a` | 座標・tune ≤5e-11、smear・振幅 相対 ≤7e-9 | `pysad/maintrack.py` (docs/08 §8、docs/history/work_main_track.md) | FOURIE ON、`ResultsOfTracking`、TURNS<0、NPARA>1 |
| 27B | 同、粒子が落ちるケース (損失周、落ちた粒子の x,y,z、落ちた粒子を除いた smear/tune) | `bench_level27b.sad`、`bench_level27b_main.sad`、`bench_level27a_pysad.py b` | 損失 (粒子, 周) 一致、他は 27a と同じ許容 | 同上 | 同上 |
| 28C/D/F | 不安定な光学 (V: QD×1.19、H: QF×0.9) の EMIT 報告の全行と `Emittance[]` の戻り値 (28e: 両方不安定は SAD 自身が NaN で回帰に入れない) | `bench_level28c.sad` `28d` `28e` `28f` | 28c: 行単位一致 (雑音の桁のみ)、28d: ≤3e-5 (SAD の T の悪条件)、28f ≤1.6e-6 | 不安定時の temit の経路 | 28e (SAD が NaN) |
| 29A | `RADINT` コマンド (INTGRL.f): ESRF の 6 通り (CELL CALC 後・EMIT 後・NORFSW・結合 (DY/ROTATE)・DP0=0.003・四極 DX)。出力ブロック (見出し〜`LINEAR CHROMATICITY GY`) の行単位比較 | `bench_level29a.sad` / `cmp_radint.py` | 1.8e-7 (= SAD の表示 7 桁; 許容 1e-6。0 同然の V 面の量 <1e-10 は雑音) | INTGRL.f の直訳 | — |
| 29A2 | RADINT: 小リング (E1/E2 付き BEND、RF なし (`NOT ENOUGH RF VOLTS`)・あり、ロール付き四極 + DY) | `bench_level29a2.sad` | 3.2e-5 (V 面の小さい量; 許容 1e-4) | 同上 | — |
| 29B | `DISP` の位置範囲: 数値・名前・`NAME.n±k`・パターン・周回・ALL・単独・省略・未定義・R/O/E モード 36 ブロック | `bench_level29b.sad` / `cmp_blocks.py --noise 1e-9` | 全ブロック一致 | DISP の位置指定 | — |
| 29C | `CALC6D;`/`CALC4D;` 単独 (再計算されない) と、4D の CALC が縦の列 (AZ..ZPY) を直前の 6D の値のまま残す癖 | `bench_level29c.sad` | 3.4e-12 (許容 1e-9) | SAD の癖 (4D CALC の縦の列) を再現 | — |
| 29D | FFS 内の `ON/OFF flag` (無効、`Out[1]:= USE`、後ろの語がフラグ) | `bench_level29d.sad` | 全ブロック一致 (フラグの値と `Out[1]` の行) | FFS のフラグ構文 | — |
| 29E | FSHIFT=1e-6 の追跡 (ESRF、NORAD 1/10 周・RAD NOFLUC・RAD+FLUC、5 粒子) | `bench_level29e.sad` | 1.9e-11 (許容 5e-11) | dvfs の扱い (docs/08 §2.1) | — |
| 30A | INS + CALC6D (qcell6d の cell=.false. の枝): 先頭 MARK の 28 成分から 6D 伝搬 (ESRF、5 通り × 8 位置 × 28 列) | `bench_level30a.sad` | 1.6e-10 (許容 1e-9) | qcell6d の INS 枝 | — |
| 30B | RFSW ON + RADCOD OFF + FSHIFT≠0 の `Emittance[Orbit->True]` (tcod の 6D Newton の直訳) | `bench_level30b.sad` | 1.2e-11 (許容 1e-10) | tcod の直訳 | — |

Level 24〜30 のベンチのラティスは `pysad/tests/kek-ls_v03_12.lat` を先頭に連結した 1 ファイル。参照は `*.result`。

その他のファイル: `bench_moncod_test.sad` (MON/MCOD の確認: この SAD ビルドでは何もしない)、`error_one_elemet.sad` (1 台だけの誤差の確認)、
`run_cod_response.py`, `plot_optics.py` (補助)。

### 3.1 一致の程度についての注意 (CALC 経路と emit 経路)

- `_num_jacobian` (lattice.py) の数値微分の刻みは 1e-7 (K1=0 の BEND の式の分岐境界に 1e-6 が当たるため)。Emittance[] は Level 4 の α_c 8e-10、Level 1 の 6D 値 ≤3e-9。
- CALC 経路 (閉軌道上の一周行列) は、BEND (ずれなし) の行列を SAD の解析式で組み、結合の式は `qmat2twiss` の移植を使う。SAD の解析行列は写像の厳密な微分と 2e-8 ずれる (px, py≠0) ので、数値微分で組むと CALC の Twiss が相対 1e-7 程度ずれる。解析式で 3A のチューンは 1.5e-10、3B は全項目 4e-13 以内、14A の DP0 スキャンは SAD と 1e-10。
- 3C の EX/EPX は 3.2e-11 / 6e-12 で一致 (`history/work_emit_report.md`)。ずれ・ロール付き BEND の行列は数値微分のまま (`history/work_ins_offset.md`)。
- オンモメンタム・閉軌道 0 の CALC (Level 1, 2 の dp=0, 12, 13) は 1e-11 で一致している。

---

## 4. 回帰テスト `run_regression.py`

```bash
cd "<repo>"
PYTHONIOENCODING=utf-8 python -u run_regression.py            # 全部 (約 15〜20 分)
QUICK=1 PYTHONIOENCODING=utf-8 python -u run_regression.py    # Level 7, 12t (純 Python 追跡), 13a, 16a を省略 (約 5 分)
python -u run_regression.py 3b 12 13a                          # 指定したものだけ (key は表の Level 列。24a〜24f, 22e, 26a〜26d, 27a/27b, 28c/28d/28f, 29a〜29e, 30a/30b も key)
```

- 各 pysad スクリプトを子プロセスで実行し、出力の「差」「最大差」「相対差」「Matched/Residual」「Score」などの行を読んで判定する。
  `run_sad.py` で動かすもの (13a/13b/14a/16a) は SAD の `.result` の `TAG {...}` (Print の出力)、Matched/Residual 行、反復行、
  DP0 スキャンの行、PolynomialFit の行、DA の地図を突き合わせる。Level 1/2 は pysad スクリプトが SAD の値を印字しないので
  `bench_level1.result` / `bench_level2.result` から参照値を読む。Level 24 以降は `cmp_tags.py` 等 (§2.2) で SAD の出力とタグ行・ブロックを比較する。
- 許容: 追跡 1e-9、CALC 1e-9 (閉軌道・結合・オフモメンタム経由は相対 1e-7、絶対下限 1e-8)、Emittance[] 相対 1e-6
  (SAD の表示が 5 桁のとき 5e-5)、マッチング (K1・Twiss・残差) 1e-9、DA は地図と Score の一致。
  色収差の当てはめ係数 c1 は 1e-5 (チューンの 1e-8 のノイズを拡大するため)、c2/c3 は参考表示のみ。
- 判定: PASS / FAIL / KNOWN (既知の不一致: 14A の EMIT dz、16A の損失周数 1 個、13A-2 の TWB) / INFO。
  ベンチの合否は FAIL があれば FAIL、KNOWN だけなら PASS*。終了コードは FAIL があれば 1。
- 子プロセスの出力は `<work-regression>/<key>.txt` に保存する。出力形式を変えたときは `run_regression.py` の各 `check_*` を直す。

---

## 5. 回帰テストの結果

直近の実行結果は `<work-regression>/<key>.txt` (各ベンチの出力) と `run_regression.py` の標準出力で確認する。表 (§3) の「現在の結果」はこれらの値。
