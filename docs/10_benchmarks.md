# 10. ベンチマーク一覧 — SAD (K64) と pysad の照合

最終更新: 2026-11-17 (現状のスナップショット)。

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
  実ラティス (Level 12 以降の 試験用のラティス) では `<試験用ラティス>.lat` の内容を先頭に貼り、`FFS;` の後にスクリプトを書く。
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
- SAD (WSL 版) は FITP>1 の GO (bench_level15a) で `*** Run time Environment` からスクリプトを再実行してしまい
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
PYTHONIOENCODING=utf-8 python -u run_sad.py <試験用ラティス>.lat bench_level15a.sad   # ラティスとスクリプトを連結 (SAD の gs と同じ)
PYTHONIOENCODING=utf-8 python -u run_regression.py                # 回帰テスト (全部)。QUICK=1 で遅いものを省略
```

- `bench_levelN_pysad.py` は pysad の API (`Lattice`, `twiss_table`, `emit6d.emittance6d`, `track6`, `match6.Matcher` …) を直接呼び、
  SAD の値 (スクリプトに埋め込み、または `bench_levelN.result` から読む) との差を印字する。出力は `bench_levelN_pysad.result` に保存してある。
- `run_sad.py lattice script` は SADScript インタプリタ ([`pysad/sadscript.py`](../pysad/sadscript.py)) + FFS コマンド層 ([`pysad/ffscmd.py`](../pysad/ffscmd.py)) で
  SAD スクリプトをそのまま実行する。MAIN 部 (素子定義) は `parser.py`、`FFS;` 以降は上記。出力は `bench_levelN_run.result`。
- 遅いもの: Level 7 (K1≠0 BEND の純 Python 追跡、約 2 分)、`bench_level12_track_pysad.py` (純 Python、約 1 分)、
  `run_sad.py bench_level13a.sad` (13A-2 の解なしの反復が 3〜4 分)、`bench_level16a.sad` (DA、numba、数十秒〜)。
  `run_regression.py` は `QUICK=1` でこれらを飛ばす。
- 比較の補助スクリプト: `cmp_tags.py` (タグ行・ブロックの比較、Level 24)、`cmp_lines.py` (Print のラベル付き数値行を順に比較)、
  `cmp_emit_report.py` (EMIT 報告の行単位比較)、`cmp_radint.py` (RADINT 出力ブロック)、`cmp_blocks.py --noise 1e-9` (DISP などのブロック比較。
  SAD の 1e-17 などの雑音と pysad の 0 は同じ扱い)。マッチングの反復列の比較は NEWTON/Matched/Residual の行だけを抜き出して比べる。

---

## 3. 全レベルの一覧

「現在の結果」は `run_regression.py` の実行値と各ベンチの `.result` / `_pysad.result` による。
SAD 側の表示桁 (10 桁、追跡は 5e-11) が比較の下限。後ろの 2 列は、そのベンチが照合している pysad の実装と SAD の仕様 (癖を含む)、および未照合の項目。

| Level | 内容 (何を照合) | スクリプト (SAD / pysad) | 現在の結果 (一致度) | 関係する pysad の実装・SAD の仕様 | 未照合の項目 |
|---|---|---|---|---|---|
| 1 | 20 セル FODO (2.5 GeV, 140 m) の線形光学 (入口・QD.1 の β, α, η, Q)、COD (QD DX=1e-4)、RADINT (I1〜I5, εx, U0)、Emittance[] (εx, σE, α, U0) | bench_level1.sad / bench_level1_pysad.py | Twiss 表示 6 桁で一致、COD x −110.8245 µm (一致)、I1〜I5 7 桁一致、Emittance[] 相対 ≤3e-9 | periodic_optics は周期分散を含む。INTGRL.f (SYMPS, BEDGE/TEDGE) の直訳。`LINE["DX","QD.1"]` は CALC で全 QD に展開 (SAD の既定): `set_element` (ファミリー) と `set_error` (1 台) を区別。SAD の MomentumCompaction は Emittance[] (6D, −R56/C) の値で RADINT の I1/C とは別 | — |
| 2 | オフモメンタム Twiss (DP0 = 0, ±0.01, ±0.02)、K1 付き BEND (RINGK) の線形光学・RADINT・Emittance[] | bench_level2.sad / bench_level2_pysad.py | Twiss 表示 6 桁で一致 (β 相対 ≤3e-8)、RADINT 一致、Emittance[] 相対 ≤5e-9 | 分散は qcell.f の δ 列 (a13, a23) から周期解を出す式。四極のハードエッジフリンジ tqfrie、tsolque の分割 ndiv=1+int(\|L·K1L\|/0.1)、qcod の直訳 (収束判定・ダンピング・初期値は直前の閉軌道 twmov)。SAD の off-momentum β は初期値で相対 1e-6 動くので SAD と同じ順で呼ぶ | 不安定な dp (RINGK dp=−0.01) で SAD は 1/cosμ の値を返すが pysad は例外 (qcell の intres/halfres 未移植) |
| 3A | 6 極 SF/SD と色収差: DP0 スキャン 21 点のチューン、3 次当てはめ (補正前・SAD の補正後 K2) | bench_level3a.sad / bench_level3a_pysad.py | チューン 21 点 ≤1.5e-10、ξ (c1) 差 ≤1.5e-6 | 薄レンズ多極の COD 写像 tthine (L/6–2L/3–L/6)。SAD 組み込みの CHRO/RADINT の色収差は BEND 成分を含まないので使わない。SAD の PolynomialFit (Packages/Fit.n) は最高次以外の係数が壊れるので、ここでは自作の当てはめ関数で当てはめる | 色収差補正 (FIT) の K2 は 5e-4 の差 (評価関数の違い、ユーザー判断で許容: 同じ K2 を入れれば結果は一致) |
| 3B | 6 極を通る COD (DELX QD.1, DELY SF.1 を 1 台だけ)、COD 下の β/α/η/Q/EY/R1–R4、色収差 | bench_level3b.sad / bench_level3b_pysad.py | 全項目が許容内 (閉軌道・Twiss・Q・R1–R4・ξ、最大 4e-13) | tchge.f (ずれ・回転)、結合時の分散 (R1–R4 と c=√(1−det R) で対角化)。ドリフトなしで続く要素は上流と一緒に動くので QF–SF 間に DQS=0.1 m を置く | 結合時の途中要素の位相 |
| 3C | x-y 結合 (DTHETA: QF.3 ROTATE, B.5 DROTATE)、閉軌道・光学・R1–R4・εx, εy, σE, U0 | bench_level3c.sad / bench_level3c_pysad.py | 閉軌道 相対 ≤1e-8、β/α/Q ≤1e-8、EX/EPX 3.2e-11 / 6e-12、Emittance[] 相対 ≤7.8e-7 | 回転の向きは x′ = cosθ x − sinθ y に統一、四極のロールは「ずらす→回す」、BEND のロールは中点まわり (LGEO=L/2, tfgeo.f)。6D 行列の z 行はシンプレクティック条件で組む (数値微分のままだと 6.6e-5 非シンプレクティックで εy が 2% 狂う)。QD K1 は −1.7×0.3 (Qy=8.0095 だと整数共鳴に近く σE が NaN) | ずれ・ロール付き BEND の行列は数値微分のまま (`history/work_ins_offset.md`) |
| 3D | 強さの誤差 (DK QF.2, DK QD* の乱数, DDK B.3)、閉軌道・光学・εx, σE, U0 | bench_level3d.sad / bench_level3d_pysad.py | 閉軌道 相対 ≤4e-8、Twiss 相対 ≤2e-8、Emittance[] ≤3.3e-8 | DK は既定キーワード (QUAD は K1) への相対誤差、DDK は K0 への絶対誤差 (terror.f)。BEND の K0 (phib = ANGLE + K0, tbende.f) | — |
| 4 | RF 空洞 (RFSW ON): (1) RADCOD OFF、(2) RADCOD ON の Emittance[] (εx, εz, σE, σz, Qx/y/z, 減衰率, U0, Vc, trf0, α, バケット高さ, 放射込み閉軌道) | bench_level4.sad / bench_level4_pysad.py | (1) 相対 ≤2.1e-9、(2) 相対 ≤1.1e-7 (α_c 8e-10)、閉軌道 x の差 1.5e-10 m (SAD tcod の収束誤差内) | tcave.f (長さ 0 の空洞)、temit.f tecalc の synchm (縦モード)、photons.f tradke の radcod 分岐。RF なしでは z が閉じる δ を選ぶ。SAD は z のずれを trf0 側へ移すので z と trf0 の分け方が違う (z/β+trf0 は一致) | — |
| 5 | 輸送ライン (INS) とフリンジ: A) F1=0.1, FRINGE=3、B) F1=0、C) DISFRIN=1。終端と要素 6 の Twiss・軌道・位相進み | bench_level5.sad / bench_level5_pysad.py | 軌道 ≤1.3e-11、Twiss 絶対 ≤2.6e-7 (β 〜10 m)、NX/NY ≤2.7e-9 | ins_twiss (入口の Twiss・軌道から非線形写像で伝搬)。tqlfre (線形フリンジ)、tbendef1 (BEND の F1)、tsetfringepe。入口の Twiss・軌道は先頭 MARK のパラメータで与える。Do[] の反復子に文字列リストは使えない、USE の引数に変数は使えない | K1≠0 BEND の F1、逆向き要素の F1/F2 |
| 6 | 粒子追跡 TrackParticles (5 粒子 × 1/10/100 周): RAD OFF、RAD ON (FLUC OFF) | bench_level6.sad / bench_level6_pysad.py | ≤5e-11 (SAD の印字精度) | 追跡は CALC/emit と別ルーチン → [`pysad/track6.py`](../pysad/track6.py) に直訳 (tdrift, tquad/tsolqu, tthin+ttfrin, tbend/tbendr, tcav, tradk)。内部座標 px/(1+δ) と粒子ごとの dv | FLUC → Level 10、K1≠0 BEND・誤差 → Level 7、MULT/SOL → 8〜11 |
| 7 | 追跡: K1=−0.05 の複合機能 BEND のリング (7-1)、誤差 DELX QD.1 / DELY SF.1 / DTHETA QF.3, B.5 / DELX B.2 (7-2) | bench_level7.sad / bench_level7_pysad.py (約 2 分) | ≤1.7e-10 (z のみ、trf0 の表示桁の差)、他 ≤5e-11 | tbendi.f (module bendib: rhoe は入口で 1 回だけ決める)、tsolrot、TENT/TEXIT.inc、tbshift/tbrot。DUMP は 5 桁なので LINE["DX",…] で 10 桁を取る。誤差の乱数は既定の種 17 の GaussRandom の列そのもの。emit6d: 誤差で周長が変わると δ≠0 の閉軌道 | — |
| 8 | 追跡: 四極の F1=0.05/FRINGE=3 (8-1)、長さ 0.5 m の空洞 (8-2, QF K1=1.76×0.3)、6 極を MULT (K2) に (8-3) | bench_level8.sad / bench_level8_pysad.py | ≤5e-11 | tsetfringep の係数 (f1=−\|K1/L·F1²\|/24, f2=\|K1/L\|·F2)、放射ありでは tradk に f1in/f1out をそのまま渡す (SAD のまま)。tcav の L≠0 (vnominal=0 で分割 1)、tmulti (eps00m, ampmaxm)。emit6d も tcave L≠0・tmulte・tqlfre を持ち trf0 を自前で出す | MULT の ANGLE≠0 (tmulta)、SK |
| 9 | 追跡: ソレノイド領域 SOL BZ=2 T (BOUND=1, GEO=1) に D1, QD を含む。RAD OFF / ON | bench_level9.sad / bench_level9_pysad.py | ≤1.5e-10 (z, trf0 は SAD 表示値)、他 ≤5e-11 | tsfrin、tdrift_solenoid、tsolqu/tsolqur の bz 分岐、tfbzs (bzs = charge·BZ/(p0/c))。**SAD の癖: tsolqur (RAD ON, bz≠0) の最後の半ステップが最後に tzsetparams を呼んだ粒子の bzp を全粒子に使う** → 粒子の並びで結果が 1e-6 変わる。pysad は要素外側・粒子内側のループで再現 | SOL の DX/DY/DZ/CHI (trots)、F1、BOUND=0、SOL 内の BEND/MULT/CAVI。**(2026-11-16) 現在の gs (2026-10-19 ビルド) は RAD ON の SOL (BZ≠0) 内の四極で、この保存済み参照 (10-14) と別の値を返す (§6.3)**。emit6d の SOL は docs/29 |
| 10 | 追跡: RAD + FLUC (SeedRandom[17])、1/10/100 周、追跡後の乱数の状態 | bench_level10.sad / bench_level10_pysad.py | ≤4.9e-11、乱数の最終状態 3 回とも一致 | SAD の乱数 (random_plugin_sad.c: 32 bit 合同法 seed←48828125·seed、Box–Muller、GCUT) を [`pysad/sadrandom.py`](../pysad/sadrandom.py) に直訳 (15 桁一致)。FLUC の放射 (tradkfn/tradkf1, tdusrn/drndsr)。乱数を引く順: 要素外側・粒子内側、放射ありの BEND/四極/MULT はスライスごとに全粒子の本体→全粒子の tradk | — |
| 11 | 追跡: MULT の K0/K1/K2 とずれ・回転 (ANGLE=0)、閉軌道が δ=1.4e-3 にずれる場合 | bench_level11.sad / bench_level11_pysad.py | ≤5.9e-11、emit の trf0 2.6e-10 | akang (K1<0 で 90° 回転)、tbfrme、tsoldz、tsolque の dx0/dy0、ttfrin nord=2・ttfrins。**SAD の癖: 放射なしの MULT (tsolqum) の z の式は dv の項を括弧の中に入れており tsolqu と一致しない** (δ=1e-3 で 1 周 5e-10) → そのまま再現 (`_tsolqu(qum=True)`) | ANGLE≠0 の MULT、SK、MULT の F1/RF |
| 12 | 実ラティス 3GLS 型 (<試験用ラティス>.lat: 3 GeV, 603.12 m, 3563 要素): 読み込み (名前・型・L・K1・ANGLE・K2・K3・E1・ROTATE)、CELL CALC (`$$$` と MKC の Twiss, NX/NY)、emit (RFSW)、追跡 5 粒子 × 1/20 周 (NORAD / RAD / RAD+FLUC) | bench_level12.sad, bench_level12_track.sad / bench_level12_pysad.py, bench_level12_track_pysad.py (純 Python, 数分), bench_level12_track_fast.py (numba) | 読み込み 全要素一致 (≤5e-11)、CELL CALC ≤4e-11、emit は SAD の表示桁 (εx・α は 5 桁 → 相対 ≤1.8e-5)、追跡 ≤4.7e-11、乱数一致、track6 と track6n の差 4e-18 | Cavity.PHI の既定値は 0。parser: SEXT K2 は per-length、OCT/DECA。CALC では空洞は RFSW に依らずドリフト (4 次元 CALC)。tthin の kord≥3、tsteer。tbendi の最後の tradk は全粒子のスライス後にまとめて呼ぶ (FLUC の乱数の割り当て) | FIT/FREE/GO → 13、SADScript → 14、DA → 16、Touschek、DISP/GEO |
| 13A | マッチング FIT/FREE/GO (小リング): 13A-1 NX/NY を QF, QD で、13A-2 + ORG の BX/BY を QF, QD, QF2 で (解なし) | bench_level13a.sad / bench_level13a_pysad.py (API)、`run_sad.py bench_level13a.sad` | 反復行・Matched 行が同じ、13A-1 は残差 5.0759e-16 = SAD (解析微分; 数値微分では 5.0763e-16 = `FFS$NumericalDerivative=True` の SAD と同じ)。13A-2 は SAD も収束せず Residual 2.5977e-3 で打ち切り (pysad 同じ)。**解析微分 (既定) では 13A-2 の最終点が K1 で 3e-10 (相対)、BX@ORG で 4.7e-8 (絶対) SAD と違う** (残差の 1e-15 の大小で決まる隣の点; 数値微分にすると SAD と 10 桁一致。run_regression では KNOWN。history/work_integration.md §2) | [`pysad/match6.py`](../pysad/match6.py) = tffsmatch/tffscalc/qdcell (tdfun…)/twfit/tweigh/tfsolv/tgrad/fmincube の直訳。微分は既定で解析微分 (match6d.py。`FFS$NumericalDerivative=True` で数値微分) | FITP>1 の照合 (SAD が再起動)、INS のマッチング、R1–R4/DX の tdfun、CHI/G* |
| 13B | マッチング (試験用のラティス): NX 47.6 NY 17.6, MKC の BX/AX/BY/AY, free Q* (7 ファミリー, 28 台) | bench_level13b.sad / bench_level13b_pysad.py、`run_sad.py bench_level13b.sad` | 2 反復 Matched、残差 3.37240e-10 = SAD (解析微分; 数値微分では 3.37243e-10)、K1 は SAD の表示 10 桁と一致、Twiss ≤4e-11 | `LINE RING=(10*(CELL CELL2))` は構文エラーなので CELLP を介す | 同上 |
| 14A | SADScript インタプリタ + FFS 層: 試験用のラティス + 色収差スキャンの台本をそのまま実行 (Get, OpenWrite, Table, Do, DP0 スキャン + FFS["cell calc;"], PolynomialFit) | bench_level14a.sad (+ _dbg, _dbg2: 4 次元 CALC の空洞の扱いの確認) / `run_sad.py bench_level14a.sad` | 21 点のチューン (DP0 スキャン) SAD と 1e-10、PolynomialFit は SAD の式どおり (係数 相対 ≤5e-6)、EMIT の U0/α/εx/σE 表示桁で一致 | [`pysad/sadscript.py`](../pysad/sadscript.py) (Pratt 構文解析、autofg.f の書式、Fit.n の PolynomialFit を評価: 平均 0 のとき xl=1 になる癖も再現)、[`pysad/ffscmd.py`](../pysad/ffscmd.py)。tffsmatch の fuzz 論理 (残差悪化時の NaN を gfortran の min/max が無視) を再現。この SAD では MON/MCOD/TCOD は何もしない | **EMIT (RADCOD ON) の dz が SAD 20.618669 mm に対し pysad 20.619466 mm (8e-7 m、z と合わせても 3e-7 m)** (run_regression では KNOWN) |
| 15A | オフモメンタム条件 (FITP, "NX nx 2") による色収差補正 (実ラティス) | bench_level15a.sad / `run_sad.py bench_level15a.sad` (→ bench_level15a_run.result) | **参照値なし**: SAD (WSL) が FITP>1 の GO で再起動する | tffscalc の idp ループ・OffMomentumWeight・INTRES 等の安定判定を直訳 (未照合) | 全部 |
| 16A | DynamicApertureSurvey[{{1,400},{1,1},{-0.5,0,0.5}}, 100, Output->6] (NORAD, RFSW, Axes "ZX") の地図・Score・戻り値 | bench_level16a.sad / `run_sad.py bench_level16a.sad` | 地図 3 行 153 文字すべて一致、Score 63 = 63。損失周数 (RES/RDAP) は損失境界の 1 粒子 (SAD 14 / pysad 15、muls=2 なので地図の文字は同じ) を run_regression で KNOWN 扱い (track6n 版では一致、docs/09 §6)。EMITS 相対 ≤1e-9 | [`pysad/dapert.py`](../pysad/dapert.py) (dapert.n → tfdapert.f → trackd.f): tinip1 の初期値、tapert の損失判定 (15 要素ごと, LOSSAMPL/LOSSDZ, NaN)、rad62 の文字表示、粒子のスロット順 (i=1..n1p, j=51→1)。[`pysad/track6n.py`](../pysad/track6n.py) (numba, 176 ns/粒子・要素) で数秒 (純 Python の track6 では 6.5 分) | Damping->True、DAPWIDTH の窓 |
| 16B〜E | Emittance[] の全桁比較 (16b)、一周行列と拡張 Twiss (16c)、Emittance[Matrix->True] の要素ごとの転送行列 (16d: 先頭 180 要素, 16e: 全要素) | bench_level16b〜e.sad / (調査用。比較は対話的に実施) | εx, εz, σE, σz, νx, νy, νs, α_c, バケット高さ, trf0, U0 が相対 ≤1.3e-11 | emit 経路の薄い多極 (tthine) の中のドリフトは dv 項 (tdrife0 の z −= L·dv) を含む。数値微分の刻みは 1e-7 (K1=0 BEND の (1,6) 成分の分岐境界に 1e-6 が当たるため; §3.1) | — |
| 16F〜L | RAD+FLUC の DA と乱数の順 (16f: 16a を RAD FLUC で、16g/h: 短い range、16i/j: DA の初期値 51 粒子を TrackParticles で、16k_*: 粒子数・周回数のスキャン、16l_j*: 単粒子で損失周の切り分け、_np1: NPARA=1 版) | bench_level16f〜l*.sad / bench_level16f/g/h_run.result | **損失のない場合は最終 seed まで一致** (16g_np1, 16k_n12/32_np1)。16l の 12 件 (単粒子、損失あり) は最終 seed 一致。16f (損失あり) と 20 small は地図 (Score 63) が SAD と一致するが最終 seed は不一致 (丸め誤差のカオス: SAD 自身も 1e-16 の摂動で変わる)。再同期テストでは 30 周 seed 全一致 | NPARA の既定値 4 では粒子 9 個以上で fork し乱数列が分かれる (NPARA=1 にする)。落ちる粒子の乱数消費 (NaN になる途中の tdusrn が 10000 回ループ、gfortran の MAX/MIN の NaN 挙動、`rot1`/`gfloor`/`xsincos`) は SAD に合わせてある (history/work_loss_rng.md, history/work_loss_rng2.md、docs/08 §7) | 自由に走らせた最終 seed の一致 (カオスのため基準にしない) |
| 22 | INS の初期条件・軌道オフセット (試験用のラティス のリング): 22a CELL → INS で DX/DY/DPX を MARK に与える (RAD OFF)、22b CELL なしの INS (AXI/BXI、EX.., R1–R4、軌道)、22c RAD ON で CELL → INS に切り替えるだけ (DP0=0.003 も)、22d RAD ON の小さなオフセット、22e DP0=0.004・四極ずれ (閉軌道≠0)・DZ (`Print` の 11 行の最大相対差、SAD の `$FORM` 12 桁)、22f 全要素 k=1..60 の Twiss、22g 単独要素 (BEND/QUAD/SEXT/DRIFT) | bench_level22a〜g.sad / `run_sad.py bench_level22a.sad` ほか (`_pysad.result`) | 22a 2.4e-10、22b 8.7e-11、22c 1.5e-10、22d 2.9e-11、22e 1.4e-11 (許容 1e-10)、22f 3.4e-10、22g 2e-12 | **CELL 後の twmov (周期解を先頭 MARK に書き戻す)、qcod の初期値は MARK の DX..DDP、INS の軌道は MARK の DX..DZ、`DP0`=`LINE["DDP",1]`**。SAD の BEND 行列は解析式で、写像の厳密な微分とずれる (px, py≠0 で 2e-8) | ずれ・ロール付き BEND の解析行列 (INS + CALC6D は 30A、ずれ・ロール付き四極・6 極は 26B で照合) |
| 24a | SADScript のリスト・文字列・行列関数 (T/U/V/W、約 600 の Print) | bench_level24a.sad / `run_sad.py` + `cmp_tags.py` | 文字列一致 (`$FORM="15.10"`): 559 タグ一致 (SAD に無い pysad の拡張 10 タグを除く) | [`pysad/sadscript.py`](../pysad/sadscript.py) の組み込み関数 | — |
| 24b / 24b2 | CHROMATICITY (実ラティス、960 要素) / CHROMATICITY (MULT の K1/K2、DP0≠0) | bench_level24b.sad, bench_level24b2.sad | ブロックの全行一致 | FFS コマンド CHROMATICITY | — |
| 24c | REPEAT/UNTIL、IF/ELSEIF/ELSE/ENDIF | bench_level24c.sad | 19 タグ一致 | FFS の制御構文 | — |
| 24d / 24d2 | CalculateOptics (4 次元) / CalculateOptics (実ラティス、6 次元) | bench_level24d.sad, bench_level24d2.sad | 数値 (絶対 1e-8 + 相対 1e-7: R1.. の 1e-16 級の雑音を許す): 29 タグ / 6 タグ一致 | CalculateOptics | — |
| 24f | Element/SetElement/FLAG/BeamLine/WriteBeamLine/USE BeamLine | bench_level24f.sad | 164 タグ (WriteBeamLine は `L001 begin…end` のブロック)、`SetElement["ORG"]` の MARK の値 (G004) を除き一致 | Element 系の関数 | `SetElement["ORG"]` の MARK の値 (G004) |
| 26A | EMIT 報告 (temitprint) の全行: 実ラティス、RFSW × RADCOD の 4 通り × 37 行 | bench_level26a.sad / `cmp_emit_report.py` | 不一致行 0 (0 同然の雑音 \|値\|<1e-11 は許す) | temitprint の直訳 (docs/history/work_emit_report.md) | — |
| 26B | 四極 FRINGE=1,2,3 (tqlfre)、ずれ・ロール付き四極、6 極のずれ・ロール: CALC (A/B/B2) と Emittance[] (C) | bench_level26b.sad | 1.3e-9 (CALC は 5e-13; 許容 3e-9) | 四極の解析行列 (docs/history/work_emit_report.md) | — |
| 26C | EMIT 報告の Orbit dilation (FSHIFT=1e-6, RFSW+RADCOD) | bench_level26c.sad | 不一致行 0 | 同上 | — |
| 26D | RFSW OFF + RADCOD ON の `Emittance[]` と CALC6D (実ラティス) | bench_level26d.sad | 2.2e-11 (許容 1e-9) | 同上 | — |
| 27A | MAIN の `TRACK` (初期分布 r·振幅+codin、smear、tune (wtunem)): RFSW ON/RAD OFF、RAD ON/FLUC OFF、RAD ON/FLUC ON、6 粒子 250 周。SAD K64 は MAIN の TRACK が動かない (SIGSEGV) ので、参照は同じ手順を SADScript で書いた版 | `bench_level27a.sad` (SAD)、`bench_level27a_main.sad` (pysad)、`bench_level27a_pysad.py a` | 座標・tune ≤5e-11、smear・振幅 相対 ≤7e-9 | [`pysad/maintrack.py`](../pysad/maintrack.py) (docs/08 §8、docs/history/work_main_track.md) | FOURIE ON、`ResultsOfTracking`、TURNS<0、NPARA>1 |
| 27B | 同、粒子が落ちるケース (損失周、落ちた粒子の x,y,z、落ちた粒子を除いた smear/tune) | `bench_level27b.sad`、`bench_level27b_main.sad`、`bench_level27a_pysad.py b` | 損失 (粒子, 周) 一致、他は 27a と同じ許容 | 同上 | 同上 |
| 28A / 28B | 不安定ライン (INTRES/HALFRES): CELL CALC の続行・Unstable 表示・GO (nstab)。28A は V・H・HV の 3 通り、28B は INTRES/HALFRES の 4 通り (NOINTRES では安定扱いで tracex=4/trace) | `bench_level28a.sad`, `bench_level28b.sad` / `run_sad.py` + `cmp_lines.py` | Print の行ごとに相対 1e-8 (H 面の TWF 1e-3、K1F 1e-6 は緩め。28A の TWU_H は 1e-8) | 不安定時の tffscalc の経路 (docs/history/work_unstable.md) | — |
| 28C/D/F | 不安定な光学 (V: QD×1.19、H: QF×0.9) の EMIT 報告の全行と `Emittance[]` の戻り値 (28e: 両方不安定は SAD 自身が NaN で回帰に入れない) | `bench_level28c.sad` `28d` `28e` `28f` | 28c: 行単位一致 (雑音の桁のみ)、28d: ≤3e-5 (SAD の T の悪条件)、28f ≤1.6e-6 | 不安定時の temit の経路 | 28e (SAD が NaN) |
| 29A | `RADINT` コマンド (INTGRL.f): 実ラティスの 6 通り (CELL CALC 後・EMIT 後・NORFSW・結合 (DY/ROTATE)・DP0=0.003・四極 DX)。出力ブロック (見出し〜`LINEAR CHROMATICITY GY`) の行単位比較 | `bench_level29a.sad` / `cmp_radint.py` | 1.8e-7 (= SAD の表示 7 桁; 許容 1e-6。0 同然の V 面の量 <1e-10 は雑音) | INTGRL.f の直訳 | — |
| 29A2 | RADINT: 小リング (E1/E2 付き BEND、RF なし (`NOT ENOUGH RF VOLTS`)・あり、ロール付き四極 + DY) | `bench_level29a2.sad` | 3.2e-5 (V 面の小さい量; 許容 1e-4) | 同上 | — |
| 29B | `DISP` の位置範囲: 数値・名前・`NAME.n±k`・パターン・周回・ALL・単独・省略・未定義・R/O/E モード 36 ブロック | `bench_level29b.sad` / `cmp_blocks.py --noise 1e-9` | 全ブロック一致 | DISP の位置指定 | — |
| 29C | `CALC6D;`/`CALC4D;` 単独 (再計算されない) と、4D の CALC が縦の列 (AZ..ZPY) を直前の 6D の値のまま残す癖 | `bench_level29c.sad` | 3.4e-12 (許容 1e-9) | SAD の癖 (4D CALC の縦の列) を再現 | — |
| 29D | FFS 内の `ON/OFF flag` (無効、`Out[1]:= USE`、後ろの語がフラグ) | `bench_level29d.sad` | 全ブロック一致 (フラグの値と `Out[1]` の行) | FFS のフラグ構文 | — |
| 29E | FSHIFT=1e-6 の追跡 (実ラティス、NORAD 1/10 周・RAD NOFLUC・RAD+FLUC、5 粒子) | `bench_level29e.sad` | 1.9e-11 (許容 5e-11) | dvfs の扱い (docs/08 §2.1) | — |
| 30A | INS + CALC6D (qcell6d の cell=.false. の枝): 先頭 MARK の 28 成分から 6D 伝搬 (実ラティス、5 通り × 8 位置 × 28 列) | `bench_level30a.sad` | 1.6e-10 (許容 1e-9) | qcell6d の INS 枝 | — |
| 30B | RFSW ON + RADCOD OFF + FSHIFT≠0 の `Emittance[Orbit->True]` (tcod の 6D Newton の直訳) | `bench_level30b.sad` | 1.2e-11 (許容 1e-10) | tcod の直訳 | — |
| 31A | 誤差コマンド (DELX/DELY/DK/DTHETA/DDK, UNIFORM/COHERENT/ADD): master グループ共有・tfadjst・触ったキーの展開と CALC (実ラティス、128 タグ) | `bench_level31a.sad` / `cmp_tags.py` | 絶対 1e-10 + 相対 1e-7 で一致 | terror.f と tfadjst の直訳 | — |
| 31B | `Twiss[]` / `LINE[]` の位置パターンの返り値 (一致 1 個ならスカラー、0 個は `{}`) (実ラティス) | `bench_level31b.sad` / `cmp_tags.py` | 絶対 1e-10 + 相対 1e-8 で一致 | 位置パターンの返り値の形 | — |

Level 24〜31 のベンチのラティスは `<試験用ラティス>.lat` を先頭に連結した 1 ファイル。参照は `*.result`。

その他のファイル: `bench_moncod_test.sad` (MON/MCOD の確認: この SAD ビルドでは何もしない)、
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

### 5.1 2026-11-17 時点の現状

- 回帰 (`PYTHONIOENCODING=utf-8 py -3.12 -u run_regression.py`、QUICK なし): **50 ベンチをすべて実行して FAIL 0** (2026-11-16 夜、205 s。11/17 の追跡の更新 (APERT) のあとは QUICK=1 の 46 件と、SKIP される 4 件の単独実行で FAIL 0)。`run_regression.py` の `BENCHES` は 50 キー (1, 2, 3a〜3d, 4〜11, 12/12t/12f, 13a, 13b, 14a, 16a, 22e, 24a〜24f, 26a〜26d, 27a/27b, 28a〜28d/28f, 29a〜29e/29a2, 30a/30b, 31a/31b)。
  `QUICK=1` は遅い 4 件 (7, 12t, 13a, 16a) を SKIP して 46 件を実行する。その 4 件も 2026-11-17 に単独で実行して**すべて PASS** (最大差: 7 が 5.0e-11、12t が 3.5e-11、13a が 4.3e-21、16a が 1.4e-31)。QUICK の 46 件と合わせて 50 ベンチ全部が FAIL 0。
  `run_regression.py` 自体は 2026-11-04 から変更なし。KNOWN (PASS*) は §4 の 3 件 (14A の EMIT dz、16A の損失周数 1 個、13A-2 の TWB)。
- pytest (`py -3.12 -m pytest pysad/tests -q --ignore=pysad/tests/test_bend_feeddown.py --ignore=pysad/tests/test_coord.py --ignore=pysad/tests/test_sext_radiation.py`): **937 通過・4 skip** (941 件収集、2026-11-17。k64 モードが既定)。
  - GUI 系 ([`test_plot.py`](../pysad/tests/test_plot.py), [`test_tclbridge.py`](../pysad/tests/test_tclbridge.py)) は全体実行で数回に 1 回失敗することがある (単独なら通る)。環境変数 `PYSAD_SAD_PACKAGES=<SAD の Packages のフォルダ>` がないと GUI 系は skip する。
  - `PYSAD_SLOW=1` で遅い試験 (CSR の MakeZL、約 7 分) も走る。通常の 937 通過にはこれは含まれない。
  - SAD 本体 (gs) が無くても通る (gs の値は試験の中か [`pysad/tests/`](../pysad/tests/) のデータに埋め込み)。

## 6. 新しい gs 照合 (2026-11-13〜11-16)

SAD のレベル別ベンチ (§3) とは別に、要素・機能ごとに gs (WSL の K64 ビルド) と pysad を同じ入力で比べ、試験 (`pysad/tests/test_*.py`) に固定した。
詳しい比較と限界は各文書:

| 項目 | 一致の程度 (gs との差) | 試験 | 文書 |
|---|---|---|---|
| SOL の Twiss (CALC / INS) | チューン 1e-15、領域内の β・α・分散・結合 1e-11 | test_sol_optics.py | docs/29 |
| SOL の Emittance[]・EMIT (12 セルのリング、RFSW) | ε_x 2e-12〜4e-9、ε_y 4e-10〜2e-8、チューン 1e-11、軌道 1e-9 | test_sol_emit.py | docs/29 |
| POL: Emittance[]・EMIT のスピン (平面・ずれ・MULT・空洞・SOL まで) | 偏極時間・スピンチューン・偏極ベクトルが一致 (MULT のリングは 1e-11) | test_spin.py | docs/36 |
| POL: TrackParticles のスピン追跡 (RAD ON のときだけ) | 4e-12 以下 (NOFLUC・FLUC)。損失の処理なし | test_spin_track.py | docs/36 |
| RADLIGHT (TrackParticles の軌道記録) と RadiationField | 差 0 (2 粒子・複数周・始点 ls≠1) | test_radlight.py, test_radfield.py | docs/36, docs/history/work_unimpl.md |
| CSR (CSRTrack・CSRHaissin・CSRSetupOY・CSRMatrix・CSRConvert ほか) | 11 桁〜1e-9 (関数ごと)。MakeZL は 1e-10 | test_csr.py (MakeZL は PYSAD_SLOW=1) | docs/36 |
| SynchroBetaEmittance | 12 セルのリング + CAVI で emitx が相対 4e-7〜1e-5 (nus=0.005〜0.02, モード 5〜9)。24 セルのリング (モード 7) で 1.4e-7〜1.4e-6。emity は 1e-18 の雑音で符号も違う。結合のあるラティスは gs が SIGSEGV | test_synchrobeta.py (許容 3e-5) | docs/history/work_unimpl.md |
| COORD / TCAVI / UND / MULT ANGLE≠0 | 追跡 (RAD FLUC まで)・Twiss・Emittance[] が一致 (TCAVI の Emittance[] は L=0 のみ、UND は K64 が磁場を設定できない) | test_coord_elem.py ほか | docs/25〜26 |
| MAP (ExternalMap)・DELCOR | 追跡 10 桁、Emittance[] 1e-11 / チューン 1e-14・軌道 1e-10 | test_map.py, test_delcor.py | docs/30, 29 |
| GaussianCoulomb(Fitted)・BBBrem・特殊関数・Fit・FindRoot | 458 式の特殊関数が一致、Fit は 5e-10〜1e-9 | test_gausscoul.py, test_bbbrem.py, test_specfn.py, test_fitroot.py | docs/32, 31 |
| 幾何 fit の tdgeo (DRIFT の L、BEND の ANGLE・L) | 幾何の有限差分と 1e-5 以内 (他の変数は中心差分のまま) | test_tdgeo.py | docs/history/work_unimpl.md |
| Skip・Read の型 | gs に合わせた | test_evalfixes.py | docs/34 |
| FIT のフラグ UNIPOL・RELW | gs と全桁一致 (QF K1・QD K1 の 2 変数で NX・NY=0.2。3 通り) | test_item4.py (3 件) | docs/06, 03 |
| DUMP の書式 (DK・DDK の `1p g12.5`) | 値は一致、書式を gs に合わせた | test_ffsextra.py | docs/03 |
| LDSolve・複素行列の Dot・Inverse | 終点 (cos 1, −sin 1) と 3 通りの呼び出しが gs と一致 | test_evalfixes.py (2 件) | docs/04, 05 |
| セッション最初の TrackParticles の `l_track=0` | gs と全桁一致 (15 値、1 回目・2 回目の呼び出し) | test_small_items.py | docs/08 |
| POL のスピン: ROTATE+DX+DY 付き空洞、MULT ANGLE≠0、BEND の DROTATE+DX | 偏極時間 13〜14 桁、スピンチューン 1e-14、偏極ベクトル 1e-13〜1e-14 | (試験なし。手動比較) | docs/07 |
| Twiss・LINE の端数位置 (3.5・4.25・4.5・リスト) | 1e-14 (平面の 12 セルのリング) | (試験なし) | docs/05 |
| ListPlot・Plot のオプション 10 例 | 戻りの Graphics 式が文字列で完全一致 (見た目ではない。`Options[g]` は差あり) | (試験なし) | docs/37 |

### 6.1 限界・未確認 (上の表の裏側)

- 上の一致はどれも、特定のエネルギー・リング (12 セルのリング、小リング、603 m の 3GLS 型 など。試験のラティスごとに異なる) と、試験に固定した入力での値で、一般の入力での一致ではない。
- SynchroBetaEmittance: 2026-11-16 に 24 セルの平面リング + CAVI (AzimuthalModes->7) の 4 点で emitx が gs と相対 1.4e-7〜1.4e-6 で一致。結合のあるラティスは gs が SIGSEGV で落ちて比較できず、pysad も NaN→整数の ValueError で失敗する。不安定なチューン、`emiout`、RFSW の組み合わせは未確認。条件数 3e6 の連立方程式の解き方 (K64 の tsvd と numpy の SVD) と閉軌道の収束履歴が差の原因と考えているが、確かめていない。
- POL のスピン追跡は、SOL (BZ≠0) 内の四極のスピンを gs と照合できていない (§6.3)。SOL 領域内の MULT・BEND・CAVI は未対応。`Polarization` は K64 が NaN を返す (B55)。
- CSR の ScanTrack・OYStability などの長い流れは未照合。`CSR.n` の `TrackAnimate` は `System["import …"]` (ImageMagick) と `System["mkdir …"]` を呼ぶため Windows の pysad では画像が作れない (2026-11-16 にソースで確認)。`Wait[mem]` の `FreeMem[]` は gs がこの環境で top が無く OpenRead エラー、pysad は psutil で 1.0 を返す。並列を増やすと速くなるかは未確認 (pysad は fork しない)。
- RADLIGHT → RadiationField → FindPeak → RadiationSpectrum の流れ (B を ANGLE=0 K0=0.1 の STEER にして 1 周) は 2026-11-16 に pysad で最後まで通ることを確認 (軌道 8 点 {4,8}、f {10,8}、ピーク位置 3、スペクトル {3,1001}) が、gs は RadiationField で SIGSEGV (B49) するので通しの比較はできない。TrackParticles の戻り (軌道) の部分は gs と一致。**複数区間 (複数周) の軌道は区間の始めの点が重なり (t が同じ点が 2 つ)、pysad の RadiationField は `x must be strictly increasing sequence` で失敗する** (gs は落ちるので比較不能)。

### 6.2 見つけた未解決の差 (いずれも `docs/history/work_unimpl.md`)

1. BEND の `F1=0.1 FRINGE=3` を RAD ON で追跡すると gs と合わない (相対 0.7)。DX・DY・ROTATE だけなら 1e-13 で一致。
2. L≠0 の CAVI に DX・DY を付けると 1e-3 ずれる。
3. K3 が強くずれの大きい MULT リング (K3=300、ずれ 1e-4) で、pysad の CALC (qcod) の閉軌道が gs と違う解 (y=-0.033) になり Emittance が NaN。gs は小さい軌道で通る。qcod の初期値・後退の違いの可能性 (未調査。2026-11-16 夜の確認では取り上げていない)。
4. SOL 領域内の BEND・MULT・CAVI は未対応 (2026-11-16 夜の確認では取り上げていない)。BEND は K64 の tsol.f どおりに実装して試したが、gs が偏向をほぼ無視して合わず取り下げた (原因未確認)。
5. `emit6d._tquade` (Python 版) の解析ヤコビアンが、x がずれた軌道で δ の微分を落とし、`emit_extra` の Orbit・Matrix 規則・`calc6d`・`touschek` で T の δ 列が 1e-5 (ε_y 1e-4) ずれる。`emit6d.tturne` (elkern) は gs と 1e-12。

### 6.3 保存済みの参照 (bench_level9) と現在の gs の食い違い

- `bench_level9.sad` / `bench_level9.result` (Level 9、SOL BZ=2 T の追跡) は、2026-10-14 に取った gs の出力を保存した参照で、pysad はこれと一致する (回帰 9 は PASS)。
- ところが現在の gs (2026-10-19 ビルド) に同じ入力を通すと、**RAD ON の SOL (BZ≠0) 内の四極**の追跡で別の値を返す (例: `S R T1 P1` の x が保存済みの -0.0027867898 に対し現在の gs では -0.0014539780)。NORAD は両者で一致する。
- **2026-11-16 夜の切り分け** (単一粒子の 1 周追跡、2.5 GeV、RAD NOFLUC、emit のあとで SOL の BZ を 0 → 0.01 T に変えた。粒子の並びの癖は効かない): (1) 領域内がドリフトだけ: gs = pysad (BZ=0・0.01 とも 9 桁以上)。(2) 領域内が四極: BZ=0 は gs = pysad、BZ=0.01 で gs は x=6.33e-4・y=−5.6e-16 (横の結合なし)、pysad は x=4.04e-4・y=1.1e-8。(3) `SB D1 QD SE` BZ=0.01 の gs は、四極なしの `SB D1 SE` と 8 桁で一致する (x=1.194038e-4、px=−5.363477e-4)。つまり **gs は RAD ON かつ BZ≠0 の領域内で四極の収束・発散力を入れていない**ように振る舞い、BZ=0.01 と 0.2 でも x がほぼ同じ (2.95e-5) で BZ に連続でない。(4) 現在の gs の 2 つのビルド (10-19 の Linux6 と、11-14 に -O0 で作った `~/sad_dbg`) は 10 桁一致 (x=2.954816161E-5) で最適化の違いではない。保存済みの参照 (10-14) を出した 10-15 の古い実行ファイルは symbol lookup error で動かず、同じ版か確認できない。**結論: pysad が壊れているのではなく、現在の gs (K64 のソース) が SOL 領域内の四極 + RAD を物理的に誤って扱う疑いが強い** (pysad は物理と保存済みの参照に合う)。ソースのどこかは未特定。
- 上のとおり断定はできない (旧版の実行ファイルで確認できない)。pysad は保存済みの参照 (K64 の tsolqur の癖「最後のステップが最後に tzsetparams を呼んだ粒子の bzp を全粒子に使う」を再現したもの) に合わせてある。現在の gs でこの癖が変わったのか、gs 側の別の変更かは調べていない。
- 回帰 9 は保存済みの参照と比べているので、現在の gs で参照を取り直すと FAIL する可能性がある (取り直しは未実施)。SOL 内の四極 + RAD の追跡を使う場合は、gs のビルドの版を明記すること。
