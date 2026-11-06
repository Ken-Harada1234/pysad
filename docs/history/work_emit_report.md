# 作業ログ: EMIT 報告の全行 (temitprint) と四極・薄い多極の解析行列 (2026-10-16)

担当: `pysad/ffscmd.py` (print_emit)、`pysad/emit6d.py`、`pysad/lattice.py`。課題 A (EMIT 報告の欠けている行)、課題 B (ずれ・ロール付き四極と薄い多極、線形フリンジ付き四極の解析行列)。
SAD の参照は WSL の `gs` (K64 fork) で、ラティスを連結した 1 ファイル (`NPARA=1`) を実行して取った。

## 0. 結果の要約

| 項目 | 作業前 | 作業後 |
|---|---|---|
| EMIT 報告 (ESRF、RFSW × RADCOD の 4 通り、各 37 行) | Extended Twiss / dl / h / fs / Imag.tune / Damping time / Tune shift / 分配数 / spin の行が無い。`Real tune` の空白が 1 桁違う。RFSW OFF + RADCOD ON は SAD と値が違う (U0 7e-5、tune 1e-5) | **4 通り 148 行がすべて SAD と一致** (0 同然の雑音の桁 |値|<1e-11 の差 38 行だけ別。下記) |
| 22e CELL (DP0=0.004、四極ずれ、閉軌道≠0) | 最大 5.3e-10 | **1.4e-11** (目標 1e-10) |
| 3C の EX/EPX (結合時の分散) | 既知の不一致として 3.2e-7 / 7e-8 としていたが、再測定すると 3.2e-11 / 6e-12 で、**この作業の前から解消済み** (10/15 の qmat2twiss と z 行の組み直しの効果。`bench_level3c_pysad.result` が古かった) | 3.2e-11 / 6e-12 (run_regression の既知リストから外した) |
| 四極 FRINGE=1〜3 + ずれ・ロール + 六極のずれ・ロール (新 bench 26b) の CALC/Emittance[] | `Emittance[]` (C) の最大相対差 1.1e-5 (チューン。ε は 2.7e-7、数値微分と tqlfre の piece が数値微分) | CALC は 1e-12 級、Emittance[] は ≤1.3e-9 (SAD の tcod 収束の差) |
| 16b 全桁 | ≤1.3e-11 | 変化なし (出力が作業前と同一) |

## 1. 課題 A: EMIT 報告 (temitprint) の書式と式

SAD の `temit.f` の `temitprint` (949-1250 行) を読み、`FFSState.print_emit` を次の順に直した。値は `emit6d.tecalc` が返す辞書に新しいキー (`cd`, `omegaz`, `heff`, `omega0`, `h0`, `ext_twiss`) を足して供給する。

| 行 (SAD の format) | 式 (SAD) | pysad |
|---|---|---|
| **Extended Twiss Parameters** (9001) | `params(iptwiss..) = tfetwiss(ri, cod, .true.)` を autofg `9.6` (幅 9) で、AX BX ZX EX / PSIX ZPX EPX / R1 R2 AY BY ZY EY / R3 R4 PSIY ZPY EPY / AZ BZ / PSIZ の 6 行 + Units 行 + 空行 | `emit6d._ext_twiss` (= `calc6d.tfetwiss(ri, 最後の一周の出口 cod)`) を `ffscmd._print_emit_exttwiss` が書く。**SAD の癖: EPX の欄に EY を渡している (引数の並びのミス)** ので EPX の値は EY と同じ。docs/11 A30 |
| Orbit dilation dl (9101) | `-dleng·1e3` [mm]、`dleng = FSHIFT·周長` (tfgeo.f:134) | グローバル `FSHIFT` × 周長。下の FSHIFT 参照 |
| Effective harmonic # h | `heff = wrfeff·c/ω0` (Vc≠0 のとき、なければ 0) | tecalc の `heff` |
| Synchrotron frequency fs | `omegaz/2π`。RFSW ON: `|Im cd(6)|·ω0/2π`、RFSW OFF: `√|α_c 2π heff Vc/pgev cos φrf| ω0/2π` | tecalc の `omegaz` |
| Imag.tune / Real tune (9012: `3(f10.7,10x)`) | `Re cd(4:6)/2π` (= log|λ|/2π、安定なら 0 の雑音) / `Im cd(4:6)/2π`。行末の `10x` は出力されない | `cd[3:6]` |
| Damping per one revolution (9013: `10x,3(a,1p g14.6,3x)`) | `Re cd(1:3)` | `fmt_g(v, 14, 6)` (gfortran の G14.6: 0.1≤|x|<1e6 は F 形式、それ以外は 7 桁の E 形式)。以前は E 形式固定 |
| **Damping time (sec)** | `-2π/ω0/Re cd(i)` | 同じ式 |
| **Tune shift due to radiation** | `Im cd(1:3)/2π` | 同じ式 |
| **Damping partition number** (9014: `3(a,f10.4,7x)`) | `J_i = Re cd(i)·4/Σ Re cd` (tecalc: `sr = 4/Σdamp`) | 同じ式 |
| Emittance / Bunch / Beam tilt / Beam size | 既存 | 変更なし |
| **Nominal spin tune, Polarization time** (9102 の最後の行) | `params(ipnnup) = h0·GSPIN` (11.7)、`params(iptaup)/60` (POL OFF は 0 = `.0000000`) | `h0·GSPIN`、`0.0`。偏極 (calpol) の行は POL ON のときだけなので出さない |

- 行末の `10x` `3x` `7x` (X 編集記述子) は出力されない。SAD の出力に末尾空白が付かないのでそれに合わせた。
- NaN の報告 (閉軌道が見つからないとき): SAD は NaN を幅に右詰めで出す (`       NaN`)。`_fr`/`_ff`/`_f107`/`_line9013` で同様にした。Extended Twiss は ri=単位行列の `tfetwiss` (BX=BY=BZ=1、他 0)、`h=0`、`Beam tilt=.00000000` も SAD の報告に合わせた (Closed orbit の Entrance の値は SAD では未初期化のゴミなので再現しない)。

### RFSW OFF + RADCOD ON は別の閉軌道 (docs/11 A32)

4 通りのうち **RFSW OFF、RADCOD ON** (R3) だけ値が違っていた。`closed_orbit6d` は RFSW OFF を RADCOD の指定に関係なく RADCOD OFF として扱っていた。
SAD の `tcod` は `rfsw=.false.` の枝 (im=5: x,px,y,py だけを閉じる) を `radcod` で実行し、**一周の放射損失を含めた閉軌道** (出口の δ は −U0/P0 = −9.87e-5)、z のずれは trf0 (Equilibrium dz) に移す。
`_tcod_radcod(rfsw=False)` (これは RFSW+RADCOD の前半と同じ関数) を RFSW OFF + RADCOD ON でも呼び、`emittance6d` の最後の一周も `radcod=radcod` で追うようにした。
ESRF の R3: U0 .2960539 MV (旧 .2960749)、dz .0034267 mm (旧 0)、α_c、Q_z の tune shift、減衰率、ε、Bucket height が SAD に一致。**CALC6D と `Emittance[]` (T1〜T7) も一致** (`bench_level26d`: 最大 2.2e-11)。

### FSHIFT (Orbit dilation) (docs/11 A33)

`dl` 行は FSHIFT を反映するだけでなく、SAD は RFSW ON のとき `dvfs=FSHIFT` を全要素の `dv` に足す (tsetdvfs)。`_Phys.dvfs` を追加し、`emittance6d(fshift=)` から設定した (`run_emit` が FFS の `FSHIFT` を渡す)。
FSHIFT=1e-6 の ESRF (RFSW+RADCOD) の報告 37 行が SAD と一致 (`bench_level26c`: 閉軌道 δ=−4.328e-3、U0=.2942751 MV、dl=−.6031200 mm)。
`Emittance[]` (`emit_extra.emittance_rules(fshift=)`、`OrbitDilation` 規則 = +FSHIFT·周長) と CALC6D (`calc6d(fshift=)`) も同じ `dvfs` を使う (`bench_level26c`: F1〜F6 が ≤1.5e-10 で一致)。未対応: 追跡 (track6/track6n) の dvfs。

### 照合 (`bench_level26a.sad`)

ESRF `pysad/tests/kek-ls_v03_12.lat` を連結した 1 ファイル、`ON RAD FLUC COD`、`USE RING`、`RFSW/NORFSW × RADCOD/NORADCOD` の順に `CELL; CALC; EMIT;`。参照 `bench_level26a.result`、pysad `bench_level26a_pysad.result`。
比較は `cmp_emit_report.py` (行単位。数値以外の部分は空白を無視して完全一致、数値は |差|≤1e-11 または相対 ≤1e-9)。結果: **不一致 0 行、雑音の桁だけの差 38 行**。

雑音の桁 (再現不能・しない): Closed orbit の y/py (1e-21〜1e-37)、Extended Twiss の R1〜R4・PSI・ZPX 等 (1e-16〜1e-39)、εy (SAD 1.4e-31 / pysad 0 など)、Beam tilt (1e-13 rad)、`Imag.tune` と RFSW OFF の z の `Real tune` の 0 の符号 (`-0.0000000` と ` 0.0000000`)。
SAD の固有値分解 (teigenc) の丸めで決まる量で、論理は同じ。`<work-esrf>/disp.sad.out` (disp スクリプトの SAD 出力) の報告 1 つ目も `bench_level26a.result` の R1 と行単位で一致した。

## 2. 課題 B: 四極・薄い多極の解析行列

### 原因の確認

`lattice._one_turn_cod_matrix` (CALC の閉軌道・Twiss の一周行列) は、δ≠0 または軌道≠0 の四極と、ずれ・ロール付きの要素を `_num_jacobian` (中心差分 eps=1e-7) で作っていた。SAD は CALC で数値微分をしない (各 `t*e` が解析行列を積算する)。
`emit6d` の解析ヤコビアン (`_USE_ANALYTIC_JAC` の枠) は、ずれ・ロールなしの四極、tqlfre なし、の場合だけだった。

作業前の切り分け: 22e の CELL の差 5e-10 は、BY で 1.5e-10 (SAD の `$FORM` 12 桁表示) など、周期解の最後に評価した軌道での一周行列の数値微分の雑音 (相対 1e-10) と一致する。
3C の EX/EPX は上のとおり別途解消済みだった。

### 解析行列

| 場所 | 追加 |
|---|---|
| `emit6d._jac_tqlfre` (新) | `tqlfre.f` の `trans1` (bz=0) の直訳 (FRINGE=1,2,3 の線形フリンジ)。行 1〜4 は cod 写像の厳密な微分、z の行は SAD のシンプレクティック条件の式のまま |
| `emit6d._roll_matrices` (新) | tsolrote のロール θ (x'=c x−s y, y'=s x+c y。px,py も) の 6×6 と出口 (転置) |
| `emit6d._jac_quad_full` (拡張) | `(dx, dy, th)` を受ける: `R_出口 · [K1<0 の −π/2] · tqfrie(+) · tqlfre(入口) · tqente + tsolque 本体 ×ndiv · tqlfre(出口) · tqfrie(−) · … · R_入口`。ずれは微分に効かない。入口の tqfrie は `mfring≠2`、出口は `≠1` (tquade.f の条件どおり)。`ph=None` は dv 項なし (lattice の CALC 経路用)。FRINGE=1,2,3 を数値微分に落としていた分岐を削除 |
| `emit6d._tquade` (krad 経路) | ずれ・ロール・K1<0 の定数行列で挟む `jfun` (以前はずれ・ロールがあると数値微分)。tqlfre の 2 つの piece に `_jac_tqlfre` |
| `lattice._one_turn_cod_matrix` (四極) | `_num_jacobian(_f, cod)` → `_jac_quad_full(cod, …, ph=None, dx, dy, ROTATION)`。解析式が作れない (pz が潰れる、K1≈0) ときだけ数値微分 |
| `lattice._thin_cod_matrix` (薄い多極) | `_num_jacobian(_f, cod)` → `emit6d._jac_thin` (tthine の各段のヤコビアンの積)。ずれ・ロールは従来の `_tchge_jac` (SAD の癖を含む) |

emit6d の薄い多極 (`_tthine_gen`) は既に `_jac_thin` + `_tchge_jac` だった。`calc6d`/`touschek`/`emit_extra` は `emit6d._tthine`/`_tquade` を呼んでおり、以前の記述 (各自のコピー) は古かった。共通化は不要。

### 検証

- **ヤコビアンの単体**: K1=±、L=0.25/0.6、(ずれ、ロール、両方、なし) × FRINGE=0,1,2,3 × DISFRIN の有無 を、COD 写像の数値微分 (刻み 1e-6、z の行は SAD が組み直すので除く) と比べると最大 2e-12。z 列 (5,6) は SAD の式と数値微分で 1e-7〜3e-6 (これは元からで、SAD の `trans1(5,6)` を使うので正しい)。`pysad/tests/test_emit_report.py`。
- **bench_level22e** (SAD との行比較, `cmp_lines.py`): 最大 **5.3e-10 → 1.4e-11** (`CELL_dp` 7.5e-12、`INS1_at_100` 1.4e-11)。
- **bench_level26b** (新規。小リング、QF=FRINGE 3、QD=FRINGE 1、1 周に別名の QFX (FRINGE 2、DX/DY/ROTATE)、QDX (FRINGE 3、DX/DY)、SFX (六極に DX/DY/ROTATE)、DP0=0 と 0.003、`Emittance[]`):
  CALC の Twiss・R1〜R4・チューン・閉軌道 (A/B/B2 の 16 行) は **≤5e-13** (SAD と 10 桁一致)。`Emittance[]` (C): ε・σE・U0 の行 1.3e-10、チューン 1.3e-9、TuneShiftByRadiation 3e-11、減衰率 1e-12、OrbitAtExit 9e-10 (絶対)。後者 2 つ (チューン・閉軌道) は SAD の tcod の収束 (RFSW OFF は codw が 100 分の 1) の違い。
  数値微分に戻して (`probes_emit/run_old.py`) 比べると、CALC の B/B2 が 2.4e-10、**C のチューンが 1.1e-5、ε・σE・U0 の行が 2.7e-7、TuneShiftByRadiation が 3.9e-8、減衰率が 1.4e-7**。tqlfre の解析ヤコビアンが無かった emit 経路の差で、この bench がそれを検出する。
- **回帰 (作業前 base_w0 と作業後 base_w3 (どちらも `base_w0/run.sh` で作る。13a 13b 14a 16a 16b 18a 22a 22e 23a 25a〜e 3c) の diff)**: 13a/13b/14a/16a/16b は出力が同一 (14a/16a は報告の新しい行が増えただけ)、18a は時間行だけ。22a 1.2e-11、23a 2.8e-10 (SAD との差は 2.7e-10 → 1.7e-11 に縮んだ)、25a 1.5e-10、25b 8e-14、25c 9e-11、25d 9.6e-10 (SAD との差 9.7e-10 → 4.5e-11)、25e 3.7e-9 (SAD との差は 4.5e-8 のまま、SAD の δ の収束差)。3C 1.8e-10。
  25a〜e の SAD との最大差は悪化していない (改善または同じ)。

## 3. テスト・道具

- `pysad/tests/test_emit_report.py` (新規 37 本): ヤコビアン (四極 K1±、FRINGE 0〜3、ずれ・ロール)、`_jac_tqlfre`、薄い多極のずれ・ロール付き行列、**`bench_level26a.result` との EMIT 報告 4 通りの行単位比較**、`bench_level26c.result` (FSHIFT)。
- ベンチ (SAD の参照 `.result` 付き): `bench_level26a.sad` (報告 4 通り)、`26b` (FRINGE/ずれ/ロール)、`26c` (FSHIFT)、`26d` (RFSW OFF + RADCOD ON の `Emittance[]` と CALC6D)。`run_regression.py` に 22e、26a〜26d を追加 (`check_lines` = ラベル付き行の最大相対差、`check_emit_report` = 報告の行単位)。
- `cmp_emit_report.py` (EMIT 報告の行比較)、`cmp_lines.py` (Print のラベル付き数値行の比較)、`probes_emit/run_old.py` (解析ヤコビアンを無効にして実行)、`probes_emit/test_jac_quad.py` (ヤコビアンの感度)、`probes_emit/eps_test.py` (数値微分の刻み依存)。
- SAD の実行 (Git Bash から): `MSYS_NO_PATHCONV=1 wsl -d Ubuntu -- <SAD>/bin/gs file.sad > file.result` (パス変換を止める。ファイルを現在のフォルダに置いて実行すると `cd` の警告が出るが動く)。

## 4. 回帰

- `QUICK=1 PYTHONIOENCODING=utf-8 python run_regression.py`: FAIL 0 (作業中の途中、28 ベンチ)。全部入りの最終結果は `00_RESUME.md` の 2026-10-16 追記。
- `python -m pytest -q pysad/tests --ignore=…` (3 本除外): 101 pass → 138 pass (新規 37)。

## 5. 未解決・次の作業

1. **解析行列の無い piece (数値微分のまま)**: MULT (`_mult_cod_map`、CALC・emit とも)、ANGLE=0 の BEND (tsteee / tmulte)、F1 フリンジ付きの BEND (`_tblfre_cod`)、DROTATE≠0 や dz≠0 の `tchge` (`_tbrote_cod`)、空洞 (tcave)。正しいが遅く、MULT 入りの CELL は SAD と 1e-9〜 の差が出うる。
2. FSHIFT (dvfs) は EMIT・`Emittance[]`・CALC6D には入れたが、追跡 (track6/track6n) には未対応。
3. RFSW OFF の `Real tune` / `Imag.tune` の 0 の符号 (teigenc の丸め) は再現しない。0 同然の桁 (εy、Beam tilt、R1〜R4) は SAD の雑音とは一致しない (論理は同じ)。
4. `CALC6D;` コマンド単独 (その後に `CELL; CALC;` なし) の直後に SAD の Twiss 配列が 6D の値に変わる件 (bench_level26d の最初の試行で見えた。`CALC6D; CELL; CALC;` なら両者一致) は pysad が再計算しない。
5. tcod が収束しないときの SAD の挙動 (報告を出して続行) は NaN の報告の書式までは合わせた (docs/history/work_item6.md §6 の 4)。Entrance の値 (未初期化) は再現不能。
