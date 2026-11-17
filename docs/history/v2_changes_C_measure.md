# v2 の変更 (担当 C): MEASURE / TRCOD / LOGMEASURE、MAIN の TRACK TURNS=-1、B41 (2026-11-17)

対象は FFS のコマンド `MEASURE`・`TRCOD`・`LOGMEASURE` と、MAIN の `TRACK … TURNS=-1` (輸送路の 1 回通過)。
実装は `pysad/measure.py`、試験は `pysad/tests/test_v2_C_measure.py` (19 件)。

**参照について**: K64 の実行ファイル (gs) は MEASURE が何も出力せず、TRCOD は SIGSEGV で落ち (B42)、MAIN の TRACK も落ちる。
gs との数値比較は 1 件もできていない。実装は K64 のソース (`tffsa.f`、`trackb.f`、`tracke.f`、`trcoda.f`、`ttinit.f`、`ttstat.f`、
`tbdecoup.f`、`tracka.f`、`track.f`) の読みで、値は解析式と pysad の独立な経路 (CALC の Twiss、EMIT、TrackParticles) で確かめた。
表示の書式は Fortran の FORMAT 文を gfortran (WSL) で同じ書式を出して合わせた (桁の位置まで。gs の実出力との比較ではない)。

## 1. 一覧

| ID | K64 の動作 | v2 の動作 | 根拠 | 試験 | k64 モード |
|---|---|---|---|---|---|
| B42 (MEASURE) | gs は何も出力しない。ソースの意図: 先頭 MARK の入口条件 (Twiss、EMITX、EMITY、DP、SIGZ、AZ) の分布の NP 個の粒子を指定の位置の入口まで追い、重心・6×6 の分散行列・射影/非結合の emittance と Twiss・ビームサイズを表示する | ソースの意図どおりに動く。`FFS["MEASURE …"]` は `{np, {重心 6}, {分散行列 6×6}}` を返す | §3 | `test_measure_*` 5 件 | 従来どおり `[pysad] FFS command MEASURE は未対応 (無視)` (変えていない) |
| B42 (TRCOD) | gs は SIGSEGV。ソースの意図: 分布を要素ごとに追い、各要素の入口での重心を Twiss の DX DPX DY DPY に、dp との相関を EX EPX EY EPY に書く | ソースの意図どおりに動く | §3 | `test_trcod_*` 3 件 | 従来どおり「未対応 (無視)」 |
| LOGMEASURE | 直前の MEASURE の xa ya xxa xya yya を `1p,5g15.7` で表示 | 同じ (v2 だけ) | 書式 9771 | `test_measure_output_and_beam_size_vs_twiss` | 従来どおり (コマンドとして認識しない) |
| B41 | `ttinit.f` が軌道の揺らぎを `rlist(kx+26..29)` から読む。MARK の添字 26..29 は ZPX ZY ZPY DP で、JDX JDPX JDY JDPY は 34..37。py の揺らぎが DP·tgauss になる | JDX JDPX JDY JDPY を読む | MARK のキーワード表 (`sadkw.py` の MARK の並び、26 番目 ZPX … 29 番目 DP、34 番目 JDX) | `test_b41_orbit_jitter_reads_jdx_in_v2_and_zpx_in_k64` | コマンドは動かないので効かない。`measure.ttinit(b41=True)` (k64 モードの既定) で K64 の読み方を再現できる。**ソースの読みのみ** (gs で未確認)。乱数を引く回数は K64 も v2 も 4 回で同じ |
| 新規 (要採番) | `ttinit.f` の JITTER OFF (QUIET) の重心合わせの基準が、分布の作り方と合わない。px の基準は `epx·dp0 − ax/bx·ex·dp0 + dpx` だが、分布の px の平均は `epx·dp0 + dpx` (ax/bx の項は分散を含まない xi に掛かる)。結合 R1..R4 も基準に入っていない。DP0≠0 で EX≠0 かつ AX≠0 (または結合あり) のとき、QUIET の重心が設計値からずれる | 分布の作り方と同じ写像で求めた設計の重心に合わせる | n=200000 の JITTER ON の平均は v2 の基準と統計誤差の 4 倍以内、K64 の基準からは 20 倍以上ずれる (MK の例: px の差 2.5e-5、統計誤差 3e-7) | `test_ttinit_quiet_centroid_v2_is_consistent_k64_is_not`、`test_ttinit_quiet_centroid_with_coupling` | `measure.ttinit(center_fix=False)` (k64 モードの既定) で K64 の式 |
| MAIN TRACK TURNS=-1 | `track.f` は TURNS=-1 を trpt=.true.・rfsw=.true.・1 回通過にし、`tracka.f` が ttinit の分布を作って `Statistics at Entrance` / `Exit` を出す。gs は MAIN の TRACK が落ちる | そのとおりに動く (`measure.track_trpt`) | §3 | `test_main_track_trpt_entrance_exit` | 従来どおりリングの経路 (`test_main_track_turns_minus1_k64_mode_unchanged`) |

## 2. 仕様 (K64 のソースの読み)

### 初期分布 (`ttinit.f` → `measure.ttinit`)
- 先頭の要素が MARK でなければ ` Missing entrance condition of beam.` を出して何もしない。
- 使う値: 先頭 MARK の AX BX AY BY EX EPX EY EPY R1 R2 R3 R4 DETR DX DPX DY DPY AZ、大域変数 EMITX EMITY SIGZ DP (FFS では `tfsetparam` が MARK に写す)、
  DP0 (= MARK の DDP)、NP (粒子数、既定 1000)。リングで CELL の CALC をした後は、周期解が先頭 MARK に入っているので整合した分布になる。
- 粒子ごとの乱数の順: pk (z)、dp、x、px、y、py。
  - z = SIGZ·pk、dp = DP0 + DP·(g − AZ·pk)/√(1+AZ²)。GAUSS OFF (既定) なら g は ±1 の一様乱数 (rms は DP/√3)、GAUSS ON なら正規乱数。
  - x = σx·g + EX·dp、px = σpx·g − AX/BX·(σx·g) + EPX·dp (σx = √(BX·EMITX)、σpx = √(EMITX/BX))。y も同じ。
    TRGAUSS OFF (TRUNI) なら単位円内の一様乱数 (rms は σ/2)。
  - 結合: R1..R4 と DETR で (x, px, y, py) を混ぜ (DETR > 1 なら x と y のモードを入れ替え)、DX DPX DY DPY を足す。
  - 軌道の揺らぎ: DX += JDX·正規乱数 (DPX DY DPY も同様。B41)。
- JITTER OFF (QUIET) なら、標本の重心の統計的なずれを引く。JITTER の既定は ON。
- TWAKE / LWAKE が ON のとき K64 は pk=0 (全粒子 z=0) にする。そのための wake 用の分布はコメントアウトされていて残っていない。
  v2 の MEASURE は z の分布を普通に作る (TrackParticles の wake がそのまま効く)。`measure.ttinit(wake=True)` で K64 の z=0 を再現できる。

### 追跡 (`tracke.f` → `measure.Tracker`)
- 追跡は FFS の TrackParticles (`fn_track`) と同じ経路。フラグ RAD FLUC RFSW などはそのまま効く。落ちた粒子 (APERT、LOSSAMPL) は統計から除く。
- 位置が MARK ならその 1 つ先、直前が MARK なら 1 つ手前まで追う (MARK は恒等なので分布は同じ)。
- FIXSEED (既定) なら乱数の種を MEASURE / TRCOD の前の値に戻す (同じコマンドは同じ分布)。MOVESEED なら進める。
- RFSW OFF では区間の追跡の後で z=0 になる (K64 の tturn と同じ。TrackParticles と共通)。
- 生き残りが 1 個以下なら統計は出さず、xa などは 1e20。

### 統計と表示 (`ttstat.f`、`tbdecoup.f`)
- 重心 (x, px/p0, y, py/p0, z, dp/p0) と分散行列 (1/n で割る。下三角を表示)。
- `x-y projected(coupled) parameters`: 分散 (dp との相関 <u·dp>/<dp²>) の寄与を引いた 2×2 から emit、β、α。ex epx ey epy は相関。
- `x-y decoupled parameters`: x-y を非結合にした emitu emitv と β α、R1..R4、detr。emitu・emitv は固有エミッタンス。
  **emitx と emity がほぼ等しいとき、非結合化は不定になる** (tbdecoup の式の分母 detp−detq が 0 に近い)。R1..R4 が大きな値になるが、K64 の式どおり。
- sigx sigy tilt (x-y の楕円の傾き)、sigpx sigpy、sigp/p、sigz、dp/p/z (= <z·dp>/sigz)。
- 実効面積 (tlum) と垂直の 4・6・8 次のモーメントは、K64 の MEASURE・TRACK からは呼ばれない (`lum=.false.`)。
  `measure.ttstat_lines(…, lum=True)` で出せる (奇数番と偶数番の粒子を二つのビームとみなす)。

### コマンド
- `MEASURE [位置]` (省略形 `MEA`): 位置は要素名・`NAME.n`・`$$$`・`^^^`。省略すると前回の MEASURE / TRCOD の位置 (最初は終端)。位置の入口での統計を表示。
- `TRCOD [位置]` (省略形 `TRC`): JITTER を OFF 扱いにして分布を作り、位置 1 から指定の位置まで、各要素の入口での重心を Twiss の DX DPX DY DPY に、
  <u·dp>/<dp²> を EX EPX EY EPY に書く (<dp²>=0 なら EX.. は書かない)。最後に位置での統計を表示。位置より後ろと、他の Twiss (BX など) は変えない。
  次の CALC で Twiss は計算し直される (K64 と同じ)。
- `LOGMEASURE` (省略形 `LOGM`): 直前の MEASURE / TRCOD の xa ya xxa xya yya (重心 x、y と <xx> <xy> <yy>)。
- `FFS["MEASURE 位置"]`: `{np, {重心 6}, {分散行列 6×6}}` (trackb.f の ret)。出力ユニット 0 (既定) では表示しない。
- MAIN の `TRACK USE=line NPART=n TURNS=-1;`: MARK の定義のキーワード EMITX EMITY DP SIGMAZ と Twiss から分布を作り、`Statistics at Entrance`、
  1 回通過、`Statistics at Exit`。戻り値 (`ffs.track_results`) に entrance / exit の (重心, 分散行列)、final、alive。乱数の種は戻さない。

## 3. 根拠 (数値)

種は既定の 17。すべて v2 モード。

| 確かめたこと | 独立な基準 | 結果 |
|---|---|---|
| 初期分布の 2 次モーメント (n=100000、GAUSS、分散と AZ あり) | <x²>=εβ+(ηδ)²、<x px>=−αε+ηη'δ²、<px²>=εγ+η'²δ²、<zδ>=−σz·DP·AZ/√(1+AZ²)、<δ²>=DP² | 統計誤差 4/√n = 1.3% 以内 |
| 一様分布 | dp の rms = DP/√3、TRUNI の横方向の rms = σ/2 | 2% 以内 (n=60000) |
| tbdecoup の非結合の emittance (R1..R4 = 0.1, 0.3, −0.05, 0.2 を入れた分布、n=200000) | 4×4 行列 Σ·J の固有値 | 相対 1e-10 以内で一致。入れた R1..R4 (±5e-3)、β α (2%)、η (±2e-3) が戻る。emitu = 9.976e-9 (入力 1e-8)、emitv = 2.003e-9 (入力 2e-9) |
| MEASURE の出口の分散行列 (輸送路 M0 → MM: QUAD 2、BEND 1、ステアリング 1。ε = 1e-12・2e-13、DP = 1e-6、NP = 1500) | M Σ0 Mᵀ。M は CALC の Twiss (β α ψ η) から作った行列、Σ0 は同じ種の初期標本 | 規格化した最大の差 1.9e-7。emittance の保存 相対 3e-8 |
| MEASURE のビームサイズ (輸送路の終端、ε = 1e-8・2e-9、DP = 1e-3、NP = 8000) | σx = √(εx·BX + (EX·DP)²)、σy = √(εy·BY) (CALC の Twiss) | 3% 以内 (統計誤差 1/√(2N) = 0.8%、色収差を含む) |
| リング (24 セルの FODO + CAVI、2 GeV、周長 139.2 m、偏向角 0.1309 rad/台) の半周 | EMIT の EMITX = 2.564e-7 m、SIGE = 4.892e-4 と Twiss: σx = 1.833 mm、σy = 0.1046 mm | MEASURE: σx = 1.803 mm (−1.6%)、σy = 0.1052 mm (+0.6%)。NP = 5000 |
| TRCOD の軌道 (入口 DX = 1 mm、DPX = −0.2 mrad、DY = 0.5 mm、DPY = 0.1 mrad、ステアリング K0 = 1e-3。軌道は最大 7.7 mm) | CALC (INS) の DX DPX DY DPY | 全要素で差 5.4e-11 以下 |
| TRCOD の分散 (ε = 1e-16、DP = 1e-5、NP = 2000) | CALC の EX EPX (EX は最大 0.59 m) | 差 5e-5 以下 (統計誤差) |
| TRCOD の重心 (DP = 0) | 入口の重心の 1 粒子の TrackParticles | 差 1e-10 以下 |
| MAIN の TRACK TURNS=-1 (NPART = 3000) | emittance の保存、出口の β = 入口の標本の (β, α) を Twiss の行列で運んだ値 | emittance 相対 9e-7、β 2e-3 以内 |
| 実効面積 (lum=True、n=40000 の正規分布) | σx·σy (同じ正規分布 2 つの重なり)、垂直のモーメント 3・15・105 | 10% 以内、モーメント 5%・15%・40% 以内 |

リングの TRCOD では、六極がなくても重心は 0 にならない (上のリングで最大 1.4 µm)。運動量の広がりの 2 次の分散 (η₂·<δ²>) による実際の重心のずれで、
誤りではない。TRCOD は「有限の広がりを持つビームの重心」を返し、1 粒子の閉軌道とは違う。

## 4. 使える範囲と条件

- v2 モード (既定) の FFS で、先頭が MARK のビームライン (輸送路、リング)。NP ≥ 2。MOMENTUM を設定すること。
- 輸送路: 先頭 MARK に BX BY AX AY (と EX.. DX..) を書き、EMITX EMITY DP SIGZ を MARK か FFS の変数で与える。
  **MARK に DP を書かないと DP は 0.01 (既定) で、運動量の広がり 1% (一様) の分布になる**。0 にするには `DP=0;`。
- リング: `CELL; CALC;` の後 (周期解が MARK に入る)。`EMIT;` の後なら EMITX・SIGE・SIGZ が大域変数に入っている (DP=SIGE; と書く)。
- 追跡は TrackParticles と同じ経路なので、TrackParticles が扱える要素・フラグの範囲で使える (RAD、FLUC、RFSW、APERT の損失)。
- 速さ: NP = 4000・14 要素の輸送路の例 (MEASURE 3 回と TRCOD) が全体で 1.3 秒、NP = 1000・194 要素のリングの CALC・EMIT・TRCOD (要素ごとに追跡) が全体で 0.7 秒。

## 5. 使えない範囲・未確認

- **gs との比較は無い** (B42)。K64 と数値が合うことは確かめていない。乱数を引く順は K64 のソースどおりだが、粒子ごとの値を gs と比べたことはない。
- `TRCOD FIXCOD` / `MVCOD` / `STNCOD` / `ALIGN`: 未対応 (エラー文を出して無視)。K64 は基準軌道を運動量ずれ用の Twiss 配列 `twiss(*,ndim,*)` に退避するが、
  pysad の Twiss はこの配列を持たない。ALIGN (`nalign`) は未翻訳。
- `PUTMEASURE` / `GETMEASURE` (パイプ経由のデータの受け渡し `nfputm` / `nfgetm`) と `ORBITJITTER`: 未翻訳。
- `DRAW` の中の粒子追跡 (`tdrwdt.f` が tracke を使う) はつないでいない。
- TWAKE / LWAKE ON の MEASURE: 動くが (TrackParticles の wake)、独立な検証はしていない。K64 は z=0 の分布にするので、K64 の動作とも違う。
- TRGAUSS OFF、DETR > 1 (モードの入れ替え) の分布: ソースどおりに書いたが、DETR > 1 は試験していない。
- MAIN の TRACK TURNS=-1 は粒子種の語 (PROTON など) を見ない (MASS・CHARGE を使う。PROTON と書くと注意を出す)。
- 表示の `GCUT =` などの autofg の桁は pysad の `fmt_real` による。gs の実出力とは比べていない。

## 6. 要判断

1. **TRPT の副作用**: K64 の MEASURE は `trpt=.true.` にしたまま戻さない。v2 は追跡の間だけ ON にして戻す (後の計算に副作用を残さない)。
   K64 のとおり残す方がよければ `cmd_measure` の `finally` の 1 行を外す。
2. **QUIET の重心合わせの修正** (§1 の「新規」): K64 の式の誤りと判断して v2 で直した。docs/21 の番号は未採番。
3. **TWAKE / LWAKE の z=0**: v2 は z の分布を作る。K64 の z=0 は使えない wake 用の分布の名残と判断した。
4. `TRCOD` の途中の統計は K64 では装置番号 0 (標準エラー) に出る。v2 は出さない (最後の位置の統計だけ)。

## 7. 使い方の例 (実際に走らせたスクリプト)

輸送路 (`run_sad_string` で実行。約 1 秒):

```
MOMENTUM = 1 GEV;
MARK M0 = (BX=8 BY=3 AX=0.5 AY=-0.3 EX=0.2 EPX=0.01 EMITX=1E-8 EMITY=2E-9 DP=1E-3 SIGZ=3E-3) MM=();
DRIFT D1 = (L=1.0);
QUAD QF = (L=0.3 K1=0.25) QD=(L=0.3 K1=-0.25);
BEND B1 = (L=1 ANGLE=0.05);
LINE TL = (M0 QF D1 B1 D1 QD D1 MM B1 D1 QF D1 QD D1);
FFS;
USE TL;
INS; CALC;
NP=4000; GAUSS;
MEASURE;                       ! 終端での統計
MEASURE MM;                    ! MM の入口での統計
LOGMEASURE;                    ! xa ya xxa xya yya
r=FFS["MEASURE $$$"];          ! {np, 重心, 分散行列}
Print[Sqrt[r[[3,1,1]]]];       ! σx = 7.647e-4 (√(εβ+(ηδ)²) = 7.626e-4)
TRCOD;                         ! 重心の軌道を Twiss の DX.. に、分散を EX.. に書く
Print[Twiss["DX","$$$"]," ",Twiss["EX","$$$"]];
```

出力 (最初の MEASURE):

```
 Statistics at $$$:   particles = 4000
 RAD: F, RFSW: T, GAUSS: T, DP =   .001000, DP0 =   .000000, GCUT = 1.000E300
                    x     px/p0         y     py/p0         z     dp/p0 
   C of M :  7.247E-06 1.454E-06 5.253E-06-3.744E-07-5.750E-05 1.840E-05
        x :  5.847E-07
    px/p0 :  1.151E-07 2.298E-08
        y :  3.167E-09 5.855E-10 3.375E-08
    py/p0 : -2.765E-10-4.597E-11-2.952E-09 3.763E-10
        z : -2.124E-09-1.216E-10 3.525E-09-5.811E-10 8.979E-06
    dp/p0 :  5.156E-07 1.016E-07-6.291E-10 1.948E-10-5.474E-08 9.879E-07
 x-y projected(coupled) parameters:
 emitx: 1.0091E-08 bx:  31.28     ax: -6.149     ex: 0.5219     epx: 0.1029    
 emity: 1.9959E-09 by:  16.91     ay:  1.479     ey:-6.3680E-04 epy: 1.9713E-04
 x-y decoupled parameters:
 emitu: 1.0093E-08 bu:  31.28     au: -6.149     eu: 0.5219     epu: 0.1029    
 emitv: 1.9934E-09 bv:  16.91     av:  1.478     ev:-7.5432E-03 epv: 8.9707E-04
    r1:-3.1873E-02 r2: 9.4578E-02 r3: 7.0391E-03 r4:-2.8911E-02detr: 2.5573E-04
   sigx: 7.6467E-04  sigy: 1.8371E-04   tilt:-5.7486E-03
  sigpx: 1.5158E-04 sigpy: 1.9400E-05
 sigp/p: 9.9395E-04  sigz: 2.9965E-03 dp/p/z:-1.8268E-05/sigz
```

(CALC の Twiss は BX = 32.16、AX = −6.284、EX = 0.510。標本 4000 個の統計誤差と DP = 1e-3 の色収差の範囲。)

リング (EMIT の平衡値の分布を半周):

```
MOMENTUM = 2 GEV;
MARK M0 = () MID=();
DRIFT D1 = (L=0.5);
QUAD QF = (L=0.4 K1=0.32) QD=(L=0.4 K1=-0.32);
BEND B1 = (L=1.5 ANGLE=0.1308996938995747);
CAVI CA = (VOLT=5E5 HARM=40);
LINE CELL1 = (QF D1 B1 D1 QD D1 B1 D1);
LINE RING = (M0 12*CELL1 MID 12*CELL1 CA);
FFS;
USE RING;
CELL; CALC;
EMIT;
EMITY=EMITX*0.01; DP=SIGE;
NP=5000; GAUSS;
MEASURE MID;                   ! sigx: 1.8032E-03 (解析値 1.833e-3)、sigy: 1.0517E-04 (1.046e-4)
```

MAIN (輸送路の 1 回通過。Entrance と Exit の統計):

```
MOMENTUM = 1 GEV;
MARK M0 = (BX=8 BY=3 AX=0.5 AY=-0.3 EX=0.2 EMITX=1E-10 EMITY=2E-11 DP=1E-5 SIGMAZ=3E-3);
DRIFT D1 = (L=1.0);
QUAD QF = (L=0.3 K1=0.25) QD=(L=0.3 K1=-0.25);
BEND B1 = (L=1 ANGLE=0.05);
LINE TL = (M0 QF D1 B1 D1 QD D1 B1 D1 QF D1 QD D1);
ON GAUSS;
TRACK USE=TL NPART=3000 TURNS=-1;
```

Python から部品を使う例:

```python
from pysad import measure as ms
from pysad.sadrandom import SadRandom
mk = dict(AX=0.5, BX=8, AY=-0.3, BY=3, EX=0.2, EMITX=1e-8, EMITY=2e-9, DP=1e-3, SIGZ=3e-3)
x = ms.ttinit(mk, 10000, SadRandom(17), gauss=True)        # (10000, 6) の正準変数
lines, sa, ss, es = ms.ttstat_lines(x, 'Entrance', lum=True)
print('\n'.join(lines))
```

## 8. 変更したファイル

- `pysad/measure.py`: 全面的に追加 (tlum / tluma は元のまま)。tbdecoup、ttinit、ttstat、ttstat_lines、mark_values、Tracker、cmd_measure、cmd_trcod、
  cmd_logmeasure、ffs_command、track_trpt。
- `pysad/ffscmd.py` (3 か所、いずれも v2 だけで効く): `is_word` に MEA_SURE / TRC_OD / LOGM_EASURE の省略形、`_exec_word` の REF_ERENCE の分岐の後に
  `measure.ffs_command` の呼び出し、`fn_ffs` に MEASURE の戻り値 (`_measure_ret`)。
- `pysad/maintrack.py` (1 か所、v2 だけ): `run_track` の TURNS=-1 で `measure.track_trpt` へ。
- `pysad/tests/test_v2_C_measure.py`: 新規 19 件。
