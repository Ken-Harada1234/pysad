# -*- coding: utf-8 -*-
"""example/kek_ls の各例の NAME.md を作る。文章は下の INFO、結果の抜粋と gs・pysad の判定は gs・pysad の出力から自動で入れる。
使い方 (リポジトリ直下から): py -3.12 example/tools/kek_ls_md.py [NAME ...]    (先に run_kekls.py で出力を作る)"""
import os, re, sys, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import kek_ls_cmp as C

EX = C.EX
TIMES = {}


def times(name, d):
    p = os.path.join(EX, d, name + '.time.txt')
    out = {}
    if os.path.exists(p):
        for l in open(p, encoding='utf-8'):
            m = re.search(r'gs (\d+) s, pysad (\d+) s', l)
            if m:
                out = {'gs': int(m.group(1)), 'py': int(m.group(2))}
    return out


# 名前: dict(kind, what, how, cmp, limits, show (結果の抜粋の接頭辞の正規表現), nshow, png, orig)
INFO = {
 'KEKLS_cell_calc_emit': dict(
  kind='FFS の `cell calc emit`', orig='(全台本に共通の冒頭。個別の元の台本は無い)',
  what='周回ライン `RINGC` (標準セル 20 個) と `RING` (入射セル込み) の周期解 (`cell calc emit`)。周長・チューン・εx・σδ・σz を、`LINE[]`・`Twiss[]`・変数 `EMITX`・`SIGE`・`SIGZ`・`MOMENTUM` で読む。',
  how='`USE` するライン名を `FFS["USE " // nm]` で切り替えて 2 回回す。周長は `LINE["S","$$$"]` (終点の S)、チューンは `Twiss["NX","$$$"]/2/Pi`。RING は RINGC より 0.9 mm 長く (入射セルの `DRID1INJ`)、εx・σδ が少し大きい。',
  cmp='全項目が全桁一致。',
  limits=['pysad の v2 (既定) モードは εx・σδ が K64 と 0.04% 違う (RINGC εx 134.37 pm rad、σδ 6.4045e-4。K64 のバグ B75 の修正。`docs/28_v2_vs_k64.md`)。gs と同じ値は `PYSAD_COMPAT=k64`。',
          'エネルギー 3 GeV、周長 570.72 m、RF 2 MV、ハーモニック数 952 の設計値のラティス。誤差なし。'],
  show=r'(circumference|tune nx|EMITX|SIGE|SIGZ|LINE)', nshow=20),
 'KEKLS_cell_optics': dict(
  kind='FFS の `cell calc`、`Twiss[]`、`LINE[]`', orig='disp_optics (要素ごとの表)',
  what='標準セル 1 個 (`ARCA`、28.536 m) の要素ごとの種類・位置・長さ・偏向角・K1・K2・βx・βy・ηx の数表。`../lattice/KEKLS_cell_layout_optics.svg` の元データ (`example/tools/kek_ls_figs.py` が gs の出力から描く)。',
  how='先頭の数行 (セルの長さ、チューン、×20 の位相進み、`MKS`・`MKC`・`MKSFA1` のツイス) と、`E ` で始まる行 (名前、種類、S、L、ANGLE、K1、K2、βx、βy、ηx。値は要素の入口)。`ARCA` は `cell calc` (周期解、`emit` なし。閉じていないラインでは `emit` は使えない)。',
  cmp='全桁一致 (数表 160 行)。',
  limits=['図 (PNG) は gs の出力から作った。pysad の出力でも同じ図になる。'],
  show=r'(ARCA|MKS |MKC |MKSFA1)', nshow=8, png=['../lattice/KEKLS_cell_layout_optics.svg']),
 'KEKLS_match_cell_tune': dict(
  kind='FFS の `fit`・`free`・`go`', orig='matching2_01 (前半)',
  what='標準セル (`ARCA`) の位相進みを (48.58/20, 17.62/20) に合わせる、セルのチューンのマッチング。`MKC` の分散 EX の目標 exmkc を 0.004、0.002、0 と変えて 3 回行う。',
  how='条件 `fit NX`・`NY`・`MKC EX exmkc`・`MKS EX 0`、変数 `free Q*` (四極の K1)。各回の結果が `exmkc target ... result: ring tune ...` の行、そのときの四極の K1 が `K1:` の行。ring tune は セル×20 で (48.58, 17.62) に戻る。',
  cmp='△ 最初の 1 点 (exmkc = 0.004) は全桁一致。2 点目以降は残差 (チューン、EX) は同じ程度 (1e-6〜1e-7) に収束するが、変数の K1 は 1e-3 違う (変数が条件より多く、解が一意でない。最適化の経路で決まる)。',
  limits=['元の台本の `free BPI6 K1`・`free BPI7 K1` は、v03_68 の偏向磁石の名前 BPI6H・BPI7H に一致せず、K1 が動かない (gs・pysad とも)。元の挙動を残し、既定は `freebend = 0`。`freebend = 1` で BPI6H・BPI7H も変数にする。',
          '元の台本の EMITX の FitFunction を除いた (閉じていない ARCA では EMITX は NaN で何も効かない)。ただし NaN の FitFunction があると pysad は fit が動かない (元の台本のまま走らせた比較で確認。gs は動く)。pysad の FitFunction の扱いの違いの疑い (README の要確認)。',
          '整理前後の比較: 元の matching2_01 (ラティスを貼り、DAsv を除き、FitFunction を除いたもの) の gs と、この例の gs で、1 点目の値は一致、2・3 点目のチューンは 5e-6 以内で一致。'],
  show=r'(design|target|SF-SF|exmkc)', nshow=10),
 'KEKLS_match_chrom': dict(
  kind='FFS の `fit`、`FITP`、6 極のクロマティシティ補正', orig='matching2_01 (後半)',
  what='リング `RINGC` の6 極 `SDH`・`SFH` の K2 で、チューンを保ちつつ運動量 DP = 0.001 でも同じチューンにする (クロマティシティを 0 にする) マッチング。ξ を `DP0 = ±0.0005` の差分で前後に表示する。',
  how='`FIT; NX nx0 4; NY ny0 4; SD* MAX 0; SF* MIN 0; DP=0.001; FREE SD*,SF*; GO; FITP 1;`。before の ξ (+0.119, −0.684) が after で (≈ 0, ≈ 0) になる。K2 の変化は SDH −12.382 → −12.437、SFH 15.134 → 15.150。',
  cmp='△ before は全桁一致。after は ξ が gs (−0.0051, −0.0065)・pysad (−0.0073, −0.0070) で、どちらも 0 に近いが違う (目標からの残りがほぼ同じ大きさ、K2 が 2e-4 違う)。',
  limits=['元の台本は前半のセル調整に続けて走る。ここは設計値の RINGC から単独で始める (整理前の前半からの連続の結果とは比べていない)。',
          '追跡 (DAsv) は含めない (`modules/modules_general.n` の `DAsv[]`。`DAturns` を小さくして使える)。'],
  show=r'(before|after)', nshow=6),
 'KEKLS_match_injcell': dict(
  kind='FFS の `ins calc`・`fit`・`free`', orig='matching2_03_INJ',
  what='入射セルの半分 `ARCINJH` の四極 `QS0`〜`QS5` の K1 を変数にして、入射点 `MKSI` で位相進みを標準セルの半分、AX = AY = EPX = 0、βx = 30 m に合わせる。v03_68 の `QS*` はすでにこの解なので、K1 を 1.2 倍して始め、解に戻ることを見る。',
  how='`before` が 1.2 倍にしたとき (位相 1.2413/0.6770、AX 2.73、βx(MKSI) 53 m)、`iteration 1` で解 (位相 1.21450/0.44050、AX・AY・EPX ≈ 0、βx 30.0 m) に戻る。2 回目以降は動かない (元の台本は 50 回繰り返す)。末尾の `type Q*` に K1 の値 (元の設計値と 1e-5 以内)。',
  cmp='○ 1 回目の結果が 1e-10 以内で一致 (K1 も 12 桁)。',
  limits=['整理前後の比較: 元の matching2_03_INJ (50 回、設計値から) の gs の最終の `QS*` の K1 は設計値と 1e-6 違い (最適化の残り)、この例の解 (摂動から 1 回) とも 3e-6 以内。βx(MKSI) は 29.9999 (元 29.999896、この例 29.999925)。',
          '入射セルの AX・AY の目標はあくまで 0 (入射セルの両端を標準セルの対称点につなぐ条件)。入射点の設計値 βx = 30 m の理由は確認していない。'],
  show=r'(ARCA end|design|before|iteration)', nshow=8),
 'KEKLS_match_ring_geo': dict(
  kind='FFS の `fit GX`・`GY`、ジオメトリ', orig='matching2_02_geo (ユーザーの指定は matching_02_geo)',
  what='入射セルの長直線部の長さ `DRID1INJ` を変数にして、リング `RING` の周回後の位置 (GX, GY) を 0 にし、リングを幾何学的に閉じる。標準の 2.94 m から始めて、v03_68 の設計値 2.765460054 m に戻ることを見る。',
  how='`design` は設計値の閉じ方 (GX・GY ≈ 1e-7)、`start` は 2.94 m にしたとき (GX 0.282 m、GY 0.205 m)、`result` はマッチング後 (DRID1INJ 2.765460041 m、GX・GY ≈ 5e-8)、周長 570.72092 m。',
  cmp='○ 全項目 5e-10 以内で一致。',
  limits=['元の `USE RINGI` は v03_68 にはなく (RING)、元の台本のまま走らせると gs が `Missing beamline in USE` で止まる (確認済み)。`RINGI` を `RING` に直した。',
          '元の台本は設計値から始めるので、何も動かない。始めの値を 2.94 m にした手順は、この例で足した。',
          '元の台本の `disp geo $$$` は、gs では `$$$` 1 行だけの表になる (要素ごとのジオメトリの表は出ない)。'],
  show=r'(design|start|result)', nshow=4),
 'KEKLS_match_repeat': dict(
  kind='診断 (modules_cell_match.n の再現性)', orig='(元の台本には無い。tune_survey_cell の確認のために足した)',
  what='同じセルのマッチング (`InitMag[]` → `MatchingARCA[1]` ×2 → `MatchingARCA[0.001]` → `ChromCor[]`) を 3 回繰り返し、毎回同じ結果になるかを見る。',
  how='A・B・C・D の行 (回数、K1 または K2 の値)。gs は 3 回とも全桁同じ。',
  cmp='△ (再現性の診断。pysad は 1 回目と 2 回目以降で値が違う (A の QFA1H: 0.5198 → 0.5201 → 0.5201。gs は 0.51151 を 3 回)。さらに、gs と pysad の最初のマッチングの結果も違う (A のリングのチューン gs 48.58012、pysad 48.5909 (収束していない))。FitFunction (6 極の対の位相条件) を fit 条件と同時に使う `MatchingARCA` での FIT の扱いの差の疑い。pysad のコードは変えていない)。',
  limits=['原因は未確認 (pysad の FIT 状態の持ち越し、または FitFunction の重みの扱い)。`KEKLS_tune_survey_cell`・`KEKLS_scan_target` の gs と pysad の違いはこれが原因とみられる。',
          '元の台本 (matching2_01 の `FitFunction` が NaN のとき pysad が動かない) と同じ現象の疑いがある。'],
  show=r'[ABCD] ', nshow=12),
 'KEKLS_scale_mini': dict(
  kind='FFS の `fit`・`FitFunction`・`LINE[]` への代入', orig='matching23_v03_35_mini (v03_35 用)',
  what='標準セル全体の長さを一様に factor 倍して縮小し (1 → 0.5)、K1 を 1/factor 倍にしてから、四極 (と、勾配つき偏向磁石 BPI6H・BPI7H) の K1 を再調整して、セルの位相進み (NX・NY) を元の値に保つ。',
  how='`design` が縮小前、`factor ...` の行が各縮小率でのセル長と位相進み、`BX*sin-0.3`・`sin_y` が FitFunction の 2 つの項 (6 極の対の位相条件)、`K1:` の行が再調整後の K1。factor 0.5 でセルは 14.268 m、四極の K1 は 2 倍前後。',
  cmp='6 点とも位相進み (2.42900, 0.88100)・FitFunction の項・K1 が 2e-6 以内で一致。',
  limits=['元は v03_35 (1.5 GeV、別のセル) 用。v03_68 に合わせて、ライン名 `ARCA1` → `ARCA`、要素名 `BPI6`・`BPI7`・`QFA4A`・`QFA4B` → `BPI6H`・`BPI7H`・`QFA4`、`MKSFA1`・`MKSFA2` → `MKSFA1.1`・`MKSFA1.2` に直した。元の v03_35 は動かして確認していない (別のラティス)。',
          '元の `free BPI6 K1` は名前が合わず効かないので `BPI6H`・`BPI7H` にした (`freebend = 0` で元の挙動)。',
          '元の K1 の縮小は 1/factor 倍 (形のとおりなら 1/factor² 倍)。再調整で補う。',
          '例は factor を 0.1 刻みの 6 点にした (元は 0.01 刻み 51 点)。元の範囲全部は `{factor,1,0.5,-0.01}`。'],
  show=r'(design|factor)', nshow=8),
 'KEKLS_scale_ncell': dict(
  kind='FFS の `fit`、`LINE[]` への代入、`LINE` の定義', orig='matching21_v03_35 (v03_35 用。matching21・matching22・matching21_v03_35_fix は同系)',
  what='標準セルの数を j (20, 18, 16, 12, 10) に変える。セル 1 個の偏向角を 20/j 倍し、四極 (と BPI6H・BPI7H) の K1 を調整してセルの位相進みを設計値に保ち、j セルのリング RINGJ の εx・チューン・周長を表示する。',
  how='`cells j:` の行がセルの位相進みと βx(MKSFA1.1)・EX(MKC)、`ring RINGJ..:` の行がリングの周長・チューン・εx。εx は j = 20 の 134.3 pm rad から、18 で 184、16 で 261、12 で 615、10 で 1064 pm rad と、(20/j)³ で増える (偏向角の 3 乗)。',
  cmp='○ 全桁 (1e-10) 一致。',
  limits=['元は v03_35 用で、FitFunction が `(EMITX − 1*1^-9)*1E9` (`1^-9` は 1。10^-9 の誤り)、さらに ARCA では EMITX が NaN なので、元のままでは条件にならない (pysad は NaN の FitFunction で fit が動かない)。この例は、元の他の版でコメントにされていた `fit NX nxcell`・`fit NY nycell` (セルの位相進みを設計値に) を条件にした。',
          'USE し直すと ANGLE が元に戻るので、リングを USE した後に ANGLE と K1 を設定し直している (元の台本は ARCA だけで、εx は出さない)。',
          'j = 14 は、位相進みを保つとリングのチューンが 34.006 (整数に近い) になり `emit` が破綻するので外した。j = 20 → 2 の元の範囲では、整数・半整数共鳴に当たる j がほかにもある。',
          'RF 空洞は 1 台 (2 MV、ハーモニック数 952)。周長が変わるので RF 周波数・バケットは設計値と違う (εx には影響しない)。'],
  show=r'(cells|   ring)', nshow=10),
 'KEKLS_tune_survey_cell': dict(
  kind='チューンサーベイ (FFS の `fit`、6 極補正)', orig='matching05 (+ modules_tune_sv_v03_68.n)。ユーザーの指定 (matching23_v03_35_mini) は縮小の台本で、チューンサーベイに当たるのはこれ (README の「元のスクリプトとの対応」)',
  what='リングのチューン (tunex, tuney) の格子の各点で、標準セルの位相進みを (tunex/20, tuney/20) に再マッチし、6 極 SDH・SFH でクロマティシティを補正して、リング RINGC のチューン・ξ・εx・6 極の強さを表示する。格子は 48.30〜48.90 (0.15 刻み、5 点) × 17.30〜17.90 (5 点) の 25 点で、整数 (49、18) と半整数 (48.5、17.5) の共鳴線のまわりを含む。',
  how='`POINT` の行に、目標 (νx, νy)、実際のリングのチューン、ξx・ξy、εx [pm rad]、K2(SDH)、K2(SFH)。gs では、εx は νx = 48.45〜48.6 付近で最小 (145〜170 pm rad)、νx = 48.30 では 365〜411 pm rad と大きい。(48.75, 17.60)・(48.75, 17.75) は6 極補正が収束せず εx が 1.4e4 pm rad に飛ぶ (図は 450 で頭打ち)。ξ が 0 から大きく外れた点は、6 極補正が収束しなかった点。',
  cmp='△ gs と pysad で、点ごとの εx・ξ・K2 が一致しない (たとえば (48.30, 17.30) の εx は gs 411、pysad 173)。マッチングの解が一意でなく、`KEKLS_match_repeat` で見たように pysad の FIT の経路が gs と違うため。リングのチューンの目標への到達は、gs では全点で 5e-4 以内。',
  limits=['元の格子は 46.05〜46.55 (0.1 刻み 6 点) × 16.05〜19 (0.1 刻み 30 点) = 180 点で、各点に DA (5000 周) を追加する。例は 25 点・DA なし。所要時間 (gs 10 秒、pysad 約 160 秒) から、元の 180 点は gs で約 1.2 分、pysad で約 19 分 (DA を除く)。DA を足すと 1 点あたり gs で数分、pysad で十数分以上と見込む (DAsv の 5000 周 × 19 運動量 × 2 方向。測定はしていない)。',
          '元の matching05 は modules_tune_sv_v03_35.n を読み込む。modules_tune_sv_v03_68.n を読み込む台本は元のフォルダに無い。v03_68 用の目標値・6 極の名前に直した (modules/modules_cell_match.n の説明)。',
          '元の `tune_survey*` (tune_survey・_linear・_momentum・_disp・_amp_dep_coef・_ncell・_write_QS) は、別のリング (6 極 S1L・S2L・S1H… と `hsx_value.csv`、チューン 15.x/5.x) の `modules_tune_sv.n` を使い、v03_68 には適用できない。'],
  show=r'POINT', nshow=8, png=['KEKLS_tune_survey_cell_gs.svg', 'KEKLS_tune_survey_cell_pysad.svg']),
 'KEKLS_ring_params': dict(
  kind='FFS の `cell calc emit` の表、`Twiss[]`', orig='calc_emit.sad (前半)、disp、disp_emit、chromaticity_plot.sad',
  what='リング RINGC・RING のパラメータ (εx、σδ、σz、減衰、αc、U0、ξ) の表示。`cell calc emit` の表 (U0、αc、減衰時間・分配、バケット高さ) と、変数と差分で出した ξ (6 極あり・自然)。',
  how='各リングの `===== NAME` の後に emit の表 (gs・pysad とも出る)、その後の `circumference`・`tune`・`EMITX`・`SIGE`・`SIGZ`・`xi` の行。ξ の自然値 (K2 = 0) は (−57.56, −89.26)、6 極ありは (+0.119, −0.684)。',
  cmp='○ 全桁一致 (6e-9 以内)。',
  limits=['ξ は `DP0 = ±0.0005` の差分。作業メモの値 (RINGC の ξ = (+0.133, −0.658)、自然 (−57.58, −88.96)) とは違う (差分の刻み 1e-4〜1e-3 では 0.118〜0.120 で変わらない。メモの値の出し方は未確認)。',
          'pysad の v2 モードは εx・σδ が K64 と 0.04% 違う。表の値の違いは lattice/README.md。'],
  show=r'(=====|circumference|tune|EMITX|SIGE|SIGZ|xi)', nshow=18),
 'KEKLS_optics_table': dict(
  kind='FFS の `disp`、`Twiss[]`、`LINE["GX"]`', orig='disp_optics、disp (`disp $$$`、`disp geo $$$`)、draw_cad.sad の入力',
  what='リング RING の `disp $$$` (終点)、`disp MKSI` (入射点)、`disp geo $$$` と、要素ごとの S・βx・βy・ηx・位相 (ψx, ψy)・ジオメトリ (GX, GY) の表 (3226 行の数表)。図 `../display/KEKLS_ring_layout_optics.png` の元データ。',
  how='`E ` で始まる行が要素ごと (名前、S [m]、BX、BY、EX、PSIX、PSIY、GX、GY。値は要素の入口)。図の左は GX–GY の配置 (長直線部は細線)、右は S に対するオプティクス。入射セルは RING の最初の標準セルの次 (S = 28.5 m から) で、入射点 `MKSI` (S ≈ 57 m) の βx = 30 m の山が見える。',
  cmp='○ 全桁一致 (1e-10 以内)。',
  limits=['`draw_cad.sad` は元の描画モジュール `drawmagvw3d.n` (元のフォルダにも無い) で `ring.vss` を書くもので、この例では扱わない。`draw` コマンドは gs・pysad とも窓か topdrawer の出力で、図は PNG に別に作った。',
          '`disp geo $$$` は、gs では `$$$` の 1 行だけの表 (要素ごとのジオメトリの表は出ない)。'],
  show=r'(tune)', nshow=2, png=['KEKLS_ring_layout_optics.png']),
 'KEKLS_calc_H': dict(
  kind='解析式 (H 関数の積分によるエミッタンス)', orig='calc_H.sad、calc_H_2.sad、calc_emit.sad (後半)',
  what='標準セル (ARCA) の偏向磁石 (勾配つきを含む) の中点の Twiss から、各磁石の H 積分と分散の積分を解析式で求め、εx = Cq γ² ΣH / Σ(θ²/L) / Jx を、`cell calc emit` の RINGC の EMITX と比べる。',
  how='磁石ごとの H 積分の行の後に、`Jx = 1.31104`、H 積分からの εx = 131.18 pm rad、SAD の εx = 134.32 pm rad、比 0.9767。',
  cmp='○ 全桁一致 (1e-10 以内)。',
  limits=['H 積分は SAD の εx より 2.3% 小さい。差の原因は未確認 (中点の Twiss から前向き・後ろ向きに積分する近似、E1・E2 の端の効果なし、減衰分配 Jx の評価など)。解析式を SAD の結果に合わせる補正はしていない。',
          '元の calc_emit.sad は勾配なしの式で、γ = 2500/0.511 (3 GeV に合わない)、途中の `stop;` で後半は動かない。calc_H.sad の Cq は `3.832/10000` (nm の単位)。',
          '整理前後の比較: 元の calc_H.sad (ラティスを貼り、NPARA だけ 4 にして) の gs の Jx = 1.3110389414、εx = 0.1311802005 (nm rad) が、この例の 131.1802005 pm rad と一致。'],
  show=r'(Jx|emittance|EMITX|ratio)', nshow=6),
 'KEKLS_chromaticity': dict(
  kind='チューンの運動量依存', orig='chromaticity_plot.sad、amp_dep_ts_dp.sad',
  what='RINGC のチューンを δ = DP0 = −3%〜+3% (0.5% 刻み) で求める。設計の6 極あり (T1) と6 極なし (T0)。',
  how='`T1 δ νx νy`・`T0 δ νx νy` の行。6 極ありは δ ±1% で Δν ≈ ±0.0x の滑らかな曲線、6 極なしは ξ = (−57.6, −89.3) の自然クロマティシティで Δν = ±1.7 (δ = ±3%) に達する。図: `KEKLS_chromaticity.svg`。',
  cmp='○ 全桁一致 (表の 26 行)。',
  limits=['元は δ を 0.1% 刻み。整理前後の比較: RING で 0.1% 刻み (61 点) にして元の chromaticity_plot.sad の前半の表と比べ、δ ≤ +2.7% で 2e-6 以内で一致。+2.8〜+3% は共鳴の近くで 0.005〜0.15 違う (元は `DP0 = 0; cell calc` で毎回周期解を求め直す。この例は直前の解から続ける)。',
          '6 極なしの表の δ = ±3% は νx が 49 や 18 の整数を越えるので、チューンの値は連続な曲線の一部に見えるだけで、その先の共鳴は無視している。'],
  show=r'T[01] ', nshow=6, png=['KEKLS_chromaticity.svg']),
 'KEKLS_amp_dep_ts': dict(
  kind='振幅依存チューンシフト (ADTS)', orig='amp_dep_ts.sad、amp_dep_ts_x.sad、amp_dep_ts_y.sad (+ modules_amp_dep_ts.n)',
  what='リングの始点のツイスから、振幅 n σ の位相 θ を 10 点で回し、各点を `ins calc` (1 周) して得たチューンの平均を求める。水平 (n = 1〜101、y = 1σ)、垂直 (n = 1〜51、x = 1σ)。',
  how='`AX n νx νy`・`AY n νx νy` の行。n = 101σ (水平 2.6 mm) で νx は 0.15 近く増え、νy が減る。垂直 51σ (0.9 mm) で νy が 0.1 以上動く。図: `KEKLS_amp_dep_ts.svg`。',
  cmp='○ 全桁一致。',
  limits=['元の y の範囲 (1〜401σ、約 7 mm まで) は、この ring では 100σ 前後から先で `ins calc` が破綻する (gs・pysad とも、NaN や 1.59e199)。例は有効な範囲にした。',
          '整理前後の比較: 元の amp_dep_ts.sad (USE RING) の gs の垂直 6 点 (1〜101) が、この例のスクリプトを RING にして走らせた値と 5e-9 以内で一致。',
          '`ResetInit[]` (初期座標を 0 に戻す) をモジュールに足した。元の台本は、失敗した点の初期座標が次の点に残ることがある。',
          '粒子の追跡との照合 (ADTS の基準は追跡) はしていない。`ins calc` の 1 周の転送の結果。'],
  show=r'(design tune|AX|AY)', nshow=14, png=['KEKLS_amp_dep_ts.svg']),
 'KEKLS_sextupole_tune_scan': dict(
  kind='6 極を変えたチューンシフト', orig='dp_sv_sx.sad (dp_sv_sx_dp_n0.03・_res・_v02_42 は同系)',
  what='6 極 SDH・SFH の K2 を設計値から (ΔSD, ΔSF) ずらし (それぞれ −3〜+3、1.5 刻み、5×5 = 25 点)、δ = −3%〜+3% (1% 刻み 7 点) で、チューンの設計値からのずれの絶対値の和 Σ|Δνx|、Σ|Δνy| を求める。',
  how='`S dsd dsf dtux dtuy` の行。設計値 (0, 0) の周りで Σ|Δν| は x が 1.25、y が 1.44。図: 25 点の地図 (`KEKLS_sextupole_tune_scan.svg`)。SD・SF とも + 方向に動かすと (3, 3) 付近で大きく悪化する。',
  cmp='△ 25 点中 20 点が 1e-8 以内で一致。残りは ΔK2 が大きい点 (±3) で、δ = ±3% のチューンが整数 (49・18) を越える点の収束が違う。',
  limits=['元は v02_42 の 6 本の6 極 (SDA・SFA・SDC・SDB・SFB・SDD) を (+dsd, +dsf)・(−dsd, −dsf) でずらし、31×31 = 961 点 × δ 13 点。基準値も v02_42 の最適値。v03_68 の6 極は SDH・SFH だけなので、基準を設計値にした。元の範囲全部は `ds0 = -3; dds = 0.2; nds = 31; dd = 0.005; nd = 12` (見込み: 25 点 × 7 δ で gs 2 秒・pysad 18 秒から、961 点 × 13 δ で gs 約 2.5 分、pysad 約 1 時間 (外挿))。',
          'チューンの動きが小さい (Σ|Δν| が小さい) 組が、ダイナミックアパーチャの大きい組とは限らない (DA は含めていない)。'],
  show=r'(design|S )', nshow=8, png=['KEKLS_sextupole_tune_scan.svg']),
 'KEKLS_scan_target': dict(
  kind='セルのマッチングの目標値のスキャン (チューンの δ 依存、ADTS)', orig='dp_sv.sad、dx_sv.sad・dx_sv_ay/by/ex.sad、dy_sv.sad・dy_sv_ay/by/ex.sad (+ _R、_v2_11)',
  what='6 極の位置 MKSFA1 でのツイスの目標 (AY 1.19、BY 6.4、EX 0.14) を 1 つずつ −10%・0・+10% に振り、セルを再マッチ → 6 極でクロマティシティ補正 → リング RINGC で ξ、δ = ±1%・±2% のチューン、ADTS (水平 61σ・垂直 31σ) を表示する。元の dx_sv・dy_sv は「ずれ」ではなく、この目標値のスキャン。',
  how='`R 名前 値 | ξx ξy | νx νy (δ = −2%, −1%, +1%, +2%) | 水平 ADTS Δν | 垂直 ADTS Δν` の行。目標を変えると ADTS の符号・大きさが変わる (gs で水平 61σ の Δνx は 0.034〜0.48)。目標 0 の 3 行 (AY 1.19、BY 6.4、EX 0.14) は同じ設定で、gs は 3 行とも同じ値 (再現性の確認)。',
  cmp='△ gs と pysad で、行ごとの ξ・Δν が一致しない。設計値の目標の 3 行は、gs は同じ値、pysad は行ごとに違う (`KEKLS_match_repeat` の pysad の再現性の問題)。',
  limits=['元の dx_sv_*・dy_sv_* は AY・BY・EX の 1 つを 3〜7 などの広い範囲で 0.1 刻みに振る (BY 3〜7 は v02 の値)。元の完全な設定は、`facs` を細かくする (1 点あたり gs で約 0.6 秒、pysad で約 13 秒)。',
          '元の振幅は x 1〜202σ、y 1〜401σ (21 点)。この ring では 100σ 前後から先が破綻するので、61σ・31σ の 1 点にした。',
          'ξ の目標を 0 にしたので、ξ ≠ 0 の行 (BY 5.76 で 2.9 など) は6 極補正が収束しなかった点。'],
  show=r'R ', nshow=9),
 'KEKLS_errors_cod': dict(
  kind='誤差 (COD・分散・チューン)', orig='print_all_error.sad、print_dk_effect.sad、print_dt_error.sad、print_cod_magnitude.sad (他の print_*.sad は補正後の表示違い)',
  what='誤差の種類ごとの影響。乱数の種 3 本 (17・111・3) × 4 通り (全部、位置ずれだけ、強さだけ、回転だけ) で、モニター PM.* (201 台) の COD の rms・最大値、分散のずれ、チューンのずれを表示する。位置ずれ ae = 40 μm、強さ fe = 4e-4、回転 te = 100 μrad (GCUT = 1σ で打ち切り)。',
  how='`R seed type | Δνx Δνy | COD rms x y [μm] | COD max x y [μm] | η のずれ rms x y [mm]` の行。位置ずれが COD の主因 (rms x 0.2〜0.4 mm、y 0.6〜1 mm)、回転は垂直の COD (rms 50〜100 μm)、強さの誤差は COD に効かずチューン・分散にだけ効く。',
  cmp='△ 種 111・3 は 1e-9 で一致。種 17 の「全部」「位置ずれ」の 2 行は違う (COD rms gs 184 / pysad 301 μm)。この種は垂直の COD が 2 mm 近い大きな場合で、閉軌道の解が前の状態に依り、gs・pysad で別の解に落ちる (乱数列は GCUT = 1 で全桁一致を確認)。',
  limits=['元は 115 個の種 (SeedRandom[187] と Random[1000] で作る)、結果を datas_e.txt に書く。例は 3 本。完全な設定は `seeds` を 115 個に (1 つの種 (4 通り) あたり gs で約 0.3 秒、pysad で約 2 秒)。',
          '位置ずれは `DELX`・`DELY` ではなく成分ごとの `LINE["DX"/"DY", ...]` (K64 の `DELX`・`DELY` は閉軌道に効かない。`modules/modules_errors.n`)。`USE RING NOEXPAND` と `cell calc noexpand` が必要 (`docs/29_sad_quirks.md` A56)。',
          '元の `FFS["MON 0 0 0 0 PM*;"]` は K64 では何もしない (エラーになる) ので除いた。',
          '`GCUT = 0` のとき、gs の `GaussRandom[]` は 1e-9 程度の値を返し、pysad は打ち切りなしのガウス乱数を返す (差の確認のみ。pysad は変えていない)。'],
  show=r'(design|R )', nshow=13),
 'KEKLS_correct_cod': dict(
  kind='誤差入りの COD・分散の補正', orig='correct_cod.sad (correct_cod_best/worst/noSX/GCUT1/v03_70 は設定違い)',
  what='乱数の種 (17・111) の誤差入りのリング RING で、COD 応答行列の特異値分解 (固有値の大きい ICX/2・ICY/2 個の疑似逆行列) により、ステアリング HD (161 台)・VD (121 台) の K0 を 0.5 倍ずつ 6 回かけて、モニター 201 台の COD を補正する。',
  how='各種で `before correction` (チューン、COD の最大・rms、η のずれ)、`iteration l:` (行列で予測した rms と SAD の rms)、`after correction`。種 111 で COD rms (x, y) が 430 / 640 μm から 42 / 23 μm (gs・pysad とも)、SAD の rms は 6 回で数 μm に収束する。',
  cmp='△ 種 111 は補正の各回で gs と pysad が 1e-2 (相対) 以内。種 17 は補正前から違う (位置ずれの閉軌道の解が違う、`KEKLS_errors_cod` と同じ)。',
  limits=['元は種 115 本 (k = 1〜115)、補正 40 回、毎回 DA (5000 周) を計算して datas.txt に書く。例は種 2 本、補正 6 回 (rms は数回で収束)、DA なし (`doDA = 1` で、周回数 200 の DA)。完全な設定は `seeds`・`niter`・`doDA`・`DAturns` (DA の周回数は 5000 のまま) を変える。見込み: 種 1 個あたり pysad で約 10 秒 (応答行列の準備に約 17 秒)、115 個で gs 約 1 分・pysad 約 20 分 (外挿、DA は別)。',
          '補正磁石 HD・VD は同じ要素名の繰り返しなので、既定の `cell calc` では K0 が全部同じ値に揃えられ、個別に効かない (元の台本が前提にしたラティスには、個別の名前があったと思われる)。`USE RING NOEXPAND` と `cell calc noexpand` で成分ごとの K0 が効く (gs・pysad とも確認)。',
          '位置ずれは成分ごとの DX・DY に直接入れる (K64 の `DELX`・`DELY` は閉軌道に効かない)。元の `MON ...` と `s!top;` (打ち間違い) は除いた。',
          '元の modules_correctoptics.n・modules_write_correctdata.n は correct_cod.sad から読み込まれていない (コメントアウトされた Get のみ)。取り込んでいない。',
          '補正後のチューンは設計値から 0.005 動く (補正磁石の K0 の寄与)。チューンの補正は含めていない (元の correct_cod.sad にも無い)。'],
  show=r'(=====|before|after)', nshow=10),
 'KEKLS_cod_response': dict(
  kind='COD 応答行列', orig='write_cod_response.sad',
  what='ステアリング (水平 HD*.*、垂直 VD*.*) の K0 → モニター PM.* の COD の応答行列を、`Response[j,i] = 0.5 sqrt(βm βc) cos(πν − |Δμ|) / sin(πν)` でリング RING の Twiss から作り、水平ステアリング 3 個に K0 = 1e-6 を与えた SAD の COD と比べる。',
  how='`horizontal steerings 161 vertical steerings 121 monitors 201`、`Responsex[...]`・`Responsey[...]` の代表値、`steering HD.n ...` の行 (最大の |COD| と、SAD の COD + K0 × Response の最大の差)。SAD の COD は −K0 × Response と 1e-4 μm (最大振幅 1.7〜3.0 μm の 5e-5) 以内で一致。',
  cmp='○ 最大の COD は 1e-8、解析式との差 (1e-4 μm) の違いは 1e-7 μm 以内で一致。',
  limits=['元の WriteCODResponce[] はモニターを `BPM.*` で探し (v03_68 の BPM は `PM`) 行列が空になり、結果を fort.92 に書く。`PM.*` に直し、書き出しは省いて、代表値と SAD の COD の照合に変えた。',
          '補正磁石 HD・VD は同じ要素名の繰り返しなので、既定の `cell calc` は K0 を全部同じ値に揃える (HD.1 に与えると全 161 個が同じ値になり、COD が 5 倍以上大きくなる)。`USE RING NOEXPAND` と `cell calc noexpand` で 1 個だけに効く (docs/29 A56)。',
          'SAD の K0 の符号は COD = −K0 × Response (PF の例 `pf/optics/calc_respx_PF` と同じ)。'],
  show=r'(horizontal|Responsex|Responsey|steering)', nshow=8),
}

ORDER = ['KEKLS_cell_calc_emit', 'KEKLS_cell_optics', 'KEKLS_match_cell_tune', 'KEKLS_match_chrom', 'KEKLS_match_injcell', 'KEKLS_match_ring_geo', 'KEKLS_match_repeat',
         'KEKLS_scale_mini', 'KEKLS_scale_ncell', 'KEKLS_tune_survey_cell', 'KEKLS_ring_params', 'KEKLS_optics_table', 'KEKLS_calc_H', 'KEKLS_chromaticity',
         'KEKLS_amp_dep_ts', 'KEKLS_sextupole_tune_scan', 'KEKLS_scan_target', 'KEKLS_errors_cod', 'KEKLS_correct_cod', 'KEKLS_cod_response']


def excerpt(name, d, pat, n):
    p = os.path.join(EX, d, name + '.gs.txt')
    rx = re.compile('^' + pat)
    out = []
    for l in open(p, encoding='utf-8', errors='replace'):
        l = l.rstrip('\n')
        if rx.match(l):
            out.append(re.sub(r'(\.\d*?)0{4,}(?=\s|$)', lambda m: m.group(1) if len(m.group(1)) > 1 else '.0', l)[:200])
    if len(out) > n:
        out = out[:n] + ['... (計 %d 行)' % len(out)]
    return '\n'.join(out)


def build(name):
    d = C.SPEC[name][0]
    i = INFO[name]
    r = C.compare(name)
    t = times(name, d)
    title = '# %s' % name
    s = [title]
    s.append('種別: %s。フォルダ: `%s`。元のスクリプト: %s。ラティス: `../lattice/kek_ls_v03_68.sad` (3 GeV、周長 570.72 m、標準セル 20 個)。フラグ: `ON RAD FLUC COD RFSW RADCOD`。' % (i['kind'], d, i['orig']))
    s.append('\n## 何を見る例か\n' + i['what'])
    s.append('\n## スクリプト\n`%s.sad` (先頭に `../lattice/kek_ls_v03_68.sad` を貼った自己完結のファイル。台本の本体は `FFS;` から。説明のコメントは台本の先頭にある)。%s' % (name, ('共通モジュールは `Get["../modules/..."]` で読む (台本のあるフォルダをカレントにして走らせる)。' if 'Get["../modules' in open(os.path.join(EX, d, name + '.sad'), encoding='utf-8').read() else '')))
    s.append('\n## 結果の読み方\n' + i['how'])
    ex = excerpt(name, d, i['show'], i['nshow'])
    s.append('\ngs の出力の抜粋 (指定の行だけ。全体は `%s.gs.txt`):\n```\n%s\n```' % (name, ex))
    for p in i.get('png', []):
        s.append('\n![%s](%s)' % (os.path.basename(p), p))
    s.append('\n## gs と pysad の比較\n判定: **%s** (比べた行 %s、最大の相対差 %s、k64 モードの pysad)。%s' % (
        r['mark'], r['n'], ('%.1e' % r['maxrel']) if r['maxrel'] is not None else '-', re.sub(r'^[○△× ]+', '', i['cmp'])))
    if t:
        s.append('\n所要時間: gs %d 秒、pysad %d 秒 (この PC。k64 モード)。' % (t['gs'], t['py']))
    s.append('\n## 未確認・限界\n' + '\n'.join('- ' + x for x in i['limits']))
    open(os.path.join(EX, d, name + '.md'), 'w', encoding='utf-8', newline='\n').write('\n'.join(s) + '\n')
    return r


if __name__ == '__main__':
    names = sys.argv[1:] or ORDER
    res = {}
    for n in names:
        res[n] = build(n)
        print(n, res[n]['mark'], res[n]['maxrel'])
    json.dump({k: {'mark': v['mark'], 'maxrel': v['maxrel'], 'n': v['n']} for k, v in res.items()}, open(os.path.join(EX, '_marks.json'), 'w'), ensure_ascii=False)
