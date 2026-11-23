# 11_sad_quirks.md から取り除いた履歴の記述

11_sad_quirks.md (現在の仕様の解説) から、版・経過に結びついた記述を取り除いた。下は取り除く前の元の行を、節の見出しつきでそのまま並べたもの (行番号は整理前のもの)。

## # 11. SAD の癖・仕様・ビルド依存の挙動 (pysad での扱い)

整理前の 3 行目:

````text
最終更新: 2026-11-18 (v2.1.0 の注記を足した)
````

整理前の 7 行目:

````text
> pysad は「再導出せず直訳」が方針なので、下記の癖は (バグに見えるものも) そのまま再現している (**k64 モード**。v2.0.0 から既定の v2 モードは、K64 のバグと分かったものの一部を直した。直したものには各行に「**v2 で修正**」と B 番号を付けた。モードによらない pysad 自身の不具合の修正は docs/47 の v2.0.0 の (d) と v2.1.0 の (c)。一覧は docs/47。v2.1.0 で直したものも同じ書き方で注記した)。
````

## ## A. 数値計算の癖 (結果の値に効くもの)

整理前の 44 行目:

````text
| A32 | **RFSW OFF でも RADCOD ON なら tcod の rfsw=.false. の枝 (im=5) が放射込みの閉軌道を求める**: x,px,y,py を閉じ、z のずれは trf0 (Equilibrium dz) に移る。出口の δ は −U0/P0 (RF がないので閉じない)。U0 (実ラティスで 2.96054e5 vs 2.96075e5 eV)、α_c、チューンが RADCOD OFF と 1e-5 違う | tcod.f, bench_level26a R3 | `emit6d.closed_orbit6d` の新しい枝 (`_tcod_radcod(rfsw=False)`)、`emittance6d` の最後の一周も radcod=True。**v2.1.0 (D1)**: 空洞があって RFSW OFF のとき、pysad は z の戻りを δ で閉じていた (K64 の tcod は x・px・y・py だけを閉じ、δ も z も動かさない。z は trf0 に移る) 不具合を直した (両モード)。軌道のずれ (四極の DX は経路長を変える) があると εx が gs と 1.8〜4% 違っていた。空洞のないリングは従来どおり δ で閉じる (回帰 3c・3d・26b)。tradke と tradkn の反作用の角度の符号が逆なのは K64 の誤り (B111) |
````

整理前の 58 行目:

````text
| A46 | **POL のスピン追跡は RAD ON のときだけスピンが動く** (スピンの歳差が `photons.f` の tradk の中の `sprot` で計算されるため。RAD OFF では POL でもスピンは変わらない)。TrackParticles の入力は 9〜11 列 {x,px,y,py,z,dp, 角度, sy/\|s\|, flag}、戻りの 7・8 列は atan(sz, sx) と sy/\|s\| (tftrack.f)。損失 (flag) の処理は高速追跡にしかなく、POL のスピン追跡では粒子が落ちない (pysad の制限) | photons.f sprot、tftrack.f の ktfresetparticles。gs と 4e-12 以下 ([`test_spin_track.py`](../pysad/tests/test_spin_track.py)) | `track6.py` の `sprot`・`bsi`、`fn_track` が POL かつ RAD ON で純 Python の追跡に切り替える (k64 モードは損失の処理なし)。**v2 で変更**: 損失の処理を足した (LOSSAMPL・LOSSMAP は両モードで gs と一致)。歳差の角を厳密な θ に (B83)、SOL 領域内・UND の縦磁場 (B81・B82) を直した (docs/43 §8a) |
````

整理前の 59 行目:

````text
| A47 | **`Polarization` / `Polarizations` (平衡偏極度) は K64 では NaN** (spdepol の `exp(cb²)·cerfc(cb)` が \|cb\|>18 で 0·∞ になる)。PolarizationTime・SpinTune・PolarizationVector は有限で gs と一致する | photons.f:spintrm。docs/21 B55 | K64 の式のまま (NaN も写す)。**v2 で修正** (B55): 有限値 (平面リングで 92.376%。スピン軌道結合の脱偏極率から。v2.1.0 でシンクロトロン側帯の項を足した。docs/43 §8a)。**v2.1.0 で誤りの位置を特定** (B125): 桁あふれを直しても K64 の率は emittance/減衰率に比例し、独立な率の約 1000 倍 (4.2 GeV・12 MV の 12 セル FODO)。2 次の零位相の項の二重数え + `spintrm` のガウスの項の 2π が余分。逆向きの BEND (挿入光源の極) では緩和率が符号つきで相殺する (B116。v2 は Σ r (b·n)²) |
````

整理前の 61 行目:

````text
| A49 | **クラスのメンバーでない名前の `obj@name[...]` は式のまま残る** (メッセージなし)。MessageName の関数形 `MessageName[f,"usage"]` は gs が f を無視して `Undefined MessageName General::"usage"` を出す (`f::usage` の形は一致)。pysad は関数形でも `f::usage` と同じ表を引くので差があり、見送り | gs の出力、work_unimpl.md | `sadclass.py` (メンバーでなければ None で式のまま)。**v2 で修正** (B70): `MessageName[f,"usage"]` の関数形は f の本文を返す (k64 モードの pysad も以前から f を見る) |
````

整理前の 71 行目:

````text
| A59 | **BEND (K1≠0) の K1 は、中央面の磁場 B_y = B0 + (K1/h) ln(r/ρ0) (r = ρ0 + x)、すなわち曲線座標の勾配 (1 + hx)∂b/∂x が一様な定義** (B''(ρ0) = −h B'(ρ0)、x y² の項なし)。B' が x に一様な定義 (elegant の CSBEND の K2 = 0・MAD の SBEND) とは、同じ K1 でも 6 極 K2 = +hK1 の差がある。**SAD → elegant・MAD の変換には本体に K2 = −hK1 を足す**。δ = 0 の光学・エミッタンスには効かず、δ≠0 のチューン (ξ・∂²ν/∂δ²)・LMA・DA・クロマティック ADTS に効く (n = 0.4 の扇形磁石 16 台、3 GeV の dν/dδ は SAD (−0.17213、+0.21082)、B' 一様 (−0.60246、+0.73786)) | `tbend.f` の tbendiinit (41-116 行)・tbendicorr (118-134 行。132 行の kick は (1 + u) ln(1 + u) − u の u³ までの級数)・tbendibody (136-160 行)、`tbendi.f` 89-158 行、CALC は `tbende.f` の tbendecorr (61-121 行)。SADHelp の BEND の K1 は B'L/Bρ とだけ書き B'' を決めていない | pysad も同じ (k64・v2 とも。追跡 `track6.py` `_tbendicorr`・numba `track6n.py` `k_bendi`、CALC `perfkern_bend.py` `bendecorr`)。**v2 のオプション (v2.3.0 U2): `compat.set_bend_k1_def('lin')` (環境変数 `PYSAD_BEND_K1_DEF=lin`) で B' 一様の定義 (= 'log' + 6 極 K2 = +hK1 を本体に足したもの、kick のポテンシャル K1 h (x³/3 − x y²/2)) にする。既定は 'log' (K64。k64・v2 とも)。k64 モードでは無視する。追跡 (`track6.py`・`track6n.py` `k_bendi`)・CALC・`Emittance[]` (`perfkern_bend.py` `bendecorr`)・MULT (ANGLE≠0) の K1 に効く。δ = 0 の設計軌道の光学・行列は変わらない。n = 0.4 の扇形磁石 16 台の dν/dδ は 'lin' で (−0.60246、+0.73786) (解析式・独立な積分・elegant K2 = 0 と一致)。条件・限界は `history/v23_U2.md`**。B160 (docs/21)、`history/v23_U1.md`、`history/v23_U2.md`、[`test_v23_U1_bend_k1_field.py`](../pysad/tests/test_v23_U1_bend_k1_field.py)・[`test_v23_U2_bend_k1_def.py`](../pysad/tests/test_v23_U2_bend_k1_def.py) |
````

## ## B. 乱数と並列 (再現には順序が重要)

整理前の 82 行目:

````text
| B5 | **NPARA の既定は Packages/init.n で 4** (initb1.f の 1 ではない)。TrackParticles は粒子 9 個以上 (npparamin=9) で fork し乱数列が分かれる。NPARA>1 だと子プロセスが出力を重複して書く | Packages/init.n:219、tftrack.f。Level 16 | 参照値は `NPARA=1;` を入れて取る。pysad の既定は `NPARA` 1。**v2.1.0 で `NPARA=n;` を FFS で設定すると gs の分け方を再現する** (子の乱数の種 = 親の種 + 2k、np ≥ 13 で分かれる。有効な並列数 = min(NPARA, np/9 + 1, ne·np/3000 + 1)。RAD + FLUC の追跡だけ)。「δ の rms が gs 4.21e-4、pysad 4.52e-4」の差の原因はこれ (NPARA=4 と 1 の比較)。NPARA=1 と 4 とも gs と全桁一致 (FODO・PF リング) |
````

## ## D. ビルド依存・環境の挙動

整理前の 123 行目:

````text
| D1 | **WSL ビルドは FITP>1 の GO で再起動する** ("*** Run time Environment" からスクリプト再実行、NPARA=1 でも)。参照値が取れない | bench_level15a | match6/ffscmd の FITP>1 (idp ループ、重み、INTRES 判定、fam 点) は直訳のまま SAD と未照合 (自己検証のみ)。v2.1.0 の確認: 今の gs は再起動せず完走し最小値の手前で止まる。pysad の最終解は独立な最小化と 5e-8 で一致 (docs/03 §9) |
````

## ## E. 未解決・要確認

整理前の 134 行目:

````text
  1. (解決: BEND の `F1 FRINGE` + RAD の追跡の差は B113。tfgeo の inext・iprev の再現。v2.1.0)
````

整理前の 135 行目:

````text
  2. (解決: L≠0 の CAVI の DX・DY の差は CAVI の TENT・TEXIT が pysad に無かったため。v2.1.0)
````

整理前の 136 行目:

````text
  3. SOL (BZ≠0) 内の四極 + RAD の追跡で、現在の gs (2026-10-19 ビルド) が保存済みの参照 (bench_level9、10-14 の gs の出力) と別の値を返す (x が -0.0027867898 から -0.0014539780)。pysad は保存済みの参照 (と tsolqur の癖 A1) に合わせている。2026-11-16 の単一粒子 1 周追跡 (BZ を 0→0.01 T) で、領域内がドリフトだけなら gs=pysad、`SB D1 QD SE` の gs は `SB D1 SE` (四極なし) と 8 桁一致するので、**現在の gs が RAD ON・BZ≠0 の領域内の四極の力を入れていない疑いが強い** (B68。10-19 と 11-14 -O0 の 2 ビルドは同値で最適化の差ではない。ソースの箇所は未特定。NORAD は両方で一致)。
````

整理前の 137 行目:

````text
  4. SOL 領域内の BEND・MULT・CAVI は k64 モードでは未対応 (v2 モードは実装した: docs/44 §7a)。BEND は K64 の tsol.f どおり (`tdrift(al,bzs,phiy,phix,…)`) に実装したが、gs が偏向をほぼ無視して (ANGLE=0.1 でも 1e-5) 合わず k64 モードの実装は取り下げた。原因は未確認。
````

整理前の 138 行目:

````text
  5. K3 が強くずれの大きい MULT リング (K3=300、ずれ 1e-4) で pysad の CALC (qcod) の閉軌道が gs と違う解 (y=-0.033) になり Emittance が NaN。gs は小さい軌道で通る。
````

整理前の 139 行目:

````text
  6. (解決しなかった・再現しない: `emit6d._tquade` の解析ヤコビアンの δ の微分の欠落は再現せず、Python 版・numba 版・数値微分が 1e-12 で一致した。v2.1.0。z の行の関数単体の欠落を補った)
````

整理前の 140 行目:

````text
  7. (解決: EMIT のあとの `USE` で RF の位相の基準 `trf0` が持ち越されるのは K64 の仕様。pysad が USE で 0 に戻していた。v2.1.0)
````

整理前の 141 行目:

````text
  8. `Check[Return[1],2]`・`Thread[a==b]` (a が値、b が未定義)・`Out[n]` の繰り返し表示など SADScript の細部 (docs/04 §6)。
````

整理前の 142 行目:

````text
- 未翻訳のまま (理由つき): MAIN の TRACK の Fourier 解析 (wjfit・wfres・wfbin・wiord・wjdraw)、(光学・Emittance への空間電荷 (WSPAC) は v2.1.0 で実装済み。追跡の WSPAC のビーム行列は空間電荷なしの平衡値のまま)。MEASURE/TRCOD (gs が動かない: B42)・`tspac.f` の追跡側 (空間電荷)・PHSROT (K64 で動作しない: B38) は v2 モードで実装した (k64 モードは未対応)。一覧は `46_unimplemented_and_unverified.md` と `16_fortran_status.md` 末尾。
````

