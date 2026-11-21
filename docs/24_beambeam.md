# BEAMBEAM (強弱ビームのビーム・ビーム相互作用)

最終更新: 2026-11-17 (現状のスナップショット。経緯は `history/` にある)

K64 の `bbstrhl.f` (beambeam, bbinit)、`bbtool.f` (cerrf, gauinv)、`matrix.f` (包絡行列) の直訳。[`pysad/beambeam.py`](../pysad/beambeam.py)。SADHelp に BEAMBEAM の記述はなく、根拠はソースだけである。

## 1. 使い方

```
BEAMBEAM BB = (BX=0.5 BY=0.01 EMIX=1e-8 EMIY=1e-10 SIGZ=0.005 DP=1e-3 SLICE=5 NP=1e10 XANGLE=0.01 ...);
```

追跡 (TrackParticles、`track6n`) の BB の位置で、弱ビーム (追跡する粒子) が強ビームから z スライスごとにキックを受ける。光学 (Twiss、CALC、EMIT) には効かない (恒等。K64 の qtwiss に分岐がない)。

| キーワード | 意味 |
| --- | --- |
| AX BX AY BY, R1..R4, EX EPX EY EPY, ZX ZPX ZY ZPY | 強ビームの Twiss、結合、分散 (IP での値) |
| EMIX EMIY, DP, SIGZ, AZ | エミッタンス、エネルギー広がり、バンチ長、⟨z pz⟩ 相関 |
| R11..R66 | 6 次元の包絡行列 (R11 R22 R33 R44 R55 R66 がすべて正のとき、Twiss より優先) |
| XANGLE | 交差角 [rad] |
| DX DPX DY DPY DZ | 強ビームの軌道 |
| SLICE | z スライス数 (0 以下で無効、最大 500) |
| NP | 強ビームの粒子数 (ce = NP re / γ) |
| BSTRL | ビームストラールング (>0 で有効、§2b) |

弱ビームの γ は MOMENTUM と電子質量から決める (bbinit。陽子でも電子質量を使うのは K64 のまま)。

## 2. 計算の流れ

1. bbinit: 包絡 Σ6 = HRB·Emit·(HRB)ᵀ (Twiss 経路) または R11..R66、交差角の変換 C·Σ·Cᵀ。Σ⁻¹ の (x px y py pz) 部分を A、z との結合を b として、スライスの 5 次元包絡 Σ5 = A⁻¹、重心 v = −Σ5 b (z に比例)。x-y 結合を消す傾き角で回す。
2. スライスの位置: Gauss の分位 (gauinv。6.4e-4 の近似のまま)。
3. 追跡: 交差角・COD の補正、傾きの回転 → 各スライスで x,y を drift (sz=(z−zs)/2)、`cerrf` を使った Bassetti-Erskine のキック (丸い・水平に広い・垂直に広い の 3 通り)、pz のキック (gxy) → 戻す。最後に limitnan (NaN は +1e4、上限 1e4、z は 1e10)、dv を更新。

## 2b. ビームストラールング (BSTRL>0)

各スライスのキック (fx0, fy0) から曲率 ρ = |f|/dsz、放射の確率 rnbs = cnbs·ρ·dsz (cnbs = 5√3 α γ/6)。`tran()` が rnbs より小さいとき 1 個の光子を放射し、`xx/rnbs` を表 (K ベッセル関数 K5/3 の積分の逆引き、`mkbesseltab`) で光子のエネルギー cubs·ρ·exp(ξ) (cubs = 1.5 ħ c γ²/m_e) に直して δ から引く。rnbs>0.5 は K64 が stop するので、pysad も止める ("Nbs is too large")。

- 乱数は pysad の `tran` (SAD と同じ LCG)。粒子 → スライスの順に 1 回ずつ引く。
- 表の定数は単精度のまま (`hbar=6.58211928e-16`、`1.1**40`、`1.02`、`5./3.`、`EPS=1e-10`)。
- PHOTONS が ON なら、光子を PhotonList に入れる (§5。BEAMBEAM の光子だけ)。ルミノシティは `LuminosityMonitor` に渡す (シンボルの値と、定義された関数 `LuminosityMonitor[x_]` の呼び出し。Σ dLum / スライス数 / 粒子数)。

## 3. 検証 (gs との比較)

台本は b1..b12 の 12 本。`TrackParticles[{1, ...}, 10, 1, 5]` で**全周を 5 周** (200 粒子) 追跡して gs と比べる。(`TrackParticles[..., 5, 1]` の 5 は終点の要素番号で、周回数ではない。最初に取った比較はこの誤りで BB を 1 回通るだけだった。全周・複数周で取り直し、同じ結果を得た。)

[`pysad/tests/test_beambeam.py`](../pysad/tests/test_beambeam.py) (24 本) の内訳:

- 相互作用 (b1〜b12): 水平に広い / 垂直に広い / 丸い、SLICE=1,3,4,5、交差角、軌道オフセット、分散と AZ、振幅が大きい場合 (cerrf の漸近の枝)、x=y=0 の粒子 (丸いビーム) は相対 0〜1e-9 で一致。
- 傾き (b6, b9: R3 R4 で x-y 結合) は k64 モードでは K64 の癖 (B34) まで再現: 座標がほぼ 0 (1e-14 程度) に潰れる場合も相対 1e-6 で一致 (v2 モードは回転を直す。§4)。
- **BSTRL** (bs1, bs10, 放射なし; bsr1, bsr10, RAD FLUC あり): `SeedRandom[7]` で 200・300 粒子、5 周、全粒子が乱数列まで gs と一致 (相対 3e-8 以内。光子のエネルギーの表の誤差 1e-9)。表の作り方、逆引き、式が合っている根拠。
- **乱数列が合う粒子数**: 5 周では 330 粒子までは全粒子が一致し、340 粒子以上では一部 (3 割ほど) が一致しなくなる (gs は同じ入力で再現する)。原因は未調査 (粒子が多いと gs が分けて処理し、別の乱数列を使うと推測するが確認していない)。このときは分布で比べた (20000 粒子、5 周、NP=3e11、放射なし): 放射した粒子の割合 gs 0.0664 / pysad 0.0650、平均 δ のずれ gs −2.14e-8 / pysad −2.21e-8、δ のずれの分布の KS 検定 p=1.00。(以前の 1 周ぶん・NP=2e12 でも、割合 0.1458 / 0.1455、KS p=0.07〜0.54 で統計誤差の範囲だった。) 個々の粒子の一致は求めない (乱数列が違うため)。
- **LuminosityMonitor** (lm2, lm3): シンボルの値と、利用者が定義した関数 `LuminosityMonitor[x_]` の呼び出し (回数と値) が gs と一致 (相対 1e-9)。
- **PHOTONS** (PhotonList): 7 個の光子の位置・方向・偏光・粒子番号・要素番号が一致、エネルギーは相対 3e-9。
- 乱数の状態: RAD FLUC の追跡のあと `Random[]` が gs と同じ値 (放射の分割が一致する根拠、§4)。
- 光学 (Twiss): gs も BEAMBEAM を恒等 (MARK と同じ) として扱い、pysad の Twiss は gs と 10 桁一致 (BEND を含む 9 行の BX BY NX EX DX で確認)。
- 光学の Emittance (2026-11-16 夜): BEAMBEAM を含むリングの CALC・Emittance[] は gs と Twiss のチューン 1e-14、運動量収縮 4e-12 で一致 (BEAMBEAM は Marker 扱い)。1 つのリングでの照合。
- キーワード STURN・COUPLE (BEAMBEAM): K64 のソースでキーワードの登録 (initb1.f) 以外に参照がなく、使われない (ソースの確認。gs での動作は試していない)。
- PhotonList の 12 列は `photons.f` の `tphotonlist` の定義と一致 (1 エネルギー=|p|·amass·gammab、2-4 位置、5-7 方向、8-10 偏光、11 粒子番号、12 要素番号。ソースとの照合)。
- wake: BEAMBEAM に WakeFunction を付けた台本 (縦と横、5 周) が gs と相対 0 で一致。

## 4. K64 の癖と、見つけた pysad の差

- **傾き回転** (B34): `sint=xangle(5)` (傾き角そのもの)、`cost=xangle(6)` (sin)。bbinit は (6)=sin、(7)=cos を書くので、ずれている。強ビームに x-y 結合があると、回転が直交でなくなり粒子の座標が潰れる。**k64 モードはこれを再現し、v2 モードは (cos t, sin t) の正しい回転にする**。強ビームに x-y 結合 (R1..R4 か R13 など) があるときだけ結果が変わる。根拠: 回転に対する共変性 (包絡を t 回した強ビームの蹴り = 粒子を −t 回し、回していない包絡で蹴り、+t 回す) が t = 0.3 で 1e-9 (蹴りに対する相対)。試験 [`test_v2_C_beambeam.py`](../pysad/tests/test_v2_C_beambeam.py) (2 件)。gs との比較は無い (K64 は座標が潰れる)。確かめたのは t = 0.3 の 1 例で、実機の衝突点の規模では未確認。
- `cerrf` の `3.2` は単精度の定数。
- 丸いビームで x=y=0 の粒子は 0/0 → NaN → px=+1e4。
- スライスの境界 (`dsz`) はビームストラールングだけが使う。
- **TRPT の放射の分割** : gs の `ndivrad` (tffs.f) は TRPT のとき `int(1+h0/100·aka·10)`、そうでなければ `int(1+h0/100·aka)`。pysad は `track6n.turn` で TRPT のとき ndivrad に渡す h0 を 10 倍にして同じ分割数にする (四極 1 分割に対し gs は 4。10 倍しないと乱数を引く数が gs と違い、放射の効果が小さくてもあとの乱数列がずれる)。

## 5. 未対応・限界

- **PHOTONS (PhotonList) は BEAMBEAM の光子だけ**。他の要素 (四極・BEND・MULT など) が放射する光子 (tphrec) は含まれず、その乱数 (光子ごとに `tgauss` を 2 回) も引かないので、RAD FLUC で他の要素が光子を出すと、そのあとの乱数列が gs とずれる (1 周の台本 ph1 は他の光子が出ず一致、3 周の台本では四極の光子が出て 43 個 (gs) と 12 個 (pysad) になる)。BEND/ANGLE≠0 の MULT を含む並びは SADError にする。PHOTONS は RAD と FLUC が必要 (gs と同じメッセージ)。
- 光子の放射角 (`radangle1`) は BEAMBEAM では曲率半径 0 (直線) なので、角度の広がりなし・偏光 (0, 0, 1) を、幾何の傾き chi で回した値になる。`tgauss` の Box–Muller のキャッシュが前から残っている場合 (利用者が GaussRandom を奇数回呼んだ後など) は未対応 (2 回引く前提)。
- K64 は光子ごとに `write(*,*) x,px,y,py,delg` を標準出力に出す。pysad も同じ形式で出す。
- **BEAMBEAM がソレノイド領域 (SOL BOUND=1 の間) にある配置**: 元の SAD K64 が SIGSEGV で落ちる (tsol.f:798 付近。四極と BEAMBEAM 以外の SEXT や COORD でも同じ。診断メッセージの書き込みではなく、メッセージを消しても落ちる。原因は未特定、docs/21 B35)。pysad は NotImplementedError にする。
- BBCUT (旧 bb.f の分岐のみ、bbstrhl.f は使わない)。
- 加速のある CAVI の後で MOMENTUM が変わる場合: bbinit は更新時の MOMENTUM を使う (pysad は Prog 作成時の運動量)。
- track6 (低速の追跡) は BEAMBEAM を扱わず NotImplementedError。

解説: docs/39_beambeam_guide.md (BEAM-BEAM)
