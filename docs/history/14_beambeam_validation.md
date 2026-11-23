# 14_beambeam の検証の詳細 (gs との比較)

docs/14_beambeam.md の「検証」の詳細。旧 docs/39_beambeam_guide.md §10 と旧 docs/24_beambeam.md §3 を元の文脈つきで移した。

## 旧 39 §10 gs との一致の度合い


| 項目 | 一致 | 備考 |
| --- | --- | --- |
| 相互作用 (b1〜b12: 水平に広い・垂直に広い・丸い、SLICE=1,3,4,5、交差角、軌道オフセット、分散と AZ、大振幅) | 相対 0〜2e-9 | 200 粒子、全周 5 周、1・77・200 番目と成分和 |
| 傾き (b6, b9: R3 R4 の x-y 結合) | 相対 1e-6 (最大値比) | k64 モード: K64 の癖 B34 まで再現。座標がほぼ 0 に潰れる (v2 モードは回転を正しくする。§11) |
| BSTRL (bs1, bs10、RAD FLUC の bsr1, bsr10) | 相対 3e-8 | SeedRandom[7]、200・300 粒子、乱数列まで一致 |
| 乱数が合う粒子数 | 5 周・330 粒子以下 | 340 粒子以上では一部 (3 割ほど) の粒子が個別には合わない。原因は未調査。分布で比較 (下) |
| LuminosityMonitor (lm2, lm3) | 相対 1e-9 | シンボルと定義関数の呼び出し |
| PHOTONS (ph4) | 位置・方向・偏光・番号は一致、エネルギーは相対 3e-9 | 7 個の光子 |
| RAD FLUC 後の `Random[]` | 一致 | TRPT の ndivrad を 10 倍にしているため |
| Twiss | 10 桁 | gs も恒等 |
| wake を併用 (縦・横、5 周) | 相対 0 | |
| BBBrem1 | 相対 1e-15 | 乱数列も |

340 粒子以上の分布比較 (20000 粒子、5 周、NP=3e11、放射なし): 放射した粒子の割合 gs 0.0664 / pysad 0.0650、平均 δ のずれ gs −2.14e-8 / pysad −2.21e-8、δ のずれの分布の KS 検定 p=1.00。1 周ぶん・NP=2e12 でも割合 0.1458 / 0.1455、KS p=0.07〜0.54 で統計誤差の範囲。個々の粒子は一致しない。

試験は [`pysad/tests/test_beambeam.py`](../pysad/tests/test_beambeam.py) (parametrize を含め BEAMBEAM 関連)、[`test_bbbrem.py`](../pysad/tests/test_bbbrem.py)、[`test_measure.py`](../pysad/tests/test_measure.py)。台本は b1..b12 の 12 本。


## 旧 24 §3 検証 (gs との比較)


台本は b1..b12 の 12 本。`TrackParticles[{1, ...}, 10, 1, 5]` で**全周を 5 周** (200 粒子) 追跡して gs と比べる。(`TrackParticles[..., 5, 1]` の 5 は終点の要素番号で、周回数ではない。最初に取った比較はこの誤りで BB を 1 回通るだけだった。全周・複数周で取り直し、同じ結果を得た。)

[`pysad/tests/test_beambeam.py`](../pysad/tests/test_beambeam.py) (24 本) の内訳:

- 相互作用 (b1〜b12): 水平に広い / 垂直に広い / 丸い、SLICE=1,3,4,5、交差角、軌道オフセット、分散と AZ、振幅が大きい場合 (cerrf の漸近の枝)、x=y=0 の粒子 (丸いビーム) は相対 0〜1e-9 で一致。
- 傾き (b6, b9: R3 R4 で x-y 結合) は k64 モードでは K64 の癖 (B34) まで再現: 座標がほぼ 0 (1e-14 程度) に潰れる場合も相対 1e-6 で一致 (v2 モードは回転を直す。§4)。
- **BSTRL** (bs1, bs10, 放射なし; bsr1, bsr10, RAD FLUC あり): `SeedRandom[7]` で 200・300 粒子、5 周、全粒子が乱数列まで gs と一致 (相対 3e-8 以内。光子のエネルギーの表の誤差 1e-9)。表の作り方、逆引き、式が合っている根拠。
- **乱数列が合う粒子数**: 5 周では 330 粒子までは全粒子が一致し、340 粒子以上では一部 (3 割ほど) が一致しなくなる (gs は同じ入力で再現する)。原因は未調査 (粒子が多いと gs が分けて処理し、別の乱数列を使うと推測するが確認していない)。このときは分布で比べた (20000 粒子、5 周、NP=3e11、放射なし): 放射した粒子の割合 gs 0.0664 / pysad 0.0650、平均 δ のずれ gs −2.14e-8 / pysad −2.21e-8、δ のずれの分布の KS 検定 p=1.00。(1 周ぶん・NP=2e12 でも、割合 0.1458 / 0.1455、KS p=0.07〜0.54 で統計誤差の範囲だった。) 個々の粒子の一致は求めない (乱数列が違うため)。
- **LuminosityMonitor** (lm2, lm3): シンボルの値と、利用者が定義した関数 `LuminosityMonitor[x_]` の呼び出し (回数と値) が gs と一致 (相対 1e-9)。
- **PHOTONS** (PhotonList): 7 個の光子の位置・方向・偏光・粒子番号・要素番号が一致、エネルギーは相対 3e-9。
- 乱数の状態: RAD FLUC の追跡のあと `Random[]` が gs と同じ値 (放射の分割が一致する根拠、§4)。
- 光学 (Twiss): gs も BEAMBEAM を恒等 (MARK と同じ) として扱い、pysad の Twiss は gs と 10 桁一致 (BEND を含む 9 行の BX BY NX EX DX で確認)。
- 光学の Emittance (2026-11-16 夜): BEAMBEAM を含むリングの CALC・Emittance[] は gs と Twiss のチューン 1e-14、運動量収縮 4e-12 で一致 (BEAMBEAM は Marker 扱い)。1 つのリングでの照合。
- キーワード STURN・COUPLE (BEAMBEAM): K64 のソースでキーワードの登録 (initb1.f) 以外に参照がなく、使われない (ソースの確認。gs での動作は試していない)。
- PhotonList の 12 列は `photons.f` の `tphotonlist` の定義と一致 (1 エネルギー=|p|·amass·gammab、2-4 位置、5-7 方向、8-10 偏光、11 粒子番号、12 要素番号。ソースとの照合)。
- wake: BEAMBEAM に WakeFunction を付けた台本 (縦と横、5 周) が gs と相対 0 で一致。


## 旧 32 BBBrem1 の検証 (旧 docs/32_gaussiancoulomb_bbbrem.md)

検証 (`test_bbbrem.py`): `SeedRandom[3]` のあと 6 事象、sqrt(s)=10,20,30,40 GeV の 4 事象で、乱数の列も含めて相対 1e-15 (最後の桁の丸め) で一致 (試験の許容は相対 1e-11)。
