# example/05_tracking の md から取り除いた履歴の記述

例の md を「現在の動作と結果の説明」にするため、版・経過・以前の状態の記述を取り除いた。取り除いた元の行を、ファイルごと・元の行番号つきで残す (行番号は整理前のもの)。この節は、整理前の md (リポジトリの先頭の版) と現在の md の差分から再構成した。

## README.md

元の 12 行目:
> | [elem_APERT_LOSSAMPL](elem_APERT_LOSSAMPL.md) | APERT 要素と LOSSAMPL | ○ (APERT 要素・LOSSAMPL とも一致。SOL 領域内の APERT は v2 で対応、gs と未比較) |

元の 18 行目:
> | [flag_SPAC](flag_SPAC.md) | 空間電荷 (1 m のドリフトの広がりとリングのチューンシフトを解析式と比べる。RADIUS<0 は絞りだけ) | △ (RADIUS<0 の絞りは一致。RADIUS>0 の空間電荷は v2 で実装し解析式に合うが、gs (K64) は力が 1/γ 倍に弱いので意図して違う。k64 モードは空間電荷も絞りも入らず gs と合わない: 2026-11-17 確認) |

元の 19 行目:
> | [flag_LOSSMAP](flag_LOSSMAP.md) | 損失の要素・周を戻りの 8・9 列目に | ○ (v2 で対応) |

元の 25 行目:
> - `FFS` の中で `MOMENTUM = ...` を設定すると、gs と同じく追跡・Emittance[] に反映され、`USE` で MAIN 側の値に戻らない (2026-11-17 に対応。以前の pysad は USE で戻した)。

## elem_APERT_LOSSAMPL.md

元の 34 行目:
> (A)(B)(C) とも gs と一致 (○)。pysad の追跡 (numba 版 `track6n.py` の `k_apert_el`、Python 版 `track6.py` の `_apert_pass`) は、APERT 要素について `tapert.f` / `tapert1` の判定をそのまま行う。LOSSAMPL ((B)、(C)) は従来から一致している。

元の 39 行目:
> - ソレノイド領域 (SOL BOUND≠0 の間) の中にある APERT: この例の最初の版では pysad が拒否していた (gs の tsol は tapert1 を呼ぶ)。v2 で追跡に対応した (`docs/47_v2_changes.md` の B69・B80。k64 モードは NotImplementedError)。この例のリング (SOL なし) では使っておらず、gs との比較はしていない。領域内の APERT のスピン追跡は未対応。

## flag_LOSSMAP.md

元の 32 行目:
> **○ 一致 (2026-11-17、v2.0.0 の再実行)。** pysad の戻りは `LOSSMAP` で 9 列になり、flag {1,1,0,0,0}、損失の要素 {0,0,242,242,242}、損失の周 {31,31,4,1,1}、`LINE["NAME", ...]` の名前 {"-","-","D1.24","D1.24","D1.24"} が gs と同じ (`NOLOSSMAP` の 7 列も同じ)。この例の最初の版では、pysad は LOSSMAP を ON にしても 7 列のままで `r[[2,8]]` が `Part: index 8 out of range (length 7)` になっていた (v1 の未対応)。v2 で実装した (`docs/47_v2_changes.md`。互換モードによらない、gs と同じ結果の機能)。

## flag_RAD_FLUC.md

元の 46 行目:
> **v2 の注記 (2026-11-17、v2.0.0 の再実行)**: `.pysad.txt` は既定の v2 モードの出力に作り直した。v2 は `Emittance[]`・EMIT・追跡の Jx の端の項 (B75) と `MomentumCompaction` (B77) を直したので、この例の pysad の値のうち εx (約 2.4% 大きい)・εz・減衰率の X と Z・σδ・バンチ長・α・バケット高さ・放射によるチューンシフトが gs と違う (K64 のバグを直したため。`docs/47_v2_changes.md`)。以下の表と本文の「pysad」の値と「一致」の判定は、v1 と同じ `PYSAD_COMPAT=k64` での結果 (旧い出力と同じことを 2026-11-17 に確認)。gs と同じ値が要るときは `PYSAD_COMPAT=k64` を付ける。

## flag_SPAC.md

元の 27 行目:
> - **フラグは最上位に書く。** `Do[..., NOSPAC; ...; SPAC; ...]` のように式の中に書くと、フラグの命令にならず切り替わらない (この例の以前の版の誤り。そのため「gs でも SPAC の効果が出ない」と書いていた)。式の中で切り替えるなら `FFS["SPAC"]`。

元の 29 行目:
> - 粒子は乱数を使わない決定的な分布。(A) は半径 1 mm の一様な円 × 長さ 0.1 m の一様なバンチ (2000 個、角度 0)、(B) は整合したガウス分布 (εx = εy = 1e-6 m、σz = 1 cm、4000 個)、(C) は以前の版と同じ 200 個。

元の 30 行目:
> - `PBUNCH = 1e13` (以前の版) は桁外れ (チューンシフトが 100 を超える) なので、効果が見えて破綻しない値にした: (A) 1e12 (20 MeV で半径が 1.44 倍)、(B) 1.5e9 (Δν ≈ −0.01)。

元の 31 行目:
> - `SPACSTEP` (キックの間隔 [m]) は pysad v2 で足した変数 (既定 0.05)。gs では普通の変数で、意味を持たない (5 cm 固定)。

元の 77 行目:
> - **gs は SPAC の粒子数が少ないと落ちる**: 以前の版に書いた「小さい RADIUS>0 (0.4 mm) で gs が不正終了する」の原因。`spkick.f` の `spmesh` は縦を 10 枚のスライスに分け、最後のスライスが空 (粒子数 np が 9·ceil(np/10)+1 より小さい) だと配列の外を読み、1 粒子だと厚さ 0 で割る。gs で確かめた: np = 2、15、33 は SIGSEGV、5 は Bus error、9・10・11・19・25・28・81・82 は NaN で全粒子が失われ、20・30・40・…・100・200 と 1 は正常。小さい絞りで粒子が減っていくと、これに当たる。pysad v2 は粒子数によらず走る (試験 `test_spac_small_particle_number`)。

元の 84 行目:
> - 以前の版にあった「別件」(FFS の中の `MOMENTUM = …` の再設定を pysad が追跡に反映しなかった) は 2026-11-17 に直っている。この例は運動量をラティスの先頭 (MAIN レベル) に置いている。

## fn_DynamicApertureSurvey_flag_DAPERT.md

元の 36 行目:
> **v2 の注記 (2026-11-17、v2.0.0 の再実行)**: `.pysad.txt` は既定の v2 モードの出力に作り直した。v2 は `Emittance[]`・EMIT・追跡の Jx の端の項 (B75) と `MomentumCompaction` (B77) を直したので、この例の pysad の値のうち εx (約 2.4% 大きい)・εz・減衰率の X と Z・σδ・バンチ長・α・バケット高さ・放射によるチューンシフトが gs と違う (K64 のバグを直したため。`docs/47_v2_changes.md`)。以下の表と本文の「pysad」の値と「一致」の判定は、v1 と同じ `PYSAD_COMPAT=k64` での結果 (旧い出力と同じことを 2026-11-17 に確認)。gs と同じ値が要るときは `PYSAD_COMPAT=k64` を付ける。

## fn_SurvivedParticles.md

元の 26 行目:
> - 損失の位置と周は戻りに入らない (LOSSMAP を使う: [flag_LOSSMAP](flag_LOSSMAP.md)。pysad も v2 で LOSSMAP に対応し、gs と同じ値が出る)。

## multi-turn-tracking.md

元の 45 行目:
> **v2 の注記 (2026-11-17、v2.0.0 の再実行)**: `.pysad.txt` は既定の v2 モードの出力に作り直した。v2 は `Emittance[]`・EMIT・追跡の Jx の端の項 (B75) と `MomentumCompaction` (B77) を直したので、この例の pysad の値のうち εx (約 2.4% 大きい)・εz・減衰率の X と Z・σδ・バンチ長・α・バケット高さ・放射によるチューンシフトが gs と違う (K64 のバグを直したため。`docs/47_v2_changes.md`)。以下の表と本文の「pysad」の値と「一致」の判定は、v1 と同じ `PYSAD_COMPAT=k64` での結果 (旧い出力と同じことを 2026-11-17 に確認)。gs と同じ値が要るときは `PYSAD_COMPAT=k64` を付ける。


# 台本 (`.sad`・`.py`) の冒頭コメントから、機能別フォルダへの移動時に取り除いた履歴の記述

v2_next・v2_C・v2_B の作業フォルダから機能別のフォルダへ移したとき、台本 (`.sad`・`.py`) の冒頭コメントにあった作業符号 (「v2next A1 項目 4」など) と、旧版の動作を語る文 (「旧 pysad は…」など) を取り除いた。取り除いた行を、元の行と直したあとの行の対で残す。

- `flag_PSPAC_resolution.sad` 1 行目:
  - 元: `! v2next A1 項目 8: PSPAC (SPCH 要素) のメッシュの分解能。gs (K64) との差の原因はメッシュ (縦: λ(z) の入れ方、横: 電位の補間)。`
  - 後: `! PSPAC (SPCH 要素) のメッシュの分解能。gs (K64) との差の原因はメッシュ (縦: λ(z) の入れ方、横: 電位の補間)。`
- `flag_WSPAC_emittance.sad` 1 行目:
  - 元: `! v2next B1 項目 3: 空間電荷 (WSPAC) の Emittance[]・EMIT・CALC6D・CALC。K64 の tintrb の twspace (ビーム行列 σ のガウス分布の線形レンズ`
  - 後: `! 空間電荷 (WSPAC) の Emittance[]・EMIT・CALC6D・CALC。K64 の tintrb の twspace (ビーム行列 σ のガウス分布の線形レンズ`
- `fn_NAFF_tunes.py` 2 行目:
  - 元: `"""v2next B3 項目 5: 追跡の座標列からチューンを取り出す NAFF (pysad/naff.py)。小さなリングと大きなリングで CALC のチューンと比べ、`
  - 後: `"""追跡の座標列からチューンを取り出す NAFF (pysad/naff.py)。小さなリングと大きなリングで CALC のチューンと比べ、`
- `flag_FLUC_nparallel.sad` 1 行目:
  - 元: `! v2next C3 項目 7: FLUC つきの追跡の乱数列。gs の既定は NPARA=4 (tftrack.f: 粒子を子プロセスに分け、子の乱数の種は 親の種 + 2k)。`
  - 後: `! FLUC つきの追跡の乱数列。gs の既定は NPARA=4 (tftrack.f: 粒子を子プロセスに分け、子の乱数の種は 親の種 + 2k)。`
- `fn_DynamicApertureSurvey_dapwidth_window.sad` 1 行目:
  - 元: `! v2next P4 項目 5: DynamicApertureSurvey の DAPWIDTH の窓 (trackd.f の粒子スロットの直訳)。gs (K64) と: NORAD の地図 (DAPWIDTH = 7 (既定)・3・12・51) と Score が一致`
  - 後: `! DynamicApertureSurvey の DAPWIDTH の窓 (trackd.f の粒子スロットの直訳)。gs (K64) と: NORAD の地図 (DAPWIDTH = 7 (既定)・3・12・51) と Score が一致`
- `flag_PSPAC_astrid2.sad` 1 行目:
  - 元: `! v2next P4 項目 3 (大きなリングでの確認): ASTRID2 (0.58 GeV、周長 45.7 m、ASTRID2_v01_01.sad) に SPCH を 3 か所 (始点・RF の前・終点) 置き、PBUNCH = 1e14 (非現実的に大きい。K64 の PSPAC の数値の一致を見るため)、`
  - 後: `! (大きなリングでの確認): ASTRID2 (0.58 GeV、周長 45.7 m、ASTRID2_v01_01.sad) に SPCH を 3 か所 (始点・RF の前・終点) 置き、PBUNCH = 1e14 (非現実的に大きい。K64 の PSPAC の数値の一致を見るため)、`
- `flag_PSPAC_k64_shape.sad` 1 行目:
  - 元: `! v2next P4 項目 3: k64 モードの PSPAC (K64 の spch.c の形状関数を直訳: 横 CIC・縦「ほぼ最近接」・双 3 次の電位・3 次スプラインの線密度)。gs と全桁一致 (k64)。v2 は別の形状関数で、粗いメッシュ (縦 4 セル) では gs と 4.5% 違う。`
  - 後: `! k64 モードの PSPAC (K64 の spch.c の形状関数を直訳: 横 CIC・縦「ほぼ最近接」・双 3 次の電位・3 次スプラインの線密度)。gs と全桁一致 (k64)。v2 は別の形状関数で、粗いメッシュ (縦 4 セル) では gs と 4.5% 違う。`
- `flag_WSPAC_tracking_beam.sad` 1 行目:
  - 元: `! v2next P4 項目 4: 追跡の WSPAC のビーム行列。gs (K64) は `CODPLOT; WSPAC; EMIT;` のあとの追跡で EMIT の最後の一周のビーム行列 (beamsize) を使う。`
  - 後: `! 追跡の WSPAC のビーム行列。gs (K64) は `CODPLOT; WSPAC; EMIT;` のあとの追跡で EMIT の最後の一周のビーム行列 (beamsize) を使う。`
- `flag_SPAC_MULT_radius.sad` 先頭付近:
  - 元: `! v2 (担当 C2): MULT`
  - 後: `! MULT`
- `flag_WSPAC_tune_shift.sad` 先頭付近:
  - 元: `! v2 (担当 C2): WSPAC`
  - 後: `! WSPAC`

## 手で書き直した記述
- 「以前の記録の…」「以前の pysad は…」「記録は docs/history/v2next_*.md §…」「(.pysad_k64.txt)」の言及など、複数行にまたがる記述は個別に書き直した (`.sad` 側に残るのは現在の動作の説明だけ)。

## 削除した k64 モードの出力
- この移動で、旧 `v2_next` の `*.pysad_k64.txt` のうち、このフォルダに移る例の分を削除した (k64 モードの出力は `PYSAD_COMPAT=k64` で再現できる。`flag_PSPAC_k64_shape` だけは k64 オプションの例として残した)。
