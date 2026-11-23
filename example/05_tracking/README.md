# 05_tracking (追跡)

使ったラティスは `00_lattice/` の `FODO_ring.sad` (24 セル、3 GeV、周長 105.6 m)、`FODO_ring_V4MV.sad` (RF 4 MV)、新規の `apert_ring.sad`、`linac_line.sad`、`spac_ring.sad`。
各例は 4 ファイル (`NAME.sad`、`NAME.gs.txt`、`NAME.pysad.txt`、`NAME.md`)。一致の記号: ○ 一致、△ 差あり、× gs 失敗または pysad 未対応。

| 名前 | 1 行説明 | gs と pysad |
|---|---|---|
| [fn_TrackParticles](fn_TrackParticles.md) | 1 周・N 周・周ごとの座標、転送行列・チューン・不変量との照合 | ○ |
| [fn_TrackParticles_loss](fn_TrackParticles_loss.md) | LOSSAMPL と flag 列 (損失) | ○ |
| [fn_SurvivedParticles](fn_SurvivedParticles.md) | 生存粒子だけを取り出す | ○ |
| [fn_DynamicApertureSurvey_flag_DAPERT](fn_DynamicApertureSurvey_flag_DAPERT.md) | 動的口径の地図 (50 周) | △ (Score 一致、地図の数文字が 1 周違う) |
| [elem_APERT_LOSSAMPL](elem_APERT_LOSSAMPL.md) | APERT 要素と LOSSAMPL | ○ (APERT 要素・LOSSAMPL とも一致。SOL 領域内の APERT は v2 が対応 (k64 は未対応)、gs と未比較) |
| [flag_TRPT](flag_TRPT.md) | 輸送ライン (CAVI で加速)、断熱減衰 | △ (粒子の座標は一致。gs は CAVI の PHI を変えても Twiss が更新されない) |
| [fn_AccelerateParticles](fn_AccelerateParticles.md) | 周ごとに運動量を上げる追跡 | ○ |
| [flag_RAD_FLUC](flag_RAD_FLUC.md) | 放射減衰と量子ゆらぎ (乱数固定) | ○ |
| [flag_JITTER](flag_JITTER.md) | JITTER/QUIET (TrackParticles には効かない) | ○ |
| [multi-turn-tracking](multi-turn-tracking.md) | EMIT のあとの多周回追跡、振幅依存チューン | ○ |
| [flag_SPAC](flag_SPAC.md) | 空間電荷 (1 m のドリフトの広がりとリングのチューンシフトを解析式と比べる。RADIUS<0 は絞りだけ) | △ (RADIUS<0 の絞りは一致。RADIUS>0 の空間電荷は v2 が実装していて解析式に合うが、gs (K64) は力が 1/γ 倍に弱いので意図して違う。k64 モードは空間電荷も絞りも入らず gs と合わない: 2026-11-17 確認) |
| [flag_LOSSMAP](flag_LOSSMAP.md) | 損失の要素・周を戻りの 8・9 列目に | ○ |
| [flag_SPAC_MULT_radius](flag_SPAC_MULT_radius.sad) | MULT の RADIUS の空間電荷 (SPAC)。四極を MULT で書くと全周に入る。四極の中の KV ビームの包絡 (包絡方程式と 0.5%)、リングのチューンシフト (SPAC・PSPAC・WSPAC と Laslett の式。SPAC・PSPAC は粒子の雑音で 5%) | v2 専用 (K64 は MULT の中の spkick が壊れている。k64 モードは空間電荷なしで走る) |
| [flag_WSPAC_tune_shift](flag_WSPAC_tune_shift.sad) | WSPAC (ビーム行列のガウス分布の力) と SELFCOD、`WSPACBEAM`。追跡のチューンシフトが Laslett の式と 1% で一致。平らなリングはエラーで知らせる (弱いスキュー四極で結合を入れた 0.5 GeV の 8 セルリング) | △ (gs (K64) とも 1%。k64 モードは WSPAC を追跡に入れない) |
| [flag_WSPAC_tracking_beam](flag_WSPAC_tracking_beam.sad) | 追跡の WSPAC のビーム行列。既定 (ゼロ電流) と v2 のオプション `WSPACSC = 1` (空間電荷つき・MINCOUP の下限) の試験粒子のチューンシフト | △ (gs (−2.938e-5、−1.210e-3) に既定が x 1.5%・y 0.7%。`WSPACSC = 1` は y が 3 倍違う: K64 の追跡には下限が入らない。k64 は追跡に WSPAC を入れない) |
| [flag_WSPAC_emittance](flag_WSPAC_emittance.sad) | WSPAC の Emittance[]・EMIT・CALC6D・CALC | △ (Emittance[] のシフトは gs より 4% (x)・3% (y) 弱い (原因未特定、一次のラスレットの式とは 0.15%)。gs の CALC・CALC6D は NaN) |
| [flag_PSPAC_resolution](flag_PSPAC_resolution.sad) | PSPAC のメッシュの分解能 (縦: λ(z) の入れ方、横: 電位の補間) | 台本の先頭に別セッションで取った gs の値と解析式を記載 (gs の 2 行目以降は K64 の static 変数のため無効) |
| [flag_PSPAC_astrid2](flag_PSPAC_astrid2.sad) | ASTRID2 (0.58 GeV、45.7 m) に SPCH 3 か所、PBUNCH = 1e14、10 周 | ○ (k64 モードは 11 桁以上一致。v2 は x が 8% 違う) |
| [flag_PSPAC_k64_shape](flag_PSPAC_k64_shape.sad) | k64 モードの PSPAC (K64 の形状関数: 横 CIC・縦「ほぼ最近接」・双 3 次の電位・3 次スプラインの線密度)。5 MeV・1e11・縦 4 セルの粗いメッシュ。`.pysad_k64.txt` は k64 オプションの出力 | ○ (k64 モードは全桁一致。v2 は別の形状関数で rms x が 4.5% 違う) |
| [flag_FLUC_nparallel](flag_FLUC_nparallel.sad) | FLUC つきの追跡の乱数列 (NPARA=1 と 4)。gs の既定は NPARA=4 (粒子を子プロセスに分け、子の種は親の種 + 2k) | ○ (k64 モードは rms dp・最終の種まで全桁一致。v2 は 6e-6 (1 周)〜1e-3 (100 周) 違う) |
| [flag_LOSSMAP_RADTAPER](flag_LOSSMAP_RADTAPER.sad) | LOSSMAP の 2 行 (落ちた要素と周)、RADTAPER。RADTAPER で閉軌道が −79 μm から −0.2 μm に戻り U0 が 0.1% 増える | (gs と一致) |
| [fn_DynamicApertureSurvey_dapwidth_window](fn_DynamicApertureSurvey_dapwidth_window.sad) | `DynamicApertureSurvey` の DAPWIDTH の窓 (7・3・12・51)。RAD + FLUC の最終 seed (DAPWIDTH = 2) | ○ (NORAD の地図・Score が一致 (境界の 1 点 1 周差)。seed は k64 で 4/5 が一致。v2 は放射の式が違うので合わない) |
| [fn_DynamicApertureSurvey_dapwidth_ring_a](fn_DynamicApertureSurvey_dapwidth_ring_a.out.txt) | 大きなリング A (kek-ls_v03_12.lat、3 GeV、603 m) の DAPWIDTH の窓 (7 と 51) の地図 (台本なし、出力のみ) | ○ (地図・Score 63 が一致。境界の 1 点 1〜2 周差) |
| [fn_NAFF_tunes](fn_NAFF_tunes.py) | 追跡の座標列からチューンを取り出す NAFF。4 つのリング (12 セル FODO、大きなリング A、PF、ASTRID2) の追跡のチューンと CALC の比較 (1024 周、小さな振幅。結果は `.out.txt`) | 追跡のチューンが CALC と 1e-9 以内。wtunem (FFT) は 1e-4。gs は使わない |

上の表のうち、拡張子つきのリンクの行 (`.sad`・`.py`) は説明の `.md` を持たない。内容と結果の要点は表と、台本・スクリプトの先頭のコメントに書いた。`.py` の例は、リポジトリ直下から `py -3.12 example/05_tracking/NAME.py` で走らせる (出力は `NAME.out.txt`)。gs と同じ動作で走らせるときは `PYSAD_COMPAT=k64`。

図: `fn_TrackParticles_phase.png`、`multi-turn-tracking_tune_amp.png`、`flag_RAD_FLUC_damping.png` (pysad の出力数値から matplotlib で描いた。gs の画面表示は撮れない)。

## このフォルダで分かったこと
- SAD の `E1`/`E2` は ANGLE に対する比として読まれる (`FODO_ring.sad` の E1=0.1309 の端の角度は 0.0343 rad)。この前提で行列の積がチューン・β・α と 10 桁一致する (fn_TrackParticles.md)。
- `FFS` の中で `MOMENTUM = ...` を設定すると、gs と同じく追跡・Emittance[] に反映され、`USE` で MAIN 側の値に戻らない (2026-11-17 確認)。
- `Do` の中のフラグの切り替えは `FFS["NAME"]` と書く (`NAME;` は Do の中では実行されない)。
- 未解決: LOSSMAP、SPAC。gs の SPAC は PBUNCH=1e13 でも効果を確認できなかった。

## 未確認・限界
- 大きなリング A (603 m) で確かめたのは DAPWIDTH だけ。PSPAC・WSPAC の追跡は ASTRID2 (0.58 GeV) までの確認で、PF リング・8 GeV のリングは使っていない。
- DAPWIDTH: `NPARA > 1`・`Damping->True` は未対応。損失の境界の粒子は gs と 1〜数周違うことがある。gs の追跡の例は `NPARA = 1` で取った (gs の既定は粒子を子プロセスに分け、PSPAC は各プロセスが自分の粒子だけで電位を作る)。
- PSPAC (k64 モード): gs の `NPARA > 1` の再現と、1 セッションの 2 回目の `TrackParticles` (gs は PSPAC が働かない) は再現しない。`PSPACNZ < 3` はエラー。
- WSPAC の追跡: 強い空間電荷 (0.2 GeV で PBUNCH ≳ 1e9) は反復が収束しない。
