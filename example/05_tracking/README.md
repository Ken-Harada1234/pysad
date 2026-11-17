# 05_tracking (追跡)

使ったラティスは `00_lattice/` の `FODO_ring.sad` (24 セル、3 GeV、周長 105.6 m)、`FODO_ring_V4MV.sad` (RF 4 MV)、新規の `apert_ring.sad`、`linac_line.sad`、`spac_ring.sad`。
各例は 4 ファイル (`NAME.sad`、`NAME.gs.txt`、`NAME.pysad.txt`、`NAME.md`)。一致の記号: ○ 一致、△ 差あり、× gs 失敗または pysad 未対応。

| 名前 | 1 行説明 | gs と pysad |
|---|---|---|
| [fn_TrackParticles](fn_TrackParticles.md) | 1 周・N 周・周ごとの座標、転送行列・チューン・不変量との照合 | ○ |
| [fn_TrackParticles_loss](fn_TrackParticles_loss.md) | LOSSAMPL と flag 列 (損失) | ○ |
| [fn_SurvivedParticles](fn_SurvivedParticles.md) | 生存粒子だけを取り出す | ○ |
| [fn_DynamicApertureSurvey_flag_DAPERT](fn_DynamicApertureSurvey_flag_DAPERT.md) | 動的口径の地図 (50 周) | △ (Score 一致、地図の数文字が 1 周違う) |
| [elem_APERT_LOSSAMPL](elem_APERT_LOSSAMPL.md) | APERT 要素と LOSSAMPL | ○ (APERT 要素・LOSSAMPL とも一致。SOL 領域内の APERT は v2 で対応、gs と未比較) |
| [flag_TRPT](flag_TRPT.md) | 輸送ライン (CAVI で加速)、断熱減衰 | △ (粒子の座標は一致。gs は CAVI の PHI を変えても Twiss が更新されない) |
| [fn_AccelerateParticles](fn_AccelerateParticles.md) | 周ごとに運動量を上げる追跡 | ○ |
| [flag_RAD_FLUC](flag_RAD_FLUC.md) | 放射減衰と量子ゆらぎ (乱数固定) | ○ |
| [flag_JITTER](flag_JITTER.md) | JITTER/QUIET (TrackParticles には効かない) | ○ |
| [multi-turn-tracking](multi-turn-tracking.md) | EMIT のあとの多周回追跡、振幅依存チューン | ○ |
| [flag_SPAC](flag_SPAC.md) | 空間電荷 (1 m のドリフトの広がりとリングのチューンシフトを解析式と比べる。RADIUS<0 は絞りだけ) | △ (RADIUS<0 の絞りは一致。RADIUS>0 の空間電荷は v2 で実装し解析式に合うが、gs (K64) は力が 1/γ 倍に弱いので意図して違う。k64 モードは空間電荷も絞りも入らず gs と合わない: 2026-11-17 確認) |
| [flag_LOSSMAP](flag_LOSSMAP.md) | 損失の要素・周を戻りの 8・9 列目に | ○ (v2 で対応) |

図: `fn_TrackParticles_phase.png`、`multi-turn-tracking_tune_amp.png`、`flag_RAD_FLUC_damping.png` (pysad の出力数値から matplotlib で描いた。gs の画面表示は撮れない)。

## このフォルダで分かったこと
- SAD の `E1`/`E2` は ANGLE に対する比として読まれる (`FODO_ring.sad` の E1=0.1309 の端の角度は 0.0343 rad)。この前提で行列の積がチューン・β・α と 10 桁一致する (fn_TrackParticles.md)。
- `FFS` の中で `MOMENTUM = ...` を設定すると、gs と同じく追跡・Emittance[] に反映され、`USE` で MAIN 側の値に戻らない (2026-11-17 に対応。以前の pysad は USE で戻した)。
- `Do` の中のフラグの切り替えは `FFS["NAME"]` と書く (`NAME;` は Do の中では実行されない)。
- 未解決: LOSSMAP、SPAC。gs の SPAC は PBUNCH=1e13 でも効果を確認できなかった。
