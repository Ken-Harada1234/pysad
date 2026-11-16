# T2_samples (PF の話題): サンプル集
元: KEK SAD チュートリアル「SAD サンプル集」(Tutorial2) の PF のスクリプト。PF 2.5 GeV リング (`../lattice/v9_16_asym.sad`、入射バンプだけ `v9_16_inj_noSF.sad`、周長 187.035 m、チューン 9.5975 / 5.2790)。ASTRID2 の話題 (チューンサーベイ・MPW など) は別の担当の `*_ASTRID2.*` (このフォルダの別のファイル)。

| ファイル (.sad/.gs.txt/.pysad.txt/.md) | サンプル集の節 | 内容 | gs と pysad |
|---|---|---|---|
| `calc_respx_PF` | calc_respx / calc_resp_analy | 軌道応答 (K0 → モニター) と解析式、符号と RFSW の確認 | ○ (応答は解析式と 0.08% で一致) |
| `calc_dispersion_PF` | calc_dispersion | 分散 + 差分と α_c の確認 | ○ |
| `calc_chromaticity_PF` | calc_chromaticity | `PolynomialFit` の 6 次、自然色収差、四極・六極の和の式 | ○ |
| `FFS_FIT_injbump_PF` | ローカルバンプ | 入射バンプ (`ins`、`fit INJ DX`、`free KK* K0`) と解析式 | ○ |
| `fn_DynamicApertureSurvey_DAsv_PF` | sample_DAsv | DA サーベイ (x・y 方向、運動量 9 点、2000 周) | △ (y の 1 行だけ gs 11 / pysad 13) |
| `fn_DynamicApertureSurvey_da_xy_PF` | da_xy | X-Y 面の DA (24 方向、1000 周) と図 | ○ |
| `FFS_mom_dep_ts_PF` | mom_dep_ts | 運動量依存チューンシフト (六極あり・なし) と図 | ○ |
| `FFS_amp_dep_ts_PF` | amp_dep_ts | 振幅依存チューンシフト (12 点平均) と図 | ○ |
| `calc_RadInt_PF` | calc_RadInt | 放射積分、εx・U0・σE・Jx の式と数値の照合 | ○ (εx は SAD と 0.24%) |
| `phase_PF` | phase | 位相空間 (60 周) と図 | ○ |
| `FFS_draw_topdrawer_PF` | topdrawer 形式の図 | topdrawer ファイルの手書き。`draw` は gs で動かず (出力なし/止まる)、pysad は未対応 | ○ (手書き部分) / `draw` は × |

- 走らせ方: `py -3.12 tools/runpf.py T2_samples/NAME.sad` (全部を回すときは `*_PF.sad` だけを指定すること。このフォルダは ASTRID2 の担当と共有)。`.sad` は `tools/build_T2.py`、図は `tools/make_figs.py`。
- 長い計算 (DA) は 5000 周 → 1000〜2000 周、方向 161 → 24 などに短縮した (各 `.md` に書いてある)。
- 未対応・書けなかったもの: `calc_respy` (垂直の応答)、`drawmagvw3d` (CAD 出力、CAD に読み込む `.vss` は作っていない)、`modules_*` の GUI 版。
