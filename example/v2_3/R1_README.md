# R1: ADTS の gs と pysad の不一致の原因 (台本の RFSW)

原因と切り分けは `docs/history/v23_R1.md`。FFS に入ると RFSW は ON なので、gs の TrackParticles は RF つきのラティスで空洞を働かせ、1 周内の z に応じた δ が (自然クロマティシティ) × δ(J) の偽の ADTS を作る。台本に `NORFSW;` を入れると解消する。

| ファイル | 内容 |
|---|---|
| `Q2_adts_gs.py` | gs の追跡の ADTS (既定で `NORFSW;`。`--rfsw` で RFSW ON のまま) |
| `R1_adts_gs_{fodo,astrid,kekls,kekls_oct0}_norfsw.txt` | gs `NORFSW;`、1024 周 (Q2 の C・B が解析と一致する) |
| `R1_adts_gs_{astrid,kekls_oct0}_rfsw.txt`、`R1_adts_gs_{astrid,kekls_oct0}_rfsw_512.txt` | 対照: RFSW ON、1024 周と 512 周 (×0 の ADTS が周数に依る。旧版 `Q2_adts_gs_*.txt` はこの条件) |

再現: `py -3.12 example/v2_3/Q2_adts_gs.py astrid --turns 1024` (gs は WSL、約 20 s)。大きなリング A は `kekls --cases 1:0,-1:0,0:0`。
