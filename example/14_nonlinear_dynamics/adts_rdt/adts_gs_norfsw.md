# SAD (gs) で ADTS を追跡で求めるときの RF の扱い (`NORFSW;`)

**この文書で分かること**: gs (SAD 本家) の TrackParticles で ADTS を求めるとき、台本に `NORFSW;` (RF を止める) を入れる必要がある理由と、その確認用の出力ファイル。pysad の追跡 (RF 止め) と同じ条件で比べるための注意です。
**前提**: ADTS が何かは `rdt_adts_basic.md`。

FFS に入ると RFSW (RF 空洞を働かせるスイッチ) は ON になる。そのため gs の TrackParticles は RF つきのラティスで空洞を働かせ、1 周内の z に応じた δ が (自然クロマティシティ) × δ(J) の偽の ADTS を作る。6 極 ×1 では ξ ≈ 0 なので見えないが、6 極 ×0 や ×−1 の ADTS が食い違う。台本に `NORFSW;` を入れると解消する (pysad の追跡は RF を止めた条件)。RF が付いたまま ADTS を求めると、ADTS はシンクロトロン振動の位相と周回数に依存し、1 周のずれで 2e-4 の差も出る。

| ファイル | 内容 |
|---|---|
| `adts_gs.py` | gs の追跡の ADTS (既定で `NORFSW;`。`--rfsw` で RFSW ON のまま) |
| `adts_gs_{fodo,astrid,kekls,kekls_oct0}_norfsw.txt` | gs `NORFSW;`、1024 周 (6 極 2 次の C と 8 極 1 次の B が解析と一致する) |
| `adts_gs_{astrid,kekls_oct0}_rfsw.txt`、`adts_gs_{astrid,kekls_oct0}_rfsw_512.txt` | 対照: RFSW ON、1024 周と 512 周 (×0 の ADTS が周数に依る。`adts_gs_*.txt` (`_norfsw` なし) は RF ON のままの出力) |

再現: `py -3.12 example/14_nonlinear_dynamics/adts_rdt/adts_gs.py astrid --turns 1024` (gs は WSL、約 20 s)。大きなリング A は `kekls --cases 1:0,-1:0,0:0`。
