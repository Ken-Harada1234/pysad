# 04_emittance_radiation 一覧
ラティスは `00_lattice/FODO_ring_V4MV.sad` (FODO_ring の RF 電圧だけ 4 MV に変更)。`RINGQ` を使う (RADINT の NaN を避けるため)。式の検算は `formula_check.py` → `formula_check.txt`。

| ファイル | 説明 | gs と pysad |
|---|---|---|
| FFS_EMITTANCE | EMIT の報告と放射積分の式による検算 (U0, τ, εx, σδ, νs, σz) | ○ |
| FFS_RADINT | RADINT の I1〜I5 と自前の積分 | ○ |
| fn_Emittance | Emittance[] の規則 | ○ |
| flag_RFSW_RAD_RADCOD_FLUC | 4 つのフラグの効き方 | △ (Emittance[] は一致、FLUC の追跡は乱数列が違う) |
| flag_RADTAPER | RADTAPER の前後 | ○ (v2 で対応。k64 モードで gs と一致。v2 の既定では εx が 2.4% 違う: Jx の修正 B75) |
| flag_CALC6D | CALC6D の縦の光学 | ○ |
| fn_SynchroBetaEmittance | シンクロ・ベータ結合の emitx (平面リング) | ○ (結合のあるラティスは gs が落ちる) |

## このフォルダで分かったこと
- RADINT の $J_x$ (0.9415) と EMIT の $J_x$ (0.9651) は 2.4% 違う。$I_4$ の端の項を 2 倍にすると EMIT に近づく (FFS_EMITTANCE.md)。**v2.0.0 で原因が分かり直した**: K64 の `Emittance[]`・EMIT が端の角の蹴りの行路を入れず $I_4$ の端の項が 2 倍になっていた (B75。RADINT が正しい)。v2 の既定では EMIT の $J_x$ は 0.9429 で RADINT に近い (gs と違う値になる)。`PYSAD_COMPAT=k64` では gs と同じ 0.9651。この節の他の例の pysad の出力 (`.pysad.txt`) は v2 のもの。
- 例 (`FFS_EMITTANCE`、`FFS_RADINT`、`fn_Emittance`、`flag_RFSW_RAD_RADCOD_FLUC`、`fn_SynchroBetaEmittance`) の判定は k64 モードでの gs との一致。v2 の既定では εx・減衰率・α などが B75・B77 の修正の分だけ gs と違う。
- `FODO_ring` の RF 1.5 MV は $U_0$ より小さく、バケットがない。
