# example/v2_next: 担当 P2 の例 (v2next、2026-11-19)

記録は `docs/history/v2next_P2.md`。実行は `py -3.12 run_sad.py P2_xxx.sad` (既定は v2。`PYSAD_COMPAT=k64` で K64 の動作)、`.py` は `py -3.12 example/v2_next/P2_xxx.py`。
gs の出力は WSL の gs (10-19 ビルド) で取った (`.gs.txt`)。`.pysad_k64.txt` は k64 モードの出力。

| 例 | 項目 | pysad と gs の関係 |
|---|---|---|
| `P2_sb_dispersion_start.sad` | 1: 始点を変えた 24 セルの FODO の SynchroBetaEmittance | gs と k64 は SB が Emittance[] の 6.54 倍 (始点の分散 η ≠ 0。1e-7 で一致)。v2 は 1.0003〜1.0004 倍 |
| `P2_sb_rings.py` (`.out.txt`、`P2_sb_ring_a.gs.txt`) | 1: 大きなリング A (kek-ls_v03_12.lat)・PF・ASTRID2 の SB と Emittance[] の比 | A の k64 は 39.88 倍 (gs と 2e-7 で一致: 6.178294e-9)、v2 は 1.0084 倍 |
| `P2_emitdiv_rect.py` (`.out.txt`) | 2: Emittance[] の EMITDIV の収束 (矩形・扇形の FODO) | SP8 (大きなリング B) の数値はファイルの冒頭の説明 (名前・パスは書かない) |
| `P2_tradke_exact.py` (`.out.txt`) | 3: 放射反作用の ddpx の符号 | 放射反作用つきのローレンツ力の積分と、v2 の tradke・追跡 (tradkn) が一致、k64 の tradke は既定の刻みで 9% ずれる。gs は走らない (Python の内部関数を使う) |
| `P2_radcod_derivative.py` (`.out.txt`) | 4: RADCOD の解析 T + D と数値微分 | v2 は BEND 1 台で 1.7e-10 (k64 は 5e-7) |
| `P2_bend_k1_face.py` (`.out.txt`) | 5: K1≠0 の BEND の端の模型 (オプション `face`) | CALC・追跡が厳密な写像と 1e-10〜1e-12 (それぞれの模型) |
