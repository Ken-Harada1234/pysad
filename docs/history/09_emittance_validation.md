# 09_emittance: gs (K64) との比較・検証の詳細 (本文から移したもの)

本文 `docs/09_emittance.md` から、gs との比較の細部を元の文脈つきでそのまま移した。

## 元の節: 2-4. 報告 (temitprint / `FFS.print_emit`)

実ラティスの RFSW × RADCOD 4 通り (`bench_level26a`) の 148 行が SAD と一致 (0 同然の雑音の桁を除く)。

## 元の節: POL のスピン (docs/16_spin_pol.md)

- 一致 (gs): 平面リングで偏極時間は全桁、スピンチューン 13 桁、偏極ベクトル 12 桁。縦ずれ・傾き・6 極のあるリングで 13/12/10 桁以上。POL 付きの Emittances・DampingRate は 11〜12 桁。エネルギー 3.0〜3.08 GeV の格子 (aγ=6.8〜6.99) で確認した範囲。

## 元の節: SynchroBetaEmittance (`pysad/synchrobeta.py`)

- gs との一致 (12 セルのリング + CAVI VOLT=2e6 HARM=64): nus=0.01・AzimuthalModes=9 で emitx 相対 6e-6 (conv は 2e-4)。nus=0.005〜0.02・AzimuthalModes=5 の 4 点で emitx 4e-7〜1e-5。emity は 1e-18 の雑音で符号も違う。試験は [`test_synchrobeta.py`](../pysad/tests/test_synchrobeta.py) (許容 3e-5)。

## 元の節: SynchroBetaEmittance (`pysad/synchrobeta.py`)

- 2026-11-16 の確認: 24 セルの平面リング + CAVI VOLT=2e6、AzimuthalModes->7 の 4 点で emitx は gs と相対 1.4e-7〜1.4e-6 (3.7e-7、4.0e-7、1.4e-6、1.4e-7) で一致、emity は雑音 (gs −2.5e-19 など、pysad 1.3e-20 など)。

## 元の節: SOL・COORD・TCAVI・UND・MAP・MULT (ANGLE≠0) と Emittance[]

| 要素 | Emittance[]・EMIT の扱い | gs との一致 | 文書 |
|---|---|---|---|
| SOL (BOUND≠0 の領域) | `solemit.py` (領域内を K64 の `tsole1` の順に、写像は追跡カーネルの数値微分、放射は bz 付きの `_State.tradke(al, 0, bzh)`) | ε_x 2e-12〜4e-9、ε_y 4e-10〜2e-8、ε_z 1e-7〜3e-6、チューン 1e-13〜1e-11 (BZ=2 と BZ=5 の 12 セルのリング、RFSW) | docs/17 |
| COORD | `emit6d._tcoord` (`piece`、z 行のシンプレクティック組み直しなし) | ε・チューン・COD・運動量収縮因子 5e-10 以下 | docs/19 |
| TCAVI | `emit6d._ttcave` (tchge → 半ドリフト → キック → 半ドリフト [tradke] → tchge)。RFSW OFF はドリフト | **L=0 のみ gs と一致** (ε 1e-11)。L≠0 は gs が壊れる (B36: EnergyLossU0 が NaN など) ので pysad は NaN を再現せず自然な値を出す | docs/19 |
| UND | k64 モード: 恒等 (K64 の `tturne` は `case default`)。K64 は磁場を設定できず磁場なしのドリフト (B37)。v2 モード: 磁場がなくても長さ L のドリフト、`FBX`・`FBY` つきは追跡のカーネルの数値微分と 1 周期 64 分割以上の `tradke` | k64: 1e-12 台。v2: gs は磁場を設定できないので比較なし (解析式で確認) | docs/19 |
| MAP | `emit6d._tmap` (`ExternalMap["EMIT", l, cod]` の `{cod1, trans[, dtrans, dbeam]}`)。定義がなければ恒等・長さ 0 | ε・チューン・減衰率・軌道 1e-11 | docs/19 |
| MULT (ANGLE≠0) | `emit6d._tmulte` → `multae.emit` | ε 1e-11、チューン 1e-13 (12 セルのリング) | docs/19 |

## 元の節: v2 モードで変わった放射・エミッタンス (RADINT・Jx・MomentumCompaction・EMITDIV ほか)

gs (2026-10-19 ビルド) の Emittance[] の DampingRate は Jx = 1.0974 (pysad の k64 と 12 桁一致)、同じ gs の RADINT は 1.00575100: K64 の中で EMIT と RADINT が食い違っており、正しいのは RADINT。

## 元の節: v2 モードで変わった放射・エミッタンス (RADINT・Jx・MomentumCompaction・EMITDIV ほか)

gs: `EMITDIV=20` で DampingRate −.000286831320204 (pysad k64 と 12 桁一致)。

## 元の節: F1・FRINGE・K0 の極・WSPAC・SynchroBeta の警告)

k64 は gs と全桁一致: 12 セル FODO (3 GeV、F1 = 0.1 FRINGE = 3) で U0 739673.742271 eV (F1 なし 781633)、Tunes {−.2933271700、−.1395704757、−.0179539450}。大きなリング A (全 BEND に `F1=0.03 FRINGE=3`) の U0 285476.566791 eV (F1 なし 296074.86)・運動量圧縮・減衰率・εx が 1e-7、PF・ASTRID2・MPW 入りも 1e-9。

## 元の節: F1・FRINGE・K0 の極・WSPAC・SynchroBeta の警告)

k64 は gs と全桁一致 (12 セル FODO の 14 通りの F1・FB1・FB2・AE1・AE2・FRINGE・DISFRIN・K0・K1・K2・DX・DY・ROTATE の組み合わせ。大きなリング A・PF・ASTRID2・MPW 入りにステアリングまたは MULT を 2 台挿入)。

## 元の節: F1・FRINGE・K0 の極・WSPAC・SynchroBeta の警告)

**k64 は gs と全桁一致** (K64 の最後のスライスの 8 倍のレンズ: B127 を再現。小リング、大きなリング A・PF・ASTRID2 で PBUNCH 1e9〜1e10)。

## 元の節: 6. 照合の要約

## 6. 照合の要約

| Level | 内容 | 一致 |
|---|---|---|
| 2 | 20 セル FODO, combined BEND, RF なし | εx 表示桁まで一致 (≤4e-10), U0 1e-15 |
| 3C/3D | 結合 (ロール) / 強さ誤差 | εx, εy, σE, U0 ≤7.8e-7 / ≤1.6e-9 |
| 4 | RFSW ON (RADCOD OFF/ON) | εx, εz, σE, σz, Q, 減衰率, U0, trf0, α, bh ≤1e-7 / ≤6e-8 |
| 8 | 四極 F1/FRINGE, 空洞 L=0.5, MULT K2 | εx, σE, σz, α, U0 表示桁 (≤2.3e-6) |
| 11 | MULT K0/K1/K2 + DX/DY/ROTATE | trf0 2.6e-10 |
| 28c/28d/28f | 不安定な光学 (V: QD×1.19、H: QF×0.9) の EMIT 報告・Emittance[] | V は行単位で一致 (雑音の桁のみ)、H は SAD の T の悪条件の範囲 (≤3e-5) |
| 30b | RFSW ON + RADCOD OFF + FSHIFT = 1e-6, 3e-7, −2e-6, 0 | ≤1.2e-11 (tcod の直訳) |
| 25a/25b/25c/25d/25h | ROTATE: 垂直ステアリング (ANGLE=0, ROTATE=−90 DEG, K0)、skew 四極、ROTATE 付き BEND (5/90 DEG, K1 付き)、6 極・8 極・十極・MULT の ROTATE、L=0 の SQ | εx/εy/U0 相対 ≤1e-8 (εy 2e-12〜1e-8)、TUN ≤4e-9 (25h は δ の収束差で 3e-8) |
| 25e | ROTATE/DX/DY/DROTATE/K0 の組み合わせ (癖のない組) | εx/εy 相対 1.6e-6 (δ の収束差)、TUN 4.5e-8 |
| 25f | th≠0 かつずれ≠0 の要素 (tchge の癖) | CALC・追跡は一致、Emittance[] の δ・チューンは SAD の tcod が閉じないため不一致 (2-1 節) |
| 12 | 試験用のラティス (3563 要素) | 表示桁 (相対 ≤2e-5) |
| 16b/16e | 同、全桁 | εx, εz, σE, σz, νx, νy, νs, α_c, bh, trf0, U0 **≤1.3e-11** |
| item6 (`bench_item6_emit.sad`) | 試験用のラティス, RFSW+**RADCOD**, `Emittance[Matrix->True]` の全桁 | εx 1.1e-12, εz 9e-12, σE 2e-14, σz 9e-12, α_c 2e-11, ν ≤9e-12, 減衰率 ≤8e-13, ビーム行列 ≤4e-11, U0・trf0 (dz) は全桁 |
| item6 (`bench_item6_tiny.sad`) | 実ラティス 1 セルのリング (ν_s=2e-4): RFSW OFF / ON / RADCOD | 横方向・U0・trf0 ≤1.3e-11、縦方向 (RFSW ON) は 1e-9 (νs が小さく 1e4 倍に拡大)。テスト [`test_emit_jac.py`](../pysad/tests/test_emit_jac.py) |
| 26a/26c/26d | EMIT 報告の全行 (実ラティスの RFSW × RADCOD 4 通り、FSHIFT=1e-6)、RFSW OFF + RADCOD ON の Emittance[] と CALC6D | 報告 148 行 + 37 行すべて一致 (雑音の桁を除く)、26d ≤2.2e-11 |
| 26b | FRINGE=1,2,3 付き四極、ずれ・ロール付き四極・6 極 (小リング): CALC と Emittance[] | CALC ≤5e-13、Emittance[] ≤1.3e-9 (数値微分では 1.1e-5) |


## 追記: EMITDIV の既定の記述 (docs/09_emittance.md)

旧い記述: 「既定の分割 (EMITDIV=1) では v2 の Emittance[] の Jx が解析値より 1.5e-3 大きい」「Jx は解析値と 1e-3 (既定)」「既定の分割 (EMITDIV=1) では低エミッタンスのリング (大きなリング A) で εx が 2.1% 大きい」。
いまの v2 の既定は EMITDIV=4 (pysad/compat.py)。上の数値は EMITDIV=1 で測ったもので、本文では条件を EMITDIV=1 と明記した。
