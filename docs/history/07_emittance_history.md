# 07_emittance.md から取り除いた履歴の記述

07_emittance.md (現在の仕様の解説) から、版・経過に結びついた記述を取り除いた。下は取り除く前の元の行を、節の見出しつきでそのまま並べたもの (行番号は整理前のもの)。

## # 07. Emittance[] / emit — SAD の temit と pysad の emit6d.py

整理前の 3 行目:

````text
> 最終更新: 2026-11-18 (v2.1.0 の後) (現状のスナップショット)。行番号は 2026-10 時点。
````

## ### 解析ヤコビアンの対象と理由

整理前の 123 行目:

````text
- まだ数値微分: `MULT`、ANGLE=0 の BEND、DROTATE≠0 の BEND、空洞 (history/work_emit_report.md §5)。F1 付き BEND の本体は v2.1.0 で `emit6d._tbende_core` (tblfre・tbendal・断片のループ) として実装した (微分の方式は未確認。F1 の無い BEND は従来と同じ演算)。
````

## ## 4b. 追加した要素・機能の範囲

整理前の 129 行目:

````text
## 4b. 追加した要素・機能の範囲
````

## ### SynchroBetaEmittance (`pysad/synchrobeta.py`)

整理前の 152 行目:

````text
- 2026-11-16 の確認: 24 セルの平面リング + CAVI VOLT=2e6、AzimuthalModes->7 の 4 点で emitx は gs と相対 1.4e-7〜1.4e-6 (3.7e-7、4.0e-7、1.4e-6、1.4e-7) で一致、emity は雑音 (gs −2.5e-19 など、pysad 1.3e-20 など)。**結合のあるラティス** (QT・SF・DY のずれ入り 8+4 セルのリング) は gs が SIGSEGV で落ちて比較できず、pysad も `ValueError: cannot convert float NaN to integer` で失敗する (エラーの出し方が粗い)。**v2.1.0 の確認 (C)**: 24 セルの平面 FODO + skew 四極 1 台 (3 GeV、RF 4 MV、6 極なし) の結合のあるリングでは gs は落ちず、k64 の pysad は gs と 7 桁以上一致した (εx 1.3603366e-6 対 1.3603365e-6、εy 4.7516654e-8 対 4.7516654e-8)。Emittance[] と、独立な 6x6 の平衡エミッタンス (自前の 1 周写像 + 放射の減衰 + 量子励起の拡散 + Lyapunov 方程式) とも一致 (ASTRID2 + skew 四極 ks = 0.05 で ε₁ 6e-5、ε₂ 1.3e-4)。SB と Emittance[] の差 (ASTRID2 +0.8%、PF +1.1%・+2.9%) は結合の有無に依らず、シンクロ・ベータ結合と色収差の分。SOL・強い結合で gs が落ちるかは未確認。
````

## ### SOL・COORD・TCAVI・UND・MAP・MULT (ANGLE≠0) と Emittance[]

整理前の 172 行目:

````text
- v2 モードの SOL 領域: 領域内の BEND・MULT・CAVI (RF つき MULT を含む)・F1 つきの端の放射・ずれて傾いたソレノイド・TRPT の加速・AUTOPHI が Emittance[]・EMIT に入る (docs/29・44)。**POL の量 (SpinTune・Polarization など)**: k64 は領域の中が K64 の式 (縦磁場を z 軸の向きに足す: B130) のままで、ずれて傾いたソレノイドではスピンチューンが追跡と 1.5e-4〜1.1e-3、n0 の向きが 2.5e-3〜1.5e-2 rad 違う。**v2 は v2.1.0 で修正** (B⊥ と B∥ を運動量の向きに分けて足す。`spin.tradke_spin`・`solemit._trk`): 追跡との差は SpinTune で 1e-13〜5e-7、n0 で 0〜1e-5 rad (小さな FODO・大きなリング A・PF・ASTRID2・MPW 入りで確認)。回転した縦磁場の歳差の閉じた式とも 5e-7。条件 1〜3 T・0.3〜3 GeV。bz ≥ 3 1/m の QUAD・MULT の D は 1e-4 台 (docs/44)。
````

## ### v2 モードで変わった放射・エミッタンス (RADINT・Jx・MomentumCompaction・EMITDIV ほか)

整理前の 201 行目:

````text
- **EMITDIV**: 大域変数 `EMITDIV`・`EMITDIVB`・`EMITDIVQ` (放射の分割の倍率) が Emittance[] に効く。**両モードで効く** (以前の pysad は無視していた。MAIN 部の `EMITDIV = 20;` も効く)。gs: `EMITDIV=20` で DampingRate −.000286831320204 (pysad k64 と 12 桁一致)。
````

整理前の 217 行目:

````text
| B115 | RADINT が ANGLE=0・K0≠0 の BEND (挿入光源の極) を無視 | K0 を曲げ角として積分に入れる (端の角なし) | ASTRID2 + MPW で U0 6.17 keV → 7.81 keV (Emittance[] 7.80 keV)。v2.4.0: 極を直線の基準軌道の場として積分 + CALC の分散の補正で、εx・Jx は RF なしの Emittance[] と 1e-4・2e-4 (docs/history/v24_RAD.md) |
````

整理前の 225 行目:

````text
| B63・B36・B21 | (バグではない / pysad は有限 / z を閉じる δ) | 変更なし | B63: SOL に L のキーワードは無い (L を書いても軌道・光学に入らない)。pysad は S 座標 (LENG) にだけ入れていたのを直した。B36: L≠0 の TCAVI のリングの Emittance[] は gs が NaN、pysad は有限 (K0=0 で同じ長さのドリフトと 12 桁一致)。B21: RFSW OFF・RADCOD OFF の閉軌道は z が閉じる δ (1e-15) |
````

整理前の 234 行目:

````text
| RADINT | BEND (K1 つきも)・QUAD・SEXT のリングと輸送路 (INS)。先頭が SEXT でもよい。空洞が複数あってもよい (同じ周波数)。GX・GY は一次の色収差 dν/dδ | MULT・SOL は積分に入らない (表示が出る)。ANGLE≠0 の BEND の K0 は入らない (ANGLE=0 の K0 の極は v2 で入る: B115。v2.4.0 で εx も RF なしの Emittance[] と 1e-4)。結合のあるリングの V 面 (I4V) は v2.1.0 で確認 (I4U + I4V が独立な I4 と 1e-7、DU・DV・Je が 1 次の摂動論と 3e-3。K1≠0 の BEND は B107 の修正後)。軌道のずれのある四極は入る (B126)。SEXT・BEND のずれは入らない。空洞は VOLT を単純に足す (PHI を見ない。位相の違う空洞のベクトル和にするかは仕様の判断)。I1〜I5 は独立な積分 (大きなリング A・B・PF・ASTRID2) と 1e-6 以下 |
````

整理前の 235 行目:

````text
| Emittance[]・EMIT | 端の角のある偏向磁石 (E1・E2)、勾配つき偏向磁石。Jx は解析値と 1e-3 (既定)、`EMITDIV=20` で 1e-4 | ANGLE=0 の BEND (補正磁石・挿入光源の極) の AE1・AE2 の端の角は K64 が端の項を約 2 倍にする (B108。v2 は直した)。K1≠0・K0≠0 の ANGLE=0 BEND と MULT の端の角 (AE) は pysad が読まず未確認。既定の分割 (EMITDIV=1) では低エミッタンスのリング (大きなリング A) で εx が 2.1% 大きい (`EMITDIV=5` で 8e-5)。v2 の Emittance[] と RADINT が矩形磁石で 1.4e-3 違う (SP8 (大きなリング B)。未調査) |
````

整理前の 240 行目:

````text
**要判断・直していない**: K1=0 の BEND の追跡のスピン (tbend + tradk) が Emittance[]・tbendi (v2) と 2.2e-5 rad 違う (放射の分割が粗く、スライスの前のエネルギーで歳差を計算する。K64 の精度の範囲)、tthine (長さのある6 極など) の bsir0 (B110。次元が合わない式だが放射なしで呼ばれ入口と出口で打ち消すので影響なし)、RADCOD の 1 周の解析 T と写像の数値微分の 2e-3 の差 (K64 の tradke の微分が近似)、RADINT の SEXT・BEND のずれ、tradke の符号 (B111) と RADTAPER の重みの選択 (「閉じる」を根拠にした判断で、gs の Emittance[] (RADCOD) と 1e-4〜1e-3 違う)。詳細は docs/47 の v2.1.0 の (d) と docs/46。
````

## ### v2.1.0 の追加 (F1・FRINGE・K0 の極・WSPAC・SynchroBeta の警告)

整理前の 242 行目:

````text
### v2.1.0 の追加 (F1・FRINGE・K0 の極・WSPAC・SynchroBeta の警告)
````

整理前の 248 行目:

````text
| **BEND の F1・FRINGE・FB1・FB2** の Emittance[]・EMIT・CALC6D | 両 (K64 の動作の追加) | `F1=0.1 FRINGE=3` など (FRINGE = 1・2・3・−1 (入口側だけ)・−2 (出口側だけ)、F1 + FB1 + FB2、K1≠0 (tbendebody の分割)、K0、RADCOD、DX・DY・ROTATE・DROTATE、POL ON の減衰率)。K64 の `tbende` の長さ補正 `alcorr`・入口の `tblfre`・断片 (K1=0 は n = 0 と ndiv+1 に `tbendef1`、K1≠0 は n = −1, 0 と ndiv+1, ndiv+2 の `tbendal`) の直訳。**周長に F1 の長さ補正が入る** (ω0・運動量圧縮・周回周波数)。k64 は gs と全桁一致: 12 セル FODO (3 GeV、F1 = 0.1 FRINGE = 3) で U0 739673.742271 eV (F1 なし 781633)、Tunes {−.2933271700、−.1395704757、−.0179539450}。大きなリング A (全 BEND に `F1=0.03 FRINGE=3`) の U0 285476.566791 eV (F1 なし 296074.86)・運動量圧縮・減衰率・εx が 1e-7、PF・ASTRID2・MPW 入りも 1e-9。放射・RF なしの Emittance[] のチューン = CALC のチューン (1e-9)。F1 の無い BEND は従来とビット単位で同じ | CALC の解析行列と `tbende` は DISFRIN=1 でも Maxwell のフリンジを入れる (追跡の tbendi は DISFRIN に従う。B134)。逆向きの要素 (`direlc<0`) の F1 の入れ替えは未対応。隣り合う BEND の `bradprev` は tbende では使われない。**TrackParticles の F1 は ANGLE=0 の BEND・MULT で未対応** |
````

