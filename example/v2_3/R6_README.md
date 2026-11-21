# R6: 解析式が入れていなかった物理の項の追加と、誤りの修正 (a)〜(e) (v2.3.0、2026-11-21。V2 で `adts()` の既定を物理的な項だけ (`fourth='r3'`・`optics='calc'`) に戻した)

式・導出・検証の表は `docs/49` §11、入門は `docs/50` §7.6、記録は `docs/history/v23_R6.md`、試験は `pysad/tests/test_v23_R6_bend_orbit.py` (14 件)、コードは `pysad/bendnl.py` (閉軌道・δ のまわりの分解、K0≠0、閉形式)・`pysad/resdrive.py` (`adts` の `fourth`・`optics`)・`pysad/rdt2nd.py` (非対角の正規形の補正、4 次の ADTS)・`pysad/rdttrack.py` (`tracked_rows`、`element_h4` の軌道・δ)・`pysad/nform.py` (δ の ADTS と δ 微分)。
ラティス: FODO (`example/00_lattice/FODO_ring.sad`、SF・SD の K2L = 4, −6、3 GeV、105.6 m)、ASTRID2 (`example/astrid/lattice/ASTRID2_v01_01.sad`、0.58 GeV、45.7 m)、ASTRID2 MPW (`ASTRID2_v01_02_id.sad` の `RINGMPW`)、KEK-LS (`pysad/tests/kek-ls_v03_12.lat`、大きなリング A。ESRF 型の HMBA、3 GeV、603 m)。
条件 (すべて): pysad k64 の追跡 (RF 空洞止め、放射なし、結合なし、閉軌道のまわり)、6 極 ×1 (指定以外)、直接の追跡は J の格子 6×6 (±位相の対 = 70 粒子)・1024 周・NAFF Hann p=1・order 2 のフィット (Jmax は FODO 1e-7、ASTRID2・MPW 2e-8、KEK-LS 2e-10)、3 次のテイラー展開の振幅 amp3 = 1e-4 √m (KEK-LS 1e-5)、偏向磁石の分割 16。実行は `PYTHONUTF8=1`。**スクリプトは出力を `example/v2_3/` へ上書きする**。

| ファイル | 内容 |
|---|---|
| `R6_resolve.py` | `mpw`・`delta <fodo\|astrid\|kekls>`・`bend`・`nf`・`adts`: (a)〜(e) の表を作る |
| `R6_mpw.out.txt` | (a) MPW: h・f4・ξ・ADTS (λ = 0, 1, 2、追跡の Jmax・order) |
| `R6_delta_<ring>.out.txt`・`.track.json` | (b) δ = ±0.002 の ADTS と ∂(ADTS)/∂δ: 解析 (精密版 `fourth='r6'`・`optics='tracked'`、既定 `fourth='r3'`・CALC の光学)・数値の正規形・直接の追跡 (追跡は `.track.json` に保存) |
| `R6_bend.out.txt` | (c) 偏向磁石の閉形式 対 取り出し、スライス数への依存 |
| `R6_nf.out.txt` | (d) 非対角の正規形の補正の独立な検証 (数値のリー変換) |
| `R6_adts.out.txt` | (e) `adts()` の 4 次の項: 既定 (`fourth='r3'`) 対 精密版 (`fourth='r6'`) 対 追跡 |

## 前後の表 (解析と追跡の相対差)

以下の「R6」の列は、精密版 (`fourth='r6'`、δ≠0 は `optics='tracked'` も) の値 (数値積分器の 0.1〜0.3% の端数を含む)。既定 (`fourth='r3'`・`optics='calc'`) は物理的な項だけで、追跡との差は R3 の値 (e の表の「R5 まで」の列) と δ≠0 の表 (docs/49 §11.4)。追跡との一致は入れた効果の検証の結果で、追跡に合わせて項を足したのではない。

| 項目 | R5 まで | R6 | 原因・直したもの |
|---|---|---|---|
| (a) MPW の h10020・h10200 | 4.6e-5・2.0e-5 | 2.5e-9・5.5e-10 | K0≠0 の偏向磁石 (内部の軌道が ±0.032 rad 傾く) を入れていなかった。閉軌道のまわりで取り出す |
| (a) MPW の f4 (6 極 ×0、19 項の最大) | 8.7e-3 | 2.6e-6 | 同上 (4 次の項は要素の写像から、軌道のまわり) |
| (a) MPW の ξy (CALC − 解析) | −0.173 | +3.5e-4 | K0 の磁石の δ の項を要素の追跡から (スロープ座標) |
| (a) MPW λ=2 の ADTS (xx, xy, yy) | −2.1e-3, −2.0e-3, −5.8e-3 (R3) | −1.5e-3, −1.8e-3, −5.3e-3 (Jmax 2e-8)、−1e-4, −3e-5, −4.9e-4 (Jmax 5e-9)、+2.5e-5, +2.6e-4, +2.9e-4 (Jmax 1e-8、order 3) | 残りは物理でなく追跡側のフィットの偏り (4νy = 8.997 が整数 9 の近く) |
| (b) ASTRID2 の ADTS (δ = ±0.002、xx, xy, yy) | 6e-3 程度 | 精密版 −9e-5/−9e-5, +2e-5/−5e-5, +3e-4/−7e-4。既定 −1.5e-3/−3.4e-3, +8e-4/−2.5e-3, −2e-4/−3.7e-3 | kick の強さの 1/(1+δ) の誤り (既定にも入る)、偏向磁石・四極・ドリフトの項、追跡の光学 (精密版) |
| (b) ASTRID2 の ∂(ADTS)/∂δ の解析/追跡 | 0.95, 0.28, 4.6 | 精密版 1.000, 0.997, 1.285。既定 0.962, 0.847, 1.971 (数値の正規形: 1.003, 0.999, 1.050) | 同上。既定の差は 4 次の積分器の分が 2 点の差で増幅される |
| (b) FODO の ∂(ADTS)/∂δ の解析/追跡 | (R2 の記録なし) | 1.049, 0.794, 1.045 (数値の正規形: 1.000, 1.001, 1.000) | 解析の残りは未解明 |
| (b) KEK-LS の ADTS (δ = ±0.002、xx, xy, yy) と ∂(ADTS)/∂δ | R2 の記録なし (8 極のフィードダウンなしで 1.4e-2, 2.0e-2, 7.8e-2) | −1.8e-5/+1.1e-5, −7.8e-5/+6.6e-5, −4.0e-4/+3.7e-4。微分の比 1.000, 1.001, 1.004 (数値の正規形 1.000) | 軌道がずれた 8 極の 6 極へのフィードダウン K2 = K3·x0 |
| (c) 閉形式と取り出し (E1・E2≠0 の FODO) | 閉形式なし (1 個のスライスで 1e-5〜8e-5) | 1.6e-8 (K1=0)、ASTRID2 1.8e-6、KEK-LS K1≠0 の 2 台 9e-5・3e-4 | 全体の和 (積分) と比べる。端の閉形式を導出 |
| (c) スライス数への依存 (n = 4 と 64) | 1e-6〜6e-5 | 6e-11 (FODO)、4e-8 (ASTRID2) | 4 次の項を要素の写像から取り出して内部の対の項を引く |
| (d) 非対角 35 調和 | 補正なし | 独立な数値のリー変換と 1.2e-8 | F4' = F4 + ½{R⁻¹W, W} |
| (e) `adts()` と追跡 (FODO, xx xy yy) | −3.6e-4, +4.1e-4, −3.4e-4 (これが既定) | 精密版 −8.8e-5, −1.7e-5, −1.2e-4 | 四極の台形則 + 偏向磁石の 4 次の項 (積分器の端数) |
| (e) ASTRID2 | −2.4e-3, −8.2e-4, −1.9e-3 | −5.8e-5, +2.0e-5, −1.8e-4 | 同上 |
| (e) ASTRID2 MPW | −2.4e-3, −8.4e-4, −2.2e-3 | −1.0e-4, −5.3e-5, −4.7e-4 | 同上 |
| (e) KEK-LS | −1.4e-4, +1.7e-4, +2.6e-3 | −6.9e-7, +1.8e-6, +8.3e-6 | 同上 |

## コマンド

```
set PYTHONUTF8=1
py -3.12 example/v2_3/R6_resolve.py mpw          (数分。追跡 6 回)
py -3.12 example/v2_3/R6_resolve.py delta astrid (追跡 3 回 × 35 s。fodo は 7 s ずつ、kekls は大きなリングの追跡で数十分)
py -3.12 example/v2_3/R6_resolve.py bend         (数秒)
py -3.12 example/v2_3/R6_resolve.py nf           (数秒)
py -3.12 example/v2_3/R6_resolve.py adts         (KEK-LS を含めて 10 s、追跡は R5 の保存値)
```

Python から: `R.adts(lat)` (既定: `fourth='r3'`・`optics='calc'`)、`R.adts(lat, fourth='r6')` (精密版)、`R.adts(lat, delta=0.002, fourth='r6', optics='tracked')` (δ≠0 の精密版)、`nform.adts_delta_derivative(lat, p, 2e-3, method='nform')` (δ 微分)、`rdt2nd.rdt4_normal_form(lat)` (4 次の 35 調和の正規形の値)。

## 使えない範囲・未確認

docs/49 §11.8・§11.9 のとおり。主なもの: δ 微分の解析式の yy (精密版で ASTRID2 28%、FODO で xy が 21% 外れる。原因は未特定。数値の正規形を使う)、K1≠0 の本体の閉形式の 1e-4 (K64 の積分器と推定)、K0≠0 の磁石に E・ROTATE がある場合、共鳴線の上での非対角の補正、`other_atoms` は四極を分割しない、KEK-LS の δ≠0 の解析は 20 秒程度。
