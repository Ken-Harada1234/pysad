# 13_wake の検証の詳細 (gs との比較)

docs/13_wake.md の「検証」の詳細。旧 docs/23_wake.md §4・§7 の検証の表と、旧 docs/38_wake_guide.md §8 を、元の文脈つきで移した。

## 追跡 (旧 23 §4)


SAD K64 (WSL の gs) と同じ台本を走らせて比べた。粒子 200〜400 個の成分ごとの総和と 3 粒子の値が、10 桁一致する。

| ケース | 内容 |
| --- | --- |
| 縦だけ、横だけ、縦と横 | マーカー 1 個、PBUNCH = 1E10 と 1E12 |
| 2 か所のウェイク | QUAD (DX≠0、同名 2 個なので `QF.1`) と MARK に別の表 |
| NBUNCH = 3 | 3 つのバンチ列 (z を 0.5 m ずつずらす) |
| 粒子損失 | LOSSAMPL = 0.0025 で 107 個が落ちる (np が途中で減っても一致) |
| flag=0 の粒子 | 7 の倍数を最初から落とす。NP = 全粒子数 (400) が最良。wake なしでも出る差 (初めから落ちた粒子の出力座標、3e-4) 以外の差はない |
| 同じ z の粒子 | 5 個ずつ同じ z (`spsort` を写して 0 差) |
| ROTATE≠0 の要素 | QF ROTATE=0.3、QD ROTATE=−0.2 (DX、DY も) |
| CAVI | VOLT=0 で L=0 と L=0.4 (フリンジあり)、DPHI≠0 |
| 加速のある CAVI | VOLT=1E5 (L=0)、3E7 (L=0.6、フリンジあり・なし)。wake なしと wake あり (厚い空洞と薄い空洞)。ndiv>1 になる |
| RF のある MULT | TRPT の加速、薄い MULT、K0・K1・K2・K3、DX・DY・ROTATE、DVOLT、DPHI、FREQ、フリンジ。wake なしと wake あり |
| MULT (VOLT=0) の wake | MULT に wake |
| SOL の領域 | 内側の QUAD (ROTATE=0.2) と終端の SOL にウェイク。wake なしの領域の追跡も一致 |
| 表の最後の行 | 最後の行を変えても結果が変わらない |


## 光学 (旧 23 §7 の検証)


gs と同じ台本の CALC の出力を比べた ([`test_wake_optics.py`](../pysad/tests/test_wake_optics.py)): 表は文字列で完全一致、中心の軌道 (`Twiss["DX"..]`) は相対 2e-8 以内 (z は絶対 3e-14 以内)。

| ケース | 内容 |
| --- | --- |
| 2 軌道 (縦+横) | `FIT $$$ DX DDP` の表。軌道ごとの AX BX NX .. 全部の列が一致 |
| 3 軌道 (縦だけ、横だけ、縦+横) | 奇数 (nfam1=1−nfam) の縦のキックの添字のずれ |
| 4、5 軌道 | 偶数・奇数、QUAD (DX・DY) の wake |
| NBUNCH=3、SIGZ=0.01、2 か所の wake | 表の使い回し (z がほぼ同じ) の分岐 |
| 28 個形式の軌道 | Twiss (BX>0) と NX の引き継ぎ |
| WAKEOPT なし | InitialOrbits だけ (独立な伝搬) |


## 旧 38 §8 gs との一致の度合い


数値は [`docs/13_wake.md`](13_wake.md) §4・§7 と試験から。使った条件: 追跡はエネルギー 1.7 GeV 前後の小さなライン (INS、NORAD NOFLUC)、粒子 200〜400 個、PBUNCH=1E10〜1E12。実際のリングでの検証はしていない。

- 追跡: 成分ごとの総和と 3 粒子の値が 10 桁一致。ケースは縦だけ・横だけ・両方、2 か所のウェイク、NBUNCH=3、粒子損失 (LOSSAMPL=0.0025 で 107 個)、flag=0 の粒子、同じ z の粒子、ROTATE≠0、CAVI (VOLT=0 と加速あり、フリンジあり・なし)、RF のある MULT (K0〜K3、DX・DY・ROTATE、DVOLT、DPHI、FREQ)、MULT (VOLT=0)、SOL 領域。
- 光学: 条件の表は文字列で完全一致、中心の軌道 (`Twiss["DX"..]`) は相対 2e-8 以内 (z は絶対 3e-14 以内)。2〜5 軌道、NBUNCH=3 と表の使い回し、28 個形式、WAKEOPT なしを含む。
- 試験: `py -3.12 -m pytest pysad/tests/test_wake.py pysad/tests/test_wake_optics.py -q` で 32 通過 (本文書の作成時に確認)。SAD が無くても通る (参照値は試験に書き込み済み)。


## 旧式 WAKE コマンドの一致 (旧 38 §5、旧 23 §2c)

追跡は `WakeFunction[…] = …` で同じ表を与えた場合と同じ結果 (gs の参照値と 2e-9、400 点・2 台の CAVI で 1e-12。6 点・5 点の表で gs の参照値と 2e-9、400 点の表・2 台の CAVI で WakeFunction と 1e-12)。
