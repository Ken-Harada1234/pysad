# fn_SurvivedParticles
種別: 関数 / フォルダ: 05_tracking / ラティス: `00_lattice/FODO_ring.sad` (3 GeV) / フラグ: `LOSSAMPL=0.05`

## 何を見る例か
`SurvivedParticles[x]` は、`TrackParticles` の戻りのビーム (7 列) から、flag=1 の粒子の列だけを取り出す。全滅のときは長さ 0 のリスト 7 本。

## スクリプト
```
r = TrackParticles[{1, beam}, 1, 1, 20];
s = SurvivedParticles[r[[2]]];     ! 7 x (生存数) の表
r3 = TrackParticles[{1, {{0.06},{0},{0},{0},{0},{0},{1}}}, 1, 1, 2];
SurvivedParticles[r3[[2]]]         ! 全滅: {{},{},{},{},{},{},{}}
```

## 結果の読み方
- 5 粒子 (うち 1 つは flag=0 で入力) のうち 2 粒子が生き残り、`Dimensions[s] = {7, 2}`、$x=\{0.001852, -0.045458\}$、flag={1, 1}。
- $x_0=0.06$ m ($>$ LOSSAMPL=0.05) の 1 粒子ではすぐ落ちるので、戻りは `{{},{},{},{},{},{},{}}`。

## 物理・解析式
なし (リストの操作)。生存数の見積もりは [fn_TrackParticles_loss](fn_TrackParticles_loss.md) の「物理・解析式」を参照。

## gs と pysad の比較
○ 一致 (全 3 行が同じ)。

## 未確認・限界
- 損失の位置と周は戻りに入らない (LOSSMAP を使う: [flag_LOSSMAP](flag_LOSSMAP.md)。pysad は未対応)。
- POL ON の 9〜11 列のビームでの動作は見ていない。
