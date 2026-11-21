# KEKLS_correct_cod
種別: 誤差入りの COD・分散の補正。フォルダ: `errors`。元のスクリプト: correct_cod.sad (correct_cod_best/worst/noSX/GCUT1/v03_70 は設定違い)。ラティス: `../lattice/kek_ls_v03_68.sad` (3 GeV、周長 570.72 m、標準セル 20 個)。フラグ: `ON RAD FLUC COD RFSW RADCOD`。

## 何を見る例か
乱数の種 (17・111) の誤差入りのリング RING で、COD 応答行列の特異値分解 (固有値の大きい ICX/2・ICY/2 個の疑似逆行列) により、ステアリング HD (161 台)・VD (121 台) の K0 を 0.5 倍ずつ 6 回かけて、モニター 201 台の COD を補正する。

## スクリプト
`KEKLS_correct_cod.sad` (先頭に `../lattice/kek_ls_v03_68.sad` を貼った自己完結のファイル。台本の本体は `FFS;` から。説明のコメントは台本の先頭にある)。共通モジュールは `Get["../modules/..."]` で読む (台本のあるフォルダをカレントにして走らせる)。

## 結果の読み方
各種で `before correction` (チューン、COD の最大・rms、η のずれ)、`iteration l:` (行列で予測した rms と SAD の rms)、`after correction`。種 111 で COD rms (x, y) が 430 / 640 μm から 42 / 23 μm (gs・pysad とも)、SAD の rms は 6 回で数 μm に収束する。

gs の出力の抜粋 (指定の行だけ。全体は `KEKLS_correct_cod.gs.txt`):
```
===== seed   17.0
before correction: tune   48.5532991311   17.7244474451
after correction: tune   48.5753723079   17.6297862396  EMITX [pm rad]  137.5274131138
===== seed  111.0
before correction: tune   48.5762095890   17.6492502682
after correction: tune   48.5738016146   17.6192299699  EMITX [pm rad]  136.6160939264
```

## gs と pysad の比較
判定: **△** (比べた行 27、最大の相対差 4.9e-01、k64 モードの pysad)。種 111 は補正の各回で gs と pysad が 1e-2 (相対) 以内。種 17 は補正前から違う (位置ずれの閉軌道の解が違う、`KEKLS_errors_cod` と同じ)。

所要時間: gs 1 秒、pysad 37 秒 (この PC。k64 モード)。

## 未確認・限界
- 元は種 115 本 (k = 1〜115)、補正 40 回、毎回 DA (5000 周) を計算して datas.txt に書く。例は種 2 本、補正 6 回 (rms は数回で収束)、DA なし (`doDA = 1` で、周回数 200 の DA)。完全な設定は `seeds`・`niter`・`doDA`・`DAturns` (DA の周回数は 5000 のまま) を変える。見込み: 種 1 個あたり pysad で約 10 秒 (応答行列の準備に約 17 秒)、115 個で gs 約 1 分・pysad 約 20 分 (外挿、DA は別)。
- 補正磁石 HD・VD は同じ要素名の繰り返しなので、既定の `cell calc` では K0 が全部同じ値に揃えられ、個別に効かない (元の台本が前提にしたラティスには、個別の名前があったと思われる)。`USE RING NOEXPAND` と `cell calc noexpand` で成分ごとの K0 が効く (gs・pysad とも確認)。
- 位置ずれは成分ごとの DX・DY に直接入れる (K64 の `DELX`・`DELY` は閉軌道に効かない)。元の `MON ...` と `s!top;` (打ち間違い) は除いた。
- 元の modules_correctoptics.n・modules_write_correctdata.n は correct_cod.sad から読み込まれていない (コメントアウトされた Get のみ)。取り込んでいない。
- 補正後のチューンは設計値から 0.005 動く (補正磁石の K0 の寄与)。チューンの補正は含めていない (元の correct_cod.sad にも無い)。
