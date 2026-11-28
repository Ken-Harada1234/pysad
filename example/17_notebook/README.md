# 17. Jupyter の notebook から SAD を実行する

SAD のラティスファイルとコマンドファイルの中身を、notebook のセルにそのまま貼って実行する例です。説明と限界は [`docs/12_pysad_api.md`](../../docs/12_pysad_api.md) §11b。

| ファイル | 内容 |
|---|---|
| [`sad_in_notebook.ipynb`](sad_in_notebook.ipynb) | 5 つの使い方の例 (出力つきで保存してあります): (1) `%%sad lattice` でラティスを貼る、(2) `%%sad` でコマンドを貼って実行する、(3) `sad.value()` で結果を Python に取り出して matplotlib で描く、(4) `%%sad --continue` で前の状態の続きを流す、(5) マジックを使わずに `pysad.run()` に文字列で渡す |

## 動かし方

```
pip install pysad-k64 ipykernel jupyterlab
```

notebook の最初のセルで `%load_ext pysad.notebook` を実行します。そのあと、セルの先頭に `%%sad lattice` (ラティス) または `%%sad` (コマンド) を書いて、中身を貼ります。

## 確認した範囲と限界

- 出力は、Windows・Python 3.12・IPython 9.19.1・ipykernel 7.4.0・nbclient 0.11.0 で、本物のカーネルを起動して notebook を実行して作りました。JupyterLab・VS Code の notebook・Colab・他の OS は未確認です。
- 例のラティスは `example/00_lattice/FODO_ring.sad` の中身 (3 GeV、24 セル、周長 105.6 m) です。(4) は QF の K1 を 0.7 から 0.75 に変えた計算で、水平チューンは 42.07 から 46.19 に変わりました。この値の物理的な意味は、この例では見ていません (使い方の例です)。
- (5) の格子は、(1) と別の簡単な FODO リング (1 セルに偏向磁石 1 台、24 セル) です。
- `--continue` で流せるのは FFS 部の文だけです。要素・LINE の定義は、新しい `%%sad` のセルに書きます。
- Tk の GUI は notebook では動きません。SAD の描画関数 (Plot 系) の notebook 内の表示は未確認です。
- notebook を開き直したときは、`%load_ext pysad.notebook` と、ラティスのセルを、もう一度実行してください (保存したラティスは、カーネルの中にあります)。
