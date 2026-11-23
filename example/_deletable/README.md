# 削除してよい作業ファイルの置き場

整理のあとに不要になった作業ファイルを、削除してよい形で分けて置いた場所です。例の説明・試験はここを参照していません (例外は下の表)。

| フォルダ | 内容 | 参照するもの |
|---|---|---|
| `lma_corrected_kekls/` | KEK-LS に誤差を与え、COD を補正し、補正後の LMA を求めた途中の出力 (種・補正の有無・周回数・放射の有無を変えた 40 本。`lma_corrected_run.py` の標準出力を保存したもの)。サマリ `15_lma_touschek/lma/lma_corrected_kekls_v03_12_t982_summary.out.txt` だけを例のフォルダに残した | `15_lma_touschek/lma/lma_corrected_summary.py` が、v03_12・982 周の 9 位置の出力を読んでサマリを作り直す (消すとサマリは作り直せない) |
