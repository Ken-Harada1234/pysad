# pysad/data: 同梱のデータファイル

pysad が持つ小さなデータファイル (pip で入る。例 (example/) や文書は pip に入れない)。

| ファイル | 内容 |
|---|---|
| `kek-ls_v03_12.lat` | KEK-LS (3 GeV、20 セル。ESRF 型の HMBA) のラティス (SAD の入力形式)。例と試験が使う (「大きなリング A」)。`pysad/tests/kek-ls_v03_12.lat` と同じ内容 |

場所は `pysad.datafiles.data_path('kek-ls_v03_12.lat')` で取れる。例の台本 (`example/` の `*.sad`) は、このラティスを `pysad/tests/kek-ls_v03_12.lat` という相対パスで読むので、リポジトリの根 (clone) をカレントにして走らせるとそのまま動く (例は pip には入れていない)。
