# example/astrid の履歴 (整理で md から取り除いた記述)

元の文書の文脈つきで、取り除いた文章をそのまま入れてある。

## 共通: 「gs と pysad の比較」の冒頭の v2 の注記 (2 文書)

各文書の「gs と pysad の比較」の冒頭にあった段落 (全文書で同文)。現在の md では「**v2 モードと k64 モード**」の見出しの段落に言い直した (日付・再実行・出力の作り直しの記述を除いた)。

```
**v2 の注記 (2026-11-17、v2.0.0 の再実行)**: `.pysad.txt` は既定の v2 モードの出力に作り直した。保存してあった旧い出力は `PYSAD_COMPAT=k64` の出力と同じ (確認済み)。v2 は K64 のバグを直したので (`docs/47_v2_changes.md`。主に `Emittance[]`・EMIT・追跡の Jx の端の項 B75 と `MomentumCompaction` B77)、この例の出力のうち `Emittance[]` の表示 (εx が数 % 大きい、εz・減衰率・α・バケット高さなど) と、それに依る量が gs と違う。以下の「一致」の判定は k64 モードの結果。gs と同じ値が要るときは `PYSAD_COMPAT=k64` を付ける。
```

対象の文書:

- `astrid/tune_survey/tune_survey_ncell_ideal_ASTRID2.md`
- `astrid/tune_survey/tune_survey_ncell_smallerror_ASTRID2.md`

