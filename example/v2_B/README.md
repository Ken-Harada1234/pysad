# example/v2_B: v2.0.0 (担当 B) の走る例

pysad v2.0.0 で直した・足した機能のうち、特殊関数、SADScript、FFS の小項目の例です。説明と根拠は [`docs/history/v2_changes_B.md`](../../docs/history/v2_changes_B.md)。

実行 (リポジトリのルートから。既定の互換モードは v2):

```
cd example/v2_B
py -3.12 ../../run_sad.py B1_special_functions.sad
```

K64 (gs) と同じ動作で走らせるときは、環境変数 `PYSAD_COMPAT=k64` を付けます。

| ファイル | 内容 | 見どころ |
|---|---|---|
| `B1_special_functions.sad` | K64 が誤値・NaN・INF・無応答になる特殊関数 (B43〜B47、B53、B54)、非整数次の PolyLog の大きい引数、複素引数の Bessel・超幾何 | 各行のコメントに真値と K64 の値。`InverseErf[.3]` だけはモードで変わる (v2 .272462714727、k64 .272462711811) |
| `B2_package_context.sad` | `BeginPackage`・`Begin` の文脈 (名前空間)。パッケージのファイルを書き出して `Get` で読む | `Stat`Mean2[…]` で呼ぶ。`` `Sigma `` で大域に公開。大域の `Counter` はパッケージに壊されない |
| `B3_shared_system.sad` | 共有メモリ、`Parallelize`、`System`、`Fork` | `System` は既定で実行しない (`PYSAD_ALLOW_SYSTEM=1` で許可)。`Fork` は Linux・WSL・macOS だけ |
| `B4_ffs_v2.sad` | MAIN の `OFF RFSW`、`INDEX` 属性、`MON`、`ORG` の負数、`NOKEEPEXP`、MOMENTUM 未設定の警告 | 小さな試験用リング (周長 7 m。MOMENTUM をわざと書いていないので 30 GeV で追跡し、警告が出る) |

`*.pysad.txt` は v2 モードの pysad の出力 (Windows、2026-11-17)。`B3` は Windows の出力なので `Fork[]` は使えないと表示されます
(WSL では `child wrote: written by child`、`PYSAD_ALLOW_SYSTEM=1` を付けると `System の戻り値 = 0`)。

## 未確認・限界
- 例の数値は上の小さなリングと単体の式のものです。実機のラティスでは確かめていません。
- `B3` の `System` の結果は OS のシェル (Windows は cmd.exe) に依ります。
- 文脈の K64 との差 (System` の範囲、表示) は [`docs/history/v2_changes_B.md`](../../docs/history/v2_changes_B.md) §2.4。
