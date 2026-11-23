# fn_Fork_Shared (Fork, OpenShared, Shared, SharedSize)
種別: 関数 / フォルダ: 09_io_system / FFS 内の SADScript (ラティスなし)。gs は SAD 1.4.4.2k64、pysad は `run_sad_string`。ファイルは `TemporaryName[]` か実行ディレクトリの作業ファイルで、最後に `DeleteFile` で消す。

## 何を見る例か
`Fork[]` はプロセスを 2 つに分け、子で 0、親で子の pid を返す。`OpenShared[n]` は共有メモリ、`Shared[s]` で読み書き、`SharedSize[x]` は値の大きさ。子が共有メモリに書いた値を親が読む。

## スクリプト
```
s = OpenShared[64];   Shared[s] = {1, 2, 3};   SharedSize[{1,2,3}];
pid = Fork[];
If[pid == 0, Shared[s] = {10, 20, 30}; Exit[0], Pause[0.5]; Print["parent ..."]];
Print[Shared[s]];   Close[s];
```

## 結果の読み方
gs: 親が `{1,2,3}` を書き、子が `{10,20,30}` に書き換えて終了し、0.5 秒待った親が `{10,20,30}` を読む。`SharedSize[{1,2,3}]`=32 バイト。親の `Fork[]` は実数 (pid、正)。

## gs と pysad の比較
△ **Windows の pysad は子プロセスを作らない (`.pysad.txt` はこの結果)。** `SharedSize[{1,2,3}]` は 32 になり gs と一致する。`Fork[]` は Windows に fork が無いので、v2 は `???System::error: Fork[] is not available on this platform` を出し (k64 モードは実数 1 を返す)、子の分岐を通らず、共有メモリは親のものだけで `{1,2,3}` のまま (gs は `{10,20,30}`)。Linux・WSL・macOS では v2 の `Fork[]` が `os.fork` で子を作る (`docs/28_v2_vs_k64.md`。記録は WSL で確認)。**この例を WSL の pysad では再実行していない** (WSL の Python に numpy が無い)。要確認。

## 未確認・限界
- 子の標準出力は親と共有されるので、子で `Print` すると行が重なる。子では何も表示していない。
- 待ちは `Wait` ではなく `Pause[0.5]` (SADHelp の `Wait` は未使用)。
- `Wait[pid]` を使うと gs が 120 秒止まるので、この例では使わない。
