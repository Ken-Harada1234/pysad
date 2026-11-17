# fn_TouschekLifetime
種別: 関数 / フォルダ: 06_collective / ラティス: `00_lattice/FODO_ring_V4MV.sad` (RF 4 MV、3 GeV) / フラグ: `RFSW; RAD; FLUC; INTRA; MINCOUP=0.01; PBUNCH=1e9`

## 何を見る例か
`INTRA` を ON にして `EMIT` を実行すると、バンチ内散乱 (IBS) の拡散と同時に、Touschek 散乱の損失率の表 (`TouschekTable`) が作られる。`TouschekLifetime[nx, ny, nz]` は、その表から寿命 [秒] を補間して返す。`TouschekLifetime[Infinity, Infinity, nz]` は運動量アクセプタンス $nz\cdot\mathrm{SIGE}$ での寿命。

## スクリプト
```
RFSW; RAD; FLUC; INTRA;
MINCOUP = 0.01;        ! 垂直エミッタンスの下限 = MINCOUP (EMITX+EMITY) (理想リングは EMITY≈0)
PBUNCH = 1e9;
EMIT;                  ! EMITX, EMITY, SIGE, SIGZ と TouschekTable
TouschekLifetime[Infinity, Infinity, nz]       ! nz = 4, 6, 8
```
- `Emittance[]` でも表は作られるが、`SIGE` などの変数が更新されるのは `EMIT` (コマンド)。
- 運動量アクセプタンスは SIGE ($\sigma_\delta=1.1914\times10^{-3}$) の `nz` 倍。RF バケットの高さは約 1% なので、`nz=8` ($9.5\times10^{-3}$) が RF の限界に近い。

## 結果の読み方
EMITX=1.735775e-7、EMITY=2.16e-21 (MINCOUP による下限が IBS の計算に使われる)、SIGE=1.1914e-3、SIGZ=16.06 mm。

| nz | 運動量アクセプタンス | 寿命 $\tau$ [s] |
|---|---|---|
| 4 | $4.766\times10^{-3}$ | $1.381730\times10^6$ |
| 6 | $7.148\times10^{-3}$ | $3.908993\times10^6$ |
| 8 | $9.531\times10^{-3}$ | $6.981234\times10^6$ |

粒子数を 2 倍 (`PBUNCH=2e9`) にして `EMIT` し直した nz=6 の寿命: $1.954497\times10^6$ s、比 `tau(1e9)/tau(2e9)=2.000000`。

## 物理・解析式
Touschek 損失率は $1/\tau\propto N/(\sigma_x\sigma_y\sigma_z)\cdot\delta_{acc}^{-2}\times C(\varepsilon)$ の形 (Bruck の式。$\delta_{acc}$ は運動量アクセプタンス、$C(\varepsilon)$ は $\varepsilon=\delta_{acc}^2\beta^2\ldots$ の緩い関数)。式の係数は検証していないが、次の 2 つの性質は確認できる。

| 性質 | 式 | 結果 (gs = pysad) |
|---|---|---|
| $1/\tau\propto N$ | $\tau(N)/\tau(2N)=2$ | 2.000000 |
| $\tau\propto\delta_{acc}^2\times$(緩い関数) | $\tau(6)/\tau(4)\gtrsim(6/4)^2=2.25$、$\tau(8)/\tau(4)\gtrsim4$ | 2.829、5.052 (ずれは $C(\varepsilon)$ の変化で 1.26 倍) |

2 つめは桁と傾向の確認だけ (係数は確かめていない)。

## gs と pysad の比較
**v2 の注記 (2026-11-17、v2.0.0 の再実行)**: `.pysad.txt` は既定の v2 モードの出力に作り直した。v2 は `Emittance[]`・EMIT・追跡の Jx の端の項 (B75) と `MomentumCompaction` (B77) を直したので、この例の pysad の値のうち εx (約 2.4% 大きい)・εz・減衰率・σδ・バンチ長・α・バケット高さ・放射によるチューンシフトが gs と違う (K64 のバグを直したため。`docs/47_v2_changes.md`)。以下の表と本文の「pysad」の値と「一致」の判定は、v1 と同じ `PYSAD_COMPAT=k64` での結果 (旧い出力と同じことを 2026-11-17 に確認)。gs と同じ値が要るときは `PYSAD_COMPAT=k64` を付ける。
Touschek 寿命は EMIT のエミッタンスに依るので、v2 では全体に約 2% 長い (例: 表の先頭の寿命 1023084 → 1044275 s)。

○ 一致 (全 10 行が印字桁まで同じ)。

## 未確認・限界
- 絶対値 ($10^6$ s 台) は、このラティス (ビーム電流・ビームサイズは人工) での値で、実機との比較はしていない。`nx`・`ny` が有限の場合 (振幅アクセプタンス) と、`TouschekTable` の中身も見ていない。
- 垂直エミッタンスは MINCOUP の下限 (0.01) で決まるので、寿命はこの値に依存する。
