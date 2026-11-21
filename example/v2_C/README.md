# example/v2_C: v2.0.0 (担当 C) の走る例

pysad v2.0.0 で直した・足した機能のうち、追跡・集団効果・要素の例です。説明と根拠は [`docs/history/v2_changes_C.md`](../../docs/history/v2_changes_C.md)
(MEASURE・TRCOD は `v2_changes_C_measure.md`、空間電荷は `v2_changes_C_spch.md` に例がある)。

実行 (既定の互換モードは v2):

```
cd example/v2_C
py -3.12 ../../run_sad.py C1_undulator.sad
```

K64 (gs) と同じ動作で走らせるときは、環境変数 `PYSAD_COMPAT=k64` を付けます (C1・C2・C3 は k64 では動きません: K64 に無い・動かない機能)。

| ファイル | 内容 | 見どころ |
|---|---|---|
| `C1_undulator.sad` | UND の磁場 `FBX` (B37)。リングにアンジュレータ (1 m、20 周期、1.8 T) を入れる | 垂直チューンのずれ、放射損失の増分が解析式 (2/3) r_e γ³ (F0²/2) L E と 0.04% で一致、減衰率の和 = 2 U0/E |
| `C2_phsrot.sad` | PHSROT (B38)。リング 1 周を位相の回転・減衰・拡散で表す | CALC の Twiss・チューンが入力どおり、Emittance[] が b/(2d − d²)、FLUC つきの追跡の分布がそれに合う (1000 粒子の統計誤差で 5〜10%。2 万粒子では 1%) |
| `C3_sol_region.sad` | ソレノイド領域の中の BEND・MULT・CAVI・APERT (B69) | CALC (INS) の軌道と追跡が全桁一致、APERT で落ちた要素が LOSSMAP に出る、リングの Emittance[] |
| `C4_loss_and_taper.sad` | LOSSMAP の 2 行、最初の tapert (B33)、RADTAPER | 落ちた要素と周。RADTAPER で閉軌道が −79 μm から −0.2 μm に戻り U0 が 0.1% 増える (gs と一致) |
| `C5_radlight.sad` | RADLIGHT の軌道 (B65) と RadiationField | 軌道の終点が幾何 (GX GY GZ) と一致、偏向磁石の放射場のピークが Liénard-Wiechert の式と 0.3% で一致 |
| `C6_misc.sad` | DELCOR (B40・B62)、ResponseMatrix の TRPT、旧式の WAKE コマンド (B66)、MULT のスキュー成分 | CAVI・MARK にも DX が入る、応答に sqrt(p(補正磁石)/p(モニタ)) が掛かる、`WAKE L 名前 s w … ENDWAKE` |

`*.pysad.txt` は v2 モードの pysad の出力 (Windows、2026-11-17)。

## 未確認・限界
- 例の数値は上の小さなリング・輸送路のものです (3 GeV、周長 31〜106 m。C5 は 10 MeV)。実機のラティスでは確かめていません。
- C1: `SLICE` は `POLE` の 20 倍以上にします (刻みが周期の整数倍だと磁場が見えません)。
- C2: 追跡の統計は 1000 粒子・3000 周 (減衰時間の 3〜9 倍) です。
- C3: 領域内の BEND は直線の座標のままの一様な双極磁場です (ANGLE でも K0 でも同じ。幾何は曲がりません)。
- C5: RADLIGHTDIV (1/γ の角度あたりの点数、既定 16) を増やすとピークの誤差が減ります (64 で 5e-5)。

## 担当 C2 の追加 (2026-11-17)

説明と根拠は [`docs/history/v2_changes_C2.md`](../../docs/history/v2_changes_C2.md) (空間電荷は `v2_changes_C2_spch.md`、MEASURE 系は `v2_changes_C2_measure.md`、
WAKE の表・DELCOR・台本の読み方は `v2_changes_C2_misc.md`)。`*.pysad.txt` はどれも v2 モードの pysad の出力。

| ファイル | 内容 | 見どころ |
|---|---|---|
| `C2a_sol_tilt.sad` | ずれて傾いたソレノイド (GEO=1 の SOL の `DX DY DPX DPY`)、端の放射 `F1`、幾何の fit | CALC の閉軌道は領域の外で 0・中はソレノイドの系、設計軌道の粒子は 1 周して 0 に戻る、F1 による U0 の増分 2.11 eV = (2/3) r_e γ³ (bz r/2)²/F1 の和、領域の中の四極で GZ を合わせる fit |
| `C2b_sol_rf.sad` | SOL 領域 (BZ = 1 T) の中の CAVI・RF つき MULT の TRPT の加速、AUTOPHI、wake | CALC (TRPT) の軌道 = 追跡、AUTOPHI で設計軌道 (DZ = 3 mm) が公称の加速を受ける (DDP 7.4e-3 → 2e-10) |
| `C2c_spin.sad` | スピン追跡 (POL + RAD): 平面リング + ソレノイドのスピンチューン、領域内の BEND、UND | cos μ = 2 (cos π aγ · cos φ/2)² − 1 と 5e-6 で一致、スピンと運動量の角 = aγ × 偏向角 |
| `C2d_taper_optics.sad` | RADTAPER の閉軌道 (収束)、CALC6D の RADTAPER、CalculateOptics の TRPT の加速、加速区間の分散 | 軌道 −2.1732e-7 (gs は既定の CODCONV で −2.1977e-7、CODCONV=1e-4 で −2.1732e-7)、CalculateOptics = CALC |
| `C2s_mult_radius.sad` | MULT の RADIUS の空間電荷 (SPAC)。四極を MULT で書くと全周に入る | 四極の中の KV ビームの包絡が包絡方程式と 0.5% で一致、リングのチューンシフトが SPAC・PSPAC・WSPAC・Laslett の式で一致 (SPAC・PSPAC は粒子の雑音で 5%) |
| `C2s_wspac.sad` | WSPAC (ビーム行列のガウス分布の力) と SELFCOD、`WSPACBEAM` | 追跡のチューンシフトが Laslett の式と 1% で一致、gs (K64) とも 1%。平らなリングはエラーで知らせる |
| `C2m_trcod_reference.sad` | TRCOD FIXCOD / MVCOD / STNCOD (基準の軌道と差 `Twiss["DDX"]`)、ORBITJITTER、TRCOD ALIGN | v2 だけ (k64 では「未対応」) |
| `C2m_putget_measure.sad` | PUTMEASURE / GETMEASURE (測定値をファイルで受け渡し、fit の目標値にする) と LWAKE ON の MEASURE | v2 だけ。実行すると `C2m_measure.dat` ができる (消してよい) |
| `C2w_wake_table.sad` | 旧式の WAKE コマンドの表を行をまたいで書く | WakeFunction と同じ結果、エラーの例つき |
| `C2w_delcor_rms.sad` | DELCOR の rms = delx (v2 の規格化。k64 は 1/√2) | 相関長 5 m と 5 km |
| `C2w_main_script.sad` | MAIN の代入 (式・単位・関数) と `FFS USE RING;` | gs の出力と全行同じ (モードによらない) |

### 未確認・限界 (C2 の分)
- C2a〜C2d は 0.5〜3 GeV の小さなリング (周長 91 m)・輸送路 (6 m)、BZ = 1〜3 T、空洞 20〜200 MV の例です。実機のラティスでは確かめていません。
- C2a: 領域の中の四極・BEND・MULT の軸はソレノイドの軸です (ソレノイドを傾けると一緒に傾きます)。SOL の DX DY DPX DPY は GEO=1 の SOL に書きます。
- C2b: wake は TRPT のときだけ働きます (従来どおり)。AUTOPHI の意味 (位相の基準 = 設計軌道の到着時刻) は K64 のソースからの読みで、SADHelp に説明がありません。
- C2c: スピンは RAD が ON のときだけ動きます (K64 と同じ)。Emittance[] の POL の量は、軸を外れたソレノイドの中では K64 の式のまま (スピンチューンで 2e-4 の差) です。
