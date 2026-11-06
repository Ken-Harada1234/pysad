# 作業ログ: 粒子が落ちる周の乱数消費 (RAD+FLUC) を SAD と一致させる

担当ファイル: `pysad/track6n.py` (dapert.py は変更不要だった)。作業用スクリプトは `<work-scratch>/`。
日付: 2026-10-15

## 問題

RAD+FLUC の追跡で、粒子が落ちない限り pysad (track6n) の乱数列は SAD (NPARA=1) と最終 seed まで一致していたが、
粒子が落ちる周の中で乱数の消費数がずれていた (bench_level16l_j22: t1 一致, t2 で SAD 1515505817 に対し pysad 3673544217)。

## 調べ方

1. SAD で `TrackParticles[beam, dest, 1, nturn]` の dest を変えながら、各 dest での `SeedRandom[][[2]]` と座標を印字
   (`<work-scratch>/j22_dest.sad` 1 周目 20 要素刻み、`j22_dest2.sad` 2 周目 20 要素刻み、`j22_dest3.sad` 940〜990 を 1 要素ずつ、
   `j22_dest4.sad` は `LOSSAMPL=1e100` にして tapert のクランプ (limitnan: |x|>3·LOSSAMPL → ±3·LOSSAMPL, NaN → +3·LOSSAMPL) を
   外した実値)。`TrackParticles[beam, dest, 1, nturn]` は tturn(1, nlat) を nturn−1 回 + tturn(1, dest) で、dest の手前までの
   座標を返す (末尾でも tapert が呼ばれるのでクランプは入るが乱数は消費しない)。
2. pysad 側は `<work-scratch>/j22_drv.py` で track6n の kernel を要素ごとに呼び、要素ごとの seed と座標を印字して突き合わせた。
3. 要素の中は `<work-scratch>/j22_quad.py` (QFA6H.22 = tsolqur) と `<work-scratch>/j22_bendi.py` (BPI6.22 = tbendi) で、
   track6n と同じ式を純 Python で 1 スライスずつ追った。
4. gfortran (-Ofast: SAD の sad.conf の FOPT) の NaN/Inf の扱いは `<work-scratch>/ftest/nantest2.f`, `normtest.f` を
   SAD と同じフラグでコンパイルして確認した (コンパイル時定数畳み込みを避けるため値はコマンド引数から読む)。

## 分かったこと (j22, 2 周目)

粒子は 2 周目の要素 963 (QFA6H.22) で発散し始め、973 (DR08H.33 の先頭の tapert) で落ちる。ずれの原因は 4 つあった。

1. **sqrt1 / pxy2dpz の NaN**: SAD の mathfun は `sqrt1(x)=x/(1+sqrt(max(1e-200, 1+x)))`、
   `pxy2dpz=-x/(1+sqrt(max(1e-200, 1-x)))` (小さい x は級数) で、px²+py²>1 でも NaN にならず有限値 (≈ −x) を返す。
   track6n は `a/(1+sqrt(1+a))` だったので px²+py²>1 の瞬間に NaN になり、以後の tdusrn が anp=NaN で 10000 回ループ
   (乱数 20001 個) していた。→ SAD と同じ式に置き換えた (級数の分岐も含む。既存ベンチの差は変わらず ≤4.7e-11)。
2. **tapert のタイミング (la の数え方)**: tturn.f では `la` は要素ループの末尾 `1020 la=la-1` で減るが、DRIFT は `go to 1011`、
   MARK/MONI は `exit` でここを通らないので la を減らさない。つまり tapert は「DRIFT/MARK/MONI 以外の要素 15 個ごと」。
   track6n の turn() は全要素で減らしていた。さらに parser が L=0 の BEND (VD/HD, ROTATE 付き) を Drift にするが SAD では
   icBEND なので数える → `Prog.la_count` (Drift は `_sad_type` が 'BEND' のときだけ数える) を追加して turn() で使う。
3. **asinz(NaN)**: SAD の asinz は `asin(min(1, max(-1, x)))`。-Ofast の gfortran では max(-1, NaN) が −1 になり
   asin(−1) = −π/2 が返る (ftest で確認)。tradkfn では theta<0 → anp<0 → p0=exp(−anp)>1 → 光子なし (乱数 1 個)。
   track6n の fmax/fmin (a if a>b else b) では NaN が返り anp=NaN → 10000 回ループになっていた。→ asinz に NaN の特別扱いを追加。
4. **norm2 の NaN**: tradkfn の `ppa=norm2([ppx,ppy,ppz])` は gfortran のインライン展開 (scale/ssq 法) で、成分に ±inf が
   2 つあると (inf/inf)² で NaN になる (normtest.f で確認: norm2([inf,-inf,x]) = NaN、norm2([inf,1,2]) = inf)。
   track6n は `sqrt(ppx²+ppy²+ppz²)` = inf → asinz(inf)=+π/2 → anp=+95 で光子が 102 個出ていた (BPI6.22 の 1 スライス目)。
   SAD は NaN → −π/2 → 光子なし。→ `norm2_3()` で gfortran の算法を再現。

gfortran -Ofast の MAX/MIN の NaN の扱いはサイトごとに違う (maxsd のオペランド順がレジスタ割付で決まる: ftest で
max(-1,NaN)=−1、max(smin,NaN)=smin だが max(gmin,NaN)=NaN、min(0.9999,NaN)=NaN) ので一般には再現できない。
ただし上の 1 の修正で SAD と同様に NaN がほとんど発生しなくなり、残る NaN は asinz/norm2 の 2 箇所の経路に限られた。
比較演算も -Ofast では NaN で真になるものがある (`v==v` が T、`v/=0` が F) が、tdusrn のループ回数には影響しない
(`p<=pn`, `pn>cutoff` は NaN で偽、ftest で確認: tdusrn(NaN) は 10000 回, tdusrn(+Inf) も 10000 回, tdusrn(−Inf) は 0 回)。

修正後の要素ごとの seed (2 周目): 963 QFA6H.22 3977592969, 965 BPI6.21 2431828505, 967 BPI6.22 2410509385 (28 個 = 28 スライス × 1),
970 QFA8H.21 4182837073, 972 QFA8H.22 1515505817, 973 の tapert で損失 — すべて SAD と同じ。
座標も SAD と同じ (QFA6H.22 の出口で x=2.6100825 (SAD 2.61008255947741, pysad 2.61008254167: 発散の増幅で 8 桁目から違うが
seed と落ちる場所は同じ)。

## 結果

- `<work-scratch>/l16l_check.py`: bench_level16l_j22/j25/j27 × t1/t2/t4/t10 の 12 件すべてで最終 seed が SAD と一致
  (j25 は 4 周目、j27 は 3 周目で損失。"NG" と出る 2 件は SAD の表示 3073029016.999999 を int() で切り捨てただけで一致)。
- `bench_level12_track_fast.py` (NORAD/RAD/FLUC の既存ベンチ): 変化なし (≤4.7e-11, 乱数一致)。
- bench_level16f (RAD FLUC の DA, 3 行 × 51 粒子 × 100 周): 下記。

## track6n.py の変更点

- `sqrt1`, `pxy2dpz`: mathfun の直訳 (級数 + max(1e-200, ·))。
- `asinz`: mathfun asinz の直訳 + NaN → −π/2。`asinz_e`, `asinz_l` は asinz の別名に。
- `norm2_3`: gfortran norm2 の再現。`tradk_core1` の ppa に使用。
- `Prog.la_count`、`turn()` の la の減らし方。
