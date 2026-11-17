"""v2 の SADScript まわり (担当 B): B18 PolynomialFit、B71 For、B70 MessageName、B61 複素 FindRoot、B48 FixedPointList、B52 ベクトル独立変数の Fit、
BeginPackage/Begin の文脈、Fork/Shared、System、Sleep。

根拠: 解析解・numpy の最小二乗・定義式。文脈は gs (K64) の実測 (2026-11-17。local の ctx1.sad) に合わせた。
k64 モードで K64 の動作が残ることも確かめる。"""
import io
import math
import os
import re
import sys
import tempfile

import numpy as np
import pytest

_NUM = r'[-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?'


def _run(body, cwd=None):
    from pysad.ffscmd import run_sad_string
    buf = io.StringIO()
    run_sad_string('FFS;\nPageWidth=1999;\n' + body + '\nExit;\n', out=buf, cwd=cwd)
    return [l.rstrip() for l in buf.getvalue().splitlines()]


def _nums(s):
    return [float(t) for t in re.findall(_NUM, s)]


def _line(lines, tag):
    for l in lines:
        if l.startswith(tag):
            return l[len(tag):]
    raise AssertionError('no line %r in %r' % (tag, lines))


# ---------------------------------------------------------------- B18 PolynomialFit
_PF = 'Print["P ",PolynomialFit[{{-1,1},{0,0},{1,1},{2,4.1},{-2,4}},2][[1]]];'


def test_b18_polynomialfit_mean_zero_v2(v2_mode):
    c = _nums(_line(_run(_PF), 'P '))
    ref = np.polyfit([-1, 0, 1, 2, -2], [1, 0, 1, 4.1, 4], 2)[::-1]              # 独立な最小二乗 (numpy)
    assert np.allclose(c, ref, rtol=0, atol=1e-13)
    # DP0 を対称に振ったチューンのスキャン (クロマティシティ) の形: x の平均が厳密に 0
    xs = [-0.002, -0.001, 0.0, 0.001, 0.002]
    ys = [0.31 + 2.5 * x - 40.0 * x * x + 3e3 * x ** 3 for x in xs]
    c = _nums(_line(_run('Print["P ",PolynomialFit[{%s},3][[1]]];' % ','.join('{%r,%r}' % t for t in zip(xs, ys))), 'P '))
    assert np.allclose(c, [0.31, 2.5, -40.0, 3e3], rtol=1e-9, atol=1e-12)


def test_b18_polynomialfit_k64_and_nonzero_mean(k64_mode):
    c = _nums(_line(_run(_PF), 'P '))
    assert np.allclose(c, [1.0257142857142862, 2.0485714285714285, 1.0142857142857142], rtol=1e-13)     # gs の値 (二項係数の和で壊れる)
    # 平均が 0 でないデータは K64 も正しく、モードで変わらない
    body = 'Print["P ",PolynomialFit[{{1,1},{2,0},{3,1},{4,4.1},{5,4}},2][[1]]];'
    c = _nums(_line(_run(body), 'P '))
    assert np.allclose(c, np.polyfit([1, 2, 3, 4, 5], [1, 0, 1, 4.1, 4], 2)[::-1], atol=1e-12)
    from pysad import compat
    compat.set_mode('v2')
    assert _nums(_line(_run(body), 'P ')) == c


# ---------------------------------------------------------------- B71 For
def test_b71_for_body_v2(v2_mode):
    r = _run('t=0;For[i=0,i<5,i++,t+=i^2];Print["T ",t," ",i];For[k=0,k<3,k++,If[k==1,Continue[]];Print["K ",k]];'
             'For[j=0,j<10,j++,If[j==2,Break[]]];Print["J ",j];For[m=0,m<4,m++];Print["M ",m];')
    assert _nums(_line(r, 'T ')) == [30.0, 5.0]
    assert [l for l in r if l.startswith('K ')] == ['K 0', 'K 2']
    assert _nums(_line(r, 'J ')) == [2.0]
    assert _nums(_line(r, 'M ')) == [4.0]                                       # v2: 3 引数は本体なしのループ


def test_b71_for_body_k64(k64_mode):
    r = _run('t=0;For[i=0,i<5,i++,t+=i^2];Print["T ",t," ",i];For[k=0,k<3,k++,Print["K ",k]];m=7;For[m=0,m<4,m++];Print["M ",m];')
    assert _nums(_line(r, 'T ')) == [0.0, 5.0]                                  # gs: 本体は実行されず、incr だけ進む
    assert not [l for l in r if l.startswith('K ')]
    assert _nums(_line(r, 'M ')) == [7.0]                                       # gs: 3 引数は未評価のまま


# ---------------------------------------------------------------- B70 MessageName
@pytest.mark.parametrize('mode', ['v2', 'k64'])
def test_b70_messagename(mode):
    from pysad import compat
    old = compat._mode
    compat.set_mode(mode)
    try:
        r = _run('f::usage="f does things";Print["A ",f::usage];Print["B ",MessageName[f,"usage"]];'
                 'General::mytag="my `1`";Print["C ",h::mytag];Print["D ",MessageName[h,"mytag"]];Print["E ",g::nosuch];')
    finally:
        compat.set_mode(old)
    assert _line(r, 'A ') == 'f does things'
    assert _line(r, 'B ') == 'f does things'                                    # 関数形でも f の本文 (K64 は f を無視して Undefined)
    assert _line(r, 'C ') == 'my `1`' and _line(r, 'D ') == 'my `1`'            # 無ければ General::tag (gs の h::mytag と同じ)
    assert _line(r, 'E ').startswith('Undefined MessageName')


# ---------------------------------------------------------------- B61 複素 FindRoot
def test_b61_findroot_complex_v2(v2_mode):
    r = _run('pc[t_,v_]:=Print[t,Re[v]," ",Im[v]];pc["A ",z/.FindRoot[z^2+1==0,{z,1.+1.*I}]];pc["B ",z/.FindRoot[Exp[z]==2,{z,1.+1.*I}]];'
             's=FindRoot[{z*w==I,z+w==2},{z,1.+1.*I},{w,.5}];Print["C ",Re[z/.s]," ",Im[z/.s]," ",Re[w/.s]," ",Im[w/.s]," ",Residual/.s];'
             'pc["D ",z/.FindRoot[z^3==1,{z,-1.+1.*I}]];Print["E ",x/.FindRoot[x^2-2==0,{x,1.}]];')
    a = _nums(_line(r, 'A '))
    assert abs(complex(a[0], a[1]) - 1j) < 1e-12                                # z² = -1 の根 i
    b = _nums(_line(r, 'B '))
    assert abs(b[0] - math.log(2.0)) < 1e-12 and (len(b) == 1 or abs(b[1]) < 1e-12)
    c = _nums(_line(r, 'C '))
    z, w = complex(c[0], c[1]), complex(c[2], c[3])
    assert abs(z * w - 1j) < 1e-12 and abs(z + w - 2.0) < 1e-12 and c[4] < 1e-24      # 連立: zw = i, z + w = 2
    d = _nums(_line(r, 'D '))
    assert abs(complex(d[0], d[1]) - complex(-0.5, math.sqrt(3.0) / 2.0)) < 1e-10      # 1 の 3 乗根
    assert abs(_nums(_line(r, 'E '))[0] - math.sqrt(2.0)) < 1e-10                # 実数の経路は変わらない (K64 の打ち切り条件)


def test_b61_findroot_complex_k64_is_unchanged(k64_mode):
    r = _run('Print["A ",FindRoot[z^2+1==0,{z,1.+1.*I}]];Print["E ",x/.FindRoot[x^2-2==0,{x,1.}]];')
    assert any('TypeError' in l for l in r)                                     # 従来どおり (gs は壊れた値)
    assert abs(_nums(_line(r, 'E '))[0] - math.sqrt(2.0)) < 1e-10


# ---------------------------------------------------------------- B48 FixedPointList、B52 Fit
@pytest.mark.parametrize('mode', ['v2', 'k64'])
def test_b48_fixedpointlist(mode):
    from pysad import compat
    old = compat._mode
    compat.set_mode(mode)
    try:
        r = _run('Print["A ",FixedPointList[(#+2/#)/2&,1.]];Print["B ",FixedPointList[Cos,1.,5]];')
    finally:
        compat.set_mode(old)
    a = _nums(_line(r, 'A '))
    x = 1.0
    ref = []
    for _ in range(len(a)):
        x = (x + 2.0 / x) / 2.0
        ref.append(x)
    assert np.allclose(a, ref, rtol=1e-15) and abs(a[-1] - math.sqrt(2.0)) < 1e-15
    b = _nums(_line(r, 'B '))
    assert len(b) == 5 and abs(b[0] - math.cos(1.0)) < 1e-15


@pytest.mark.parametrize('mode', ['v2', 'k64'])
def test_b52_fit_vector_independent_variable(mode):
    """Fit[data, fun, {x1,x2}, ...] は最小二乗解になる (K64 は narg を出して最小にならない値を返す)。根拠は numpy の lstsq"""
    from pysad import compat
    rng = np.random.default_rng(7)
    x1 = rng.uniform(-1, 1, 25)
    x2 = rng.uniform(-1, 1, 25)
    y = 1.5 * x1 - 0.7 * x2 + 0.3 + rng.normal(0, 0.05, 25)
    data = ','.join('{{%r,%r},%r}' % (float(a), float(b), float(c)) for a, b, c in zip(x1, x2, y))
    old = compat._mode
    compat.set_mode(mode)
    try:
        r = _run('s=Fit[{%s},a*u+b*v+c,{u,v},{a,1},{b,1},{c,0}];Print["F ",a/.s," ",b/.s," ",c/.s," ",ChiSquare/.s];' % data)
    finally:
        compat.set_mode(old)
    got = _nums(_line(r, 'F '))
    A = np.stack([x1, x2, np.ones(25)], axis=1)
    sol, res, *_ = np.linalg.lstsq(A, y, rcond=None)
    assert np.allclose(got[:3], sol, rtol=1e-6, atol=1e-8)
    assert abs(got[3] - float(res[0])) < 1e-8 * max(1.0, float(res[0]))


# ---------------------------------------------------------------- 文脈 (BeginPackage / Begin)
_PKG = '''BeginPackage[Stat`];
Begin[Stat`];
Counter = 0;
sq[x_] := x^2;
Mean2[l_List] := (Counter++; Plus @@ l/Length[l]);
Var2[l_List] := Plus @@ (sq /@ (l - Mean2[l]))/Length[l];
Calls[] := Counter;
`Exported[l_List] := Sqrt[Var2[l]];
End[];
EndPackage[];
'''


@pytest.mark.parametrize('mode', ['v2', 'k64'])
def test_context_package_file(mode):
    """.n のパッケージを Get で読む (example/08_scripting/fn_Get_Needs_BeginPackage と同じ形。gs の出力と同じ値)"""
    from pysad import compat
    d = tempfile.mkdtemp()
    with open(os.path.join(d, 'ex_stat_pkg.n'), 'w') as f:
        f.write(_PKG)
    old = compat._mode
    compat.set_mode(mode)
    try:
        r = _run('Counter=77;sq[x_]:=-x;\nGet["ex_stat_pkg.n"];\n'
                 'Print["A ",Stat`Mean2[{1,2,3,4}]," ",Stat`Var2[{1,2,3,4}]," ",Stat`Calls[]," ",Stat`Counter];\n'
                 'Print["B ",Mean2[{1,2}]];\nPrint["C ",Counter," ",sq[3]];\nPrint["D ",Context[Stat`sq]," ",Context[Counter]," ",Context[Sin]," ",Context[]];\n'
                 'Print["E ",Exported[{1,2,3,4}]];\nPrint["F ",$Context," ",$ContextPath];\nPrint["G ",Names["Stat`*"]];', cwd=d)
    finally:
        compat.set_mode(old)
        os.remove(os.path.join(d, 'ex_stat_pkg.n'))
        os.rmdir(d)
    assert _nums(_line(r, 'A ')) == [2.5, 1.25, 2.0, 2.0]                       # gs: 2.5 1.25 2 2
    assert _line(r, 'B ') == 'Mean2[{1,2}]'                                     # 文脈なしの名前は見えない (gs)
    assert _nums(_line(r, 'C ')) == [77.0, -3.0]                                # 大域の同名の変数・関数はパッケージに壊されない
    assert _line(r, 'D ') == 'Stat` ` System` System`'
    assert abs(_nums(_line(r, 'E '))[0] - math.sqrt(1.25)) < 1e-15              # `名前 で大域に公開した関数
    assert _line(r, 'F ') == '` {`,System`}'
    g = _line(r, 'G ')
    if mode == 'k64':
        assert g == '{"Stat`"}'                                                 # gs
    else:
        assert '"Stat`Mean2"' in g and '"Stat`Counter"' in g                    # v2: 文脈の記号の一覧


def test_context_resolution_rules(k64_mode):
    """gs で確かめた名前解決 (現在の文脈 → $ContextPath の順、無ければ現在の文脈に作る。` は大域)"""
    r = _run('''ug = 5; uf[x_] := x + 100; vg[x_] := x + 200; Counter = 77;
r1 = BeginPackage[Stat`];
Print["1 ", $ContextPath, " ", r1];
Begin[Stat`];
Print["2 ", $Context];
Counter = 0;
sq[x_Real] := x^2;
Mean2[l_List] := (Counter++; Plus @@ l/Length[l]);
Print["3 ", ug, " ", uf[1], " ", `ug, " ", Head[MOMENTUM], " ", sq[3.], " ", Context[sq], " ", Context[ug]];
tt = {Head[1.5], Head["a"], Head[{1}], True, Pi};
`gl = 9;
Stat2`zz = 4;
End[];
Print["4 ", sq[3.], " ", Mean2[{1, 3}], " ", Counter];
EndPackage[];
Print["5 ", $Context, " ", $ContextPath];
Print["6 ", sq[3.], " | ", Stat`sq[3.], " | ", Counter, " | ", Stat`Counter, " | ", gl, " | ", Stat2`zz, " | ", Stat`gl];
Print["7 ", Stat`tt];
BeginPackage[Stat`, `]; Begin[Stat`];
Print["9 ", sq[2.], " ", uf[1], " ", vg[1], " ", $ContextPath];
newone = 3;
End[]; EndPackage[];
Print["10 ", newone, " ", Stat`newone];
Begin[Stat`]; qq = 8; End[];
Print["12 ", qq, " ", Stat`qq];
BeginPackage[Stat`]; pp = 6; EndPackage[];
Print["14 ", pp, " ", Stat`pp];''')
    assert _line(r, '1 ') == '{`Stat`,`System`} {`Stat`,`System`}'               # gs: ` が $ContextPath に無いと、文脈の記号は `Stat` と表示される
    assert _line(r, '2 ') == '`Stat`'
    assert _line(r, '3 ') == 'ug uf[1] 5 Real 9 `Stat` `Stat`'                    # 大域の ug・uf は ` が $ContextPath に無いと見えない。`ug は大域
    assert _nums(_line(r, '4 ')) == [9.0, 2.0, 77.0]                            # End の後: 現在の文脈 (大域) の Counter が先
    assert _line(r, '5 ') == '` {`,System`}'
    assert _line(r, '6 ') == 'sq[3] | 9 | 77 | 1 | 9 | 4 | Stat`gl'
    assert _line(r, '7 ').startswith('{Real,String,List,1,3.14159')
    # ` を $ContextPath に入れると大域の vg が見える。uf は先に Stat` の中で読んだので Stat`uf (未定義) が先に見つかる (gs も同じ)
    assert _line(r, '9 ') == '4 uf[1] 201 {Stat`,`,System`}'
    assert _line(r, '10 ') == 'newone 3'
    assert _line(r, '12 ') == 'qq 8'
    assert _line(r, '14 ') == '6 Stat`pp'


def test_context_package_using_ffs_and_builtins(v2_mode):
    """文脈の中から、組み込み関数・SAD のパッケージ (自動読み込み)・大域の変数 (` を付ける) が使える"""
    r = _run('''gscale = 2;
BeginPackage[My`];
Begin[My`];
helper[x_] := Module[{t = x}, t*`gscale];
fit[d_] := PolynomialFit[d, 1][[1]];
`MyTotal[l_List] := helper[Plus @@ l];
`MyFit = fit;
End[];
EndPackage[];
Print["A ", MyTotal[{1, 2, 3}]];
Print["B ", MyFit[{{1, 3}, {2, 5}, {3, 7}}]];
Print["C ", helper[1], " ", My`helper[1]];''')
    assert _nums(_line(r, 'A ')) == [12.0]
    assert np.allclose(_nums(_line(r, 'B ')), [1.0, 2.0], atol=1e-12)
    assert _line(r, 'C ').startswith('helper[1] ') and _nums(_line(r, 'C '))[-1] == 2.0


# ---------------------------------------------------------------- Shared / Fork / System / Sleep
@pytest.mark.parametrize('mode', ['v2', 'k64'])
def test_shared_memory_in_process(mode):
    from pysad import compat
    old = compat._mode
    compat.set_mode(mode)
    try:
        r = _run('s=OpenShared[128];Print["U ",Shared[s]];Shared[s]={1,2,3};Print["A ",Shared[s]," ",SharedSize[{1,2,3}]," ",SharedSize["abcdefghi"]," ",SharedSize[1.5]];'
                 'Shared[s]=5;Shared[s]++;Shared[s]+=10;Print["B ",Shared[s]];Shared[s]={x->1,"str",{2,Null}};Print["C ",Shared[s]];'
                 'Shared[s]=Table[i,{i,100}];Print["D ",Shared[s]];Close[s];Shared[s]=1;')
    finally:
        compat.set_mode(old)
    assert _line(r, 'U ') == 'Undefined'
    assert _nums(_line(r, 'A ')) == [1.0, 2.0, 3.0, 32.0, 24.0, 0.0]            # gs: SharedSize[{1,2,3}] = 32 ((3+3-2)*8)
    assert _nums(_line(r, 'B ')) == [16.0]
    assert _line(r, 'C ') == '{x->1,"str",{2,}}'
    assert any('Shared::toolarge' in l for l in r)                              # OpenShared[128] に 100 個の実数は入らない (K64 と同じ判定)
    assert _line(r, 'D ') == '{x->1,"str",{2,}}'                                # 入らなかった値は書かれない
    assert any('Shared::notopen' in l for l in r)


@pytest.mark.skipif(hasattr(os, 'fork'), reason='Windows (fork が無い環境) の動作')
def test_fork_on_windows():
    from pysad import compat
    old = compat._mode
    try:
        compat.set_mode('v2')
        r = _run('pid=Fork[];Print["A ",Head[pid]];')
        assert any('Fork[] is not available' in l and 'Parallelize' in l for l in r)    # 分かるエラー (子は作れない)
        compat.set_mode('k64')
        r = _run('pid=Fork[];Print["A ",Head[pid]," ",pid>0];')
        assert _line(r, 'A ') == 'Real 1'                                       # 従来どおり (親の側)
    finally:
        compat.set_mode(old)


@pytest.mark.skipif(not hasattr(os, 'fork'), reason='POSIX (Linux・WSL・macOS) だけ')
def test_fork_and_shared_on_posix(v2_mode):
    """子が共有メモリに書いた値を親が読む (gs の example/09_io_system/fn_Fork_Shared と同じ結果)"""
    r = _run('n=3;sh=Table[OpenShared[256],{n}];pids=Table[0,{n}];'
             'Do[If[(pids[[k]]=Fork[])==0,Shared[sh[[k]]]={k,k^2,"child"};Exit[k]],{k,n}];'
             'Do[w=Wait4[pids[[k]]];Print["C ",k," ",w[[2]]," ",Exited/.Drop[w,2]," ",Shared[sh[[k]]]],{k,n}];')
    for k in (1, 2, 3):
        assert 'C %d %d 1 {%d,%d,"child"}' % (k, k, k, k * k) in r


@pytest.mark.parametrize('mode', ['v2', 'k64'])
def test_parallelize_serial(mode):
    from pysad import compat
    old = compat._mode
    compat.set_mode(mode)
    try:
        r = _run('f[k_]:=k^2+1;Print["P ",Parallelize[f,4,64,2]];')
    finally:
        compat.set_mode(old)
    assert _nums(_line(r, 'P ')) == [2.0, 5.0, 10.0, 17.0]


def test_system_is_disabled_by_default(monkeypatch):
    from pysad import compat, miscfn
    monkeypatch.delenv('PYSAD_ALLOW_SYSTEM', raising=False)
    miscfn.allow_system(None)
    d = tempfile.mkdtemp()
    mark = os.path.join(d, 'made_by_system.txt').replace(os.sep, '/')
    old = compat._mode
    try:
        compat.set_mode('v2')
        r = _run('r=System["echo x > %s"];Print["A ",r];Print["B ",System["[ -d %s ]"]," ",System["[ -e %s ]"]," ",System["test -f %s"]];' % (mark, d.replace(os.sep, '/'), mark, mark))
        assert not os.path.exists(mark)                                          # 実行されていない
        assert _line(r, 'A ') == '$Failed'
        assert any('PYSAD_ALLOW_SYSTEM' in l for l in r)                         # 許可の仕方を知らせる
        assert _nums(_line(r, 'B ')) == [0.0, 1.0, 1.0]                          # 存在確認だけはシェルなしで判定 (0: 真)
        compat.set_mode('k64')
        r = _run('r=System["echo x > %s"];Print["A ",r];' % mark)
        assert not os.path.exists(mark) and _nums(_line(r, 'A ')) == [0.0]       # 従来どおり
    finally:
        compat.set_mode(old)
        os.rmdir(d)


def test_system_runs_only_with_permission(monkeypatch):
    from pysad import compat, miscfn
    monkeypatch.delenv('PYSAD_ALLOW_SYSTEM', raising=False)
    monkeypatch.delenv('PYSAD_SYSTEM_SHELL', raising=False)
    d = tempfile.mkdtemp()
    old = compat._mode
    try:
        miscfn.allow_system(True)
        compat.set_mode('v2')
        r = _run('Print["A ",System["echo from-shell > made.txt"]];Print["B ",System["exit 3"]];Print["C ",System["echo hello-out"]];', cwd=d)
        assert os.path.isfile(os.path.join(d, 'made.txt'))
        assert open(os.path.join(d, 'made.txt')).read().strip() == 'from-shell'
        assert _nums(_line(r, 'A ')) == [0.0] and _nums(_line(r, 'B ')) == [3.0]      # v2: 終了コード (SADHelp の "system return code")
        assert any(l.strip() == 'hello-out' for l in r)                          # コマンドの標準出力は出力先に出る
        compat.set_mode('k64')
        r = _run('Print["B ",System["exit 3"]];', cwd=d)
        assert _nums(_line(r, 'B ')) == [1.0]                                    # K64: 起動できれば終了コードに関係なく 1
        # SADScript からは許可を変えられない
        miscfn.allow_system(False)
        compat.set_mode('v2')
        r = _run('SetEnv["PYSAD_ALLOW_SYSTEM","1"];Environment["PYSAD_ALLOW_SYSTEM","1"];Print["D ",System["echo no"]];', cwd=d)
        assert _line(r, 'D ') == '$Failed'
        miscfn.allow_system(None)                                                # 環境変数に従う設定でも、スクリプトの SetEnv は効かない (起動時の値だけ)
        r = _run('SetEnv["PYSAD_ALLOW_SYSTEM","1"];Print["D ",System["echo no"]];', cwd=d)
        assert _line(r, 'D ') == '$Failed'
    finally:
        miscfn.allow_system(None)
        compat.set_mode(old)
        os.environ.pop('PYSAD_ALLOW_SYSTEM', None)
        for f in os.listdir(d):
            os.remove(os.path.join(d, f))
        os.rmdir(d)


def test_sleep_waits():
    import time
    t0 = time.time()
    r = _run('Print["S ",Sleep[0.2]];')
    assert time.time() - t0 >= 0.19 and _nums(_line(r, 'S ')) == [1.0]
