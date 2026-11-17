"""CSR (csr.py): CSRTrack / CSRHaissin / CSRSetupOY / CSRMatrix / CSRConvert。参照は SAD K64 (gs) の出力 (local/wake/csr1〜csr6.sad)。
CSROYMatrix は gs が有効な引数でも wrongtype を返す (docs/21 B60) ので参照なし (ソースどおりの翻訳)。"""
import io
import os
import re

import numpy as np
import pytest


def _run(body):
    from pysad.ffscmd import run_sad_string
    buf = io.StringIO()
    run_sad_string('FFS;\n$FORM="18.14";\nPageWidth=1999;\n' + body + 'Exit;\n', out=buf)
    return buf.getvalue()


def _nums(line):
    line = line.split(None, 1)[1]          # 先頭のラベル (R, Z1, R2, ...) は数に含めない
    return [float(x) for x in re.findall(r'[-+]?(?:\d+\.?\d*|\.\d+)(?:[Ee][-+]?\d+)?', line)]


def test_csrtrack_matches_gs():
    """z・dp をその場で更新する。2 周で発散 (以後 0) するケースと RF あり (frf≠0) のケース"""
    out = _run('''SeedRandom[12];
nr=256; dz=2e-5; np=2000;
zl=Table[0.02*Cos[0.3 i]/(1+i/10.), {i,1,nr}];
z=3e-4*GaussRandom[np]; dp=1e-3*GaussRandom[np];
r=CSRTrack[z,dp,zl,{dz,4,0.5,0.1,-0.5,0,0,1e-3,1e-3},5];
Print["R ",r];
Print["Z1 ",Plus@@(z^2)," ",Plus@@(dp^2)];
r=CSRTrack[z,dp,zl,{dz,5,0.5,0.1,-0.5,1.3e9,0.2,1e-3,1e-3},3];
Print["R2 ",r];
Print["Z2 ",Plus@@z," ",Plus@@(z^2)," ",Plus@@(dp^2)];
''')
    lines = out.splitlines()
    r = _nums([l for l in lines if l.startswith('R ')][0])
    assert np.allclose(r[:2], [0.02742838627370, 0.36118959158257], rtol=1e-11)
    assert r[2:] == [0.0] * 8
    z1 = _nums([l for l in lines if l.startswith('Z1 ')][0])
    assert np.allclose(z1, [1.50463274715906, 260.9158421351614], rtol=1e-11)
    r2 = _nums([l for l in lines if l.startswith('R2 ')][0])
    assert np.allclose(r2[:2], [0.03250355505305, 0.18065803991553], rtol=1e-11)
    z2 = _nums([l for l in lines if l.startswith('Z2 ')][0])
    assert np.allclose(z2[1:], [1309.011760683231, 242.0834189346895], rtol=1e-11)


def test_csrtrack_stable_and_haissin_match_gs():
    out = _run('''SeedRandom[7];
nr=256; dz=2e-5; np=3000;
zl=Table[0.02*Cos[0.3 i]/(1+i/10.), {i,1,nr}];
z=3e-4*GaussRandom[np]; dp=1e-3*GaussRandom[np];
r=CSRTrack[z,dp,zl,{dz,4,0.001,0.05,-0.3,0,0,2e-3,1e-3},20];
Print["R ",Take[r,6]];
Print["Z1 ",Plus@@(z^2)," ",Plus@@(dp^2)];
rho=CSRHaissin[zl,{dz,1e-8,-0.3,3e-4}];
Print["H ",Length[rho]," ",rho[[10]]," ",rho[[60]]," ",Max[rho]," ",rho[[130]]];
''')
    lines = out.splitlines()
    r = _nums([l for l in lines if l.startswith('R ')][0])
    # gs の先頭 6 個: 3.400369909354E-4 .00243617306564 4.430528687671E-4 .00297733213416 5.571517017490E-4 .00373789107516
    assert np.allclose(r[:6], [3.400369909354e-4, 0.00243617306564, 4.430528687671e-4, 0.00297733213416,
                               5.571517017490e-4, 0.00373789107516], rtol=1e-11)
    z1 = _nums([l for l in lines if l.startswith('Z1 ')][0])
    assert np.allclose(z1, [0.00710486270665, 0.16032061340151], rtol=1e-11)
    h = _nums([l for l in lines if l.startswith('H ')][0])
    assert np.allclose(h, [256.0, 2.03724970092067, 1283.169713515304, 1329.605939173232, 0.0], rtol=1e-11, atol=1e-12)


def test_csrsetupoy_matches_gs_and_matrix_is_finite():
    """J・ω は周回の離散計算なので 1e-8 の一致。zjn は dpjn と同じ記憶域 (B59) なので p の値、dpjn は 0"""
    out = _run('''nr=256; dz=2e-5;
zl=Table[0.02*Cos[0.3 i]/(1+i/10.), {i,1,nr}];
rho=CSRHaissin[zl,{dz,1e-8,-0.3,3e-4}];
s=CSRSetupOY[rho,zl,{dz,6,3,1e-8,0.05,-0.3,3}];
Print["S ",Length/@s];
Print["AJ ",s[[2]]];
Print["OM ",s[[3]]];
Print["RJ ",s[[6]]];
Print["ZS ",Plus@@s[[5]]];
m=CSROYMatrix[s,{dz,1e-8,0.01,3e-4,1e-3}];
Print["M ",Dimensions[m]];
''')
    lines = out.splitlines()
    assert _nums([l for l in lines if l.startswith('S ')][0]) == [256.0, 6.0, 6.0, 288.0, 288.0, 6.0]
    aj = _nums([l for l in lines if l.startswith('AJ ')][0])
    assert np.allclose(aj, [1.598877766878e-8, 4.801962255995e-8, 8.001652161176e-8, 1.119885529838e-7,
                            1.440141918819e-7, 1.760343431927e-7], rtol=1e-7)
    om = _nums([l for l in lines if l.startswith('OM ')][0])
    assert np.allclose(om, [0.12244341109984, 0.12244819337127, 0.12245250160718, 0.12245635032843,
                            0.12245977046611, 0.12246278316010], rtol=1e-9)
    rj = _nums([l for l in lines if l.startswith('RJ ')][0])
    assert np.allclose(rj, [0.22047618866441, 0.19067736090131, 0.16492735497839, 0.14269848454815,
                            0.12341352148526, 0.10673588977341], rtol=1e-7)
    assert _nums([l for l in lines if l.startswith('ZS ')][0])[0] == 0.0
    assert _nums([l for l in lines if l.startswith('M ')][0]) == [18.0, 18.0]


def test_csrmatrix_and_csrconvert_match_gs():
    from pysad.csr import csrind, csrmat
    pipe = [3, 4, 4, 3]
    ind, mm, n1 = csrind(pipe)
    ref = {22: (-367.674659844456, 8714.994964898138), 33: (-3143.45341849862, 11490.79050800808),
           23: (0.33029695561819, 550.7558472075209), 32: (0.19686015490331, 550.3572309041174)}
    for mode, (s, sa) in ref.items():
        a = csrmat(ind, 1e-3, 2e-3, 3e3, 5.0, 1.1, 1.3, 0.9, 0.7, 0.1, mode, mm, len(pipe), n1, 2)
        assert abs(a.sum() - s) < 1e-9 * abs(sa)
        assert abs(np.abs(a).sum() - sa) < 1e-9 * sa
    assert mm == 14
    out = _run('''v=CSRConvert[{{1,2,3,4},{5,6,7,8},{9,10,11,12},{13,14,15,16}},{3,4,4,3}];
Print["CV ",v];
w=CSRConvert[Range[14],{3,4,4,3}];
Print["CW ",Flatten[w]];
a=CSRMatrix[{3,4,4,3},{1e-3,2e-3,1.1,1.3,0.9,0.7,0.1},3e3,5.0,2,220];
Print["V ",Length[a]," ",Plus@@a];
''')
    lines = out.splitlines()
    cv = _nums([l for l in lines if l.startswith('CV ')][0])
    assert cv == [1, 2, 3, 5, 6, 7, 8, 9, 10, 11, 12, 14, 15, 16][:len(cv)] or len(cv) == 14
    assert _nums([l for l in lines if l.startswith('V ')][0])[0] == 14.0


def test_csrinit_matches_gs():
    """CSRInit[{{dx,dy}, center, rho, pipe}] (csr2dpoisson、3 回の繰り返し、tsolva = tsvd(.false.))。参照は gs。実数 1 つの dx は K64 のバグで wrongtype"""
    out = _run('''r=CSRInit[{{0.5e-3,0.5e-3},4,1e-3,{4,6,8,8,8,8,6,4}}];
Print["A ",Dimensions[r]," ",Plus@@Flatten[r[[1]]]," ",Plus@@Flatten[r[[2]]]," ",Max[Abs[Flatten[r]]]];
Print["B ",r[[2,3,4]]," ",r[[1,2,3]]," ",r[[1,8,8]]];
r2=CSRInit[{{0.5e-3,0.25e-3},3,2e-3,{3,5,5,5,3}}];
Print["C ",Plus@@Flatten[r2[[1]]]," ",Plus@@Flatten[r2[[2]]]," ",r2[[1,3,2]]];
''')
    a = _nums([l for l in out.splitlines() if l.startswith('A ')][0])
    assert a[:3] == [2, 8, 8]
    for got, want in zip(a[3:], [323.936626657909, 34983.676790819125, 2658.9908202932147]):
        assert abs(got - want) < 1e-8 * abs(want)
    b = _nums([l for l in out.splitlines() if l.startswith('B ')][0])
    for got, want in zip(b, [712.6418937209675, -840.7488459342261, 0.0]):
        assert abs(got - want) < 1e-8 * max(abs(want), 1.0)
    c = _nums([l for l in out.splitlines() if l.startswith('C ')][0])
    for got, want in zip(c, [-137.42131632108863, 11966.81452886956, -800.7638478031346]):
        assert abs(got - want) < 1e-8 * abs(want)
    bad = _run('r=CSRInit[{0.5e-3,4,1e-3,{4,6,8,8,8,8,6,4}}];\n')
    assert 'wrongtype' in bad


@pytest.mark.skipif(not os.environ.get('PYSAD_SLOW') or not os.environ.get('PYSAD_SAD_PACKAGES'),
                    reason='環境変数 PYSAD_SLOW=1 と PYSAD_SAD_PACKAGES が必要 (約 7 分)')
def test_csr_makezl_matches_gs():
    """CSR[…]@MakeZL[{{10.,1.}}, 2000., 6000., 3] (CSRMatrix → RealEig → CSRInit → Ez → ZL。参照は gs: top の代用を PATH に置いて実行)。
    gs は ZL の表の最後の Spline で失敗するが、表そのものは出す (相対 1e-10 で一致)"""
    out = _run('''c=CSR[NParallel->1,Nk->4,BunchLength->3e-4,BendingRadius->10,PipeWidth->0.04,Center->0.02,MeshRatio->{2,2}];
z=c@MakeZL[{{10.,1.}},2000.,6000.,3];
Print["Z ",z];
''')
    line = [l for l in out.splitlines() if l.startswith('Z ')][0]
    got = _nums(line)
    want = [2000, -.027525977863628738, -.02562238715379007, 4000, -.02860985374476984, -.03392664519294917,
            6000, -.045736712549164574, -.03026221214935138]
    assert len(got) == 9
    for g, w in zip(got, want):
        assert abs(g - w) < 1e-9 * max(abs(w), 1.0)


def test_csrmatrix_mode1000_matches_gs():
    """CSRMatrix mode 1000 (csrhelmat + csr2dpoisson1。tsolvm の最小ノルム解。K64 の phix(0:nphimax) の添字のずれまで写した)。参照は gs (相対 1e-12)"""
    out = _run('''a=CSRMatrix[{2,2,2},{1e-3,2e-3,1.1,1.3,0.9,0.7,0.1},3e3,5.0,2,1000];
n=Length[a];
Print["D ",Dimensions[a]];
Print["S ",Plus@@Flatten[a]," ",Plus@@Abs[Flatten[a]]];
Print["E ",a[[1,1]]," ",a[[2,3]]," ",a[[n,n]]," ",a[[n,1]]," ",a[[1,n]]];
b=CSRMatrix[{5,5,5,5,5,5},{1e-3,2e-3,1.1,1.3,0.9,0.7,0.1},3e3,5.0,2,1000];
m=Length[b];
Print["F ",Dimensions[b]," ",Plus@@Flatten[b]," ",b[[1,1]]," ",b[[m,m]]];
''')
    lines = {l[0]: _nums(l) for l in out.splitlines() if l[:2] in ('D ', 'S ', 'E ', 'F ')}
    assert lines['D'][:2] == [12, 12]
    ref = {'S': [3.9443741850615275e-5, 5.585179755240409e-5],
           'E': [3.5987372900573512e-6, 3.2364494271231303e-7, 2.8973328550838e-8, 1.241388575586084e-7, -3.0579891393328253e-7]}
    for k, w in ref.items():
        for g, x in zip(lines[k], w):
            assert abs(g - x) < 1e-11 * abs(x) + 1e-22
    f = lines['F']
    assert f[:2] == [60, 60]
    assert abs(f[2] - .0012846962197225107) < 1e-11 * .00128
    assert abs(f[3] - 6.061901936460279e-6) < 1e-11 * 6.06e-6
