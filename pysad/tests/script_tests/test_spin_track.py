"""POL 時の TrackParticles (粒子ごとのスピン追跡: track6 の tradk の sprot、bsi)。参照は SAD K64 (gs) の出力 (local/wake/pt*.sad)。
RAD ON (NOFLUC / FLUC) のときだけスピンが動く (K64 は tradk の中で歳差を計算する)。入力は 9〜11 列 {x,px,y,py,z,dp, 角度, sy/|s|, flag}、
戻りの 7・8 列は atan(sz, sx) と sy/|s|。"""
import io
import re

CASES = {}
CASES['pt1'] = ('MOMENTUM = 3GEV;\nDRIFT D1 = (L=1) D2 = (L=0.5);\nQUAD QF = (L=0.3 K1=0.7) QD = (L=0.3 K1=-0.6);\nBEND B = (L=1.2 ANGLE=0.2617993877991494 E1=0.1308996938995747 E2=0.1308996938995747);\nCAVI RFC = (VOLT=1.5e6 HARM=64);\nMARK ST = ();\nLINE CELL1 = (QF D2 B D2 QD D1);\nLINE RING = (ST 12*CELL1 RFC);\nFFS;\n$FORM="18.14";\nPageWidth = 1999;\nUSE RING;\nCELL;\nCALC;\nRFSW; RAD; NOFLUC;\nPOL;\nem = Emittance[];\nPrint["EM ",{SpinTune,NominalSpinTune}/.em];\nx0={1e-4,0};px0={0,0};y0={0,1e-4};py0={0,0};z0={0,0};dp0={0,1e-3};\nr=TrackParticles[{1, {x0,px0,y0,py0,z0,dp0,{0,0},{0.1,0.2},{1,1},{1,1}}}, 20, 1, 3];\nPrint["R ",Length[r]," ",Dimensions[r[[2]]]];\nPrint["S ",r[[2]]];\nExit;\n',
    [0.0001302713362823, 0.000375905169844, -1.30817025201e-05, 0.0002065853632128, 3.03905945049e-20, 0.0001651632081391, 2.3097998852e-20, -0.000104740471749, 0.0003781806765284, -0.00288869797463, -6.2614386424e-05, 0.0009188061726734, -2.14774794795588, -2.10190448442089, 0.10000000225364, 0.20083028975489, 1.0, 1.0])
CASES['pt3'] = ('MOMENTUM = 3GEV;\nDRIFT D1 = (L=1) D2 = (L=0.5);\nQUAD QF = (L=0.3 K1=0.7) QD = (L=0.3 K1=-0.6);\nBEND B = (L=1.2 ANGLE=0.2617993877991494 E1=0.1308996938995747 E2=0.1308996938995747);\nCAVI RFC = (VOLT=1.5e6 HARM=64);\nMARK ST = ();\nLINE CELL1 = (QF D2 B D2 QD D1);\nLINE RING = (ST 12*CELL1 RFC);\nFFS;\n$FORM="18.14";\nPageWidth = 1999;\nUSE RING;\nCELL;\nCALC;\nRFSW; RAD; FLUC;\nPOL;\nem = Emittance[];\nPrint["EM ",{SpinTune,NominalSpinTune}/.em];\nx0={1e-4,0};px0={0,0};y0={0,1e-4};py0={0,0};z0={0,0};dp0={0,1e-3};\nr=TrackParticles[{1, {x0,px0,y0,py0,z0,dp0,{0,0},{0.1,0.2},{1,1},{1,1}}}, 20, 1, 3];\nPrint["R ",Length[r]," ",Dimensions[r[[2]]]];\nPrint["S ",r[[2]]];\nExit;\n',
    [0.0001693990676281, 0.000438344758072, -2.36624905398e-06, 0.000231715009685, 3.80603156669e-20, 0.0001652467535584, 2.28512696119e-20, -0.000104687713312, 0.0004419843400547, -0.00302480379272, -4.68611831562e-05, 0.00100987820921, -2.14877240071237, -2.10016117744606, 0.10000000219484, 0.20083082556555, 1.0, 1.0])
CASES['pt4'] = ('MOMENTUM = 3GEV;\nDRIFT D1 = (L=1) D2 = (L=0.5) D3=(L=0.2);\nQUAD QF = (L=0.3 K1=0.7) QD = (L=0.3 K1=-0.6);\nQUAD QFE = (L=0.3 K1=0.7 DY=1e-4) QDE = (L=0.3 K1=-0.6 DY=-1.5e-4 ROTATE=0.001);\nQUAD QT = (L=0.2 K1=0.02 ROTATE=0.5);\nSEXT SF = (L=0.1 K2=2.0 DY=1e-4);\nBEND B = (L=1.2 ANGLE=0.2617993877991494 E1=0.1308996938995747 E2=0.1308996938995747);\nBEND BR = (L=1.2 ANGLE=0.2617993877991494 E1=0.1308996938995747 E2=0.1308996938995747 ROTATE=0.0015);\nCAVI RFC = (VOLT=1.5e6 HARM=64);\nMARK ST = ();\nLINE CELL1 = (QF D2 B D2 QD D1);\nLINE CELL2 = (QFE D3 SF D3 B D2 QDE D1);\nLINE CELL3 = (QF D2 BR D2 QT D1 QD D1);\nLINE RING = (ST 8*CELL1 CELL2 CELL3 CELL2 CELL3 RFC);\nFFS;\n$FORM="18.14";\nPageWidth = 1999;\nUSE RING;\nCELL;\nCALC;\nRFSW; RAD; NOFLUC;\nPOL;\nem = Emittance[];\nPrint["EM ",{SpinTune,NominalSpinTune}/.em];\nx0={1e-4,0};px0={0,0};y0={0,1e-4};py0={0,0};z0={0,0};dp0={0,1e-3};\nr=TrackParticles[{1, {x0,px0,y0,py0,z0,dp0,{0,0},{0.1,0.2},{1,1},{1,1}}}, 20, 1, 5];\nPrint["R ",Length[r]," ",Dimensions[r[[2]]]];\nPrint["S ",r[[2]]];\nExit;\n',
    [0.0009815965106435, 0.0006477291081407, 0.0002421362510943, 0.0003199206023513, 0.00203888308936, 0.00219455099708, -0.00128700883856, -0.00140580045512, 0.0007401474213401, -0.00535615833209, -5.70191897792e-05, 0.0008859909686444, 2.92091002248385, 3.0076631820702, 0.09173213264127, 0.19242358616221, 1.0, 1.0])
CASES['pt5'] = ('MOMENTUM = 3GEV;\nDRIFT D1 = (L=1) D2 = (L=0.5);\nQUAD QF = (L=0.3 K1=0.7) QD = (L=0.3 K1=-0.6);\nBEND B = (L=1.2 ANGLE=0.2617993877991494 E1=0.1308996938995747 E2=0.1308996938995747);\nCAVI RFC = (VOLT=1.5e6 HARM=64);\nMARK ST = ();\nMULT MM = (L=0.1 K1=0.05 K2=2.0 DY=1e-4);\nMULT MQ = (L=0.2 K1=0.02 K3=20 DY=-3e-5 DX=2e-5);\nLINE CELL1 = (QF D2 B D2 MM QD D1 MQ);\nLINE RING = (ST 12*CELL1 RFC);\nFFS;\n$FORM="18.14";\nPageWidth = 1999;\nUSE RING;\nCELL;\nCALC;\nRFSW; RAD; NOFLUC;\nPOL;\nem = Emittance[];\nPrint["EM ",{SpinTune,NominalSpinTune}/.em];\nx0={1e-4,0};px0={0,0};y0={0,1e-4};py0={0,0};z0={0,0};dp0={0,1e-3};\nr=TrackParticles[{1, {x0,px0,y0,py0,z0,dp0,{0,0},{0.1,0.2},{1,1},{1,1}}}, 20, 1, 5];\nPrint["R ",Length[r]," ",Dimensions[r[[2]]]];\nPrint["S ",r[[2]]];\nExit;\n',
    [-0.00033234746885, 0.000184081091859, 7.564606709837e-06, -5.39379562111e-05, -1.33621629576e-05, 0.0001797387437109, -9.38444590722e-06, 4.221302552685e-05, 0.0005429204656226, -0.00429068856584, -3.73355207485e-05, 0.0009190395784645, 1.14157814490935, 1.22999193806627, 0.09995579568639, 0.19986294825146, 1.0, 1.0])
CASES['pt6'] = ('MOMENTUM = 3GEV;\nDRIFT D1 = (L=1) D2 = (L=0.5);\nQUAD QF = (L=0.3 K1=0.7) QD = (L=0.3 K1=-0.6);\nBEND B = (L=1.2 ANGLE=0.2617993877991494 E1=0.1308996938995747 E2=0.1308996938995747);\nCAVI RFC = (VOLT=1.5e6 HARM=64);\nMARK ST = ();\nBEND BK = (L=1.2 ANGLE=0.2617993877991494 K1=0.1 E1=0.1308996938995747 E2=0.1308996938995747 DX=1e-4);\nLINE CELL1 = (QF D2 BK D2 QD D1);\nLINE RING = (ST 12*CELL1 RFC);\nFFS;\n$FORM="18.14";\nPageWidth = 1999;\nUSE RING;\nCELL;\nCALC;\nRFSW; RAD; NOFLUC;\nPOL;\nem = Emittance[];\nPrint["EM ",{SpinTune,NominalSpinTune}/.em];\nx0={1e-4,0};px0={0,0};y0={0,1e-4};py0={0,0};z0={0,0};dp0={0,1e-3};\nr=TrackParticles[{1, {x0,px0,y0,py0,z0,dp0,{0,0},{0.1,0.2},{1,1},{1,1}}}, 20, 1, 5];\nPrint["R ",Length[r]," ",Dimensions[r[[2]]]];\nPrint["S ",r[[2]]];\nExit;\n',
    [2.611896522049e-05, 0.0008178070444642, -9.58883224659e-05, 3.91530834294e-05, 2.54116997826e-20, 0.0001658711434302, 1.79436023088e-21, -4.96951349496e-05, 0.0002653802312743, -0.00511811871732, -6.13096938468e-05, 0.0008877675251837, 2.92576743068356, 3.01388208008631, 0.1000000042561, 0.19999572319429, 1.0, 1.0])
CASES['pt7'] = ('MOMENTUM = 3GEV;\nDRIFT D1 = (L=1) D2 = (L=0.5);\nQUAD QF = (L=0.3 K1=0.7) QD = (L=0.3 K1=-0.6);\nBEND B = (L=1.2 ANGLE=0.2617993877991494 E1=0.1308996938995747 E2=0.1308996938995747);\nCAVI RFC = (VOLT=1.5e6 HARM=64);\nMARK ST = ();\nSEXT SF = (L=0.1 K2=2.0 DY=1e-4);\nMULT TH = (L=0 K2=0.8 DY=5e-5);\nLINE CELL1 = (QF D2 SF B D2 TH QD D1);\nLINE RING = (ST 12*CELL1 RFC);\nFFS;\n$FORM="18.14";\nPageWidth = 1999;\nUSE RING;\nCELL;\nCALC;\nRFSW; RAD; NOFLUC;\nPOL;\nem = Emittance[];\nPrint["EM ",{SpinTune,NominalSpinTune}/.em];\nx0={1e-4,0};px0={0,0};y0={0,1e-4};py0={0,0};z0={0,0};dp0={0,1e-3};\nr=TrackParticles[{1, {x0,px0,y0,py0,z0,dp0,{0,0},{0.1,0.2},{1,1},{1,1}}}, 20, 1, 5];\nPrint["R ",Length[r]," ",Dimensions[r[[2]]]];\nPrint["S ",r[[2]]];\nExit;\n',
    [-8.82045089561e-05, 0.00101961131947, 8.508932952526e-05, -0.000450082074055, -8.5559334739e-08, 0.0001094685901498, -1.50624074082e-07, 7.564174064823e-05, 0.000778736735474, -0.00508370914882, -3.6374763284e-05, 0.0009092859598788, 1.14101050811342, 1.23269586254823, 0.09999877672396, 0.19992473240412, 1.0, 1.0])
CASES['pt7b'] = ('MOMENTUM = 3GEV;\nDRIFT D1 = (L=1) D2 = (L=0.5);\nQUAD QF = (L=0.3 K1=0.7) QD = (L=0.3 K1=-0.6);\nBEND B = (L=1.2 ANGLE=0.2617993877991494 E1=0.1308996938995747 E2=0.1308996938995747);\nCAVI RFC = (VOLT=1.5e6 HARM=64);\nMARK ST = ();\nSEXT SF = (L=0.1 K2=2.0 DY=1e-4);\nMULT TH = (L=0 K2=0.8 DY=5e-5);\nLINE CELL1 = (QF D2 B D2 TH QD D1);\nLINE RING = (ST 12*CELL1 RFC);\nFFS;\n$FORM="18.14";\nPageWidth = 1999;\nUSE RING;\nCELL;\nCALC;\nRFSW; RAD; NOFLUC;\nPOL;\nem = Emittance[];\nPrint["EM ",{SpinTune,NominalSpinTune}/.em];\nx0={1e-4,0};px0={0,0};y0={0,1e-4};py0={0,0};z0={0,0};dp0={0,1e-3};\nr=TrackParticles[{1, {x0,px0,y0,py0,z0,dp0,{0,0},{0.1,0.2},{1,1},{1,1}}}, 20, 1, 5];\nPrint["R ",Length[r]," ",Dimensions[r[[2]]]];\nPrint["S ",r[[2]]];\nExit;\n',
    [-0.000110919014271, 0.000798583261666, 1.098099837006e-05, -6.01915151646e-05, -6.31570214964e-08, -5.02468170226e-05, -3.35423074943e-08, 2.832270326749e-05, 0.0007961924513628, -0.00546786211618, -5.74132237886e-05, 0.0008836279001867, 2.9238646841851, 3.01345972439317, 0.10000000732172, 0.20014840635491, 1.0, 1.0])
CASES['pt8'] = ('MOMENTUM = 3GEV;\nDRIFT D1 = (L=1) D2 = (L=0.5);\nQUAD QF = (L=0.3 K1=0.7) QD = (L=0.3 K1=-0.6);\nBEND B = (L=1.2 ANGLE=0.2617993877991494 E1=0.1308996938995747 E2=0.1308996938995747);\nCAVI RFC = (VOLT=1.5e6 HARM=64);\nMARK ST = ();\nSOL SO = (L=0.3 BZ=0.05 BOUND=1 GEO=1);\nSOL SO2 = (L=0 BZ=0 BOUND=1);\nLINE CELL1 = (QF D2 B D2 QD D1);\nLINE RING = (ST 6*CELL1 SO D1 SO2 6*CELL1 RFC);\nFFS;\n$FORM="18.14";\nPageWidth = 1999;\nUSE RING;\nCELL;\nCALC;\nRFSW; RAD; NOFLUC;\nPOL;\nem = Emittance[];\nPrint["EM ",{SpinTune,NominalSpinTune}/.em];\nx0={1e-4,0};px0={0,0};y0={0,1e-4};py0={0,0};z0={0,0};dp0={0,1e-3};\nr=TrackParticles[{1, {x0,px0,y0,py0,z0,dp0,{0,0},{0.1,0.2},{1,1},{1,1}}}, 20, 1, 5];\nPrint["R ",Length[r]," ",Dimensions[r[[2]]]];\nPrint["S ",r[[2]]];\nExit;\n',
    [0.0004249162095486, 0.0003786208160105, 0.0001505578963919, 0.0001152721126939, 8.208999827024e-06, 0.0006137805350307, -7.26174663743e-06, -0.000580734479323, 0.0007955700987898, -0.00500516137404, -5.81431085112e-05, 0.0008887947163192, 2.92321944283329, 3.01210171490603, 0.09823359818979, 0.19488147117007, 1.0, 1.0])
CASES['pt13'] = ('MOMENTUM = 3GEV;\nDRIFT D1 = (L=1) D2 = (L=0.5);\nQUAD QF = (L=0.3 K1=0.7) QD = (L=0.3 K1=-0.6);\nBEND B = (L=1.2 ANGLE=0.2617993877991494 E1=0.1308996938995747 E2=0.1308996938995747);\nCAVI RFC = (VOLT=1.5e6 HARM=64);\nMARK ST = ();\nQUAD QT = (L=0.2 K1=0.02 ROTATE=0.5);\nSEXT SF = (L=0.1 K2=2.0 DY=1e-4);\nLINE CELL1 = (QF D2 B D2 SF QT QD D1);\nLINE RING = (ST 12*CELL1 RFC);\nFFS;\n$FORM="18.14";\nPageWidth = 1999;\nUSE RING;\nCELL;\nCALC;\nRFSW; RAD; FLUC;\nPOL;\nem = Emittance[];\nPrint["EM ",{SpinTune,NominalSpinTune}/.em];\nx0={1e-4,0};px0={0,0};y0={0,1e-4};py0={0,0};z0={0,0};dp0={0,1e-3};\nr=TrackParticles[{1, {x0,px0,y0,py0,z0,dp0,{0,0},{0.1,0.2},{1,1},{1,1}}}, 20, 1, 5];\nPrint["R ",Length[r]," ",Dimensions[r[[2]]]];\nPrint["S ",r[[2]]];\nExit;\n',
    [0.0004277990692368, 0.0009427696407098, -0.000123165961362, -0.000156807666783, -2.36655713711e-06, 0.0001228764976836, 1.349610306189e-05, 7.445825329957e-05, 0.0008747281552658, -0.00534818068285, -1.98455130668e-05, 0.0009744125477602, 1.14175586796584, 1.23348579649629, 0.10019398876307, 0.20009254017568, 1.0, 1.0])


def _run(src):
    from pysad.ffscmd import run_sad_string
    buf = io.StringIO()
    run_sad_string(src, out=buf)
    line = [l for l in buf.getvalue().splitlines() if l.startswith('S ')][0]
    return [float(x) for x in re.findall(r'[-+]?\d*\.?\d+(?:[Ee][-+]?\d+)?', line[1:].replace('{', ' ').replace('}', ' '))]


def _check(name, tol=1e-8):
    src, ref = CASES[name]
    got = _run(src)
    assert len(got) == len(ref)
    for a, b in zip(got, ref):
        assert abs(a - b) <= tol * max(abs(b), 1e-6) + 1e-14, (name, a, b)


def test_planar_ring_spin_tracking_matches_gs():
    _check('pt1')


def test_coupled_ring_spin_tracking_matches_gs():
    """縦ずれ・回転した四極と偏向・6 極のあるリング"""
    _check('pt4')


def test_spin_tracking_with_fluc_matches_gs():
    """RAD FLUC: 乱数の順と、光子が 0 個のときのスピン歳差 (K64 の tradkfn の elseif)"""
    _check('pt3')
    _check('pt13')


def test_spin_tracking_mult_bend_k1_sext_match_gs():
    """MULT (K1・K2・K3、tmulti の bsi)、K1 のある BEND (tbendi の bsi)、SEXT・薄い MULT"""
    for n in ('pt5', 'pt6', 'pt7', 'pt7b'):
        _check(n)


def test_spin_tracking_in_solenoid_region_matches_gs():
    """SOL 領域 (BOUND=1) 内のドリフト (tdrift_solenoid の bsi=bz·al)。K64 は SOL の長さを周長に入れない"""
    _check('pt8')
