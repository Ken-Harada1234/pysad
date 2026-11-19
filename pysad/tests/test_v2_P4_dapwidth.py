# -*- coding: utf-8 -*-
"""v2next 担当 P4 項目 5: DynamicApertureSurvey の DAPWIDTH の窓 (trackd.f の粒子スロットの管理の直訳)。参照は gs (K64、2026-10-19 ビルド、2026-11-19)。

gs の trackd は、行 (ivar1 の値) ごとに、損失が見つかった最小の j (nxm) から下へ DAPWIDTH 点の窓の中だけを追跡し、同時に追うスロットは
npmax = n1p·DAPWIDTH 個 (NPARA=1)。損失・完了でスロットが空くと配置ループに戻り、空きがあるときは 16 周ごとに戻る。窓の外の点は追跡せず、
ntloss = 2^29 のまま (地図は `*`)。以前の pysad は常に全 51 点を追い (DAPWIDTH 既定 51)、`Maximum number of particles` も違っていた (gs の既定は 7)。
どのモードでも gs の動作 (k64 の再現)。使った条件: 3 GeV、24 セルの FODO (周長 36 m、SF・SD 付き、RF 4 MV)、NORAD と RAD + FLUC (NPARA=1)。
NPARA > 1 (子プロセスに分ける) は未対応。損失の境界の粒子 (カオス的) は gs と 1 周ほど違うことがある (下の許容)。
"""
import io
import re

import pytest

_HEAD = '''MOMENTUM = 3GEV;
DRIFT D1 = (L=1.0) D2 = (L=0.5) DS = (L=0.1);
QUAD  QF = (L=0.3 K1=0.7) QD = (L=0.3 K1=-0.6);
SEXT  SF = (L=0.2 K2=4) SD = (L=0.2 K2=-7);
BEND  B  = (L=1.2 ANGLE=0.2617993877991494 E1=0.1308996938995747 E2=0.1308996938995747);
MARK  ORIGIN = ();
CAVI RF = (L=0 VOLT=4.0E6 FREQ=5.0E8 PHI=0);
LINE CELL = (SF DS QF D2 B D2 QD DS SD D1);
LINE RING = (ORIGIN RF 24*CELL);
FFS;
USE RING; CELL; CALC;
RFSW; RAD; NOFLUC;
EMIT;
NPARA = 1;
$FORM = "12.6";
'''


def _run(body):
    from pysad.ffscmd import run_sad_string
    buf = io.StringIO()
    run_sad_string(_HEAD + body + '\nExit;\n', out=buf)
    return buf.getvalue()


def _maps(out):
    """地図の行 (a1, nsc, 51 文字) の並び。呼び出しごとの塊に分ける"""
    blocks, cur = [], None
    for l in out.splitlines():
        if l.startswith(' Turns ='):
            cur = []
            blocks.append(cur)
        m = re.match(r'^\s+(-?\d+\.\d+)\s{0,2}(\d+) (.{51})\s*$', l)
        if m and cur is not None:
            cur.append((float(m.group(1)), int(m.group(2)), m.group(3)))
    return blocks


def _seeds(out):
    return [float(m) for m in re.findall(r'SEED[^,]*, ([0-9.]+)', out)]


_DA = 'res = DynamicApertureSurvey[{{1,300},{1,1},{-0.5,0,0.5}}, 100, Output->6];\n'


@pytest.mark.parametrize('dw,npmax,head', [(None, 42, '*UUUUUUU'), (3, 18, '*****UUU'), (12, 72, 'UUUUUUUU'), (51, 306, 'UUUUUUUU')])
def test_dapwidth_window_maps_match_gs(dw, npmax, head):
    """NORAD の 3 行 (z = −0.5・0・0.5)、100 周。gs: スコア 24 (行ごと nsc = 8)。窓の外 (j = 1 から) は `*`。DAPWIDTH の既定は 7 で
    `Maximum number of particles` は n1·min(DAPWIDTH, 51)·2 (gs: 42・18・72・306)。損失の境界の 1 点 (j = 14) は gs と 1 周違いうるので先頭 13 文字で比べる"""
    body = ('NORAD;\n' + ('DAPWIDTH = %d;\n' % dw if dw is not None else '') + _DA)
    out = _run(body)
    assert ' Turns =  100.000000 Maximum number of particles =  %s' % ('%d.000000' % npmax).rjust(10) in out or \
        re.search(r'Maximum number of particles =\s+%d\.0' % npmax, out)
    m = _maps(out)[0]
    assert [r[1] for r in m] == [8, 8, 8]
    gs_rows = {'-0.50': 'BAnUg', '0.00': 'L9UUj', '0.50': 'b87UJ'}
    for a1, nsc, row in m:
        assert row[:8] == head
        assert row[8:13] == gs_rows['%.2f' % a1]


# gs (NPARA=1、RAD + FLUC、DAPWIDTH = 2、ntloss の最終 seed): SeedRandom[200 + nt] から。y = 1、x = {1, 150}、z = {−0.5, 0, 0.5, 1}。
# 窓の動き (スロットの空き・詰め替え・16 周ごとの再配置) が gs と同じなら乱数の列が同じで、最終の seed まで一致する
_GS_SEEDS_DW2 = {5: (111.0, 43421809.0), 10: (108.0, 1130603987.0), 20: (108.0, 3849420681.0), 25: (108.0, 2237612085.0)}


@pytest.mark.parametrize('nt', sorted(_GS_SEEDS_DW2))
def test_dapwidth_window_slot_dynamics_reproduce_gs_random_stream(nt):
    body = ('RAD; FLUC;\nDAPWIDTH = 2;\nSeedRandom[%d];\n' % (200 + nt)
            + 'res = DynamicApertureSurvey[{{1,150},{1,1},{-0.5,0,0.5,1.0}}, %d, Output->0];\n' % nt
            + 'Print["SC ", res[[1]], " SEED ", SeedRandom[]];\n')
    out = _run(body)
    m = re.search(r'SC\s+([0-9.]+) SEED \{"SAD",\s+([0-9.]+)', out)
    score, seed = float(m.group(1)), float(m.group(2))
    assert (score, seed) == _GS_SEEDS_DW2[nt]


def test_dapwidth_one_slot_and_single_row_seed_like_gs():
    """DAPWIDTH = 1 (スロット 1 つ。点を 1 つずつ追う)、x = {1, xm}、z = 0、2 周: gs の最終 seed と一致 (xm = 60・100。ほかの xm も 130 まで一致、
    xm = 140 だけ損失の境界の粒子が 1 周違って seed が違う)"""
    ref = {60: (51.0, 1684391341.0), 100: (43.0, 2339435309.0)}
    for xm, (score, seed) in ref.items():
        body = ('RAD; FLUC;\nDAPWIDTH = 1;\nSeedRandom[300];\n'
                + 'res = DynamicApertureSurvey[{{1,%d},{1,1},{0.0}}, 2, Output->0];\n' % xm
                + 'Print["SC ", res[[1]], " SEED ", SeedRandom[]];\n')
        m = re.search(r'SC\s+([0-9.]+) SEED \{"SAD",\s+([0-9.]+)', _run(body))
        assert (float(m.group(1)), float(m.group(2))) == (score, seed), xm


def test_dapwidth_result_of_dapert_marks_unlaunched_points_with_2_pow_29():
    out = _run('NORAD;\nDAPWIDTH = 3;\n' + _DA + 'Print["R1 ", res[[2,2,1,3]]];\n')
    row = [float(v) for v in re.search(r'R1 \{([^}]*)\}', out.replace('\n', ' ')).group(1).split(',')]
    assert len(row) == 51
    assert row[:5] == [2.0 ** 29] * 5 and row[5] == 100.0 and row[-1] == 0.0


_RING_A_GS = {
    7: ['**************ZZZZZZZ1112225ZZZZ212121', '**************ZZZZZZZ1117135VZZZ222111', '**************ZZZZZZZ23 24234ZZZ32E31 '],
    51: ['ZZZZZZZZZZZZZZZZZZZZZ1112225ZZZZ212121', 'ZZZZZZZZZZZZZZZZZZZZZ1117135VZZZ222111', 'ZZZZZZZZZZZZZZZZZZZZZ23 24234ZZZ32D31 '],
}


def test_dapwidth_window_on_large_ring_a_matches_gs():
    """大きなリング A (`kek-ls_v03_12.lat`、3 GeV、周長 603 m、3563 要素、RFSW + NORAD): `{{1,400},{1,1},{-0.5,0,0.5}}`、60 周で、DAPWIDTH = 7 (既定: j = 1〜14 は窓の外で `*`) と 51 の
    地図 (21 文字ずつ。それ以降は 0 周で空白) と Score 63 が gs と一致 (損失の境界の 1 点 (z = 0.5 の j = 34) だけ 1〜2 周違い、文字が違う)"""
    from pathlib import Path
    p = Path(__file__).resolve().parent / 'kek-ls_v03_12.lat'
    if not p.exists():
        pytest.skip('大きなリング A のラティスがない')
    from pysad.ffscmd import run_sad_string
    tail = '''
FFS;
NPARA = 1;
USE RING; CELL; CALC;
RFSW; RAD; NOFLUC;
EMIT;
NORAD;
$FORM = "12.6";
res = DynamicApertureSurvey[{{1,400},{1,1},{-0.5,0,0.5}}, 60, Output->6];
Print["SCORE ", res[[1]]];
DAPWIDTH = 51;
res = DynamicApertureSurvey[{{1,400},{1,1},{-0.5,0,0.5}}, 60, Output->6];
Print["SCORE51 ", res[[1]]];
Exit;
'''
    buf = io.StringIO()
    run_sad_string(p.read_text(encoding='utf-8', errors='replace') + tail, out=buf)
    out = buf.getvalue()
    maps = _maps(out)
    assert len(maps) == 2
    for blk, dw in zip(maps, (7, 51)):
        assert [r[1] for r in blk] == [21, 21, 21]
        mism = 0
        for (a1, nsc, row), ref in zip(blk, _RING_A_GS[dw]):
            row = row[:len(ref)]
            mism += sum(1 for a, b in zip(row, ref) if a != b)
        assert mism <= 1, (dw, blk)
    assert [r[2][:14] for r in maps[0]] == ['*' * 14] * 3 and all('*' not in r[2] for r in maps[1])
    assert 'SCORE    63.000000' in out and 'SCORE51    63.000000' in out
