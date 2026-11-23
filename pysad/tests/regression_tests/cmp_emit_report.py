# -*- coding: utf-8 -*-
"""EMIT 報告の行単位比較: python cmp_emit_report.py SAD.result PYSAD.result
'Closed orbit:' から 'Nominal spin tune' 行までのブロックを順に取り出して行ごとに比べる。
完全一致、または数値の差 |差|<=1e-11 (0 同然の雑音) か相対差<=1e-9 で、数値以外の部分 (空白は無視) が同じなら一致。
違えば行を並べて表示 (数値は最大相対差も)。run_regression.py からは blocks()/same() を使う。"""
import re
import sys

NUM = re.compile(r'[-+]?(?:\d+\.\d*|\.\d+|\d+)(?:[eE][-+]?\d+)?')


def blocks(lines):
    """行の列 → EMIT 報告のブロック (行のリスト) のリスト"""
    out, cur = [], None
    for l in lines:
        l = l.rstrip('\n')
        if l.strip() == 'Closed orbit:':
            cur = [l]
            continue
        if cur is not None:
            cur.append(l)
            if l.startswith('Nominal spin tune'):
                out.append(cur)
                cur = None
    return out


def load_blocks(path):
    return blocks(open(path, encoding='utf-8', errors='replace').readlines())


def same(lx, ly, noise=1e-11, rel=1e-9):
    """True なら一致 (0 を雑音として許す)。戻り値は (一致か, 雑音だけの差か)"""
    if lx == ly:
        return True, False
    nx, ny = NUM.findall(lx), NUM.findall(ly)
    # 数値以外の並びは空白を全部除いて比べる ('tune:-0.0000000' と 'tune: 0.0000000' のように負の 0 で空白の位置が変わる)
    if len(nx) != len(ny) or ''.join(NUM.sub('#', lx).split()) != ''.join(NUM.sub('#', ly).split()):
        return False, False
    # 0 同然の量 (Beam tilt など、どちらも |v|<1e-9) は丸め雑音として同じとみなす (方針: 1e-10 以下の差は追わない)
    ok = all(abs(float(p) - float(q)) <= noise or abs(float(p) - float(q)) <= rel * max(abs(float(p)), abs(float(q)))
             or (abs(float(p)) < 1e-9 and abs(float(q)) < 1e-9)
             for p, q in zip(nx, ny))
    return ok, ok


def main():
    a = load_blocks(sys.argv[1])
    b = load_blocks(sys.argv[2])
    print('blocks', len(a), len(b))
    nbad = nnoise = 0
    for k, (x, y) in enumerate(zip(a, b)):
        print(f'--- block {k + 1}: lines {len(x)} / {len(y)}')
        for i in range(max(len(x), len(y))):
            lx = x[i] if i < len(x) else '<none>'
            ly = y[i] if i < len(y) else '<none>'
            ok, noisy = same(lx, ly)
            if ok:
                nnoise += 1 if noisy else 0
                continue
            nbad += 1
            nx, ny = NUM.findall(lx), NUM.findall(ly)
            r = ''
            if len(nx) == len(ny) and nx:
                d = [abs(float(p) - float(q)) / max(abs(float(p)), abs(float(q)), 1e-300) for p, q in zip(nx, ny)]
                r = ' maxrel=%.1e' % max(d)
            print(f'  L{i}{r}\n   SAD  : {lx}\n   pysad: {ly}')
    print('mismatched lines:', nbad, ' noise-only lines (|diff|<=1e-11 or rel<=1e-9):', nnoise)
    return 1 if nbad else 0


if __name__ == '__main__':
    sys.exit(main())
