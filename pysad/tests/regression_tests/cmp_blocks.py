# -*- coding: utf-8 -*-
"""cmp_blocks.py SAD.result PYSAD.out [TAG正規表現] [--num tol] [--noise x] — Print["Dn ..."] で区切られたブロックごとに、出力行を突き合わせる。

ブロック: 行頭 (空白なし) が TAG正規表現 (既定 `[A-Z]+[0-9]+[a-z]* `) で始まる行から次のタグ行まで。
SAD の結果ファイルにはコマンドのエコー行 (先頭 1 空白で始まり ';' で終わる行) があるので、それは両方から除く。
既定は完全一致 (末尾の空白は無視)。--num tol を付けると、文字部分 (数を '#' にした形) が一致し数値が相対 tol 以内なら一致扱い。
--noise x を付けると、両方の絶対値が x 未満の数 (SAD の 1e-17 などの丸め雑音と pysad の 0) は一致扱い。(その場合、数と数の間の空白の数は問わない)。不一致のブロックは最初の差の行を表示。終了コード: 不一致があれば 1。"""
import re
import sys

NUM = re.compile(r'[-+]?(?:\d+\.\d*|\.\d+|\d+)(?:[eEdD][-+]?\d+)?')


def read_blocks(path, tag):
    lines = open(path, encoding='utf-8', errors='replace').read().split('\n')
    blocks, cur, name = [], None, None
    for l in lines:
        if tag.match(l):
            name = l.split()[0] if l.split() else l
            cur = []
            blocks.append((l.strip(), cur))
            continue
        if cur is None:
            continue
        s = l.rstrip()
        if s.startswith(' ') and s.endswith(';') and not s.startswith('  '):
            continue                       # コマンドのエコー
        if not s.strip():
            continue
        cur.append(s)
    return blocks


def same(a, b, tol, noise=None):
    """a, b が一致するか。数を除いた文字部分は完全一致が必要。数は、文字列が同じ / 両方 |x|<noise (0 同然の雑音) / 相対 tol 以内 のどれか"""
    if a == b:
        return True
    if tol is None and noise is None:
        return False
    ta, tb = NUM.split(a), NUM.split(b)
    if noise is not None:                  # 雑音の桁は幅が違う (.0 と .00000) ので空白の数は問わない
        ta, tb = [' '.join(x.split()) for x in ta], [' '.join(x.split()) for x in tb]
    na, nb = NUM.findall(a), NUM.findall(b)
    if ta != tb or len(na) != len(nb):
        return False
    for x, y in zip(na, nb):
        if x == y:
            continue
        fx, fy = float(x.replace('D', 'e')), float(y.replace('D', 'e'))
        if noise is not None and abs(fx) < noise and abs(fy) < noise:
            continue
        if tol is not None and abs(fx - fy) <= tol * max(abs(fx), abs(fy)):
            continue
        return False
    return True


def compare(sad_path, py_text, tag, tol=None, noise=None, skip=()):
    """run_regression 用: ブロックごとに比べて (ブロック数, 不一致のタグのリスト) を返す。skip = 比べないタグの先頭文字列"""
    import tempfile
    A = read_blocks(sad_path, tag)
    with tempfile.NamedTemporaryFile('w', suffix='.txt', delete=False, encoding='utf-8') as f:
        f.write(py_text)
    try:
        B = read_blocks(f.name, tag)
    finally:
        import os
        os.unlink(f.name)
    pb = {t: ls for t, ls in B}
    bad = []
    n = 0
    for t, la in A:
        if any(t.startswith(s) for s in skip):
            continue
        n += 1
        lb = pb.get(t)
        if lb is None or len(la) != len(lb) or any(not same(x, y, tol, noise) for x, y in zip(la, lb)):
            bad.append(t.split()[0])
    return n, bad


def main(argv):
    tol = None
    noise = None
    args = []
    i = 1
    while i < len(argv):
        if argv[i] == '--num':
            tol = float(argv[i + 1])
            i += 2
        elif argv[i] == '--noise':
            noise = float(argv[i + 1])
            i += 2
        else:
            args.append(argv[i])
            i += 1
    sad, py = args[0], args[1]
    tag = re.compile(args[2] if len(args) > 2 else r'[A-Z]+[0-9]+[a-z]* ')
    A, B = read_blocks(sad, tag), read_blocks(py, tag)
    bad = 0
    pb = {t: ls for t, ls in B}
    for t, la in A:
        if t not in pb:
            print(f'{t}: pysad にブロックがない')
            bad += 1
            continue
        lb = pb[t]
        n = max(len(la), len(lb))
        diff = [(k, la[k] if k < len(la) else '<なし>', lb[k] if k < len(lb) else '<なし>') for k in range(n)
                if k >= len(la) or k >= len(lb) or not same(la[k], lb[k], tol, noise)]
        if diff:
            bad += 1
            k, x, y = diff[0]
            print(f'{t}: {len(diff)}/{n} 行が不一致 (行数 SAD {len(la)} / pysad {len(lb)}); 最初 (行 {k + 1}):')
            print('   SAD  :', x[:170])
            print('   pysad:', y[:170])
        else:
            print(f'{t}: 一致 ({len(la)} 行)')
    print(f'ブロック {len(A)} 個, 不一致 {bad}')
    return 1 if bad else 0


if __name__ == '__main__':
    sys.exit(main(sys.argv))
