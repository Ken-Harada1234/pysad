# -*- coding: utf-8 -*-
"""PF の履歴ラティス 4 本 (example/pf/lattice/pf_*.sad) を、元の履歴フォルダのファイルから作る。
使い方: py -3.12 example/tools/pf_history/mk_lattices.py <履歴フォルダ>
元のファイルは読むだけ。内容は改行を LF にして、先頭に新しいヘッダを付け、元のヘッダはその下に「元のヘッダ」として残す。
本体の変更は、2005 の `INJ` を MARK に足す 1 行だけ (元は未定義で、gs が警告を出す)。"""
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__)); EX = os.path.dirname(os.path.dirname(HERE))
OUT = os.path.join(EX, 'pf', 'lattice')

COMMON = """!  エネルギー 2.5 GeV (MOMENTUM=2.5 GEV)、偏向磁石 28 台 (矩形、E1=E2=0.5)、空洞 CAVI 2 台 (各 VOLT 850 kV)。
!  元は設計値のラティス (誤差なし)。数値は ON RAD FLUC COD RFSW RADCOD; のあとの `cell calc emit` の値
!  (gs と pysad の k64 モードで全桁一致。example/pf/lattice/PF_history.md の表)。
"""
SPECS = [
 dict(src='pf_460nmrad', dst='pf_1981_460nmrad.sad', name='pf_1981_460nmrad', era='1981 年',
  what='PF リング建設時のオプティクス (中エネルギー・低 beta・高 nuy のオプティクス)',
  numbers="""!  周長      187.0741 m (RING の 2*HALF = リング全体。HALF は半周)
!  チューン  nux = 5.2987, nuy = 4.1638 (cell calc の周期解。cell calc emit は水平安定)
!  エミッタンス  名前は 460 nm rad、計算値は 4.034e-7 m rad (403 nm rad)。名前と合わない (要確認)
!            元のヘッダの NUX=8.454 NUY=3.300 も合わない (実測 5.2987 / 4.1638。この文面は 1986 年のラティス
!            pf_128nmrad のヘッダと同じで、そちらでは一致する。要確認)""",
  tidy="""!  整理: なし (台本の中身は元のまま。改行を LF にしただけ)。
!        RING = (IMARK III 2*HALF RMARK)、HALF = (HMARK -RQUAD LQUAD) で周回のラインになっている (Get・USE・計算コマンドは無い)。
!        -RQUAD の中の -RQUADM (二重の反転) を pysad が正しく展開しなかったので、pysad/parser.py を直した (gs は元から正しい)。"""),
 dict(src='pf_128nmrad', dst='pf_1986_128nmrad.sad', name='pf_1986_128nmrad', era='1986 年',
  what='最初の高輝度化のオプティクス (PF 履歴フォルダの pf0 と同じ内容)',
  numbers="""!  周長      187.0741 m (RING の 2*HALF = リング全体。HALF は半周)
!  チューン  nux = 8.45400, nuy = 3.30000 (元のヘッダの NUX=8.454 NUY=3.300 と一致)
!  エミッタンス  名前は 128 nm rad、計算値は 1.213e-7 m rad (121 nm rad)。名前と合わない (要確認)""",
  tidy="""!  整理: なし (台本の中身は元のまま。改行を LF にしただけ)。
!        RING = (IMARK III 2*HALF RMARK)、HALF = (HMARK -RQUAD LQUAD) で周回のラインになっている (Get・USE・計算コマンドは無い)。
!        -RQUAD の中の -RQUADM (二重の反転) を pysad が正しく展開しなかったので、pysad/parser.py を直した (gs は元から正しい)。"""),
 dict(src='pf_036nmrad', dst='pf_1997_036nmrad.sad', name='pf_1997_036nmrad', era='1997 年',
  what='弧部改造・高輝度化 (超低エミッタンス) のオプティクス (標準セルの位相進み 105 度、Version 3.6.1、2000-02-22)',
  numbers="""!  周長      187.0741 m
!  チューン  nux = 9.5972, nuy = 4.2902
!  エミッタンス  名前は 36 nm rad、計算値は 3.473e-8 m rad (34.7 nm rad)。近い (差 3.5%)""",
  tidy="""!  整理: なし (台本の中身は元のまま。改行を LF にしただけ)。Get・USE・計算コマンドは無い。
!        周回のライン RING = (RMARK HALFE HALFW EMARK)。四極 QUAD と六極 SEXT の定義が何度も出てくるのは機械研究ごとの
!        値の履歴 (2000-05-29、2002-11-25、`105 deg Present Optics`) で、あとの定義が有効 (最後が `105 deg Present Optics`)。"""),
 dict(src='pf_036nmrad_2005', dst='pf_2005_036nmrad.sad', name='pf_2005_036nmrad', era='2005 年',
  what='直線部増強のオプティクス (標準セルの位相進み 115 度。PM と STR 付きの最初の版 = First edition)',
  numbers="""!  周長      187.0355 m
!  チューン  nux = 9.5975, nuy = 5.2790
!  エミッタンス  名前は 36 nm rad、計算値は 3.454e-8 m rad (34.5 nm rad)。近い (差 4%)
!  参考: 同じ PF の `v9_16_asym.sad` (Second edition) と周長・チューンは全桁同じ、エミッタンスも相対差 4e-7 で同じ""",
  tidy="""!  整理: MARK の一覧に `INJ =()` を 1 行足した (元は STRNE の中で `INJ` を使うのに定義が無く、gs が
!        `INJ is not defined yet` の警告を出す。pysad は黙って飛ばす。光学には効かない)。それ以外は元のまま。
!        TEST・NCELL (NW・SW が未定義) など使われない定義は元のまま残した。周回のライン RING = (HALFW HALFE RMARK)。"""),
]

def main(srcdir):
    for sp in SPECS:
        t = open(os.path.join(srcdir, sp['src']), 'rb').read().decode('ascii').replace('\r\n', '\n').replace('\r', '\n')
        if not t.endswith('\n'): t += '\n'
        # 元のヘッダ = 先頭の ! で始まる連続した行
        lines = t.split('\n'); k = 0
        while k < len(lines) and lines[k].startswith('!'): k += 1
        hdr, body = lines[:k], '\n'.join(lines[k:])
        if sp['src'] == 'pf_036nmrad_2005':
            old = '        INJ' 
            assert 'INJ   =' not in body and ' MARK   MKD   =()' in body
            body = body.replace(' MARK   MKD   =()\n', ' MARK   MKD   =()\n        INJ   =()\n', 1)
        head = ['!' + '=' * 76,
                '!  %s.sad : PF リング (2.5 GeV) の %s のオプティクス' % (sp['name'], sp['era']),
                '!  %s' % sp['what'],
                '!' + '-' * 76,
                '!  元のファイル: %s' % sp['src']]
        txt = '\n'.join(head) + '\n' + sp['numbers'] + '\n' + COMMON.rstrip('\n') + '\n' + sp['tidy'] + '\n' + \
              '!' + '=' * 76 + '\n!  ---- 元のヘッダ (%s の先頭) ----\n' % sp['src'] + '\n'.join(hdr) + '\n' + body
        open(os.path.join(OUT, sp['dst']), 'w', encoding='utf-8', newline='\n').write(txt)
        print(sp['dst'], len(txt))
if __name__ == '__main__':
    main(sys.argv[1])
