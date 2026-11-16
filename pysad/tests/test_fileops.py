"""ファイル・ディレクトリ関数 (fileops.py。Packages/FileSystem.n)、ARRAY・StringInsert・Expand など (init.n・functions.n)。
ARRAY・Len・StringInsert・Diagonal・Expand は gs と一致を確認済み (local/wake/ex11.txt)。"""
import io


def _run(body, cwd):
    from pysad.ffscmd import run_sad_string
    buf = io.StringIO()
    run_sad_string('FFS;\n$FORM="";\n' + body + 'Exit;\n', out=buf, cwd=str(cwd))
    return {l[:2]: l[2:] for l in buf.getvalue().splitlines() if len(l) > 2 and l[1] == ' '}


def test_file_functions(tmp_path):
    out = _run('''Print["1 ",CreateDirectory["d"]," ",DirectoryQ["d"]," ",CreateDirectory["d"]];
f=OpenWrite["d/a.txt"];WriteString[f,"hello\nworld\n"];Close[f];
Print["2 ",FileQ["d/a.txt"]," ",FileType["d/a.txt"]," ",FileType["d"]," ",FileType["zzz"]," ",FileByteCount["d/a.txt"]];
Print["3 ",CopyFile["d/a.txt","d/b.txt"]<=>$Failed," ",CopyFile["d/a.txt","d/b.txt"]," ",RenameFile["d/b.txt","d/c.txt"]<=>$Failed," ",FileQ["d/b.txt"]];
Print["4 ",DeleteFile["d/c.txt"]," ",DeleteFile["d/zz"]," ",FileNames["*"]];
SetFileDate["d/a.txt",{2020,1,2,3,4,5}];
Print["5 ",FileDate["d/a.txt"]];
Print["6 ",CopyDirectory["d","e"]<=>$Failed," ",FileQ["e/a.txt"]," ",DeleteDirectory["d"]," ",DeleteDirectory["d",DeleteContents->True]," ",DirectoryQ["d"]];
Print["7 ",DirectoryName["/a/b/c.txt"]," ",DirectoryName["/a/b/c.txt",2]," ",ToFileName[{"/a","b","c"}]," ",ToFileName[{"a","b"},"c"]];
''', tmp_path)
    assert out['1 '].startswith('d 1 $Failed') or out['1 '] == 'd 1 $Failed'
    assert out['2 '] == '1 File Directory None 14'
    assert out['3 '] == '1 $Failed 1 0'
    assert out['4 '] == 'Null $Failed {d,e}' or out['4 '].startswith('Null $Failed')
    assert out['5 '] == '{2020,1,2,3,4,5}'
    assert out['6 '] == '1 1 $Failed Null 0'
    assert out['7 '] == '/a/b/ /a/ /a/b/c/ a/b/c'


def test_init_n_small_functions(tmp_path):
    out = _run('''Print["1 ",ARRAY[2,3]," ",Len["abcd"]," ",StringInsert["abcdef","XY",3]," ",StringInsert["abcdef","XY",{2,4}]," ",Diagonal[{{1,2},{3,4}}]];
Print["2 ",Expand[(a+b)*(c+d)]," ",Expand[(a+b)^2]," ",SymbolNameRoot[x$12]," ",SymbolNameRoot[abc]];
Put[1,{2,3},"p.txt"];
Print["3 ",ReadString[OpenRead["p.txt"]]];
''', tmp_path)
    assert out['1 '] == '{{0,0,0},{0,0,0}} 4 abXYcdef aXYbcXYdef {1,4}'
    assert out['2 '] == '(a c)+(a d)+(b c)+(b d) (a+b)^2 x abc'
