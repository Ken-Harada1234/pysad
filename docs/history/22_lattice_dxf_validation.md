# ラティスの CAD 図 (DXF) の検証の詳細 (gs との比較)

本文 [../22_lattice_dxf.md](../22_lattice_dxf.md) から移した、K64 (gs) との一致の記述。元の文書: 51_lattice_dxf.md §4 (本文の §7)。

- **gs との一致**: 同じ台本 `DrawMagnetsDXF[]` を gs (SAD 1.4.4.2k64、WSL) と pysad で走らせ、8 本の DXF (FODO 3、SOL 1、ASTRID2 2、KEK-LS 2) が**全行 (文字列として) 一致**した。pysad の Python 関数 `write_dxf` の出力も同じ。再現は `py -3.12 example/dxf_lattice/run_dxf_examples.py`。台本と Python の一致は pysad だけで走る試験 (`test_script_matches_python`、3 通りのオプション) にもある。
  - **gs との一致**: 同じ台本 (`DXF_FODO_3d.sad`・`DXF_KEKLS_3d.sad`) を gs (SAD 1.4.4.2k64、WSL。`local/wake/gsrun.sh` の呼び出し) と pysad で走らせ、5 本の 3D の DXF (FODO の 1 セル `body`・`split`・全周、KEK-LS の 1 セル `body`・`split`) が**全行一致** (gs は行末が LF、pysad は CRLF の違いだけ)。pysad の台本と Python 関数 `write_dxf(..., mode='3d')` も全行一致 (`run_dxf_examples_3d.py run`)。試験にも pysad の台本と Python 関数の一致 (`body`・`split`・`MagnetSpec`・`Scale`・`Overlap`) がある。
  - **2D の既定の出力**: 現在の `DrawMagnetsDXF.n` と `write_dxf` で、既存の DXF 8 本と v2 の 5 本 (台本のとおりに作り直したもの) が gs・pysad の台本とも**全行一致**。
  - **gs との一致**: KEK-LS の標準セルの 5 本 (`MagnetSpec`・`Overlap`・`UseMagnetSpec` を組み合わせた `DXF_KEKLS_cell_v2_*.dxf`) と全周 1 本 (`Overlap->"remove"`) が、gs・pysad の台本・Python 関数 (台本から読む場合と dict で与える場合。全周は Python 関数と pysad の台本) で全行一致。pysad の台本と Python 関数の一致は 4 通りのオプションの試験にもある。再現は `py -3.12 example/dxf_lattice/run_dxf_examples_v2.py` (既存の DXF・PNG は上書きしない)。
  - 既定のまま (`overlap="keep"`、`MagnetSpec` なし) の出力は、既存の `DXF_KEKLS_cell.dxf` と全行一致 (`test_example_v2_files`)。
