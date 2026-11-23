# 09_dynamic_aperture.md から取り除いた履歴の記述

09_dynamic_aperture.md (現在の仕様の解説) から、版・経過に結びついた記述を取り除いた。下は取り除く前の元の行を、節の見出しつきでそのまま並べたもの (行番号は整理前のもの)。

## # 09. DynamicApertureSurvey — dapert.n → tfdapert → trackd と pysad の dapert.py

整理前の 5 行目:

````text
確認状況 (2026-11-17 時点): [`pysad/dapert.py`](../pysad/dapert.py) の最終更新は 2026-10-16、`bench_level16a.sad` は 2026-10-14。`Damping->True` は SADError、`DampEnough` は読まない (`dapert.py` の冒頭と `ffscmd.fn_dasurvey` で確認)。`ExpandElementValues`・`ReferenceOrbit` は v2.4.0 で読む (下記)。
````

## ## 5. pysad の対応 (dapert.py / ffscmd.fn_dasurvey)

整理前の 161 行目:

````text
**v2.4.0 で対応**: `ExpandElementValues->False` (成分ごとの値を揃えない。既定 True は揃える。gs の help: 'Machine errors may be reset') と `ReferenceOrbit->{x0,px0,y0,py0,z0,dp0}` (その軌道のまわりで追跡。`Default` は emit の閉軌道)。以前は常に揃えて ReferenceOrbit は読まなかったので、補正した COD・K0 を持つリングで gs と DA が違った (DDBA v07_14 の誤差の例: gs 34σ・29σ、以前の pysad は補正が消えた状態の 0〜38σ。`docs/history/v24_A.md` §A5、試験 `test_v24_A5_dapert_options.py`)。
````

