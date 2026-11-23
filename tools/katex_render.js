// 標準入力の JSON [{tex, display}] を KaTeX で HTML にして、標準出力に JSON [{html, err}] で返す。
// 第 1 引数 (省略可): output ('html' | 'htmlAndMathml')。変換できない数式は、赤字の TeX で残し err に理由を入れる。
const fs = require('fs');
const path = require('path');
const katex = require(path.join(__dirname, 'katex', 'katex.min.js'));

const output = process.argv[2] || 'html';
// 物理でよく使う記法は KaTeX の標準で足りる。足りないものはここに書く: '\\name': '定義'
const macros = { '\\lambdabar': '\\bar{\\lambda}' };   // 換算コンプトン波長の ƛ (KaTeX に無い。λ に横棒で代用)
const base = { output, macros, strict: 'ignore', trust: false, maxExpand: 1000, errorColor: '#cc0000' };

const items = JSON.parse(fs.readFileSync(0, 'utf8'));
const out = items.map((it) => {
  const opt = Object.assign({}, base, { displayMode: !!it.display });
  try {
    return { html: katex.renderToString(it.tex, Object.assign({}, opt, { throwOnError: true })), err: null };
  } catch (e) {
    let html = '';
    try { html = katex.renderToString(it.tex, Object.assign({}, opt, { throwOnError: false })); } catch (e2) { html = ''; }
    return { html: html, err: String(e && (e.rawMessage || e.message) || e) };
  }
});
process.stdout.write(JSON.stringify(out));
