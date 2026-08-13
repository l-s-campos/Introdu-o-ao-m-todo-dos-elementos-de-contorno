/**
 * Build multipage GitHub Pages site from Typst HTML exports.
 */
import { execFileSync } from "child_process";
import {
  mkdirSync,
  readFileSync,
  writeFileSync,
  existsSync,
  readdirSync,
  rmSync,
} from "fs";
import { join, dirname } from "path";
import { fileURLToPath } from "url";

const __dirname = dirname(fileURLToPath(import.meta.url));
const ROOT = __dirname;
const DOCS = join(ROOT, "docs");
const RAW = join(ROOT, "docs", "_raw");

const PAGES = [
  { id: "index", title: "Início", file: "index.html", src: "site/pages/index.typ" },
  { id: "apresentacao", title: "Apresentação", file: "apresentacao.html", src: "site/pages/apresentacao.typ" },
  { id: "glossario", title: "Glossário", file: "glossario.html", src: "site/pages/glossario.typ" },
  { id: "interpolacao", title: "Interpolação", file: "interpolacao.html", src: "site/pages/interpolacao.typ" },
  { id: "equacoes-diferenciais", title: "Equações diferenciais", file: "equacoes-diferenciais.html", src: "site/pages/equacoes-diferenciais.typ" },
  { id: "indo-para-2d", title: "Indo para 2D", file: "indo-para-2d.html", src: "site/pages/indo-para-2d.typ" },
  { id: "gmsh", title: "Gmsh e malhas", file: "gmsh.html", src: "site/pages/gmsh.typ" },
  { id: "laplace-2d", title: "Laplace 2D", file: "laplace-2d.html", src: "site/pages/laplace-2d.typ" },
  { id: "poisson-2d", title: "Poisson 2D", file: "poisson-2d.html", src: "site/pages/poisson-2d.typ" },
  { id: "elasticidade-2d", title: "Elasticidade 2D", file: "elasticidade-2d.html", src: "site/pages/elasticidade-2d.typ" },
  { id: "propgeo-3d", title: "Propgeo 3D", file: "propgeo-3d.html", src: "site/pages/propgeo-3d.typ" },
  { id: "trabalhos-finais", title: "Trabalhos finais", file: "trabalhos-finais.html", src: "site/pages/trabalhos-finais.typ" },
  { id: "viga-euler", title: "Viga de Euler (extra)", file: "viga-euler.html", src: "site/pages/viga-euler.typ" },
  { id: "extra", title: "Extra", file: "extra.html", src: "site/pages/extra.typ" },
];

const CSS = `/* site chrome */
:root {
  --bg: #f7f4ef;
  --surface: #ffffff;
  --ink: #1c1917;
  --muted: #57534e;
  --accent: #0f766e;
  --accent-soft: #ccfbf1;
  --line: #e7e5e4;
  --code-bg: #f5f5f4;
  --nav-w: 280px;
  --shadow: 0 10px 30px rgba(28, 25, 23, 0.06);
}
* { box-sizing: border-box; }
html { scroll-behavior: smooth; }
body {
  margin: 0;
  font-family: "Segoe UI", "Helvetica Neue", Arial, sans-serif;
  color: var(--ink);
  background:
    radial-gradient(circle at top left, #ecfeff 0, transparent 28%),
    radial-gradient(circle at top right, #fff7ed 0, transparent 24%),
    var(--bg);
  line-height: 1.65;
}
a { color: var(--accent); }
a:hover { color: #115e59; }
.layout {
  display: grid;
  grid-template-columns: var(--nav-w) 1fr;
  min-height: 100vh;
}
.sidebar {
  position: sticky;
  top: 0;
  height: 100vh;
  overflow: auto;
  padding: 1.5rem 1rem 2rem;
  background: rgba(28, 25, 23, 0.94);
  color: #fafaf9;
  border-right: 1px solid rgba(255,255,255,0.06);
}
.brand {
  display: block;
  color: #fff !important;
  text-decoration: none !important;
  font-weight: 700;
  font-size: 1.05rem;
  line-height: 1.35;
  margin-bottom: 0.35rem;
}
.brand-sub {
  color: #a8a29e;
  font-size: 0.85rem;
  margin-bottom: 1.25rem;
}
.nav-label {
  text-transform: uppercase;
  letter-spacing: 0.08em;
  font-size: 0.72rem;
  color: #a8a29e;
  margin: 1rem 0 0.5rem 0.55rem;
}
.nav a {
  display: block;
  color: #e7e5e4 !important;
  text-decoration: none !important;
  padding: 0.55rem 0.75rem;
  border-radius: 0.65rem;
  margin-bottom: 0.2rem;
  font-size: 0.95rem;
}
.nav a:hover { background: rgba(255,255,255,0.06); color: #fff !important; }
.nav a.active {
  background: var(--accent);
  color: #fff !important;
  font-weight: 600;
}
.content-wrap {
  padding: 1.5rem clamp(1rem, 3vw, 2.5rem) 3rem;
}
.topbar {
  display: none;
  gap: 0.75rem;
  align-items: center;
  margin-bottom: 1rem;
}
.menu-btn {
  border: 1px solid var(--line);
  background: var(--surface);
  border-radius: 0.7rem;
  padding: 0.55rem 0.8rem;
  font: inherit;
  cursor: pointer;
}
.article {
  max-width: 920px;
  margin: 0 auto;
  background: var(--surface);
  border: 1px solid var(--line);
  border-radius: 1.1rem;
  box-shadow: var(--shadow);
  padding: clamp(1.25rem, 3vw, 2.4rem);
}
.article h1, .article h2, .article h3, .article h4 {
  line-height: 1.25;
  color: #0c0a09;
}
.article h1 { font-size: 2rem; margin-top: 0; }
.article h2 {
  margin-top: 2rem;
  padding-bottom: 0.35rem;
  border-bottom: 1px solid var(--line);
}
.article pre {
  background: var(--code-bg) !important;
  border: 1px solid var(--line);
  border-radius: 0.75rem;
  overflow-x: auto;
  margin: 0;
  padding: 0.9rem 1rem;
}
.code-block {
  position: relative;
  margin: 1.1rem 0;
}
.code-block pre {
  padding-right: 4.2rem;
}
.copy-btn {
  position: absolute;
  top: 0.45rem;
  right: 0.45rem;
  z-index: 2;
  display: inline-flex;
  align-items: center;
  gap: 0.3rem;
  border: 1px solid var(--line);
  background: rgba(255,255,255,0.92);
  color: var(--muted);
  border-radius: 0.5rem;
  padding: 0.28rem 0.55rem;
  font: inherit;
  font-size: 0.78rem;
  line-height: 1;
  cursor: pointer;
  box-shadow: 0 1px 2px rgba(0,0,0,0.04);
  transition: background .15s ease, color .15s ease, border-color .15s ease;
}
.copy-btn:hover {
  background: #fff;
  color: var(--ink);
  border-color: #a8a29e;
}
.copy-btn:focus-visible {
  outline: 2px solid var(--accent);
  outline-offset: 2px;
}
.copy-btn.copied {
  background: var(--accent-soft);
  color: var(--accent);
  border-color: #99f6e4;
}
.copy-btn svg {
  width: 14px;
  height: 14px;
  flex: 0 0 auto;
}
.article img {
  max-width: 100%;
  height: auto;
  border-radius: 0.6rem;
}
.article figure { margin: 1.4rem auto; text-align: center; }
.pager {
  display: flex;
  justify-content: space-between;
  gap: 1rem;
  margin-top: 2rem;
  padding-top: 1rem;
  border-top: 1px solid var(--line);
  flex-wrap: wrap;
}
.pager a {
  text-decoration: none !important;
  border: 1px solid var(--line);
  border-radius: 0.8rem;
  padding: 0.75rem 1rem;
  background: #fafaf9;
  color: var(--ink) !important;
  min-width: min(100%, 220px);
}
.pager a:hover { border-color: var(--accent); background: var(--accent-soft); }
.pager .label { display:block; color: var(--muted); font-size: 0.8rem; margin-bottom: 0.15rem; }
.card-grid {
  display: grid;
  grid-template-columns: repeat(auto-fit, minmax(220px, 1fr));
  gap: 0.9rem;
  margin: 1.5rem 0 0.5rem;
}
.card {
  display: block;
  text-decoration: none !important;
  color: inherit !important;
  border: 1px solid var(--line);
  border-radius: 0.9rem;
  padding: 1rem 1.05rem;
  background: linear-gradient(180deg, #fff, #fafaf9);
  transition: transform .15s ease, box-shadow .15s ease, border-color .15s ease;
}
.card:hover {
  transform: translateY(-2px);
  border-color: #99f6e4;
  box-shadow: var(--shadow);
}
.card .num {
  display: inline-flex;
  width: 1.7rem;
  height: 1.7rem;
  align-items: center;
  justify-content: center;
  border-radius: 999px;
  background: var(--accent-soft);
  color: var(--accent);
  font-weight: 700;
  font-size: 0.85rem;
  margin-bottom: 0.55rem;
}
.card strong { display: block; font-size: 1.02rem; margin-bottom: 0.2rem; }
.card span { color: var(--muted); font-size: 0.9rem; }
.links {
  display: flex;
  flex-wrap: wrap;
  gap: 0.6rem;
  justify-content: center;
  margin: 1.2rem 0 0.2rem;
}
.links a {
  text-decoration: none !important;
  background: var(--accent);
  color: white !important;
  padding: 0.55rem 0.9rem;
  border-radius: 999px;
  font-size: 0.92rem;
}
.backdrop {
  display: none;
  position: fixed;
  inset: 0;
  background: rgba(0,0,0,.35);
  z-index: 20;
}
/* soften typst default full-page styles inside article */
.article > style { display: none; }
@media (max-width: 920px) {
  .layout { grid-template-columns: 1fr; }
  .topbar { display: flex; }
  .sidebar {
    position: fixed;
    left: 0; top: 0; bottom: 0;
    width: min(86vw, 320px);
    z-index: 30;
    transform: translateX(-105%);
    transition: transform .2s ease;
  }
  .sidebar.open { transform: translateX(0); }
  .backdrop.show { display: block; }
}
`;

function findTypst() {
  try {
    execFileSync("typst", ["--version"], { stdio: "ignore" });
    return "typst";
  } catch {
    const win = "C:/Users/lucas.s.campos/AppData/Local/Microsoft/WinGet/Packages/Typst.Typst_Microsoft.Winget.Source_8wekyb3d8bbwe/typst-x86_64-pc-windows-msvc/typst.exe";
    if (existsSync(win)) return win;
    throw new Error("typst not found");
  }
}

function extractBodyAndAssets(html) {
  // Typst HTML: full document. Keep head styles/math fonts and body inner HTML.
  const titleMatch = html.match(/<title>(.*?)<\/title>/i);
  const title = titleMatch ? titleMatch[1] : "Notas";
  const styleMatches = [...html.matchAll(/<style[\s\S]*?<\/style>/gi)].map((m) => m[0]);
  const bodyMatch = html.match(/<body[^>]*>([\s\S]*)<\/body>/i);
  let body = bodyMatch ? bodyMatch[1] : html;
  // remove possible empty wrappers noise
  body = body.trim();
  return { title, styles: styleMatches.join("\n"), body };
}

function navHtml(currentId) {
  return PAGES.map((p) => {
    const cls = p.id === currentId ? ' class="active"' : "";
    return `<a href="${p.file}"${cls}>${p.title}</a>`;
  }).join("\n");
}

function pagerHtml(currentId) {
  const idx = PAGES.findIndex((p) => p.id === currentId);
  const prev = idx > 0 ? PAGES[idx - 1] : null;
  const next = idx >= 0 && idx < PAGES.length - 1 ? PAGES[idx + 1] : null;
  const left = prev
    ? `<a href="${prev.file}"><span class="label">Anterior</span><strong>${prev.title}</strong></a>`
    : "<span></span>";
  const right = next
    ? `<a href="${next.file}" style="text-align:right"><span class="label">Próximo</span><strong>${next.title}</strong></a>`
    : "<span></span>";
  return `${left}${right}`;
}

function enhanceIndexBody(body) {
  // Replace plain numbered list of chapter links with cards if present; else prepend cards.
  const cards = PAGES.filter((p) => p.id !== "index")
    .map((p, i) => {
      const blurb =
        p.id === "extra"
          ? "Materiais complementares e medidas de erro"
          : "Abrir capítulo";
      return `<a class="card" href="${p.file}"><div class="num">${i + 1}</div><strong>${p.title}</strong><span>${blurb}</span></a>`;
    })
    .join("\n");

  const heroExtras = `
<div class="links">
  <a href="https://1drv.ms/b/s!AmfyGvdmTYong45aJ5g2TBxKCkygcQ?e=DJQ9oC" target="_blank" rel="noopener">Apostila</a>
  <a href="https://forms.gle/7gKy3k1TqHaCUkVD8" target="_blank" rel="noopener">Entregas</a>
  <a href="https://www.youtube.com/playlist?list=PLajnQa6HBzEIrJXXrQfUAeYwdMk1ygAhr" target="_blank" rel="noopener">Playlist YouTube</a>
  <a href="https://github.com/l-s-campos/BEM" target="_blank" rel="noopener">Código BEM</a>
</div>
<div class="card-grid">
${cards}
</div>`;

  // Keep typst body, append card grid at top after first heading block
  return body + heroExtras;
}

function wrapPage(page, rawHtml) {
  const { title, styles, body } = extractBodyAndAssets(rawHtml);
  const content = page.id === "index" ? enhanceIndexBody(body) : body;
  const pageTitle =
    page.id === "index"
      ? "Introdução ao método dos elementos de contorno"
      : `${page.title} · Introdução ao método dos elementos de contorno`;

  return `<!DOCTYPE html>
<html lang="pt">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>${pageTitle}</title>
<meta name="authors" content="Lucas S. Campos">
${styles}
<style>${CSS}</style>
</head>
<body>
<div class="backdrop" id="backdrop"></div>
<div class="layout">
  <aside class="sidebar" id="sidebar">
    <a class="brand" href="index.html">Introdução ao método dos elementos de contorno</a>
    <div class="brand-sub">Notas de aula · Lucas S. Campos</div>
    <div class="nav-label">Capítulos</div>
    <nav class="nav">
      ${navHtml(page.id)}
    </nav>
  </aside>
  <div class="content-wrap">
    <div class="topbar">
      <button class="menu-btn" id="menuBtn" type="button">Menu</button>
      <strong>${page.title}</strong>
    </div>
    <main class="article">
      ${content}
      <div class="pager">
        ${pagerHtml(page.id)}
      </div>
    </main>
  </div>
</div>
<script>
const btn = document.getElementById('menuBtn');
const sidebar = document.getElementById('sidebar');
const backdrop = document.getElementById('backdrop');
function closeMenu(){ sidebar.classList.remove('open'); backdrop.classList.remove('show'); }
function openMenu(){ sidebar.classList.add('open'); backdrop.classList.add('show'); }
btn?.addEventListener('click', () => {
  if (sidebar.classList.contains('open')) closeMenu(); else openMenu();
});
backdrop?.addEventListener('click', closeMenu);

(function enhanceCodeBlocks() {
  const iconCopy = '<svg viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round" aria-hidden="true"><rect x="9" y="9" width="13" height="13" rx="2"></rect><path d="M5 15H4a2 2 0 0 1-2-2V4a2 2 0 0 1 2-2h9a2 2 0 0 1 2 2v1"></path></svg>';
  const iconCheck = '<svg viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round" aria-hidden="true"><path d="M20 6 9 17l-5-5"></path></svg>';

  async function copyText(text) {
    if (navigator.clipboard && window.isSecureContext) {
      await navigator.clipboard.writeText(text);
      return;
    }
    const ta = document.createElement('textarea');
    ta.value = text;
    ta.setAttribute('readonly', '');
    ta.style.position = 'fixed';
    ta.style.left = '-9999px';
    document.body.appendChild(ta);
    ta.select();
    document.execCommand('copy');
    document.body.removeChild(ta);
  }

  document.querySelectorAll('main.article pre').forEach((pre) => {
    if (pre.closest('.code-block')) return;
    const wrap = document.createElement('div');
    wrap.className = 'code-block';
    pre.parentNode.insertBefore(wrap, pre);
    wrap.appendChild(pre);

    const button = document.createElement('button');
    button.type = 'button';
    button.className = 'copy-btn';
    button.setAttribute('aria-label', 'Copiar c\u00f3digo');
    button.innerHTML = iconCopy + '<span>Copiar</span>';
    wrap.appendChild(button);

    let resetTimer;
    button.addEventListener('click', async () => {
      const code = pre.innerText.replace(/\\n$/, '');
      try {
        await copyText(code);
        button.classList.add('copied');
        button.innerHTML = iconCheck + '<span>Copiado</span>';
        clearTimeout(resetTimer);
        resetTimer = setTimeout(() => {
          button.classList.remove('copied');
          button.innerHTML = iconCopy + '<span>Copiar</span>';
        }, 1600);
      } catch (err) {
        button.innerHTML = iconCopy + '<span>Erro</span>';
        clearTimeout(resetTimer);
        resetTimer = setTimeout(() => {
          button.innerHTML = iconCopy + '<span>Copiar</span>';
        }, 1600);
      }
    });
  });
})();
</script>
</body>
</html>
`;
}

function compilePage(typst, page) {
  const out = join(RAW, page.file);
  console.log("compile", page.id);
  execFileSync(
    typst,
    ["compile", "--features", "html", "--format", "html", "--root", ROOT, page.src, out],
    { stdio: "inherit", cwd: ROOT }
  );
  return out;
}

function main() {
  const typst = findTypst();
  mkdirSync(DOCS, { recursive: true });
  mkdirSync(RAW, { recursive: true });

  // clean previous html pages (keep folder)
  for (const f of readdirSync(DOCS)) {
    if (f.endsWith(".html") || f === ".nojekyll") {
      // remove old single-page index etc later rewritten
    }
  }

  for (const page of PAGES) {
    const rawPath = compilePage(typst, page);
    const raw = readFileSync(rawPath, "utf8");
    const wrapped = wrapPage(page, raw);
    writeFileSync(join(DOCS, page.file), wrapped, "utf8");
  }

  writeFileSync(join(DOCS, ".nojekyll"), "", "utf8");
  // remove raw intermediates
  rmSync(RAW, { recursive: true, force: true });

  console.log("\nBuilt multipage site in docs/:");
  for (const p of PAGES) console.log(" -", p.file);
}

main();
