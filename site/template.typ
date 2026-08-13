// Shared template for multipage HTML site

#let site-title = "Introdução ao método dos elementos de contorno"

#let pages = (
  (id: "index", title: "Início", file: "index.html"),
  (id: "apresentacao", title: "Apresentação", file: "apresentacao.html"),
  (id: "glossario", title: "Glossário", file: "glossario.html"),
  (id: "interpolacao", title: "Interpolação", file: "interpolacao.html"),
  (id: "equacoes-diferenciais", title: "Equações diferenciais", file: "equacoes-diferenciais.html"),
  (id: "indo-para-2d", title: "Indo para 2D", file: "indo-para-2d.html"),
  (id: "gmsh", title: "Gmsh e malhas", file: "gmsh.html"),
  (id: "laplace-2d", title: "Laplace 2D", file: "laplace-2d.html"),
  (id: "poisson-2d", title: "Poisson 2D", file: "poisson-2d.html"),
  (id: "elasticidade-2d", title: "Elasticidade 2D", file: "elasticidade-2d.html"),
  (id: "propgeo-3d", title: "Propgeo 3D", file: "propgeo-3d.html"),
  (id: "trabalhos-finais", title: "Trabalhos finais", file: "trabalhos-finais.html"),
  (id: "viga-euler", title: "Viga de Euler (extra)", file: "viga-euler.html"),
  (id: "extra", title: "Extra", file: "extra.html"),
)

#let site-css = `
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
  text-decoration: none;
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
  color: #e7e5e4;
  text-decoration: none;
  padding: 0.55rem 0.75rem;
  border-radius: 0.65rem;
  margin-bottom: 0.2rem;
  font-size: 0.95rem;
}
.nav a:hover { background: rgba(255,255,255,0.06); color: #fff; }
.nav a.active {
  background: var(--accent);
  color: #fff;
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
.article pre, .article code {
  font-family: ui-monospace, "Cascadia Code", "Consolas", monospace;
}
.article pre {
  background: var(--code-bg) !important;
  border: 1px solid var(--line);
  border-radius: 0.75rem;
  overflow-x: auto;
  padding: 0.9rem 1rem !important;
}
.article img {
  max-width: 100%;
  height: auto;
  border-radius: 0.6rem;
}
.article figure { margin: 1.4rem auto; text-align: center; }
.article figcaption { color: var(--muted); font-size: 0.92rem; margin-top: 0.4rem; }
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
  text-decoration: none;
  border: 1px solid var(--line);
  border-radius: 0.8rem;
  padding: 0.75rem 1rem;
  background: #fafaf9;
  color: var(--ink);
  min-width: min(100%, 220px);
}
.pager a:hover { border-color: var(--accent); background: var(--accent-soft); }
.pager .label { display:block; color: var(--muted); font-size: 0.8rem; margin-bottom: 0.15rem; }
.hero {
  text-align: center;
  padding: 1rem 0 0.5rem;
}
.hero h1 { font-size: clamp(1.8rem, 4vw, 2.5rem); margin-bottom: 0.4rem; }
.hero p { color: var(--muted); margin: 0.35rem 0; }
.card-grid {
  display: grid;
  grid-template-columns: repeat(auto-fit, minmax(220px, 1fr));
  gap: 0.9rem;
  margin: 1.5rem 0 0.5rem;
}
.card {
  display: block;
  text-decoration: none;
  color: inherit;
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
  text-decoration: none;
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
`

#let nav-html(current) = {
  let items = pages.map(p => {
    let cls = if p.id == current { " active" } else { "" }
    `<a class="nav-link` + cls + `" href="` + p.file + `">` + p.title + `</a>`
  }).join("\n")
  items
}

#let page-index(current) = {
  let i = 0
  for p in pages {
    if p.id == current { return i }
    i += 1
  }
  0
}

#let pager-html(current) = {
  let idx = page-index(current)
  let prev = if idx > 0 { pages.at(idx - 1) } else { none }
  let next = if idx + 1 < pages.len() { pages.at(idx + 1) } else { none }
  let left = if prev != none {
    `<a href="` + prev.file + `"><span class="label">Anterior</span><strong>` + prev.title + `</strong></a>`
  } else { `<span></span>` }
  let right = if next != none {
    `<a href="` + next.file + `" style="text-align:right"><span class="label">Próximo</span><strong>` + next.title + `</strong></a>`
  } else { `<span></span>` }
  left + right
}

#let site-page(current: "index", title: site-title, body) = {
  set document(title: title, author: "Lucas S. Campos")
  set text(lang: "pt", size: 11pt)
  set heading(numbering: none)
  set par(justify: true)
  set math.equation(numbering: none)

  show raw.where(block: true): it => block(
    width: 100%,
    fill: luma(245),
    inset: 10pt,
    radius: 4pt,
    stroke: 0.5pt + luma(220),
    it,
  )
  show link: underline

  // HTML shell via raw HTML export hooks when available; fallback is content only.
  // We inject chrome by wrapping content and using a post-build script for full shell.
  body
}

// For PDF-friendly single chapter compile
#let chapter-page(title, body) = {
  set document(title: title + " — " + site-title, author: "Lucas S. Campos")
  set text(lang: "pt", size: 11pt)
  set heading(numbering: "1.1")
  set par(justify: true)
  set math.equation(numbering: none)
  show raw.where(block: true): it => block(
    width: 100%,
    fill: luma(245),
    inset: 10pt,
    radius: 4pt,
    stroke: 0.5pt + luma(220),
    it,
  )
  show link: underline
  body
}
