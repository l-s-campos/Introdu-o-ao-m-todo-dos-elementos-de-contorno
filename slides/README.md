# Apresentações (Touying)

Slides **com o conteúdo completo de cada aula**, gerados a partir de `chapters/*.typ`.
Pacote: [Touying](https://typst.app/universe/package/touying/) · tema `university` · 16:9.

## Compilar um

```bash
typst compile --root . slides/06-laplace-2d.typ slides/06-laplace-2d.pdf
```

## Compilar todos

```bash
for f in slides/[0-9]*.typ; do
  b=$(basename "$f" .typ)
  typst compile --root . "$f" "slides/${b}.pdf"
done
```

## Estrutura

| Item | Papel |
|------|--------|
| `_theme.typ` | tema comum (`@preview/touying:0.7.4`) |
| `0x-*.typ` / `9x-*.typ` | um deck por capítulo |
| `slides/*.pdf` | PDFs compilados (junto das fontes) |

Cada `==` das notas vira um slide (seções longas usam fonte menor). Fórmulas, código, tabelas e imagens das notas entram no deck.

## Capítulos

| Fonte | PDF |
|-------|-----|
| `01-apresentacao.typ` | [pdf](01-apresentacao.pdf) |
| `02-glossario.typ` | [pdf](02-glossario.pdf) |
| `03-interpolacao.typ` | [pdf](03-interpolacao.pdf) |
| `04-equacoes-diferenciais.typ` | [pdf](04-equacoes-diferenciais.pdf) |
| `05-indo-para-2d.typ` | [pdf](05-indo-para-2d.pdf) |
| `06-laplace-2d.typ` | [pdf](06-laplace-2d.pdf) |
| `08-erros.typ` | [pdf](08-erros.pdf) |
| `09-poisson-2d.typ` | [pdf](09-poisson-2d.pdf) |
| `10-elasticidade-2d.typ` | [pdf](10-elasticidade-2d.pdf) |
| `11-indo-para-3d.typ` | [pdf](11-indo-para-3d.pdf) |
| `12-trabalhos-finais.typ` | [pdf](12-trabalhos-finais.pdf) |
| `90-viga-euler.typ` | [pdf](90-viga-euler.pdf) |
| `91-contato-halfspace.typ` | [pdf](91-contato-halfspace.pdf) |
