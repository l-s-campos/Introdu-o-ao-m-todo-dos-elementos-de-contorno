# Apresentações (Touying)

Slides **com o conteúdo completo de cada aula**, gerados a partir de `chapters/*.typ`.
Pacote: [Touying](https://typst.app/universe/package/touying/) · tema `university` · 16:9.

## Compilar um

```bash
typst compile --root . slides/06-laplace-2d.typ docs/slides/06-laplace-2d.pdf
```

## Compilar todos

```bash
mkdir -p docs/slides
for f in slides/[0-9]*.typ; do
  b=$(basename "$f" .typ)
  typst compile --root . "$f" "docs/slides/${b}.pdf"
done
```

## Estrutura

| Item | Papel |
|------|--------|
| `_theme.typ` | tema comum (`@preview/touying:0.7.4`) |
| `0x-*.typ` / `9x-*.typ` | um deck por capítulo |
| `docs/slides/*.pdf` | PDFs compilados |

Cada `==` das notas vira um slide (seções longas usam fonte menor). Fórmulas, código, tabelas e imagens das notas entram no deck.

## Capítulos

| Fonte | PDF |
|-------|-----|
| `01-apresentacao.typ` | [pdf](../docs/slides/01-apresentacao.pdf) |
| `02-glossario.typ` | [pdf](../docs/slides/02-glossario.pdf) |
| `03-interpolacao.typ` | [pdf](../docs/slides/03-interpolacao.pdf) |
| `04-equacoes-diferenciais.typ` | [pdf](../docs/slides/04-equacoes-diferenciais.pdf) |
| `05-indo-para-2d.typ` | [pdf](../docs/slides/05-indo-para-2d.pdf) |
| `06-laplace-2d.typ` | [pdf](../docs/slides/06-laplace-2d.pdf) |
| `08-erros.typ` | [pdf](../docs/slides/08-erros.pdf) |
| `09-poisson-2d.typ` | [pdf](../docs/slides/09-poisson-2d.pdf) |
| `10-elasticidade-2d.typ` | [pdf](../docs/slides/10-elasticidade-2d.pdf) |
| `11-indo-para-3d.typ` | [pdf](../docs/slides/11-indo-para-3d.pdf) |
| `12-trabalhos-finais.typ` | [pdf](../docs/slides/12-trabalhos-finais.pdf) |
| `90-viga-euler.typ` | [pdf](../docs/slides/90-viga-euler.pdf) |
| `91-contato-halfspace.typ` | [pdf](../docs/slides/91-contato-halfspace.pdf) |
