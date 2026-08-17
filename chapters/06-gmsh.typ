// Gmsh e malhas — ponte (geometria já em "Indo para 2D")
// CDC detalhada: capítulo Laplace 2D

= Gmsh e malhas para o BEM
<gmsh>

A *geometria* com Gmsh + `BEM_gmsh` (`.geo`, orientação, `format2d`, `shapefun`,
`geometric_props`, elementos descontínuos) está no capítulo *Indo para 2D*.

Este capítulo só fixa o que *ainda não* entra lá: o papel dos grupos físicos como
*portadores de CDC* e o caminho até o solver. O dicionário completo
(`"0;T"`, `"1;q"`, elasticidade, `attach_analytical!`, `solve`) é da aula *Laplace 2D*.

== Dois papéis do Gmsh no curso

#table(
  columns: (auto, auto),
  inset: 7pt,
  stroke: 0.5pt + luma(200),
  [*Aula*], [*Uso do Gmsh / `BEM_gmsh`*],
  [*Indo para 2D*], [malha, $N_k$, $J$, $upright(bold(n))$, perímetro\/área\/centróide],
  [*Laplace 2D*], [CDC nos grupos físicos, montagem $H,G$, `solve`, erro],
)

== Lembrete de instalação

```julia
using DrWatson
@quickactivate :BEM
include(datadir("Laplace", "Laplace_dad.jl"))
```

GUI do sistema: #link("https://gmsh.info/#Download")[gmsh.info].

== Grupos físicos → CDC (adiantado)

A partir de Laplace 2D, os *nomes* das curvas físicas carregam tipo e valor:

#table(
  columns: (auto, auto),
  inset: 8pt,
  stroke: 0.5pt + luma(200),
  [*Nome*], [*Significado (Laplace)*],
  [`"0;T"`], [Dirichlet: potencial\/temperatura $T$],
  [`"1;q"`], [Neumann: fluxo $q = -k partial T \/ partial n$],
)

Exemplos: `"0;100"`, `"1;0"` (isolado). Elasticidade usa `tx;ux;ty;uy` (ver Laplace \/ elasticidade).

Até lá, nos `.geo` de geometria pode usar só nomes descritivos (`"contorno"`, `"Domain"`),
como na aula anterior.

== Fluxo completo (preview)

```julia
props = Laplace(1.0)
msh   = quadrado(ndiv=20, show=false)
dad   = format2d(msh, props)   # lê CDC se os grupos estiverem no formato do repo

# a partir de Laplace 2D:
# attach_analytical!(dad, ana_laplace_linear(; direction=SA[1.0, 0.0]))
# H_G_full_direct(dad, 20)
# solve(dad)
# @show rel_error(dad)

plot_geo(dad)
gp = geometric_props(dad)
@show gp.perimeter gp.area gp.centroid
```

Montagem hierárquica (`H_G_Hmat`), ONELAB e geradores avançados: *Laplace 2D* e docs do `BEM_gmsh`.

== Boas práticas (repete o essencial)

1. Contorno externo anti-horário; furos horários.
2. Elementos descontínuos no campo (`format2d`) — cantos sem nó compartilhado de CDC.
3. `ordem` da malha alinhada a `tipo`.
4. Estudo de convergência: variar `ndiv` e medir erro (capítulo de erros + Laplace).

== Exercícios

1. Releia o capítulo *Indo para 2D* e rode `data/examples/geo_unit_square.jl`.
2. Abra um `.geo` de `data/Laplace/` e identifique o que é *geometria* e o que já é *CDC* (strings `"0;..."` \/ `"1;..."`).
3. (Opcional) Prepare um `.geo` só com geometria e outro com CDC no formato do repo, para usar na próxima aula.

== Leitura

- Capítulo *Indo para 2D*
- `BEM_gmsh/docs` (getting started, ONELAB)
- Manual Gmsh: #link("https://gmsh.info/doc/texinfo/gmsh.html")[documentação]
