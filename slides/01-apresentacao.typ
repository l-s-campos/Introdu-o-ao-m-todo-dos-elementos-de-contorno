#import "_theme.typ": *
#show: bem-slides.with(
  title: [Apresentação],
  subtitle: [Por que BEM · curso 30 h],
)

= Por que BEM?

== Ideia em uma frase

#keybox[
  Discretizar só o *contorno* $Gamma = partial Omega$ — não o volume.
]

- Menos dimensão na malha
- Bom para exterior, fratura, contato, SF conhecida
- Preço: matrizes cheias e integrais singulares

== Por que menos usado que MEF?

+ Formulação matemática mais densa
+ Poucos códigos didáticos simples
+ Poucos cursos na graduação
+ Singularidades e CDC no contorno

*Este curso ataca sobretudo esses pontos.*

== Vantagens e limites

#grid(columns: (1fr, 1fr), gutter: 1.2em)[
  *Vantagens*
  - Malha só em $Gamma$
  - Precisão no contorno
  - Exterior natural
][
  *Limites*
  - $H, G$ densas
  - SF necessária
  - Não linearidade no domínio mais difícil
]

== Trilha

Notação → interpolação → ED → 2D → Laplace → Poisson → elasticidade → 3D → trabalhos

Código: `BEM_gmsh` (Julia + Gmsh)

== Laboratório

- Julia + DrWatson + Gmsh
- Smoke test: quadrado $T = x$
- Site, entregas e playlist no README
