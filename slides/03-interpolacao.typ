#import "_theme.typ": *
#show: bem-slides.with(
  title: [Interpolação],
  subtitle: [Nós, formas, Gauss],
)

= No BEM

== Aproximação no elemento

$ T(xi) approx sum_k N_k (xi) T_k $

- $N_k$ — funções de forma
- $xi in [-1, 1]$
- Base de `shapefun` no pacote

== Lagrange

- $N_i (xi_j) = delta_(i j)$
- $sum_i N_i = 1$
- Grau $p$ → $p+1$ nós

== Gauss–Legendre

- Quadratura em $[-1,1]$
- Nós descontínuos no `BEM_gmsh`

== Mensagem

#keybox[
  Boa interpolação + boa quadratura → $H$ e $G$ confiáveis.
]
