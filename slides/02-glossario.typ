#import "_theme.typ": *
#show: bem-slides.with(
  title: [Glossário e notação],
  subtitle: [Símbolos do curso],
)

= Domínio

== Símbolos

- $Omega$ — domínio
- $Gamma = partial Omega$ — contorno (normal saindo)
- $upright(bold(n))$ — normal exterior
- $upright(bold(x))$ — campo · $upright(bold(x))_d$ — fonte
- $r = |upright(bold(x)) - upright(bold(x))_d|$

== Potencial

#keybox[
  $ q := - k (partial T) / (partial n) $
]

- $T$ — potencial / temperatura
- $T^*, q^*$ — solução fundamental
- $H T = G q$ · após CDC: $A x = b$

== CDC e Gmsh

- Dirichlet / Neumann
- Laplace: `"0;T"` e `"1;q"`
- Elasticidade: `tipo;valor` por direção

== Siglas

BEM/MEC · SF · CDC · DIBEM · RBF · H-matriz
