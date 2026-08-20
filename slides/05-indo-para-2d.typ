#import "_theme.typ": *
#show: bem-slides.with(
  title: [Indo para 2D],
  subtitle: [Malha · N · J · n · geo],
)

= O salto

== 1D → 2D

- Contorno vira curva $Gamma$
- Elementos em $xi in [-1,1]$
- Precisamos de $N_k$, $J$, $upright(bold(n))$

CDC → capítulo Laplace 2D

== Pipeline

```text
Gmsh → format2d → dad → geometric_props / plot_geo
```

== Descontínuo

- Nós de campo no interior do elemento
- Evita conflito de CDC no canto

== Geo

$ J = |d upright(bold(x))/d xi| $

Radial ou divergência → $P$, $A$, centróide

#keybox[
  Geometria sólida *antes* de montar $H,G$.
]
