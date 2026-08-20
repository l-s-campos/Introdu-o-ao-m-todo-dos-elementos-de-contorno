#import "_theme.typ": *
#show: bem-slides.with(
  title: [Poisson 2D],
  subtitle: [DIBEM · domínio],
)

= Fonte de domínio

== PDE

$ nabla^2 T = f $

BIE ganha $integral_Omega T^* f$

== DIBEM

#keybox[
  $ H T - G q = M f $
]

- RBF nos nós (contorno + internos)
- `DIBEM(dad)` monta a matriz `M`

== Estacionário

```julia
H_G_full_direct(dad, 16)
M = DIBEM(dad)
d = M * fvec
applyBC(dad)
dad.b .+= d
```

== Tempo

$M$ como massa: `solve_transient` / Houbolt
