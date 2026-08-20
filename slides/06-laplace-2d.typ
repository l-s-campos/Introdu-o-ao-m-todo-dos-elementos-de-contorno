#import "_theme.typ": *
#show: bem-slides.with(
  title: [Laplace 2D],
  subtitle: [Primeiro BEM 2D completo],
)

= Pipeline

== Sete passos

PDE → Green → SF → $H,G$ → CDC → solve → internos

== Código mínimo

```julia
dad = format2d(msh, Laplace(1.0))
H_G_full_direct(dad, 20)
solve(dad)
rel_error(dad)
```

== CDC

`"0;T"` Dirichlet · `"1;q"` Neumann

== Integral

#keybox[
  $ c T = integral T q^* - integral T^* q $
]

$c = 1/2$ no contorno liso

== Diagonal de $H$

$T equiv 1$, $q equiv 0$ → $H_(i i) = - sum_(j != i) H_(i j)$

== Extra

H-matriz só para $N$ grande
