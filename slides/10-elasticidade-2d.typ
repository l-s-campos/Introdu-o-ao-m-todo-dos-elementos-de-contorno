#import "_theme.typ": *
#show: bem-slides.with(
  title: [Elasticidade 2D],
  subtitle: [Kelvin · 2 DOFs],
)

= Mesmo pipeline

== Troca de papéis

- Campo: $T$ → $upright(bold(u))$
- Neumann: $q$ → $upright(bold(t))$
- DOFs/nó: 1 → 2
- SF: potencial → Kelvin

== Código

```julia
dad = format2d(msh, Elasticity(E, nu, rho))
H_G_full_direct(dad, 16)
solve(dad)
```

== CDC

`tipo;val;tipo;val` em $x$ e $y$

0 = deslocamento · 1 = tração

== Labs

Patch $upright(bold(u)) = bold(epsilon) upright(bold(x))$

Tubo / Kirsch / viga em `data/elastico/iso/`

== Plane stress

`Elasticity(...; plane_stress=true)`
