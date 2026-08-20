#import "_theme.typ": *
#show: bem-slides.with(
  title: [Equações diferenciais],
  subtitle: [PVI · MDF · BEM 1D],
)

= Dois mundos

== PVI e contorno

- PVI: marcha no tempo / passo a passo
- Contorno: sistema algébrico nas pontas (BEM 1D)

== MDF

$ D_(x x) upright(bold(T)) = upright(bold(b)) $

Nós no domínio · matriz esparsa

== BEM 1D

#keybox[
  $ H T = G Q $ só nos extremos.
]

SF 1D · no tempo: método das linhas + $M$ nos Gauss

== Takeaway

SF → integral no contorno → $H,G$ → CDC → solve

(o mesmo padrão sobe para 2D/3D)
