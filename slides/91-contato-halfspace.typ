#import "_theme.typ": *
#show: bem-slides.with(
  title: [Contato half-space],
  subtitle: [Extra · Hertz · Pohrt–Li],
)

= Extra

== Não é contorno fechado

#keybox[
  $ u = K * p $ na grade do half-space / half-plane
]

$p >= 0$, gap $>= 0$ (zona ativa)

== 2D e 3D

- Linha: Flamant · `solve_line_contact`
- Esfera: Boussinesq · `solve_normal_contact`

== Labs

`hertz_line_2d.jl` · `contact_pohrt_li.jl`

== Trabalhos

Mapa mínimo da proposta E
