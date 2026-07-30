#import "../template.typ": chapter-page
#chapter-page("Extra")[
  = Extra
  <extra>

  == Materiais complementares

  - #link("https://1drv.ms/b/s!AmfyGvdmTYong45aJ5g2TBxKCkygcQ?e=DJQ9oC")[Apostila]
  - #link("https://forms.gle/7gKy3k1TqHaCUkVD8")[Entregas]
  - #link("https://www.youtube.com/playlist?list=PLajnQa6HBzEIrJXXrQfUAeYwdMk1ygAhr")[Playlist YouTube]
  - #link("https://github.com/l-s-campos/BEM")[Código BEM (GitHub)]

  == Medidas de erro

  erro médio

  $bar.v.double upright(bold(u)) bar.v.double_"med" = sum | u_i |$

  $epsilon_u = (bar.v.double upright(bold(u)) - upright(bold(u))^"exact" bar.v.double_"med") / (bar.v.double upright(bold(u))^"exact" bar.v.double_"med")$

  erro máximo

  $bar.v.double upright(bold(u)) bar.v.double_"max" = max | u_i |$

  $epsilon_u = (bar.v.double upright(bold(u)) - upright(bold(u))^"exact" bar.v.double_"max") / (bar.v.double upright(bold(u))^"exact" bar.v.double_"max")$

  norma $L_2$ do erro

  $bar.v.double upright(bold(u)) bar.v.double_(L_2) = sqrt(sum_(i = 1)^(N_e) (u_i)^2)$

  $epsilon_u = (bar.v.double upright(bold(u)) - upright(bold(u))^"exact" bar.v.double_(L_2)) / (bar.v.double upright(bold(u))^"exact" bar.v.double_(L_2))$
]
