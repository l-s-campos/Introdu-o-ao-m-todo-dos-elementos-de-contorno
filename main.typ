// Introdução ao método dos elementos de contorno
// markdown2typst + tex2typst → Typst nativo (sem TeX)

#set document(
  title: "Introdução ao método dos elementos de contorno",
  author: "Lucas S. Campos",
)

#set text(lang: "pt", size: 11pt)
#set heading(numbering: "1.1")
#set par(justify: true)
#set math.equation(numbering: none)

#show raw.where(block: true): it => block(
  width: 100%,
  fill: luma(245),
  inset: 10pt,
  radius: 4pt,
  stroke: 0.5pt + luma(220),
  it,
)

#show link: underline

#align(center)[
  #v(2cm)
  #text(24pt, weight: "bold")[Introdução ao método dos elementos de contorno]
  #v(0.8cm)
  #text(14pt)[Notas de aula]
  #v(0.4cm)
  #text(12pt)[Lucas S. Campos]
  #v(1cm)
  #text(11pt)[
    Materiais complementares: \
    #link("https://1drv.ms/b/s!AmfyGvdmTYong45aJ5g2TBxKCkygcQ?e=DJQ9oC")[Apostila] ·
    #link("https://forms.gle/7gKy3k1TqHaCUkVD8")[Entregas] ·
    #link("https://www.youtube.com/playlist?list=PLajnQa6HBzEIrJXXrQfUAeYwdMk1ygAhr")[Playlist YouTube]
  ]
]

#pagebreak()
#outline(title: "Sumário", indent: auto, depth: 2)
#pagebreak()

#include "chapters/apresentacao.typ"
#include "chapters/interpolacao.typ"
#include "chapters/equacoes-diferenciais.typ"
#include "chapters/viga-euler.typ"
#include "chapters/indo-para-2d.typ"
#include "chapters/laplace-2d.typ"
#include "chapters/erros.typ"
#include "chapters/poisson-2d.typ"
#include "chapters/elasticidade-2d.typ"
#include "chapters/propgeo-3d.typ"
