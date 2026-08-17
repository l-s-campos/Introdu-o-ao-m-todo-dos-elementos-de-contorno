// Introdução ao método dos elementos de contorno
// markdown2typst + tex2typst → Typst nativo (sem TeX)
// Curso 30 h — código: BEM_gmsh

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
  #text(14pt)[Notas de aula — curso 30 h]
  #v(0.4cm)
  #text(12pt)[Lucas S. Campos]
  #v(1cm)
  #text(11pt)[
    Materiais complementares: \
    #link("https://1drv.ms/b/s!AmfyGvdmTYong45aJ5g2TBxKCkygcQ?e=DJQ9oC")[Apostila] ·
    #link("https://forms.gle/7gKy3k1TqHaCUkVD8")[Entregas] ·
    #link("https://www.youtube.com/playlist?list=PLajnQa6HBzEIrJXXrQfUAeYwdMk1ygAhr")[Playlist YouTube] \
    Código: #link("https://github.com/l-s-campos/BEM_gmsh")[BEM_gmsh]
  ]
]

#pagebreak()
#outline(title: "Sumário", indent: auto, depth: 2)
#pagebreak()

// ---- Núcleo do curso (30 h) ----
#include "chapters/01-apresentacao.typ"
#include "chapters/02-glossario.typ"
#include "chapters/03-interpolacao.typ"
#include "chapters/04-equacoes-diferenciais.typ"
#include "chapters/05-indo-para-2d.typ"
#include "chapters/06-laplace-2d.typ"
#include "chapters/09-poisson-2d.typ"
#include "chapters/10-elasticidade-2d.typ"
#include "chapters/11-propgeo-3d.typ"
#include "chapters/12-trabalhos-finais.typ"

// ---- Apêndice (referência) ----
#include "chapters/08-erros.typ"

// ---- Extra (fora da trilha principal de 30 h) ----
#pagebreak()
#align(center)[
  #v(2cm)
  #text(20pt, weight: "bold")[Material extra]
  #v(0.5cm)
  #text(12pt)[Tópicos avançados / aprofundamento opcional]
]
#pagebreak()

#include "chapters/90-viga-euler.typ"
