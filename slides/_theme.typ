// Tema compartilhado — slides com conteúdo completo da aula
#import "@preview/touying:0.7.4": *
#import themes.university: *

#let bem-slides(
  title: [],
  subtitle: [Curso 30 h · BEM_gmsh],
  body,
) = {
  show: university-theme.with(
    aspect-ratio: "16-9",
    config-info(
      title: title,
      subtitle: subtitle,
      author: [Lucas S. Campos],
      date: none,
      institution: [Introdução ao MEC / BEM],
    ),
    config-common(slide-level: 2),
  )

  set text(size: 14pt, lang: "pt")
  set par(justify: true, leading: 0.55em)
  show raw.where(block: true): it => {
    set text(size: 10pt)
    block(
      width: 100%,
      fill: luma(245),
      inset: 7pt,
      radius: 3pt,
      stroke: 0.4pt + luma(220),
      it,
    )
  }
  show raw.where(block: false): set text(size: 0.9em)
  // imagens um pouco menores
  show image: it => {
    set align(center)
    it
  }

  title-slide()
  body
}

#let keybox(body) = block(
  width: 100%,
  fill: rgb("#ecfdf5"),
  inset: 9pt,
  radius: 4pt,
  stroke: 0.5pt + rgb("#99f6e4"),
  body,
)
