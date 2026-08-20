// Tema compartilhado — apresentações do curso BEM
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
    config-common(
      slide-level: 2,
    ),
  )

  set text(size: 20pt, lang: "pt")
  set par(justify: false)

  title-slide()
  body
}

// Atalhos de layout
#let keybox(body) = block(
  width: 100%,
  fill: rgb("#ecfdf5"),
  inset: 14pt,
  radius: 6pt,
  stroke: 0.6pt + rgb("#99f6e4"),
  body,
)

#let muted(body) = text(fill: luma(80), size: 0.9em, body)
