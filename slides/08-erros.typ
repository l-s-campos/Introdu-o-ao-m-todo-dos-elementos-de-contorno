#import "_theme.typ": *
#show: bem-slides.with(
  title: [Medidas de erro],
  subtitle: [Apêndice],
)

= Como reportar

== Sempre diga

+ Onde mediu (contorno / internos)
+ Qual norma
+ Erro *relativo*
+ $N$ ou $h$

== Normas

- $epsilon_infinity$
- $epsilon_1$ (média)
- $epsilon_2$ (euclidiana)

#keybox[
  `rel_error(dad)` = $epsilon_2$ na primária do `dad`.
]

== Convergência

Um parâmetro por vez: `ndiv`, ordem, `npg`.

Tabela $N$ × erros + taxa aparente.
