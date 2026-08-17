// Apêndice — Medidas de erro
// Referência para todos os labs do curso (BEM_gmsh)

= Apêndice: medidas de erro
<erros>

Este apêndice fixa *como* reportar erro nos exercícios e trabalhos.
Não é uma aula da trilha: use-o quando um capítulo pedir “erro médio / $L_2$ / `rel_error`”.

== Princípio

Sempre que houver referência $upright(bold(u))^"ex"$ nos *mesmos* pontos que o numérico $upright(bold(u))$:

1. diga *onde* mediu (contorno, internos, ambos);
2. diga *qual* norma;
3. prefira *erro relativo* (normalizado pela referência);
4. em estudos de convergência, reporte *pelo menos duas* normas (ex.: $L_2$ e $infinity$).

Seja $e_i = u_i - u_i^"ex"$ o erro nodal (componente a componente; em elasticidade empilhe $u_x,u_y$).

== Normas em pontos discretos

*Máximo ($infinity$):*

$
norm(e)_infinity = max_i |e_i| ,
quad
epsilon_infinity = norm(e)_infinity \/ max(norm(u^"ex")_infinity, epsilon_"floor") .
$

*Média ($L_1$ relativa):*

$
norm(e)_1 = sum_i |e_i| ,
quad
epsilon_1 = norm(e)_1 \/ max(sum_i |u_i^"ex"|, epsilon_"floor") .
$

*$L_2$ (euclidiana nos graus de liberdade amostrados):*

$
norm(e)_2 = (sum_i e_i^2)^(1\/2) ,
quad
epsilon_2 = norm(e)_2 \/ max(norm(u^"ex")_2, epsilon_"floor") .
$

Use $epsilon_"floor" > 0$ pequeno (ex. $10^(-16)$) só para evitar divisão por zero quando a referência é nula.

#block(
  width: 100%, fill: luma(248), inset: 10pt, radius: 4pt, stroke: 0.5pt + luma(200),
)[
  *Atenção tipográfica.* A “norma $L_2$ do *erro*” é $norm(u - u^"ex")_2$, *não*
  $norm(u)_2 - norm(u^"ex")_2$. Nas notas antigas, a barra $bar.v.double u - u^"ex" bar.v.double$
  significava aplicar a norma ao *vetor diferença*.
]

== O que o `BEM_gmsh` calcula

Com solução analítica no cache (`attach_analytical!` ou `apply_analytical_bc!`) e após `solve`:

```julia
err = rel_error(dad)   # ||T_num - T_ana||_2 / ||T_ana||_2  (mesmos DOFs)
@show err
```

Essência (`src/Core/Analytical.jl`):

```julia
function rel_error(dad::BEMdata; t=0.0)
    Tana = analytical(dad; t=t)
    Tnum = dad.T   # ou dad.u em elasticidade
    n = min(length(Tnum), length(Tana))
    num = norm(Tnum[1:n] .- Tana[1:n])
    den = norm(Tana[1:n])
    return den > 0 ? num / den : num
end
```

Ou seja: `rel_error` $equiv$ $epsilon_2$ no vetor de solução primária armazenado
(potencial $T$ ou deslocamentos empilhados). *Não* substitui erro em fluxo $q$,
tensões ou sensores internos — calcule à parte quando o exercício pedir.

== Receita mínima de convergência

```julia
rows = []
for ndiv in (8, 16, 32)
    msh = quadrado(ndiv=ndiv, show=false)
    dad = format2d(msh, Laplace(1.0))
    attach_analytical!(dad, ana_laplace_linear(; direction=SA[1.0, 0.0]))
    H_G_full_direct(dad, 16)
    solve(dad)
    push!(rows, (; ndiv, N=dad.n, e2=rel_error(dad)))
end
@show rows
# ordem aparente: log(e_coarse/e_fine) / log(h_coarse/h_fine)
```

Boas práticas:
- refine *um* parâmetro de cada vez (`ndiv`, `tipo`\/ordem, `npg`);
- anote $N$ (DOFs), não só `ndiv`;
- se a solução for polinomial e o elemento a reproduz, espere erro de máquina no patch test;
- em cantos com CDC mista, compare também pontos *longe* da singularidade geométrica.

== Onde medir

#table(
  columns: (auto, auto),
  inset: 7pt,
  stroke: 0.5pt + luma(200),
  [*Local*], [*Uso típico*],
  [Contorno (nós BEM)], [padrão do `rel_error`; CDCs e $H T = G q$],
  [Pontos internos], [pós-processamento; mapas; exercícios de Poisson],
  [Sensores \/ linha], [transiente, ondas, trabalhos finais],
  [Fluxo $q$ ou tração $t$], [sempre que a CDC de Neumann for o alvo],
)

== Checklist de relatório

+ Definição da norma e do conjunto de pontos.
+ Tabela $N$ (ou $h$) × erros.
+ Uma frase sobre a taxa observada vs. a esperada (ordem do elemento).
+ Se usou só `rel_error`, declare que é $epsilon_2$ na primária do `dad`.
