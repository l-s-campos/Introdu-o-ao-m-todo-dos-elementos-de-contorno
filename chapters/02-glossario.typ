// Glossário e notação do curso
// Curso 30 h — BEM_gmsh

= Glossário e notação
<glossario>

Este capítulo fixa a *notação padrão* das notas. Quando um capítulo legado usar outro símbolo, a tabela abaixo prevalece.

== Domínio e contorno

#table(
  columns: (auto, auto),
  inset: 7pt,
  stroke: 0.5pt + luma(200),
  [*Símbolo*], [*Significado*],
  [$Omega$], [domínio (área em 2D, volume em 3D)],
  [$Gamma = partial Omega$], [contorno (orientado; normal saindo de $Omega$)],
  [$upright(bold(n)) = (n_x, n_y)$], [normal unitária exterior],
  [$upright(bold(x)) = (x, y)$], [ponto de campo (integração)],
  [$upright(bold(x))_d$ ou $p$], [ponto fonte (colocação)],
  [$r = | upright(bold(x)) - upright(bold(x))_d |$], [distância fonte–campo],
  [$N_k (xi)$], [função de forma no elemento de contorno],
  [$J$], [jacobiano da transformação para $xi in [-1, 1]$],
)

== Problema de potencial (Laplace / Poisson)

Usamos *temperatura/potencial* $T$ e *fluxo de contorno* $q$ com a convenção do `BEM_gmsh`:

$ q := - k (partial T)/(partial n) . $

#table(
  columns: (auto, auto),
  inset: 7pt,
  stroke: 0.5pt + luma(200),
  [*Símbolo*], [*Significado*],
  [$T$], [potencial / temperatura],
  [$q$], [fluxo no contorno, $q = -k partial T \/ partial n$],
  [$k$], [condutividade (Laplace: `Laplace(k)`)],
  [$f$], [fonte de domínio em Poisson: $nabla^2 T = f$],
  [$T^* , q^*$], [solução fundamental e sua derivada normal],
  [$H , G$], [matrizes de influência: $H T = G q$ (+ termo de domínio)],
  [$M$], [matriz de domínio (DIBEM / “massa”)],
  [$A , upright(bold(x)) , upright(bold(b))$], [sistema final após CDCs: $A upright(bold(x)) = upright(bold(b))$],
)

Em textos de mecânica dos fluidos o potencial costuma ser $phi$ ou $u$; *nestas notas de potencial térmico* preferimos $T$. Nos exercícios de membrana, a deflexão é $w$ e a equação é a de Poisson em $w$.

== Elasticidade 2D

#table(
  columns: (auto, auto),
  inset: 7pt,
  stroke: 0.5pt + luma(200),
  [*Símbolo*], [*Significado*],
  [$u_i$ ou $upright(bold(u)) = (u_x, u_y)$], [deslocamento],
  [$t_i$ ou $upright(bold(t))$], [tração no contorno, $t_i = sigma_(i j) n_j$],
  [$sigma_(i j) , epsilon_(i j)$], [tensão e deformação],
  [$E , nu$], [módulo de Young e coeficiente de Poisson],
  [$mu = E \/ (2(1 + nu))$], [módulo de cisalhamento],
  [$U_(i j) , T_(i j)$], [SF de deslocamento e de tração (Kelvin)],
  [$H , G$], [blocos $2 times 2$ por nó: análogo vetorial de $H T = G q$],
)

*Não* use $nu$ (nu) e $v$ (velocidade ou função peso) no mesmo parágrafo sem deixar claro. A função peso dos resíduos é $v$ ou $T^*$; o Poisson do material é sempre $nu$.

== Condições de contorno (CDC)

#table(
  columns: (auto, auto),
  inset: 7pt,
  stroke: 0.5pt + luma(200),
  [*Tipo*], [*O que é preescrito*],
  [Dirichlet], [$T$ (potencial) ou $u_i$ (deslocamento)],
  [Neumann], [$q$ (fluxo) ou $t_i$ (tração)],
  [Mista], [Dirichlet em parte de $Gamma$, Neumann no resto],
  [Robin / convecção], [$q = h (T - T_infinity)$ — combinação linear],
)

No Gmsh / `BEM_gmsh`: Laplace `"0;T"` ou `"1;q"`; elasticidade `"tx;ux;ty;uy"` (ver capítulo *Laplace 2D*).

== Integração e erro

#table(
  columns: (auto, auto),
  inset: 7pt,
  stroke: 0.5pt + luma(200),
  [*Símbolo*], [*Significado*],
  [$n$ ou $N$], [número de nós / pontos de colocação no contorno],
  [$n_"elem"$], [número de elementos de contorno],
  [$n_"pg"$ / `npg`], [pontos de Gauss por elemento],
  [$epsilon_u$], [erro relativo (médio, máximo ou $L_2$ — capítulo *Medidas de erro*)],
  [`rel_error(dad)`], [erro relativo agregado no `BEM_gmsh`],
)

== Siglas

#table(
  columns: (auto, auto),
  inset: 7pt,
  stroke: 0.5pt + luma(200),
  [*Sigla*], [*Significado*],
  [BEM / MEC], [Boundary Element Method / Método dos Elementos de Contorno],
  [MEF / FEM], [Método dos Elementos Finitos],
  [SF], [solução fundamental],
  [CDC], [condição de contorno],
  [DIBEM], [Dual Integration Boundary Element Method — termo de domínio por RBF],
  [RBF], [função de base radial],
  [H-matriz], [matriz hierárquica (compressão de blocos de $H$ e $G$)],
  [PVI], [problema de valor inicial],
)

== Gráficos no curso

Salvo menção em contrário, os scripts usam *Plots.jl*:

```julia
using Plots
plot(x, y; xlabel="x", ylabel="T", label="numérico", lw=2)
```

Funções de visualização do repositório de código (quando usadas) podem ter backend próprio; nestas notas o padrão didático é `Plots`.
