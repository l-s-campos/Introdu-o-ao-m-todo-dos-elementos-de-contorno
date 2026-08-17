// Laplace 2D — CDC + H,G + solve (BEM_gmsh)
// Geometria de contorno: capítulo "Indo para 2D"

= Laplace 2D
<laplace-2d>

Geometria, $N_k$, $J$ e $upright(bold(n))$ já vieram de *Indo para 2D*.
Aqui o Gmsh deixa de ser só malha e passa a carregar *CDC*; montamos $H,G$ e resolvemos o primeiro BEM 2D completo.

== Objetivos

+ Rodar o fluxo mínimo do quadrado ($T = x$) no `BEM_gmsh`.
+ Ler CDCs nos *grupos físicos* (`"0;T"`, `"1;q"`) e ligar a `dad.BC` \/ `dad.BV`.
+ Seguir o pipeline PDE $arrow.r$ equação integral $arrow.r$ $H T = G q$ $arrow.r$ $A x = b$ $arrow.r$ internos.
+ Entender a diagonal de $H$ (corpo rígido) e o papel de `npg` na montagem.

== Setup

Código: #link("https://github.com/l-s-campos/BEM_gmsh")[`BEM_gmsh`].
Na pasta do projeto:

```julia
using Pkg
Pkg.activate(".")
Pkg.instantiate()

using DrWatson
@quickactivate :BEM
include(datadir("Laplace", "Laplace_dad.jl"))
```

Se o capítulo *Indo para 2D* ainda não rodou, comece por `data/examples/geo_unit_square.jl`.

== Fluxo mínimo (quadrado, $T = x$)

Rode *antes* da álgebra — o resto do capítulo só nomeia cada linha.

```julia
props = Laplace(1.0)
msh   = quadrado(ndiv=20, show=false)
dad   = format2d(msh, props)

attach_analytical!(dad, ana_laplace_linear(; direction=SA[1.0, 0.0]))
H_G_full_direct(dad, 20)   # npg = 20 (Gauss na montagem)
solve(dad)

@show rel_error(dad)
plot_geo(dad)
```

As CDCs já vêm dos grupos físicos do `quadrado` (ver abaixo). `attach_analytical!` *não* sobrescreve CDC: só guarda a solução exata para `rel_error`.

== Por que Laplace \/ Poisson importam

A *mesma* PDE descreve vários problemas: mudam o nome do campo, a lei constitutiva e o significado de $f$ e do fluxo. Com $f = 0$ é Laplace; com fonte de domínio, Poisson (capítulo seguinte).

#figure(
  kind: table,
  caption: [Exemplos de problemas físicos descritos por Laplace ou Poisson.],
)[
#set text(size: 9.5pt)
#table(
  columns: (1.15fr, 1.15fr, 1.35fr, 1.35fr),
  inset: 6pt,
  stroke: 0.5pt + luma(200),
  align: (left, left, left, left),
  table.header(
    [*Equação*],
    [*Problema físico*],
    [*Grandezas*],
    [*Lei \/ fluxo*],
  ),
  [$nabla^2 phi = 0$],
  [Torção de Saint-Venant],
  [$phi$: função de empenamento],
  [Hooke (tensões a partir de $phi$)],

  [$nabla · (S nabla u) + f = 0$],
  [Deflexão de membranas],
  [$u$: deflexão; $f$: carga; $S$: tensão],
  [$bold(D) = S bold(I)$],

  [$nabla · (bold(D) nabla T) + f = 0$],
  [Condução de calor],
  [$T$: temperatura; $bold(D)$: condutividade; $f$: geração],
  [Fourier: $bold(q) = - bold(D) nabla T$],

  [$nabla · (nabla phi) + f = 0$],
  [Escoamento potencial (irrot., invíscido)],
  [$phi$: potencial de velocidade; $f$: fonte],
  [$bold(v) = nabla phi$],

  [$nabla · (bold(D) nabla phi) + f = 0$],
  [Meio poroso],
  [$phi$: carga piezométrica; $bold(D)$: permeabilidade],
  [Darcy: $bold(q) = - bold(D) nabla phi$],

  [$nabla · (bold(D) nabla V) + f = 0$],
  [Potencial elétrico],
  [$V$: potencial; $bold(D)$: condutividade elétrica],
  [Ohm: $bold(q) = - bold(D) nabla V$],
)
]

#block(
  width: 100%,
  fill: luma(248),
  inset: 10pt,
  radius: 4pt,
  stroke: 0.5pt + luma(200),
)[
  *No curso.* O protótipo do `BEM_gmsh` é potencial $T$ com
  $q = -k partial T \/ partial n$ (glossário). Anisotropia ou $f != 0$ mudam a SF
  ou pedem termo de domínio (Poisson \/ DIBEM). As matrizes $H,G$ seguem a *mesma*
  lógica: campo no contorno + fluxo conjugado.
]

== Grupos físicos $arrow.r$ CDC

No Gmsh, o *nome* da curva física carrega tipo e valor. O `format2d` lê isso e preenche `dad.BC` (0 = Dirichlet, 1 = Neumann) e `dad.BV` (valor).

#table(
  columns: (auto, auto),
  inset: 8pt,
  stroke: 0.5pt + luma(200),
  [*Nome*], [*Significado (Laplace)*],
  [`"0;T"`], [Dirichlet: $T$ prescrito],
  [`"1;q"`], [Neumann: $q = -k partial T \/ partial n$ prescrito],
)

Exemplos: `"0;100"`, `"1;0"` (isolado). Elasticidade: `"tx;ux;ty;uy"` (capítulo próprio).

=== O `quadrado` do pacote (CDCs de $T = x$)

Em `data/Laplace/Laplace_dad.jl`, com $k = 1$ e solução $T = x$
($q = - partial T \/ partial n$):

```julia
# Laplace_dad.jl — trecho de quadrado(...)
# l1 bottom, l2 right, l3 top, l4 left
gmsh.model.addPhysicalGroup(1, [l1, l3], -1, "1;0")   # isolado (q = 0)
gmsh.model.addPhysicalGroup(1, [l2], -1, "1;-1")      # direita: q = -1
gmsh.model.addPhysicalGroup(1, [l4], -1, "0;0")       # esquerda: T = 0
gmsh.model.addPhysicalGroup(2, [s1], -1, "Domain")
```

Conferência rápida: normal saindo à direita é $+upright(bold(e))_x$, $partial T\/partial n = 1$, logo $q = -1$; topo\/base $partial T\/partial n = 0$; esquerda $T = 0$.

=== Cantos e elemento descontínuo

Em um vértice, duas arestas podem pedir CDCs *incompatíveis* no mesmo nó
(ex.: $T$ de um lado e $q$ do outro, ou dois $T$ diferentes).
Por isso o `format2d` coloca nós de *campo* no interior do elemento
(`discontinuous_nodes_weights` \/ Gauss–Legendre) — cada elemento tem seus $T_i,q_i$,
sem nó compartilhado no canto. A geometria pode continuar isoparamétrica nos vértices
(`Equispaced`), como no capítulo *Indo para 2D*.

=== `attach_analytical!` vs CDC da malha

```julia
# Analytical.jl (essência)
function attach_analytical!(dad::BEMdata, ana::AnalyticalSolution)
    set_cache!(dad; analytical=ana)
    return dad
end
```

Só registra a referência para erro \/ gráficos.
Para *forçar* Dirichlet analítico em todos os nós (patch test), o pacote tem
`apply_analytical_bc!(dad, ana)` — isso *sim* sobrescreve `BC` \/ `BV`.
No fluxo didático padrão do quadrado, *não* chame `apply_analytical_bc!`:
as CDCs do `.msh` já batem com $T = x$.

== Mapa do BEM (leia antes da álgebra)

Cada seção seguinte preenche *um* passo. A diagonal de $H$ é detalhe do passo 4.

#block(
  width: 100%,
  fill: luma(248),
  inset: 12pt,
  radius: 6pt,
  stroke: 0.6pt + luma(200),
)[
  #set par(justify: false)
  #set text(size: 10.5pt)
  *1. PDE* $quad nabla^2 T = 0$ em $Omega$, CDCs em $Gamma$

  $arrow.b.double$

  *2. Resíduos + Green* $quad$ o operador passa de $T$ para o peso $v$

  $arrow.b.double$

  *3. Solução fundamental* $quad$ $v = T^*$ com $-nabla^2 T^* = delta(x - x_d)$;
  o domínio some $arrow.r$ equação integral só no contorno

  $arrow.b.double$

  *4. Discretização* $quad$ $T,q approx N_k$; integrais $arrow.r$ $H,G$
  (quadratura `npg`, singularidades, *diagonal de* $H$)

  $arrow.b.double$

  *5. CDC* $quad$ em cada nó, $T$ *ou* $q$ conhecido $arrow.r$ $H T = G q$ vira $A x = b$

  $arrow.b.double$

  *6. Solve* $quad x = A^(-1) b$ $arrow.r$ $(T,q)$ completo no contorno

  $arrow.b.double$

  *7. Internos* $quad$ com contorno conhecido, $T(x_d)$ interior é *só* integral
  — sem novo sistema linear
]

#table(
  columns: (auto, auto),
  inset: 7pt,
  stroke: 0.5pt + luma(200),
  [*Passo*], [*No `BEM_gmsh`*],
  [1–3 formulação], [`fundamental` \/ núcleos na montagem],
  [4 $H,G$], [`H_G_full_direct(dad, npg)`],
  [5 CDC], [grupos Gmsh + `format2d`; troca de colunas em `applyBC`],
  [6 solve], [`solve(dad)` $arrow.r$ `applyBC` + `bem_linsolve`],
  [7 internos], [`pontointerno=true` em `format2d`; `plot_geo`],
)

#block(
  width: 100%,
  fill: rgb("#ecfdf5"),
  inset: 10pt,
  radius: 4pt,
  stroke: 0.5pt + rgb("#99f6e4"),
)[
  *Ideia-chave.* No MEF montamos “rigidez” no *volume*. No BEM montamos $H$ e $G$ no *contorno*;
  o preço é matrizes cheias e integrais singulares quando a fonte $x_d$ cai no elemento integrado.
]

== Formulação

=== Passos 1–3: da PDE à equação integral

Problema canônico do capítulo:

$ nabla^2 T = 0 quad "em" Omega, $

com $q = -k partial T \/ partial n$ no contorno (quando $f != 0$, ver *Poisson 2D*).

*Em cinco movimentos:*

1. Resíduo ponderado: $integral_Omega v nabla^2 T dif Omega = 0$ para pesos $v$ adequados.
2. Integração por partes \/ 2ª identidade de Green ($u = T$, peso $v$):

$ integral_Omega (v nabla^2 T - T nabla^2 v) dif Omega
  = integral_Gamma (v partial_n T - T partial_n v) dif s . $

3. Escolhe-se a *solução fundamental* $v = T^*$ tal que
   $-nabla^2 T^* = delta(x - x_d)$ (em sentido distribucional).
4. O termo de domínio com $nabla^2 T^*$ vira avaliação em $x_d$; com $nabla^2 T = 0$,
   sobra uma *equação integral só em* $Gamma$.
5. No contorno liso o coeficiente de $T(x_d)$ é $1\/2$ (ângulo sólido); no interior é $1$;
   em cantos, $c(x_d) != 1\/2$ — por isso o código prefere elementos descontínuos ou
   calcula a diagonal de $H$ de forma *indireta* (passo 4).

Equação integral no contorno liso:

$ c(x_d) T(x_d)
  = integral_Gamma T q^* dif s - integral_Gamma T^* q dif s ,
  quad c = 1\/2 " (liso)" . $

No pacote, os núcleos 2D ($U hat(=) T^*$, $T hat(=) q^*$ no código) são:

```julia
# Laplace/Fundamental.jl
function fundamental(props::Laplace, r::SVector{2}, n::SVector{2})
    k = props.k
    R = _R(r)                                    # |r|
    G = -log(R) / (2π * k)                       # T*
    H = dot(r, n) / (R^2 * 2π)                   # q*
    return KernelPair(G, H)
end
```

Ou seja (notação das notas):

$ T^* = - (ln r)\/(2 pi k) ,
  quad
  q^* = (upright(bold(r)) · upright(bold(n))) \/ (2 pi r^2) . $

=== Passo 4: discretização $arrow.r$ $H$ e $G$

Em cada elemento descontínuo com $m$ nós:

$ T(xi) = sum_(k=1)^m N_k (xi) T_k ,
  quad
  q(xi) = sum_(k=1)^m N_k (xi) q_k . $

Colocando a fonte em cada nó $i$ e somando elementos $j = 1,...,n_"elem"$:

$ sum_j sum_k
  ( integral_(Gamma_j) q^* N_k dif Gamma ) T_k^(j)
  -
  sum_j sum_k
  ( integral_(Gamma_j) T^* N_k dif Gamma ) q_k^(j)
  = 0
  quad "(+ termo livre na diagonal de" H")" . $

Em forma matricial global:

$ H T = G q . $

Definições operacionais (linha $i$ = fonte no nó $i$; coluna ligada ao GL $k$):

$
H_(i k) = integral_Gamma q^* (x, x_i) N_k (x) dif Gamma
  + c_i delta_(i k) ,
quad
G_(i k) = integral_Gamma T^* (x, x_i) N_k (x) dif Gamma .
$

- $H$: núcleo *mais* singular ($q^* ~ 1\/r$ em 2D no elemento da fonte).
- $G$: singularidade *fraca* ($ln r$), em geral integrável com quadratura adequada.
- `npg` em `H_G_full_direct(dad, npg)` fixa a quadratura de montagem (eco do Gauss do cap. 05).

Trecho real da montagem densa (`src/Core/Assembly_full.jl`):

```julia
function H_G_full_direct(dad::BEMdata{<:Scalar}; npg=20, threaded=true)
    _init_quadrature!(dad, npg)
    H = zeros(dad.nt, dad.nt)
    G = zeros(dad.nt, dad.n)
    set_cache!(dad; H, G)
    # ... para cada fonte i e elemento el:
    #     integrate_element(...)  ou  _far_nodal_scalar!(...)
    #     H[i, j] += ...;  G[i, j] += ...
    @inbounds for i in 1:dad.nt
        H[i, i] = 0.0
        H[i, i] = -sum(H[i, :])     # diagonal por soma nula (corpo rígido)
    end
    return H, G
end
```

No elemento singular o pacote usa Guiggiani; no quase-singular, regra $sinh$; longe, lumping nodal — detalhes em `integrate_element`. Para a aula, o essencial é: *integra o que for regular e fecha a diagonal de $H$ por identidade*.

==== Diagonal de $H$: corpo a temperatura constante

Quando a fonte está no mesmo elemento, integrar $H_(i i)$ “na marra” é delicado.
Em vez disso, usa-se a solução exata trivial de Laplace

$ T(x) equiv 1 , quad q(x) equiv 0 . $

Ela implica $H {1} = G {0} = {0}$, logo *cada linha de* $H$ *soma zero*.
Com os $H_(i j)$ ($i != j$) já calculados,

$ H_(i i) = - sum_(j != i) H_(i j) . $

Isso *já inclui* o fator de ângulo sólido ($1\/2$ no liso, outro valor em cantos).
No código: exatamente o laço `H[i,i] = -sum(H[i,:])` acima.

#block(
  width: 100%,
  fill: luma(248),
  inset: 10pt,
  radius: 4pt,
  stroke: 0.5pt + luma(200),
)[
  *Ordem prática (passo 4):*
  1. zere $H$ e $G$;
  2. para cada par (fonte $i$, elemento $j$), some contribuições regulares \/ tratadas;
  3. $H_(i i) = -sum_(j != i) H_(i j)$;
  4. diagonal de $G$: em geral *integre* (singularidade fraca) — o truque indireto raramente é necessário.
]

=== Passos 5–6: CDC e $A x = b$

Depois de $H T = G q$ ainda falta, em cada nó, *uma* informação.
A CDC decide qual coluna vai para a esquerda (incógnita) e o que alimenta $b$.

Regra por nó $i$:
- *Dirichlet* ($T_i$ conhecido, `BC[i]==0`): incógnita $q_i$ $arrow.r$ coluna $i$ de $-G$ em $A$; $H_(:i) T_i$ vai a $b$;
- *Neumann* ($q_i$ conhecido, `BC[i]==1`): incógnita $T_i$ $arrow.r$ coluna $i$ de $H$ em $A$; $G_(:i) q_i$ vai a $b$.

No pacote (`src/Core/Boundary_conditions.jl`):

```julia
function applyBC(dad::BEMdata{<:Union{Laplace,OrthotropicLaplace}}, A, B, b)
    n = dad.n
    for bc in 1:n
        if dad.BC[bc] == 0          # Dirichlet: incógnita = q
            A[:, bc] .= .-view(dad.G, :, bc)
            B[:, bc] .= .-view(dad.H, :, bc)
        end
    end
    fill!(b, 0)
    for j in 1:n
        if dad.BC[j] == 0
            b .-= view(dad.H, :, j) .* dad.BV[j]
        else
            b .+= view(dad.G, :, j) .* dad.BV[j]
        end
    end
end
```

e o solve estacionário:

```julia
# Core/Solver.jl (essência)
function solve(dad::BEMdata{<:Union{Laplace,OrthotropicLaplace}}; ...)
    applyBC(dad; ...)
    x = bem_linsolve(dad.A, dad.b)
    # split_sol! → dad.T, dad.q
    return dad.T
end
```

==== Exemplo didático (1 elemento por lado)

Condução unidimensional com um elemento por aresta do retângulo:

#image("../assets/laplace-2d/exemplo-unidirecional.png", width: 80%)

Esquema (macron = conhecido):

$ H mat(macron(T_1); T_2; macron(T_3); T_4)
  = G mat(q_1; macron(q_2); q_3; macron(q_4)) . $

Reorganizando colunas conhecidas \/ desconhecidas obtém-se $A x = b$ com
$x = (q_1, T_2, q_3, T_4)$. No `quadrado` real as CDCs misturam Dirichlet e Neumann
como na seção dos grupos físicos — o mecanismo é o mesmo `applyBC`.

=== Passo 7: pontos internos

Para $x_d in Omega$ (fora de $Gamma$), $c = 1$:

$ T(x_d) = integral_Gamma T q^* dif s - integral_Gamma T^* q dif s . $

Com $(T,q)$ *já conhecidos em todo o contorno*, isso é pós-processamento:
*não* se monta nem se resolve outro $A x = b$.
Ative com `format2d(..., pontointerno=true)`.

== Lab guiado: quadrado $T = x$

Um único script (relatório em tela). CDCs = grupos do `quadrado`; referência = `ana_laplace_linear`.

```julia
using DrWatson
@quickactivate :BEM
using LinearAlgebra, Printf
include(datadir("Laplace", "Laplace_dad.jl"))

props = Laplace(1.0)
msh   = quadrado(ndiv=16, show=false)
dad   = format2d(msh, props; pontointerno=true)

println("N = ", dad.n, "  internos = ", length(dad.internalNodes))

attach_analytical!(dad, ana_laplace_linear(; direction=SA[1.0, 0.0]))
H_G_full_direct(dad, 16)
solve(dad)

@printf "rel_error = %.6e\n" rel_error(dad)

println(" i |    x |    y |   T_num | T_exato |    q_num")
for i in 1:min(8, dad.n)
    x, y = dad.Nodes[i]
    @printf "%2d | %5.3f | %5.3f | %7.5f | %7.5f | %8.5f\n" i x y dad.T[i] x dad.q[i]
end

eT = maximum(abs(dad.T[i] - dad.Nodes[i][1]) for i in 1:dad.n)
@printf "max |T_num - x| = %.6e\n" eT
plot_geo(dad)
```

Observar:
- `rel_error` cai com `ndiv in {8,16,32}` (tabela no exercício 0);
- lados verticais: $T approx x$; horizontais: $q approx 0$;
- internos seguem $T approx x$ sem novo sistema.

== Exercícios

Notação: $T$, $q = -k partial T \/ partial n$. Erros: *Apêndice: medidas de erro* e\/ou `rel_error(dad)`.
Ordem do elemento: `tipo` \/ `ordem` no gerador e em `format2d` (linear $p=1$, quadrático $p=2$, …).

*Escada*

+ *E0 — convergência no quadrado.* Para `ndiv in {8,16,32}` (e, se puder, `tipo in {1,2}`),
  rode o lab $T = x$, monte a tabela `ndiv, tipo, N, rel_error, max|T-x|` e comente a taxa aparente.

+ *E1 — leia a malha.* Abra o trecho de grupos físicos de `quadrado` (ou outro `.geo` em `data/Laplace/`)
  e classifique cada aresta como Dirichlet ou Neumann. Confira `dad.BC` \/ `dad.BV` após `format2d`.

+ *E2 — setor circular.* Elementos lineares, quadráticos e cúbicos; erro médio, máximo e $L_2$ no contorno
  para várias malhas; gráfico de convergência. Analítico:

$ T(theta) = theta\/pi ,
  quad
  q(x) = -1\/(pi x) quad "em" quad y = 0 . $

#image("../assets/laplace-2d/setor-circular.png", width: 80%)

+ *E3 — placa mista.*

- $T(x = 0) = 0$
- $T(x = 1) = cos(pi y)$
- $q(y = 0) = q(y = 1) = 0$

#image("../assets/laplace-2d/placa-mista.png", width: 80%)

$ T^"an" = sinh(pi x) cos(pi y) \/ sinh(pi) . $

Varie o número de elementos; tabela do erro percentual no ponto interno
$(sqrt(2)\/2, sqrt(2)\/2)$.

+ *E4 — mapa de cores.* Placa da figura: obtenha $T$ no contorno (+ internos se quiser) e
  produza um mapa (`plot_geo` e\/ou export para o Gmsh). Critério: figura legível + CDCs usadas no texto.

#image("../assets/laplace-2d/placa-mapa.png", width: 80%)

== Desafio

Com base neste #link("https://onlinelibrary.wiley.com/doi/epdf/10.1002/fld.1650030504")[artigo],
estime o coeficiente de sustentação de um perfil NACA via BEM.
#link("https://youtu.be/_G4yNayAPPE")[gravação]

== Extra: montagem hierárquica (H-matriz)

*Fora da trilha obrigatória* (também aparece em *Trabalhos finais*, proposta D).
Até $N tilde.eq 10^3$–$5 times 10^3$, a montagem *densa* `H_G_full_direct` basta para o curso.

Quando $N$ cresce, memória $O(N^2)$ e fatoração direta doem. O pacote oferece

```julia
msh = quadrado(ndiv=80, show=false, nome="quad_fine")
dad = format2d(msh, Laplace(1.0); pontointerno=false)
attach_analytical!(dad, ana_laplace_linear(; direction=SA[1.0, 0.0]))

H_G_Hmat(dad; atol=1e-6, nmax=32)
@show compression_ratio(dad.H)
solve(dad)                 # caminho iterativo quando A é hierárquica
@show rel_error(dad)
```

#table(
  columns: (auto, auto, auto),
  inset: 7pt,
  stroke: 0.5pt + luma(200),
  [*Critério*], [*Densa*], [*H-matriz*],
  [API], [`H_G_full_direct`], [`H_G_Hmat`],
  [Memória], [$O(N^2)$], [tipicamente $O(N log N)$],
  [Solver], [LU (`bem_linsolve`)], [GMRES \/ blocos (`solve_Hmat`)],
  [Curso 30 h], [padrão], [opcional \/ monografia],
)

Regra: comece denso e grosso; só compre quando o denso não couber ou demorar demais.
Detalhes: `docs/src/pt-br/performance.md` no repositório.

== Leituras e código

- `data/Laplace/Laplace_dad.jl` — `quadrado`, CDCs
- `src/Laplace/Fundamental.jl` — $T^*$, $q^*$
- `src/Core/Assembly_full.jl` — `H_G_full_direct`, diagonal de $H$
- `src/Core/Boundary_conditions.jl` — `applyBC`
- `src/Core/Solver.jl` — `solve`
- `src/Core/Analytical.jl` — `attach_analytical!`, `apply_analytical_bc!`
- Capítulo *Indo para 2D* (geometria) · *Poisson 2D* (fonte $f$) · *Medidas de erro*
