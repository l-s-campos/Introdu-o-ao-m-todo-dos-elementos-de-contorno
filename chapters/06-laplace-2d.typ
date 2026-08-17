// Laplace 2D — Gmsh + CDC + BEM_gmsh
// Geometria de contorno: capítulo "Indo para 2D"

= Laplace 2D
<laplace-2d>

Geometria, $N_k$, $J$ e $upright(bold(n))$ já vieram de *Indo para 2D*.
Aqui o Gmsh deixa de ser só malha e passa a carregar *CDC*; montamos $H,G$ e resolvemos o primeiro BEM 2D completo.

== Objetivos

+ Ativar o `BEM_gmsh` e rodar o fluxo mínimo do quadrado ($T = x$).
+ Ler CDCs nos *grupos físicos* do `.geo` \/ `.msh` (`"0;T"`, `"1;q"`).
+ Seguir o pipeline PDE $arrow.r$ equação integral $arrow.r$ $H T = G q$ $arrow.r$ $A x = b$ $arrow.r$ internos.
+ Entender a diagonal de $H$ (corpo rígido) e quando usar H-matriz.

== Setup

Código: #link("https://github.com/l-s-campos/BEM_gmsh")[`BEM_gmsh`].
Clone, abra Julia *na pasta do projeto* e prepare o ambiente:

```julia
using Pkg
Pkg.activate(".")
Pkg.instantiate()

using DrWatson
@quickactivate :BEM
include(datadir("Laplace", "Laplace_dad.jl"))
```

GUI do sistema (opcional): #link("https://gmsh.info/#Download")[gmsh.info].

== Fluxo mínimo (quadrado unitário, $T = x$)

```julia
props = Laplace(1.0)
msh   = quadrado(ndiv=20, show=false)
dad   = format2d(msh, props)

attach_analytical!(dad, ana_laplace_linear(; direction=SA[1.0, 0.0]))
H_G_full_direct(dad, 20)   # densa; H_G_Hmat(dad) se N for grande
solve(dad)

@show rel_error(dad)
plot_geo(dad)
# export_results_to_gmsh(dad, msh, :T; viewer=false)
gp = geometric_props(dad)
@show gp.perimeter gp.area gp.centroid
```

Os grupos físicos da malha já carregam as CDCs. A geometria foi preparada em *Indo para 2D*; a convenção de CDC é desenvolvida *neste* capítulo.

== Grupos físicos $arrow.r$ CDC

No `BEM_gmsh`, o *nome* da curva física no Gmsh carrega tipo e valor da condição de contorno.

#table(
  columns: (auto, auto),
  inset: 8pt,
  stroke: 0.5pt + luma(200),
  [*Nome*], [*Significado (Laplace)*],
  [`"0;T"`], [Dirichlet: potencial \/ temperatura $T$ prescrita],
  [`"1;q"`], [Neumann: fluxo $q = -k partial T \/ partial n$ prescrito],
)

Exemplos: `"0;100"`, `"1;0"` (isolado). Elasticidade usa `"tx;ux;ty;uy"` (capítulo de elasticidade).

#table(
  columns: (auto, auto),
  inset: 7pt,
  stroke: 0.5pt + luma(200),
  [*Onde*], [*Papel do Gmsh \/ `BEM_gmsh`*],
  [*Indo para 2D*], [malha, $N_k$, $J$, $upright(bold(n))$, perímetro \/ área \/ centróide],
  [*Este capítulo*], [CDC nos grupos, montagem $H,G$, `solve`, erro],
)

Boas práticas (repete o essencial da geometria):
1. Contorno externo anti-horário; furos horários.
2. Elementos descontínuos no campo (`format2d`) — cantos sem nó compartilhado de CDC.
3. `ordem` da malha alinhada a `tipo`.
4. Convergência: variar `ndiv` e medir erro (`rel_error` \/ capítulo de erros).

== Mapa do BEM (leia isto antes da álgebra)

Antes das fórmulas, fixe o *pipeline* inteiro. Cada bloco das seções seguintes preenche *um* passo desta cadeia — inclusive o cálculo da diagonal de $H$, que é só um detalhe técnico do passo 4.

#block(
  width: 100%,
  fill: luma(248),
  inset: 12pt,
  radius: 6pt,
  stroke: 0.6pt + luma(200),
)[
  #set par(justify: false)
  #set text(size: 10.5pt)
  *1. PDE* $quad nabla^2 T = 0$ em $Omega$, com CDCs em $Gamma = partial Omega$

  $arrow.b.double$

  *2. Resíduos ponderados + integração por partes (Green)* $quad$
  o operador sai de $T$ e vai para a função peso $v$

  $arrow.b.double$

  *3. Solução fundamental (SF)* $quad$
  escolhe-se $v = T^*$ tal que $-nabla^2 T^* = delta(x - x_d)$;
  o domínio some e sobra uma *equação integral só no contorno*

  $arrow.b.double$

  *4. Discretização* $quad$
  $T$ e $q$ ≈ funções de forma nos elementos; integrais → matrizes $H$ e $G$
  (aqui entram quadratura, singularidades e a *diagonal de* $H$)

  $arrow.b.double$

  *5. Condições de contorno (CDC)* $quad$
  em cada nó, $T$ *ou* $q$ é conhecido → reorganiza-se $H T = G q$ em $A x = b$

  $arrow.b.double$

  *6. Solve* $quad x = A^(-1) b$ $arrow.r$ contorno completo $(T, q)$ em todos os nós

  $arrow.b.double$

  *7. Pontos internos (pós-processamento)* $quad$
  com $(T, q)$ no contorno já conhecidos, $T(x_d)$ no interior é só uma integral
  — *sem* novo sistema linear
]

Correspondência com o `BEM_gmsh`:

#table(
  columns: (auto, auto),
  inset: 7pt,
  stroke: 0.5pt + luma(200),
  [*Passo*], [*No código*],
  [1–3 formulação], [já embutida em `fundamental` / montagem],
  [4 $H, G$], [`H_G_full_direct(dad, npg)` ou `H_G_Hmat(dad)`],
  [5 CDC], [grupos físicos Gmsh + `format2d` / `attach_analytical!`],
  [6 solve], [`solve(dad)`],
  [7 internos], [automático se `pontointerno=true`; ver também `plot_geo`],
)

#block(
  width: 100%,
  fill: rgb("#ecfdf5"),
  inset: 10pt,
  radius: 4pt,
  stroke: 0.5pt + rgb("#99f6e4"),
)[
  *Ideia-chave.* No MEF montamos “rigidez” no *volume*. No BEM montamos $H$ e $G$ no *contorno*; o preço é que essas matrizes são cheias e as integrais podem ser singulares quando o ponto fonte $x_d$ cai no elemento que está sendo integrado.
]

== Formulação

=== Passos 1–3: da PDE à equação integral

A equação de Laplace é dada por:

$ nabla^2 T = 0, $

usando o método dos resíduos ponderados e aplicando a segunda identidade de Green (com $u = T$ e peso $v$):

$ integral_Omega (v nabla^2 T - T nabla^2 v) d Omega = integral_Gamma (v (partial T)/(partial n) - T (partial v)/(partial n)) d s $

obtém-se a equação integral de contorno:

$ 0.5 T(x_d , y_d) = integral_Gamma T q^* d s - integral_Gamma T^* q d s, $

onde $T^*$ e $q^*$ são as soluções fundamentais:

$ T^* = (-1)/(2 pi k) ln r $

$ q^* = 1/(2 pi r^2) [(x - x_d) n_x + (y - y_d) n_y] . $

O coeficiente $0.5$ à esquerda vale para ponto fonte $x_d$ *sobre* um contorno liso. É o termo de corpo rígido / ângulo sólido: geometricamente, metade do “salto” da solução fundamental fica de cada lado do contorno. (Em cantos o fator não é $1/2$; por isso o código prefere elementos *descontínuos* ou calcula a diagonal de $H$ de forma indireta — passo 4.)

=== Passo 4: discretização → matrizes $H$ e $G$

Podemos representar a temperatura e fluxo em um elemento descontinuo com  $m$  nós como:
$T = N_1 T_1 + N_2 T_2 + N_3 T_3 + ... + N_m T_m$
$q = N_1 q_1 + N_2 q_2 + N_3 q_3 + ... + N_m q_m$
onde estamos aproximando nossa distribuição de temperatura e fluxo no elemento por uma função polinomial de  grau  $(m - 1)$ .

Discretizando em $n$ elementos de contorno descontínuos:

$0.5 T(x_d , y_d) = sum_(j = 1)^(n_(e l e m)) [integral_(Gamma_j) T q^* d Gamma] - sum_(j = 1)^(n_(e l e m)) [integral_(Gamma_j) T^* q d Gamma]$

Usando a representação de   $T$ e $q$  na equação acima temos:

$ 0.5 T(x_d , y_d) = sum_(j = 1)^(n_(e l e m)) {integral_(Gamma_j) [mat(delim: #none, N_1, N_2, N_3, dots.c, N_m)] [mat(delim: #none, T_1; T_2; T_3; dots.v; T_m)]_j q^* d Gamma} \ -
sum_(j = 1)^(n_(e l e m)) {integral_(Gamma_j) T^* [mat(delim: #none, N_1, N_2, N_3, dots.c, N_m)] [mat(delim: #none, q_1; q_2; q_3; dots.v; q_m)]_j d Gamma} $

Repetindo essa equação para cada diferente nó-fonte $x_d = x_i$ montamos um sistema matricial:

$ H T = G q $

- Coluna ligada a $T$: integrais de $q^* N_k$ → entradas de $H$ (núcleo *mais* singular).
- Coluna ligada a $q$: integrais de $T^* N_k$ → entradas de $G$ (singularidade fraca, integrável).
- A linha $i$ é a equação integral escrita com ponto fonte no nó $i$.

==== Diagonal de $H$: por que é especial?

Quando o ponto fonte $x_i$ está *no mesmo elemento* que o ponto de integração, $r -> 0$ e $q^* ~ 1/r$ (em 2D) deixa de ser uma integral comum. Integrar $H_(i i)$ “na marra” é possível, mas delicado (transformações, partes finitas, etc.).

==== Método indireto (corpo a temperatura constante)

Em vez de atacar a singularidade forte de frente, usamos uma solução *exata* trivial do problema de Laplace:

$ T(x) equiv 1, quad q(x) equiv 0 . $

Ela satisfaz $nabla^2 T = 0$ e a equação integral. No sistema discreto isso vira

$ H {1} = G {0} = {0} $.

Ou seja, *cada linha de* $H$ *soma zero*. Como os termos fora da diagonal $H_(i j)$ ($i != j$) são integrais regulares (já calculadas por Gauss), a diagonal sai de graça:

$ H_(i i) = - sum_(j = 1, j != i)^N H_(i j) , quad i = 1, ..., N . $

Isso *já inclui* o fator de ângulo sólido (o “$0.5$” no contorno liso, ou outro valor em cantos).

#block(
  width: 100%,
  fill: luma(248),
  inset: 10pt,
  radius: 4pt,
  stroke: 0.5pt + luma(200),
)[
  *Ordem prática na montagem (passo 4):*
  1. zere $H$ e $G$;
  2. para cada par (nó-fonte $i$, elemento $j$), integre e *some* nas colunas certas de $H$ e $G$ — *pule* ou trate à parte o caso singular $i in$ elemento $j$ em $H$;
  3. corrija a diagonal: $H_(i i) = -sum_(j != i) H_(i j)$;
  4. (opcional) a diagonal de $G$ tem só singularidade fraca e em geral *é* integrada com quadratura adequada — o truque indireto raramente é necessário.
]

No `BEM_gmsh`, os passos 2–3 estão dentro de `H_G_full_direct`.

=== Passos 5–6: CDCs e o sistema $A x = b$

Depois de montar $H T = G q$, *ainda não se resolve nada*: em cada nó falta uma informação. A CDC diz qual coluna vai para o lado esquerdo (incógnita) e qual contribui para o lado direito (dado).

==== Exemplo (1 elemento por lado)

A fim de ilustrar como se aplica as condições de contorno e se calcula as variáveis desconhecidas será analisado um problema de condução de calor unidirecional com uma discretização de um elemento por lado.

#image("../assets/laplace-2d/exemplo-unidirecional.png", width: 80%)

As equações obtidas podem ser escritas na forma matricial, como:

$ (mat(delim: #none, H_11, H_12, H_13, H_14; H_21, H_22, H_23, H_24; H_31, H_32, H_33, H_34; H_41, H_42, H_43, H_44; ))
(mat(delim: #none, macron(T_1); T_2; macron(T_3); T_4))
=
(mat(delim: #none, G_11, G_12, G_13, G_14; G_21, G_22, G_23, G_24; G_31, G_32, G_33, G_34; G_41, G_42, G_43, G_44; ))
(mat(delim: #none, q_1; macron(q_2); q_3; macron(q_4))) $

onde $macron(T)$ e $macron(q)$ são termos *conhecidos* (CDC).

Separando os termos conhecidos dos desconhecidos:

$ (mat(delim: #none, -G_11, H_12, -G_13, H_14; -G_21, H_22, -G_23, H_24; -G_31, H_32, -G_33, H_34; -G_41, H_42, -G_43, H_44; ))(mat(delim: #none, q_1; T_2; q_3; T_4)) = (mat(delim: #none, -H_11, G_12, -H_13, G_14; -H_21, G_22, -H_23, G_24; -H_31, G_32, -H_33, G_34; -H_41, G_42, -H_43, G_44; ))(mat(delim: #none, macron(T_1); macron(q_2); macron(T_3); macron(q_4))) $

Assim, pode-se escrever $A x = b$. Resolvendo o sistema linear obtém-se as variáveis que faltavam no contorno. No `BEM_gmsh` isso é o `solve(dad)`.

Regra prática por nó $i$:
- se $T_i$ é preescrito (Dirichlet): a incógnita é $q_i$ → coluna $i$ de $-G$ entra em $A$, e $H_(: i) T_i$ vai para $b$;
- se $q_i$ é preescrito (Neumann): a incógnita é $T_i$ → coluna $i$ de $H$ entra em $A$, e $G_(: i) q_i$ vai para $b$.

=== Passo 7: pontos internos

A equação integral para pontos *interiores* ($x_d in Omega$, fora de $Gamma$) é ligeiramente modificada — o coeficiente à esquerda vira $1$, não $0.5$:

$ T(x_d , y_d) = integral_Gamma T q^* d s - integral_Gamma T^* q d s $

Detalhe importante: com $T$ e $q$ *já conhecidos em todo o contorno* (passos 5–6), a temperatura em qualquer ponto interno é só avaliar essas integrais. *Não* se monta nem se resolve um novo $A x = b$.

== Exemplo resolvido ponta a ponta (`BEM_gmsh`)

Problema: quadrado unitário, $k = 1$, CDCs tais que a solução exata é $T(x,y) = x$ (e portanto $q = - partial T / partial n$). Malha Gmsh com `ndiv=16`, montagem densa, impressão de números e gráfico.

```julia
using DrWatson
@quickactivate :BEM
using LinearAlgebra, Printf
include(datadir("Laplace", "Laplace_dad.jl"))

# --- 1. Física e malha ---
props = Laplace(1.0)                       # k = 1
msh   = quadrado(ndiv=16, show=false)      # grava .msh e devolve o caminho
dad   = format2d(msh, props; pontointerno=true)

println("Arquivo de malha: ", msh)
println("Nós de contorno N = ", dad.n)
println("Pontos internos   = ", length(dad.internalNodes))

# --- 2. Campo analítico T = x  (bate com as CDCs padrão do quadrado) ---
ana = ana_laplace_linear(; direction=SA[1.0, 0.0])
attach_analytical!(dad, ana)

# --- 3. Montagem H, G (densa) e solução ---
H_G_full_direct(dad, 16)                   # npg = 16
solve(dad)

# --- 4. Números na tela ---
err = rel_error(dad)
@printf "Erro relativo (rel_error) = %.6e\n" err

println("\nPrimeiros nós do contorno:")
println(" i |      x |      y |    T_num |  T_exato |     q_num")
for i in 1:min(8, dad.n)
    x, y = dad.Nodes[i]
    T_ex = x                               # T = x
    @printf "%2d | %6.3f | %6.3f | %8.5f | %8.5f | %9.5f\n" i x y dad.T[i] T_ex dad.q[i]
end

# Erro pontual máximo em T no contorno
eT = maximum(abs(dad.T[i] - dad.Nodes[i][1]) for i in 1:dad.n)
@printf "\nmax |T_num - x| no contorno = %.6e\n" eT

# --- 5. Gráfico ---
fig = plot_geo(dad)                        # visualização do solver (backend do projeto)
# save("laplace_quadrado.png", fig)        # descomente para gravar PNG
fig
```

O que você deve observar:
- `rel_error` cai ao aumentar `ndiv` (faça `ndiv = 8, 16, 32` e anote numa tabela);
- nos lados verticais, $T approx x$ (constante em cada lado); nos horizontais, $q approx 0$;
- pontos internos, se houver, seguem $T approx x$ sem novo sistema linear.

== Montagem densa vs. H-matriz

Até $N tilde.eq 3 dot.op 10^3$–$5 dot.op 10^3$ nós, a montagem *densa* (`H_G_full_direct`) costuma ser a escolha didática e prática: implementação simples, fatoração LU direta, ótima para estudos de convergência.

Limitações do BEM denso:
- memória $O(N^2)$ para armazenar $H$ e $G$;
- tempo de montagem e de fatoração também $O(N^2)$–$O(N^3)$;
- em malhas finas (furos, camadas limite geométricas, 3D) o custo explode antes da física ficar interessante.

Quando $N$ cresce, o `BEM_gmsh` oferece montagem *hierárquica*:

```julia
msh = quadrado(ndiv=80, show=false, nome="quad_fine")
dad = format2d(msh, Laplace(1.0); pontointerno=false)
attach_analytical!(dad, ana_laplace_linear(; direction=SA[1.0, 0.0]))

H_G_Hmat(dad; atol=1e-6, nmax=32)   # blocos comprimidos
@show compression_ratio(dad.H)
solve(dad)                          # GMRES sobre operador misto
@show rel_error(dad)
```

#table(
  columns: (auto, auto, auto),
  inset: 7pt,
  stroke: 0.5pt + luma(200),
  [*Critério*], [*Densa* `H_G_full_direct`], [*H-matriz* `H_G_Hmat`],
  [Memória], [$O(N^2)$], [tipicamente $O(N log N)$],
  [Uso no curso], [padrão; $N tilde.eq 10^2$–$10^3$], [malhas grandes; projeto avançado],
  [Solver], [LU direta], [iterativo (GMRES)],
  [Controle], [`npg`], [`atol`, `nmax`, …],
)

Regra prática: comece *sempre* denso e com malha grossa; só mude para H-matriz quando o denso não couber na RAM ou demorar demais. Detalhes em `docs/src/pt-br/performance.md` do repositório.

== Exercícios

Use as medidas de erro do capítulo *Medidas de erro* (e/ou `rel_error(dad)` do `BEM_gmsh`). Notação: potencial $T$, fluxo $q = -k partial T / partial n$.

*Malha e CDC*
1. Rode `data/examples/geo_unit_square.jl` e confira `geometric_props`.
2. Abra um `.geo` de `data/Laplace/` e separe geometria de CDC (`"0;..."` \/ `"1;..."`).

*Problemas de Laplace*

+ Resolva esse problema com elementos lineares, quadráticos e cúbicos e calcule o erro médio, o erro máximo e a norma $L_2$ do erro no contorno para diferentes discretizações. Faça um gráfico (Plots) comparando a convergência dos três tipos de elementos. A solução analítica é dada por:

$ T(theta) &= theta/pi \ q(x) &= - 1/(pi x) quad "em" quad y = 0 $

#image("../assets/laplace-2d/setor-circular.png", width: 80%)

+ Analise o seguinte problema de placa com condições de contorno mistas:

- $T(x = 0) = 0$
- $T(x = 1) = cos(pi y)$
- $q(y = 0) = q(y = 1) = 0$

#image("../assets/laplace-2d/placa-mista.png", width: 80%)

A solução analítica é:

$ T^"an" = sinh(pi x) cos(pi y) / sinh(pi) $

Varie o número de elementos e faça uma tabela com o erro percentual no ponto interno $(x, y) = (sqrt(2) / 2, sqrt(2) / 2)$.

+ Faça um mapa de cor da distribuição de temperatura em uma placa com dimensões e condições de contorno mostradas na figura.

#image("../assets/laplace-2d/placa-mapa.png", width: 80%)

== Desafio

Baseado nesse #link("https://onlinelibrary.wiley.com/doi/epdf/10.1002/fld.1650030504")[artigo] calcule o coeficiente de sustentação de um perfil NACA usando o BEM.

#link("https://youtu.be/_G4yNayAPPE")[gravação]
