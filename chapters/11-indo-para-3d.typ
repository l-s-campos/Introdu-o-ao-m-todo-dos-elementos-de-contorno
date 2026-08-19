// Indo para 3D — geometria de superfície + exemplos BEM 3D simples
// Extensão do cap. Indo para 2D; código BEM_gmsh (format3d, geometric_props, Laplace 3D)
// Nome legado: Propgeo 3D

= Indo para 3D
<indo-para-3d>

No capítulo *Indo para 2D* o laboratório aprendeu a malhar só o contorno $Gamma$,
calcular $N_k$, $J$, $upright(bold(n))$ e reduzir área\/centróide a integrais em $Gamma$.

Esta é a *última aula* da trilha: o domínio vira um *sólido* $Omega subset RR^3$ e o
“contorno” uma *superfície* fechada — mas o *roteiro* é o mesmo do 2D.

#block(
  width: 100%,
  fill: rgb("#ecfdf5"),
  inset: 10pt,
  radius: 4pt,
  stroke: 0.5pt + rgb("#99f6e4"),
)[
  *Mensagem da aula.* 2D e 3D são *semelhantes* no `BEM_gmsh`:
  `format*` $arrow.r$ `geometric_props` $arrow.r$ `H_G_full_direct` $arrow.r$ `solve`.
  Mudam a dimensão de $Gamma$, a SF e o número de DOFs — não a lógica.
]

Nesta aula você:
1. gera superfície com Gmsh (`mesh_cube`, …);
2. monta o `dad` com `format3d`;
3. valida a malha com `geometric_props` ($A$, $V$, centróide);
4. repete o pipeline em *Laplace 3D* e, no exercício final, em *elasticidade 3D* simples.

== Objetivos

+ Ver o que muda (e o que *não* muda) de 2D para 3D.
+ Usar `format3d` e inspecionar nós, normais e faces.
+ Fechar geometria com `geometric_props` 3D (orientação inclusa).
+ Rodar Laplace 3D mínimo no cubo ($T = z$) com a *mesma* API do 2D.
+ Nos exercícios: um de geo, um de Laplace, um de elasticidade — fechamento do curso.

== Mapa

+ Por que “indo para 3D” (tabela 2D $arrow.r$ 3D)
+ Ambiente e `format3d`
+ Elemento de superfície: $J$ e $upright(bold(n))$
+ Propgeo: divergência e radial ($A$, $V$, $upright(bold(c))$)
+ Lab G — cubo (geometria)
+ Lab L — Laplace 3D ($T = z$)
+ Orientação e cavidades
+ Armadilhas
+ Exercícios finais (3)
+ Leituras

== Por que “indo para 3D”

#table(
  columns: (auto, auto, auto),
  inset: 7pt,
  stroke: 0.5pt + luma(200),
  [*Peça*], [*2D*], [*3D*],
  [Domínio $Omega$], [área no plano], [volume no espaço],
  [Contorno $Gamma$], [curva(s) fechada(s)], [superfície(s) fechada(s)],
  [Elementos], [1D em $xi in [-1,1]$], [2D na face ($xi,eta$); tri\/quad no Gmsh],
  [API malha], [`format2d`], [`format3d`],
  [Geo rápida], [`geometric_props` → $P,A,upright(bold(c))$], [`geometric_props` → $A,V,upright(bold(c))$],
  [Laplace], [`Laplace(k)` + SF $ln r$], [mesmo tipo + SF $1\/r$ em 3D],
  [DOFs escalares por nó], [$1$], [$1$ (potencial); elasticidade 3D teria $3$],
)

#block(
  width: 100%,
  fill: luma(248),
  inset: 10pt,
  radius: 4pt,
  stroke: 0.5pt + luma(180),
)[
  *Ideia-chave.* Malhar só $Gamma$ em 3D *não* elimina geometria: $N_k$, $J$ de superfície
  e $upright(bold(n))$ são os tijolos de `geometric_props` e, em seguida, de $H$ e $G$.
  Se $V$ e as normais estiverem errados, o BEM 3D herda o erro.
]

== Ambiente

```julia
using DrWatson
@quickactivate :BEM
include(datadir("Laplace", "cube_mesh.jl"))
```

(Na primeira vez: `Pkg.activate` no clone + `Pkg.instantiate`.)
GUI Gmsh: #link("https://gmsh.info/#Download")[gmsh.info] — útil para ver faces e grupos.

== `format3d`: do `.msh` ao `dad`

O Gmsh gera *elementos de dimensão 2* (faces). O `format3d` (`Input.jl`) lê esses
elementos, coloca nós de campo (descontínuos por face, grau `tipo`), calcula
jacobiano de superfície e normal:

```julia
# format3d — essência
# - getElements(2): quads 4 nós (tipo Gmsh 3) ou 9 nós (tipo 10)
# - NOS = N * nós_geométricos
# - n = (∂_ξ x) × (∂_η x)  normalizado
# - Physical Group nas *faces* (dim 2): CDC "0;T" / "1;q" se for Laplace
dad = format3d(msh, Laplace(1.0); tipo=1, pontointerno=false)
@show dad.dimension   # 3
@show dad.n           # nós de colocação de superfície
```

#table(
  columns: (auto, auto),
  inset: 7pt,
  stroke: 0.5pt + luma(200),
  [*Campo útil*], [*Papel*],
  [`dad.Nodes`], [pontos 3D de colocação],
  [`dad.Normal`], [normal unitária na face (orientação da malha)],
  [`dad.elements`], [conectividade + `Jacobian` por ponto do elemento],
  [`dad.elem_weight`], [pesos de quadratura reutilizados no propgeo default],
)

CDC nesta aula: os grupos do `mesh_cube` já trazem Dirichlet no topo\/base e
Neumann nos lados (para o Lab L). Para *só* geometria (Lab G), os nomes importam
pouco — importam $J$ e $upright(bold(n))$.

== Elemento de superfície: $J$ e $upright(bold(n))$

Em cada face, com parâmetros $(xi, eta)$:

$
upright(bold(x))(xi, eta) = sum_k N_k (xi, eta)\, upright(bold(X))_k ,
quad
partial_xi upright(bold(x)) , quad partial_eta upright(bold(x)) .
$

$
upright(bold(n)) parallel partial_xi upright(bold(x)) times partial_eta upright(bold(x)) ,
quad
J = | partial_xi upright(bold(x)) times partial_eta upright(bold(x)) | ,
quad
dif Gamma = J\, dif xi\, dif eta .
$

É o análogo 3D de $J = |d upright(bold(x))\/d xi|$ no contorno 2D. O sinal de
$upright(bold(n))$ fixa fora\/dentro — veja a seção de orientação.

== Propgeo: geometria só com $Gamma$

=== Divergência

$
integral_Omega nabla · upright(bold(F))\, dif V
=
integral_Gamma upright(bold(F)) · upright(bold(n))\, dif Gamma .
$

#table(
  columns: (auto, auto, auto),
  inset: 7pt,
  stroke: 0.5pt + luma(200),
  [*$upright(bold(F))$*], [*$nabla · upright(bold(F))$*], [*Uso*],
  [$upright(bold(x))\/3$], [$1$], [$V = (1\/3) integral_Gamma upright(bold(x)) · upright(bold(n))\, dif Gamma$],
  [—], [—], [$A = integral_Gamma dif Gamma$],
  [campos com $x_i$], [constantes], [momentos $arrow.r$ centróide],
)

=== Radial no pacote

Como no 2D, o código integra primitivas ao longo do raio $0 arrow.r upright(bold(x))$
e monta com o fator

$ j_"fac" = (upright(bold(n)) · upright(bold(x))) \/ R^3 ,
  quad R = |upright(bold(x))| . $

```julia
# GeometricProperties.jl — essência geometric_props_3d
# Area += wJ
# Fv,Fx,Fy,Fz = _calc_F_3d(x, ...)   # ∫₀ᴿ ρ² dρ e momentos
# Vol += Fv * dot(n,x)/R^3 * wJ
```

API:

```julia
gp = geometric_props(dad)
# GeometricProps3D(surface_area, volume, centroid)

gp2 = geometric_props(dad; npg_radial=12, npg_boundary=8)
# npg_boundary=nothing → pesos já no dad (default)
```

Pólo na *origem* dos nós: corpos longe de $upright(bold(0))$ ou não estrelados podem
pedir translação da malha ou mais cuidado na orientação.

== Lab G — Cubo unitário (só geometria)

Analítico: $A = 6$, $V = 1$, $upright(bold(c)) = (1\/2,1\/2,1\/2)$ se o cubo é $[0,1]^3$.

```julia
using DrWatson
@quickactivate :BEM
using LinearAlgebra
include(datadir("Laplace", "cube_mesh.jl"))

rows = []
for ndiv in (2, 3, 4, 6)
    msh = mesh_cube(L=1.0, ndiv=ndiv, show=false, nome="cube_G_$ndiv")
    dad = format3d(msh, Laplace(1.0); pontointerno=false)
    gp  = geometric_props(dad)
    gpB = geometric_props(dad; npg_boundary=8)
    push!(rows, (;
        ndiv, n=dad.n,
        A=gp.surface_area, V=gp.volume, c=gp.centroid,
        eA=abs(gp.surface_area - 6), eV=abs(gp.volume - 1),
        eV_b=abs(gpB.volume - 1),
    ))
end
display(rows)

# normais para fora? (cubo centrado em c)
dad = format3d(mesh_cube(L=1.0, ndiv=4, show=false), Laplace(1.0); pontointerno=false)
gp = geometric_props(dad)
c = gp.centroid
mean_out = sum(dot(dad.Normal[i], dad.Nodes[i] - c) for i in 1:dad.n) / dad.n
@show mean_out   # > 0 se n aponta para fora do sólido
```

*Esperado:* erros de $A$ e $V$ caindo com `ndiv` (ordens de grandeza do teste do repo:
poucos % em malha moderada). Se $V < 0$, inverta orientação global.

== Lab L — Laplace 3D mínimo ($T = z$ no cubo)

O mesmo cubo, agora com o solver. Solução exata compatível com as CDC do
`mesh_cube`: $T=0$ em $z=0$, $T=1$ em $z=L$, lados isolados $arrow.r$ $T = z\/L$
(com $L=1$, $T=z$).

Espelho de `scripts/cube_3d_laplace.jl`:

```julia
using DrWatson
@quickactivate :BEM
using LinearAlgebra
include(datadir("Laplace", "cube_mesh.jl"))

msh   = mesh_cube(L=1.0, ndiv=3, show=false)
props = Laplace(1.0)
dad   = format3d(msh, props; pontointerno=false)
println(dad)

# geometria primeiro (hábito)
gp = geometric_props(dad)
@show gp.surface_area gp.volume gp.centroid

# BEM 3D — mesma API do 2D
H_G_full_direct(dad; npg=8, threaded=true)
solve(dad)

# erro vs T = z
Terr = [dad.T[i] - dad.Nodes[i][3] for i in 1:dad.n]
err = norm(Terr) / max(norm(getindex.(dad.Nodes, 3)), eps())
@show err, extrema(dad.T)
```

*Variantes:* `ndiv in {2,3,4,5}`; `npg in {4,8,12}`; tabela `n`, `err`, tempo se quiser.
SF 3D no pacote: $T^* = 1\/(4 pi k R)$ (`Fundamental.jl`) — você não monta na mão.

#block(
  width: 100%,
  fill: rgb("#ecfdf5"),
  inset: 10pt,
  radius: 4pt,
  stroke: 0.5pt + rgb("#99f6e4"),
)[
  *Ordem de trabalho em 3D.* (1) `geometric_props` ok → (2) `H_G_full_direct` + `solve` →
  (3) erro de campo. Se (1) falha, não debugue singularidade de $H$ ainda.
]

=== CDC do cubo (referência)

Em `cube_mesh.jl` (faces = Physical Group dim 2):

```julia
# bottom z=0 → "0;0"   (Dirichlet T=0)
# top    z=L → "0;1"   (Dirichlet T=1)
# sides      → "1;0"   (Neumann q=0)
```

Mesma gramática do Laplace 2D (`"0;T"` \/ `"1;q"`), agora em *superfícies*.

== Orientação e cavidades

#table(
  columns: (auto, auto),
  inset: 7pt,
  stroke: 0.5pt + luma(200),
  [*Situação*], [*Efeito*],
  [$upright(bold(n))$ saindo de $Omega$ em toda $partial Omega$], [$V>0$; BEM com normal exterior padrão],
  [Normal global invertida], [$V$ muda de sinal; fluxos com sinal trocado],
  [Cavidade com $upright(bold(n))$ para o vazio (fora do material)], [Volume da cavidade subtraído — correto],
  [Cavidade com a mesma orientação do exterior], [Volume somado — errado],
)

2D: exterior anti-horário, furo horário. 3D: mão direita nas faces para $upright(bold(n))$
*para fora do material* em todas as componentes de $partial Omega$.

== Armadilhas

#table(
  columns: (auto, auto),
  inset: 7pt,
  stroke: 0.5pt + luma(200),
  [*Sintoma*], [*Causa típica*],
  [$V approx 0$ ou absurdo], [Superfície aberta; face faltando; $J=0$],
  [$V < 0$], [Orientação invertida],
  [$A$ ok, $V$ ruim], [Radial\/pólo; tente `npg_boundary`; transladar malha],
  [`format3d` “Unsupported element”], [Face não é quad 4\/9 no Gmsh atual],
  [Laplace com erro enorme, geo ok], [CDC de face; `npg` baixo; confunda $n$ nós],
  [Misturar `format2d` em `.msh` 3D], [Use `format3d`],
)

== Exercícios finais (3)

Última entrega da trilha. Apêndice de erros. Em *cada* item: reporte `dad.n` (e `dimension`),
o que mediu, e *uma linha* explicitando o paralelo 2D $arrow.r$ 3D
(ex.: “`format2d` $arrow.r$ `format3d`”; “$P,A$ $arrow.r$ $A,V$”; “$T$ $arrow.r$ $upright(bold(u))$ com 3 DOFs”).

=== E1 — Geometria (cubo ou paralelepípedo)

*2D análogo:* Lab de `geometric_props` no quadrado \/ coroa (*Indo para 2D*).

1. Gere o cubo unitário (`mesh_cube`, $L=1$) com `ndiv in {2,3,4,6}`.
2. `dad = format3d(msh, Laplace(1.0); pontointerno=false)`.
3. `gp = geometric_props(dad)` e, numa malha, também `npg_boundary=8`.
4. Tabela: `ndiv`, $n$, $A$, $V$, $upright(bold(c))$, $e_A=|A-6|$, $e_V=|V-1|$,
   $|upright(bold(c))-(1\/2,1\/2,1\/2)|$.
5. Confira `mean(n · (x - c)) > 0` (normal para fora).
6. *Opcional (+):* paralelepípedo $L_x times L_y times L_z$ e fórmulas analíticas de $A,V$.

*Critério:* $e_A$ e $e_V$ caindo com `ndiv`; sinal de $V$ correto; texto do paralelo 2D\/3D.

=== E2 — Laplace 3D (potencial no cubo)

*2D análogo:* quadrado com $T=x$ (*Laplace 2D*).

1. Mesmo cubo; CDC do `mesh_cube` ($T=0$ em $z=0$, $T=1$ em $z=1$, lados $q=0$).
2. Pipeline *idêntico* ao 2D:

```julia
H_G_full_direct(dad; npg=8)
solve(dad)
```

3. Erro vs $T = z$ nos nós (como no Lab L): $epsilon_2$ e $epsilon_infinity$ (apêndice).
4. Tabela `ndiv` (e, se puder, `npg`) × $n$ × erros.
5. *Antes* do solve: uma linha com `gp.volume` da mesma malha (geo ok $arrow.r$ campo).

*Critério:* erro de $T$ caindo com refino; menção explícita de que a API é a do Laplace 2D.

=== E3 — Elasticidade 3D (patch de Dirichlet no cubo)

*2D análogo:* patch $upright(bold(u)) = bold(epsilon)\, upright(bold(x))$ (*Elasticidade 2D*).

Objetivo: ver que o pipeline *continua o mesmo* com 3 DOFs por nó.

1. Superfície do cubo unitário (`mesh_cube` + `format3d`), agora com

```julia
dad = format3d(msh, Elasticity(1.0, 0.3, 1.0); tipo=1, pontointerno=false)
```

2. Patch de deformação uniforme (ex. $epsilon_(x x)=0.01$, demais nulos):
   em *todo* nó de superfície e em *cada* direção $d=1,2,3$,

```julia
# u = ε · x  (Dirichlet em todos os DOFs de contorno)
for i in 1:dad.n
    x = dad.Nodes[i]
    u = (εxx * x[1], εxy * x[1] + εyy * x[2], ...)  # escolha um ε simples
    for d in 1:3
        dad.BC[3*(i-1)+d] = 0          # Dirichlet
        dad.BV[3*(i-1)+d] = u[d]
    end
end
```

   Sugestão mínima: $epsilon_(x x)=0.01$, $epsilon_(y y)=epsilon_(z z)=epsilon_(i j)=0$
   $arrow.r$ $upright(bold(u)) = (0.01\, x,\, 0,\, 0)$.

3. De novo a *mesma* sequência:

```julia
H_G_full_direct(dad; npg=8)
solve(dad)
```

4. Erro em deslocamento: $epsilon_2$ e $epsilon_infinity$ empilhando os 3 componentes
   (ou `rel_error` se anexar analítico compatível). Compare 2 refinamentos de `ndiv`.

5. No relatório, complete a tabela mental do curso:

#table(
  columns: (auto, auto, auto),
  inset: 7pt,
  stroke: 0.5pt + luma(200),
  [*Etapa*], [*2D (aulas anteriores)*], [*3D (esta aula)*],
  [Malha], [`format2d`], [`format3d`],
  [Geo], [$P,A,upright(bold(c))$], [$A,V,upright(bold(c))$],
  [Laplace], [$T$ \/ $q$; `"0;T"`], [$T$ \/ $q$ em *faces*],
  [Elasticidade], [$upright(bold(u)),upright(bold(t))$ 2 DOFs], [$upright(bold(u)),upright(bold(t))$ *3* DOFs],
  [Montagem \/ solve], [`H_G_full_direct` + `solve`], [*iguais*],
)

*Critério:* patch com erro pequeno ou decrescente; CDC 3 DOFs documentadas; tabela 2D\/3D preenchida com suas palavras.

#block(
  width: 100%,
  fill: luma(248),
  inset: 10pt,
  radius: 4pt,
  stroke: 0.5pt + luma(200),
)[
  *Fechamento.* Se E1–E3 rodaram, você repetiu em 3D o arco inteiro do curso:
  geometria de contorno $arrow.r$ potencial $arrow.r$ elasticidade — com o mesmo código de montagem.
  Trabalhos finais aprofundam um tema; a trilha de 30 h encerra aqui.
]

== Depois do curso

- Trabalhos finais (propostas A–E no `BEM_gmsh`).
- H-matrizes, trinca, contato, multirregião — quando a geometria 2D\/3D já estiver sólida.
- Internos (`pontointerno=true`) e pós-processamento volumétrico.

== Leituras e código

- Cap. *Indo para 2D* — radial, divergência, `format2d` (espelho desta aula)
- Cap. *Laplace 2D* — `H_G_full_direct` + `solve` (mesma API no Lab L \/ E2)
- Cap. *Elasticidade 2D* — patch Dirichlet (espelho do E3 com 2 DOFs)
- *Apêndice: medidas de erro*
- `src/Core/Input.jl` — `format3d`
- `src/Core/GeometricProperties.jl` — `geometric_props_3d`, `_calc_F_3d`
- `data/Laplace/cube_mesh.jl` — `mesh_cube`
- `scripts/cube_3d_laplace.jl` — Lab L
- `src/Laplace/Fundamental.jl` — SF 3D $1\/(4 pi k R)$
- `test/test_new_features.jl` — $A approx 6$, $V approx 1$
- #link("https://1drv.ms/f/s!AmfyGvdmTYongqYn5kjjlZaMHr9h2w?e=z0sXvU")[Arquivos legados propgeo]
- #link("https://youtu.be/Uc-rxXDBU6I")[gravação]
