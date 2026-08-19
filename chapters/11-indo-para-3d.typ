// Indo para 3D — geometria de superfície + exemplos BEM 3D simples
// Extensão do cap. Indo para 2D; código BEM_gmsh (format3d, geometric_props, Laplace 3D)
// Nome legado: Indo para 3D

= Indo para 3D
<indo-para-3d>

No capítulo *Indo para 2D* o laboratório aprendeu a malhar só o contorno $Gamma$,
calcular $N_k$, $J$, $upright(bold(n))$ e reduzir área\/centróide a integrais em $Gamma$.

Aqui o domínio é um *sólido* $Omega subset RR^3$. O “contorno” passa a ser uma
*superfície* fechada (malha de elementos 2D em $RR^3$). O foco continua sendo
*geometria + pipeline do `BEM_gmsh`*:

1. gerar superfície com Gmsh (`mesh_cube`, …);
2. montar o `dad` com `format3d`;
3. validar malha com `geometric_props` (área, volume, centróide — o antigo “propgeo”);
4. rodar um *Laplace 3D mínimo* (cubo, $T = z$) com o mesmo `H_G_full_direct` + `solve` do 2D.

CDC detalhada e problemas 3D avançados ficam para projetos; o essencial é *sair do plano*
sem perder o fio do curso.

== Objetivos

+ Ver o que muda de 2D para 3D na malha e no `BEMdata`.
+ Usar `format3d` e inspecionar nós, normais e elementos de face.
+ Obter $A$, $V$ e centróide só com a superfície (`geometric_props` 3D).
+ Conferir orientação (normal exterior) e o papel de cavidades.
+ Rodar dois exemplos simples: (G) só geometria no cubo; (L) Laplace $T=z$ no cubo.
+ Comparar erros de geometria e de potencial (apêndice de erros).

== Mapa

+ Por que “indo para 3D”
+ Ambiente e `format3d`
+ Elemento de superfície: $J$ e $upright(bold(n))$
+ Indo para 3D: divergência e radial
+ Lab G — cubo unitário (só geometria)
+ Lab L — Laplace 3D mínimo ($T = z$)
+ Orientação e cavidades
+ Armadilhas
+ Exercícios
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

== Indo para 3D (geometria só com $Gamma$)

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

== Exercícios

Apêndice de erros. Sempre reporte `dad.n` e o que mediu ($A,V$ ou $T$).

=== E0 — Labs G e L

(i) Tabela do Lab G (`ndiv`, $e_A$, $e_V$, default vs `npg_boundary`).
(ii) Tabela do Lab L (`ndiv`, `npg`, erro $T=z$).
(iii) Uma frase: a geo já era boa quando o potencial convergiu?

=== E1 — Paralelepípedo

$L_x, L_y, L_z$ distintos. Geo: $A$ e $V$ analíticos.
Opcional: Laplace com $T = z\/L_z$ (CDC análogas ao cubo).

=== E2 — Esfera facetada

Só geometria: $A$ e $V$ vs $4 pi R^2$ e $(4\/3) pi R^3$; erro vs número de faces.

=== E3 — Cubo com cavidade

(a) $V approx V_"ext" - V_"cav"$ com orientação correta;
(b) inverta só a cavidade e relate o sinal de $V$;
(c) (opcional) um potencial simples se montar CDC coerentes.

=== E4 — Divergência na mão

Com `Nodes`, `Normal`, pesos do `dad`:

$ V_3 = (1\/3) sum_i (upright(bold(x))_i · upright(bold(n))_i)\, (w J)_i $

Compare com `gp.volume`.

== O que fica para depois

- Elasticidade 3D, H-matrizes 3D, problemas industriais — trabalhos \/ pesquisa.
- Internos (`pontointerno=true`) e pós-processamento volumétrico.
- SF e singularidades 3D em profundidade (além do uso via `H_G_full_direct`).

== Leituras e código

- Cap. *Indo para 2D* — radial, divergência, `format2d`, descontínuo
- Cap. *Laplace 2D* — pipeline $H,G$, CDC, `solve` (idêntico na API)
- *Apêndice: medidas de erro*
- `src/Core/Input.jl` — `format3d`
- `src/Core/GeometricProperties.jl` — `geometric_props_3d`, `_calc_F_3d`
- `data/Laplace/cube_mesh.jl` — `mesh_cube`
- `scripts/cube_3d_laplace.jl` — Lab L
- `src/Laplace/Fundamental.jl` — SF 3D $1\/(4 pi k R)$
- `test/test_new_features.jl` — $A approx 6$, $V approx 1$
- #link("https://1drv.ms/f/s!AmfyGvdmTYongqYn5kjjlZaMHr9h2w?e=z0sXvU")[Arquivos legados propgeo]
- #link("https://youtu.be/Uc-rxXDBU6I")[gravação]
