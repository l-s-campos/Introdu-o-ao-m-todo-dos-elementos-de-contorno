// Propgeo 3D — RASCUNHO (não incluído em main.typ / site até aprovação)
// Extensão do cap. Indo para 2D: área, volume e centróide só com a superfície
// Código: BEM_gmsh GeometricProperties.jl + format3d

= Propgeo 3D
<propgeo-3d>

No capítulo *Indo para 2D* perímetro, área e centróide saíram de integrais *apenas*
no contorno $Gamma = partial Omega$ (radial ou divergência), com o `dad` 2D.

Aqui o domínio é um *sólido* $Omega subset RR^3$. A fronteira é uma *superfície*
fechada (malha de elementos 2D em 3D). Volume, área superficial e centróide
continuam sendo integrais *só em* $Gamma$ — sem tetraedros de volume.

Isso é o “propgeo 3D”: a mesma filosofia do BEM (informação no contorno), agora
só geometria.

== Objetivos

+ Ligar as fórmulas de divergência 3D à API `geometric_props`.
+ Entender a redução *radial* usada no pacote (`_calc_F_3d` + fator $upright(bold(n))·upright(bold(x))\/R^3$).
+ Rodar `mesh_cube` + `format3d` + `geometric_props` e conferir $A=6$, $V=1$ no cubo unitário.
+ Controlar orientação (normal exterior) e o efeito de cavidades.
+ Comparar quadratura do `dad` vs `npg_boundary`.

== Mapa

+ Por que só a superfície
+ Fórmulas (divergência e radial)
+ Malha 3D no `BEM_gmsh` (`format3d`)
+ API `geometric_props` \/ `geometric_props_3d`
+ Lab: cubo unitário
+ Orientação e furos
+ Armadilhas
+ Exercícios
+ Leituras

== Por que só a superfície

Identidade da divergência em 3D: para um campo vetorial $upright(bold(F))$ regular em $Omega$,

$
integral_Omega nabla · upright(bold(F))\, dif V
=
integral_Gamma upright(bold(F)) · upright(bold(n))\, dif Gamma
$

($upright(bold(n))$ exterior). Escolhas clássicas:

#table(
  columns: (auto, auto, auto),
  inset: 7pt,
  stroke: 0.5pt + luma(200),
  [*$upright(bold(F))$*], [*$nabla · upright(bold(F))$*], [*Resultado em $Gamma$*],
  [$upright(bold(x))\/3$], [$1$], [$V = (1\/3) integral_Gamma upright(bold(x)) · upright(bold(n))\, dif Gamma$],
  [$x_i upright(bold(e))_i$ (sem soma)], [constante], [momentos $arrow.r$ centróide],
  [$upright(bold(0))$ trivial], [—], [área: $A = integral_Gamma dif Gamma$ (só medida de superfície)],
)

#block(
  width: 100%,
  fill: rgb("#ecfdf5"),
  inset: 10pt,
  radius: 4pt,
  stroke: 0.5pt + rgb("#99f6e4"),
)[
  *Elo com o BEM.* Em potencial 3D, $H$ e $G$ também vivem em $Gamma$. Propgeo 3D
  treina a mesma infraestrutura: nós de superfície, normais, jacobianos de
  elementos 2D em $RR^3$, orientação — *sem* ainda montar SF.
]

== Duas leituras da mesma integral

=== Forma “divergência pura” (didática)

$
A = integral_Gamma dif Gamma ,
quad
V = 1/3 integral_Gamma upright(bold(x)) · upright(bold(n))\, dif Gamma ,
quad
upright(bold(x))_G = (1\/V) integral_Omega upright(bold(x))\, dif V .
$

Os momentos de volume $integral_Omega x_i\, dif V$ também reduzem a $Gamma$ com
$upright(bold(F))$ adequado (ex. $F = (x_i^2 \/ 2)\, upright(bold(e))_i$).

=== Forma “radial” (o que o pacote implementa)

Como no 2D, fixa-se um pólo (na prática a *origem* das coordenadas dos nós) e
integra-se ao longo do raio até a superfície. Em cada ponto de Gauss superficial
$upright(bold(x))$ com $R = |upright(bold(x))|$:

1. Primitivas radiais ao longo do segmento $0 arrow.r upright(bold(x))$:

```julia
# GeometricProperties.jl — _calc_F_3d (essência)
# ρ = x * q,  q ∈ [0,1];  fac ~ ρ² dρ
# Fv ~ ∫₀ᴿ ρ² dρ   (→ volume quando montado com o jacobiano angular)
# Fx ~ ∫₀ᴿ ρ_x ρ² dρ   etc. (momentos)
```

2. Fator de “ângulo sólido \/ fluxo radial” na superfície:

$
j_"fac" = (upright(bold(n)) · upright(bold(x))) \/ R^3 .
$

3. Soma:

$
V approx sum (F_V\, j_"fac"\, w J) ,
quad
integral x_i dif V approx sum (F_i\, j_"fac"\, w J) ,
quad
A approx sum w J .
$

4. Centróide: $upright(bold(x))_G = (integral upright(bold(x))\, dif V) \/ V$ se $V != 0$.

#block(
  width: 100%,
  fill: luma(248),
  inset: 10pt,
  radius: 4pt,
  stroke: 0.5pt + luma(200),
)[
  *Pólo na origem.* A fórmula radial assume que os raios de $upright(bold(0))$ a $Gamma$
  cobrem $Omega$ de forma consistente (corpo *estrelado* em relação à origem, ou
  pelo menos que o sinal de $upright(bold(n))·upright(bold(x))$ marque dentro\/fora).
  *Traduza a malha* se o sólido estiver longe de $upright(bold(0))$ e o erro crescer
  sem motivo — ou confie na forma de divergência equivalente quando o código a
  reproduz via as mesmas primitivas.
]

Para $f equiv 1$ e geometria simples, radial e divergência *coincidem* (a menos de
quadratura). Em 2D isso já foi o par “área radial vs Green”.

== Malha de superfície e `format3d`

Em 3D o BEM (e o propgeo) leem *elementos de superfície* do Gmsh (`getElements(2)`):
quads 4 nós (tipo 3) ou 9 nós (tipo 10) no `format3d` atual.

```julia
# Input.jl — format3d (essência)
function format3d(filename, properties; tipo=1, pontointerno=true, ...)
    # nós 3D; laço nos elementos de superfície
    # NOS, normal = shapefun2D × geometria; n = ∂_ξx × ∂_ηx normalizado
    # CDC nos Physical Group *de face* (dim 2): "0;T" / "1;q" se for Laplace
    return BEMdata(...; dimension=3, ...)
end
```

Gerador pronto de cubo (`data/Laplace/cube_mesh.jl`):

```julia
include(datadir("Laplace", "cube_mesh.jl"))
msh = mesh_cube(L=1.0, ndiv=4, show=false)
# faces com grupos físicos de CDC (para Laplace 3D);
# para *só* propgeo, os nomes de CDC não importam — importam nós, J, n
```

`dad.dimension == 3` faz `geometric_props(dad)` despachar para `geometric_props_3d`.

== API no `BEM_gmsh`

```julia
using DrWatson
@quickactivate :BEM
include(datadir("Laplace", "cube_mesh.jl"))

msh = mesh_cube(L=1.0, ndiv=4, show=false)
dad = format3d(msh, Laplace(1.0); pontointerno=false)

gp = geometric_props(dad)
# == geometric_props_3d(dad; npg_radial=12, npg_boundary=nothing)

@show gp.surface_area   # ≈ 6
@show gp.volume         # ≈ 1
@show gp.centroid       # ≈ (0.5, 0.5, 0.5)
```

Retorno:

```julia
struct GeometricProps3D
    surface_area::Float64
    volume::Float64
    centroid::SVector{3,Float64}
end
```

=== Palavras-chave

#table(
  columns: (auto, auto),
  inset: 7pt,
  stroke: 0.5pt + luma(200),
  [*Keyword*], [*Efeito*],
  [`npg_radial=12`], [Gauss na primitiva radial $0 arrow.r R$ (`_calc_F_3d`)],
  [`npg_boundary=nothing`], [Usa `elem.Jacobian` × `dad.elem_weight` já no `dad` (default)],
  [`npg_boundary=N`], [Reintegra cada face com $N times N$ Gauss (tri via mapa Duffy; quad via `shapefun2D`)],
)

Trecho do laço default (`npg_boundary=nothing`):

```julia
# geometric_props_3d — essencial
for elem in dad.elements
    for k in eachindex(elem.index)
        x = dad.Nodes[i];  n = dad.Normal[i]
        wJ = elem.Jacobian[k] * w_el[k]
        Area += wJ
        R = norm(x);  R < 1e-14 && continue
        jfac = dot(n, x) / R^3
        Fv, Fx, Fy, Fz = _calc_F_3d(x, qsi_r, w_r)
        Vol += Fv * jfac * wJ
        # xdV, ydV, zdV += ...
    end
end
c = Vol == 0 ? 0 : (xdV, ydV, zdV) ./ Vol
```

== Lab — cubo unitário

Alvo analítico: $A = 6$, $V = 1$, $upright(bold(x))_G = (1\/2,1\/2,1\/2)$.

```julia
using DrWatson
@quickactivate :BEM
using LinearAlgebra
include(datadir("Laplace", "cube_mesh.jl"))

rows = []
for ndiv in (2, 3, 4, 6)
    msh = mesh_cube(L=1.0, ndiv=ndiv, show=false, nome="cube_geo_$ndiv")
    dad = format3d(msh, Laplace(1.0); pontointerno=false)
    gp  = geometric_props(dad)
    gp2 = geometric_props(dad; npg_boundary=8)
    push!(rows, (;
        ndiv,
        n = dad.n,
        A = gp.surface_area,
        V = gp.volume,
        c = gp.centroid,
        eA = abs(gp.surface_area - 6),
        eV = abs(gp.volume - 1),
        eA2 = abs(gp2.surface_area - 6),
        eV2 = abs(gp2.volume - 1),
    ))
end
display(rows)
```

*Esperado (ordem de grandeza do teste do repo):* $A$ com erro relativo de poucos %
já em malha grossa; $V$ um pouco mais sensível à quadratura\/orientação
(`rtol` ~ 0.05–0.15 no teste sintético). Refine `ndiv` e compare `npg_boundary`.

Sanidade extra:

```julia
# normais apontando para fora? média de n · (x - center) > 0
center = gp.centroid
mean_out = sum(dot(dad.Normal[i], dad.Nodes[i] - center) for i in 1:dad.n) / dad.n
@show mean_out   # deve ser > 0 se n é exterior
```

== Orientação e cavidades

#table(
  columns: (auto, auto),
  inset: 7pt,
  stroke: 0.5pt + luma(200),
  [*Situação*], [*Efeito em $V$*],
  [Todas as faces com $upright(bold(n))$ saindo de $Omega$], [$V > 0$ correto],
  [Malha inteira com normal invertida], [$V$ muda de *sinal* (e centróide estranho)],
  [Cavidade com normal *saindo do material* (para o furo)], [Volume do furo *subtraído* — correto],
  [Cavidade com a mesma orientação do exterior], [Volume do furo *somado* — errado],
)

Em 2D: contorno externo anti-horário, furo horário. Em 3D: regra da mão direita
nas faces de modo que $upright(bold(n))$ aponte *para fora do material* em *todas*
as componentes de $partial Omega$ (exterior e paredes de vazios).

== Armadilhas

#table(
  columns: (auto, auto),
  inset: 7pt,
  stroke: 0.5pt + luma(200),
  [*Sintoma*], [*Causa típica*],
  [$V approx 0$ ou absurdo], [Superfície aberta; furos sem tampa; normal nula],
  [$V < 0$], [Orientação global invertida],
  [$A$ ok, $V$ ruim], [Pólo \/ radial; refine `npg_radial` ou `npg_boundary`],
  [Centróide longe do óbvio], [$V$ errado propaga; cheque sinal de $V$ primeiro],
  [`format3d` erro de tipo], [Elemento de superfície não quad 4\/9 no Gmsh atual],
  [CDC aparecem no `.geo`], [Irrelevante para propgeo puro; não confundir com volume],
)

== Exercícios

Use erros relativos $|A-A_"ex"|\/A_"ex"$, $|V-V_"ex"|\/V_"ex"$ e
$|upright(bold(c))-upright(bold(c))_"ex"|$ (apêndice de erros, normas em $RR^3$).

=== E0 — Lab do cubo

Tabela `ndiv` × $A,V,upright(bold(c)), e_A, e_V$ com default e com `npg_boundary=8`
(ou 12). Uma frase: o que converge mais rápido, $A$ ou $V$?

=== E1 — Paralelepípedo

$L_x, L_y, L_z$ distintos (adapte `mesh_cube` ou `.geo` próprio).
Analítico: $A = 2(L_x L_y + L_y L_z + L_z L_x)$, $V = L_x L_y L_z$,
centróide no centro geométrico se a origem for o canto $(0,0,0)$.

=== E2 — Esfera facetada

Gmsh: esfera (ou icosaedro refinado). Compare $A$ e $V$ com $4 pi R^2$ e
$(4\/3) pi R^3$. Estude erro vs número de faces.

=== E3 — Cubo com cavidade

Cubo externo + cavidade cúbica ou cilíndrica *fechada*.
(a) $V approx V_"ext" - V_"cav"$; (b) inverta só a orientação da cavidade e relate o sinal;
(c) área superficial = soma das áreas (externa + parede da cavidade).

=== E4 — (Opcional) Divergência “na mão”

Com os `Nodes`, `Normal`, `Jacobian` do `dad`, implemente

$ V_3 = (1\/3) sum_i (upright(bold(x))_i · upright(bold(n))_i)\, w J_i $

e compare com `gp.volume` do radial. Discrepâncias: quadratura e definição de $w J$.

== O que *não* é esta aula

- Montagem $H,G$ 3D e Laplace 3D no cubo (`scripts/cube_3d_laplace.jl`) — próximo
  passo natural, mas *depois* da geometria estar correta.
- Propgeo é verificação de malha\/orientação *antes* de confiar no BEM 3D.

== Leituras e código

- Cap. *Indo para 2D* — radial, divergência, `geometric_props` 2D, orientação de furos
- `src/Core/GeometricProperties.jl` — `geometric_props_3d`, `_calc_F_3d`
- `src/Core/Input.jl` — `format3d`
- `data/Laplace/cube_mesh.jl` — `mesh_cube`
- `scripts/cube_3d_laplace.jl` — uso de `format3d` (solver, não só geo)
- `test/test_new_features.jl` — teste $A approx 6$, $V approx 1$ no cubo
- #link("https://1drv.ms/f/s!AmfyGvdmTYongqYn5kjjlZaMHr9h2w?e=z0sXvU")[Arquivos legados]
- #link("https://youtu.be/Uc-rxXDBU6I")[gravação]
