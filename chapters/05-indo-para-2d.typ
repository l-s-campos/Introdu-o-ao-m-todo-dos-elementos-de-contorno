// Indo para 2D + Gmsh (geometria)
// Pacote do curso: codes/geom_gmsh (GeomGmsh). CDC ficam para Laplace 2D.
//

= Indo para 2D
<indo-para-2d>

Nessa aula usamos o código
#link("https://1drv.ms/u/c/278a4d66f71af267/IQCDqMd75dfHSJ3QtZrz22ikASAEDJ8N_AtmrycBd_sFGtM?e=CM4iSb")[`GeomGmsh`] (Julia + Gmsh; recorte de geometria do `BEM_gmsh`).
O foco é *geometria de contorno*: malha, elementos, jacobiano, normal e
integrais só em $Gamma$.
*Condições de contorno* (o que cada aresta “vale” em $T$ ou $q$) ficam para o capítulo *Laplace 2D*.

Os trechos Julia abaixo são *pedaços do próprio pacote* (`src/`) ou dos
exemplos em `scripts/` e `_solu/`. Os mesmos nomes aparecem no código-fonte.

== Objetivos

+ Ativar o ambiente `GeomGmsh` e gerar uma malha de contorno com Gmsh.
+ Entender o elemento descontínuo de ordem $p$ via `discontinuous_nodes_weights` + `shapefun` (Lagrange baricêntrico).
+ Calcular $J$, tangente e normal como em `format2d`.
+ Obter perímetro, área e centróide com `geometric_props` nas *duas* estratégias: `:radial` (default) e `:divergence`.
+ Rodar os exemplos `scripts/circulo.jl` (2D) e `scripts/toro.jl` (3D).

== Mapa

+ Por que 2D e o que muda em relação ao BEM 1D
+ Ambiente `GeomGmsh` e Gmsh (geometria, sem CDC)
+ Elemento descontínuo + Lagrange baricêntrico (`Interpolation.jl` / `Input.jl`)
+ Jacobiano, normal (como em `format2d`)
+ `geometric_props` — radial e divergência (`GeometricProperties.jl`)
+ Exemplo: disco unitário
+ Exemplo: toro 3D (`format3d`)
+ Exercícios

== Por que “indo para 2D”

No BEM 1D o “contorno” eram dois pontos. Em 2D:

- $Omega subset RR^2$ tem fronteira $Gamma = partial Omega$ *unidimensional* (curva fechada, eventualmente com furos);
- a geometria de $Gamma$ é aproximada por *elementos de contorno* $Gamma_e$;
- campo e fluxo no contorno usam interpolação em $xi in [-1,1]$ sobre cada $Gamma_e$;
- integrais de área em $Omega$ que admitem redução a $Gamma$ (radial ou divergência) preparam as identidades do BEM.

#image("../assets/indo-para-2d/perimetro.png", width: 75%)

#block(
  width: 100%,
  inset: 10pt,
  radius: 4pt,
  stroke: 0.5pt + luma(180),
  fill: luma(248),
)[
  *Ideia-chave.* Malhar só $Gamma$ não elimina a necessidade de boa geometria:
  $N_k$, $J$ e $upright(bold(n))$ em cada elemento são os tijolos de `geometric_props` e, depois, de $H$ e $G$.
]

== Ambiente: `GeomGmsh` + Gmsh

Na pasta `codes/geom_gmsh`:

```bash
julia --project=. -e "using Pkg; Pkg.instantiate()"
julia --project=. scripts/circulo.jl
```

O pacote `Gmsh.jl` já é dependência. A GUI do sistema
(#link("https://gmsh.info/#Download")[gmsh.info]) ajuda a inspecionar `.msh`.

Geradores (`circulo`, `anel`, `placa_com_furo`, `cubo`, `toro`, …) em `src/Meshes.jl`
seguem o mesmo modo de sessão: `gmsh.initialize()` → gerar → gravar
`datadir(nome * ".msh")` → `show && gmsh.fltk.run()` → `gmsh.finalize()`.
Depois `format2d` / `format3d` e `geometric_props`.

=== Por que Gmsh no BEM (só geometria)

1. Curvas, arcos, splines, furos e sólidos OCC (toro) com orientação controlável.
2. Malha de contorno (curvas) + superfície 3D; volume opcional (pontos internos / orientação).
3. Ordem geométrica alinhada ao grau do elemento (`ordem` / `tipo`).
4. O mesmo gerador alimenta Laplace, Poisson e elasticidade no `BEM_gmsh` completo.

*Nesta aula* os grupos físicos *nomeiam* pedaços de curva e a superfície do domínio.
A convenção `"0;T"` / `"1;q"` entra em *Laplace 2D*.

=== Orientação (crítico)

- Contorno *externo*: *anti-horário* → normal *para fora* de $Omega$.
- Contorno de *furo*: *horário* → normal ainda saindo de $Omega$.
- $A < 0$ (ou $V < 0$ em 3D) nas fórmulas abaixo denuncia orientação invertida.

`format2d(; orient=true)` (default) corrige o sinal de $upright(bold(n))$ com as
células 2D da malha de superfície. `placa_com_furo(; reverse_hole=true)` já
inverte o laço do furo no Gmsh.

== O que `format2d` faz na geometria

Trecho essencial de `src/Input.jl`:

```julia
# nós de colocação descontínuos = Gauss–Legendre (p+1 pontos)
qsi, wi = discontinuous_nodes_weights(p)   # = gausslegendre(p+1)

# geometria isoparamétrica nos nós Equispaced de grau p
Ngeo, dNgeo = shapefun(Equispaced(p), qsi)

# para cada aresta da malha Gmsh, com vértices X:
nodes[idx]  .= Ngeo * X          # posições dos nós de campo
dx           = dNgeo * X
J            = norm.(dx)
normal[idx]  = tan2normal.(dx ./ J)
L            = abs(dot(J, wi))   # comprimento do elemento
```

Ou seja: *os mesmos* `shapefun` e pesos que você verá na montagem de $H$ e $G$.

Fluxo geométrico mínimo — o mesmo de `scripts/circulo.jl`:

```julia
using GeomGmsh

msh = circulo(; ndiv=16, show=false, nome="ex_circulo", ordem=2)
dad = format2d(msh, Laplace(1.0); tipo=2, pontointerno=false)
g = geometric_props(dad)          # strategy=:radial
println("  perimeter = ", g.perimeter, "  (exact 2π)")
println("  area      = ", g.area, "  (exact π)")
println("  centroid  = ", g.centroid, "  (exact 0,0)")
plot_geo(dad; show_normals=true)
```

`format2d` lê a malha e monta `BEMdata` (elementos, nós de colocação, normais, pesos).
*Nesta aula* use `plot_geo` e `geometric_props`, não `solve` / `H_G_*`.

#image("../assets/indo-para-2d/circulo-geo.png", width: 80%)

Saída típica (`tipo=2`, `ndiv=16`): $P approx 6.2827$ (exato $2 pi$),
$A approx 3.1411$ (exato $pi$), centróide na origem.

== Elemento de contorno descontínuo

No `GeomGmsh` o padrão é *elemento descontínuo*: graus de liberdade de campo
*não* são compartilhados nos vértices entre elementos vizinhos
(evita conflito de CDC em cantos — detalhe em Laplace 2D).

#image("../assets/indo-para-2d/elemento-parabolico.png", width: 70%)

=== Mapa de referência

$ Gamma_e $ é a imagem de $xi in [-1,1]$:

$
  upright(bold(x))(xi) = sum_(k=1)^(p+1) N_k (xi)\, upright(bold(x))_k^((e)) .
$

- *Geometria* da aresta Gmsh: nós `Equispaced(p)` (vértices + meios, se ordem 2).
- *Campo* (colocação): nós `Legendre(p)` = Gauss–Legendre em $(-1,1)$ via
  `discontinuous_nodes_weights(p)`.

```julia
# src/Input.jl
function discontinuous_nodes_weights(p::Integer)
    p >= 1 || error("degree must be ≥ 1, got $p")
    return gausslegendre(p + 1)
end
```

`tipo=1,2,3` em `format2d` escolhe o grau $p$ do campo (e alinha a ordem da malha se preciso).

=== Lagrange baricêntrico (como no pacote)

O núcleo está em `src/Interpolation.jl`: pesos baricêntricos + matriz de interpolação
(Berrut–Trefethen) e `shapefun` (= $N$ e $d N\/d xi$ nos pontos pedidos).

Pesos (para nós $x_0$ do polinômio):

$ w_i = 1 \/ product_(j != i) (x_i - x_j) . $

Funções de forma no ponto $xi$ (fora dos nós):

$
  N_i (xi) = (w_i \/ (xi - x_i)) \/ sum_j (w_j \/ (xi - x_j)) .
$
#image("../assets/indo-para-2d/funcoes-forma-d.png", width: 70%)

No código:

```julia
# Interpolation.jl — interpolation_matrix (forma baricêntrica)
# M[j,i] = w[i] / (xx - x0[i]), normalizado pela soma da linha
# se xx == nó: linha de Kronecker

function shapefun(poly::AbstractPolynomial, x)
    L = interpolation_matrix(poly, x)   # N(x) nos pontos pedidos
    L, L * poly.Dmat                    # N e dN/dξ
end
```

=== `Dmat` (resumo)

Com valores nodais $u_i = u(xi_i)$, o interpolante $u(xi)=sum_i N_i(xi) u_i$ tem

$ u'(xi_k) = sum_i N_i'(xi_k) u_i quad arrow.r.double quad upright(bold(u))' = D upright(bold(u)) , $

onde $D_(k i) = N_i'(xi_k)$. No pacote, $D$ = `poly.Dmat`, montada *uma vez* por `diff_matrix`.

Com os pesos baricêntricos $w_i$ e nós $x_i$:

$
  D_(k i) = (w_i \/ w_k) \/ (x_k - x_i) quad (k != i),
  wide
  D_(k k) = - sum_(i != k) D_(k i)
$

(a diagonal impõe $sum_i N_i' = 0$: derivada de constante é zero).

`shapefun` devolve $N$ e $d N\/d xi$ em pontos *quaisquer* $xi$ via

```julia
L = interpolation_matrix(poly, x)   # N
dN = L * poly.Dmat                  # N_i'(ξ) = Σ_k N_k(ξ) D_ki
```

`Dmat` *não* é o jacobiano da malha: age no interpolante 1D em $xi$. O jacobiano geométrico vem depois,
$dif upright(bold(x)) \/ dif xi = sum_i (dif N_i \/ dif xi) upright(bold(X))_i$.

Uso didático (os mesmos tipos do `dad`):

```julia
p = 2
poly = Legendre(p)                 # campo descontínuo (nós de Gauss)
ξ = nodes(poly)                    # nós de colocação
wN = weights(poly)                 # pesos *baricêntricos* (não de Gauss)
ξg, wg = gausslegendre(6)
N, dN = shapefun(poly, ξg)         # size (6 × 3) se p=2

poly_geo = Equispaced(p)           # geometria contínua da aresta
Ngeo, dNgeo = shapefun(poly_geo, ξg)
```

Partição da unidade (teste rápido, exercício E1):

```julia
poly = Legendre(3)
for ξ in range(-1, 1; length=21)
    N, dN = shapefun(poly, ξ)
    @assert abs(sum(N) - 1) < 1e-11
    @assert abs(sum(dN)) < 1e-9
end
```

Caso contínuo clássico $p=2$, nós $xi in {-1,0,1}$ (`Equispaced(2)`):

$
  N_1 = xi(xi-1)\/2 , quad N_2 = 1-xi^2 , quad N_3 = xi(xi+1)\/2 .
$

#image("../assets/indo-para-2d/funcoes-forma.png", width: 70%)

#block(
  width: 100%,
  inset: 10pt,
  radius: 4pt,
  stroke: 0.5pt + luma(180),
  fill: luma(248),
)[
  *No `BEMdata`.* Depois de `format2d`, `dad.element_type` é o polinômio do campo
  (`Legendre(p)`), `dad.elem_weight` são os pesos de Gauss da colocação, e cada
  `dad.elements[e]` guarda índices dos nós, jacobianos nos nós de colocação e comprimento.
  `dad.Nodes` / `dad.Normal` são vistas de `collocation` (nós de contorno; internos se `pontointerno=true`).
]

== Jacobiano, comprimento e normal

Com nós geométricos $upright(bold(X))_k$ da aresta e $N, N'$ em $xi$:

$
  upright(bold(x))' (xi) = sum_k N_k' (xi)\, upright(bold(X))_k ,
  quad
  J = |upright(bold(x))'| ,
  quad
  upright(bold(n)) = (y', -x') \/ J
$

(normal “à esquerda” do sentido de percurso — contorno externo anti-horário ⇒ exterior).

No pacote (`format2d` + `tan2normal`):

```julia
dx = dNgeo * X
J  = norm.(dx)
normal = tan2normal.(dx ./ J)   # (dx,dy) → (dy, -dx)/|·|  (unitário)
L  = abs(dot(J, wi))            # ∫ J dξ ≈ Σ J_k w_k
```

Perímetro global: $P = sum_e L_e$ (é o `g.perimeter`).

#image("../assets/indo-para-2d/normal-1.png", width: 32%)
#image("../assets/indo-para-2d/normal-2.png", width: 32%)

== Propriedades geométricas: duas estratégias

`geometric_props` (`src/GeometricProperties.jl`) aceita

```julia
g  = geometric_props(dad)                         # strategy=:radial
g2 = geometric_props(dad; strategy=:divergence)
g3 = geometric_props(dad; npg_radial=20)
g4 = geometric_props(dad; npg_boundary=16)
```

`npg_boundary=nothing` reutiliza `elem.Jacobian` $times$ `dad.elem_weight`
(colocação). Um inteiro reintegra cada elemento com `shapefun(dad.element_type, ξ)`
e Gauss.

=== Integração radial (`strategy=:radial`, default)

Fixe o pólo na *origem*. Para $upright(bold(x)) in Gamma$,
$upright(bold(r))=upright(bold(x))$, $r=|upright(bold(r))|$, $hat(upright(bold(r)))=upright(bold(r))\/r$.

$
  dif A = rho \, dif rho \, dif theta ,
  quad
  dif theta = (upright(bold(n)) · hat(upright(bold(r)))) (dif Gamma) \/ r .
$

#image("../assets/indo-para-2d/area-polar.png", width: 70%)
#image("../assets/indo-para-2d/area-contorno.png", width: 70%)

$
  I = integral_Omega f \, dif A
  = integral_Gamma F(upright(bold(x))) thin (upright(bold(n)) · upright(bold(r))) / r^2 dif Gamma ,
$

$
  F = integral_0^r f(rho hat(upright(bold(r)))) \, rho \, dif rho .
$

No código, $F$ para $f=1$, $f=x$, $f=y$ é numérico na direção radial
(`_calc_F_2d` com `npg_radial` pontos de Gauss em $rho$):

```julia
# GeometricProperties.jl — _geom_accum_2d_radial
# r = |x|,  r̂ = x/r,  nr = n · r̂
# Fa, Fx, Fy = ∫_0^r {1, x, y} ρ dρ     (_calc_F_2d)
# A  += Fa * nr / r * wJ
# Sx += Fx * nr / r * wJ
# Sy += Fy * nr / r * wJ
# centróide = (Sx, Sy) / A
```

```julia
function _calc_F_2d(r, theta, qsi, w)
    dro = r / 2
    Fa = Fx = Fy = 0.0
    for i in eachindex(qsi)
        ρ = r / 2 * (qsi[i] + 1)      # ρ ∈ [0, r]
        x, y = ρ * cos(theta), ρ * sin(theta)
        Fa += ρ * dro * w[i]          # ∫ ρ dρ  (f = 1)
        Fx += x * ρ * dro * w[i]      # ∫ x ρ dρ
        Fy += y * ρ * dro * w[i]
    end
    return Fa, Fx, Fy
end
```

Para $f=1$ a forma fechada é $F = r^2\/2$, e recupera-se
$A = 1/2 integral_Gamma (upright(bold(n)) · upright(bold(x))) dif Gamma$.
Com `npg_radial=12` (default) essa primitiva polinomial já é exata.

=== Teorema da divergência (`strategy=:divergence`)

$
  integral_Omega nabla · upright(bold(F)) \, dif A
  = integral_Gamma upright(bold(F)) · upright(bold(n)) \, dif Gamma .
$

#table(
  columns: (auto, auto, auto),
  inset: 7pt,
  stroke: 0.5pt + luma(200),
  [*$upright(bold(F))$*], [*$nabla · upright(bold(F))$*], [*Resultado*],
  [$(x,0)$ ou $(0,y)$], [$1$], [$A = integral_Gamma x n_x dif Gamma$],
  [$(x,y)\/2$], [$1$], [$A = 1/2 integral_Gamma upright(bold(x)) · upright(bold(n)) dif Gamma$],
  [$(x^2\/2, 0)$], [$x$], [$integral_Omega x dif A = integral_Gamma (x^2\/2) n_x dif Gamma$],
  [$(0, y^2\/2)$], [$y$], [$integral_Omega y dif A = integral_Gamma (y^2\/2) n_y dif Gamma$],
  [$(0, y^3\/3)$], [$y^2$], [$I_x = integral_Omega y^2 dif A = integral_Gamma (y^3\/3) n_y dif Gamma$],
)

É o acumulador `_geom_accum_2d_divergence` (forma fechada; `npg_radial` é ignorado):

```julia
P += wJ
A += 0.5 * dot(x, n) * wJ
Sx += 0.5 * x[1]^2 * n[1] * wJ
Sy += 0.5 * x[2]^2 * n[2] * wJ
```

A outra fórmula de área da tabela ($A = integral_Gamma x n_x dif Gamma$) está em `area_divergence`:

```julia
function area_divergence(dad::BEMdata)
    # mesmos (x, n, wJ) da colocação
    A = 0.0
    w_el = dad.elem_weight
    for elem in dad.elements
        for k in eachindex(elem.index)
            i = elem.index[k]
            wJ = elem.Jacobian[k] * w_el[k]
            A += dad.Nodes[i][1] * dad.Normal[i][1] * wJ
        end
    end
    return A
end
```

Comparar no disco (ou na coroa, exercício E2):

```julia
g = geometric_props(dad)                          # :radial
gd = geometric_props(dad; strategy=:divergence)
Adiv = area_divergence(dad)
@show g.area gd.area Adiv
```

No `scripts/circulo.jl` as três áreas coincidem até o erro de quadratura de $Gamma$
($approx 3.14110$).

#block(
  width: 100%,
  inset: 10pt,
  radius: 4pt,
  stroke: 0.5pt + luma(180),
  fill: luma(248),
)[
  *Radial × divergência.* Para $f=1,x,y$ coincidem (primitiva polinomial + mesma
  quad. em $Gamma$). Radial generaliza $f$ via $F$. Divergência é imediata quando
  existe $upright(bold(F))$ polinomial simples ($I_x$ via `moment_Ix`, etc.).
]

== Exemplo 2D: disco unitário

Arquivo `scripts/circulo.jl`. Malha de quatro arcos (`circulo` em `Meshes.jl`),
ordem 2, 16 elementos 1D.

```julia
msh = circulo(; ndiv=16, show=false, nome="ex_circulo", ordem=2)
dad = format2d(msh, Laplace(1.0); tipo=2, pontointerno=false)
g = geometric_props(dad)
plot_geo(dad; show_normals=true, title="unit disk tipo=2 ndiv=16")
```

#table(
  columns: (auto, auto, auto),
  inset: 7pt,
  stroke: 0.5pt + luma(200),
  [*Grandeza*], [*Numérico*], [*Exato*],
  [perímetro], [$6.2827$], [$2 pi approx 6.2832$],
  [área], [$3.1411$], [$pi approx 3.1416$],
  [centróide], [$approx (0,0)$], [$(0,0)$],
)

A malha linear (`tipo=1`) é um $n$-ágono regular:
$P = 2 n sin(pi\/n) = 2 pi - O(n^(-2))$,
$A = (n\/2) sin(2 pi\/n) = pi - O(n^(-2))$.
Com `tipo=2` os arcos parabólicos aproximam o círculo de verdade e o erro cai bem mais rápido.

== Exemplo 3D: toro

A *mesma* API vale em 3D: gerador → `format3d` → `geometric_props`.
Arquivo `scripts/toro.jl`. O gerador `toro` (`Meshes.jl`) chama
`gmsh.model.occ.addTorus` (raio maior $R$, raio do tubo $r$, eixo $z$, centro na origem),
malha a superfície em *quads* (recombine) e, se possível, o volume (para orientar $upright(bold(n))$).

Analítico:

$
  S = 4 pi^2 R r , quad V = 2 pi^2 R r^2 , quad upright(bold(c)) = upright(bold(0)) .
$

```julia
using GeomGmsh

R, r = 2.0, 0.5
msh = toro(; nome="ex_toro", R=R, r=r, lc=0.25, show=false)
dad = format3d(msh, Laplace(1.0); tipo=1, pontointerno=false)
g  = geometric_props(dad)                         # :radial
gd = geometric_props(dad; strategy=:divergence)
plot_geo(dad; show_normals=false, show_nodes=false)
```

Em 3D o acumulador radial usa $dif V = rho^2 dif rho\, dif Omega$ e

$
  j_"fac" = (upright(bold(n)) · upright(bold(x))) \/ R^3 ,
  quad R = |upright(bold(x))| ,
$

com `_calc_F_3d` para $F = integral_0^R f rho^2 dif rho$. A divergência é
$V = (1\/3) integral_Gamma upright(bold(x)) · upright(bold(n))\, dif Gamma$.

#image("../assets/indo-para-2d/toro-geo.png", width: 80%)

#table(
  columns: (auto, auto, auto),
  inset: 7pt,
  stroke: 0.5pt + luma(200),
  [*Grandeza*], [*Numérico (`lc=0.25`)*], [*Exato ($R=2$, $r=1\/2$)*],
  [área de superfície], [$39.18$], [$4 pi^2 R r approx 39.48$],
  [volume], [$9.57$], [$2 pi^2 R r^2 approx 9.87$],
  [centróide], [$approx (0,0,0)$], [$(0,0,0)$],
)

Radial e divergência coincidem até o roundoff (o volume difere na 15ª casa).
O cubo unitário (`cubo` + `format3d`) é o outro teste 3D do pacote: $S=6$, $V=1$,
$upright(bold(c))=(1\/2,1\/2,1\/2)$ — aparece no exercício E5 e no capítulo *Indo para 3D*.

`format3d` lê faces Gmsh *quad* (tipo 3 linear ou tipo 10 quadrático).
A normal de superfície é

$
  upright(bold(n)) parallel partial_xi upright(bold(x)) times partial_eta upright(bold(x)) ,
  quad
  J = | partial_xi upright(bold(x)) times partial_eta upright(bold(x)) | .
$

```julia
# Input.jl — format3d (essência)
N, dNx, dNy = shapefun2D(Equispaced(1), qsi)
nodes[idx]  .= N * X
dx1, dx2     = dNx * X, dNy * X
J            = norm.(cross.(dx1, dx2))
normal[idx]  = cross.(dx1, dx2) ./ J
```

== Boas práticas de malha (geometria)

1. Contorno externo anti-horário; furos horários.
2. Elementos descontínuos no campo (padrão `format2d`) — cantos sem nó compartilhado de CDC.
3. `ordem` da malha Gmsh alinhada a `tipo` em `format2d`.
4. Refino perto de cantos e furos (`placa_com_furo` como modelo).
5. Nesta aula: `pontointerno=false` basta para prop. de $Gamma$.
6. Em 3D, volume no `.msh` ajuda `format3d` a orientar $upright(bold(n))$ (`toro` tenta `generate(3)`).

== Exercícios

Entrega: scripts no ambiente `GeomGmsh` (`julia --project=codes/geom_gmsh`),
uso de `plot_geo` quando fizer sentido, números com erro vs analítico e 2–3 frases
de interpretação.
Fontes: `src/{Input,Interpolation,GeometricProperties,Meshes}.jl`,
`scripts/circulo.jl`, `scripts/toro.jl`. Soluções de referência em `_solu/`.

+ *E1 — `shapefun` e partição da unidade.*
  Para $p in {1,2,3,4}$:
  (a) crie `poly = Legendre(p)` e, em 21 abscissas de $[-1,1]$, verifique
  $sum_i N_i = 1$ e $sum_i N_i' = 0$ via `shapefun`;
  (b) plote as colunas de $N(xi)$ (`Plots`) para $p=2$ e $p=3$;
  (c) com `Equispaced(2)` e $xi in {-1,0,1}$, confira $N$ com
  $xi(xi-1)\/2$, $1-xi^2$, $xi(xi+1)\/2$;
  (d) compare `nodes(Legendre(p))` com `discontinuous_nodes_weights(p)[1]` — são o quê?
  (`weights(poly)` são baricêntricos; o segundo retorno de `discontinuous_nodes_weights` é a quad. de Gauss.)

+ *E2 — Coroa circular: geometria, ordem e duas fórmulas de área.*
  Domínio $1 <= r <= 2$ (`anel`). Analítico: $P = 6 pi$, $A = 3 pi$, centróide na origem.

  (a) Malhas `ndiv in {8,16,32}` e `tipo in {1,2}`. Para cada par, `geometric_props(dad)`
  e erros relativos de $P$ e $A$.

  (b) No *mesmo* `dad`, compare
  `geometric_props(dad; strategy=:divergence).area` e `area_divergence(dad)` com `g.area`
  (`:radial`). Há diferença sistemática ao mudar `tipo`?

  (c) `geometric_props(dad; npg_boundary=16)` e `geometric_props(dad; npg_radial=20)`
  vs o default. Quando a reintegração de $Gamma$ muda $A$? E `npg_radial`?

  (d) `plot_geo` na malha mais grossa e na mais fina; as normais no furo apontam para fora de $Omega$?

  (e) Tabela `ndiv, tipo, n_col, e_P, e_A, e_Adiv`. Qual combinação atinge $e_A < 10^(-3)$ com menos nós?

+ *E3 — Círculo unitário e ordem em $n$.*
  Disco unitário com `circulo` (`ndiv in {8,16,32,64}`, `tipo=1` e `tipo=2`).
  Estude $P$ e $A$ vs $2 pi$ e $pi$. Qual ordem aparente em $n$ (ou em $h ~ 1\/n$)?
  Lembrete: `tipo=1` é o $n$-ágono regular.

+ *E4 — Momento $I_x$.*
  #image("../assets/indo-para-2d/exercicio-1.png", width: 55%)
  #image("../assets/indo-para-2d/exercicio-2.png", width: 55%)
  Calcule $I_x = integral_Omega y^2 dif A$ nas duas figuras
  (`moment_Ix`, divergência com $upright(bold(F))=(0,y^3\/3)$) e compare com
  $I_x = (a^4)/96 (9 sqrt(3) - 2 pi)$ e
  $I_x = 2 · 10^4 pi - (20^2 pi)/2 (80\/(3 pi))^2 + ((20^2 pi)/2)(15 + 80\/(3 pi))^2$.
  Geradores: `triangulo_equilatero`, `triangulo_com_incirculo`, `setor60`, `semicirculo`.

+ *E5 — Furo e orientação.*
  `placa_com_furo`.
  (a) `geometric_props`: confira $A = A_"ret" - pi R^2$ (e o perímetro $P_"ret"+2 pi R$).
  (b) `reverse_hole=false` e `format2d(; orient=false)`: relate o sinal de $A$.
  (c) `plot_geo`: normais no contorno externo vs furo.
  Extra 3D: `cubo` + `format3d` (S = 6, V = 1) ou o toro de `scripts/toro.jl`.

== O que fica para Laplace 2D

- Convenção `"0;T"`, `"1;q"` (e elasticidade `tx;ux;ty;uy`).
- `attach_analytical!`, `H_G_*`, `solve`, `rel_error`.
- CDC em cantos e o papel do elemento descontínuo no conflito de nós.

== Leituras e código

- `codes/geom_gmsh/scripts/circulo.jl`, `scripts/toro.jl`
- `codes/geom_gmsh/src/Input.jl` — `format2d`, `format3d`, `discontinuous_nodes_weights`
- `codes/geom_gmsh/src/Interpolation.jl` — pesos baricêntricos, `shapefun`, `shapefun2D`
- `codes/geom_gmsh/src/GeometricProperties.jl` — `:radial`, `:divergence`, `moment_Ix`
- `codes/geom_gmsh/src/Meshes.jl` — `circulo`, `anel`, `placa_com_furo`, `toro`, `cubo`
- `codes/geom_gmsh/_solu/` — exercícios E1 a E5 resolvidos
- Manual Gmsh: #link("https://gmsh.info/doc/texinfo/gmsh.html")[documentação]
