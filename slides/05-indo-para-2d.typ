#import "_theme.typ": *
#show: bem-slides.with(
  title: [Indo para 2D],
  subtitle: [Malha · N · J · n · propgeo · GeomGmsh],
)

// Conteúdo completo da aula (chapters/05-indo-para-2d.typ)

= Indo para 2D

== Indo para 2D

A partir desta aula o laboratório usa o pacote
#link("https://github.com/l-s-campos/BEM_gmsh")[`GeomGmsh`]
em `codes/geom_gmsh` (Julia + Gmsh; recorte de geometria do `BEM_gmsh`).

Foco: *geometria de contorno* — malha, elementos, jacobiano, normal e
integrais só em $Gamma$.

*Condições de contorno* (o que cada aresta “vale” em $T$ ou $q$) ficam para o capítulo *Laplace 2D*.

Os trechos Julia são *pedaços do próprio pacote* (`src/`) ou dos exemplos em `scripts/` e `_solu/`.

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

#set text(size: 16pt)

No BEM 1D o “contorno” eram dois pontos. Em 2D:

- $Omega subset RR^2$ tem fronteira $Gamma = partial Omega$ *unidimensional* (curva fechada, eventualmente com furos);
- a geometria de $Gamma$ é aproximada por *elementos de contorno* $Gamma_e$;
- campo e fluxo no contorno usam interpolação em $xi in [-1,1]$ sobre cada $Gamma_e$;
- integrais de área em $Omega$ que admitem redução a $Gamma$ (radial ou divergência) preparam as identidades do BEM.

#image("../assets/indo-para-2d/perimetro.png", width: 52%)

== Ideia-chave

#keybox[
  *Ideia-chave.* Malhar só $Gamma$ não elimina a necessidade de boa geometria:
  $N_k$, $J$ e $upright(bold(n))$ em cada elemento são os tijolos de `geometric_props` e, depois, de $H$ e $G$.
]

== Ambiente: GeomGmsh + Gmsh

Na pasta `codes/geom_gmsh`:

```bash
julia --project=. -e "using Pkg; Pkg.instantiate()"
julia --project=. scripts/circulo.jl
julia --project=. scripts/toro.jl
```

O pacote `Gmsh.jl` já é dependência. GUI: #link("https://gmsh.info/#Download")[gmsh.info].

Geradores (`circulo`, `anel`, `placa_com_furo`, `cubo`, `toro`, …) em `src/Meshes.jl`:
`gmsh.initialize()` → gerar → `datadir(nome * ".msh")` → `finalize()`.
Depois `format2d` / `format3d` e `geometric_props`.

== Por que Gmsh no BEM (só geometria)

1. Curvas, arcos, splines, furos e sólidos OCC (toro) com orientação controlável.
2. Malha de contorno (curvas) + superfície 3D; volume opcional (pontos internos / orientação).
3. Ordem geométrica alinhada ao grau do elemento (`ordem` / `tipo`).
4. O mesmo gerador alimenta Laplace, Poisson e elasticidade no `BEM_gmsh` completo.

*Nesta aula* os grupos físicos *nomeiam* pedaços de curva e a superfície do domínio.
A convenção `"0;T"` / `"1;q"` entra em *Laplace 2D*.

== Orientação (crítico)

- Contorno *externo*: *anti-horário* → normal *para fora* de $Omega$.
- Contorno de *furo*: *horário* → normal ainda saindo de $Omega$.
- $A < 0$ (ou $V < 0$ em 3D) denuncia orientação invertida.

`format2d(; orient=true)` (default) corrige o sinal de $upright(bold(n))$ com as
células 2D da malha de superfície. `placa_com_furo(; reverse_hole=true)` já
inverte o laço do furo no Gmsh.

== O que format2d faz (Input.jl)

#set text(size: 16pt)

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

Os *mesmos* `shapefun` e pesos entram depois em $H$ e $G$.

== Fluxo mínimo — scripts/circulo.jl

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

`format2d` monta `BEMdata`. *Nesta aula:* `plot_geo` e `geometric_props`, não `solve` / `H_G_*`.

== Disco unitário (plot_geo)

#image("../assets/indo-para-2d/circulo-geo.png", width: 62%)

== Disco: números (`tipo=2`, `ndiv=16`)

#table(
  columns: (auto, auto, auto),
  inset: 6pt,
  stroke: 0.5pt + luma(200),
  [*Grandeza*], [*Numérico*], [*Exato*],
  [perímetro], [$6.2827$], [$2 pi approx 6.2832$],
  [área], [$3.1411$], [$pi approx 3.1416$],
  [centróide], [$approx (0,0)$], [$(0,0)$],
)

A malha linear (`tipo=1`) é um $n$-ágono regular:
$P = 2 n sin(pi\/n) = 2 pi - O(n^(-2))$,
$A = (n\/2) sin(2 pi\/n) = pi - O(n^(-2))$.
Com `tipo=2` os arcos parabólicos aproximam o círculo e o erro cai bem mais rápido.

== Elemento de contorno descontínuo

No `GeomGmsh` o padrão é *elemento descontínuo*: graus de liberdade de campo
*não* são compartilhados nos vértices entre elementos vizinhos
(evita conflito de CDC em cantos — detalhe em Laplace 2D).

#image("../assets/indo-para-2d/elemento-parabolico.png", width: 52%)

== Mapa de referência

$ Gamma_e $ é a imagem de $xi in [-1,1]$:

$
  upright(bold(x))(xi) = sum_(k=1)^(p+1) N_k (xi)\, upright(bold(x))_k^((e)) .
$

- *Geometria* da aresta Gmsh: nós `Equispaced(p)` (vértices + meios, se ordem 2).
- *Campo* (colocação): nós `Legendre(p)` = Gauss–Legendre em $(-1,1)$ via
  `discontinuous_nodes_weights(p)`.

```julia
function discontinuous_nodes_weights(p::Integer)
    p >= 1 || error("degree must be ≥ 1, got $p")
    return gausslegendre(p + 1)
end
```

`tipo=1,2,3` em `format2d` = grau $p$ do campo.

== Lagrange baricêntrico (Interpolation.jl)

Pesos (nós $x_0$): $ w_i = 1 \/ product_(j != i) (x_i - x_j) . $

$
  N_i (xi) = (w_i \/ (xi - x_i)) \/ sum_j (w_j \/ (xi - x_j)) .
$

```julia
function shapefun(poly::AbstractPolynomial, x)
    L = interpolation_matrix(poly, x)   # N(x)
    L, L * poly.Dmat                    # N e dN/dξ
end
```

#image("../assets/indo-para-2d/funcoes-forma-d.png", width: 42%)

== Dmat (resumo)

$ u'(xi_k) = sum_i N_i'(xi_k) u_i quad arrow.r.double quad upright(bold(u))' = D upright(bold(u)) $

$
D_(k i) = (w_i \/ w_k) \/ (x_k - x_i) quad (k != i),
wide
D_(k k) = - sum_(i != k) D_(k i)
$

`poly.Dmat` é montada *uma vez* por `diff_matrix`. A diagonal impõe $sum_i N_i' = 0$.

```julia
L = interpolation_matrix(poly, x)   # N
dN = L * poly.Dmat                  # N_i'(ξ) = Σ_k N_k(ξ) D_ki
```

`Dmat` *não* é o jacobiano da malha.

== Uso didático (tipos do dad)

#set text(size: 16pt)

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

Partição da unidade (E1):

```julia
poly = Legendre(3)
N, dN = shapefun(poly, ξ)
@assert abs(sum(N) - 1) < 1e-11
@assert abs(sum(dN)) < 1e-9
```

== Equispaced(2) = Lagrange clássico

$
  N_1 = xi(xi-1)\/2 , quad N_2 = 1-xi^2 , quad N_3 = xi(xi+1)\/2 .
$

#image("../assets/indo-para-2d/funcoes-forma.png", width: 52%)

Usado na *geometria* da aresta Gmsh de ordem 2.

== No BEMdata

#keybox[
  Depois de `format2d`, `dad.element_type` é o polinômio do campo (`Legendre(p)`),
  `dad.elem_weight` são os pesos de Gauss da colocação, e cada `dad.elements[e]`
  guarda índices, jacobianos nos nós de colocação e comprimento.
  `dad.Nodes` / `dad.Normal` são vistas de `collocation`.
]

== Jacobiano, comprimento e normal

$
  upright(bold(x))' (xi) = sum_k N_k' (xi)\, upright(bold(X))_k ,
  quad
  J = |upright(bold(x))'| ,
  quad
  upright(bold(n)) = (y', -x') \/ J
$

```julia
dx = dNgeo * X
J  = norm.(dx)
normal = tan2normal.(dx ./ J)   # (dx,dy) → (dy, -dx)/|·|
L  = abs(dot(J, wi))            # ∫ J dξ ≈ Σ J_k w_k
```

Perímetro global: $P = sum_e L_e$ (`g.perimeter`).

#image("../assets/indo-para-2d/normal-1.png", width: 32%)
#image("../assets/indo-para-2d/normal-2.png", width: 32%)

== geometric_props — duas estratégias

```julia
g  = geometric_props(dad)                         # strategy=:radial
g2 = geometric_props(dad; strategy=:divergence)
g3 = geometric_props(dad; npg_radial=20)
g4 = geometric_props(dad; npg_boundary=16)
```

`npg_boundary=nothing` reutiliza `elem.Jacobian` $times$ `dad.elem_weight`.
Um inteiro reintegra cada elemento com `shapefun(dad.element_type, ξ)` e Gauss.

== Integração radial (`:radial`, default)

#set text(size: 16pt)

Pólo na *origem*. $upright(bold(r))=upright(bold(x))$, $r=|upright(bold(r))|$, $hat(upright(bold(r)))=upright(bold(r))\/r$.

$
  dif A = rho \, dif rho \, dif theta ,
  quad
  dif theta = (upright(bold(n)) · hat(upright(bold(r)))) (dif Gamma) \/ r .
$

#image("../assets/indo-para-2d/area-polar.png", width: 42%)
#image("../assets/indo-para-2d/area-contorno.png", width: 42%)

== Radial: F e o loop

$
  I = integral_Omega f \, dif A
  = integral_Gamma F(upright(bold(x))) thin (upright(bold(n)) · upright(bold(x))) / r^2 dif Gamma
$

$
  F = integral_0^r f(rho hat(upright(bold(r)))) \, rho \, dif rho .
$

```julia
# GeometricProperties.jl — _geom_accum_2d_radial
# Fa, Fx, Fy = ∫_0^r {1, x, y} ρ dρ     (_calc_F_2d)
# A  += Fa * nr / r * wJ
# Sx += Fx * nr / r * wJ
# Sy += Fy * nr / r * wJ
# centróide = (Sx, Sy) / A
```

== `_calc_F_2d`

```julia
function _calc_F_2d(r, theta, qsi, w)
    dro = r / 2
    Fa = Fx = Fy = 0.0
    for i in eachindex(qsi)
        ρ = r / 2 * (qsi[i] + 1)      # ρ ∈ [0, r]
        x, y = ρ * cos(theta), ρ * sin(theta)
        Fa += ρ * dro * w[i]          # ∫ ρ dρ  (f = 1)
        Fx += x * ρ * dro * w[i]
        Fy += y * ρ * dro * w[i]
    end
    return Fa, Fx, Fy
end
```

Para $f=1$: $F = r^2\/2$ $arrow.r$ $A = 1/2 integral_Gamma (n · x) dif Gamma$.
Com `npg_radial=12` a primitiva polinomial já é exata.

== Teorema da divergência (`:divergence`)

$
  integral_Omega nabla · upright(bold(F)) \, dif A
  = integral_Gamma upright(bold(F)) · upright(bold(n)) \, dif Gamma .
$

#table(
  columns: (auto, auto, auto),
  inset: 5pt,
  stroke: 0.5pt + luma(200),
  [*$upright(bold(F))$*], [*$nabla · upright(bold(F))$*], [*Resultado*],
  [$(x,0)$ ou $(0,y)$], [$1$], [$A = integral_Gamma x n_x$],
  [$(x,y)\/2$], [$1$], [$A = 1/2 integral_Gamma x · n$],
  [$(x^2\/2, 0)$], [$x$], [$integral_Omega x = integral_Gamma (x^2\/2) n_x$],
  [$(0, y^2\/2)$], [$y$], [$integral_Omega y = integral_Gamma (y^2\/2) n_y$],
  [$(0, y^3\/3)$], [$y^2$], [$I_x$ (`moment_Ix`)],
)

== Acumulador de divergência

```julia
P += wJ
A += 0.5 * dot(x, n) * wJ
Sx += 0.5 * x[1]^2 * n[1] * wJ
Sy += 0.5 * x[2]^2 * n[2] * wJ
```

`npg_radial` é ignorado. A outra fórmula ($A = integral_Gamma x n_x$) está em `area_divergence`:

```julia
function area_divergence(dad::BEMdata)
    A = 0.0
    w_el = dad.elem_weight
    for elem in dad.elements, k in eachindex(elem.index)
        i = elem.index[k]
        wJ = elem.Jacobian[k] * w_el[k]
        A += dad.Nodes[i][1] * dad.Normal[i][1] * wJ
    end
    return A
end
```

== Comparar no disco

```julia
g = geometric_props(dad)                          # :radial
gd = geometric_props(dad; strategy=:divergence)
Adiv = area_divergence(dad)
@show g.area gd.area Adiv
```

No `scripts/circulo.jl` as três áreas coincidem até o erro de quadratura de $Gamma$ ($approx 3.14110$).

#keybox[
  *Radial × divergência.* Para $f=1,x,y$ coincidem (primitiva polinomial + mesma quad. em $Gamma$).
  Radial generaliza $f$ via $F$. Divergência é imediata quando existe $upright(bold(F))$ polinomial simples (`moment_Ix`, etc.).
]

== Exemplo 2D: disco unitário

Arquivo `scripts/circulo.jl`. Quatro arcos, ordem 2, 16 elementos 1D.

```julia
msh = circulo(; ndiv=16, show=false, nome="ex_circulo", ordem=2)
dad = format2d(msh, Laplace(1.0); tipo=2, pontointerno=false)
g = geometric_props(dad)
plot_geo(dad; show_normals=true, title="unit disk tipo=2 ndiv=16")
```

#image("../assets/indo-para-2d/circulo-geo.png", width: 48%)

== Exemplo 3D: toro — scripts/toro.jl

A *mesma* API vale em 3D: gerador → `format3d` → `geometric_props`.

`toro` chama `gmsh.model.occ.addTorus` (raio maior $R$, tubo $r$, eixo $z$, origem).

$
  S = 4 pi^2 R r , quad V = 2 pi^2 R r^2 , quad upright(bold(c)) = upright(bold(0)) .
$

```julia
R, r = 2.0, 0.5
msh = toro(; nome="ex_toro", R=R, r=r, lc=0.25, show=false)
dad = format3d(msh, Laplace(1.0); tipo=1, pontointerno=false)
g  = geometric_props(dad)                         # :radial
gd = geometric_props(dad; strategy=:divergence)
plot_geo(dad; show_normals=false, show_nodes=false)
```

== Toro (malha de superfície)

#image("../assets/indo-para-2d/toro-geo.png", width: 72%)

== Toro: números (`R=2`, $r=1\/2$, `lc=0.25`)

#table(
  columns: (auto, auto, auto),
  inset: 6pt,
  stroke: 0.5pt + luma(200),
  [*Grandeza*], [*Numérico*], [*Exato*],
  [área de superfície], [$39.18$], [$4 pi^2 R r approx 39.48$],
  [volume], [$9.57$], [$2 pi^2 R r^2 approx 9.87$],
  [centróide], [$approx (0,0,0)$], [$(0,0,0)$],
)

Radial = divergência até roundoff.
3D: $V = (1\/3) integral_Gamma x · n$; radial usa `_calc_F_3d` e
$j_"fac" = (n · x)\/R^3$.

== format3d (essência)

$
  upright(bold(n)) parallel partial_xi upright(bold(x)) times partial_eta upright(bold(x)) ,
  quad
  J = | partial_xi upright(bold(x)) times partial_eta upright(bold(x)) | .
$

```julia
N, dNx, dNy = shapefun2D(Equispaced(1), qsi)
nodes[idx]  .= N * X
dx1, dx2     = dNx * X, dNy * X
J            = norm.(cross.(dx1, dx2))
normal[idx]  = cross.(dx1, dx2) ./ J
```

Faces Gmsh: só quads (tipo 3 linear ou tipo 10 quadrático).
`cubo` + `format3d`: $S=6$, $V=1$, $c=(1\/2,1\/2,1\/2)$.

== Boas práticas de malha

1. Contorno externo anti-horário; furos horários.
2. Elementos descontínuos no campo (padrão `format2d`).
3. `ordem` da malha Gmsh alinhada a `tipo` em `format2d`.
4. Refino perto de cantos e furos (`placa_com_furo`).
5. Nesta aula: `pontointerno=false` basta para prop. de $Gamma$.
6. Em 3D, volume no `.msh` ajuda a orientar $upright(bold(n))$ (`toro` tenta `generate(3)`).

== Exercícios — entrega

Scripts no ambiente `GeomGmsh` (`julia --project=codes/geom_gmsh`),
`plot_geo` quando fizer sentido, erro vs analítico e 2–3 frases.

Fontes: `src/{Input,Interpolation,GeometricProperties,Meshes}.jl`,
`scripts/circulo.jl`, `scripts/toro.jl`. Soluções em `_solu/`.

== E1 — shapefun e partição da unidade

Para $p in {1,2,3,4}$:

(a) `poly = Legendre(p)`: em 21 abscissas, $sum_i N_i = 1$ e $sum_i N_i' = 0$.

(b) Plote as colunas de $N(xi)$ para $p=2$ e $p=3$.

(c) `Equispaced(2)` em $xi in {-1,0,1}$ vs Lagrange clássico.

(d) `nodes(Legendre(p))` vs `discontinuous_nodes_weights(p)[1]` — são o quê?
(`weights(poly)` são baricêntricos; o segundo retorno é a quad. de Gauss.)

== E2 — coroa circular (`anel`)

$1 <= r <= 2$. Exato: $P = 6 pi$, $A = 3 pi$, centróide na origem.

(a) `ndiv in {8,16,32}`, `tipo in {1,2}`: erros de $P$ e $A$.

(b) No *mesmo* `dad`: `:divergence`, `area_divergence` vs `:radial`.

(c) `npg_boundary=16` e `npg_radial=20` vs default.

(d) `plot_geo`: normais no furo para fora de $Omega$?

(e) Tabela `ndiv, tipo, n_col, e_P, e_A, e_Adiv`. Quem atinge $e_A < 10^(-3)$ com menos nós?

== E3 — círculo unitário, ordem em $n$

`circulo` com `ndiv in {8,16,32,64}`, `tipo=1` e `tipo=2`.

$P$ e $A$ vs $2 pi$ e $pi$. Qual ordem aparente em $n$ (ou $h ~ 1\/n$)?

Lembrete: `tipo=1` é o $n$-ágono regular.

== E4 — momento $I_x$

#image("../assets/indo-para-2d/exercicio-1.png", width: 38%)
#image("../assets/indo-para-2d/exercicio-2.png", width: 38%)

$I_x = integral_Omega y^2 dif A$ via `moment_Ix` ($F=(0,y^3\/3)$).

Fórmulas da lista; geradores: `triangulo_equilatero`, `triangulo_com_incirculo`, `setor60`, `semicirculo`.

== E5 — furo e orientação

`placa_com_furo`.

(a) $A = A_"ret" - pi R^2$, $P = P_"ret" + 2 pi R$.

(b) `reverse_hole=false` e `format2d(; orient=false)`: sinal de $A$.

(c) `plot_geo`: normais no externo vs furo.

Extra 3D: `cubo` ($S=6$, $V=1$) ou `scripts/toro.jl`.

== O que fica para Laplace 2D

- Convenção `"0;T"`, `"1;q"` (e elasticidade `tx;ux;ty;uy`).
- `attach_analytical!`, `H_G_*`, `solve`, `rel_error`.
- CDC em cantos e o papel do elemento descontínuo.

== Leituras e código

- `codes/geom_gmsh/scripts/circulo.jl`, `scripts/toro.jl`
- `src/Input.jl` — `format2d`, `format3d`, `discontinuous_nodes_weights`
- `src/Interpolation.jl` — `shapefun`, `shapefun2D`
- `src/GeometricProperties.jl` — `:radial`, `:divergence`, `moment_Ix`
- `src/Meshes.jl` — `circulo`, `anel`, `placa_com_furo`, `toro`, `cubo`
- `_solu/` — E1 a E5 resolvidos
