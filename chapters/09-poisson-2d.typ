// Poisson 2D — DIBEM (fonte de domínio)
// Pré-requisito: Laplace 2D

= Poisson 2D
<poisson-2d>

No capítulo *Laplace 2D* o problema era $nabla^2 T = 0$ e o sistema ficou

$ H T = G q $

só com integrais em $Gamma$. Agora admitimos *fonte de domínio*:

$ nabla^2 T = f(x) quad "em" Omega . $

A equação integral *ganha um termo em* $Omega$. Em vez de malha volumétrica de EF, o
`BEM_gmsh` monta o operador *DIBEM* (Direct Interpolation BEM): uma matriz $M$ tal que

$
  bold(d) approx M bold(f) ,
  quad
  d_i approx integral_Omega T^* (x, x_i)\, f(x)\, dif Omega ,
$

usando RBF nos nós de colocação (contorno + internos) e redução das integrais de
volume ao contorno — o mesmo espírito do `geometric_props` no cap. *Indo para 2D*.

== Objetivos

+ Derivar *por que* aparece $integral_Omega T^* f$ a partir do Laplace.
+ Entender DIBEM: RBF $arrow.r$ $F$ $arrow.r$ primitivas no contorno $arrow.r$ $M$.
+ Montar `H_G_full_direct` + `DIBEM` e usar o cache `dad.M`.
+ Resolver Poisson estacionário $H T - G q = M f$ (RHS de domínio).
+ Ver $M$ como “massa” no transiente (ponte).
+ Exercícios clássicos com API DIBEM e erros do apêndice.

== Mapa

+ De Laplace a Poisson (PDE $arrow.r$ BIE)
+ DIBEM em detalhe (ideia, fórmulas, código)
+ Onde $M$ entra no sistema
+ Lab A: sanidade ($f=0$ e tamanho de $M$)
+ Lab B: Poisson manufaturado $u = x^2+y^2$ ($f=4$)
+ Transiente (ponte)
+ Armadilhas
+ Exercícios
+ Leituras

== De Laplace a Poisson

=== PDE e notação

$
  nabla^2 T = f quad "em" Omega ,
  quad
  q = - k (partial T)\/(partial n) quad "em" Gamma .
$

- $f = 0$: Laplace (capítulo anterior) — $H T = G q$.
- $f$ conhecida: Poisson estacionário — $H T - G q = M f$.
- no tempo: $M$ multiplica $dot(T)$ ou $accent(T, dot.double)$ (calor\/onda).

Física típica (tabela de aplicações do Laplace): geração de calor, membrana
($S nabla^2 w = -p$), etc. Mudam os *nomes* de $T$ e $f$.

=== Identidade integral

Mesma identidade de Green do Laplace, peso = SF $T^*$ com $-nabla^2 T^* = delta(x-x_d)$:

$
  integral_Omega (T^* nabla^2 T - T nabla^2 T^*) dif Omega
  =
  integral_Gamma (T^* partial_n T - T partial_n T^*) dif s .
$

Com $nabla^2 T = f$ e o salto da SF:

$
  c(x_d) T(x_d)
  =
  integral_Gamma T q^* dif s
  -
  integral_Gamma T^* q dif s
  +
  integral_Omega T^* f dif Omega .
$

#table(
  columns: (auto, auto),
  inset: 7pt,
  stroke: 0.5pt + luma(200),
  [*Onde*], [*$c(x_d)$*],
  [Interior], [$1$],
  [Contorno liso], [$1\/2$],
  [Canto], [$c != 1\/2$; no código, diagonal de $H$ indireta],
)

$H$ e $G$ são *os mesmos* do Laplace. O bloco novo é só

$ d(x_d) := integral_Omega T^* (x, x_d)\, f(x)\, dif Omega . $

Discretamente, nas colocações $x_i$:

$
  H T - G q = bold(d) ,
  quad
  bold(d) approx M bold(f) .
$

== DIBEM em detalhe

=== Ideia em uma frase

Interpola $f$ (ou outra densidade de domínio) por RBF nos $N =$ `dad.nt` pontos
(contorno + internos) e transforma $integral_Omega T^* f$ em combinações de
*integrais só em* $Gamma$, montando a matriz $M$ uma vez.

=== Passo a passo

Sejam $x_1,...,x_N$ as colocações (`point(dad,i)`, $i=1..N$).

*1. Interpolação*

$
  f(x) approx sum_(j=1)^N phi.alt(|x - x_j|) alpha_j
$
(mais monômios de baixa ordem se `poly_deg` $>= 0$).


Nos nós:

$
  F bold(alpha) = bold(f) ,
  quad
  F_(i j) = phi.alt(|x_i-x_j|)
  \ (i != j),\
  F_(i i)\ "com regularização" .
$

No pacote o default é `PHS()` (polyharmonic spline; ver `Radial_Basis_Functions.jl`).

*2. Integral contra a SF*

$
  d(x_i)
  =
  integral_Omega T^* (x,x_i) f(x) dif Omega
  approx
  sum_j alpha_j integral_Omega T^* (x,x_i) phi.alt_j (x) dif Omega .
$

*3. Redução ao contorno*

Integrais radiais $integral_Omega g(r) dif Omega$ com $r=|x-x_i|$ viram integrais em
$Gamma$ via primitiva + fator $upright(bold(n))·upright(bold(r))\/r^2$ em 2D
(mesma geometria da área no cap. *Indo para 2D*).

O código acumula, por fonte $i$:

- `IF[i]` — contribuição de contorno ligada à RBF;
- `ID[i]` — contribuição ligada à SF (caso $f equiv 1$).

Longe do elemento: lumping nodal; perto: Gauss regular (`_dibem_accumulate_IF_ID!`).

*4. Matriz $M$*

Forma densa típica (RBF genérica em `DIBEM_dense`):

```julia
M = IF' / F .* D          # D_ij ~ T*(xi, xj), i≠j
for i in 1:dad.nt
    M[i, i] = 0
    M[i, i] = -sum(M[i, :]) + ID[i]
end
```

A diagonal impõe a identidade do caso constante:

$
  M bold(1) = bold(I D)
  quad arrow.r.double.long quad
  M_(i i) = I D_i - sum_(j != i) M_(i j) .
$

#block(
  width: 100%,
  fill: rgb("#ecfdf5"),
  inset: 10pt,
  radius: 4pt,
  stroke: 0.5pt + rgb("#99f6e4"),
)[
  *Analogia com $H$.*
  Laplace: $T equiv 1$, $q equiv 0$ $arrow.r$ $H bold(1)=0$ $arrow.r$
  $H_(i i)=-sum_(j!=i) H_(i j)$.
  DIBEM: $f equiv 1$ $arrow.r$ $M bold(1)=bold(I D)$ $arrow.r$ diagonal de $M$ por linha.
  Em ambos: *não integre o singular “na marra”*.
]

Depois disso:

$ bold(d) = M bold(f) . $

=== O que DIBEM *não* é

- Não é malha de volume MEF: internos são *centros de RBF* e sensores, não elementos de $Omega$.
- Não substitui $H,G$: só constrói $M$.
- Poucos internos, aglomerados ou fora de $Omega$ $arrow.r$ $F$ mal-condicionada e $M$ ruim.

=== Código no pacote

```julia
# Laplace/Domain.jl — essência DIBEM_dense
# F_ij = φ(|xi-xj|),  D_ij = T*(xi,xj)
# IF, ID = primitivas no contorno
# M = (IF'/F) .* D  +  diagonal via ID

function DIBEM_dense(dad::BEMdata{<:Laplace}; rbf=PHS())
    # ... monta F, D ...
    _dibem_accumulate_IF_ID!(IF, ID, dad, rbf; mon=mon, IP=IP)
    M = IF' / F .* D
    for i in 1:dad.nt
        M[i, i] = 0
        M[i, i] = -sum(M[i, :]) + ID[i]
    end
    set_cache!(dad; M, dibem_F=F, dibem_rbf=rbf, dibem_method=:dense)
    return M
end
```

API:

```julia
DIBEM(dad)                         # :dense, rbf=PHS()
DIBEM(dad; rbf=PHS(3; poly_deg=0))
DIBEM(dad; method=:hmatrix)        # N grande — extra / trabalhos
```

Exige colocações de domínio: `format2d(..., pontointerno=true)` (ou lista de internos no `dad`).
Então `dad.nt = dad.n + n_"int"`.

== Onde $M$ entra no sistema

=== Estacionário (Poisson)

Equação discreta *antes* das CDC:

$ H T - G q = M f . $

1. `H_G_full_direct(dad, npg)` — monta $H,G$ (como no Laplace).
2. `DIBEM(dad)` — monta $M$.
3. Amostra $f$ em *todas* as colocações $i=1..N$ $arrow.r$ vetor `fvec`.
4. `d = M * fvec`.
5. `applyBC(dad)` — reorganiza $H T - G q$ em $A x = b_"CDC"$ (igual ao Laplace).
6. Soma o domínio no RHS: `dad.b .+= d` (mesmo tamanho que as linhas de $H$, `dad.nt`).
7. Resolve $A x = b$ e espalha $T,q$ (`bem_linsolve` + `split_sol!`, ou equivalente).

#block(
  width: 100%,
  fill: luma(248),
  inset: 10pt,
  radius: 4pt,
  stroke: 0.5pt + luma(200),
)[
  *Por que somar $d$ em `b` depois do `applyBC`?*
  As CDC só movem colunas de $H$ e $G$. O termo $M f$ já está no *lado direito* da
  identidade integral; na forma $A x = b$ ele permanece como contribuição conhecida
  em todas as equações de colocação (contorno e internos).
]

=== Transiente \/ outros

#table(
  columns: (auto, auto),
  inset: 7pt,
  stroke: 0.5pt + luma(200),
  [*Modelo*], [*Papel de $M$*],
  [Calor $dot(T) ~ kappa nabla^2 T$], [`solve_transient` usa $M$ como massa],
  [Onda], [`solve_Houbolt` \/ `solve_transient_o2`],
  [Helmholtz $nabla^2 T + kappa^2 T = 0$], [desloca $H - kappa^2 M$ (`κ2` em `applyBC`)],
)



== Poisson manufaturado $u = x^2+y^2$ ($f = 4$)

No quadrado unitário, $u = x^2+y^2$ $arrow.r$ $nabla^2 u = 4$.
Dirichlet com o valor exato em todo o contorno.

```julia
using DrWatson
@quickactivate :BEM
using LinearAlgebra, Statistics
include(datadir("Laplace", "Laplace_dad.jl"))

ufun(p) = p[1]^2 + p[2]^2
fval = 4.0

msh = quadrado(ndiv=12, show=false, nome="pois_dibem")
dad = format2d(msh, Laplace(1.0); pontointerno=true)

for i in 1:dad.n
    dad.BC[i] = 0
    dad.BV[i] = ufun(dad.Nodes[i])
end

H_G_full_direct(dad, 16)
M = DIBEM(dad; rbf=PHS(3; poly_deg=1))

# f em todas as colocações (contorno + internos)
fvec = fill(fval, dad.nt)
# se f for campo: fvec = [f(point(dad, i)) for i in 1:dad.nt]
d = M * fvec

applyBC(dad)                 # A, b só com CDC (H,G)
dad.b .+= d                  # domínio

x = dad.A \ dad.b
Tfull = zeros(dad.nt)
qfull = zeros(dad.n)
Tfull[1:length(x)] .= x
split_sol!(dad, Tfull, qfull)
set_cache!(dad; T=Tfull, q=qfull)

pts = [dad.Nodes; dad.internalNodes]
uex = ufun.(pts)
u   = dad.T
rmse = sqrt(mean(abs2, u .- uex))
einf = maximum(abs, u .- uex)
@show rmse, einf
plot_geo(dad)
```

*Variantes (tabela):* `ndiv in {8,12,16}`; `PHS(3; poly_deg=1)` vs `PHS(5; poly_deg=2)`;
misto Dirichlet\/Neumann com $q = -2(x n_x + y n_y)$ nos lados horizontais.

Normas: *Apêndice: medidas de erro*. Aqui RMSE e $epsilon_infinity$ em `pts`.

#block(
  width: 100%,
  fill: luma(248),
  inset: 10pt,
  radius: 4pt,
  stroke: 0.5pt + luma(200),
)[
  *Implementação.* O caminho `applyBC` $arrow.r$ `b .+= M*f` $arrow.r$ `A\\b` $arrow.r$ `split_sol!`
  deixa o papel de $M$ *explícito*. Se a versão do pacote expuser um helper de Poisson
  estacionário com DIBEM, use-o — a matemática é a mesma: $H T - G q = M f$.
]

== Transiente (ponte)

Com $M$ no cache:

```julia
H_G_full_direct(dad, 16)
DIBEM(dad)
# sol = solve_transient(dad, 0.01, 1.0)     # calor
# solve_Houbolt(dad, 0.05, 2.0)             # 2ª ordem
```

Estabilidade Δt ↔ malha ↔ qualidade de $M$: monografia (proposta C dos trabalhos).
Nos exercícios E3–E4, o foco é sensor no tempo + erro, não a teoria completa de CFL.

== Armadilhas

#table(
  columns: (auto, auto),
  inset: 7pt,
  stroke: 0.5pt + luma(200),
  [*Sintoma*], [*Causa típica*],
  [`DIBEM` falha ou $M$ absurda], [Sem `pontointerno`; pontos ruins; RBF inadequada],
  [Erro grande com $f$ simples], [Malha grossa; `poly_deg` baixo; CDC errada],
  [Solução “desloca” com $f=0$], [Somou `d` não nulo por engano; `fvec` sujo],
  [Dimensão de `b` ≠ `d`], [Misturou `n` e `nt`; internos desligados],
  [Neumann manufaturado explode], [Sinal de $q=-k partial_n u$],
  [Transiente explode], [Δt grande; $M$ pobre],
)

== Exercícios

Apêndice de erros. Sempre reporte $N$ (`dad.n`, `dad.nt`) e a norma usada.

=== E1 — Membrana triangular

$ S nabla^2 w = -f $, $a=5$, $f=10$, $S=1$ (unidades do enunciado); contorno $w=0$.

#image("../assets/poisson-2d/membrana.png", width: 80%)

$
  w = -f/(2 S) [
    1/2 (x^2+y^2)
    - 1/(a sqrt(3)) (y^3 - 3 x^2 y)
    - a^2\/18
  ] .
$

No código a PDE é $nabla^2 T = f_"bem"$ com $T=w$: use $f_"bem" = -f\/S$.
DIBEM + RHS; erros em ≥ 100 internos; mapa de $w$.

=== E2 — Elipse

$ nabla^2 u = 4 - x^2 $ no domínio da figura.

$
  u = [
    1.6 - 1/246 (50 x^2 - 8 y^2 + 33.6)
  ] (x^2\/4 + y^2 - 1) .
$

#image("../assets/poisson-2d/elipse.png", width: 80%)

Amostra $f(x,y)=4-x^2$ em cada colocação para montar `fvec`.
Erros de $u$ (internos) e de $q$ no contorno ($q=-partial_n u$, normal exterior);
≥ 3 malhas.

=== E3 — Transiente (placa \/ cubo)

$T_0=0$; uma face em $T=1$; propriedades unitárias. Série:

$
  T(Y,t)
  =
  1 - (4\/pi) sum_(n=0)^infinity
  ((-1)^n)\/(2n+1)
  exp{ -((2n+1)^2 pi^2 kappa t)\/(4 L^2) }
  cos(((2n+1) pi Y)\/(2 L)) .
$

#image("../assets/poisson-2d/cubo-transiente.png", width: 80%)

`DIBEM` + `solve_transient` (ou Houbolt); sensores vs $t$; erro em 2–3 instantes.
