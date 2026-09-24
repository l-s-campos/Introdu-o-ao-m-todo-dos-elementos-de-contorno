#import "_theme.typ": *
#show: bem-slides.with(
  title: [Poisson 2D],
  subtitle: [DIBEM · termo de domínio],
)

// Conteudo completo da aula (chapters/09-poisson-2d.typ)

= Poisson 2D

== Poisson 2D - intro

#set text(size: 12.5pt)
No capítulo *Laplace 2D* o problema era $nabla^2 T = 0$ e o sistema ficou

$ H T = G q $

só com integrais em $Gamma$. Agora admitimos *fonte de domínio*:

$ nabla^2 T = f(x) quad "em" Omega . $

A equação integral *ganha um termo em* $Omega$. Em vez de malha volumétrica de EF, o
`BEM_gmsh` monta o operador *DIBEM* (Dual Integration Boundary Element Method): uma matriz $M$ tal que

$
  bold(d) approx M bold(f) ,
  quad
  d_i approx integral_Omega T^* (x, x_i)\, f(x)\, dif Omega .
$

O produto $T^* f$ é aproximado por RBF mais um polinômio de primeira ordem nas
colocações (contorno + internos). As integrais das bases e de $T^*$ sozinha caem no
contorno — o mesmo espírito do `geometric_props` no cap. *Indo para 2D*.

#set text(size: 18pt)

== Objetivos

+ Derivar *por que* aparece $integral_Omega T^* f$ a partir do Laplace.
+ Entender DIBEM: produto $T^* f$ $arrow.r$ pesos $S$ $arrow.r$ $M$, diagonal por $I_s$.
+ Montar `H_G_full_direct` + `DIBEM` e usar o cache `dad.M`.
+ Resolver Poisson estacionário $H T - G q = M f$ (RHS de domínio).
+ Ver $M$ como “massa” no transiente (ponte).
+ Exercícios clássicos com API DIBEM e erros do apêndice.

#set text(size: 18pt)

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

#set text(size: 18pt)

== De Laplace a Poisson

== PDE e notação

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

#set text(size: 18pt)

== Identidade integral

#set text(size: 14pt)
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

#set text(size: 18pt)

== DIBEM em detalhe

== Ideia em uma frase

#set text(size: 16pt)
Aproxima-se o produto $T^*(x, y)\, f(y)$ por RBF mais um polinômio de primeira ordem.
A integral em $Omega$ fica $S$ vezes os valores nodais desse produto. Fora da diagonal,
$M_(i j) = S_j T^*(x_i, x_j)$; a diagonal sai de $integral T^*$, sem avaliar a solução
fundamental na singularidade.

$x$ é o ponto fonte; $y$ percorre $Omega$. Os centros das bases são as $N =$ `dad.nt`
colocações (contorno + internos).

#set text(size: 18pt)

== 1. Aproximação

#set text(size: 15pt)
$
  T^*(x, y) f(y)
  approx
  sum_(j=1)^N phi_j(y) alpha_j
  + c_0 + c_1 y_1 + c_2 y_2
  =
  Phi(y) alpha + P(y) c.
$

$P(y) = mat(1, y_1, y_2)$ e $c = mat(c_0; c_1; c_2)$.
$alpha$ pesa as bases; $c$ pesa o polinômio.
Default do pacote: `PHS()` (`poly_deg = 2`; o caso linear é `poly_deg = 1`).

#set text(size: 18pt)

== 2. Coeficientes

#set text(size: 15pt)
Nos pontos $x_i$, com $F_(i j) = phi_j(x_i)$:

$
  F alpha + P^T c = T^* f,
  quad
  P alpha = 0.
$

$
  A = mat(F, P^T; P, 0),
  quad
  mat(alpha; c) = A^(-1) mat(T^* f; 0).
$

A segunda equação do bloco é a reprodução polinomial.

#set text(size: 18pt)

== 3. Integração sobre o domínio

#set text(size: 14pt)
$
  integral_Omega T^* f dif y
  = I alpha + I_p c,
  quad
  I = integral_Omega Phi dif y,
  quad
  I_p = integral_Omega P dif y.
$

$
  mat(S, S_p) = mat(I, I_p) A^(-1),
  quad
  integral_Omega T^* f dif y = S (T^* f).
$

$S$ não depende de $f$ nem do ponto fonte. O bloco $S_p$ multiplica o zero do lado direito.

`IF` guarda $I$, `IP` guarda $I_p$, `ID[i]` guarda $I_s(x_i) = integral T^* dif y$.
Essas integrais caem em $Gamma$ pela primitiva vezes $upright(bold(n)) dot.op upright(bold(r)) \/ r^2$.
$I_s$ não usa a RBF. Longe do elemento: lumping; perto: Gauss (`_dibem_accumulate_IF_ID!`).

#set text(size: 18pt)

== 4. Matriz $M$

#set text(size: 16pt)
O mesmo $S$ vale em todo ponto fonte. Para $i != j$,

$
  M_(i j) = S_j T^*(x_i, x_j),
  quad
  (M f)_i = sum_(j=1)^N M_(i j) f_j.
$

$T^*(x_i, x_i)$ é singular: $M_(i i)$ não sai dessa fórmula.

No pacote, $S$ é `dibem_c` e $D_(i j) = T^*(x_i, x_j)$ ($i != j$). Fora da diagonal:
`M = D .* S'`.

#set text(size: 18pt)

== 5. Diagonal indireta

#set text(size: 15pt)
Com $f equiv 1$, a integral é $I_s$, sem RBF. A discretização reproduz esse caso:

$
  M bold(1) = bold(I_s),
  quad
  M_(i i) = (I_s)_i - sum_(j != i) M_(i j).
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
  DIBEM: $f equiv 1$ $arrow.r$ $M bold(1)=bold(I_s)$ (`ID` no código).
  Em ambos, a singularidade fica fora da quadratura.
]

$ bold(d) = M bold(f). $

#set text(size: 18pt)

== 6. Regularização constante

#set text(size: 14pt)
A diagonal indireta é a decomposição

$
  integral_Omega T^* f dif y
  =
  integral_Omega T^* (f - f_0) dif y
  + f_0 I_s.
$

Com $f_0 = f_i$, a diferença se anula em $y = x_i$. Só entram $j != i$:

$
  integral_Omega T^* f dif y
  approx
  sum_(j != i) M_(i j) (f_j - f_i) + f_i I_s.
$

O fator de $f_i$ é $M_(i i)$ do passo 5. A condição $M bold(1) = bold(I_s)$ é essa regularização em forma matricial.

#set text(size: 18pt)

== 7. Regularização linear

#set text(size: 13pt)
Com $r = y - x$ e $f_0 = f_i$,

$
  integral T^* f
  =
  integral T^* (f - f_0 - f_0' dot.op r)
  +
  integral T^* (f_0 + f_0' dot.op r).
$

$F alpha = f$ e $(F')_(i j) = nabla phi_j(x_i)$, logo $f_0' = F' (F backslash f)$ na linha $x_i$.
A parcela linear integra sem RBF:

$
  integral T^* (f_0 + f_0' dot.op r) dif y
  = f_0 I_s + f_0' dot.op I_r,
  quad
  I_r = integral T^* r dif y.
$

Em $y = x_i$, $r = 0$: o resto e a sua derivada se anulam. Para $i != j$, com $r_(i j) = x_j - x_i$,

$
  M_(i j) = S_j (T^*(x_i, x_j) - (F' F^(-1)) r_(i j)).
$

$(F' F^(-1)) r_(i j) = F' (F backslash r)$ é o mesmo operador de $f_0'$, na posição relativa.
$f_0 + f_0' dot.op r$ segue em $I_s$ e $I_r$. A diagonal continua a do passo 5.

#set text(size: 18pt)

== O que DIBEM *não* é

- Não é malha de volume MEF: internos são *centros de RBF* e sensores, não elementos de $Omega$.
- Não substitui $H,G$: só constrói $M$.
- Poucos internos, aglomerados ou fora de $Omega$ $arrow.r$ $F$ mal-condicionada e $M$ ruim.

#set text(size: 18pt)

== Código no pacote

#set text(size: 14pt)
`DIBEM_dense` monta a $M$ dos passos 4 e 5 (ação = regularização constante).
No código, `P` tem uma linha por nó: `K = [F P; P' 0]`.

```julia
# F_ij = φ_j(x_i)
# [S; λ] = [F P; P' 0] \ [IF; IP]   # S = dibem_c
# D_ij = T*(xi, xj), i ≠ j
# M = D .* S' ;  M_ii = ID[i] - Σ_{j≠i} M_ij

S = _dibem_poly_c(F, IF, pts, rbf; IP=IP)
M = D .* S'
for i in 1:dad.nt
    M[i, i] = 0
    M[i, i] = -sum(M[i, :]) + ID[i]
end
```

```julia
DIBEM(dad)                         # :dense, rbf=PHS()
DIBEM(dad; rbf=PHS(3; poly_deg=0))
DIBEM(dad; method=:hmatrix)        # N grande — extra / trabalhos
```

Exige `format2d(..., pontointerno=true)`. Então `dad.nt = dad.n + n_"int"`.

#set text(size: 18pt)

== Onde entra no sistema

== Estacionário (Poisson)

#set text(size: 15pt)
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

#set text(size: 18pt)

== Transiente / outros

#table(
  columns: (auto, auto),
  inset: 7pt,
  stroke: 0.5pt + luma(200),
  [*Modelo*], [*Papel de $M$*],
  [Calor $dot(T) ~ kappa nabla^2 T$], [`solve_transient` usa $M$ como massa],
  [Onda], [`solve_Houbolt` \/ `solve_transient_o2`],
  [Helmholtz $nabla^2 T + kappa^2 T = 0$], [desloca $H - kappa^2 M$ (`κ2` em `applyBC`)],
)

#set text(size: 18pt)

== Poisson manufaturado ()

#set text(size: 14pt)
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

#set text(size: 18pt)

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

#set text(size: 18pt)

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

#set text(size: 18pt)

== Exercícios

Apêndice de erros. Sempre reporte $N$ (`dad.n`, `dad.nt`) e a norma usada.

#set text(size: 18pt)

== E1 - Membrana triangular

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

#set text(size: 18pt)

== E2 - Elipse

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

#set text(size: 18pt)

== E3 - Transiente (placa / cubo)

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

#set text(size: 18pt)
