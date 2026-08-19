// Poisson 2D — RASCUNHO (não incluído em main.typ / site até aprovação)
// Fonte de domínio: DIBEM + particular/DRM no BEM_gmsh
// Pré-requisito: Laplace 2D (H, G, CDC, solve)

= Poisson 2D
<poisson-2d>

No capítulo *Laplace 2D* o problema era $nabla^2 T = 0$ e o sistema ficou

$ H T = G q $

só com integrais em $Gamma$. Agora admitimos *fonte de domínio*:

$ nabla^2 T = f(x) quad "em" Omega . $

A equação integral *ganha um termo em* $Omega$. O BEM clássico não quer malha de
volume tipo MEF: o `BEM_gmsh` oferece duas famílias de tratamento —

1. *DIBEM* — monta um operador $M$ tal que a contribuição de domínio é $M bold(beta)$
   (fonte estacionária, massa em transiente, etc.);
2. *Particular + BEM homogêneo* (`solve_poisson_rbf_bem!`) — aproxima uma particular
   $u_p$ com $nabla^2 u_p approx f$, resolve Laplace em $u_h = u - u_p$, soma de volta.

Os dois usam RBF e pontos internos; mudam *onde* a interpolação entra.

== Objetivos

+ Derivar *por que* aparece $integral_Omega T^* f$ a partir do Laplace.
+ Distinguir o papel de $M$ (DIBEM) e o caminho particular\/DRM.
+ Montar `H_G_full_direct` + `DIBEM` e entender o cache `dad.M`.
+ Resolver um Poisson manufaturado com `solve_poisson_rbf_bem!` e medir erro.
+ Ler os exercícios clássicos (membrana, elipse, transiente) com API explícita.

== Mapa

+ De Laplace a Poisson (PDE $arrow.r$ BIE)
+ Duas estratégias no pacote (quadro)
+ DIBEM em detalhe (ideias + fórmulas + código)
+ Particular \/ DRM (`solve_poisson_rbf_bem!`)
+ Lab A: $M$ e sanidade ($f=0$)
+ Lab B: Poisson manufaturado $u = x^2+y^2$
+ Transiente (ponte)
+ Armadilhas
+ Exercícios em escada
+ Leituras

== De Laplace a Poisson

=== PDE e notação

$ nabla^2 T = f quad "em" Omega,
  quad
  q = - k (partial T)\/(partial n) quad "em" Gamma . $

- $f = 0$: Laplace (capítulo anterior).
- $f$ conhecida: Poisson estacionário (este capítulo).
- $f$ vira operador no tempo ($dot(T)$, $accent(T, dot.double)$): calor\/onda — DIBEM como “massa”.

Física típica (lembre a tabela de aplicações do Laplace): geração de calor, carga em membrana
($S nabla^2 w = -p$), etc. Só mudam os nomes de $T$ e $f$.

=== Identidade integral

Parte-se da *mesma* identidade de Green do Laplace, com peso = SF $T^*$
( $-nabla^2 T^* = delta(x - x_d)$ ):

$
integral_Omega (
  T^* nabla^2 T - T nabla^2 T^*
) dif Omega
=
integral_Gamma (
  T^* partial_n T - T partial_n T^*
) dif s .
$

Com $nabla^2 T = f$ e o salto da SF em $x_d$:

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
  [Interior $x_d in Omega$], [$1$],
  [Contorno liso], [$1\/2$],
  [Canto], [$c != 1\/2$ (ângulo sólido); diagonal de $H$ indireta no código],
)

Tudo que já valia para $H$, $G$, CDC e `solve` *continua*.
O único bloco novo é aproximar

$ d(x_d) := integral_Omega T^* (x, x_d)\, f(x)\, dif Omega $

sem integrar em volume “na marra”.

=== Sistema discreto (visão)

Com as mesmas colocações do Laplace,

$ H T - G q = bold(d) , $

onde $bold(d)_i approx d(x_i)$. As estratégias diferem em *como* obter $bold(d)$
(ou em *eliminar* $bold(d)$ via particular).

== Duas estratégias no `BEM_gmsh`

#table(
  columns: (1.1fr, 1.2fr, 1.2fr),
  inset: 7pt,
  stroke: 0.5pt + luma(200),
  [*Ideia*], [*DIBEM*], [*Particular + BEM (DRM\/RBF)*],
  [Objeto central], [Operador $M$: $bold(d) approx M bold(f)$], [$u = u_h + u_p$, $nabla^2 u_h = 0$],
  [API], [`DIBEM(dad)` $arrow.r$ `dad.M`], [`solve_poisson_rbf_bem!(dad, f)`],
  [Arquivo], [`Laplace/Domain.jl`, `Domain_fast.jl`], [`Laplace/ParticularSolution.jl`],
  [Uso forte], [Transiente (massa); Helmholtz com $kappa^2$; fontes genéricas], [Poisson estacionário didático],
  [Pontos internos], [Centros RBF + linhas de $M$], [Centros do fit de $f$ e de $u_p$],
)

#block(
  width: 100%,
  fill: luma(248),
  inset: 10pt,
  radius: 4pt,
  stroke: 0.5pt + luma(200),
)[
  *No curso.* Para *entender* o termo de domínio e o cache $M$, use DIBEM (Labs A).
  Para *fechar um Poisson estacionário* com poucas linhas e analítico manufaturado,
  use `solve_poisson_rbf_bem!` (Lab B). Os exercícios “clássicos” (membrana, elipse)
  podem seguir o Lab B; transiente puxa de volta o $M$ do DIBEM (ou DRM temporal).
]

== DIBEM em detalhe

*DIBEM* = *Direct Interpolation Boundary Element Method*:
interpola o campo de domínio (fonte $f$, ou $dot(T)$, …) por RBF nos nós de colocação
(contorno + internos) e reduz as integrais de volume a *integrais de contorno*
via primitivas radiais — o mesmo espírito do `geometric_props` no cap. *Indo para 2D*.

=== Passo a passo conceitual

Sejam $x_1,...,x_N$ os pontos de colocação ($N =$ `dad.nt` = nós de contorno + internos).

*1. Interpolação da densidade de domínio* $beta$ (no Poisson estacionário, $beta = f$):

$
beta(x) approx sum_(j=1)^N phi.alt_j (x)\, alpha_j ,
quad
phi.alt_j (x) = phi.alt(|x - x_j|)
$

(mais polinômios de baixa ordem se a RBF pedir — `poly_deg` no pacote).

Nos nós:

$ F bold(alpha) = bold(beta) ,
  quad
  F_(i j) = phi.alt(|x_i - x_j|)
  quad (i != j),\ 
  F_(i i)\ "regularizado" .
$

*2. Integral contra a SF.* Queremos

$
d(x_i) = integral_Omega T^* (x, x_i)\, beta(x)\, dif Omega
approx sum_j alpha_j integral_Omega T^* (x, x_i)\, phi.alt_j (x)\, dif Omega .
$

*3. Redução ao contorno.* Integrais do tipo $integral_Omega g(r) dif Omega$ com $r = |x - x_i|$
viram integrais em $Gamma$ usando divergência \/ primitiva radial
(fator $upright(bold(n)) · upright(bold(r)) \/ r^2$ em 2D — igual à ideia de área).

O código acumula, por fonte $i$:

- `IF[i]` — primitiva radial da RBF $phi.alt$ “vista” do contorno;
- `ID[i]` — primitiva radial da própria SF $T^*$ (caso $beta equiv 1$).

*4. Montagem de $M$.* Forma densa típica no pacote (RBF genérica):

$
M approx (upright(bold(I F))^top F^(-1)) compose D
$

com $D_(i j) ~ T^* (x_i, x_j)$ fora da diagonal, e depois *correção de diagonal*
para que $M bold(1) = upright(bold(I D))$ (identidade do caso constante — análogo à
soma nula de $H$):

$
M_(i i) = I D_i - sum_(j != i) M_(i j) .
$

Assim, a contribuição de domínio na colocação $i$ é

$ d_i approx (M bold(f))_i . $

#block(
  width: 100%,
  fill: rgb("#ecfdf5"),
  inset: 10pt,
  radius: 4pt,
  stroke: 0.5pt + rgb("#99f6e4"),
)[
  *Analogia com a diagonal de $H$.*
  Em Laplace, $T equiv 1$, $q equiv 0$ $arrow.r$ $H bold(1) = 0$ $arrow.r$ $H_(i i) = -sum_(j != i) H_(i j)$.
  Em DIBEM, $f equiv 1$ $arrow.r$ $M bold(1) = bold(I)_1$ (integrais de $T^*$ no domínio)
  $arrow.r$ diagonal de $M$ por linha. Mesma filosofia: *não integre o singular “na marra”*.
]

=== O que *não* é DIBEM

- *Não* é malha volumétrica de elementos finitos: os internos não carregam “rigidez” de volume.
- *Não* substitui $H$ e $G$: só constrói $M$ (ou alimenta marchadores que usam $M$).
- Pontos internos ruins (aglomerados, fora de $Omega$, poucos) $arrow.r$ $F$ mal-condicionada e $M$ ruim.

=== Snippet real (`Domain.jl`)

```julia
# DIBEM_dense — essência (Laplace/Domain.jl)
# 1) F_ij = φ(|xi-xj|),  D_ij = T*(xi,xj)   (i≠j)
# 2) IF, ID por integração no contorno (primitivas radiais)
# 3) M = (IF' / F) .* D
# 4) diagonal: M[i,i] = -sum(M[i,:]) + ID[i]

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

API de alto nível:

```julia
DIBEM(dad)                              # default :dense, rbf=PHS()
DIBEM(dad; rbf=PHS(3; poly_deg=0))
DIBEM(dad; method=:hmatrix)             # grandes N — extra / trabalhos
```

`dad.nt` inclui internos se `format2d(..., pontointerno=true)`.

=== Onde $M$ entra depois

#table(
  columns: (auto, auto),
  inset: 7pt,
  stroke: 0.5pt + luma(200),
  [*Problema*], [*Papel de $M$*],
  [Poisson estacionário “manual”], [$H T - G q = M f$ (após CDC, $M f$ no RHS)],
  [Calor $dot(T) = kappa nabla^2 T$], [Massa: marchadores `solve_transient`, …],
  [Onda], [`solve_Houbolt` \/ `solve_transient_o2` com $M accent(T, dot.double)$],
  [Helmholtz $nabla^2 T + kappa^2 T = 0$], [Deslocamento $H - kappa^2 M$ (via `κ2` em `applyBC`)],
)

Para o aluno: depois de `DIBEM(dad)`, `has_cache(dad, :M)` é verdadeiro e `dad.M` é o operador.

== Particular + BEM homogêneo (DRM \/ RBF)

Ideia clássica de *solução particular*:

$
T = T_h + T_p ,
quad
nabla^2 T_p approx f ,
quad
nabla^2 T_h = 0 .
$

1. Ajusta-se $T_p$ (RBF global\/local sobre contorno+internos) com $nabla^2 T_p approx f$.
2. Transfere-se a particular para as CDCs:

$
T_h = T - T_p quad "(Dirichlet)" ,
quad
q_h = q - q_p quad "(Neumann)" .
$

3. Resolve-se *Laplace* em $T_h$ com o pipeline já conhecido (`H_G_full_direct` + `solve`).
4. Recupera-se $T = T_h + T_p$ (e $q$ análogo).

No pacote isso está empacotado:

```julia
# ParticularSolution.jl — essência
function solve_poisson_rbf_bem!(dad, f; method=:global,
        basis=PHS(3; poly_deg=1), npg=16)
    # fit f → particular up, qp
    # BV ← BV - (Tp ou qp) conforme BC
    H_G_full_direct(dad, npg)
    solve(dad)                 # homogêneo
    # T ← Th + Tp ;  q ← qh + qp
    return dad.T
end
```

Forma DRM equivalente (matrizes): `build_drm_matrices` monta $M$ tal que

$ H u - G q = M b quad "quando" nabla^2 u = b , $

com $M = (H Psi - G eta) F^(-1)$ (particular $Psi$ das RBF). É o primo “dual reciprocity”
do DIBEM; no curso, trate-os como *duas implementações da mesma necessidade*
(representar domínio só com contorno + pontos).

== Lab A — Sanidade: $M$ existe e $f = 0$ não estraga Laplace

```julia
using DrWatson
@quickactivate :BEM
using LinearAlgebra
include(datadir("Laplace", "Laplace_dad.jl"))

props = Laplace(1.0)
msh   = quadrado(ndiv=16, show=false)
dad   = format2d(msh, props; pontointerno=true)

attach_analytical!(dad, ana_laplace_linear(; direction=SA[1.0, 0.0]))
H_G_full_direct(dad, 16)
@show rel_error(dad)   # ainda sem solve? → solve primeiro
solve(dad)
e0 = rel_error(dad)

DIBEM(dad)             # monta M; não altera H,G já prontos
@show size(dad.M), dad.nt, dad.n
# com f=0 a contribuição de domínio é nula: o solve de Laplace permanece
solve(dad)
@show rel_error(dad), e0
plot_geo(dad)
```

*Esperado:* `dad.nt > dad.n` (há internos); `rel_error` do patch $T=x$ continua pequeno.
Se `DIBEM` falhar: confira `pontointerno=true` e malha válida.

== Lab B — Poisson manufaturado $u = x^2 + y^2$ ($f = 4$)

Solução exata no quadrado unitário: $u = x^2 + y^2$, logo $nabla^2 u = 4$.
CDC Dirichlet em todo o contorno com o valor exato (como no teste do repo).

```julia
using DrWatson
@quickactivate :BEM
using LinearAlgebra, Statistics
include(datadir("Laplace", "Laplace_dad.jl"))

ufun(p) = p[1]^2 + p[2]^2

msh = quadrado(ndiv=12, show=false, nome="pois_xy2")
dad = format2d(msh, Laplace(1.0); pontointerno=true)

# Dirichlet exato em todos os nós de contorno
for i in 1:dad.n
    dad.BC[i] = 0
    dad.BV[i] = ufun(dad.Nodes[i])
end

u = solve_poisson_rbf_bem!(dad, 4.0;
        method=:global, basis=PHS(3; poly_deg=1), npg=14)

pts = [dad.Nodes; dad.internalNodes]
uex = ufun.(pts)
rmse = sqrt(mean(abs2, u .- uex))
einf = maximum(abs, u .- uex)
@show rmse, einf, maximum(abs, u)

# opcional: comparar métodos de fit
# res = compare_poisson_rbf_bem(dad, 4.0, ufun; methods=(:global, :local), npg=12)
```

*Variantes (faça tabela):*
- `ndiv in {6,8,12,16}`;
- `basis = PHS(3; poly_deg=1)` vs `PHS(5; poly_deg=2)`;
- misto: lados verticais Dirichlet, horizontais Neumann com
  $q = - partial_n u = -2 (x n_x + y n_y)$ (ver `test/test_poisson_drm.jl`).

Normas: *Apêndice: medidas de erro*. Aqui RMSE e $epsilon_infinity$ nos pontos `pts` bastam;
`rel_error(dad)` só após `attach_analytical!` compatível — no manufaturado, calcule direto
como acima.

== Transiente (ponte, não é o núcleo da aula)

Com operador de domínio $M$:

#table(
  columns: (auto, auto, auto),
  inset: 7pt,
  stroke: 0.5pt + luma(200),
  [*Modelo*], [*API típica*], [*Notas*],
  [Calor (1ª ordem)], [`DIBEM` + `solve_transient`], [Δt e malha acoplados],
  [Onda \/ 2ª ordem], [`solve_Houbolt` \/ `solve_transient_o2`], [Proposta C dos trabalhos],
  [DRM no tempo], [`solve_transient_drm!`], [ParticularSolution.jl],
)

```julia
H_G_full_direct(dad, 16)
DIBEM(dad)
# sol = solve_transient(dad, 0.01, 1.0)
# solve_Houbolt(dad, 0.05, 2.0)
```

Exercícios de série de Fourier\/Bessel no final treinam *sensor no tempo* + erro
(apêndice), não a teoria completa de estabilidade — isso é monografia (proposta C).

== Armadilhas frequentes

#table(
  columns: (auto, auto),
  inset: 7pt,
  stroke: 0.5pt + luma(200),
  [*Sintoma*], [*Causa típica*],
  [`DIBEM` \/ fit falha ou $M$ absurda], [Sem internos; pontos fora de $Omega$; RBF inadequada],
  [Erro grande com $f$ polinomial simples], [`poly_deg` baixo demais; malha grossa; CDC errada],
  [Patch Laplace piora depois de `DIBEM`], [Reusou `solve` com RHS de domínio não nulo sem querer],
  [Neumann manufaturado explode], [Esqueceu $q = -k partial_n u$ (sinal \/ $k$)],
  [Internos bons, contorno ruim (ou vice-versa)], [CDC vs particular mal transferida; meça os dois],
  [Transiente instável], [Δt grande; $M$ pobre; ver proposta C],
)

== Exercícios

Use o *Apêndice: medidas de erro*. Sempre varie malha (`ndiv` \/ ordem) e reporte $N$.

=== E0 — Labs A e B

Entregue: (i) `size(M)`, `nt`, `n` e `rel_error` do patch $T=x$ após DIBEM;
(ii) tabela `ndiv` × RMSE\/$epsilon_infinity$ para $u=x^2+y^2$;
(iii) *uma* frase sobre o efeito de `poly_deg` ou do método `:global` vs `:local`.

=== E1 — Membrana triangular

Contorno fixo ($w=0$), carga uniforme e tensão constante:

$ S nabla^2 w = - f ,
  quad a = 5.0,\ f = 10,\ S = 1
  quad "(unidades do enunciado)" . $

#image("../assets/poisson-2d/membrana.png", width: 80%)

Analítico:

$
w = - f/(2 S) [
  1/2 (x^2 + y^2)
  - 1/(a sqrt(3)) (y^3 - 3 x^2 y)
  - a^2 \/ 18
] .
$

*Fazer:* malha triangular (Gmsh) com Dirichlet nulo na borda; resolva com
`solve_poisson_rbf_bem!` (ou DIBEM + RHS, se for o fluxo que a equipe adotar);
erro médio e $L_2$ em ≥ 100 internos; mapa de $w$.
Atenção: a PDE do código é $nabla^2 T = f_"bem"$ — ajuste o sinal\/escala
($f_"bem" = -f\/S$ se $T=w$).

=== E2 — Elipse com fonte polinomial

$ nabla^2 u = 4 - x^2 $

no domínio elíptico da figura. Solução analítica (material do curso):

$
u = [
  1.6 - 1/246 (50 x^2 - 8 y^2 + 33.6)
] (
  x^2 \/ 4 + y^2 - 1
) .
$

Fluxo no contorno (com $k=1$, $q = - partial_n u$): derive de $u$ ou use a expressão
legada do capítulo oficial — *confira o sinal* com a normal exterior.

#image("../assets/poisson-2d/elipse.png", width: 80%)

*Fazer:* erros de potencial nos internos e de fluxo no contorno; ≥ 3 malhas.

=== E3 — Transiente (placa \/ cubo unitário)

Temperatura inicial nula; uma face sobe e permanece em $T=1$; propriedades unitárias.
Série (1D efetiva na coordenada normal à face aquecida):

$
T(Y,t)
=
1 - 4/pi sum_(n=0)^infinity
  ((-1)^n)/(2n+1)
  exp{ - ((2n+1)^2 pi^2 kappa t)/(4 L^2) }
  cos( ((2n+1) pi Y)/(2 L) ) .
$

#image("../assets/poisson-2d/cubo-transiente.png", width: 80%)

*Fazer:* `DIBEM` + `solve_transient` (ou Houbolt); sensores em $Y$ fixos vs $t$;
tabela de erro em 2–3 instantes; comente Δt.

=== E4 — Extra: cilindro oco transiente

$a=1$, $b=2$, $T(a,t)=1$, simetria; série com $J_0,Y_0$ e raízes de

$ J_0 (a x) Y_0 (b x) - J_0 (b x) Y_0 (a x) = 0 $

(`Roots.jl` \/ `fzeros`). DIBEM + marchador; erro em $r$ médios.

== O que fica para depois

- H-matriz \/ FMM em `DIBEM(; method=:hmatrix)` — custo (trabalhos, proposta D).
- Multirregião com fonte por subdomínio — proposta B.
- Elasticidade com força de corpo — DIBEM elástico no pacote (fora do núcleo 30 h).

== Leituras e código

- Cap. *Laplace 2D* — $H,G$, CDC, diagonal, `solve`
- Cap. *Indo para 2D* — integração radial (mesma geometria da redução de domínio)
- *Apêndice: medidas de erro*
- `src/Laplace/Domain.jl`, `Domain_fast.jl` — `DIBEM`
- `src/Laplace/ParticularSolution.jl` — `solve_poisson_rbf_bem!`, `build_drm_matrices`
- `test/test_poisson_drm.jl`, `scripts/poisson_rbf_bem_compare.jl`
- `src/Laplace/Solver.jl` — transiente
- #link("https://youtu.be/uSREar_ejnM")[gravação]
