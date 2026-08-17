// Elasticidade 2D — RASCUNHO (não incluído em main até aprovação)

= Elasticidade 2D
<elasticidade-2d>

Mesmo pipeline do Laplace: malha de contorno, SF, $H$\/$G$ (agora *blocos* 2×2 por nó),
CDC, `solve`. A primária deixa de ser o escalar $T$ e passa a ser o deslocamento
$upright(bold(u)) = (u_x, u_y)$; o conjugado de Neumann é a tração $upright(bold(t))$.

Notação (glossário): Poisson do *material* $nu$; função peso da formulação $v$ — não misturar.

== Objetivos

+ Ligar equilíbrio + Hooke à BIE de elasticidade (Kelvin).
+ Rodar patch de Dirichlet e um problema clássico no `BEM_gmsh`.
+ Ler CDC `"tx;ux;ty;uy"`.
+ Reportar erro em $upright(bold(u))$ (e tensões, se o exercício pedir).

== Mapa

+ Problema de contorno elástico 2D
+ SF de Kelvin (o que muda vs Laplace)
+ Código: patch + tubo
+ Exercícios (cilindro, Kirsch, viga)

== Problema contínuo (resumo)

Equilíbrio: $partial_j sigma_(i j) + f_i = 0$.

Cinemática: $epsilon_(i j) = (partial_i u_j + partial_j u_i)\/2$.

Hooke isotrópico (com $E$, $nu$, $mu = E\/(2(1+nu))$).

Eliminando tensões chega-se ao sistema de Navier em $upright(bold(u))$ — a SF de Kelvin é a resposta a força pontual nesse operador.

No contorno, para cada nó e direção, prescreve-se *ou* deslocamento *ou* tração
(analogamente a $T$ vs $q$).

== Equação integral (forma operacional)

Com núcleos tensoriais $U_(i j)$ (deslocamento fundamental) e $T_(i j)$ (tração fundamental),

$
c_(i k) u_k (x_d)
  + integral_Gamma T_(i j) (x_d, x) u_j (x) dif Gamma
  = integral_Gamma U_(i j) (x_d, x) t_j (x) dif Gamma
$

(mais termos de domínio se houver $f_i$ — DIBEM elástico no pacote, fora do núcleo da aula).

Discretização $arrow.r$ blocos $H,G$ com 2 DOFs por nó; `H_G_full_direct` e `solve` já conhecidos.

No pacote (`src/Elasticity/Fundamental.jl`), Kelvin 2D usa $nu$ efetivo
(plane strain default; plane stress via flag em `Elasticity(...)`):

```julia
# essência: U ~ (3-4ν) log(1/R) δ + r̂⊗r̂ ;  T ~ 1/R × ...
function fundamental(props::Elasticity, r::SVector{2}, n::SVector{2})
    ν = effective_nu(props)
    μ = props.mu
    # ... monta Mat 2×2 U e T
    return KernelPair(U, Tker)
end
```

Tensões em pontos internos: pós-processamento com núcleos $D$ e $S$ (fórmulas longas no material legado);
no curso, prefira sensores do `dad` \/ scripts de referência quando existirem.

== CDC no Gmsh

Quatro campos por curva: `"tx;ux;ty;uy"` — em cada direção, o flag 0\/1 escolhe valor prescrito de tração ou deslocamento (mesmo espírito de `"0;T"` \/ `"1;q"` do Laplace; detalhes no cap. *Laplace 2D*).

== Lab 1 — patch de Dirichlet

```julia
using DrWatson
@quickactivate :BEM
include(datadir("Laplace", "Laplace_dad.jl"))

msh = quadrado_elasticity(ndiv=15, show=false)
dad = format2d(msh, Elasticity(1.0, 0.3, 1.0); pontointerno=false)

ana = ana_elasticity_patch(; E=1.0, ν=0.3, εxx=0.01)
apply_analytical_bc!(dad, ana)   # u = ε · x em todo o contorno

H_G_full_direct(dad, 16)
solve(dad)
@show rel_error(dad)
plot_geo(dad)
```

Espere erro muito pequeno (patch polinomial). Se falhar: orientação, `plane_strain`\/stress, ou BC.

== Lab 2 — tubo pressurizado (referência no repo)

```julia
include(datadir("elastico", "iso", "pressurized_tube.jl"))
msh = mesh_pressurized_tube(; ndiv=16, show=false)
dad = format2d(msh, Elasticity(E, ν, 1.0); pontointerno=true)
H_G_full_direct(dad, 16)
solve(dad)
# ur_tube / σr_tube / σθ_tube no mesmo arquivo
```

== Exercícios

Resolva com o `BEM_gmsh`; ≥ 3 refinamentos; erros em deslocamento (e tensão se houver analítico).
Apêndice de erros: use $epsilon_2$ (`rel_error`) + $epsilon_infinity$ em sensores.

=== Cilindro pressurizado

$R_a = 50$, $R_b = 100$ (mm); $E = 200$ GPa; $nu = 0.32$; $P = 100$ N\/mm.

#image("../assets/elasticidade-2d/cilindro.png", width: 80%)

$
u_r = ((1 + nu) P R_a^2)/((R_b^2 - R_a^2) E) [(1 - 2 nu) r + R_b^2\/r] \
sigma_r = (P R_a^2)/(R_b^2 - R_a^2) (1 - R_b^2\/r^2) \
sigma_theta = (P R_a^2)/(R_b^2 - R_a^2) (1 + R_b^2\/r^2)
$

=== Placa com furo (Kirsch)

$R = 50$ mm; $E = 100$ GPa; $nu = 0.25$; tração remota $P = 1$ N\/mm.
Estado plano de *tensão*: ajuste $E',nu'$ se a SF do código estiver em deformação plana
(`Elasticity(..., plane_stress=true)` ou fórmulas do texto legado).

#image("../assets/elasticidade-2d/placa-furo.png", width: 80%)

$
nu' = nu\/(1+nu) ,
quad
E' = E (1 - (nu')^2\/(1+nu')^2) .
$

Campos $sigma_r, sigma_theta, tau_(r theta)$: solução de Kirsch (material do capítulo).

=== Viga com cisalhamento parabólico

$L = 48$, $D = 12$ (mm); $E = 300$ GPa; $nu = 0.3$; $P = 1000$ N.
Analítico de $u_1,u_2,sigma_x x, tau_x y$ no texto legado \/ figura.

#image("../assets/elasticidade-2d/viga.png", width: 80%)

Malhas: `data/elastico/iso/` ou Gmsh próprio.

== Leituras

- Cap. *Laplace 2D* (pipeline idêntico em espírito)
- `src/Elasticity/Fundamental.jl`, `data/elastico/iso/`
- #link("https://youtu.be/R6-_ECEQXRk")[gravação]
