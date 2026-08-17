// Poisson 2D — RASCUNHO (não incluído em main.typ / site até aprovação)
// DIBEM + fonte de domínio; base: Laplace 2D

= Poisson 2D
<poisson-2d>

Laplace montou $H T = G q$ *sem* fonte. Aqui

$ nabla^2 T = f quad "em" Omega $

e a equação integral ganha um termo de *domínio*:

$ c(x_d) T(x_d)
  = integral_Gamma T q^* dif s
  - integral_Gamma T^* q dif s
  + integral_Omega T^* f dif Omega . $

O BEM clássico não quer malha volumétrica de EF. O `BEM_gmsh` trata o termo extra por *DIBEM*
(Direct Interpolation BEM): RBF + redução a integrais de contorno.

== Objetivos

+ Entender por que aparece $integral_Omega T^* f$.
+ Seguir a ideia DIBEM (RBF $arrow.r$ matriz $M$) sem se afogar em matrizes densas na página.
+ Rodar `H_G_full_direct` + `DIBEM` + `solve` (e, se quiser, um transiente).
+ Resolver um problema estacionário com $f$ conhecida e medir erro (apêndice de erros).

== Mapa

+ PDE Poisson e o termo de domínio
+ Ideia DIBEM (o que é $M$)
+ Snippets do pacote
+ Lab mínimo
+ Transiente (opcional curto)
+ Exercícios

== Da PDE ao termo de domínio

Com a mesma SF de Laplace, a identidade de Green + $-nabla^2 T^* = delta$ deixa

$ integral_Omega T^* (nabla^2 T) dif Omega
  arrow.r
  integral_Omega T^* f dif Omega $

no lado da equação integral. Tudo que já sabemos ($H,G$, CDC, diagonal de $H$) permanece;
só entra um vetor (ou operador) de domínio no sistema final.

== Ideia DIBEM (sem parede de matrizes)

1. Amostramos $f$ (e a interpolação) em pontos do contorno *e* internos
   (`pontointerno=true` — os internos da superfície Gmsh \/ lista do `dad`).
2. Escrevemos uma interpolação por RBF $phi.alt(r)$ (no pacote, default `PHS`, etc.):

$ f(x) approx sum_j phi.alt(|x - x_j|) alpha_j
  quad arrow.r.double quad
  bold(alpha) = F^(-1) bold(f)_"nos" . $

3. A integral $integral_Omega T^* f$ vira combinação de integrais *radiais* de $phi.alt$ e de $T^*$,
   reduzidas ao contorno (mesmo espírito da integração radial do cap. *Indo para 2D*).
4. O resultado global é um operador $M$ tal que a contribuição de domínio é da forma $M bold(f)$
   (ou $M accent(T, dot)$ no transiente).

Diagonal de $M$: como em $H$, o termo singular se fecha por identidade
(caso $f equiv 1$ $arrow.r$ integrais $I_1$ conhecidas \/ montadas):

$ M_(i i) = I_(1 i) - sum_(j != i) M_(i j) . $

#block(
  width: 100%, fill: luma(248), inset: 10pt, radius: 4pt, stroke: 0.5pt + luma(200),
)[
  *Leitura.* Você *não* precisa rederivar todas as matrizes $F,D,s_f$ na mão para usar o curso.
  Precisa saber: *existe* $M$; ela vive no cache após `DIBEM(dad)`; entra no RHS ou no marchador temporal.
]

== No `BEM_gmsh`

```julia
# Laplace/Domain.jl — essência DIBEM_dense
# F_ij = φ(|x_i-x_j|),  D_ij ~ T*(x_i,x_j)
# IF, ID = integrais radiais no contorno
# M ≈ (IF' / F) .* D   e depois
# M[i,i] = -sum(M[i,:]) + ID[i]
```

API de alto nível:

```julia
using DrWatson
@quickactivate :BEM
include(datadir("Laplace", "Laplace_dad.jl"))

props = Laplace(1.0)
msh   = quadrado(ndiv=24, show=false)
dad   = format2d(msh, props; pontointerno=true)  # internos → centros RBF

H_G_full_direct(dad, 16)
DIBEM(dad)                    # cache M; default RBF do pacote
# DIBEM(dad; rbf=PHS(3; poly_deg=0))

# Estacionário *sem* fonte ainda é solve(dad) após H,G.
# Com fonte de domínio estacionária, a contribuição M*f entra no RHS
# (ver Domain.jl / docs do repo e exemplos em data/).
solve(dad)
plot_geo(dad)
```

Pontos internos da malha *não* são elementos finitos de volume: são centros de interpolação
para o DIBEM e sensores de pós-processamento.

== Transiente (ponte curta)

Com $M$ no papel de “massa” de domínio:

- calor: `solve_transient(dad, Δt, t_f)` (1ª ordem);
- onda \/ 2ª ordem: `solve_Houbolt` ou `solve_transient_o2`.

Detalhes e CFL: proposta C dos *Trabalhos finais*. Aqui basta saber a ordem:

```julia
H_G_full_direct(dad, 16)
DIBEM(dad)
# sol = solve_transient(dad, 0.01, 1.0)
```

== Lab mínimo sugerido

1. Repita o quadrado $T = x$ (Laplace) com `pontointerno=true` e confira que, *sem* fonte, o DIBEM não estraga o patch se $f=0$.
2. Monte um caso com $f$ constante ou polinomial baixa (membrana ou elipse dos exercícios) e compare $T$ em internos com o analítico.
3. Tabela `ndiv` × `rel_error` e, se pedir o enunciado, $epsilon_1$ \/ $epsilon_infinity$ no apêndice.

== Exercícios

Use o apêndice de erros. Varie a discretização.

+ *Membrana triangular.* Contorno fixo, $S nabla^2 w = -f$, $a = 5$, $f = 10$, $S = 1$ (unidades do enunciado).

#image("../assets/poisson-2d/membrana.png", width: 80%)

$ w = - f/(2 S) [1/2 (x^2 + y^2) - 1/(a sqrt(3)) (y^3 - 3 x^2 y) - 1/18 a^2] $

Erro médio e $L_2$ em ≥ 100 internos; mapa de $w$.

+ *Elipse estacionária.* $nabla^2 u = 4 - x^2$ com analítico do material legado; erros de potencial (internos) e fluxo (contorno).

#image("../assets/poisson-2d/elipse.png", width: 80%)

+ *Transiente 1D efetivo (placa\/cubo unitário).* Face aquecida a $T=1$, resto conforme figura; propriedades unitárias; série de Fourier do material. Use `solve_transient` (ou Houbolt) e compare sensores no tempo.

#image("../assets/poisson-2d/cubo-transiente.png", width: 80%)

+ *Extra — cilindro oco transiente.* $a=1$, $b=2$, $T(a,t)=1$; série com Bessel; raízes via `Roots.jl`. Simetria + DIBEM + marchador.

== Leituras e código

- Cap. *Laplace 2D* (pipeline $H,G$, CDC)
- `src/Laplace/Domain.jl`, `Domain_fast.jl` — `DIBEM`
- `src/Laplace/Solver.jl` — transiente
- Apêndice: medidas de erro
- #link("https://youtu.be/uSREar_ejnM")[gravação]
