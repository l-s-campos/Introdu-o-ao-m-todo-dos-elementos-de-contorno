// Elasticidade 2D — Kelvin + BEM_gmsh
// Pré-requisito: Laplace 2D

= Elasticidade 2D
<elasticidade-2d>

O pipeline é o *mesmo* do Laplace: malha de contorno $arrow.r$ SF $arrow.r$ $H,G$
$arrow.r$ CDC $arrow.r$ $A x = b$ $arrow.r$ `solve`. O que muda:

#table(
  columns: (auto, auto, auto),
  inset: 7pt,
  stroke: 0.5pt + luma(200),
  [*Peça*], [*Laplace*], [*Elasticidade 2D*],
  [Campo primário], [$T$ (escalar)], [$upright(bold(u))=(u_x,u_y)$],
  [Conjugado de Neumann], [$q = -k partial_n T$], [tração $upright(bold(t))$, $t_i = sigma_(i j) n_j$],
  [SF], [$T^*$, $q^*$], [tensores de Kelvin $U_(i j)$, $T_(i j)$],
  [DOFs por nó], [$1$], [$2$ ($H,G$ em blocos $2 times 2$)],
  [CDC no Gmsh], [`"0;T"` \/ `"1;q"`], [`"tipo_x;val_x;tipo_y;val_y"`],
  [Props], [`Laplace(k)`], [`Elasticity(E, ν, ρ; plane_strain\/stress)`],
)

Notação (glossário): Poisson do *material* $nu$; *não* use $v$ para o Poisson
(reserve $v$ à função peso da formulação, se aparecer).

== Objetivos

+ Ligar equilíbrio + Hooke à BIE de contorno (Kelvin).
+ Entender plane strain vs plane stress no `Elasticity(...)`.
+ Ler e montar CDC vetoriais no Gmsh.
+ Rodar o *patch test* de Dirichlet e um problema clássico do repo.
+ Reportar erro em $upright(bold(u))$ (e tensões nos exercícios) — apêndice de erros.

== Mapa

+ Problema contínuo (resumo operacional)
+ De Navier à equação integral
+ Kelvin no `BEM_gmsh` + diagonal de $H$
+ CDC `"tipo;val;tipo;val"`
+ Lab 1 — patch $upright(bold(u)) = bold(epsilon) upright(bold(x))$
+ Lab 2 — tubo pressurizado
+ Tensões (pós-processamento, sem muro de fórmulas)
+ Armadilhas
+ Exercícios (cilindro, Kirsch, viga)
+ Leituras

== Problema contínuo (o que importa para o BEM)

=== Equilíbrio, cinemática, Hooke

$
partial_j sigma_(i j) + f_i = 0 ,
quad
epsilon_(i j) = 1/2 (partial_i u_j + partial_j u_i) .
$

Isotrópico linear: $E$, $nu$, módulo de cisalhamento

$ mu = E \/ (2(1+nu)) . $

Eliminando $sigma$ e $epsilon$ chega-se às equações de *Navier* em $upright(bold(u))$.
A SF de Kelvin é a resposta a uma *força pontual* unitária nesse operador
(análogo a $T^*$ ser a resposta a uma fonte pontual em Laplace).

No contorno, *por nó e por direção* ($x$ e $y$): ou se prescreve deslocamento $u_i$,
ou se prescreve tração $t_i$ — nunca os dois no mesmo DOF (como $T$ vs $q$).

=== Estado plano

#table(
  columns: (auto, auto, auto),
  inset: 7pt,
  stroke: 0.5pt + luma(200),
  [*Modelo*], [*Hipótese*], [*No pacote*],
  [Plane *strain*], [$epsilon_z = 0$ (seção longa)], [`Elasticity(E,ν,ρ)` default],
  [Plane *stress*], [$sigma_z = 0$ (chapa fina)], [`Elasticity(E,ν,ρ; plane_stress=true)`],
)

O código usa um Poisson *efetivo* $tilde(nu)$ nos núcleos 2D:

- strain: $tilde(nu) = nu$;
- stress: $tilde(nu) = nu\/(1+nu)$  
  (e $lambda, mu$ em cache — `lame_constants` em `Structures.jl`).

#block(
  width: 100%,
  fill: luma(248),
  inset: 10pt,
  radius: 4pt,
  stroke: 0.5pt + luma(200),
)[
  *Não* misture: se o problema físico é chapa fina (Kirsch, viga fina do repo),
  use `plane_stress=true`. Se for cilindro “longo” \/ tubo do repo, use o default
  strain. Ajustar $E',nu'$ *na mão* e ainda pedir strain no construtor é receita
  de inconsistência — prefira a flag do `Elasticity`.
]

== Da PDE à equação integral

Mesma lógica do Laplace (resíduos $arrow.r$ identidades de Betti\/Green $arrow.r$ SF):

$
c_(i k)(x_d)\, u_k (x_d)
+ integral_Gamma T_(i j)(x_d,x)\, u_j (x)\, dif Gamma
=
integral_Gamma U_(i j)(x_d,x)\, t_j (x)\, dif Gamma
+
"termo de domínio se" f_i != 0 .
$

- $U_(i j)$: deslocamento fundamental (Kelvin) — singularidade fraca ($ln r$ em 2D).
- $T_(i j)$: tração fundamental — singularidade forte ($1\/r$).
- $c_(i k)$: termo livre (ângulo sólido tensorial); no contorno liso $c = (1\/2) I$.

Discretização com as *mesmas* $N_k$ do cap. *Indo para 2D* \/ Laplace:

$
upright(bold(u)) = sum_k N_k upright(bold(u))_k ,
quad
upright(bold(t)) = sum_k N_k upright(bold(t))_k
$

$arrow.r$ sistema global (blocos $2 times 2$ por par de nós)

$ H upright(bold(U)) = G upright(bold(T)) $

(vetores empilhados $u_x,u_y$ por nó). CDC reorganizam colunas $arrow.r$ $A x = b$;
`solve(dad)` faz `applyBC` + solve + espalha `dad.T`\/`dad.u` e trações.

Força de corpo $f_i$: DIBEM elástico existe no pacote — *fora* do núcleo desta aula
(análogo ao Poisson; ver trabalhos se precisar).

== Kelvin no `BEM_gmsh`

```julia
# Elasticity/Fundamental.jl — essência 2D
function fundamental(props::Elasticity, r::SVector{2}, n::SVector{2})
    ν = effective_nu(props)   # tilde(ν)
    μ = props.mu
    R = norm(r)
    # U ~ [(3-4ν) log(1/R) I + r̂⊗r̂] / (const · μ)
    # T ~ (1/R) · [termos com dr/dn, (1-2ν), ...]
    return KernelPair(U, T)   # Mat 2×2 cada
end
```

Montagem: o *mesmo* `H_G_full_direct(dad, npg)` do Laplace, no branch *vectorial*
(`Assembly_full.jl`): laços por fonte, `integrate_element`, e diagonal de $H$ por
*corpo rígido* (translação: linhas de cada bloco somam zero de forma análoga).

```julia
# após montar contribuições regulares (vectorial):
# H[ii, ii] corrigido para que translação rígida ⇒ tração nula
```

Você *não* precisa reprogramar Kelvin: confira que `dad.properties isa Elasticity`
antes de montar.

== CDC no Gmsh — formato real do pacote

O parser (`parse_pairs` em `Input.jl`) lê o nome do grupo físico como

`tipo_1; valor_1; tipo_2; valor_2`

para os dois DOFs $(x,y)$:

#table(
  columns: (auto, auto),
  inset: 8pt,
  stroke: 0.5pt + luma(200),
  [*`tipo`*], [*Significado*],
  [`0`], [Dirichlet: o `valor` é $u$ naquela direção],
  [`1`], [Neumann: o `valor` é $t$ (tração) naquela direção],
)

Exemplos:

#table(
  columns: (auto, auto),
  inset: 7pt,
  stroke: 0.5pt + luma(200),
  [*String*], [*Leitura*],
  [`"0;0;0;0"`], [$u_x=0$, $u_y=0$ (engaste)],
  [`"1;0;1;0"`], [$t_x=0$, $t_y=0$ (livre)],
  [`"1;0;0;0"`], [$t_x=0$, $u_y=0$ (rolo horizontal)],
  [`"0;0;1;0"`], [$u_x=0$, $t_y=0$ (rolo vertical)],
  [`"1;P;1;0"`], [$t_x=P$, $t_y=0$ (tração uniforme em $x$)],
)

#block(
  width: 100%,
  fill: luma(248),
  inset: 10pt,
  radius: 4pt,
  stroke: 0.5pt + luma(200),
)[
  Textos antigos às vezes escrevem `"tx;ux;ty;uy"`. No `BEM_gmsh` o que vale é
  *pares* `(tipo, valor)` por direção — igual ao espírito de `"0;T"` \/ `"1;q"`,
  só que *duas* vezes. Veja `quadrado_elasticity`, `pressurized_tube.jl`,
  `plate_with_hole.jl`.
]

Cantos: elementos *descontínuos* no campo (como no Laplace) evitam um nó com duas
CDCs incompatíveis. `format2d` + `tipo`\/ordem alinhados à malha.

== Lab 1 — Patch test (deformação uniforme)

Solução exata $upright(bold(u)) = bold(epsilon)\, upright(bold(x))$ (polinômio linear) —
o BEM com elementos adequados deve reproduzi-la com erro muito pequeno.

```julia
using DrWatson
@quickactivate :BEM
include(datadir("Laplace", "Laplace_dad.jl"))
# espelho: data/examples/elasticity_patch.jl

E, ν, εxx = 1.0, 0.3, 0.01
msh = quadrado_elasticity(ndiv=10, show=false, nome="ex_elast_patch", ordem=2)
dad = format2d(msh, Elasticity(E, ν, 1.0; plane_strain=true);
               tipo=2, pontointerno=false)

ana = ana_elasticity_patch(; E=E, ν=ν, εxx=εxx)
apply_analytical_bc!(dad, ana)   # sobrescreve BC/BV com u = ε·x (e t se misto)

H_G_full_direct(dad; npg=12)
solve(dad)
@show rel_error(dad)             # ‖u - u_ana‖₂ / ‖u_ana‖₂
plot_geo(dad)
```

`ana_elasticity_patch` (`Analytical.jl`): $u = epsilon upright(bold(x))$ e tração de
Hooke em plane strain a partir de $sigma$ constante.

*Esperado:* `rel_error` pequeno (no exemplo do repo, tipicamente $<< 10^(-1)$ já em malha grossa; refine se precisar).
*Se falhar:* malha\/orientação, `plane_strain` inconsistente com o analítico, ou `tipo`≠ordem.

Default do `quadrado_elasticity` *sem* `apply_analytical_bc!`: esquerda engastada
`"0;0;0;0"`, resto livre `"1;0;1;0"` — útil para outros testes; o patch *exige*
o `apply_analytical_bc!` (Dirichlet em todo o contorno).

== Lab 2 — Tubo pressurizado (repo)

```julia
include(datadir("elastico", "iso", "pressurized_tube.jl"))
# define Ra, Rb, E, ν, P, ur_tube, σr_tube, σθ_tube, mesh_pressurized_tube

msh = mesh_pressurized_tube(; ndiv=16, show=false)
dad = format2d(msh, Elasticity(E, ν, 1.0; plane_strain=true); pontointerno=true)

# CDCs do .geo: simetria nos eixos; confira se a pressão interna no arco
# está como você espera (ajuste BV/BC nos nós do furo se o gerador deixar livre).
# Tração de pressão no sólido: t = -P n  (n exterior ao material).

H_G_full_direct(dad, 16)
solve(dad)
plot_geo(dad)

# sensores: compare u_r numérico com ur_tube(r) nos nós / internos
```

Analítico (Lamé, pressão interna $P$, exterior livre) — como no arquivo do repo:

$
sigma_r (r) = c (1 - R_b^2\/r^2) ,
quad
sigma_theta (r) = c (1 + R_b^2\/r^2) ,
quad
c = R_a^2 P \/ (R_b^2 - R_a^2) ,
$

$
u_r (r) = ((1+nu)\/E)\, c\, [ (1-2 nu) r + R_b^2\/r ]
quad "(plane strain do arquivo)" .
$

== Tensões (pós-processamento)

Com $upright(bold(u)), upright(bold(t))$ no contorno, tensões em pontos *interiores*
vêm de integrais com núcleos $D_(k i j)$ e $S_(k i j)$ (hipersingulares no contorno).

#block(
  width: 100%,
  fill: rgb("#ecfdf5"),
  inset: 10pt,
  radius: 4pt,
  stroke: 0.5pt + rgb("#99f6e4"),
)[
  *Nesta aula.* Não decore as páginas de $S$ e $D$. Use:
  (1) deslocamentos e `rel_error` no patch;
  (2) nos exercícios, sensores onde o repo já dá $sigma$ analítico
  (`σr_tube`, Kirsch, viga) e compare componentes em pontos internos \/ contorno
  com a API de pós-processamento disponível na sua versão do `BEM_gmsh`.
  Fórmulas completas: material legado \/ livros de BEM elástico.
]

== Armadilhas

#table(
  columns: (auto, auto),
  inset: 7pt,
  stroke: 0.5pt + luma(200),
  [*Sintoma*], [*Causa típica*],
  [Patch com erro grande], [`plane_strain` ≠ analítico; ordem\/`tipo`; BC não aplicadas],
  [Corpo “voa” (modo rígido)], [Só Neumann em todo o contorno — prenda o corpo rígido],
  [Tração com sinal trocado], [Normal para dentro; pressão $t = -P n$ com $n$ exterior],
  [Kirsch longe do analítico], [Usou strain em problema de stress (ou o contrário)],
  [Cantos ruins], [Elemento contínuo com CDC mista — use descontínuo],
  [`parse` estranho da CDC], [String com número ímpar de campos; label sem pares],
)

== Exercícios

`BEM_gmsh` + malhas em `data/elastico/iso/` (ou Gmsh próprio).
≥ 3 refinamentos. Apêndice de erros: $epsilon_2$ (`rel_error` quando houver ana)
e $epsilon_infinity$ em sensores de $upright(bold(u))$ e $sigma$.

=== E0 — Patch

Rode o Lab 1 com `ndiv in {6,10,16}` e `tipo in {1,2}`. Tabela `N`, `rel_error`.
Comente a taxa (deve saturar em erro de integração\/máquina se o patch for exato).

=== E1 — Cilindro \/ tubo pressurizado

$R_a = 50$, $R_b = 100$ (mm); $E = 200$ GPa; $nu = 0.32$; $P = 100$ N\/mm
(valores alinhados a `pressurized_tube.jl` em MPa·mm).

#image("../assets/elasticidade-2d/cilindro.png", width: 80%)

Compare $u_r$ e $sigma_r, sigma_theta$ com Lamé; simetria de quarto de círculo.
Documente as CDCs de simetria e a pressão no raio interno.

=== E2 — Placa com furo (Kirsch)

$R = 50$ mm; tração remota $P$ em $x$; $E$, $nu$ como em `plate_with_hole.jl`.
*Plane stress:* `Elasticity(..., plane_stress=true)`.

#image("../assets/elasticidade-2d/placa-furo.png", width: 80%)

$
sigma_(r r) = (P\/2)(1 - a^2\/r^2)
  + (P\/2)(1 - 4 a^2\/r^2 + 3 a^4\/r^4) cos 2 theta ,
$
$
sigma_(theta theta) = (P\/2)(1 + a^2\/r^2)
  - (P\/2)(1 + 3 a^4\/r^4) cos 2 theta ,
$
$
sigma_(r theta) = -(P\/2)(1 + 2 a^2\/r^2 - 3 a^4\/r^4) sin 2 theta .
$

(com $a=R$; confira a forma exata no arquivo do repo se os coeficientes diferirem.)
Malha: `mesh_plate_with_hole`. Erros em sensores perto do furo e longe.

=== E3 — Viga em balanço (cisalhamento parabólico)

$L=48$, $D=12$; $E$, $nu$, $P$ como `cantilever_beam.jl` (plane stress).
Analítico:

$
u_1 = - (P y)\/(6 E I) [ (6L-3x)x + (2+nu)(y^2 - D^2\/4) ] ,
$
$
u_2 = (P)\/(6 E I) [ 3 nu y^2 (L-x) + (4+5 nu) D^2 x \/ 4 + (3L-x) x^2 ] ,
$
$
sigma_(x x) = - P(L-x) y \/ I ,
quad
tau_(x y) = - P\/(2 I) (D^2\/4 - y^2) ,
quad
I = D^3\/12 .
$

#image("../assets/elasticidade-2d/viga.png", width: 80%)

O gerador deixa o extremo livre com `"1;0;1;0"`: para o exercício, imponha a tração
parabólica $t_y (y)$ no extremo $x=L$ (via BV nos nós ou grupo físico adequado) e
engaste em $x=0$. Compare $u$ e $sigma$ em $x=L\/2$.

== O que fica para depois

- Trinca dual \/ $K_I$ — trabalhos, proposta A (`src/Crack/`).
- Anisotropia Lekhnitskii — `data/elastico/aniso/`.
- Contato Hertz — proposta E.
- DIBEM com força de corpo — espelho do Poisson.

== Leituras e código

- Cap. *Laplace 2D* — pipeline $H,G$, CDC, diagonal, `solve`
- Cap. *Indo para 2D* — $N_k$, $J$, $upright(bold(n))$, descontínuo
- *Apêndice: medidas de erro*
- `src/Elasticity/Fundamental.jl` — Kelvin
- `src/Core/Structures.jl` — `Elasticity`, `plane_stress`
- `src/Core/Analytical.jl` — `ana_elasticity_patch`
- `src/Core/Input.jl` — `parse_pairs`, `format2d`
- `data/examples/elasticity_patch.jl`
- `data/elastico/iso/pressurized_tube.jl`, `plate_with_hole.jl`, `cantilever_beam.jl`
- #link("https://youtu.be/R6-_ECEQXRk")[gravação]
