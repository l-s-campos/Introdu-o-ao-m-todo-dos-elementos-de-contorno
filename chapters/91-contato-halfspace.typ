// Contato em semi-espaço / semi-plano (extra)
// BEM_gmsh: ContactHalfPlane2D, ContactHalfSpace (Pohrt–Li)
// Fora da trilha obrigatória de 30 h — afim da proposta E dos trabalhos finais

= Contato por BEM de half-space (extra)
<contato-halfspace>

#block(
  width: 100%,
  fill: luma(245),
  inset: 10pt,
  radius: 4pt,
  stroke: 0.5pt + luma(180),
)[
  *Material extra / opcional.* Não faz parte da trilha obrigatória de 30 h.
  Usa o módulo de *contato em semi-espaço* do
  #link("https://github.com/l-s-campos/BEM_gmsh")[`BEM_gmsh`]
  (`src/Contact/`), *não* o pipeline `format2d` + `H_G_full_direct` das aulas de
  Laplace\/elasticidade em domínio limitado.

  Pré-requisitos úteis: *Elasticidade 2D* (Hooke, tração) e a ideia de solução
  fundamental. Aprofunda a *proposta E* dos trabalhos finais.
]

== Objetivos

+ Entender o modelo de *half-plane* (2D) e *half-space* (3D) para contato.
+ Ligar núcleos de influência (Flamant \/ Boussinesq–Cerruti) à convolução na grade.
+ Rodar Hertz *linha* (`ContactHalfPlane2D`) e Hertz *esfera* (`ContactHalfSpace`).
+ Ver o esboço de deslizamento parcial (Coulomb) no half-space.
+ Comparar com analítico e reportar erro (apêndice).

== Mapa

+ Contato ≠ BEM de contorno fechado
+ Modelo half-plane / half-space
+ Núcleos e convolução (FFT)
+ Problema normal (zona ativa)
+ Lab 1 — Hertz linha (2D)
+ Lab 2 — Hertz esfera (3D, grade)
+ Lab 3 — partial slip (opcional)
+ Armadilhas
+ Exercícios
+ Leituras

== Contato no curso vs. este módulo

#table(
  columns: (auto, auto, auto),
  inset: 7pt,
  stroke: 0.5pt + luma(200),
  [*Aspecto*], [*Laplace \/ elasticidade (trilha)*], [*Contato half-space*],
  [Geometria], [Contorno fechado $Gamma = partial Omega$], [Superfície plana infinita (grade)],
  [Incógnita típica], [$T$ ou $upright(bold(u))$ no contorno], [Pressão $p$ (e tração tangencial $tau$) na grade],
  [Operador], [$H,G$ densos ou H-matriz], [Convolução com núcleo de influência $K$],
  [API], [`format2d`\/`3d`, `solve`], [`solve_line_contact`, `solve_normal_contact`, …],
  [Referência clássica], [Green \/ Kelvin], [Hertz; Pohrt & Li (2014)],
)

#block(
  width: 100%,
  fill: rgb("#ecfdf5"),
  inset: 10pt,
  radius: 4pt,
  stroke: 0.5pt + rgb("#99f6e4"),
)[
  *Mesma filosofia BEM.* A resposta elástica a cargas na superfície é dada por uma
  *SF* (Flamant em 2D, Boussinesq\/Cerruti em 3D). O “sistema” vira
  $u = K * p$ (convolução) na grade, com restrições de contato ($p >= 0$, gap $>= 0$).
]

== Modelo físico (resumo)

=== Half-plane (linha de contato, 2D)

- Sólido: semi-plano elástico $y <= 0$ (plane strain).
- Indentador rígido (ex. cilindro de raio $R$) ou gap inicial $g_0(x)$.
- Pressão normal $p(x) >= 0$ só na zona de contato $C$.
- Deflexão normal $u(x)$ *para dentro* do sólido (convenção do pacote).

Módulo de contato (plane strain):

$ E^* = (2 G)\/(1 - nu) = E \/ (1 - nu^2) . $

No código: `ElasticHalfPlane2D(G, ν; h=…)` e `contact_modulus(hp)`.

=== Half-space (contato 3D na grade)

- Semi-espaço $z <= 0$; grade retangular na superfície.
- Gap $g_0(x,y)$ (ex. esfera: $g_0 = (x^2+y^2)\/(2R)$).
- Indentação rígida $delta > 0$ “para dentro”.
- `ElasticHalfSpace(G, ν; hx, hy)`.

=== Hertz (referência)

*Linha* (cilindro \/ half-plane), força $F$ por unidade de comprimento:

$ a = sqrt((4 F R)\/(pi E^*)) ,
  quad
  p_0 = (2 F)\/(pi a) ,
  quad
  p(x) = p_0 sqrt(1 - (x\/a)^2)_+ . $

*Esfera* no half-space, indentação $delta$:

$ a = sqrt(R delta) ,
  quad
  F = (4\/3) E^* sqrt(R)\, delta^(3\/2) ,
  quad
  p_0 = (3\/2) F \/ (pi a^2) . $

(O pacote expõe `hertz_line`, `hertz_line_pressure`, `hertz_halfwidth`, `hertz_pressure`.)

== Núcleos de influência e convolução

=== 2D — Flamant integrado no segmento

Deflexão no centro da célula $i$ por pressão unitária na célula $i+d_i$:

```julia
# ContactHalfPlane2D.jl — essência
# K ~ -((1-ν)/(π G)) ∫ ln|x-ξ| dξ  no segmento de comprimento h
influence_coeff_2d(di, hp)
```

A matriz densa seria Toeplitz; na prática usa-se FFT (`precompute_kernel_2d`,
`fc_forward_2d`, `fc_inverse_2d`).

=== 3D — Love \/ Cerruti em retângulos (Pohrt–Li)

Coeficientes $K_(a b)$ para pares de componentes (normal–normal $K_(z z)$, etc.),
integrados na célula $h_x times h_y$:

```julia
influence_coeff(Kzz, di, dj, hs)
precompute_kernels(nx, ny, hs; components=(Kzz,))
u = fc_forward(p, Kzz, prep)       # u = K * p
p = fc_inverse(u_target, mask, Kzz, prep)
```

Sinal: pressão e $u_z$ positivos *para dentro* do half-space.

== Problema normal (zona ativa)

Condições de Karush–Kuhn–Tucker (fricção nula):

$
p >= 0 ,
quad
g = g_0 - delta + u >= 0 ,
quad
p\, g = 0 .
$

No contato $C$: $u = delta - g_0$ (fechamento do gap); fora: $p = 0$.

Algoritmo típico (Polonsky–Keer \/ active set no pacote):

1. chute de $C$ (interferência geométrica $delta > g_0$);
2. inverter $u = K * p$ só em $C$ com $p$ livre;
3. cortar $p < 0$ (sair de $C$);
4. onde ainda penetra, entrar em $C$;
5. repetir até estabilizar.

```julia
sol = solve_line_contact(gap0, δ, hp)           # 2D
sol = solve_normal_contact(gap0, δ, hs)         # 3D grade
# sol.p, sol.contact, sol.force, sol.u
```

== Lab 1 — Hertz linha (half-plane)

Espelho de `scripts/hertz_line_2d.jl`:

```julia
using DrWatson
@quickactivate :BEM

G, ν = 1.0, 0.3
N = 256
L = 2.0
x = collect(range(-L, L; length=N))
h = x[2] - x[1]
hp = ElasticHalfPlane2D(G, ν; h=h)

R = 1.0
F_target = 0.05
hz = hertz_line(F_target, R, hp)   # a, p0, Estar analíticos
@show hz.a, hz.p0, hz.Estar

gap0 = @. x^2 / (2R)

# bisseção em δ até a zona de contato ~ a_Hertz
# (no half-plane o log torna δ absoluto ambíguo — casa-se a largura)
δ_lo, δ_hi = 1e-6, 0.5
sol = nothing
for _ in 1:40
    δ_mid = 0.5 * (δ_lo + δ_hi)
    sol = solve_line_contact(gap0, δ_mid, hp; tol=1e-12)
    a_num = begin
        idx = findall(sol.contact)
        isempty(idx) ? 0.0 : 0.5 * (x[maximum(idx)] - x[minimum(idx)])
    end
    if a_num < hz.a
        δ_lo = δ_mid
    else
        δ_hi = δ_mid
    end
end

p0_num = maximum(sol.p)
@show sol.force, p0_num, p0_num / hz.p0

p_hz, _ = hertz_line_pressure(sol.force, R, hp, x)
# plot: sol.p vs p_hz
```

*Esperado:* $p(x)\/p_0$ próximo da semi-elipse de Hertz; $a_"num"\/a$ e $p_0$ com erro
caindo ao refinar $N$ (grade maior e domínio $L$ adequado).

== Lab 2 — Hertz esfera (half-space)

Espelho de `scripts/contact_pohrt_li.jl` (parte normal):

```julia
using DrWatson
@quickactivate :BEM

E, ν = 1.0, 0.3
G = E / (2(1 + ν))
N = 64
L = 2.0
x = range(-L, L; length=N)
y = range(-L, L; length=N)
hx, hy = x[2] - x[1], y[2] - y[1]
hs = ElasticHalfSpace(G, ν; hx=hx, hy=hy)

R = 1.0
gap0 = [(x[i]^2 + y[j]^2) / (2R) for i in 1:N, j in 1:N]
δ = 0.05

sol = solve_normal_contact(gap0, δ, hs; tol=1e-7)
Estar = contact_modulus(hs)
a_hz = sqrt(R * δ)
F_hz = (4/3) * Estar * sqrt(R) * δ^(3/2)
p0_hz = 1.5 * F_hz / (π * a_hz^2)

@show count(sol.contact), sol.force, sol.force / F_hz, maximum(sol.p) / p0_hz
# heatmap(sol.p)
```

*Esperado:* $F_"num"\/F_"Hertz"$ e $max p \/ p_0$ próximos de 1 em grade fina
($N >= 64$ já dá ideia; $N=128$ melhora a borda do contato).

== Lab 3 — Deslizamento parcial (opcional)

Com $p$ e a máscara de contato do Lab 2:

```julia
μf = 0.3
d = 0.2 * μf * δ
ps = solve_partial_slip(sol.p, sol.contact, d, μf, hs; direction=:x, tol=1e-6)
@show count(ps.stick), count(ps.slip), ps.force_t
@show abs(ps.force_t) / (μf * sol.force)   # ≤ 1
```

Ideia Cattaneo–Mindlin: núcleo em *stick*, anel em *slip* com $|tau| = mu p$.
O solver do half-space é a versão em grade (Pohrt–Li §5); analíticos 2D de
cilindro estão em `CattaneoMindlin.jl` (`cattaneo_shear`, …).

== Armadilhas

#table(
  columns: (auto, auto),
  inset: 7pt,
  stroke: 0.5pt + luma(200),
  [*Sintoma*], [*Causa típica*],
  [Contato “cola” na borda do domínio], [$L$ pequeno demais — alargue a grade],
  [$p$ oscila \/ ruído], [Grade grossa; `tol` frouxo; plot sem suavizar],
  [Força longe de Hertz], [$delta$ ou $R$ inconsistentes; $E^*$ errado ($G,nu$)],
  [Half-plane: $delta$ não casa com teoria], [Kernel log — normalize por $a$ ou $F$, não por $delta$ absoluto],
  [Sinal de $u$\/$p$ invertido], [Convenção “para dentro do sólido”],
  [Partial slip estranho], [$nu != 1\/2$: acoplamento N–T negligenciado no algoritmo simples],
)

== Exercícios

Apêndice de erros. Sempre declare $N$, $L$, $E^*$ e o que é analítico.

=== E1 — Convergência Hertz linha

Lab 1 com $N in {64,128,256,512}$ (e $L$ fixo generoso).
Tabela: $a_"num"\/a$, $p_0$"num"\/$p_0$, $F$"num"\/$F$, RMSE de $p$ na zona de contato
(amostrado vs `hertz_line_pressure`). Uma figura $p(x)$.

=== E2 — Convergência Hertz esfera

Lab 2 com $N in {32,64,128}$. Tabela $F\/F_"hz"$, $max p\/p_0$, fração de nós em contato.
Comente o custo (FFT) ao dobrar $N$.

=== E3 — Partial slip ou Cattaneo 2D

*Ou* complete o Lab 3 e discuta stick\/slip e $|F_t|\/(mu F_n)$;
*ou* use `solve_cattaneo_halfplane` \/ `cattaneo_shear` e compare $q(x)$ ao analítico
para um $Q\/(f P)$ fixo.

=== E4 — (Desafio) Rugosidade ou indentador não circular

Altere $g_0$ (ex. superposição de senos de pequena amplitude, ou gap elíptico) e
compare força–indentação com o caso liso. Sem analítico completo: use refino e
conservação $F = sum p\, h_x h_y$.

== Ligação com a trilha e trabalhos

- *Elasticidade 2D:* mesmo $E$, $nu$, noção de tração; aqui o “contorno” é a grade plana.
- *Trabalhos, proposta E:* este capítulo é o mapa mínimo; a monografia exige
  Pareto de refino, mais de um algoritmo e discussão de limites do half-space.
- Contato entre *dois corpos malhados* (NTN\/NTS, `MortarContact2D`) é outro nível —
  não misturar com half-space na primeira leitura.

== Leituras e código

- Pohrt & Li (2014), *Physical Mesomechanics* — formulação completa N+T na grade
- Johnson, *Contact Mechanics* — Hertz, Cattaneo–Mindlin
- `src/Contact/ContactHalfPlane2D.jl` — Flamant, Hertz linha
- `src/Contact/ContactHalfSpace.jl` — Pohrt–Li, `solve_normal_contact`, `solve_partial_slip`
- `src/Contact/CattaneoMindlin.jl` — analíticos e solvers 2D de partial slip
- `scripts/hertz_line_2d.jl`, `scripts/contact_pohrt_li.jl`
- `scripts/compare_contact_acceleration.jl` — variantes de aceleração
- Cap. *Elasticidade 2D* · *Trabalhos finais* (proposta E) · *Apêndice: medidas de erro*
