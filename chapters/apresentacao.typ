// Apresentação (versão didática para quem não conhece o BEM)
// Aula 1 — sem BEM_gmsh; gráficos com Plots.jl

= Apresentação
<apresentacao>

== Objetivos desta aula

Ao final, você deve ser capaz de:

+ Explicar, em uma frase, o que o método dos elementos de contorno (BEM / MEC) faz de diferente em relação a métodos de domínio (diferenças finitas, elementos finitos).
+ Citar duas situações em que o BEM costuma ser atrativo e duas em que ele é desconfortável.
+ Partir de $T'' = 0$ em um intervalo, chegar ao sistema $H T = G Q$ nas pontas e aplicar condições de contorno.
+ Calcular $T$ em pontos *internos* a partir só dos dados de contorno já resolvidos.
+ Resolver no Julia (com `Plots`) os dois casos-modelo e o exercício de convecção.

Não é necessário, nesta aula, clonar repositório de produção nem gerar malha com Gmsh: o laboratório é o BEM *em uma dimensão*, onde o “contorno” são só dois pontos.

== Gancho: um problema que só “mora” nas pontas

Considere uma barra (ou a espessura de uma grande placa) no intervalo

$ Omega = (x_0, x_f), quad Gamma = {x_0, x_f}. $

Em regime permanente, sem geração de calor, a temperatura $T(x)$ obedece à equação de Laplace 1D

$ (d^2 T)/(d x^2) = 0, $

com *duas* condições de contorno escolhidas entre temperatura e fluxo

$ Q(x) := (d T)/(d x) $

nas extremidades. Exemplos:

- Dirichlet: $T(x_0)$ e $T(x_f)$ dados;
- Neumann / mista: uma temperatura e um fluxo, ou dois fluxos compatíveis com o balanço.


Em 1D o contorno $Gamma$ *é* o par de pontas. A pergunta do BEM fica literal:

#align(center)[
  _Se a física no interior está toda amarrada pelo operador diferencial,
  será que basta conhecer o que acontece em $Gamma$?_
]

Para Laplace 1D a resposta analítica é óbvia ($T$ é reta). O valor da aula é outro: montar a *mesma* lógica que, em 2D/3D, transforma um problema de domínio em um problema só de contorno.

== O BEM em uma figura

#image("../assets/apresentacao/bem-overview.jpeg", width: 85%)

Ideia-chave (guarde esta frase):

#align(center)[
  *Identidade integral + solução fundamental
  $=>$ equações cujas incógnitas vivem no contorno $Gamma$.*
]

Comparação rápida com métodos de domínio:

#table(
  columns: (auto, auto, auto),
  inset: 7pt,
  stroke: 0.5pt + luma(200),
  [*Aspecto*], [*FEM / diferenças (domínio)*], [*BEM (contorno)*],
  [O que se malha], [todo o $Omega$], [só $Gamma = partial Omega$],
  [Incógnitas primárias], [nos nós do volume/área], [nos nós do contorno],
  [Interior], [obtido junto com a solução], [pós-processamento (avaliar a identidade num ponto interno)],
  [Matriz típica], [esparsa (interação local)], [cheia e, em geral, não simétrica],
  [Ingrediente extra], [funções de forma no domínio], [*solução fundamental* do operador],
)

#image("../assets/image.png", width: 85%)

Em 2D, malhar só a curva que cerca a peça reduz uma dimensão da discretização. Em troca, cada ponto de contorno “enxerga” todos os outros: a matriz deixa de ser esparsa. Não existe almoço grátis — só um trade-off diferente.

== Quando o BEM costuma valer a pena

*Vantagens típicas*

+ Só o contorno é discretizado: menos geometria para preparar quando $Omega$ é grande e a física é linear e homogênea.
+ Domínios infinitos ou semi-infinitos (solo, escoamento externo, espaço aberto) entram com naturalidade via solução fundamental adequada — sem “truncar caixas” enormes.
+ Valores internos são calculados *depois*, só onde interessam (útil em laços de otimização).
+ Boa resolução de concentrações de tensão / gradientes em bordas, fendas e cargas concentradas, quando a formulação está bem posta.

*Limitações típicas*

+ Matrizes cheias e não simétricas: custo e memória crescem rápido se a implementação for ingênua.
+ É preciso conhecer (ou saber derivar) a *solução fundamental* do operador. Nem todo material / não-linearidade tem SF simples.
+ Estruturas muito delgadas e alguns acoplamentos são chatas de tratar só com BEM clássico.
+ Não-linearidade ou heterogeneidade no *domínio* em geral reintroduz integrais de volume (ou truques equivalentes).

*Por que ainda se ensina pouco na graduação*

+ Formulação matemática mais densa na entrada (integrais, SF, singularidades).
+ Menos códigos didáticos curtos do que no universo do FEM.
+ Singularidades nas integrais exigem cuidado numérico.
+ Mudar a ideia padrão dos métodos de domínio tem um custo envolvido.

Este curso  ataca sobretudo formulação passo a passo, integração de termos delicados e a passagem 1D $->$ 2D.

== Notação mínima (aula 1)

Alinhada ao glossário do curso, com um atalho só para o 1D:

#table(
  columns: (auto, auto),
  inset: 7pt,
  stroke: 0.5pt + luma(200),
  [*Símbolo*], [*Significado aqui*],
  [$Omega = (x_0, x_f)$], [domínio 1D],
  [$Gamma = {x_0, x_f}$], [contorno],
  [$n$], [normal *exterior* a $Omega$: $n(x_0) = -1$, $n(x_f) = +1$],
  [$T$], [temperatura / potencial],
  [$Q = d T \/ d x$], [derivada espacial (atalho 1D)],
  [$x_d$], [ponto fonte (onde “colocamos” a SF)],
  [$T^* (x, x_d)$, $Q^* = partial T^* \/ partial x$], [solução fundamental e sua derivada],
  [$H$, $G$], [matrizes de influência: $H T = G Q$],
)

*Armadilha de notação.* No restante do curso o fluxo de contorno do código é

$ q := - k (partial T)/(partial n). $

Em 1D, $partial T \/ partial n = n \, Q$. Com $k = 1$, $q = - n Q$. Nesta aula trabalhamos com $Q$ para a álgebra ficar transparente; quando formos ao 2D, a incógnita de contorno volta a ser $q$.

== Da integração por partes ao contorno

Sejam $u$ e $v$ funções regulares em $[x_1, x_2]$. Integração por partes:

$
  integral_(x_1)^(x_2) u(x) v'(x) dif x
  = [u(x) v(x)]_(x_1)^(x_2)
  - integral_(x_1)^(x_2) v(x) u'(x) dif x.
$

O termo avaliado nas pontas *é* o contorno 1D. Com a normal exterior $n(x_1)=-1$, $n(x_2)=+1$,

$
  [u v]_(x_1)^(x_2)
  = u(x_2) v(x_2) n(x_2) + u(x_1) v(x_1) n(x_1)
  = sum_(x in Gamma) u(x) v(x) n(x),
$

que em dimensão maior se escreve $integral_Gamma u v n dif Gamma$ (ou $integral_Gamma u v upright(bold(n)) · dif upright(bold(Gamma))$).

#image("../assets/apresentacao/contorno-1d.png", width: 60%)


Leitura operacional:

+ Integração por partes *transfere* derivadas de uma função para a outra.
+ Cada transferência produz termos de contorno.
+ No BEM repetimos o processo até o operador diferencial cair inteiro sobre a função peso (a SF). O número de integrações acompanha a ordem do operador (Laplace: 2; biarmônico/viga: 4; …).

== Laplace 1D, em camadas

=== Camada 0 — problema forte

$ (d^2 T)/(d x^2) = 0 quad "em" quad (x_0, x_f), $

mais duas CDC entre ${T(x_0), T(x_f), Q(x_0), Q(x_f)}$.

=== Camada 1 — resíduo ponderado (ainda sem escolher o peso)

Para uma função peso $T^*$ suficientemente regular,

$ integral_(x_0)^(x_f) T^* (d^2 T)/(d x^2) dif x = 0. $

Duas integrações por partes levam a

$
  integral_(x_0)^(x_f) T^* T'' dif x
  = [T^* Q - T Q^*]_(x_0)^(x_f)
  - integral_(x_0)^(x_f) T (T^*)'' dif x.
$


=== Camada 2 — escolha da função peso: solução fundamental

Em vez de um peso polinomial arbitrário, o BEM escolhe $T^*$ especial:

$ - (d^2 T^* (x, x_d))/(d x^2) = delta(x - x_d). $

Aqui $delta$ é a delta de Dirac: a integral de $f(x) delta(x-x_d)$ devolve $f(x_d)$ se o intervalo cobre $x_d$, e $0$ caso contrário. Intuição: $T^*(·, x_d)$ é a resposta em todo o eixo de uma fonte unitária em $x_d$, para o operador $-d^2\/dif x^2$.

*Por que isso ajuda?* O termo de domínio vira amostragem do campo:

$
  integral_(x_0)^(x_f) T (T^*)'' dif x
  = - integral_(x_0)^(x_f) T(x) delta(x - x_d) dif x
  = - T(x_d) quad ("se" x_d in (x_0, x_f)).
$

Com $T'' = 0$, a identidade colapsa para uma relação *só com valores de contorno* e com o valor $T(x_d)$.

=== Camada 3 — SF explícita em 1D

Uma SF conveniente é

$
  T^*(x, x_d) = - 1/2 |x - x_d|,
  quad
  Q^*(x, x_d) = (partial T^*)/(partial x) = - 1/2 "sign"(x - x_d).
$

Confira: longe de $x_d$, $T^*$ é linear por partes, logo $(T^*)'' = 0$; o “salto” da derivada em $x_d$ produz a delta. Em 1D, $T^*(x_d, x_d) = 0$ — a singularidade é branda. Em 2D a SF de Laplace envolve $log r$ e a integração exige mais cuidado (aulas seguintes).

Substituindo e reorganizando, para $x_d$ no intervalo obtém-se a equação integral de contorno

$
  - T(x_d)
  = T(x_f) Q^*(x_f, x_d) - T^*(x_f, x_d) Q(x_f)
  - T(x_0) Q^*(x_0, x_d) + T^*(x_0, x_d) Q(x_0).
$

=== Camada 4 — colocação no contorno ($H T = G Q$)

O contorno só tem dois pontos. Colocamos a fonte em cada um deles: $x_d = x_0$ e $x_d = x_f$.

Valores da SF:

#table(
  columns: (auto, auto, auto),
  inset: 7pt,
  stroke: 0.5pt + luma(200),
  [*Quantidade*], [*$x_d = x_0$*], [*$x_d = x_f$*],
  [$T^*(x_0, x_d)$], [$0$], [$-(x_f-x_0)/2$],
  [$T^*(x_f, x_d)$], [$-(x_f-x_0)/2$], [$0$],
  [$Q^*(x_0, x_d)$], [$+1\/2$], [$+1\/2$],
  [$Q^*(x_f, x_d)$], [$-1\/2$], [$-1\/2$],
)

Equação em $x_d = x_0$:

$
  - T(x_0)
  = - 1/2 T(x_f) + (x_f - x_0)/2 Q(x_f) - 1/2 T(x_0).
$

Equação em $x_d = x_f$:

$
  - T(x_f)
  = - 1/2 T(x_f) - 1/2 T(x_0) - (x_f - x_0)/2 Q(x_0).
$

Reorganizando em matriz (com $T_0 = T(x_0)$, etc.):

$
  mat(0.5, -0.5; -0.5, 0.5)
  mat(T_0; T_f)
  =
  (x_f - x_0)
  mat(0, -0.5; 0.5, 0)
  mat(Q_0; Q_f).
$

Ou seja,

$
  H T = G Q,
  quad
  H = mat(0.5, -0.5; -0.5, 0.5),
  quad
  G = (x_f - x_0) mat(0, -0.5; 0.5, 0).
$

Leitura que você vai reencontrar em 2D:

- cada *linha* de $H$ e $G$ corresponde a uma colocação (um $x_d$);
- cada *coluna* corresponde à influência de um grau de liberdade de contorno ($T$ ou $Q$ naquele nó);
- $H$ multiplica temperaturas; $G$ multiplica fluxos $Q$.

*Condições de contorno.* Das quatro quantidades $(T_0, T_f, Q_0, Q_f)$, *duas* são dados e duas são incógnitas. Move-se para a esquerda tudo o que é desconhecido e para a direita tudo o que é conhecido, até obter $A x = b$ com $A$ $2 times 2$.

=== Camada 5 — pontos internos (pós-processamento)

Com $T$ e $Q$ já conhecidos *nas pontas*, a mesma identidade com $x_d in (x_0, x_f)$ devolve $T(x_d)$ sem resolver outro sistema. Para a SF deste capítulo:

$
  T(x_d)
  = 1/2 (T_0 + T_f)
  + (x_d - x_0)/2 Q_0
  + (x_d - x_f)/2 Q_f.
$

Isso materializa a vantagem “interior sob demanda”: o sistema linear é só de contorno; o campo interno é avaliação.

== Exemplos numéricos (Julia + Plots)

Ambiente mínimo desta aula:

```julia
using Pkg
Pkg.add("Plots")   # uma vez por ambiente
using Plots, LinearAlgebra
```

Matrizes $H$ e $G$ no intervalo $[x_0, x_f]$:

```julia
function bem1d_matrices(x0, xf)
    L = xf - x0
    H = [0.5 -0.5; -0.5 0.5]
    G = L * [0.0 -0.5; 0.5 0.0]
    return H, G, L
end
```

Montagem a partir de $H T - G Q = 0$: colunas de incógnitas vão para $A$; termos conhecidos vão para $b$.

```julia
"""Resolve H*T = G*Q com CDC mistas.
`T_data` / `Q_data`: valor numérico se conhecido, `NaN` se incógnita.
Em cada nó, prescreva exatamente uma entre T e Q."""
function solve_bem1d(H, G, T_data, Q_data)
    n = 2
    T_known = .!isnan.(T_data)
    Q_known = .!isnan.(Q_data)
    @assert count(T_known) + count(Q_known) == n
    @assert all(T_known .!= Q_known)

    A = zeros(n, n)
    b = zeros(n)
    col_of_T = zeros(Int, n)
    col_of_Q = zeros(Int, n)
    col = 0

    for i in 1:n
        if T_known[i]
            b .-= H[:, i] * T_data[i]
        else
            col += 1
            col_of_T[i] = col
            A[:, col] .+= H[:, i]
        end
    end
    for i in 1:n
        if Q_known[i]
            b .+= G[:, i] * Q_data[i]
        else
            col += 1
            col_of_Q[i] = col
            A[:, col] .-= G[:, i]   # -G * Q_desconhecido
        end
    end

    x = A \ b
    T = zeros(n)
    Q = zeros(n)
    for i in 1:n
        T[i] = T_known[i] ? T_data[i] : x[col_of_T[i]]
        Q[i] = Q_known[i] ? Q_data[i] : x[col_of_Q[i]]
    end
    return T, Q, A, b, x
end

T_interior(xd, x0, xf, T, Q) =
    0.5 * (T[1] + T[2]) + (xd - x0)/2 * Q[1] + (xd - xf)/2 * Q[2]
```

=== Caso 1 — Dirichlet

$
  T(0) = 100, quad T(1) = 0
  quad => quad T_"exata"(x) = 100 - 100 x, quad Q = -100.
$

```julia
x0, xf = 0.0, 1.0
H, G, L = bem1d_matrices(x0, xf)

T_data = [100.0, 0.0]   # T conhecido nos dois nós
Q_data = [NaN, NaN]     # Q desconhecido
T, Q, A, b, x = solve_bem1d(H, G, T_data, Q_data)
@show T Q               # Q ≈ [-100, -100]

xs = range(x0, xf; length=100)
T_num   = [T_interior(xd, x0, xf, T, Q) for xd in xs]
T_exato = @. 100 - 100 * xs
@show maximum(abs, T_num - T_exato)

plot(xs, T_exato; label="exato", xlabel="x", ylabel="T",
     title="Caso 1 — Dirichlet", lw=2)
plot!(xs, T_num; label="BEM (interior)", ls=:dash, lw=2)
```

=== Caso 2 — mista (um dado em cada ponta)

Pedagogicamente mais limpo do que dois dados no mesmo extremo:

$
  T(0) = 100, quad Q(1) = 0
  quad => quad T equiv 100, quad Q equiv 0.
$

```julia
T_data = [100.0, NaN]
Q_data = [NaN, 0.0]
T, Q, A, b, x = solve_bem1d(H, G, T_data, Q_data)
@show T Q    # T ≈ [100, 100], Q ≈ [0, 0]
```

*Caso 2b (variante).* $T(0)=100$, $Q(0)=0$ também fecha algebricamente e devolve $T_f=100$, $Q_f=0$. Serve para ver que bastam *dois dados independentes* entre os quatro slots — não necessariamente um em cada lado.

== Roteiro mental

+ Operador no domínio $->$ resíduo ponderado.
+ Integrar por partes até o operador cair na função peso.
+ Escolher o peso = SF ($-L^* T^* = delta$).
+ Domínio vira $T(x_d)$; sobram termos só em $Gamma$.
+ Colocar $x_d$ nos nós de contorno $->$ $H T = G Q$.
+ Aplicar CDC $->$ $A x = b$.
+ Interior: reavaliar a identidade com $x_d in Omega$.

== Exercícios

+ *Convecção (Robin) à direita.*
  No extremo direito, $Q_f = h (T_f - T_infinity)$ com $T_infinity = 20$ °C e $h = 3$ (unidades consistentes do modelo 1D). À esquerda, $T_0 = 100$ °C. Como $Q_f$ depende de $T_f$, substitua essa relação *antes* de montar $A x = b$ (a linha correspondente mistura colunas de $H$ e $G$).
  Resolva analítica e numericamente. Compare $T(x)$ e os fluxos nas pontas.

+ *Fonte uniforme (Poisson 1D).*
  Com geração $b$ constante, a identidade ganha um termo de domínio

  $
    - T(x_d)
    = ..."(termos de contorno)..."
    + integral_(x_0)^(x_f) T^*(x, x_d) b dif x.
  $

  Para $b$ constante a integral é analítica. Considere a placa de espessura $L = 2$ cm, $k = 1$ W/(m·K), $b = 1000$ kW/m³, faces a $T_A = 100$ °C e $T_B = 200$ °C. Atente às unidades ($L$ em metros). Compare com

  #image("../assets/apresentacao/placa-calor.png", width: 50%)

  $
    T(x)
    = [(T_B - T_A)/L + b/(2 k) (L - x)] x + T_A.
  $

  *Nota:* a equação forte passa a envolver $b$ e $k$; mantenha a mesma convenção de sinal da SF e do termo de domínio.

+ *Perguntas de checagem (sem código).*
  - Por que a matriz do BEM, em geral, é cheia?
  - O que precisa existir para o BEM “clássico” de um operador linear?
  - Se só $T_0$ e $T_f$ são dados, o sistema devolve o quê?

== Preparar o Julia desta aula

```powershell
winget install julia -s msstore
julia --version
```

```julia
using Pkg
Pkg.add("Plots")
using Plots
```

Gráficos desta aula usam *Plots.jl*. Capítulos seguintes reintroduzem malha, Gmsh e o fluxo de trabalho completo quando a geometria deixar de ser um intervalo.

== Extra (fora da trilha da aula 1) — radiação

Condição não linear de radiação entre superfícies:

$
  q_n = kappa f_s f_epsilon.alt (u^4 - u_R^4)
  = underbrace(kappa f_s f_epsilon.alt (u^2 + u_R^2)(u + u_R), h_r (u))
  (u - u_R),
$

com $kappa approx 5.699 times 10^(-8)$ W/(m²·K⁴) (Stefan–Boltzmann), $0 <= f_s <= 1$ fator de forma e $0 < f_epsilon.alt <= 1$ emissividade. Parece convecção com $h_r$ dependente da própria temperatura: resolve-se o problema linear, atualiza-se $h_r$, itera-se até a variação de $u$ ficar pequena.

Projeto opcional *depois* de dominar Robin linear — não é pré-requisito da aula 1.
