// Interpolação
// markdown2typst + tex2typst

= Interpolação
<interpolacao>

Gráficos deste capítulo usam *CairoMakie* (ver glossário).


== Interpolação polinomial

Suponha que queremos conhecer a população às vezes entre os anos do censo ou estimar populações futuras. Uma técnica é encontrar um polinômio que passa por todos os pontos de dados.

Dado $n$ pontos $(t_1 , y_1), ..., (t_n , y_n)$, onde os $t_i$ são todos distintos, o problema *de interpolação polinomial* é encontrar um polinomial $P$ de grau menor que n tal que $p(t_i) = y_i$ para todos $i$.

O problema de interpolação polinomial tem uma solução única. Uma vez encontrado o polinômio interpolador, ele pode ser avaliado em qualquer lugar para estimar ou prever valores.

=== *Interpolação como um sistema linear*

Dados os dados $(t_i , y_i)$ por $i = 1, ..., n$, buscamos um polinômio

$ p(t) = c_1 + c_2 t + c_3 t^2 + dots.c + c_n t^(n - 1) , $

de tal modo que $y_i = p(t_i)$ para todos $i$. Essas condições são usadas para determinar os coeficientes $c_1 ..., c_n$:

$ c_1 + c_2 t_1 + dots.c + c_(n - 1) t_1^(n - 2) + c_n t_1^(n - 1) = y_1 , \ c_1 + c_2 t_2 + dots.c + c_(n - 1) t_2^(n - 2) + c_n t_2^(n - 1) = y_2 , \ c_1 + c_2 t_3 + dots.c + c_(n - 1) t_3^(n - 2) + c_n t_3^(n - 1) = y_3 , \ dots.v wide \ c_1 + c_2 t_n + dots.c + c_(n - 1) t_n^(n - 2) + c_n t_n^(n - 1) = y_n . $

Essas equações formam um sistema linear para os coeficientes  $c_i :$

$ mat(delim: "[", 1, t_1, dots.c, t_1^(n - 2), t_1^(n - 1);
1, t_2, dots.c, t_2^(n - 2), t_2^(n - 1);
1, t_3, dots.c, t_3^(n - 2), t_3^(n - 1);
dots.v, dots.v, , dots.v, dots.v;
1, t_n, dots.c, t_n^(n - 2), t_n^(n - 1); )
mat(delim: "[", c_1;
c_2;
c_3;
dots.v;
c_n)
=
mat(delim: "[", y_1;
y_2;
y_3;
dots.v;
y_n), $

ou simplesmente, $upright(bold(V)) upright(bold(c)) = upright(bold(y))$. A Matrix $upright(bold(V))$ é de um tipo especial chamado de Vandermonde.

A interpolação polinomial pode, portanto, ser formulada como um sistema linear de equações com uma matriz de Vandermonde.
Criamos dois vetores para dados sobre a população da China. O primeiro tem os anos dos dados do censo e o outro tem a população, em milhões de pessoas.

```julia
year = [1982,2000,2010,2015];
pop = [1008.18, 1262.64, 1337.82, 1374.62];
t = year .- 1980.0
y = pop;
V = [ t[i]^j for i=1:4, j=0:3 ]
#4×4 Matrix{Float64}:
 #1.0   2.0     4.0      8.0
 #1.0  20.0   400.0   8000.0
 #1.0  30.0   900.0  27000.0
 #1.0  35.0  1225.0  42875.0
 c = V \ y 
 #4-element Vector{Float64}:
# 962.2387878787877
 # 24.12775468975476
  #-0.592262049062053
   #0.006843867243867301

using Polynomials
p = Polynomial(c)   
p(2005-1980)         
```

O valor oficial da população para 2005 foi 1303,72, então nosso resultado é bastante bom.

```julia
using CairoMakie
fig = Figure(size=(700, 400))
ax = Axis(fig[1, 1], xlabel="anos desde 1980", ylabel="população (milhões)",
          title="População da China")
scatter!(ax, t, y, label="real")
tt = range(0, 35; length=500)
yy = p.(tt)
lines!(ax, tt, yy, label="interpolante")
axislegend(ax, position=:lt)
fig
```

== Exercícios

+ Suponha que você queira interpolar os pontos (-1,0), (0,1), (2,0), (3,1) e (4,2) por um polinômio de grau o mais baixo possível.
  *(a)* ✍ Qual é o grau máximo necessário desse polinômio?
  *(b)* ✍ Escreva um sistema linear de equações para os coeficientes do polinômio interpolador.
  *(c)* ⌨ Use Julia para resolver numericamente o sistema em (b).
+ *(a)* ✍ Suponha que você quer encontrar um polinômio cúbico $p$ tal que $p(- 1) = - 2$, $p' (- 1) = 1$, $p(1) = 0$, e $p' (1) = - 1$. (Isso é conhecido como um _Interpolador de Hermite._) Escreva um sistema linear de equações para os coeficientes de $p$.
  *(b)* ⌨ Use Julia para resolver o sistema linear na parte (a), e faça um gráfico de $p$ sobre $-1 <= x <= 1$.

= Continuando com interpolação

Dado $n + 1$ pontos distintos $(t_0 , y_0)$, $(t_1 , y_1), ..., (t_n , y_n)$, com $t_0 < t_1 < ... < t_n$ chamados de nós, o problema de interpolação é encontrar uma função $p(x)$, chamada de interpolante, tal que $p(t_k) = y_k$  para $k = 0, ..., n$.

Aqui $t_k$ são os nós e $x$ denota a variável independente contínua.

Os nós de interpolação são numerados de 0 a $n$. Isso é conveniente para nossas declarações matemáticas, mas menos em um idioma como Julia, no qual os índices vetoriais começam com 1. Os índices em um código de computador têm o mesmo significado que os nomeados idênticos nas fórmulas matemáticas e portanto, deve ser incrementado por um sempre que usado em um contexto de indexação.

Os polinômios são o primeiro candidato óbvio a servir como funções interpolando. Eles são fáceis de trabalhar e, em vimos que um sistema linear de equações pode ser usado para determinar os coeficientes de um polinômio que passa por todos os membros de um conjunto de pontos. No entanto, não é difícil encontrar exemplos para os quais a interpolação polinomial leva a resultados inutilizáveis.

Aqui estão alguns pontos que poderíamos considerar ser observações de uma função desconhecida em \[-1,1\].

```julia
n = 5
t = range(-1,1,length=n+1)
y = @. t^2 + t + 0.05*sin(20*t)

using CairoMakie
fig = Figure(); ax = Axis(fig[1, 1]); scatter!(ax, t, y, label="dados"); fig
```

O interpolante polinomial, calculado usando o `fit`, parece muito bom.

```julia
using CairoMakie
p = Polynomials.fit(t, y, n)
xx = range(-1, 1; length=400)
fig = Figure(); ax = Axis(fig[1, 1])
scatter!(ax, t, y, label="dados")
lines!(ax, xx, p.(xx), label="interpolante")
axislegend(ax, position=:lt)
fig
```

Mas agora considere um conjunto diferente de pontos gerados quase exatamente da mesma maneira.

```julia
n = 18
t = range(-1,1,length=n+1)
y = @. t^2 + t + 0.05*sin(20*t)

using CairoMakie
fig = Figure(); ax = Axis(fig[1, 1]); scatter!(ax, t, y); fig
```

Os pontos em si não têm nada de especial. Mas observe o que acontece com o interpolante polinomial.

```julia
using CairoMakie
p = Polynomials.fit(t, y, n)
x = range(-1, 1; length=1000)
fig = Figure(); ax = Axis(fig[1, 1])
scatter!(ax, t, y, label="dados")
lines!(ax, x, p.(x), label="interpolante")
axislegend(ax, position=:lt)
fig
```

Certamente deve haver funções que são mais representativas desses pontos!

== Interpolação por polinômios por partes

Para manter pequenos graus de polinomial enquanto interpolam grandes conjuntos de dados, escolheremos interpolantes dos polinômios por partes. Especificamente, o interpolante p deve ser um polinômio em cada subintervalo $[t_(k - 1) , t_k]$ para  $k = 1, ..., n .$

Geralmente, designamos antecipadamente um grau máximo  para cada parte polinomial de p (x).

```julia
using CairoMakie, Dierckx
p1 = Spline1D(t, y; k=1)
p2 = Spline1D(t, y; k=2)
p3 = Spline1D(t, y; k=3)
xx = range(-1, 1; length=400)
fig = Figure(); ax = Axis(fig[1, 1])
scatter!(ax, t, y, label="dados")
lines!(ax, xx, p1.(xx), label="linear por partes")
lines!(ax, xx, p2.(xx), label="quadrático por partes")
lines!(ax, xx, p3.(xx), label="cúbico por partes")
axislegend(ax, position=:lt)
fig
```

== Exercício

+ Os dois vetores a seguir definem uma forma geometrica.

```julia
x = [ 0,0.51,0.96,1.06,1.29,1.55,1.73,2.13,2.61,
      2.19,1.76,1.56,1.25,1.04,0.58,0 ]
y = [ 0,0.16,0.16,0.43,0.62,0.48,0.19,0.18,0,
      -0.12,-0.12,-0.29,-0.30,-0.15,-0.16,0 ]
```

Podemos considerar x e y como funções de um parâmetro s, com os pontos sendo valores dados em  $s = 0, 1, ..., 15$.

(a) Interpole os pontos usando spline1d e plote a imagem.

(b) Uma desvantagem do resultado na parte (a) é o canto perceptível no lado esquerdo, que corresponde a s = 0 de cima e s = 15 de baixo. Teste adicionar a palavra-chave `periodic = true` à chamada de spline1d e plote o resultado.

== *Estabilidade da interpolação polinomial*

Escolhemos uma função em relação ao intervalo  $[0, 1]$ .

```julia
using CairoMakie
f = x -> sin(exp(2*x))
t = (0:6) ./ 6
y = f.(t)
p = Polynomials.fit(t, y)
xx = range(0, 1; length=400)
fig = Figure(); ax = Axis(fig[1, 1], title="Interpolante equidistante, n=6")
lines!(ax, xx, f.(xx), label="função")
scatter!(ax, t, y, label="nós")
lines!(ax, xx, p.(xx), label="interpolante")
axislegend(ax, position=:lb)
fig
```

Isso parece bom. Queremos rastrear o comportamento do erro à medida que $N$ aumenta. Estimaremos o erro no interpolante contínuo, amostrando-o em um grande número de pontos e tomando a norma máxima.

```julia
using LinearAlgebra
n = 5:5:60;  err = zeros(size(n))
x = range(0,1,length=2001)      # pontos para medir o erro
for (i,n) in enumerate(n)
	t = (0:n)/n                   
	y = f.(t)                     
	p = Polynomials.fit(t, y)
	err[i] = norm((@. f(x) - p(x)), Inf)
end
using CairoMakie
fig = Figure(); ax = Axis(fig[1, 1], xlabel="n", ylabel="max error",
          yscale=log10, title="Erro de interpolação para nós equidistantes")
scatterlines!(ax, n, err)
fig
```

O erro diminui inicialmente como seria de esperar, mas começa a crescer. Ambas as fases ocorrem a taxas exponenciais em $n$, ou seja,  $O(k^n)$, aparecendo linear em um gráfico semi-log.

== Fenômeno de Runge

A decepcionante perda de convergência  é um sinal de mau condicionamento devido ao uso de nós igualmente espaçados.

```julia
using CairoMakie, LinearAlgebra
f = x -> 1/(x^2 + 16)
xx = range(-1, 1; length=400)
fig = Figure(); ax = Axis(fig[1, 1], title="Função teste")
lines!(ax, xx, f.(xx)); fig
```

Essa função possui infinitamente muitas derivadas contínuas em toda a linha real e parece fácil de aproximar em $[- 1, 1]$. Começamos fazendo interpolação polinomial equispacada para alguns pequenos valores de $n$.

```julia
using CairoMakie
x = range(-1, 1; length=2501)
fig = Figure(); ax = Axis(fig[1, 1], xlabel="x", ylabel="|f-p|", yscale=log10,
          title="Erro para graus baixos", limits=(nothing, (1e-20, 1)))
for n in 4:4:12
    tt = range(-1, 1; length=n+1)
    pp = Polynomials.fit(tt, f.(tt))
    lines!(ax, x, abs.(f.(x) .- pp.(x)), label="grau $n")
end
axislegend(ax); fig
```

A convergência até agora parece bastante boa, embora não seja uniformemente. No entanto, observe o que acontece à medida que continuamos aumentando o grau.

```julia
using CairoMakie
fig = Figure(); ax = Axis(fig[1, 1], xlabel="x", ylabel="|f-p|", yscale=log10,
          title="Erro para graus altos", limits=(nothing, (1e-20, 1)))
for n in @. 12 + 15*(1:3)
    tt = range(-1, 1; length=n+1)
    pp = Polynomials.fit(tt, f.(tt))
    lines!(ax, x, abs.(f.(x) .- pp.(x)), label="grau $n")
end
axislegend(ax); fig
```

A convergência no meio não pode ficar melhor do que a precisão da máquina em relação aos valores da função. Portanto, a manutenção da lacuna crescente entre o centro e as extremidades empurra as curvas de erro exponencialmente rapidamente nas extremidades, destruindo a convergência.

A observação da instabilidade é conhecida como o fenômeno de Runge. O fenômeno de Runge é uma instabilidade manifestada quando os nós do interpolante são igualmente espaçados e o grau do polinomial aumenta.
Significativamente, a convergência observada  é estável dentro de uma parte do intervalo. Ao redistribuir os nós de interpolação, sacrificaremos um pouco da convergência na parte do meio, a fim de melhorá-la perto das extremidades e resgatar o processo globalmente.

Uma família de nós que fornece convergência estável para a interpolação polinomial é os pontos  Chebyshev do segundo tipo definido por:

$ t_k = - cos((k pi)/n), wide k = 0, ..., n . $

```julia
using CairoMakie
x = range(-1, 1; length=2001)
fig = Figure(); ax = Axis(fig[1, 1], xlabel="x", ylabel="|f-p|", yscale=log10,
          title="Erro com os pontos de Chebyshev", limits=(nothing, (1e-20, 1)))
for n in [4, 10, 16, 40]
    tt = [-cos(pi * k / n) for k in 0:n]
    pp = Polynomials.fit(tt, f.(tt))
    lines!(ax, x, abs.(f.(x) .- pp.(x)), label="grau $n")
end
axislegend(ax); fig
```

A partir do grau 16, o erro está dentro da precisão de máquina e ele permanece lá à medida que $n$ aumenta.

=== Exercício

+ Para cada caso, calcule o interpolante polinomial usando  $n$ nós de Chebyshev do segundo tipo em  $[- 1, 1]$  por  $n = 4, 8, 12, ..., 60$ . Em cada valor de $n$, calcule o erro (ou seja,  $max | p(x) - f(x) |$  avaliado em 4000 valores de $x$. Usando uma escala log-linear, plote o erro em função de  $n$ e, em seguida, determine uma boa aproximação à constante $k$ de $O(n^(-k))$.
(a) $f(x) = 1 \/ (25 x^2 + 1)$
(b) $f(x) = tanh(5 x + 2)$
(c)  $f(x) = cosh(sin x)$
(d)  $f(x) = sin(cosh x)$

== Integração numérica

A primitiva de $e^x$ é simples, isso torna a avaliação de $integral_(-1)^1 e^x d x$  pelo teorema fundamental trivial.

```julia
using FastGaussQuadrature, LinearAlgebra
exato = exp(1)-exp(-1)

x, w = gausslegendre(3)
#([-0.7745966692414834, 0.0, 0.7745966692414834], [0.5555555555555556, 0.8888888888888888, 0.5555555555555556])

 f(x) = exp(x)

In = dot(w, f.(x))
exato-In
```

A abordagem numérica é muito robusta. Por exemplo, $e^(sin x)$ não tem primitiva conhecida. Mas numericamente, sua integral não é mais difícil de ser calculada.

```julia
f(x) = exp(sin(x))
In = dot(w, f.(x))
```

Quando você olha para os gráficos dessas funções, o que é notável é que uma dessas áreas é muito simples, enquanto o outro é analiticamente muito difícil. Do ponto de vista numérico, eles são praticamente o mesmo problema.

```julia
using CairoMakie
xx = range(-1, 1; length=400)
fig = Figure(size=(700, 500))
ax1 = Axis(fig[1, 1], xlabel="x", ylabel="exp(x)", limits=(nothing, (0, 2.7)))
ax2 = Axis(fig[2, 1], xlabel="x", ylabel="exp(sin(x))", limits=(nothing, (0, 2.7)))
band!(ax1, xx, zero(xx), exp.(xx)); lines!(ax1, xx, exp.(xx))
band!(ax2, xx, zero(xx), exp.(sin.(xx))); lines!(ax2, xx, exp.(sin.(xx)))
fig
```

A integração numérica, que também passa pelo nome mais antigo quadratura, é executada pela combinação de valores do integrando amostrados nos nós. Nesta seção, assumiremos nós igualmente espaçados:

$ t_i = a + i h, quad h = (b - a)/n , wide i = 0, ..., n . $

A Integração numérica consiste em uma lista de pesos $w_0 , ..., w_N$ escolhidos de modo que:

$ integral_a^b f(x) thin d x approx h sum_(i = 0)^n w_i f(t_i) = h [w_0 f(t_0) + w_1 f(t_1) + dots.c w_n f(t_n)] . $

Uma maneira direta de derivar fórmulas de integração é encontrar um interpolante e operar exatamente nele.

=== Regra do trapézio

Uma das fórmulas de integração mais importantes resulta da integração do interpolante linear por partes. Geometricamente, a fórmula trapezoidal de áreas de trapezoides que se aproximam da região sob a curva y = f (x)

#image("../assets/interpolacao/trapezio.svg", width: 80%)

Usando áreas de triângulos, é trivial derivar que

$ w_i = cases(1 comma & i = 1 comma ... comma n - 1 comma, 1/2 comma & i = 0 comma n .) $

```julia
"""
    trapezoidal(f,a,b,n)

Aplique a fórmula de integração do trapezoidal no intervalo [`a`,` b`], quebrado em ppartes iguais.Retorna a estimativa, um vetor de nós e um vetor de valores do integrando nos nós.
"""

function trapezoidal(f,a,b,n)
    h = (b-a)/n
    t = range(a,b,length=n+1)
    y = f.(t)
    T = h * ( sum(y[2:n]) + 0.5*(y[1] + y[n+1]) )
    return T,t,y
end
```

Vamos aproximar $f(x) = e^(sin 7 x)$ no intervalo $[0, 2]$.

```julia
f = x -> exp(sin(7*x));
a = 0;  b = 2;

Integral = 2.6632197827615394 #exato
T,t,y = trapezoidal(f,a,b,40)
@show (T,Integral-T);

n = [ 10^n for n in 1:5 ]
err = []
for n in n
    T,t,y = trapezoidal(f,a,b,n)
    push!(err,Integral-T)
end

foreach(args -> println(args[1], "	", args[2]), zip(n, err))
```

Cada aumento em um fator de 10 em $n$ reduz o erro em um fator de cerca de 100, o que é consistente com a convergência de segunda ordem.

=== Quadratura de Gauss

Vamos considerar a fórmula de integração numérica genérica:

$ integral_(-1)^1 f(x) d x approx sum_(k = 1)^n w_k f(t_k) = Q_n [f], $

Como existem  $n$  nós e $n$  pesos disponíveis para escolher, parece plausível esperar  que a quadratura consiga representar exatamente os polinômios de grau $m = 2 n - 1$ , e essa intuição acaba sendo correta.

Vamos testar isso na integral:

$ integral_(-1)^1 1/(1 + 4 x^2) d x = arctan(2) . $

```julia
f = x->1/(1+4*x^2);
exato = atan(2);

n = 8:4:96
errT = zeros(size(n))
errG = zeros(size(n))
for (k,n) in enumerate(n)
  errT[k] = abs(exato - trapezoidal(f,-1,1,n)[1])
  x, w = gausslegendre(n)
  errG[k] = abs(exato - dot(w, f.(x)))
end

errT[iszero.(errT)] .= NaN
errG[iszero.(errG)] .= NaN
using CairoMakie
fig = Figure(); ax = Axis(fig[1, 1], xlabel="nós", ylabel="erro", yscale=log10,
          title="integração numérica", limits=(nothing, (1e-16, 1)))
scatterlines!(ax, collect(n), errT, label="trapézio")
scatterlines!(ax, collect(n), errG, label="Gauss-Legendre")
axislegend(ax); fig
    
    
```

E agora com uma integral com um integrando mais pontudo:

$integral_(-1)^1 1/(1 + 16 x^2) thin d x = 1/2 arctan(4) .$

```julia
f = x->1/(1+16*x^2);
exato = atan(4)/2;

n = 8:4:96
errT = zeros(size(n))
errG = zeros(size(n))
for (k,n) in enumerate(n)
  errT[k] = abs(exato - trapezoidal(f,-1,1,n)[1])
  x, w = gausslegendre(n)
  errG[k] = abs(exato - dot(w, f.(x)))
end

errT[iszero.(errT)] .= NaN
errG[iszero.(errG)] .= NaN
using CairoMakie
fig = Figure(); ax = Axis(fig[1, 1], xlabel="nós", ylabel="erro", yscale=log10,
          title="integração numérica", limits=(nothing, (1e-16, 1)))
scatterlines!(ax, collect(n), errT, label="trapézio")
scatterlines!(ax, collect(n), errG, label="Gauss-Legendre")
axislegend(ax); fig
```

=== Exercícios

+ Usando a mudança de variável:     $z = phi.alt(x) = a + (b - a) ((x + 1))/2$ a integral pode ser reescrita como:     $integral_a^b f(z) thin d z = (b - a)/2 integral_(-1)^1 f(phi.alt(x)) thin d x .$ Sabendo disso, use a quadratura de Gauss para integrar     $integral_(pi \/ 2)^pi x^2 sin 8 x thin d x = - (3 pi^2)/32 .$ Mostre resultados para diferentes valores de $n$ até que se obtenha uma convergência de 10 digitos.
+ Integre  numericamente a  função  $f(x) = - x log(| x - 1 \/ 2 |)$ no intervalo $[0, 1 \/ 2]$ usando  a quadratura de Gauss. Use 5, 10, 15 e 20 pontos de Gauss.
  $ I_a = integral_0^(1 \/ 2) - x log(| x - 1/2 |) d x = 1/16 (3 + 2 log 2) = 0, 274143 $

== Integrais singulares ou quasi-singulares

A última função integrada tem uma singularidade quando $x = 1 \/ 2$. Mesmo com essa singularidade, e tendendo a infinito nesse ponto, essa função pode ser integrada analiticamente. Numericamente é interessante usar técnicas especiais para tratar de integrais desses tipo.

Considere uma integral do tipo:

$ I_s = integral_(-1)^1 (g(xi))/((xi - a)^2 + b^2) dif xi $

Uma ideia para avaliar numericamente integrais singulares($b = 0$) ou quasi-singulares($b != 0$) é  encontrar uma transformação cujo Jacobiano seja zero ou bem próximo de zero no ponto singular. Por exemplo, para integrais quasi-singulares podemos usar:

$ xi = a + b sinh(mu s - eta) $

onde

$ mu = 1/2 {"arcsinh" ((1 + a)/b) + "arcsinh" ((1 - a)/b)} \ eta = 1/2 {"arcsinh" ((1 + a)/b) - "arcsinh" ((1 - a)/b)} $

e seu jacobiano é dado por:

$ (dif xi)/(dif s) = b mu cosh(mu s - eta) $

Aplicando essa transformação, obtém-se:

$ I_s = integral_(-1)^1 (g(s(xi)))/((s(xi) - a)^2 + b^2) (dif xi)/(dif s) dif s $

essa transformação suaviza o integrando e acelera a convergência da quadratura.

```julia

function sinhtrans(u, a, b)
    μ = 1 / 2 * (asinh((1 + a) / b) + asinh((1 - a) / b))
    η = 1 / 2 * (asinh((1 + a) / b) - asinh((1 - a) / b))

    x = a .+ b * sinh.(μ * u .- η)
    J = b * μ * cosh.(μ * u .- η)
    x, J
end
```

```julia
using CairoMakie
x = range(-1, 1; length=100)
xt, jac = sinhtrans(x, 0.5, 0.01)
fig = Figure(size=(700, 500))
ax1 = Axis(fig[1, 1], xlabel="s", ylabel="ξ")
ax2 = Axis(fig[2, 1], xlabel="s", ylabel="dξ/ds")
lines!(ax1, collect(x), xt); lines!(ax2, collect(x), jac)
fig
```

aplicando essa técnica para a integral do último exemplo pode-se observar uma redução significativa do erro para a mesma quantidade de pontos.

```julia
f = x->1/(1+16*x^2);
exato = atan(4)/2;

n = 8

 x, w = gausslegendre(n)
 xt,J= sinhtrans(x, 0, 1/4)
 errG = abs(exato - dot(w, f.(x)))
 errGt = abs(exato - dot(w, J.*f.(xt)))
```

=== Exercício

+ Integre, outra vez , numericamente a  função  $f(x) = - x log(| x - 1 \/ 2 |)$ no intervalo $[0, 1 \/ 2]$ usando  a quadratura de Gauss mas agora usando a transformada de Monegato com graus 3, 4 e 5. Use 5, 10, 15 e 20 pontos de Gauss e compare com os resultados.

```julia
"""
transformada de Monegato de grau q
singularidade em s0
"""
function Monegato(t, s0, q=5.0)
    δ = 2^(-q) * ((1 + s0)^(1 / q) + (1 - s0)^(1 / q))^q
    t0 = ((1 + s0)^(1 / q) - (1 - s0)^(1 / q)) / ((1 + s0)^(1 / q) + (1 - s0)^(1 / q))
    s = s0 .+ δ * (t .- t0) .^ q
    ds = q * δ * (t .- t0) .^ (q - 1)
    s, ds
end
```

+ Use 10 pontos de Gauss para calcular as integrais abaixo. Compare com a solução analítica. Quais integrais tiveram os maiores erros? Por que?

$ mat(delim: #none, integral_(-2)^(+3) (x^6 - 2 x^5 + 7) dif x, 144.04761904761904; integral_0^(+3) 2 ln(x + 1) dif x, 5.090354888959125; integral_0^(+3) ln(x) dif x, 0.2958368660043291; integral_0^pi lr(\(sin(x)) times cos(2 x)) dif x, -2/3 = - 0.66666666666667; integral_0^pi (sin(3 x) dot.op cos(2 x) + 2 x^3 + 3 x sin(x)) dif x, 59.3293234777059; integral_0^1 1/(x + 1) dif x, ln(2) = 0.6931471805599453; integral_0^1 1/(x + 0.1) dif x, ln(11) = 2.3978952727983707; ) $

== Desafio

Considere o resfriamento de uma aleta circular por meio de transferência de calor por convecção ao longo de seu comprimento. A convecção dá origem a uma perda de calor ou termo de sumidouro dependente da temperatura na equação governante. Mostrada na Figura está uma aleta cilíndrica com área de seção transversal uniforme A. A base está a uma temperatura de 100C (TB) e a extremidade direita está isolada. A aleta está exposta a uma temperatura ambiente de 20C. A transferência de calor unidimensional nesta situação é governada por$dif/(dif x) (k A (dif T)/(dif x)) - h P(T - T_infinity) = 0$

onde h é o coeficiente de transferência de calor por convecção, P o perímetro, k a
condutividade térmica do material e T— a temperatura ambiente.
Calcule usando BEM a distribuição de temperatura ao longo da aleta e compare os resultados com a solução analítica fornecida por $(T - T_infinity)/(T_B - T_infinity) = (cosh [n(L - x)])/(cosh(n L))$

#image("../assets/interpolacao/aleta.png", width: 80%)
