// Indo para 2D
// markdown2typst + tex2typst

= Indo para 2D
<indo-para-2d>

== Cálculo do perímetro de figuras planas

Uma vez que a integral ao longo do contorno no $s$ pode ser bastante difícil, ou mesmo impossível, uma estratégia para o cálculo do perímetro é a divisão do contorno em uma soma de pequenos
pedaços $s_1 , s_2 , . . ., s_n$, ou seja, $s = sum_(i = 1)^n s$, onde $n$ é o número de pedaços em que o contorno foi dividido. Uma vez
que estes pedaços podem ter uma forma qualquer, cada pedaço $s_i$ será aproximado por uma forma conhecida. Por simplicidade, esta forma é quase sempre dada por um polinômio (linha reta, parábola, etc).

#image("../assets/indo-para-2d/perimetro.png", width: 80%)

Dessa maneira, cada pedaço $s_1 , s_2 , . . ., s_n$ é aproximado por formas conhecidas $Gamma_1 , Gamma_2 , . . ., Gamma_n$, chamados Elementos de Contorno.

Considere agora que os elementos de contorno $Gamma_i$ sejam parabólicos, ou seja, são descritos por polinômios de 2ª ordem (equação de uma parábola). Desta forma, são necessários 3 pontos de $Gamma_i$ para se definir a parábola. Estes pontos são dados por $(x_1 , y_1)$, $(x_2 , y_2)$ e $(x_3 , y_3)$, que correspondem respectivamente a $xi = - 1$, $xi = 0$ e $xi = 1$.

#image("../assets/indo-para-2d/elemento-parabolico.png", width: 80%)

Criando uma função parabólica para relacionar $x$ com $xi$, tem-se: $x = a xi^2 + b xi + c$ sendo que:

$ x(xi = - 1) = x_1 = a(- 1)^2 + b(- 1) + c \ x(xi = 0) = x_2 = a(0)^2 + b(0) + c \ x(xi = + 1) = x_3 = a(+ 1)^2 + b(+ 1) + c $

resolvendo esse sistema obtem-se $x = N_1 x_1 + N_2 x_2 + N_3 x_3$ onde $N_1 = xi/2 (xi - 1)$, $N_2 = (1 - xi)(1 + xi)$ e $N_3 = xi/2 (xi + 1)$ são as funções de forma quadráticas contínuas.

#image("../assets/indo-para-2d/funcoes-forma.png", width: 80%)

A derivada de $x$ em relação a $xi$ é dada por $(d x)/(d xi) = (d [N_1 (xi)])/(d xi) x_1 + (d [N_2 (xi)])/(d xi) x_2 + (d [N_3 (xi)])/(d xi) x_3$ onde $(d N_1)/(d xi) = xi - 1/2 ,$
$(d N_2)/(d xi) = - 2 xi$
e $(d N_3)/(d xi) = xi + 1/2$
.

O comprimento $Gamma$ do perímetro da figura é então dado por:

$Gamma = sum_(i = 1)^n integral_(-1)^1 (d Gamma)/(d xi) d xi = sum_(i = 1)^n integral_(-1)^1 sqrt(((d x)/(d xi))^2 + ((d y)/(d xi))^2) d xi = sum_(i = 1)^n [integral_(-1)^1 J(xi) d xi]$ em que a função $J(xi) = sqrt(((d x)/(d xi))^2 + ((d y)/(d xi))^2)$ é chamada jacobiano da transformação. Usando quadratura de Gauss com $p$ pontos de Gauss, conclui-se que $Gamma = sum_(i = 1)^n [sum_(k = 1)^p w_k J(xi_k)]$.

#link("https://1drv.ms/f/s!AmfyGvdmTYongqYkwQK1hIXEROJbEA?e=omRWUR")[código para propriedade geométrica]

```julia
include("dad.jl")
include("format.jl")
include("calcula_PropGeom.jl")
include("calc_fforma.jl")
using FastGaussQuadrature,GLMakie

PONTOS,SEGMENTOS,MALHA=dad_1(10,1) #Arquivo de entrada de dados

NOS,ELEM=format_dad(PONTOS,SEGMENTOS,MALHA)# formata os dados (cria as
  # matrizes NOS e ELEM a partir das matrizes PONTOS, SEGMENTOS e MALHA)
display(mostra_geo(SEGMENTOS,PONTOS,ELEM,NOS))

Perimetro,Area,xbarra,ybarra=calcula_PropGeom(ELEM,NOS,4,4); # Calcula as propriedades geom�tricas de figuras plana
println("Perimetro calculado numericamente: $Perimetro")
println("Área calculada numericamente: $Area")
println("Centróide calculado numericamente: ( $xbarra, $ybarra)")
```

#image("../assets/indo-para-2d/geo-exemplo.png", width: 80%)

== Área

A integral de uma função $f(x, y)$ sobre a área  $Omega$ de uma figura plana é dada pela integral $I = integral_Omega f(x, y) d x d y .$ A integral pode ser escrita, em coordenadas polares,
como: $I = integral_Omega f(x(rho, theta), y(rho, theta)) rho d rho d theta = integral_theta integral_0^r f(x(rho, theta), y(rho, theta)) rho d rho d theta$  definindo $F$ como $F(rho, theta) = integral_0^r f(x(rho, theta), y(rho, theta)) rho d rho$, pode-se escrever:$I = integral_theta F(rho, theta) d theta$.

como $cos alpha = (r (d theta)/2)/((d Gamma)/2)$

podemos escrever: $d theta = (arrow(n) . arrow(r))/r d Gamma$ e finalmente $I = integral_Gamma F (arrow(n) . arrow(r))/r d Gamma$.

Dividindo $Gamma$ em uma soma de pequenos pedaços do contorno:

$I = sum_(i = 1)^(N E) integral_(Gamma_i) F (arrow(n) dot.op arrow(r))/r d Gamma$.

#image("../assets/indo-para-2d/area-polar.png", width: 80%)

#image("../assets/indo-para-2d/area-contorno.png", width: 80%)

=== Cálculo do vetor normal $arrow(n)$

#image("../assets/indo-para-2d/normal-1.png", width: 80%)

#image("../assets/indo-para-2d/normal-2.png", width: 80%)

$arrow(S)$= vetor tangente ao contorno $Gamma$. $arrow(S) = d x arrow(i) + d y arrow(j)$ vetor não unitário

$arrow(s)$= Vetor unitário tangente ao contorno $Gamma$

$arrow(s) = arrow(S)/(| arrow(S) |) = (d x arrow(i) + d y arrow(j))/sqrt(d x^2 + d y^2)$

Dividindo tudo por $d xi$, tem-se:

$arrow(s) = ((d x)/(d xi))/sqrt(((d x)/(d xi))^2 + ((d y)/(d xi))^2) arrow(i) + ((d y)/(d xi))/sqrt(((d x)/(d xi))^2 + ((d y)/(d xi))^2) arrow(j)$

$arrow(s) = ((d x)/(d xi))/((d Gamma)/(d xi)) arrow(i) + ((d y)/(d xi))/((d Gamma)/(d xi)) arrow(j) = s_x arrow(i) + s_y arrow(j)$

$arrow(n)$ = vetor unitário normal ao contorno $Gamma$ apontando para fora do domínio $Omega$

$arrow(n) = n_x arrow(i) + n_y arrow(j)$
Temos duas possibilidades, uma vez que $arrow(s)$ é unitário:

$n_x = - s_y quad
n_y = s_x =>(- s_y) s_x + (s_x) s_y = 0$

$n_x = s_y quad
n_y = - s_x =>(s_y) s_x + (- s_x) s_y = 0$
O vetor $arrow(s)$ é sempre definido percorrendo $Gamma$ no sentido anti-horário, se $Gamma$ for um contorno externo e horário, caso $Gamma$ seja um contorno interno. Portanto, a segunda possibilidade é escolhida de modo que a normal esteja apontada para fora do domínio.

=== Exercícios

#image("../assets/indo-para-2d/exercicio-1.png", width: 80%)

$I_x = (a^4)/96 (9 sqrt(3) - 2 pi)$

#image("../assets/indo-para-2d/exercicio-2.png", width: 80%)

$I_x = 2 dot.op 10^4 pi - (20^2 dot.op pi)/2 (80/(3 dot.op pi))^2
+ ((20^2 pi)/2) (15 + 80/(3 pi))^2$

Calcule, modificando o programa propgeo, os momentos de inércia em relação ao eixo x para as duas figuras e compare com o resultado analítico.

- Gravações
  #link("https://youtu.be/TZvHS2JzKOs")[https://youtu.be/TZvHS2JzKOs]
  #link("https://youtu.be/08vADm8JguE")[https://youtu.be/08vADm8JguE]
  #link("https://youtu.be/O14y7y9ZCLs")[https://youtu.be/O14y7y9ZCLs]
  #link("https://youtu.be/YbKmFAwBbwo")[https://youtu.be/YbKmFAwBbwo]
