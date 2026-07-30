// Laplace 2D
// markdown2typst + tex2typst

= Laplace 2D
<laplace-2d>

- Código
  #link("https://github.com/l-s-campos/BEM")[https://github.com/l-s-campos/BEM]
  faça o download  e rode de dentro da pasta:
  ```julia
  using Pkg
  
  Pkg.activate(pwd())
  
  Pkg.instantiate()
  ```
  \@which

= Formulação

A equação de Laplace é dada por:

$ nabla^2 T = 0, $

usando o método dos resíduos ponderados e aplicando a segunda identidade de Green:

$ integral_Omega (v nabla^2 u - u nabla^2 v) d Omega = integral_Gamma (v (partial u)/(partial n) - u (partial v)/(partial n)) d s $

obtém-se a equação integral de contorno:

$ 0.5 T(x_d , y_d) = integral_Gamma T q^* d s - integral_Gamma T^* q d s, $

onde $T^*$ e $q^*$ são as soluções fundamentais:

$ T^* = (-1)/(2 pi k) ln r $

$ q^* = 1/(2 pi r^2) [(x - x_d) n_x + (y - y_d) n_y] . $

Podemos representar a temperatura e fluxo em um elemento descontinuo com  $m$  nós como:
$T = N_1 T_1 + N_2 T_2 + N_3 T_3 + ... + N_m T_m$
$q = N_1 q_1 + N_2 q_2 + N_3 q_3 + ... + N_m q_m$
onde estamos aproximando nossa distribuição de temperatura e fluxo no elemento por uma função polinomial de  grau  $(m - 1)$ .

Discretizando em $n$ elementos de contorno descontínuos:

$0.5 T(x_d , y_d) = sum_(j = 1)^(n_(e l e m)) [integral_(Gamma_j) T q^* d Gamma] - sum_(j = 1)^(n_(e l e m)) [integral_(Gamma_j) T^* q d Gamma]$

Usando a representação de   $T$ e $q$  na equação acima temos:

$ 0.5 T(x_d , y_d) = sum_(j = 1)^(n_(e l e m)) {integral_(Gamma_j) [mat(delim: #none, N_1, N_2, N_3, dots.c, N_m)] [mat(delim: #none, T_1; T_2; T_3; dots.v; T_m)]_j q^* d Gamma} \ -
sum_(j = 1)^(n_(e l e m)) {integral_(Gamma_j) T^* [mat(delim: #none, N_1, N_2, N_3, dots.c, N_m)] [mat(delim: #none, q_1; q_2; q_3; dots.v; q_m)]_j d Gamma} $

Repetindo essa equação para cada diferente nó podemos montar um sistema matricial:

$ H T = G q $

== Método indireto para o cálculo da diagonal da matriz $H$

A diagonal da matriz $H$ contém uma singularidade forte que dificulta o cálculo numérico desses termos. Uma alternativa não faz a integração de maneira explícita mas usa uma propriedade da matriz $H$ decorrente da modelagem de um corpo sob temperatura constante. Sem perder a generalidade, considere que todos os nós de um corpo encontre-se com a temperatura $T = 1$. Neste caso, o fluxo será nulo em todos os nós, ou seja, $q = 0$  em todos os nós. Desta forma, a equação matricial é reescrita como $H {1} = G {0}$.

Daí, os termos da diagonal da matriz $[H]$ pode ser calculado da seguinte forma:

$H_(i i) = - sum_(j = 1)^N H_(i j) , " com " i != j, " para " i = 1, 2, . . ., N,$

uma vez que todos os termos de fora da diagonal são integrais regulares e já foram previamente calculados.

Um procedimento parecido, usando uma outra distribuição de temperatura, pode ser feito com a diagonal da matriz $G$ mas isso normalmente não é necessário por se tratar de uma singularidade fraca

== Exemplo

A fim de ilustrar como se aplica as condições de contorno e se calcula as variáveis desconhecidas será analisado um problema de condução de calor unidirecional com uma discretização de um elemento por lado.

#image("../assets/laplace-2d/exemplo-unidirecional.png", width: 80%)

As equações obtidas podem ser escritas na forma matricial, como:

$ (mat(delim: #none, H_11, H_12, H_13, H_14; H_21, H_22, H_23, H_24; H_31, H_32, H_33, H_34; H_41, H_42, H_43, H_44; ))
(mat(delim: #none, macron(T_1); T_2; macron(T_3); T_4))
=
(mat(delim: #none, G_11, G_12, G_13, G_14; G_21, G_22, G_23, G_24; G_31, G_32, G_33, G_34; G_41, G_42, G_43, G_44; ))
(mat(delim: #none, q_1; macron(q_2); q_3; macron(q_4))) $

onde $macron(T)$ e $macron(q)$ são termos conhecidos.

Separando os termos conhecidos dos desconhecidos:

$ (mat(delim: #none, -G_11, H_12, -G_13, H_14; -G_21, H_22, -G_23, H_24; -G_31, H_32, -G_33, H_34; -G_41, H_42, -G_43, H_44; ))(mat(delim: #none, q_1; T_2; q_3; T_4)) = (mat(delim: #none, -H_11, G_12, -H_13, G_14; -H_21, G_22, -H_23, G_24; -H_31, G_32, -H_33, G_34; -H_41, G_42, -H_43, G_44; ))(mat(delim: #none, macron(T_1); macron(q_2); macron(T_3); macron(q_4))) $

Assim, pode-se escrever $A x = b$ e resolvendo o sistema linear calcula-se os valores das variáveis desconhecidas.

== Pontos internos

A equação integral para pontos internos é ligeiramente modificada:

$ T(x_d , y_d) = integral_Gamma T q^* d s - integral_Gamma T^* q d s $

Aqui temos um detalhe importante, uma vez conhecida $T$ e $q$  no contorno a temperatura pode ser calculada em qualquer ponto interno sem ser necessário resolver algum sistema linear.

== Exercícios

erros

+ Resolva esse problema com elementos lineares, quadráticos e cúbicos e calcule o erro médio, erro máximo e a norma l2 do erro no contorno para diferentes discretizações. Faça um gráfico comparando a convergência dos dos 3 tipos de elementos. A solução analítica é dada por:

$ u(theta) &= theta/pi \ q(x) &= - 1/(pi x) quad "em" quad y = 0 $

#image("../assets/laplace-2d/setor-circular.png", width: 80%)

+ Analise o seguinte problema:

Analise uma placa com as seguintes condições de contorno:

- $T(x = 0) = 0$
- $T(x = 1) = cos(pi y)$
- $q(y = 0) = q(y = 1) = - k (partial u)/(partial n) = 0$

#image("../assets/laplace-2d/placa-mista.png", width: 80%)

A solução analítica para este problema é dada por:

$T^(a n) = sinh(pi x) cos(pi y) \/ sinh(pi)$

Analise o problema usando diferentes números de elementos e faça uma tabela mostrando o erro percentual em um ponto interno de coordenadas $(x, y) = (sqrt(2) \/ 2, sqrt(2) \/ 2)$  para os diversos casos analisados.

+ Faça um mapa de cor da distribuição de temperatura em uma placa com dimensões e condições de contorno mostradas na figura.

#image("../assets/laplace-2d/placa-mapa.png", width: 80%)

== Desafio

Baseado nesse #link("https://onlinelibrary.wiley.com/doi/epdf/10.1002/fld.1650030504")[artigo] calcule o coeficiente de sustentação de um perfil NACA usando o BEM.

#link("https://youtu.be/_G4yNayAPPE")[gravação]
