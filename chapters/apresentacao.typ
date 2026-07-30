// Apresentação
// markdown2typst + tex2typst

= Apresentação
<apresentacao>

== Método dos elementos de contorno

O método dos elementos de contorno é bem conhecido entre engenheiros e cientistas. Ele tem demonstrado sua superioridade em relação a outros métodos numéricos, especialmente quando usado para modelar uma aplicação apropriada. Apesar da popularidade do método dos elementos de contorno, ele não é tão comum entre os engenheiros quanto o método dos elementos finitos. As razões para isso podem ser resumidas da seguinte forma:

1- A complexidade da formulação matemática.
2- A falta de códigos simples.
3- A falta de  cursos de elementos de contorno entre estudantes de graduação.
4- A dificuldade no tratamento de alguns modelos numéricos, como singularidade.
5- A dificuldade em modificar programas de elementos de contorno em relação aos desenvolvidos usando elementos finitos.
6- A falta de versatilidade dos códigos de elementos de contorno.
7- A mudança de estratégia de modelagem de elementos finitos para elementos de contorno.

== Vantagens e desvantagens

O Método dos Elementos de Contorno (BEM), como qualquer outro método numérico, tem suas vantagens e desvantagens. As vantagens do método dos elementos de contorno são as seguintes:

1-  Apenas o contorno do problema precisa ser discretizado, o que leva a uma fácil preparação de dados e menores requisitos de computação.
2- O tratamento exato de domínios infinitos e semi-infinitos.
3- As incógnitas em locais internos são calculadas na fase de pós-processamento, o que simplifica qualquer procedimento de otimização.
4- Resultados precisos no caso de concentrações de tensões devido a fissuras ou cargas concentradas.

Por outro lado, as desvantagens do método dos elementos de contorno são as seguintes:
1- As matrizes do sistema são não simétricas e totalmente preenchidas.
2- As soluções fundamentais nem sempre são fáceis de obter.
3- A dificuldade em tratar estruturas delgadas.
4- A discretização do domínio necessária no caso de aplicações não lineares.

Desafios que precisam ser enfrentados para que o método tenha uma maior adesão:

1- Uma visão mais profunda dos aspectos matemáticos e numéricos do método.
2- Um método sistemático para a derivação das soluções fundamentais e particulares.
3- Fórmula de integração estável.
4- Programas de uso geral e pequenos programas disponíveis para engenheiros.
5- Acoplamento entre elementos de contorno e elementos finitos.

#image("../assets/apresentacao/bem-overview.jpeg", width: 80%)

== Por que JULIA?

A escolha da linguagem JULIA é justificada por várias razões. JULIA é uma linguagem de programação de alto desempenho, projetada especificamente para computação científica e análise numérica. Suas bibliotecas e funcionalidades são otimizadas para realizar cálculos complexos de forma eficiente, o que é crucial para a implementação de métodos como o Método dos Elementos de Contorno (BEM). Além disso, JULIA possui uma sintaxe simples e intuitiva, o que facilita a escrita e a leitura de códigos, tornando o processo de desenvolvimento mais ágil.

*Desempenho:* JULIA é projetada para ter um desempenho superior em cálculos numéricos e científicos, superando Python em termos de velocidade. Isso se deve à sua capacidade de compilar código em tempo de execução, o que resulta em execução mais rápida. Python, por outro lado, pode ser mais lento devido à sua natureza interpretada, embora existam bibliotecas como Cython que ajudam a melhorar o desempenho.

Normalmente, uma linguagem é utilizada para o desenvolvimento rápido e prototipagem, enquanto a outra é usada para obter desempenho otimizado. Por exemplo, um cientista pode usar Python para escrever a lógica do código devido à sua simplicidade e riqueza de bibliotecas, mas precisará reescrever partes críticas em C ou Fortran para alcançar a velocidade necessária.

Essa abordagem tem várias desvantagens:

+ *Complexidade de Manutenção:* Manter código em duas linguagens diferentes pode ser complexo e propenso a erros, especialmente quando as mudanças precisam ser sincronizadas entre as duas versões.
+ *Curva de Aprendizado:* Exige que os desenvolvedores sejam proficientes em ambas as linguagens, o que pode não ser sempre o caso.
+ *Integração Difícil:* A integração entre diferentes linguagens pode ser complicada, exigindo ferramentas e técnicas adicionais para gerenciar a comunicação entre elas.

Uma solução para esse problema é a utilização de linguagens como JULIA, que são projetadas para oferecer tanto facilidade de uso quanto alto desempenho. Isso elimina a necessidade de usar duas linguagens diferentes, simplificando o desenvolvimento, a manutenção e a execução de programas complexos.

Para instalar abra o terminal do windows e rode:

```powershell
winget install julia -s msstore
julia
```

== Formulação do BEM-1d

Considere $u$ e $v$ como duas funções da variável independente x em um espaço unidimensional. A seguinte fórmula de integração por partes é bem conhecida:

$ integral_(x = x_1)^(x = x_2) u(x) d v(x) = mat(delim: "[", u(x) v(x))_(x = x_1)^(x = x_2) - integral_(x = x_1)^(x = x_2) v(x) d u(x) $

ela pode ser rescrita de maneira simplificada como:

$ integral_(x = x_1)^(x = x_2) u(x) v' (x) d x = [u(x) v(x)]_(x = x_1)^(x = x_2) - integral_(x = x_1)^(x = x_2) v(x) u' (x) d x $

onde $'$ representa a derivada. O primeiro termo do lado direito pode ser rescrito levando em conta o *contorno*. No caso unidimensional isso consiste exatamente nos pontos $x_1$ e $x_2$.

#image("../assets/apresentacao/contorno-1d.png", width: 80%)

$ [u(x) v(x)]_(x = x_1)^(x = x_2) &= [u(x) v(x) n(x)]_(x = x_2) + [u(x) v(x) n(x)]_(x = x_1) \ &=> sum_(x = x_1 , x_2) u(x) v(x) n(x) => integral_Gamma u(x) v(x) n(x) d Gamma $

Os demais termos tratam de uma integral no domínio e pode ser escritas como:

$ integral_Omega u(x) v' (x) d Omega = integral_Gamma u(x) v(x) n(x) d Gamma - integral_Omega v(x) u' (x) d Omega $

É importante notar:
1- A ideia principal da integração por partes, é trocar o operador diferencial da função $v$ para a função $u$.
2- Ao fazer essa troca, alguns termos de contorno aparecem.
3- A integração por partes foi feita apenas uma vez. No entanto, no BEM a integração por partes pode ser realizada uma, duas ou até quatro vezes dependendo da ordem da derivada.
4- Esse processo pode ser facilmente estendido para 2 ou 3 dimensões.

== Equação de Laplace

A formulação diferencial da equação de Laplace é dada por:

$ (d^2 T(x))/(d x^2) = 0 $

Usando integração por partes duas vezes na expressão dos resíduos ponderados:

$ integral_(x_0)^(x_f) T^* (d^2 T)/(d x^2) d x = (T^* (d T)/(d x))_(x_0)^(x_f) - integral_(x_0)^(x_f) (d T^*)/(d x) (d T)/(d x) d x = (T^* (d T)/(d x) - T (d T^*)/(d x))_(x_0)^(x_f) - integral_(x_0)^(x_f) (d^2 T^*)/(d x^2) T d x . $

Reescrevendo em termos do fluxo $Q = (d T)/(d x)$:

$ integral_(x_0)^(x_f) T^* (d^2 T)/(d x^2) d x = - (T Q^* - T^* Q)_(x_0)^(x_f) - integral_(x_0)^(x_f) (d^2 T^*)/(d x^2) T d x . $

Considerando que $T^*$ é a solução fundamental do problema, ou seja:

$ -(d^2 T^* (x, x_d))/(d x^2) = delta(x - x_d) . $

onde $delta$ é a função delta de Dirac.  Pode-se observar que $T^* = - 1/2 | x - x_d |$ e $Q^* = (d T^*)/(d x) = - 1/2 "sign" (x - x_d)$. Substituindo essas funções na equação integral obtém-se:

$ -T(x_d) = T(x_f) Q^* (x_f , x_d) - T^* (x_f , x_d) Q(x_f) - T(x_0) Q^* (x_0 , x_d) + T^* (x_0 , x_d) Q(x_0) . $

Quando $x_d$ é $x_0$ tem-se:

$ T^* (x_0 , x_0) = 0 \ T^* (x_f , x_0) = - (x_f - x_0)/2 \ Q^* (x_0 , x_0) = 1/2 \ Q^* (x_f , x_0) = - 1/2 $

Quando $x_d$ é $x_f$ tem-se:

$ T^* (x_0 , x_f) = - (x_f - x_0)/2 \ T^* (x_f , x_f) = 0 \ Q^* (x_0 , x_f) = 1/2 \ Q^* (x_f , x_f) = - 1/2 $

Logo, para o ponto $x_d = x_0$ a equação se torna:

$ -T(x_0) = - T(x_f) 1/2 + (x_f - x_0)/2 Q(x_f) - T(x_0) 1/2 + 0 Q(x_0) $

e para o ponto $x_d = x_f$

$ -T(x_f) = - T(x_f) 1/2 - 0 Q(x_f) - T(x_0) 1/2 - (x_f - x_0)/2 Q(0) $

Escrevendo as duas equações em forma matricial, tem-se:

$ [mat(delim: #none, 0.5, -0.5; -0.5, 0.5)] [mat(delim: #none, T(0); T(1))] = (x_f - x_0) [mat(delim: #none, 0, -0.5; 0.5, 0)] [mat(delim: #none, Q(0); Q(1))] . $

=== Exemplos

Caso 1: Condição de Dirichlet
$T(0) = 100, T(1) = 0 quad T(x) = 100 - 100 x$

Caso 2: Condição mista
$T(0) = 100, Q(0) = 0 quad T(x) = 100$

```julia
x0=0
xf=1
l=xf-x0
A1=[0.5 -.5 0 l/2
-.5 .5 -l/2 0
1 0 0 0 
0 1 0 0]
b=[0,0,100,0]
x1=A1\b

A2=[0.5 -.5 0 l/2
-.5 .5 -l/2 0
1 0 0 0 
0 0 1 0]
b=[0,0,100,0]
x2=A2\b
```

Usando essa equação podemos calcular o valor da temperatura nos pontos internos:

```julia
x0=0
xf=1
xs=range(x0,xf,length=100)
T=0.5*(x1[1]+x1[2]).+(xs.-x0)/2*x1[3].+(xs.-xf)/2*x1[4]
using Plots
plot(xs,T,legend=false,xlabel="x",ylabel="T",marker=:c)
```

== Exercícios

1- Uma outra possível condição de contorno é a convecção:

$ q = h(T - T_infinity) $

considere o problema onde o lado direito está exposto a $T_infinity = 20 ° "C"$, no lado esquerdo $T(0) = 100 ° "C"$, e $h = 3 "W/°C"$. Resolva esse problema analicamente e numericamente usando BEM.

2- Quando existe uma fonte de calor distribuída $b$ mais um termo aparece na equação integral:

$ -T(x_d) = T(x_f) Q^* (x_f , x_d) - T^* (x_f , x_d) Q(x_f) - T(x_0) Q^* (x_0 , x_d) + T^* (x_0 , x_d) Q(x_0) + integral_(x_0)^(x_f) T^* (x, x_d) b d x . $

Considere uma carga constante(isso torna a integral restante muito fácil de ser integrada analiticamente) e resolva um problema de uma grande placa de espessura L = 2 cm com condutividade térmica constante k = 1 W/m.K e geração uniforme de calor b = 1000 kW/m3. As faces A e B estão a temperaturas de 100C e 200C, respectivamente. Compare com a solução analítica:

#image("../assets/apresentacao/placa-calor.png", width: 80%)

$ T = [(T_B - T_A)/L + b/(2 k) (L - x)] x + T_A $

=== Extra

Quando dois corpos trocam calor por radiação, o fluxo é proporcional à diferença da quarta potência de suas temperaturas absolutas: $q_n = kappa f_s f_epsilon.alt (u^4 - u_R^4)$ onde u, uR são as temperaturas absolutas dos corpos radiantes, 𝜅 = 5.699 × 10−8 W∕(m2 K4) é a constante de Stefan-Boltzmann, 0 ≤ fs ≤ 1 é o fator de forma da radiação e 0 \< f𝜖 ≤ 1 é a emissividade superficial, definida como o poder emissivo relativo de um corpo em comparação ao de um corpo negro ideal. A emissividade superficial também é igual ao coeficiente de absorção, definido como a fração da energia térmica incidente em um corpo que é absorvida. A radiação pode ser vista como uma condição de contorno convectiva onde o coeficiente de transferência de calor convectivo depende da temperatura dos corpos radiantes. Escrevendo:

$ q_n = kappa f_s f_epsilon.alt (u^4 - u_R^4) = underbrace(kappa f_s f_epsilon.alt (u^2 + u_R^2)(u + u_R), h_r (u)) (u - u_R) $

Os problemas de radiação devem ser resolvidos por iteração: primeiro, o problema linear é resolvido, então o coeficiente de transferência de calor convectivo é atualizado e a solução é repetida. O critério de parada é baseado no tamanho da variação de temperatura. Normalmente, são necessárias poucas iterações.

Implemente essa condição de contorno no BEM.
