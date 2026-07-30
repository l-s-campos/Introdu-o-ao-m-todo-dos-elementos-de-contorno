// Poisson 2D
// markdown2typst + tex2typst

= Poisson 2D
<poisson-2d>

= MECID

A equação de Poisson é dada por:

$ nabla^2 T = f(x, y), $

um novo termo aparece na equação integral de contorno:

$ 0.5 T(x_d , y_d) = integral_Gamma T q^* d s - integral_Gamma T^* q d s + integral_Omega T^* f d Omega, $

esse termo precisa de um tratamento extra. Vamos aproximar $T^* f$ por funções de base radial:

$ T^* f = sum_(j = 1)^n phi.alt(r_j) alpha_j . $

$phi.alt$ é uma função de base radial que será tomada aqui por $r^2 log(r)$. $alpha$ são coeficientes desconhecidos que podem ser obtidos aplicando essa equação $n$ vezes e montando um sistema matricial:

$ mat(delim: "[", phi.alt(upright(X)_1 - upright(X)_1), phi.alt(upright(X)_1 - upright(X)_2), dots.c, phi.alt(upright(X)_1 - upright(X)_N); phi.alt(upright(X)_2 - upright(X)_1), phi.alt(upright(X)_2 - upright(X)_2), dots.c, phi.alt(upright(X)_2 - upright(X)_N); dots.v, dots.v, dots.down, dots.v; phi.alt(upright(X)_N - upright(X)_1), phi.alt(upright(X)_N - upright(X)_2), dots.c, phi.alt(upright(X)_N - upright(X)_N)) mat(delim: "{", alpha_1; alpha_2; dots.v; alpha_N) = mat(delim: "{", upright(f) (upright(X)_1) T^* (X_d comma X_1); upright(f) (upright(X)_2) T^* (X_d comma X_2); dots.v; upright(f) (upright(X)_N) T^* (X_d comma X_N)) = mat(delim: "{", upright(T)^* (X_d comma X_1), dots.c, 0; 0, upright(T)^* (X_d comma X_2), 0; dots.v, dots.down, dots.v; 0, 0, upright(T)^* (X_d comma X_N)) mat(delim: "{", upright(f) (upright(X)_1); upright(f) (upright(X)_2); dots.v; upright(f) (upright(X)_N)) $

escrito de maneira sintética:

$ alpha = F^(-1) D f . $

Voltando a integral de domínio:

$ integral_Omega T^* f d Omega = integral_Omega sum_(j = 1)^n phi.alt(r_j) alpha_j d Omega = sum_(j = 1)^n integral_Omega phi.alt(r_j) d Omega F^(-1) D f = s F^(-1) D f = s_f D f $

$ integral_Omega T^* f d Omega = mat(delim: "{", upright(s)_1, upright(s)_2, ..., upright(s)_n) F^(-1) mat(delim: "{", upright(T)^* (X_d comma X_1), dots.c, 0; 0, upright(T)^* (X_d comma X_2), 0; dots.v, dots.down, dots.v; 0, 0, upright(T)^* (X_d comma X_N)) mat(delim: "{", upright(f) (upright(X)_1); upright(f) (upright(X)_2); dots.v; upright(f) (upright(X)_N)) = mat(delim: "{", upright(s)_(f 1) upright(T)^* (X_d comma X_1), upright(s)_(f 2) upright(T)^* (X_d comma X_2), ..., upright(s)_(f N) upright(T)^* (X_d comma X_N)) mat(delim: "{", upright(f) (upright(X)_1); upright(f) (upright(X)_2); dots.v; upright(f) (upright(X)_N)) $

Fazendo isso para os $n$ pontos fontes:

$ mat(delim: "{", upright(s)_(f 1) upright(T)^* (X_1 comma X_1), upright(s)_(f 2) upright(T)^* (X_1 comma X_2), ..., upright(s)_(f N) upright(T)^* (X_1 comma X_N); upright(s)_(f 1) upright(T)^* (X_2 comma X_1), upright(s)_(f 2) upright(T)^* (X_2 comma X_2), ..., upright(s)_(f N) upright(T)^* (X_2 comma X_N); dots.v; upright(s)_(f 1) upright(T)^* (X_N comma X_1), upright(s)_(f 2) upright(T)^* (X_N comma X_2), ..., upright(s)_(f N) upright(T)^* (X_N comma X_N)) mat(delim: "{", upright(f) (upright(X)_1); upright(f) (upright(X)_2); dots.v; upright(f) (upright(X)_N)) = M f $

a diagonal dessa matriz não é bem definida devido a singularidade da solução fundamental. Ela será calculada de maneira indireta como feito com a matriz $H$.

Esse procedimento tem de ser capaz de calcular essa integral quando a função $f$ for constante e igual a 1.

$ I_(1 d) = integral_Omega T^* (X_d , X) d Omega $

Igualando as equações:

$ M [1] = [I_1] $

a diagonal da matriz $M$  tem de ser dada por

$ M_(i i) = I_(1 i) - sum_(j = 1)^N M_(i j) , " com " i != j, " para " i = 1, 2, . . ., N, $

```julia
Ht, Gt = BEM.calc_HeGt(dad)
A, b = BEM.aplicaCDC(Ht, Gt, dad)

M = BEM.Monta_M_RIMd(dad, npg)# calc_HeG_potencial linha 310
f=ones(nc(dad)+ni(dad))*10
x = A \ (b+M*f)#carga distribuída
x = A \ (b+M*u̇)#difusão
T, q = separa(dad, x) #format 479
Ti=x[nc(dad)+1:end]
```

== Exercícios

Todos os exercícios tem de ser resolvidos com diferentes discretizações.

+ Determine a superfície de deflexão de uma membrana elástica em forma de um triângulo equilátero com comprimento lateral a = 5,0 m. A membrana está fixa ao longo de sua borda e é submetida a uma carga distribuída uniformemente f = 10 kN/m² e uma tensão S = 1 kN/m. Os eixos coordenados são tomados como mostrado na figura.

A equação que rege esse problema é: $S nabla^2 w = − f$

#image("../assets/poisson-2d/membrana.png", width: 80%)

A solução analítica é dada por:

$ w = - f/(2 S) [1/2 (x^2 + y^2) - 1/(a sqrt(3)) (y^3 - 3 x^2 y) - 1/18 a^2] $

Calcule o erro médio e a norma L2 do erro em 100 pontos internos. Faça um gráfico com a distribuição de deflexão.

+ Considere um cubo unitário inicialmente a temperatura zero e submetido a um aquecimento súbito em uma de suas faces.

#image("../assets/poisson-2d/cubo-transiente.png", width: 80%)

A face aquecida é elevada e mantida a uma temperatura unitária. Supõe-se que as propriedades do material sejam unitárias. A solução analítica para este problema de exemplo pode ser encontrada como:

$ T(Y, t) = 1 - 4/pi sum_(n = 0)^infinity ((- 1)^n)/(2 n + 1) exp {-((2 n + 1)^2 pi^2 kappa t)/(4 L^2)} cos ((2 n + 1) pi Y)/(2 L) $

Usando algum método de análise transiente da aula 2, resolva esse problema e compare com a solução analítica.

+ Considere uma elipse que é governada por: $nabla^2 u = 4 - x^2$

a solução analítica é dada por:

$u = [1.6 - 1/246 (50 x^2 - 8 y^2 + 33.6)]((x^2)/4 + y^2 - 1)$e $mat(delim: #none, q = 0.4(x^2 + 8 y^2) + 1/246 (-50 x^3 - 96 x y^2 + 83.2 x) x/2 +; 1/246 (-96 x^2 y + 32 y^3 - 83.2 y) y)$calcule os erros do potencial nos pontos internos e do fluxo no contorno.

#image("../assets/poisson-2d/elipse.png", width: 80%)

extra

Um cilindro oco com temperatura inicial zero é considerado. O raio interno a=1 e o raio externo b = 2. A superfície interna do cilindro é mantida a uma temperatura=1. Nesse caso, são impostas restrições de simetria. As propriedades do material são unitárias. A solução analítica é dada por:

$ T(r, t) = (upright(l n) (b \/ r))/(upright(l n) (b \/ a)) T_1 + pi sum_(n = 1)^infinity (J_0 (b alpha_n) J_0 (a alpha_n))/(J_0^2 (a alpha_n) - J_0^2 (b alpha_n)) {J_0 (r alpha_n) Y_0 (b alpha_n) - J_0 (b alpha_n) Y_0 (r alpha_n)} upright(e)^(-i alpha_n^2) , $

onde $T_i$ é a temperatura constante na superfície interna, $J_0$ e $Y_0$ são as funções de Bessel de primeira e segunda espécies, respectivamente, e $alpha$ é a raiz de:

$ J_0 (a x) Y_0 (b x) - J_0 (b x) Y_0 (a x) = 0 . $

você pode usar a função #link("https://juliamath.github.io/Roots.jl/stable/roots/#Searching-for-all-zeros-in-an-interval")[fzeros] para encontrar $alpha$

== Gravação

#link("https://youtu.be/uSREar_ejnM")[https://youtu.be/uSREar\_ejnM]
