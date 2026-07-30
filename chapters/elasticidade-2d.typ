// Elasticidade 2D
// markdown2typst + tex2typst

= Elasticidade 2D
<elasticidade-2d>

= Equações importantes

A relação de equilíbrio de tensão pode ser escrita como:

$ (partial sigma_(x x))/(partial x) + (partial sigma_(x y))/(partial y) + f_x = 0 \ (partial sigma_(y x))/(partial x) + (partial sigma_(y y))/(partial y) + f_y = 0 $

$ (partial sigma_(i j))/(partial x_j) + f_i = 0 $

Relação deformação-deslocamento pode ser escrita como:

$ epsilon_(x x) = (partial u_x)/(partial x) ; quad epsilon_(y y) = (partial u_y)/(partial y) ; quad epsilon_(x y) = 1/2 ((partial u_x)/(partial y) + (partial u_y)/(partial x)) $

$ epsilon_(i j) = 1/2 ((partial u_i)/(partial x_j) + (partial u_j)/(partial x_i)) $

Relação tensão-deformação pode ser escrita como:

$ & upright(bold(epsilon))_(x x) = (1/E) sigma_(x x) + ((-v)/E) sigma_(y y) \ & epsilon_(y y) = ((-v)/E) sigma_(x x) + (1/E) sigma_(y y) \ & epsilon_(x y) = 1/(2 mu) sigma_(x y) $

$ mu = E/(2(1 + v)) $

Essas três relações podem ser usadas para chegar na seguinte equação diferencial de deslocamentos:

$ (partial^2 u_x)/(partial x^2) + (partial^2 u_x)/(partial y^2) + 1/(1 - 2 v) ((partial^2 u_x)/(partial x^2) + (partial^2 u_y)/(partial x partial y)) = (-f_x)/mu \ (partial^2 u_y)/(partial x^2) + (partial^2 u_y)/(partial y^2) + 1/(1 - 2 v) ((partial^2 u_y)/(partial y^2) + (partial^2 u_x)/(partial x partial y)) = (-f_y)/mu $

$ (partial^2 u_i)/(partial x_j partial x_j) + (1/(1 - 2 v)) (partial^2 u_j)/(partial x_i partial x_j) = (-f_i)/mu $

A solução fundamental dessa equação de deslocamento é dada por:

$ U_(i j) (p, Q) = 1/(8 pi mu(1 - v)) {(3 - 4 v) ln [1/(r(p, Q))] delta_(i j) + (partial r(p, Q))/(partial x_i) (partial r(p, Q))/(partial x_j)} $

e para a tração:

$ T_(i j) (p, Q) &= (-1)/(4 pi(1 - nu) r(p, Q)) ((partial r(p, Q))/(partial n)) \
& times [(1 - 2 nu) delta_(i j) + 2 (partial r(p, Q))/(partial x_i) (partial r(p, Q))/(partial x_j)] \
& +(1 - 2 v)/(4 pi(1 - v) r(p, Q)) [(partial r(p, Q))/(partial x_j) n_i - (partial r(p, Q))/(partial x_i) n_j] $

Usando um procedimento semelhante ao apresentado para o problema potencial obtemos a equação integral de contorno:

$ mat(delim: "[", u_x (p); u_y (p)) & +integral_Gamma mat(delim: "[", T_(x x) (p comma Q), T_(x y) (p comma Q); T_(y x) (p comma Q), T_(y y) (p comma Q)) mat(delim: "[", u_x (Q); u_y (Q)) dif Gamma(Q) \ &= integral_Gamma mat(delim: "[", U_(x x) (p comma Q), U_(x y) (p comma Q); U_(y x) (p comma Q), U_(y y) (p comma Q)) mat(delim: "[", t_x (Q); t_y (Q)) dif Gamma(Q) $

$ u_i (p) + integral_Gamma T_(i j) (p, Q) u_j (Q) "d" Gamma(Q) = integral_Gamma U_(i j) (p, Q) t_j (Q) "d" Gamma(Q) $

Essa equação pode ser derivada e junto com as relações tensão deformação obtém-se a equação que pode ser usada para calcular a tensão:

$ sigma_(i j) (p) + integral_Gamma lr(\{(2 mu nu)/(1 - 2 nu) delta_(i j) (partial T_(m k) (p, Q))/(partial x_m)) \ +mu [(partial T_(i k) (p, Q))/(partial x_j) + (partial T_(j k) (p, Q))/(partial x_i)]} u_k (Q) dif Gamma(Q) \ = integral_Gamma lr(\{(2 mu nu)/(1 - 2 nu) delta_(i j) (partial U_(m k) (p, Q))/(partial x_m)) + mu [(partial U_(i k) (p, Q))/(partial x_j) + (partial U_(j k) (p, Q))/(partial x_i)]} t_k (Q) dif Gamma(Q) $

que pode ser reescrita como:

$ sigma_(i j) (p) + integral_Gamma S_(k i j) (p, Q) u_k (Q) upright(space.nobreak d) Gamma(Q) = integral_Gamma D_(k i j) (p, Q) t_k (Q) upright(space.nobreak d) Gamma(Q) $

onde os tensores $S_(k i j)$ e $D_(k i j)$ são dados por:

$ S_(k i j) (p, Q) &= mu/(2 pi(1 - nu)) (1/(r^2)) n_i [2 nu (partial r)/(partial x_j) (partial r)/(partial x_k) + (1 - 2 nu) delta_(j k)] \
& +mu/(2 pi(1 - nu)) (1/(r^2)) n_j [2 nu (partial r)/(partial x_i) (partial r)/(partial x_k) + (1 - 2 nu) delta_(i k)] \
& +mu/(2 pi(1 - v)) (1/(r^2)) n_k [2(1 - 2 v) (partial r)/(partial x_i) (partial r)/(partial x_j) -(1 - 4 v) delta_(i j)] \
& +mu/(pi(1 - nu)) (1/(r^2))((partial r)/(partial n)) [(1 - 2 nu) delta_(i j) (partial r)/(partial x_k) + nu(delta_(j k) (partial r)/(partial x_i) + delta_(i k) (partial r)/(partial x_j)) \
& -4 (partial r)/(partial x_i) (partial r)/(partial x_j) (partial r)/(partial x_k)] $

$ D_(k i j) (p, Q) = 1/(4 pi(1 - nu)) (1/r) [(1 - 2 nu)(delta_(j k) (partial r)/(partial x_i) + delta_(i k) (partial r)/(partial x_j) - delta_(i j) (partial r)/(partial x_k)) \ + 2 (partial r)/(partial x_i) (partial r)/(partial x_j) (partial r)/(partial x_k)] $

== Exercícios

=== Cilindro pressurizado

Raio interno $R_a = 50 m m$

Raio externo   $R_b = 100 m m$

Modulo de elasticidade   $E = 200 G P a$

Poisson   $nu = 0.32$

Pressão $P = 100 N \/ m m$

#image("../assets/elasticidade-2d/cilindro.png", width: 80%)

$ u_r = ((1 + v) p r_a^2)/((r_b^2 - r_a^2) E) [(1 - 2 v) r + (r_b^2)/r] \
sigma_r = (p r_a^2)/(r_b^2 - r_a^2) - (r_a^2 r_b^2 p)/(r_b^2 - r_a^2) 1/(r^2) \
sigma_h = (p r_a^2)/(r_b^2 - r_a^2) + (r_a^2 r_b^2 p)/(r_b^2 - r_a^2) 1/(r^2) $

=== Placa com furo

Raio $R = 50 m m$

Modulo de elasticidade   $E = 100 G P a$

Poisson   $nu = 0.25$

Pressão $P = 1 N \/ m m$

Esse problema é de estado plano de tensão. Para ser tratado pelas soluções fundamentais de estado plano de deformação as propriedades tem de ser ajustadas:

#image("../assets/elasticidade-2d/placa-furo.png", width: 80%)

$ nu' = nu/(1 + nu) \ E' = E [1 - (nu^(' 2))/((1 + nu')^2)] $

A solução analítica é dada por:

$ mat(delim: #none, sigma_r = sigma/2 (1 - (a^2)/(r^2)) + sigma/2 (1 + (3 a^4)/(r^2) - (4 a^2)/(r^2)) cos(2 theta); sigma_theta = sigma/2 (1 + (a^2)/(r^2)) - sigma/2 (1 + (3 a^4)/(r^4)) cos(2 theta); tau_(r theta) = - sigma/2 (1 - (3 a^4)/(r^4) + (4 a^2)/(r^2)) sin(2 theta)) $

=== Viga

Uma viga com um carregamento parabólico $t_2 (y) = - P/(2 I) ((D^2)/4 - y^2)$ tem solução analítica:

Comprimento $L = 48 m m$

Altura $D = 12 m m$

Modulo de elasticidade   $E = 300 G P a$

Poisson   $nu = 0.3$

Pressão $P = 1000 N$

#image("../assets/elasticidade-2d/viga.png", width: 80%)

$ u_1 (x, y) = - (P y)/(6 E I) [(6 L - 3 x) x + (2 + v)(y^2 - (D^2)/4)] $

$ u_2 (x, y) = P/(6 E I) [(3 v) y^2 (L - x) + (4 + 5 v) (D^2 x)/4 + (3 L - x) x^2] $

$ sigma_(x x) (x, y) = - (P(L - x) y)/I \ tau_(x y) (x, y) = - P/(2 I) ((D^2)/4 - y^2) $

#link("https://youtu.be/R6-_ECEQXRk")[gravação]
