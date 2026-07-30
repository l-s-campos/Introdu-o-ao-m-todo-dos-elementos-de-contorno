// Viga de Euler
// markdown2typst + tex2typst

= Viga de Euler
<viga-euler>

== Teoria de vigas

Considerando a viga representada a equação governante do problema pode ser escrita como :

$E I (dif^4 u)/(dif x^4) = q(x)$

A partir de $u$  as podemos definir:

$phi = (dif u)/(dif x) , quad M = - E I (dif^2 u)/(dif x^2) , quad upright(e) quad Q = - E I (dif^3 u)/(dif x^3)$ que são a rotação, momento fletor e cortante.

#image("../assets/viga-euler/viga-teoria.png", width: 80%)

usando essas definições, aplicando o método dos resíduos ponderados e integração por partes 4 vezes obtemos:

$ u(xi) &= + u^* (xi, x) Q(x) |_(x = L) - u^* (xi, x) Q(x) |_(x = 0) \ & -lr(phi^* (xi, x) M(x)|)_(x = L) + lr(phi^* (xi, x) M(x)|)_(x = 0) \ & +lr(M^* (xi, x) phi(x)|)_(x = L) - lr(M^* (xi, x) phi(x)|)_(x = 0) \ & -Q^* (xi, x) u(x) |_(x = L) + Q^* (xi, x) u(x) |_(x = 0) \ & +integral_0^L u^* (xi, x) q(x) dif x, $

onde:

$ u^* (xi, x) = (r^3)/(12 E I) , \

phi^* (xi, x) = (r^2)/(4 E I) ((partial r)/(partial x)), \
M^* (xi, x) = - r/2 ((partial r)/(partial x))^2 , \
Q^* (xi, x) = - 1/2 ((partial r)/(partial x))^3 " e " \
r = | x - xi | $

onde $(partial r)/(partial x) = - 1 upright(quad s e quad) x < xi$ e $(partial r)/(partial x) = 1 upright(quad s e quad) x > xi$.

Para esse problema precisamos de definir mais uma equação para conseguirmos montar um sistema de equação adequado. Isso é obtido derivando a equação anterior em relação à $xi$:

$ phi(xi) = + lr((partial u^* (xi, x))/(partial xi) Q(x)|)_(x = L) - lr((partial u^* (xi, x))/(partial xi) Q(x)|)_(x = 0) \ -lr((partial phi^* (xi, x))/(partial xi) M(x)|)_(x = L) + lr((partial phi^* (xi, x))/(partial xi) M(x)|)_(x = 0) \ +lr((partial M^* (xi, x))/(partial xi) phi(x)|)_(x = L) - lr((partial M^* (xi, x))/(partial xi) phi(x)|)_(x = 0) \ +integral_0^L (partial u^* (xi, x))/(partial xi) q(x) dif x, $

onde:

$ & (partial u^* (xi, x))/(partial xi) = (r^2)/(4 E I) ((partial r)/(partial xi)), \ & (partial phi^* (xi, x))/(partial xi) = r/(2 E I) ((partial r)/(partial xi))((partial r)/(partial x)), \ & (partial M^* (xi, x))/(partial xi) = 1/2 ((partial r)/(partial xi)) ((partial r)/(partial x))^2 , $

onde $(partial r)/(partial xi) = 1 upright(quad s e quad) x < xi$ e $(partial r)/(partial xi) = - 1 upright(quad s e quad) x > xi$.

```julia
function solfundviga(ξ, x,E,I,L)
	r=abs(x-ξ)
	if x==ξ
		if ξ==0
			drdx=1
		elseif ξ==L
			drdx=-1
		else
			drdx=0
		end
	else
		drdx=(x-ξ)/r
	end
	drdξ=-drdx
	
	u_star = r^3 / (12 * E * I)
	phi_star = (r^2 / (4 * E * I)) * drdx
	M_star = -r/2 * drdx^2
	Q_star = -1/2 * drdx^3
	
	du_star_dxi = (r^2 / (4 * E * I)) * drdξ
	dphi_star_dxi = (r / (2 * E * I)) * drdξ * drdx
	dM_star_dxi = 1/2 * drdξ * drdx^2
	
	u_star ,phi_star ,M_star ,Q_star ,du_star_dxi ,dphi_star_dxi ,dM_star_dxi
end
```

Considere uma viga com $L = 4 "m"$, $E = 50 "GPa"$, $I = 0.0036 "m"^2$ e uma carga pontual no centro da viga de valor $10$KN.

```julia
L=4
E=50e9
Iv=0.0036

#solfundviga(ξ, x,E,I,L)
u_00 ,phi_00 ,M_00 ,Q_00 ,du_00 ,dphi_00,dM_00= solfundviga(0,0,E,Iv,L)
u_0L ,phi_0L ,M_0L ,Q_0L ,du_0L ,dphi_0L,dM_0L= solfundviga(0,L,E,Iv,L)
u_L0 ,phi_L0 ,M_L0 ,Q_L0 ,du_L0 ,dphi_L0,dM_L0= solfundviga(L,0,E,Iv,L)
u_LL ,phi_LL ,M_LL ,Q_LL ,du_LL ,dphi_LL,dM_LL= solfundviga(L,L,E,Iv,L)
#multiplica U e phi
H=-[Q_00 -M_00 -Q_0L  M_0L
		0 dM_00 0  -dM_0L
		Q_L0 -M_L0 -Q_LL  M_LL
		0 dM_L0 0  -dM_LL]
		#multiplica V e M
G=[-u_00 phi_00 u_0L  -phi_0L
		-du_00 dphi_00 du_0L  -dphi_0L
		-u_L0 phi_L0 u_LL  -phi_LL
		-du_L0 dphi_L0 du_LL  -dphi_LL]

uc0 ,phic0 ,Mc0 ,Qc0 ,duc0 ,dphic0,dMc0=solfundviga(0,L/2,E,Iv,L)
ucL ,phicL ,McL ,QcL ,ducL ,dphicL,dMcL=solfundviga(L,L/2,E,Iv,L)
b=1e4*[uc0;duc0 ;ucL;ducL]

A=[-G[:,1] H[:,2] -G[:,3] H[:,4]]
x=A\b
```

=== Exercício

1- Generalize esse código para qualquer condição de contorno. Ele é capaz de resolver problemas hiperestáticos? Como?

2- Generalize esse código para qualquer carga distribuída e compare com a solução analítica.

= Efeitos transientes

#image("../assets/viga-euler/viga-transiente.png", width: 80%)

$ E I (partial^4 u)/(partial x^4) + rho A (partial^2 u)/(partial t^2) = f(x, t) $

Dividindo a equação por $rho A$ e definindo

$ c = sqrt((E I)/(rho A)) $

obtém-se a expressão:

$ c^2 (partial^4 u)/(partial x^4) + (partial^2 u)/(partial t^2) = (f(x, t))/(rho A) $

a solução fundamental para esse problema é dada por:

$ u^* (x, xi, t, tau) = 1/c {r/2 [upright(S) (r/sqrt(2 pi a)) - upright(C) (r/sqrt(2 pi a))] + sqrt(a)/sqrt(2 pi) [sin((r^2)/(4 a)) + cos((r^2)/(4 a))]} $

onde as funções S e C são denominadas #link("https://kiranshila.github.io/FresnelIntegrals.jl/dev/")[integrais de Fresnel], $r = | x - xi |$e $a = c(t - tau)$.

$ theta^* (x, xi, t, tau) = + 1/c {1/2 [upright(S) (r/sqrt(2 pi a)) - upright(C) (r/sqrt(2 pi a))]}((partial r)/(partial x)), \ M^* (x, xi, t, tau) = - (E I)/c {1/(2 sqrt(2 pi a)) [sin((r^2)/(4 a)) - cos((r^2)/(4 a))]}((partial r)/(partial x))^2 , \ Q^* (x, xi, t, tau) = - (E I)/c {r/(4 sqrt(2 pi a^3)) [sin((r^2)/(4 a)) + cos((r^2)/(4 a))]}((partial r)/(partial x))^3 , $

A equação integral para esse problema é dada por:

$ u(xi, t) = & -1/(rho A) integral_0^t [+ u^* Q - theta^* M + M^* theta - Q^* u]_(x = 0) dif tau \ & +1/(rho A) integral_0^t [+ u^* Q - theta^* M + M^* theta - Q^* u]_(x = L) dif tau \ & +1/(rho A) integral_0^t integral_0^L u^* f dif x dif tau . $

A solução do problema, que envolve quatro incógnitas de contorno, requer ao menos
duas equações integrais distintas. Assim, escreve-se, para a rotação:

$ theta(xi, t) = & -1/(rho A) {integral_0^t [(partial u^*)/(partial xi) Q - (partial theta^*)/(partial xi) M + (partial M^*)/(partial xi) theta - (partial Q^*)/(partial xi) u]_(x = 0) dif tau} \ & +1/(rho A) {integral_0^t [(partial u^*)/(partial xi) Q - (partial theta^*)/(partial xi) M + (partial M^*)/(partial xi) theta - (partial Q^*)/(partial xi) u]_(x = L) dif tau} \ & +1/(rho A) {integral_0^t integral_0^L (partial u^*)/(partial xi) f dif x dif tau}, $

=== Exercício

3 - Faça um gráfico 3d de $u^*$ e $d u^* \/ d xi$. Considere $xi$ e $tau$ iguais a zero, que o tempo varia de 0 a 10s e que o  $x$ varia de $0$ a $L$. #link("https://docs.juliaplots.org/stable/gallery/gr/generated/gr-ref050/#gr_ref050")[Exemplo de superfície].

= Desafio

1- Implemente um código que usando a formulação transiente e resolva um problema de vibração livre.

2-Esse problema também pode ser resolvido usando a solução fundamental estática e uma integral de domínio como feito para o problema de Poisson. Faça isso usando Houbolt e compare.

- Compare com uma solução analítica disponível no #link("https://www.pearson.com/en-us/subject-catalog/p/mechanical-vibrations/P200000003425/9780137515288")[Rao]
  #image("../assets/viga-euler/rao-tabela.png", width: 80%)

3-Implemente a formulação da #link("https://link.springer.com/article/10.1007/s40996-020-00359-z")[viga de Timoshenko] e de #link("https://link.springer.com/article/10.1007/s00366-019-00774-5")[BICKFORD-REDDY] e compare as três para diferentes tamanhos de viga.
