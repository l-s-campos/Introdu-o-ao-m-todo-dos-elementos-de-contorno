// Propgeo 3D — RASCUNHO

= Propgeo 3D
<propgeo-3d>

Extensão do capítulo *Indo para 2D*: volume, área superficial e centróide de um sólido
a partir *só* da superfície fechada $Gamma = partial Omega$.

== Objetivos

+ Reconhecer a redução volume $arrow.r$ superfície (divergência \/ radial 3D).
+ Rodar `format3d` + `geometric_props` num sólido fechado.
+ Conferir orientação (normal para fora) e o efeito de cavidades.

== Ideia

Em 2D, área e momentos saíram de integrais em $Gamma$ (radial ou Green).
Em 3D, com normal exterior $upright(bold(n))$:

$
"Área" = integral_Gamma dif Gamma ,
quad
"Vol" = (1\/3) integral_Gamma upright(bold(x)) · upright(bold(n)) dif Gamma
$

(e análogos para $integral x_i dif V$ via primitivas radiais — como em `geometric_props_3d`).

O pacote integra com as funções de forma *superficiais* já no `dad` (tri\/quad) e normais armazenadas.

== No `BEM_gmsh`

```julia
using DrWatson
@quickactivate :BEM

# malha de superfície fechada (ex. cubo em data/Laplace)
dad = format3d(datadir("Laplace", "cubo.msh"), Laplace(1.0))
gp  = geometric_props(dad)   # despacha geometric_props_3d se dimension==3
@show gp.surface_area gp.volume gp.centroid
```

Essência do laço (`GeometricProperties.jl`):

```julia
# para cada nó de contorno / ponto de Gauss superficial:
#   Area += wJ
#   jfac = dot(n, x) / R^3
#   (Fv,Fx,Fy,Fz) = _calc_F_3d(...)   # primitivas radiais
#   Vol += Fv * jfac * wJ   etc.
```

Parâmetros: `npg_radial`, `npg_boundary` (reavalia quadratura se não quiser os pesos do `dad`).

== Boas práticas

1. Superfície *fechada* e manifold orientável.
2. Normal *saindo* do material; cavidades: orientação oposta (como furos 2D).
3. Compare com valor analítico antes de confiar em malha grossa.
4. Elementos degenerados \/ normal nula: lixo em $V$ — valide `plot` da malha.

== Exercícios

1. Cubo (ou esfera facetada): `format3d` + `geometric_props` vs $A$, $V$, centróide analíticos; refine `ndiv`.
2. Cubo com furo cilíndrico: superfície interna orientada para o material; repita $V$ e $A$.
3. (Opcional) Compare `npg_boundary=nothing` (pesos do `dad`) vs `npg_boundary=16`.

== Leituras

- Cap. *Indo para 2D* (radial \/ divergência)
- `src/Core/GeometricProperties.jl`
- #link("https://1drv.ms/f/s!AmfyGvdmTYongqYn5kjjlZaMHr9h2w?e=z0sXvU")[Arquivos legados]
- #link("https://youtu.be/Uc-rxXDBU6I")[gravação]
