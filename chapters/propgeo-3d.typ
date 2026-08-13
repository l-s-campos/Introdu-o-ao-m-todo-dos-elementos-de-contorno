// Propgeo 3D
// markdown2typst + tex2typst + BEM_gmsh

= Propgeo 3D
<propgeo-3d>

A mesma ideia do capítulo *Indo para 2D* se estende a sólidos: volume, área superficial e centróide podem ser obtidos por integrais *apenas na superfície* fechada $Gamma = partial Omega$.

No `BEM_gmsh`:

```julia
using DrWatson
@quickactivate :BEM

# Requer malha de superfície fechada (ex.: cubo em data/Laplace)
# dad = format3d(datadir("Laplace", "cubo.msh"), Laplace(1.0))
# gp  = geometric_props(dad)   # GeometricProps3D
# @show gp.surface_area gp.volume gp.centroid
```

O despacho é automático: `dad.dimension == 3` seleciona `geometric_props_3d`. A integração usa as funções de forma 2D dos elementos de superfície e as normais já armazenadas no `BEMdata`.

== Exercício

1. Gere um cubo (ou esfera faceted) no Gmsh, leia com `format3d` e compare área, volume e centróide com o valor analítico.
2. Introduza um furo cilíndrico (superfície interna orientada para dentro do material) e repita o cálculo.

== Material complementar

- #link("https://1drv.ms/f/s!AmfyGvdmTYongqYn5kjjlZaMHr9h2w?e=z0sXvU")[Arquivos legados propgeo 3D]
- #link("https://youtu.be/Uc-rxXDBU6I")[gravação]
- Código atual: `src/Core/GeometricProperties.jl` no `BEM_gmsh`
