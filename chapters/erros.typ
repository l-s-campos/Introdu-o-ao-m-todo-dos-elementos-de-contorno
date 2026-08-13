// Medidas de erro
// markdown2typst + tex2typst

= Medidas de erro
<erros>

Nas comparações com solução analítica use sempre *erro relativo* normalizado. Seja $upright(bold(u))$ o vetor numérico e $upright(bold(u))^"exact"$ a referência nos mesmos pontos.

== Erro médio

$ bar.v.double upright(bold(u)) bar.v.double_"med" = sum_i | u_i | $

$ epsilon_u = (bar.v.double upright(bold(u)) - upright(bold(u))^"exact" bar.v.double_"med") / (bar.v.double upright(bold(u))^"exact" bar.v.double_"med") $

== Erro máximo

$ bar.v.double upright(bold(u)) bar.v.double_"max" = max_i | u_i | $

$ epsilon_u = (bar.v.double upright(bold(u)) - upright(bold(u))^"exact" bar.v.double_"max") / (bar.v.double upright(bold(u))^"exact" bar.v.double_"max") $

== Norma $L_2$ do erro

$ bar.v.double upright(bold(u)) bar.v.double_(L_2) = sqrt(sum_(i = 1)^(N_e) (u_i)^2) $

$ epsilon_u = (bar.v.double upright(bold(u)) - upright(bold(u))^"exact" bar.v.double_(L_2)) / (bar.v.double upright(bold(u))^"exact" bar.v.double_(L_2)) $

== No `BEM_gmsh`

Com solução analítica anexada (`attach_analytical!` ou `apply_analytical_bc!`):

```julia
err = rel_error(dad)   # erro relativo agregado no contorno
@show err
```

Reporte sempre as três medidas (médio, máximo, $L_2$) em estudos de convergência e indique se o erro é medido no contorno, em pontos internos ou em ambos.
