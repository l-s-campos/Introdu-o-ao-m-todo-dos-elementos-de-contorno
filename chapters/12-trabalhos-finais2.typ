// Trabalhos finais — RASCUNHO (polimento)
// Oficial ainda é 12-trabalhos-finais.typ até aprovação

= Trabalhos finais
<trabalhos-finais>

Curso de *30 horas* — Introdução ao Método dos Elementos de Contorno. \
Código: #link("https://github.com/l-s-campos/BEM_gmsh")[`BEM_gmsh`].

Propostas de *monografia curta*: formulação, malha/CDC, estudo quantitativo e discussão crítica.
O núcleo BEM já está no repositório — a equipe usa a biblioteca como *plataforma*, não como caixa-preta.

== Regras gerais

*Equipes:* 1–2 alunos. *Uma* proposta por equipe.

*Entregáveis (todas):*

1. *Relatório* (12–20 p., PDF): formulação BEM do problema; 3–8 refs; malha/CDC; herdado vs. próprio; números/tabelas/figuras; erro/convergência/custo (*Apêndice: medidas de erro*); limites e futuro.
2. *Pasta reproduzível:* scripts, `.geo`, `README` com comandos, seeds, versões.
3. *Contribuição própria* (≥ 2): malha nova; CDC/pós não trivial; comparação de dois métodos; estudo paramétrico físico; implementação auxiliar.
4. *Defesa* 10–15 min (formulação + código).

*Não conta:* copiar `intro.jl`, mudar `ndiv`, plotar $T=x$.

*Rubrica:* formulação 20 % · originalidade 25 % · rigor numérico 25 % · discussão 15 % · reprodutibilidade 15 %.

== Propostas (resumo + âncoras)

=== A — Fratura (Dual BEM, $K_I$, propagação)
`src/Crack/`, `center_crack.jl`, `scripts/crack_central.jl`, Feddersen. Curva $K_I(a\/W)$; convergência na ponta; ≥ 1 passo MTS.

=== B — Multirregião / contraste
`src/MultiRegion/`, `two_regions.jl`. Salto na interface; erro vs. contraste $k_1\/k_2$; geometria própria.

=== C — Dinâmica (onda \/ calor)
`DIBEM` + `solve_transient` \/ Houbolt \/ `solve_transient_o2`; `wave_propagation.jl`. CFL prático; dispersão; RBF vs. tempo.

=== D — H-matrizes (Pareto erro × custo)
`H_G_Hmat`, `src/Hmat/`, scripts de profile. Denso vs. hierárquico; ver *extra* do cap. Laplace 2D.

=== E — Contato (Hertz / fretting)
`src/Contact/`, `hertz_line_2d.jl`, Pohrt–Li. $p(x)$ vs. analítico; zona ativa.

== Cronograma sugerido

W1 plano · W2 primeiro número · W3 estudo a meio · W4 relatório + defesa.


== Escolha orientada

#table(
  columns: (auto, auto, auto),
  inset: 6pt,
  stroke: 0.5pt + luma(200),
  [*Se curte…*], [*Prefira*], [*Cuidado*],
  [Sólidos \/ $K_I$], [A], [Singularidade; malha dual],
  [Interfaces], [B], [Condicionamento],
  [Dinâmica], [C], [erro DIBEM],
  [HPC], [D], [Medir tempo/memória],
  [Contato], [E], [Zona ativa],
)

O detalhamento completo de frentes e entregas específicas de cada proposta permanece no capítulo oficial
`12-trabalhos-finais.typ` até a fusão deste rascunho (este arquivo é o *mapa* + regras limpas;
ao aprovar, fundir o detalhe A–E do oficial com este cabeçalho).
