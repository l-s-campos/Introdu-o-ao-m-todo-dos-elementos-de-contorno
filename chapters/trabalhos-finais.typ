// Trabalhos finais — propostas avançadas (curso 30 h)
// Ancoradas nas capacidades reais do BEM_gmsh (não são tutoriais)

= Trabalhos finais
<trabalhos-finais>

Curso de *30 horas* — Introdução ao Método dos Elementos de Contorno. \
Código: #link("https://github.com/l-s-campos/BEM_gmsh")[`BEM_gmsh`].

Estas propostas *não* são exercícios de “rodar o exemplo e plotar”. Cada uma exige
*leitura de formulação*, *extensão de malha/CDC*, *estudo quantitativo* e
*discussão crítica* no nível de uma monografia curta. O núcleo do BEM já está no
repositório — o trabalho é *usar a biblioteca como plataforma de pesquisa*, não
como caixa-preta.

== Regras gerais

*Equipes:* 1–2 alunos. *Uma* proposta por equipe (combinar com o docente se quiser
fusão parcial).

*Carga esperada:* 20–40 h extraclasse depois do núcleo Laplace/elasticidade.
Começar cedo.

*Entregáveis obrigatórios (todas as propostas):*

1. *Relatório* (12–20 páginas, PDF) com:
   - formulação matemática do problema *no formalismo BEM* (não só a PDE);
   - estado da arte curto (3–8 referências relevantes, não genéricas);
   - descrição da malha Gmsh e das CDCs (tipos físicos);
   - o que foi herdado do `BEM_gmsh` vs. o que a equipe implementou/adaptou;
   - resultados com *números*, tabelas e figuras;
   - análise de erro / convergência / custo;
   - limitações e trabalho futuro.
2. *Repositório ou pasta reproduzível:* scripts Julia, `.geo`/geradores, `README`
   com comandos exatos (`julia --project=... script.jl`), seeds e versões.
3. *Contribuição própria* (pelo menos *dois* itens da lista):
   - gerador de malha *novo* (geometria que não está pronta em `data/`);
   - CDC ou pós-processamento não trivial (export Gmsh, sensor em linha, SIF, …);
   - comparação sistemática entre *dois* métodos/solvers/discretizações;
   - estudo de parâmetro físico com curva e interpretação;
   - implementação auxiliar (RBF alternativa, critério de parada, remalhagem simples, …).
4. *Defesa oral* de 10–15 min (perguntas sobre formulação e sobre o código).

*O que *não* conta como trabalho final:*
copiar `scripts/intro.jl`, mudar `ndiv` e entregar um gráfico de $T=x$.

*Rubrica sugerida:*

#table(
  columns: (auto, auto),
  inset: 6pt,
  stroke: 0.5pt + luma(200),
  [*Critério*], [*Peso*],
  [Formulação e domínio do BEM no problema escolhido], [20 %],
  [Originalidade / contribuição além dos scripts prontos], [25 %],
  [Rigor numérico (convergência, verificação, analítico/literatura)], [25 %],
  [Discussão crítica (limites, custo, o que falhou)], [15 %],
  [Reprodutibilidade e clareza do relatório/código], [15 %],
)

#pagebreak()

== Proposta A — Mecânica da fratura com Dual BEM (trinca + $K_I$ + propagação)

*Nível.* Avançado. \
*Âncoras no código.* `src/Crack/`, `data/elastico/iso/center_crack.jl`,
`scripts/crack_central.jl`, `data/examples/crack_feddersen.jl`, grupos físicos
tipo `"5;…"` (faces duais).

=== Problema científico

Placa finita com trinca central (ou trinca de borda / geometria *nova* proposta
pela equipe) sob tração remota. Calcular fatores de intensidade de tensão
($K_I$, e $K_(I I)$ se houver modo misto) pelo *Dual BEM* (BIE de deslocamento em
uma face + BIE de tração na face coincidente) e confrontar com solução de
referência (Feddersen / Tada / handbook).

Em seguida, dar *pelo menos um passo de propagação* com critério MTS (e, se
couber, estimativa de vida via Paris já esboçada no script) — discutindo o que o
código faz e o que ainda é simplificado.

=== O que deve ir além do script pronto

Escolha *no mínimo duas* frentes:

- Variar $a/W$ em uma curva $K_I(a/W)/K_I^infinity$ e comparar com a correção de
  Feddersen (ou outra da literatura) em *vários* pontos — não um único $a$.
- Estudo de convergência em $n_"crack"$, ordem do elemento e `npg`; extrair taxa
  observada e discutir singularidade $r^(-1/2)$ na ponta.
- Geometria alternativa: trinca de borda, trinca inclinada (modo misto), ou
  dois furos + trinca — malha Gmsh *da equipe*.
- Pós-processamento: perfil de COD (abertura) ao longo da trinca; estimativa de
  $K_I$ por extrapolação de COD vs. por tensão; comparar.
- Um experimento numérico de *caminho*: vários passos MTS com remalhagem manual
  ou semi-automática (mesmo que tosca), documentando falhas.

=== Entregas específicas

- Tabela $K_I^"num"$ vs. $K_I^"ref"$ com erro relativo para ≥ 4 razões $a/W$ ou
  ≥ 4 malhas.
- Figura da malha com faces duais destacadas e normais opostas.
- Discussão: por que o BEM dual evita a degeneração de faces coincidentes; o que
  acontece se a face B for mal orientada.
- Seção “o que o Dual BEM ainda não faz neste trabalho” (contato entre faces,
  3D, plástico, …).

=== Pontos de partida (não são o trabalho)

```julia
using DrWatson
@quickactivate :BEM
using .Crack
include(datadir("elastico", "iso", "center_crack.jl"))
# ver scripts/crack_central.jl e a API CrackTip / propagate!
```

#pagebreak()

== Proposta B — Multirregião, interfaces e materiais contrastantes

*Nível.* Avançado. \
*Âncoras no código.* `src/MultiRegion/`, `data/Laplace/two_regions.jl`,
`scripts/two_regions_interface.jl`; elasticidade anisotrópica em
`data/elastico/aniso/` (Lekhnitskii); DIBEM se houver fonte por região.

=== Problema científico

Sistemas em que *um único* contorno exterior não basta: duas ou mais sub-regiões
com condutividades (ou módulos) diferentes, acopladas por condições de
interface (continuidade de $T$ e equilíbrio de $q$, ou analogamente em
elasticidade).

A equipe formula o bloco do sistema global, implementa/adapta um caso com
*contraste forte* de propriedades e quantifica:

- salto numérico na interface ($|T_a - T_b|$, $|q_a + q_b|$);
- erro vs. solução analítica *por região* (quando existir) ou vs. solução de
  referência monodomínio equivalente;
- sensibilidade ao contraste $k_1/k_2$ (ou $E_1/E_2$) em vários logaritmos de
  razão.

=== O que deve ir além do script pronto

Pelo menos duas frentes:

- Geometria *não* retângulo-retângulo: inclusão circular/elíptica, parede
  compostas, ou três regiões; malha Gmsh própria com grupos por região.
- Caso com solução analítica clássica (ex.: cilindro composto sob pressão /
  condução radial em coroas concêntricas com $k$ diferentes) *ou* verificação
  por refino cruzado + conservação de fluxo global.
- Comparar montagem multirregião vs. “truque” de um só domínio quando o contraste
  → 1 (regressão).
- Extensão: interface com resistência de contato térmico ($q = h(T_a-T_b)$) —
  mesmo que implementada de forma simplificada no sistema de interface — *ou*
  um caso elástico anisotrópico (`AnisotropicElasticity` / Lekhnitskii) com
  verificação de patch ou furo.
- Perfil de $q$ normal ao longo da interface e balanço integral
  $integral_Gamma_("ext") q dif s approx 0$ (estacionário sem fonte).

=== Entregas específicas

- Diagrama de blocos do sistema algébrico multirregião (incógnitas por região +
  multiplicadores/condições de interface).
- Curva erro e/ou fluxo residual vs. contraste e vs. $N$.
- Discussão de condicionamento: o que acontece com $k_1/k_2 = 10^3$ ou $10^(-3)$.

=== Pontos de partida

```julia
using DrWatson
@quickactivate :BEM
include(datadir("Laplace", "two_regions.jl"))
# scripts/two_regions_interface.jl
# data/elastico/aniso/*.jl  para a variante mecânica
```

#pagebreak()

== Proposta C — Dinâmica no contorno: ondas / calor transiente com análise de esquema

*Nível.* Avançado. \
*Âncoras no código.* `DIBEM`, `solve_transient`, `solve_transient_o2`,
`solve_Houbolt`, `data/Laplace/wave_propagation.jl`,
`scripts/wave_propagation.jl`, catálogo `:bar_sudden`, `:annulus`,
`:membrane_v0`, `:membrane_forced`, `:ricker`, …

=== Problema científico

Acoplar as matrizes de contorno $H,G$ a uma matriz de domínio tipo massa $M$
(DIBEM) e integrar no tempo um problema *hiperbólico* (onda escalar) *ou*
*parabólico* (difusão) com solução de referência (série de Fourier/Bessel ou
d’Alembert em barra).

O foco *não* é “gerar um vídeo”. É responder com números:

- qual a relação entre $Delta t$, $h_"elem"$, velocidade $c$ (ou difusividade) e
  estabilidade/precisão;
- como Houbolt, integrador de 2ª ordem (`solve_transient_o2`) e, se aplicável,
  1ª ordem se comportam em *dispersão*, *dissipação numérica* e custo;
- se a escolha de RBF / densidade de pontos internos domina o erro espacial.

=== O que deve ir além do script pronto

Pelo menos duas frentes:

- Comparação *head-to-head* Houbolt vs. `solve_transient_o2` (e/ou vs. referência
  modal) na *mesma* malha, com tabelas de erro em sensores e tempo de CPU.
- Estudo de CFL prático: malha fixa, varrer $Delta t$; depois $Delta t$ fixo,
  varrer $n_"div"$ e $n_"int"$.
- Caso forçado ou com pulso Ricker: espectro ou tempos de chegada vs. teoria.
- Um caso *novo* de domínio (L-shape, coroa, membrana com obstáculo) com malha
  própria + condição inicial/no contorno documentada.
- Extra de peso: acoplar exportação temporal para views Gmsh / animação
  reprodutível + métrica $L_2(Omega)$ aproximada via pontos internos ao longo do
  tempo.

=== Entregas específicas

- Formulação discreta: de $H T = G q + M f$ (ou $M accent(u, dot.double)$) ao marchador usado.
- Gráficos espaço-tempo ou snapshots em instantes teóricos de reflexão.
- Tabela de erros vs. $Delta t$ e vs. $N$ com ordem observada.
- Discussão honesta: onde o DIBEM polui a solução em relação ao erro de tempo.

=== Pontos de partida

```julia
using DrWatson
@quickactivate :BEM
include(datadir("Laplace", "Laplace_dad.jl"))
include(datadir("Laplace", "wave_propagation.jl"))
dad, meta = wave_problem(:bar_sudden; ndiv=12, n_int=6)
# scripts/wave_propagation.jl  (WAVE_CASE, WAVE_SOLVER=houbolt|diffeq|both)
```

#pagebreak()

== Proposta D — H-matrizes, custo e fidelidade em malhas grandes

*Nível.* Avançado. \
*Âncoras no código.* `H_G_Hmat`, `src/Hmat/`, `scripts/compare_hbs_hss.jl`,
`scripts/profile_assembly.jl`, `scripts/plate_large.jl`, docs de performance.

=== Problema científico

Para um problema *verificável* (Laplace $T=x$ ou elasticidade patch / placa com
furo com solução de Kirsch), construir a *fronteira de Pareto* precisão × custo.
Não basta um run com H-matriz: o trabalho é *instrumentar*, varrer $N$ e
responder com dados *quando* comprimir vale a pena.

=== O que deve ir além do script pronto

Pelo menos duas frentes:

- Montagem densa vs. H-matriz (e, se possível, mais de um formato: HODLR / HSS /
  parâmetros `atol`, `nmax`) na *mesma* geometria.
- Curvas de `rel_error` (ou erro em sensores) vs. $N$; tempo de montagem, tempo de
  solve, memória / `compression_ratio`; $N$ em que o denso se torna inviável na
  máquina usada.
- Geometria não trivial com $N$ alto (malha própria) *ou* script de benchmark
  automatizado com tabela/CSV reprodutível.
- Estudo de threads e/ou falha controlada (afrouxar `atol` até o erro saturar)
  com interpretação; perfil (`@time`, `TimerOutputs`) e gargalo (integração vs.
  álgebra).

=== Entregas específicas

- Gráfico único “erro × tempo” com famílias densa / H-matriz.
- Tabela $N$, memória, `compression_ratio`, tempos, erro.
- Texto explícito: *em que regime* a equipe recomenda H-matriz neste hardware.

=== Pontos de partida

```bash
julia --project=. scripts/compare_hbs_hss.jl
julia --project=. scripts/profile_assembly.jl
julia --project=. scripts/plate_large.jl
```

#pagebreak()

== Proposta E — Contato elástico por BEM de semi-espaço / semi-plano

*Nível.* Avançado. \
*Âncoras no código.* `src/Contact/`, `scripts/hertz_line_2d.jl`,
`scripts/contact_pohrt_li.jl`, `scripts/compare_contact_acceleration.jl`,
dados de fretting / Cattaneo–Mindlin em `data/elastico/iso/`.

=== Problema científico

Reproduzir e *estender* um contato clássico (Hertz linha 2D ou Cattaneo–Mindlin /
fretting) no formalismo de semi-espaço ou semi-plano do `BEM_gmsh`. O centro é a
física do contato (zona ativa, pressão, deslizamento), não só chamar um script.

=== O que deve ir além do script pronto

Pelo menos duas frentes:

- Pressão de contato vs. solução analítica de Hertz (erro em $p_0$, $a_"contato"$)
  e convergência com refino da malha de potencial contato.
- Estudo paramétrico (carga, raio, módulo efetivo) colapsado na forma adimensional
  de Hertz.
- Discussão do algoritmo de região de contato (ativo/inativo) e comparação de
  *duas* estratégias de aceleração do repositório, se aplicável.
- Extensão: superfície com rugosidade simples ou indentador não circular; *ou*
  carga em “ciclo” com análise de deslizamento parcial; *ou* visualização da zona
  de contato evolutiva + comparação de CPU entre variantes.

=== Entregas específicas

- Figura $p(x)/p_0$ vs. $x/a$ sobreposta à elipse de Hertz.
- Tabela de erros ($p_0$, semi-largura, força resultante) vs. malha.
- Seção sobre limites do modelo de semi-espaço / semi-plano.

=== Pontos de partida

```bash
julia --project=. scripts/hertz_line_2d.jl
julia --project=. scripts/contact_pohrt_li.jl
julia --project=. scripts/compare_contact_acceleration.jl
```

#pagebreak()

== Cronograma sugerido (a partir da 2ª metade do curso)

#table(
  columns: (auto, auto),
  inset: 6pt,
  stroke: 0.5pt + luma(200),
  [*Semana*], [*Marco*],
  [W1], [Escolha da proposta + leitura da API/scripts âncora + 1 página de plano],
  [W2], [Malha/CDC mínimas rodando; primeiro número comparável à referência],
  [W3], [Estudo paramétrico ou de convergência a meio caminho; texto da formulação],
  [W4], [Congelar resultados; relatório + README; ensaio da defesa],
)

== Integridade e uso de IA

- Cite o `BEM_gmsh`, as notas e *todas* as soluções analíticas/handbooks.
- Podem usar assistentes de código, mas a equipe deve explicar *qualquer* trecho
  na defesa. Código que a equipe não entende = trecho inválido na nota.
- Não entregue resultados que não conseguiu reproduzir do zero no ambiente da
  disciplina.

== Escolha orientada

#table(
  columns: (auto, auto, auto),
  inset: 6pt,
  stroke: 0.5pt + luma(200),
  [*Se você curte…*], [*Prefira*], [*Cuidado com*],
  [Mecânica dos sólidos / $K_I$], [A — fratura dual], [Singularidade na ponta; malha dual],
  [Materiais / interfaces], [B — multirregião], [Condicionamento com contraste alto],
  [Dinâmica / sinais], [C — ondas/transiente], [CFL escondido; erro do DIBEM],
  [HPC / algoritmos], [D — H-matrizes], [Medir mal tempo/memória],
  [Tribologia / contato], [E — Hertz / fretting], [Algoritmo de zona ativa],
)

Em dúvida, fale com o docente *antes* de W2: trocar de proposta depois do primeiro
marco custa caro.
