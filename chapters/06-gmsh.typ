// Gmsh e malhas para o BEM
// Capítulo de apoio ao código BEM_gmsh

= Gmsh e malhas para o BEM
<gmsh>

O #link("https://gmsh.info/")[Gmsh] é um gerador de malhas de código aberto amplamente usado em elementos finitos e, cada vez mais, em elementos de contorno. No curso usamos o repositório #link("https://github.com/l-s-campos/BEM_gmsh")[`BEM_gmsh`], que acopla o Gmsh ao fluxo de trabalho em Julia: geometria → grupos físicos (condições de contorno) → malha `.msh` → `format2d` / `format3d` → montagem e solução.

== Por que Gmsh no BEM?

No BEM clássico 2D, a malha é *apenas o contorno* (curvas). Ainda assim, o Gmsh é útil porque:

1. descreve geometrias com furos, arcos, splines e múltiplas regiões de forma robusta;
2. associa *nomes de grupos físicos* às condições de contorno (Dirichlet / Neumann);
3. gera malhas de ordem 1 ou 2 (elementos lineares ou quadráticos);
4. permite refinamento local (campos de tamanho perto de cantos e furos);
5. integra-se ao ONELAB (interface gráfica de parâmetros + pós-processamento).

Para o BEM, o domínio volumétrico da malha 2D do Gmsh serve sobretudo para *pontos internos* e para métodos de domínio (DIBEM). Os elementos de contorno são extraídos das *curvas físicas* (dimensão 1).

== Instalação rápida

No ambiente do projeto `BEM_gmsh`:

```julia
using Pkg
Pkg.activate(".")   # pasta do BEM_gmsh
Pkg.instantiate()

using DrWatson
@quickactivate :BEM
```

O pacote Julia `Gmsh.jl` já vem como dependência e embute (ou localiza) o executável do Gmsh. Para a interface gráfica completa, instale também o #link("https://gmsh.info/#Download")[Gmsh do sistema].

== Anatomia de um arquivo `.geo`

Um script Gmsh típico para o BEM 2D tem quatro blocos:

1. *Pontos e curvas* — geometria;
2. *Superfície* — domínio (para pontos internos / DIBEM);
3. *Grupos físicos* — nomes que codificam as CDCs;
4. *Geração da malha* — `Mesh 2`, opcionalmente `RecombineMesh` e ordem.

Exemplo mínimo (quadrado unitário, campo $T = x$):

```geo
// quadrado_bem.geo
lc = 0.1;

Point(1) = {0, 0, 0, lc};
Point(2) = {1, 0, 0, lc};
Point(3) = {1, 1, 0, lc};
Point(4) = {0, 1, 0, lc};

Line(1) = {1, 2}; // base
Line(2) = {2, 3}; // direita
Line(3) = {3, 4}; // topo
Line(4) = {4, 1}; // esquerda

Curve Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};

Transfinite Curve {1, 2, 3, 4} = 12;
Transfinite Surface{1};
Recombine Surface{1};

// Convenção BEM_gmsh (Laplace): "tipo;valor"
// tipo 0 = Dirichlet (T), tipo 1 = Neumann (q = -k ∂T/∂n)
Physical Curve("1;0")  = {1, 3}; // isolado
Physical Curve("1;-1") = {2};    // q = -1  (∂T/∂n = +1 se k=1)
Physical Curve("0;0")  = {4};    // T = 0
Physical Surface("Domain") = {1};

Mesh 2;
```

Gere a malha no terminal:

```bash
gmsh -2 quadrado_bem.geo -o quadrado_bem.msh
```

ou, de dentro do Julia / `BEM_gmsh`, use os geradores prontos (recomendado no curso):

```julia
using DrWatson
@quickactivate :BEM
include(datadir("Laplace", "Laplace_dad.jl"))

msh = quadrado(ndiv=20, show=false)  # grava .msh e devolve o caminho
```

== Convenção de condições de contorno

Os *nomes* dos grupos físicos carregam o tipo e o valor da CDC. Isso evita arquivos auxiliares de contorno.

=== Laplace / Poisson

#table(
  columns: (auto, auto),
  inset: 8pt,
  stroke: 0.5pt + luma(200),
  [*Nome do grupo*], [*Significado*],
  [`"0;T"`], [Dirichlet: temperatura (ou potencial) preescrita $T$],
  [`"1;q"`], [Neumann: fluxo $q = -k partial T \/ partial n$],
)

Exemplos: `"0;100"`, `"1;0"` (isolado), `"1;-1"`.

=== Elasticidade 2D

Cada aresta usa *quatro* tokens `tx;ux;ty;uy`:

- token ímpar (`tx`, `ty`): `0` = deslocamento preescrito, `1` = tração preescrita;
- token par (`ux`, `uy`): valor correspondente.

#table(
  columns: (auto, auto),
  inset: 8pt,
  stroke: 0.5pt + luma(200),
  [*Nome*], [*Significado*],
  [`"0;0;0;0"`], [engaste ($u_x = u_y = 0$)],
  [`"1;0;1;0"`], [livre de tração],
  [`"1;P;1;0"`], [tração $t_x = P$, $t_y = 0$],
  [`"0;0;1;0"`], [rolo vertical ($u_x = 0$, $t_y = 0$)],
)

== Do arquivo `.msh` ao `BEMdata`

O fluxo canônico no `BEM_gmsh` é:

```julia
using DrWatson
@quickactivate :BEM
include(datadir("Laplace", "Laplace_dad.jl"))

props = Laplace(1.0)                 # k = 1
msh   = quadrado(ndiv=20, show=false)
dad   = format2d(msh, props)         # lê grupos físicos → CDCs

attach_analytical!(dad, ana_laplace_linear(; direction=SA[1.0, 0.0]))
H_G_full_direct(dad, 20)             # monta H e G (densa)
solve(dad)

@show rel_error(dad)
plot_geo(dad)
```

Palavras-chave úteis de `format2d`:

- `tipo=1|2|3` — ordem do elemento de contorno (linear, quadrático, cúbico), quando compatível com a malha;
- `pontointerno=true|false` — usa nós internos da malha de superfície para pós-processamento / DIBEM;
- `finalize=false, reopen=false` — quando a sessão Gmsh já está aberta (pipelines ONELAB).

Para problemas grandes:

```julia
H_G_Hmat(dad; atol=1e-6, nmax=32)   # montagem hierárquica
solve(dad)
@show compression_ratio(dad.H)
```

== Gerando geometria em Julia (API Gmsh)

Além de arquivos `.geo`, o `BEM_gmsh` constrói malhas pela API Julia. Esqueleto de um anel (quarto de coroa):

```julia
using DrWatson
@quickactivate :BEM

function meu_quarto_circulo(; ndiv=12, ri=1.0, re=2.0, nome="meu_qc")
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    gmsh.model.add(nome)
    lc = (re - ri) / max(ndiv, 4)

    c  = gmsh.model.geo.addPoint(0.0, 0.0, 0.0, lc)
    p1 = gmsh.model.geo.addPoint(ri, 0.0, 0.0, lc)
    p2 = gmsh.model.geo.addPoint(re, 0.0, 0.0, lc)
    p3 = gmsh.model.geo.addPoint(0.0, re, 0.0, lc)
    p4 = gmsh.model.geo.addPoint(0.0, ri, 0.0, lc)

    lbot  = gmsh.model.geo.addLine(p1, p2)
    aout  = gmsh.model.geo.addCircleArc(p2, c, p3)
    lleft = gmsh.model.geo.addLine(p3, p4)
    ain   = gmsh.model.geo.addCircleArc(p4, c, p1)

    cl = gmsh.model.geo.addCurveLoop([lbot, aout, lleft, ain])
    s  = gmsh.model.geo.addPlaneSurface([cl])
    gmsh.model.geo.synchronize()

    for crv in (lbot, aout, lleft, ain)
        gmsh.model.mesh.setTransfiniteCurve(crv, ndiv)
    end
    gmsh.model.mesh.setTransfiniteSurface(s)

    gmsh.model.addPhysicalGroup(1, [lbot, lleft], -1, "1;0")   # isolado
    gmsh.model.addPhysicalGroup(1, [ain], -1, "0;100")         # T = 100
    gmsh.model.addPhysicalGroup(1, [aout], -1, "1;-200")       # q = -200
    gmsh.model.addPhysicalGroup(2, [s], -1, "Domain")

    gmsh.model.mesh.generate(2)
    out = datadir("Laplace", nome * ".msh")
    mkpath(dirname(out))
    gmsh.write(out)
    gmsh.finalize()
    return out
end
```

Geradores prontos no repositório (vale estudar o código-fonte):

- `data/Laplace/Laplace_dad.jl` — `quadrado`, `placa_com_furo`, `quadrado_elasticity`;
- `data/Laplace/potencial_problems.jl` — Laquini, quarto de círculo, Moulton;
- `data/elastico/iso/*.jl` — cilindro pressurizado, placa com furo, viga em balanço, trinca central.

== Propriedades geométricas via contorno

O módulo `geometric_props` do `BEM_gmsh` calcula perímetro, área e centróide *só com integrais de contorno* — o mesmo espírito do capítulo “Indo para 2D”:

```julia
using DrWatson
@quickactivate :BEM
include(datadir("Laplace", "Laplace_dad.jl"))

msh = quadrado(ndiv=30, show=false)
dad = format2d(msh, Laplace(1.0); pontointerno=false)

gp = geometric_props(dad)
@show gp.perimeter gp.area gp.centroid
```

Em 3D (`format3d` + malha de superfície fechada):

```julia
# exemplo conceitual — requer um .geo/.msh 3D (ex.: cubo)
# dad3 = format3d(datadir("Laplace", "cubo.msh"), Laplace(1.0))
# gp3  = geometric_props(dad3)   # área superficial, volume, centróide
```

== ONELAB: Gmsh como interface do solver

O `BEM_gmsh` publica um cliente ONELAB. Parâmetros (`ndiv`, $k$, tipo de montagem, …) ficam no painel do Gmsh; ao clicar em *Run*, a malha é gerada, o BEM resolve e as views de pós-processamento voltam para a GUI.

```julia
using BEM
dad = bem_onelab_run("data/Laplace/onelab_square.geo"; ndiv=12)
# interativo:
# bem_onelab_gui("data/Laplace/onelab_square.geo")
```

Na linha de comando:

```bash
julia --project=. scripts/bem_onelab.jl data/Laplace/onelab_square.geo
julia --project=. scripts/bem_onelab.jl data/Laplace/onelab_square.geo --gui
```

Para registrar como solver do Gmsh do sistema (*Tools → Options → Solver*):

`julia --project=C:/caminho/BEM_gmsh C:/caminho/BEM_gmsh/scripts/bem_onelab.jl %s`

== Boas práticas de malha para o BEM

1. *Contorno fechado e orientado* — curvas no sentido anti-horário no contorno externo (normal para fora).
2. *Furos* — contornos internos no sentido horário.
3. *Cantos com CDC mista* — elementos descontínuos (padrão do código) evitam conflito de nós compartilhados.
4. *Refino local* — use `Distance` + `Threshold` perto de furos e reentrâncias (`placa_com_furo` é um bom modelo).
5. *Ordem geométrica = ordem de interpolação* — malha de ordem 2 combina com elementos quadráticos.
6. *Estudo de convergência* — varie `ndiv` (ou `lc`) e reporte erro médio / máximo / $L_2$ (capítulo de medidas de erro).
7. *Não exagere nos pontos internos* — eles encarecem o DIBEM; comece com malha de contorno e poucos internos.

== Exercícios

1. Reproduza o quadrado unitário com $T = x$ a partir de um `.geo` próprio (não use `quadrado`). Compare `rel_error` com `ndiv = 8, 16, 32`.
2. Modele uma placa retangular com um furo circular deslocado do centro. Imponha $T = 0$ à esquerda, $T = 1$ à direita e isolamento no restante. Plote o campo com `plot_geo` e exporte com `export_results_to_gmsh`.
3. Calcule perímetro, área e centróide da figura do item 2 com `geometric_props` e confira com a fórmula analítica (retângulo menos círculo).
4. Abra `onelab_square.geo` no modo GUI e varie `ndiv`. Registre o tempo de montagem densa vs. `assembly=hmat` para o maior `ndiv` que sua máquina permitir.

== Leitura e código de referência

- Documentação do projeto: `BEM_gmsh/docs` (EN e pt-BR)
- Começando: `docs/src/pt-br/getting_started.md`
- ONELAB: `docs/src/api/onelab.md`
- Exemplos: `scripts/intro.jl`, `data/examples/`, `scripts/bem_onelab.jl`
- Manual Gmsh: #link("https://gmsh.info/doc/texinfo/gmsh.html")[gmsh.info/doc]
