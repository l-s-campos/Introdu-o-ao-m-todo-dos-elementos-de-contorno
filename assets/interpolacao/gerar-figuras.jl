# Gera as figuras dos slides 03-interpolacao.typ a partir dos códigos da aula.
# Uso: julia assets/interpolacao/gerar-figuras.jl

const OUT = @__DIR__
ENV["GKSwstype"] = "100"  # GR headless (Windows / sem display)

using Plots
using Polynomials
using Interpolations
using LinearAlgebra
using FastGaussQuadrature

gr()
default(
    size = (960, 540),
    dpi = 160,
    legendfontsize = 9,
    tickfontsize = 9,
    guidefontsize = 11,
    titlefontsize = 13,
    margin = 4Plots.mm,
    lw = 2,
)

saveplot(name) = (path = joinpath(OUT, name); savefig(path); println("  saved ", name); path)

println("=== 1. Vandermonde / população da China ===")
year = [1982, 2000, 2010, 2015]
pop = [1008.18, 1262.64, 1337.82, 1374.62]
t = year .- 1980.0
y = pop
V = [t[i]^j for i = 1:4, j = 0:3]
display(V)
c = V \ y
display(c)
p = Polynomial(c)
println("p(2005-1980) = ", p(2005 - 1980))

tt = range(0, 35; length = 500)
scatter(t, y; label = "real", xlabel = "anos desde 1980",
        ylabel = "população (milhões)", title = "População da China", legend = :topleft)
plot!(tt, p.(tt); label = "interpolante", lw = 2)
saveplot("china-pop.png")

println("=== 2. Dados n=5 ===")
n = 5
t = range(-1, 1; length = n + 1)
y = @. t^2 + t + 0.05 * sin(20 * t)
scatter(t, y; label = "dados", legend = :topleft, title = "Dados, n = 5",
        xlabel = "t", ylabel = "y")
saveplot("dados-n5.png")

println("=== 3. Interpolante n=5 ===")
p = fit(t, y, n)
xx = range(-1, 1; length = 400)
scatter(t, y; label = "dados", legend = :topleft, title = "Interpolante, n = 5",
        xlabel = "t", ylabel = "y")
plot!(xx, p.(xx); label = "interpolante", lw = 2)
saveplot("interp-n5.png")

println("=== 4. Dados n=18 ===")
n = 18
t = range(-1, 1; length = n + 1)
y = @. t^2 + t + 0.05 * sin(20 * t)
scatter(t, y; label = "dados", legend = :topleft, title = "Dados, n = 18",
        xlabel = "t", ylabel = "y")
saveplot("dados-n18.png")

println("=== 5. Interpolante n=18 (oscilação) ===")
p = fit(t, y, n)
x = range(-1, 1; length = 1000)
scatter(t, y; label = "dados", legend = :topleft, title = "Interpolante, n = 18",
        xlabel = "t", ylabel = "y")
plot!(x, p.(x); label = "interpolante", lw = 2)
saveplot("interp-n18.png")

println("=== 6. Linear vs cúbico por partes ===")
p1 = LinearInterpolation(t, y)
p3 = CubicSplineInterpolation(t, y)
xx = range(first(t), last(t); length = 400)
scatter(t, y; label = "dados", legend = :topleft, title = "Polinômios por partes, n = 18",
        xlabel = "t", ylabel = "y")
plot!(xx, p1.(xx); label = "linear por partes", lw = 2)
plot!(xx, p3.(xx); label = "cúbico por partes", lw = 2)
saveplot("spline-por-partes.png")

println("=== 7. Curva paramétrica (cúbica periódica) ===")
x = [0, 0.51, 0.96, 1.06, 1.29, 1.55, 1.73, 2.13, 2.61,
     2.19, 1.76, 1.56, 1.25, 1.04, 0.58, 0]
y = [0, 0.16, 0.16, 0.43, 0.62, 0.48, 0.19, 0.18, 0,
     -0.12, -0.12, -0.29, -0.30, -0.15, -0.16, 0]
s = 0:15
itpx = scale(interpolate(x, BSpline(Cubic(Periodic(OnGrid())))), s)
itpy = scale(interpolate(y, BSpline(Cubic(Periodic(OnGrid())))), s)
ss = range(0, 15; length = 400)
plot(itpx.(ss), itpy.(ss); label = "cúbica periódica", lw = 2, aspect_ratio = 1,
     title = "Spline cúbica periódica", xlabel = "x(s)", ylabel = "y(s)")
scatter!(x, y; label = "dados")
saveplot("curva-periodica.png")

# comparação natural (pedido no exercício)
itpx_n = scale(interpolate(x, BSpline(Cubic(Natural(OnGrid())))), s)
itpy_n = scale(interpolate(y, BSpline(Cubic(Natural(OnGrid())))), s)
plot(itpx.(ss), itpy.(ss); label = "cúbica periódica", lw = 2, aspect_ratio = 1,
     title = "Periódica vs natural", xlabel = "x(s)", ylabel = "y(s)")
plot!(itpx_n.(ss), itpy_n.(ss); label = "cúbica natural", lw = 2, ls = :dash)
scatter!(x, y; label = "dados")
saveplot("curva-periodica-vs-natural.png")

println("=== 8. Interpolante equidistante n=6 ===")
f = x -> sin(exp(2 * x))
t = (0:6) ./ 6
y = f.(t)
p = fit(t, y)
xx = range(0, 1; length = 400)
plot(xx, f.(xx); label = "função", title = "Interpolante equidistante, n=6",
     legend = :bottomleft, lw = 2, xlabel = "x", ylabel = "f(x)")
scatter!(t, y; label = "nós")
plot!(xx, p.(xx); label = "interpolante", lw = 2, ls = :dash)
saveplot("interp-equidist-n6.png")

println("=== 9. Erro vs n (nós equidistantes) ===")
n = 5:5:60
err = zeros(size(n))
x = range(0, 1, length = 2001)
for (i, n) in enumerate(n)
    t = (0:n) / n
    y = f.(t)
    p = fit(t, y)
    err[i] = norm((@. f(x) - p(x)), Inf)
end
plot(n, err; yscale = :log10, xlabel = "n", ylabel = "max error",
     title = "Erro de interpolação para nós equidistantes",
     marker = :circle, label = false, lw = 2)
saveplot("erro-equidist.png")
println("  err = ", err)

println("=== 10. Função teste Runge ===")
f = x -> 1 / (x^2 + 16)
xx = range(-1, 1; length = 400)
plot(xx, f.(xx); title = "Função teste", label = false, lw = 2,
     xlabel = "x", ylabel = "f(x)")
saveplot("runge-funcao.png")

println("=== 11. Erro graus baixos ===")
x = range(-1, 1; length = 2501)
plt = plot(xlabel = "x", ylabel = "|f-p|", yscale = :log10, title = "Erro para graus baixos",
           ylims = (1e-20, 1), legend = :bottomright)
for n in 4:4:12
    tt = range(-1, 1; length = n + 1)
    pp = fit(tt, f.(tt))
    plot!(plt, x, abs.(f.(x) .- pp.(x)); label = "grau $n", lw = 2)
end
plt
saveplot("runge-erro-baixo.png")

println("=== 12. Erro graus altos ===")
plt = plot(xlabel = "x", ylabel = "|f-p|", yscale = :log10, title = "Erro para graus altos",
           ylims = (1e-20, 1), legend = :bottomright)
for n in @. 12 + 15 * (1:3)
    tt = range(-1, 1; length = n + 1)
    pp = fit(tt, f.(tt))
    plot!(plt, x, abs.(f.(x) .- pp.(x)); label = "grau $n", lw = 2)
end
plt
saveplot("runge-erro-alto.png")

println("=== 13. Erro Chebyshev ===")
x = range(-1, 1; length = 2001)
plt = plot(xlabel = "x", ylabel = "|f-p|", yscale = :log10, title = "Erro com os pontos de Chebyshev",
           ylims = (1e-20, 1), legend = :bottomright)
for n in [4, 10, 16, 40]
    tt = [-cos(pi * k / n) for k in 0:n]
    pp = fit(tt, f.(tt))
    plot!(plt, x, abs.(f.(x) .- pp.(x)); label = "grau $n", lw = 2)
end
plt
saveplot("chebyshev-erro.png")

println("=== 14. Exercício: equidistante vs Chebyshev ===")
f = x -> 1 / (1 + 25x^2)
x_test = range(-1, 1; length = 2501)
ns = [8, 16, 32]
err_eq = Float64[]
err_ch = Float64[]
for n in ns
    t_eq = range(-1, 1; length = n + 1)
    t_ch = [-cos(pi * k / n) for k in 0:n]
    p_eq = fit(collect(t_eq), f.(t_eq))
    p_ch = fit(t_ch, f.(t_ch))
    push!(err_eq, maximum(abs, f.(x_test) .- p_eq.(x_test)))
    push!(err_ch, maximum(abs, f.(x_test) .- p_ch.(x_test)))
end
println("  n = ", ns)
println("  err_eq = ", err_eq)
println("  err_ch = ", err_ch)

n = 16
t_eq = range(-1, 1; length = n + 1)
t_ch = [-cos(pi * k / n) for k in 0:n]
p_eq = fit(collect(t_eq), f.(t_eq))
p_ch = fit(t_ch, f.(t_ch))
xx = range(-1, 1; length = 800)
plot(xx, f.(xx); label = "f", lw = 2, title = "n = 16", xlabel = "x", ylabel = "y")
plot!(xx, p_eq.(xx); label = "equidistante", ls = :dash, lw = 2)
plot!(xx, p_ch.(xx); label = "Chebyshev", ls = :dot, lw = 2)
saveplot("runge-n16-eq-vs-cheb.png")

plot(ns, err_eq; yscale = :log10, marker = :circle, label = "equidistante",
     xlabel = "n", ylabel = "‖f-p‖∞", lw = 2, title = "Erro máximo: equidistante vs Chebyshev")
plot!(ns, err_ch; marker = :square, label = "Chebyshev", lw = 2)
saveplot("runge-erro-eq-vs-cheb.png")

println("=== 15. Gauss: exp(x) e exp(sin x) ===")
exato = exp(1) - exp(-1)
xg, wg = gausslegendre(3)
println("  nós, pesos = ", (xg, wg))
fexp(x) = exp(x)
In = dot(wg, fexp.(xg))
println("  exp(x): In = ", In, "  erro = ", exato - In)

fsin(x) = exp(sin(x))
In = dot(wg, fsin.(xg))
println("  exp(sin x): In = ", In)

xx = range(-1, 1; length = 400)
pA = plot(xx, exp.(xx); fillrange = 0, fillalpha = 0.25, label = false,
          xlabel = "x", ylabel = "exp(x)", ylims = (0, 2.7), lw = 2)
pB = plot(xx, exp.(sin.(xx)); fillrange = 0, fillalpha = 0.25, label = false,
          xlabel = "x", ylabel = "exp(sin(x))", ylims = (0, 2.7), lw = 2)
plot(pA, pB; layout = (2, 1), size = (700, 500))
saveplot("integrandos-exp.png")

println("=== 16. Trapézio ===")
function trapezoidal(f, a, b, n)
    h = (b - a) / n
    t = range(a, b, length = n + 1)
    y = f.(t)
    T = h * (sum(y[2:n]) + 0.5 * (y[1] + y[n+1]))
    return T, t, y
end

f = x -> exp(sin(7 * x))
a = 0
b = 2
Integral = 2.6632197827615394
T, t, y = trapezoidal(f, a, b, 40)
@show (T, Integral - T)

ntrap = [10^k for k in 1:5]
err = []
for n in ntrap
    T, t, y = trapezoidal(f, a, b, n)
    push!(err, Integral - T)
end
foreach(args -> println(args[1], "\t", args[2]), zip(ntrap, err))

println("=== 17. Trapézio vs Gauss  1/(1+4x²) ===")
f = x -> 1 / (1 + 4 * x^2)
exato = atan(2)
n = 8:4:96
errT = zeros(size(n))
errG = zeros(size(n))
for (k, n) in enumerate(n)
    errT[k] = abs(exato - trapezoidal(f, -1, 1, n)[1])
    xg, wg = gausslegendre(n)
    errG[k] = abs(exato - dot(wg, f.(xg)))
end
errT[iszero.(errT)] .= NaN
errG[iszero.(errG)] .= NaN
plot(collect(n), errT; yscale = :log10, xlabel = "nós", ylabel = "erro",
     title = "integração numérica  1/(1+4x²)", ylims = (1e-16, 1),
     marker = :circle, label = "trapézio", lw = 2)
plot!(collect(n), errG; marker = :circle, label = "Gauss-Legendre", lw = 2)
saveplot("quad-4x2.png")

println("=== 18. Trapézio vs Gauss  1/(1+16x²) ===")
f = x -> 1 / (1 + 16 * x^2)
exato = atan(4) / 2
n = 8:4:96
errT = zeros(size(n))
errG = zeros(size(n))
for (k, n) in enumerate(n)
    errT[k] = abs(exato - trapezoidal(f, -1, 1, n)[1])
    xg, wg = gausslegendre(n)
    errG[k] = abs(exato - dot(wg, f.(xg)))
end
errT[iszero.(errT)] .= NaN
errG[iszero.(errG)] .= NaN
plot(collect(n), errT; yscale = :log10, xlabel = "nós", ylabel = "erro",
     title = "integração numérica  1/(1+16x²)", ylims = (1e-16, 1),
     marker = :circle, label = "trapézio", lw = 2)
plot!(collect(n), errG; marker = :circle, label = "Gauss-Legendre", lw = 2)
saveplot("quad-16x2.png")

println("=== 19. sinhtrans ===")
function sinhtrans(u, a, b)
    μ = 1 / 2 * (asinh((1 + a) / b) + asinh((1 - a) / b))
    η = 1 / 2 * (asinh((1 + a) / b) - asinh((1 - a) / b))
    x = a .+ b * sinh.(μ * u .- η)
    J = b * μ * cosh.(μ * u .- η)
    x, J
end

function plot_quasising_sinh(a = 0.0, b = 0.05; g = nothing, n = 800)
    g === nothing && (g = ξ -> 1.0)
    f(ξ) = g(ξ) / ((ξ - a)^2 + b^2)
    ξ = range(-1, 1; length = n)
    s = range(-1, 1; length = n)
    ξs, J = sinhtrans(collect(s), a, b)
    before = f.(ξ)
    after = @. f(ξs) * abs(J)
    p1 = plot(ξ, before; lw = 2, label = "f(ξ)",
              xlabel = "ξ", ylabel = "integrando",
              title = "Antes (quase-singular)", legend = :topright)
    vline!(p1, [a]; ls = :dash, color = :gray, label = "a = $a")
    p2 = plot(s, after; lw = 2, label = "f(ξ(s)) · |J|",
              xlabel = "s", ylabel = "integrando transformado",
              title = "Depois (× Jacobiano sinh)", legend = :topright)
    _, i = findmin(abs.(ξs .- a))
    vline!(p2, [s[i]]; ls = :dash, color = :gray, label = "ξ(s) ≈ a")
    plot(p1, p2; layout = (1, 2), size = (900, 350))
end

x = range(-1, 1; length = 100)
xt, jac = sinhtrans(x, 0.5, 0.01)
p1 = plot(collect(x), xt; xlabel = "s", ylabel = "ξ", label = false, lw = 2, title = "mapa ξ(s)")
p2 = plot(collect(x), jac; xlabel = "s", ylabel = "dξ/ds", label = false, lw = 2, title = "Jacobiano")
plot(p1, p2; layout = (2, 1), size = (700, 500))
saveplot("sinh-mapa-jac.png")

plot_quasising_sinh(0.0, 0.05)
saveplot("quasising-a0-b005.png")

plot_quasising_sinh(0.5, 0.02)
saveplot("quasising-a05-b002.png")

println("=== 20. Gauss + sinhtrans na integral 1/(1+16x²) ===")
f = x -> 1 / (1 + 16 * x^2)
exato = atan(4) / 2
n = 8
xg, wg = gausslegendre(n)
xt, J = sinhtrans(xg, 0, 1 / 4)
errG = abs(exato - dot(wg, f.(xg)))
errGt = abs(exato - dot(wg, J .* f.(xt)))
println("  errG = ", errG, "  errGt = ", errGt)

println("=== 21. Tabela b↓ + plot b=1e-3 ===")
Iex(b) = (2 / b) * atan(1 / b)
fb(ξ, b) = 1 / (ξ^2 + b^2)
bs = (1e-1, 1e-2, 1e-3)
ns = (8, 16, 32)
println("b\tn\terr_cru\terr_sinh")
for b in bs, n in ns
    s, w = gausslegendre(n)
    e_cru = abs(Iex(b) - dot(w, fb.(s, b)))
    ξ, J = sinhtrans(s, 0.0, b)
    e_sinh = abs(Iex(b) - dot(w, fb.(ξ, b) .* abs.(J)))
    println("$b\t$n\t$e_cru\t$e_sinh")
end
plot_quasising_sinh(0.0, 1e-3)
saveplot("quasising-b001.png")

println("=== 22. Monegato (definição) ===")
function Monegato(t, s0, q::Integer = 5; atol = 1e-14)
    t = float.(t)
    s0 = float(s0)
    if isapprox(s0, 1; atol = atol)
        q ≥ 2 || throw(ArgumentError("Sato no extremo +1 exige q ≥ 2"))
        u = @. 1 - t
        s = @. 1 - u^q / 2^(q - 1)
        ds = @. q * u^(q - 1) / 2^(q - 1)
        return s, ds
    elseif isapprox(s0, -1; atol = atol)
        q ≥ 2 || throw(ArgumentError("Sato no extremo -1 exige q ≥ 2"))
        u = @. 1 + t
        s = @. -1 + u^q / 2^(q - 1)
        ds = @. q * u^(q - 1) / 2^(q - 1)
        return s, ds
    end
    (-1 < s0 < 1) || throw(ArgumentError("s0 deve estar em (-1,1) ou ±1"))
    (q ≥ 3 && isodd(q)) || throw(ArgumentError(
        "Monegato–Sloan interior exige grau ímpar q = 3,5,7,… (recebeu q=$q)"))
    aq = (1 + s0)^(1 / q)
    bq = (1 - s0)^(1 / q)
    δ = (1 / 2)^q * (aq + bq)^q
    t0 = (aq - bq) / (aq + bq)
    u = @. t - t0
    s = @. s0 + δ * u^q
    ds = @. q * δ * u^(q - 1)
    return s, ds
end
s_m, ds_m = Monegato(range(-1, 1; length = 5), 1, 4)
println("  Monegato s0=1 q=4  s = ", s_m)
println("  ds = ", ds_m)

println("\nOK — figuras em ", OUT)
