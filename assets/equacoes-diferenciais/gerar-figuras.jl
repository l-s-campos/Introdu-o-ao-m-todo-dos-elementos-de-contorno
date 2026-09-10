# Gera as figuras dos slides 04-equacoes-diferenciais.typ
# Uso: julia assets/equacoes-diferenciais/gerar-figuras.jl

const OUT = @__DIR__
ENV["GKSwstype"] = "100"

using DifferentialEquations
using Plots
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

println("=== 1. PVI Tsit5 ===")
f = (u, p, t) -> sin((t + u)^2)
u0 = 1.0
tspan = (0.0, 4.0)
ivp = ODEProblem(f, u0, tspan)
sol = solve(ivp, Tsit5())
plot(sol.t, sol.u; xlabel = "t", ylabel = "u(t)", label = "solução", lw = 2)
scatter!(sol.t, sol.u; label = "nós adaptativos", ms = 3)
saveplot("pvi-tsit5.png")
println("  sol(1.0) = ", sol(1.0))

println("=== 2. Euler ===")
function euler(ivp, n)
    a, b = ivp.tspan
    h = (b - a) / n
    t = [a + i * h for i in 0:n]
    u0 = ivp.u0 isa Number ? float(ivp.u0) : float.(ivp.u0)
    U = Vector{typeof(u0)}(undef, n + 1)
    U[1] = u0
    for i in 1:n
        U[i+1] = U[i] + h * ivp.f(U[i], ivp.p, t[i])
    end
    return t, U
end

f = (u, p, t) -> sin((t + u)^2)
ivp = ODEProblem(f, -1.0, (0.0, 4.0))
t20, u20 = euler(ivp, 20)
t50, u50 = euler(ivp, 50)
u_ref = solve(ivp, Tsit5(); reltol = 1e-14, abstol = 1e-14)
plot(t20, u20; marker = :circle, label = "Euler n=20", xlabel = "t", ylabel = "u", lw = 2)
plot!(t50, u50; marker = :circle, label = "Euler n=50", lw = 2)
plot!(t50, u_ref.(t50); color = :black, lw = 2, label = "referência")
saveplot("euler-n20-n50.png")

println("=== 3. Erro Euler ===")
ns = [round(Int, 5 * 10^k) for k in 0:0.5:3]
err_E = Float64[]
for n in ns
    local t, U
    t, U = euler(ivp, n)
    push!(err_E, maximum(abs, u_ref.(t) .- U))
end
plot(ns, err_E; xscale = :log10, yscale = :log10, marker = :circle,
     label = "Euler", xlabel = "n", ylabel = "‖e‖∞", lw = 2)
plot!(ns, err_E[1] .* (ns[1] ./ ns); ls = :dash, label = "O(1/n)")
saveplot("euler-erro.png")
println("  ns = ", ns)
println("  err_E = ", err_E)

println("=== 4. RK4 vs Euler ===")
function rk4(ivp, n)
    a, b = ivp.tspan
    h = (b - a) / n
    t = [a + i * h for i in 0:n]
    u0 = ivp.u0 isa Number ? float(ivp.u0) : float.(ivp.u0)
    U = Vector{typeof(u0)}(undef, n + 1)
    U[1] = u0
    for i in 1:n
        ui, ti = U[i], t[i]
        k1 = h * ivp.f(ui, ivp.p, ti)
        k2 = h * ivp.f(ui + k1 / 2, ivp.p, ti + h / 2)
        k3 = h * ivp.f(ui + k2 / 2, ivp.p, ti + h / 2)
        k4 = h * ivp.f(ui + k3, ivp.p, ti + h)
        U[i+1] = ui + (k1 + 2k2 + 2k3 + k4) / 6
    end
    return t, U
end

err_R = Float64[]
for n in ns
    local t, U
    t, U = rk4(ivp, n)
    push!(err_R, maximum(abs, u_ref.(t) .- U))
end
plot(ns, err_E; xscale = :log10, yscale = :log10, marker = :circle, label = "Euler p≈1", lw = 2)
plot!(ns, err_R; marker = :square, label = "RK4 p≈4", lw = 2)
plot!(ns, err_E[1] .* (ns[1] ./ ns); ls = :dash, label = "O(n⁻¹)")
plot!(ns, err_R[1] .* (ns[1] ./ ns) .^ 4; ls = :dot, label = "O(n⁻⁴)")
saveplot("euler-vs-rk4.png")
println("  err_R = ", err_R)

println("=== 5. Pêndulos acoplados ===")
function couple(u, p, t)
    γ, L, k = p
    g = 9.8
    du = similar(u)
    du[1] = u[3]
    du[2] = u[4]
    du[3] = -γ * u[3] - (g / L) * sin(u[1]) + k * (u[2] - u[1])
    du[4] = -γ * u[4] - (g / L) * sin(u[2]) + k * (u[1] - u[2])
    return du
end

u0 = [1.25, -0.5, 0.0, 0.0]
tspan = (0.0, 50.0)
γ, L = 0.0, 0.5
sol0 = solve(ODEProblem(couple, u0, tspan, [γ, L, 0.0]), Tsit5())
sol1 = solve(ODEProblem(couple, u0, tspan, [γ, L, 1.0]), Tsit5())
plot(sol0.t, [u[1] for u in sol0.u]; label = "θ1, k=0", xlims = (20, 50), lw = 2)
plot!(sol0.t, [u[2] for u in sol0.u]; label = "θ2, k=0", lw = 2)
plot!(sol1.t, [u[1] for u in sol1.u]; label = "θ1, k=1", ls = :dash, lw = 2)
plot!(sol1.t, [u[2] for u in sol1.u]; label = "θ2, k=1", ls = :dash, lw = 2,
      xlabel = "t", title = "Pêndulos acoplados")
saveplot("pendulos-acoplados.png")

println("=== 6. diffmat + PVC Laplace 1D ===")
function diffmat(n, xspan)
    a, b = xspan
    h = (b - a) / n
    x = [a + i * h for i in 0:n]
    dp = fill(0.5 / h, n)
    dm = fill(-0.5 / h, n)
    Dx = diagm(-1 => dm, 1 => dp)
    Dx[1, 1:3] = [-1.5, 2, -0.5] / h
    Dx[n+1, n-1:n+1] = [0.5, -2, 1.5] / h
    d0 = fill(-2 / h^2, n + 1)
    dp2 = ones(n) / h^2
    Dxx = diagm(-1 => dp2, 0 => d0, 1 => dp2)
    Dxx[1, 1:4] = [2, -5, 4, -1] / h^2
    Dxx[n+1, n-2:n+1] = [-1, 4, -5, 2] / h^2
    return x, Dx, Dxx
end

n = 40
x, Dx, Dxx = diffmat(n, (-1.0, 1.0))
A = copy(Dxx)
A[1, :] .= 0
A[1, 1] = 1
A[end, :] .= 0
A[end, end] = 1
rhs = zeros(n + 1)
rhs[1] = 100
rhs[end] = 0
T = A \ rhs
T_ex = @. 50 * (1 - x)
plot(x, T_ex; label = "exato", lw = 2, xlabel = "x", ylabel = "T")
scatter!(x, T; label = "MDF", ms = 3)
saveplot("pvc-laplace1d.png")
@show maximum(abs, T - T_ex)

println("=== 7. Método das linhas: calor ===")
function heat_rhs!(du, u, p, t)
    Dxx, κ, ua, ub = p
    u[1] = ua
    u[end] = ub
    mul!(du, Dxx, u)
    du .*= κ
    du[1] = 0
    du[end] = 0
    return nothing
end

n = 100
x, Dx, Dxx = diffmat(n, (-1.0, 1.0))
κ = 1.0
ua, ub = 0.0, 2.0
u0 = @. 1 + sin(pi * x / 2) + 3 * (1 - x^2) * exp(-4 * x^2)
u0[1] = ua
u0[end] = ub
prob = ODEProblem(heat_rhs!, u0, (0.0, 0.75), (Dxx, κ, ua, ub))
sol = solve(prob, Tsit5())
plt = plot(xlabel = "x", ylabel = "u(x,t)", title = "Calor / difusão", legend = :topleft)
for tt in 0:0.1:0.7
    plot!(plt, x, sol(tt); label = "t=$tt", lw = 2)
end
plt
saveplot("calor-linhas.png")

println("=== 8. Onda (esqueleto, sem plot na aula) ===")
function wave_rhs!(dy, y, p, t)
    Dxx, c2, nnode = p
    u = @view y[1:nnode]
    v = @view y[nnode+1:end]
    du = @view dy[1:nnode]
    dv = @view dy[nnode+1:end]
    u[1] = 0.0
    u[end] = 0.0
    v[1] = 0.0
    v[end] = 0.0
    du .= v
    mul!(dv, Dxx, u)
    dv .*= c2
    du[1] = 0.0
    du[end] = 0.0
    dv[1] = 0.0
    dv[end] = 0.0
    return nothing
end

function solve_wave_mdf(; n = 100, c = 1.0, tspan = (0.0, 2.0))
    x, Dx, Dxx = diffmat(n, (0.0, 1.0))
    nnode = length(x)
    u0 = sin.(pi .* x)
    v0 = zeros(nnode)
    y0 = vcat(u0, v0)
    prob = ODEProblem(wave_rhs!, y0, tspan, (Dxx, c^2, nnode))
    sol = solve(prob, Tsit5(); reltol = 1e-8, abstol = 1e-8)
    return sol, x
end
solW, xW = solve_wave_mdf()
uex(x, t) = sin(pi * x) * cos(pi * t)
eW = maximum(abs, solW(1.0)[1:length(xW)] .- uex.(xW, 1.0))
println("  wave e∞(t=1) = ", eW)

println("=== 9. BEM 1D calor ===")
Tstar(x, xd) = -0.5 * abs(x - xd)
Qstar(x, xd) = x == xd ? 0.0 : -0.5 * sign(x - xd)

function bem1d_HG(x0, xf)
    L = xf - x0
    H = [0.5 -0.5; -0.5 0.5]
    G = L * [0.0 -0.5; 0.5 0.0]
    return H, G, L
end

function bem1d_M(x0, xf; κ = 1.0, nq = 6)
    nc = 2
    ξ, ŵ = gausslegendre(nq)
    jac = (xf - x0) / 2
    xg = collect(@. (x0 + xf) / 2 + jac * ξ)
    w = collect(ŵ .* jac)
    ni = nq
    xd = vcat([x0, xf], xg)
    N = nc + ni
    Mg = zeros(N, ni)
    for p in 1:N, j in 1:ni
        Mg[p, j] = w[j] / κ * Tstar(xg[j], xd[p])
    end
    H, G, _ = bem1d_HG(x0, xf)
    Hb = zeros(N, 2)
    Gb = zeros(N, 2)
    Hb[1:2, :] .= H
    Gb[1:2, :] .= G
    for p in 3:N
        xdp = xd[p]
        Hb[p, 1] = -Qstar(x0, xdp)
        Hb[p, 2] = Qstar(xf, xdp)
        Gb[p, 1] = -Tstar(x0, xdp)
        Gb[p, 2] = Tstar(xf, xdp)
    end
    return (; x0, xf, xg, w, xd, Mg, Hb, Gb, nc, ni, N, κ)
end

function bem_heat_factor(mesh)
    (; Gb, Mg, ni, N) = mesh
    rows = vcat(1, 3:N)
    nz = 1 + ni
    A = zeros(nz, nz)
    for (eq, p) in enumerate(rows)
        A[eq, 1] = -Gb[p, 1]
        A[eq, 2:end] .= Mg[p, :]
    end
    return lu(A), rows
end

function bem_heat_rhs!(dy, y, p, t)
    mesh, TL, F, rows = p
    (; Hb, nc, ni, xg) = mesh
    Tf = y[1]
    Tg = y[2:end]
    @assert length(Tg) == ni
    r = zeros(1 + ni)
    for (eq, pidx) in enumerate(rows)
        r[eq] = -(Hb[pidx, 1] * TL + Hb[pidx, 2] * Tf)
        if pidx > nc
            r[eq] -= Tg[pidx-nc]
        end
    end
    z = F \ r
    Tdotg = z[2:end]
    dy[1] = Tdotg[argmax(xg)]
    dy[2:end] .= Tdotg
    return nothing
end

function solve_bem_heat(; L = 1.0, κ = 1.0, TL = 0.0, tspan = (0.0, 0.5), nq = 8,
    T0 = nothing)
    mesh = bem1d_M(0.0, L; κ = κ, nq = nq)
    F, rows = bem_heat_factor(mesh)
    if T0 === nothing
        T0 = x -> TL + (0.0 - TL) * (x / L)
    end
    y0 = vcat(T0(L), T0.(mesh.xg))
    prob = ODEProblem(bem_heat_rhs!, y0, tspan, (mesh, TL, F, rows))
    sol = solve(prob, Tsit5(); reltol = 1e-8, abstol = 1e-8)
    return sol, mesh
end

function profile_at(solB, mesh, TL, t)
    y = solB(t)
    Tf, Tg = y[1], y[2:end]
    x = vcat(mesh.x0, mesh.xg, mesh.xf)
    T = vcat(TL, Tg, Tf)
    perm = sortperm(x)
    return x[perm], T[perm]
end

L, κ = 1.0, 1.0
λ0 = π / (2L)
T_exact(x, t) = sin(λ0 * x) * exp(-κ * λ0^2 * t)
TL = 0.0
tspan = (0.0, 0.5)
solB, mesh = solve_bem_heat(; L = L, κ = κ, TL = TL, tspan = tspan, nq = 8,
    T0 = x -> sin(λ0 * x))
ts = range(tspan...; length = 80)
TR_bem = [solB(t)[1] for t in ts]
TR_ex = [T_exact(L, t) for t in ts]
p1 = plot(ts, TR_ex; lw = 2, label = "analítico", xlabel = "t", ylabel = "T(L,t)",
          title = "Ponta direita")
plot!(p1, ts, TR_bem; lw = 2, ls = :dash, label = "BEM + linhas (nq=$(mesh.ni))")
t_snap = 0.25
xp, Tp = profile_at(solB, mesh, TL, t_snap)
xx = range(0, L; length = 200)
p2 = plot(xx, T_exact.(xx, t_snap); lw = 2, label = "analítico",
          xlabel = "x", ylabel = "T", title = "Perfil t = $t_snap")
plot!(p2, xp, Tp; marker = :circle, lw = 2, ls = :dash, label = "BEM (nós)")
plot(p1, p2; layout = (1, 2), size = (900, 350))
saveplot("bem-calor-analitico.png")
err_TR = maximum(abs, TR_bem .- TR_ex)
err_prof = maximum(abs, Tp .- T_exact.(xp, t_snap))
@show err_TR err_prof
println("nq\terr_T(L,t)\terr_perfil(t=$t_snap)")
for nq in (2, 4, 6, 8, 12)
    soln, mn = solve_bem_heat(; L = L, κ = κ, TL = TL, tspan = tspan, nq = nq,
        T0 = x -> sin(λ0 * x))
    eL = maximum(abs, [soln(t)[1] - T_exact(L, t) for t in ts])
    xq, Tq = profile_at(soln, mn, TL, t_snap)
    ep = maximum(abs, Tq .- T_exact.(xq, t_snap))
    println("$nq\t$eL\t$ep")
end

println("=== 10. AB4 + Houbolt ===")
function ab4(ivp, n)
    a, b = ivp.tspan
    h = (b - a) / n
    t = [a + i * h for i in 0:n]
    u0 = ivp.u0 isa Number ? float(ivp.u0) : float.(ivp.u0)
    U = Vector{typeof(u0)}(undef, n + 1)
    _, Us = rk4(ODEProblem(ivp.f, ivp.u0, (a, a + 3h), ivp.p), 3)
    U[1:4] .= Us
    σ = [-9, 37, -59, 55] ./ 24
    F = [ivp.f(U[i], ivp.p, t[i]) for i in 1:4]
    for i in 4:n
        U[i+1] = U[i] + h * sum(σ[j] * F[j] for j in 1:4)
        F = (F[2], F[3], F[4], ivp.f(U[i+1], ivp.p, t[i+1]))
    end
    return t, U
end

function houbolt_oscillator(ω, u0, v0, tf, n)
    h = tf / n
    t = collect(range(0, tf; length = n + 1))
    u = zeros(n + 1)
    fsys = (y, p, τ) -> [y[2], -ω^2 * y[1]]
    _, Ys = rk4(ODEProblem(fsys, [u0, v0], (0.0, 2h), nothing), 2)
    u[1] = u0
    u[2] = Ys[2][1]
    u[3] = Ys[3][1]
    for i in 3:n
        rhs = (5u[i] - 4u[i-1] + u[i-2]) / h^2
        u[i+1] = rhs / (2 / h^2 + ω^2)
    end
    return t, u
end

ω = 2π
tH, uH = houbolt_oscillator(ω, 1.0, 0.0, 10.0, 400)
tR, UR = rk4(ODEProblem((y, p, τ) -> [y[2], -ω^2 * y[1]], [1.0, 0.0], (0.0, 10.0), nothing), 400)
plot(tH, uH; label = "Houbolt", lw = 2)
plot!(tR, [y[1] for y in UR]; label = "RK4", ls = :dash, lw = 2,
      xlabel = "t", ylabel = "u", title = "u'' + ω²u = 0")
saveplot("houbolt-vs-rk4.png")

tAB, UAB = ab4(ODEProblem(f, -1.0, (0.0, 4.0)), 200)
println("  AB4 last = ", UAB[end], "  ref = ", u_ref(4.0))

println("\nOK — figuras em ", OUT)
