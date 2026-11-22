using LinearAlgebra, SparseArrays, StaticArrays, NLopt, JLD2, CairoMakie

include("mesh.jl")
include("fem.jl")
include("objectives.jl")
include("mma.jl")
include("plot.jl")

function main(θ_deg=[60.0]; abrir=true)
    lx, ly = 0.3, 0.1                        # [m]
    nx, ny = 180, 60
    E, ν, t = 2e9, 0.3, 0.01
    r, Eₘᵢₙ = 0.02, 1e-4
    F_total, Vₜ, d_max = 218.0, 0.47, 0.001  # [N], [-], [m]
    ρ₀, q = 0.6, 8
    p = [2., 2., 3., 4., 4.]
    β = [1., 4., 4., 6., 16.]
    maxiter, inner_maxeval = 100, 16

    mesh = ring_mesh(lx, ly, nx, ny)
    (; ndofs, element_dofs, free_dofs, dofs_x, dofs_y, void_elements, full_elements) = mesh
    Kₑ = plane_stress_Ke(ν)                  # unit E and t
    H = HelmholtzFilter(mesh, r*nx/lx)
    pat = AssemblyPattern(element_dofs, ndofs, free_dofs)
    prob = (; Kₑ, H, pat, element_dofs, free_dofs, dofs_x, dofs_y,
            void_elements, full_elements, Vₜ, Eₘᵢₙ)

    F = zeros(ndofs, length(θ_deg))
    F_n = F_total / length(dofs_y)
    for (i, θ) in enumerate(θ_deg)
        F[dofs_x, i] .= F_n * cosd(θ)
        F[dofs_y, i] .= F_n * sind(θ)
    end
    F ./= F_total
    u_scale = F_total/(E*t)
    simp_result = simp_mma(prob, F, d_max*(E*t)/F_total; ρ₀, p, β, q, maxiter, inner_maxeval)

    Û = solve_state(prob, last(simp_result.ρ̄), last(p), F)
    u = u_scale .* Û                         # [m]
    c_cases = reshape([dot(f, û) for (f, û) in zip(eachcol(F), eachcol(Û))], 1, :)
    c_cases .*= F_total^2/(E*t)              # [N·m]
    c = norm(vec(c_cases), q)
    r_nodes = hypot.(u[dofs_x, :], u[dofs_y, :])
    d_cases = [norm(rᵢ, 20) for rᵢ in eachcol(r_nodes)]
    max_displacement = vec(maximum(r_nodes; dims=1))

    datos = joinpath(@__DIR__, "topopt.jld2")
    pdf = joinpath(@__DIR__, "topopt.pdf")
    jldsave(datos; simp_result, θ_deg, p, β, q, ρ₀, nx, ny, lx, ly, r, E, ν, t,
            Eₘᵢₙ, F_total, Vₜ, d_max, u, c, c_cases, d_cases, max_displacement)
    plot_results(pdf, simp_result, nx, ny, p, β, Vₜ, d_max, θ_deg)
    abrir && run(`xdg-open $pdf`; wait=false)
    return simp_result
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
