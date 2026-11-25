using LinearAlgebra, SparseArrays, StaticArrays, NLopt, JLD2, CairoMakie

include("mesh.jl")
include("fem.jl")
include("objectives.jl")
include("mma.jl")
include("plot.jl")

function main(θ_deg=[-90.0]; open=true, filename=nothing)
    # mesh
    lx, ly = 0.3, 0.1                        # [m]
    nx, ny = 180, 60
    # properties
    E, ν, t = 2e9, 0.3, 0.01
    # filter and soft kill
    r, Eₘᵢₙ = 0.02, 1e-4                     # [m], [-]
    # Force and constrains
    F_total, Vₜ, d_max = 218.0, 0.47, 0.001  # [N], [-], [m]
    # SIMP and q-norm
    ρ₀, q = 0.6, 8
    p = [2., 2., 3., 4., 4.]
    β = [1., 4., 4., 6., 16.]
    # Opt parameters
    maxiter, inner_maxeval = 100, 16
    # Setup Constant Components
    mesh = ring_mesh(lx, ly, nx, ny)
    (; ndofs, element_dofs, free_dofs, dofs_x, dofs_y, void_elements, full_elements) = mesh
    Kₑ   = plane_stress_Ke(ν)                  # unit E and t
    H    = HelmholtzFilter(mesh, r*nx/lx)
    pat  = AssemblyPattern(element_dofs, ndofs, free_dofs)
    prob = (; Kₑ, H, pat, element_dofs, free_dofs, dofs_x, dofs_y,
            void_elements, full_elements, Vₜ, Eₘᵢₙ)
    F = zeros(ndofs, length(θ_deg))
    F_n = F_total / length(dofs_y)
    for (i, θ) in enumerate(θ_deg)
        F[dofs_x, i] .= F_n*cosd(θ)
        F[dofs_y, i] .= F_n*sind(θ)
    end
    # Adimentional data
    F ./= F_total
    u_scale = F_total/(E*t)
    # Optimization
    simp_result = simp_mma(prob, F, d_max*(E*t)/F_total; ρ₀, p, β, q, maxiter, inner_maxeval)
    # Post processing data
    Û = solve_state(prob, last(simp_result.ρ̄), last(p), F)
    u = u_scale .* Û                         # [m]
    c_cases = reshape([dot(f, û) for (f, û) in zip(eachcol(F), eachcol(Û))], 1, :)
    c_cases .*= F_total^2/(E*t)              # [N·m]
    c = norm(vec(c_cases), q)
    r_nodes = hypot.(u[dofs_x, :], u[dofs_y, :])
    d_cases = [norm(rᵢ, 20) for rᵢ in eachcol(r_nodes)]
    max_displacement = vec(maximum(r_nodes; dims=1))
    # Save data
    if isnothing(filename)
        angle_suffix = join((replace(string(θ), r"\.0$" => "") for θ in θ_deg), "_")
        filename = "topopt_$(angle_suffix)_.jld2"
    end
    datos = joinpath(@__DIR__, filename)
    jldsave(datos; simp_result, θ_deg, p, β, q, ρ₀, nx, ny, lx, ly, r, E, ν, t,
            Eₘᵢₙ, F_total, Vₜ, d_max, u, c, c_cases, d_cases, max_displacement)
    # Plot
    # outdir = joinpath(@__DIR__, "output", "pdf")
    # pdf = plot_results(outdir, simp_result, mesh, nx, ny, p, β, Vₜ, d_max;
    #                    compliance_scale=F_total^2/(E*t))
    # open && run(`xdg-open $pdf`; wait=false)
    return nothing
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end

function path_2links(; npoints=8 )
    l = 0.3 * 4 / 6  # distance between ring centers [m]
    x_pickup, y_floor = -0.20, -0.10
    x_place, y_surface = -0.33,  0.15

    x_samples = range(x_pickup, x_place; length=npoints)
    t_samples = range(0.0, 1.0; length=npoints)
    y_fractions = log1p.(9 .* t_samples) ./ log1p(9)
    targets = [(x, y_floor + (y_surface - y_floor) * fraction)
               for (x, fraction) in zip(x_samples, y_fractions)]

    residual_tol = 1e-10  # [m]
    opt = NLopt.Opt(:LD_SLSQP, 2)
    NLopt.lower_bounds!(opt, [-π, -π])
    NLopt.upper_bounds!(opt, [ π,  π])
    NLopt.maxeval!(opt, 100)
    NLopt.stopval!(opt, residual_tol^2 / 2)

    q = [π, π/2]
    q₁_deg = Vector{Float64}(undef, length(targets))
    q₂_deg = similar(q₁_deg)
    θ₁_deg = similar(q₁_deg)
    θ₂_deg = similar(q₁_deg)
    residuals = similar(q₁_deg)
    for (i, (x, y)) in enumerate(targets)
        NLopt.min_objective!(opt, (angles, grad) -> begin
            q₁ , q₂ = angles
            s₁ , c₁ = sincos(q₁)
            s₁₂, c₁₂ = sincos(q₁ + q₂)
            rₓ, rᵧ = l*(c₁+c₁₂)-x, l*(s₁+s₁₂)-y
            if !isempty(grad)
                grad[1] = -l*(s₁+s₁₂)*rₓ + l*(c₁+c₁₂)*rᵧ
                grad[2] = l * (-s₁₂ * rₓ + c₁₂ * rᵧ)
            end
            return (rₓ^2 + rᵧ^2)/2
        end)

        objective, q, _ = NLopt.optimize(opt, q)
        q₁_deg[i] = rad2deg(q[1])
        q₂_deg[i] = rad2deg(q[2])
        # Each FEM +x axis points from the loaded ring to the supported ring.
        θ₁_deg[i] = 90.0 - q₁_deg[i]
        θ₂_deg[i] = 90.0 - (q₁_deg[i] + q₂_deg[i])
        residuals[i] = sqrt(2 * objective)
    end
    display(θ₁_deg)
    display(θ₂_deg)
    return (; l, targets, y_floor, q₁_deg, q₂_deg, θ₁_deg, θ₂_deg, θ_deg=θ₂_deg,
            residuals, residual_tol)
end

function generate_all_data()
    (; θ₁_deg, θ₂_deg) = path_2links(npoints=16)
    # independent
    for θ in vcat(θ₁_deg, θ₂_deg)
        main([θ])
    end
    # n points same time
    for npoints in 2:16
        (; θ₁_deg, θ₂_deg) = path_2links(; npoints)
        main(θ₁_deg; filename="topopt_eslabon1_n$(npoints).jld2")
        main(θ₂_deg; filename="topopt_eslabon2_n$(npoints).jld2")
    end
    return nothing
end
