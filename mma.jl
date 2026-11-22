non_discreteness(ρ̄) = sum(4 .* ρ̄ .* (1 .- ρ̄)) / length(ρ̄)

function simp_mma(prob, F, d_max; ρ₀, p, β, q, maxiter, inner_maxeval)
    ne = length(prob.element_dofs)
    ρ = impose_passive!(fill(ρ₀, ne), prob)
    lower, upper = zeros(ne), ones(ne)
    lower[prob.full_elements] .= 1
    upper[prob.void_elements] .= 0

    c_hist = Float64[]
    c_cases = [Float64[] for _ in axes(F, 2)]
    g_hist = [Float64[] for _ in axes(F, 2)]
    Mnd_hist = Float64[]
    V_hist = Float64[]
    stage_ends = Int[]
    ρ_hist = Vector{Float64}[]
    ρ̃_hist = Vector{Float64}[]
    ρ̄_hist = Vector{Float64}[]

    for (pₛ, βₛ) in zip(p, β)
        c_ref = first(compliance_cases!(Float64[], ρ, prob, pₛ, βₛ, F, q))
        opt = NLopt.Opt(:LD_MMA, ne)
        opt.params["inner_maxeval"] = inner_maxeval
        NLopt.lower_bounds!(opt, lower)
        NLopt.upper_bounds!(opt, upper)
        NLopt.maxeval!(opt, maxiter)
        NLopt.xtol_rel!(opt, 1e-8)

        NLopt.min_objective!(opt, (ρ, ∂c∂ρ) -> begin
            local ρ̄
            c, C, ρ̄ = compliance_cases!(∂c∂ρ, ρ, prob, pₛ, βₛ, F, q)
            push!(c_hist, c)
            for i in eachindex(C)
                push!(c_cases[i], C[i])
            end
            push!(Mnd_hist, non_discreteness(ρ̄))
            push!(V_hist, volume(ρ̄))
            ∂c∂ρ ./= c_ref
            return c / c_ref
        end)

        NLopt.inequality_constraint!(opt, (ρ, ∂V∂ρ) -> begin
            g = volume_constraint!(∂V∂ρ, ρ, prob, βₛ)
            ∂V∂ρ ./= prob.Vₜ
            return g / prob.Vₜ
        end, 1e-8 / prob.Vₜ)

        for (i, f) in enumerate(eachcol(F))
            NLopt.inequality_constraint!(opt, (ρ, ∂gᵢ∂ρ) -> begin
                gᵢ = nodal_displacement_constraint!(∂gᵢ∂ρ, ρ, prob, pₛ, βₛ, f, d_max) / d_max^2
                ∂gᵢ∂ρ ./= d_max^2
                push!(g_hist[i], gᵢ)
                return gᵢ
            end, 1e-8)
        end

        c, ρ, status = NLopt.optimize(opt, ρ)
        ρ̃, ρ̄ = physical_density(prob, ρ, βₛ)
        push!(stage_ends, length(c_hist))
        push!(ρ_hist, copy(ρ))
        push!(ρ̃_hist, ρ̃)
        push!(ρ̄_hist, ρ̄)
        println("β = $βₛ, p = $pₛ, c = $(c*c_ref), V = $(volume(ρ̄)), ",
                "status = $status, evals = $(NLopt.numevals(opt))")
    end

    return (ρ=ρ_hist, ρ̃=ρ̃_hist, ρ̄=ρ̄_hist, c=c_hist, c_cases=stack(c_cases),
            Mnd=Mnd_hist, V=V_hist, g_d=stack(g_hist), stage_ends=stage_ends)
end
