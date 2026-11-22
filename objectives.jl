function compliance_cases!(∂J∂ρ, ρ, prob, p, β, F, q)
    ρ̃, ρ̄ = physical_density(prob, ρ, β)
    U = solve_state(prob, ρ̄, p, F)
    C = [dot(f, u) for (f, u) in zip(eachcol(F), eachcol(U))]
    # q-norm of compliance, big q to focus on worst case.
    J = norm(C, q)
    if !isempty(∂J∂ρ)
        # max(...,floatmin) to avoid /0
        ∂J∂C = @. (C / max(J, floatmin(typeof(J))))^(q - 1)
        ∂J∂ρ̄ = zeros(length(ρ̄))
        for (e, dofs) in enumerate(prob.element_dofs)
            # Now each compliance ad efect on sensivity
            for (i, u) in enumerate(eachcol(U))
                uₑ = u[dofs]
                ∂J∂ρ̄[e] += ∂J∂C[i] * dot(uₑ, prob.Kₑ * uₑ)
            end
            ∂J∂ρ̄[e] *= -∂E_SIMP∂ρ̄(ρ̄[e], p, prob.Eₘᵢₙ)
        end
        # To mantain voids and fill elements before and after
        impose_passive!(∂J∂ρ̄, prob; full_value=0.0)
        # _/∂ρ̄ * ∂ρ̄∂ρ̃ * ∂ρ̃/∂ρ = _/∂ρ
        ∂J∂ρ .= prob.H(∂ρ̄∂ρ̃(ρ̃, β) .* ∂J∂ρ̄)
        impose_passive!(∂J∂ρ, prob; full_value=0.0)
    end
    return J, C, ρ̄
end

volume(ρ̄) = sum(ρ̄) / length(ρ̄)

function volume_constraint!(∂V∂ρ, ρ, prob, β)
    ρ̃, ρ̄ = physical_density(prob, ρ, β)
    if !isempty(∂V∂ρ)
        ∂V∂ρ̃ = ∂ρ̄∂ρ̃(ρ̃, β) ./ length(ρ̄)
        impose_passive!(∂V∂ρ̃, prob; full_value=0.0)
        ∂V∂ρ .= prob.H(∂V∂ρ̃)
        impose_passive!(∂V∂ρ, prob; full_value=0.0)
    end
    return volume(ρ̄) - prob.Vₜ
end

function nodal_displacement_constraint!(∂g∂ρ, ρ, prob, p, β, f, d_max)
    ρ̃, ρ̄ = physical_density(prob, ρ, β)
    K = factor_stiffness(prob, ρ̄, p)
    u = solve_state(prob, ρ̄, p, f; K)
    (; dofs_x, dofs_y) = prob
    uₓ, uᵧ = u[dofs_x], u[dofs_y]
    # Constrain euclidian disntance ^2
    s = @. uₓ^2 + uᵧ^2
    D² = norm(s, 10)
    if !isempty(∂g∂ρ)
        # max(...,floatmin) to avoid /0
        ∂g∂s = @. (s / max(D², floatmin(typeof(D²))))^9
        G = zeros(length(u))
        @. G[dofs_x] = 2 * uₓ * ∂g∂s
        @. G[dofs_y] = 2 * uᵧ * ∂g∂s
        λ = solve_state(prob, ρ̄, p, G; K)
        ∂g∂ρ̄ = zeros(length(ρ̄))
        for (e, dofs) in enumerate(prob.element_dofs)
            uₑ, λₑ = u[dofs], λ[dofs]
            ∂g∂ρ̄[e] = -∂E_SIMP∂ρ̄(ρ̄[e], p, prob.Eₘᵢₙ) * dot(λₑ, prob.Kₑ * uₑ)
        end
        impose_passive!(∂g∂ρ̄, prob; full_value=0.0)
        # _/∂ρ̄ * ∂ρ̄∂ρ̃ * ∂ρ̃/∂ρ = _/∂ρ
        ∂g∂ρ .= prob.H(∂ρ̄∂ρ̃(ρ̃, β) .* ∂g∂ρ̄)
        impose_passive!(∂g∂ρ, prob; full_value=0.0)
    end
    return D² - d_max^2
end
