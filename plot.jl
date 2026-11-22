function plot_results(pdf, result, nx, ny, p, β, Vₜ, d_max, θ_deg)
    nrows = cld(length(p), 2)
    fig = Figure(size=(960, 185*nrows + 470), fontsize=14)

    for (i, ρ̄) in enumerate(result.ρ̄)
        row, col = div(i-1, 2)+1, mod(i-1, 2)+1
        Mᵈ = round(100 * (1 - non_discreteness(ρ̄)); digits=1)
        ax = Axis(fig[row, col]; title="p = $(p[i]), β = $(β[i])  |  Mᵈ = $Mᵈ %", aspect=DataAspect())
        hidedecorations!(ax)
        heatmap!(ax, 1:nx, 1:ny, reshape(ρ̄, nx, ny);
                 colormap=Reverse(:grays), colorrange=(0, 1), interpolate=false)
    end

    ρ_scale = range(0.0, 1.0; length=64)
    ax_scale = Axis(fig[1:nrows, 3]; width=16, yaxisposition=:right,
                    ylabel="Densidad ρ̄", yticks=[0.0, 0.5, 1.0])
    hidexdecorations!(ax_scale)
    heatmap!(ax_scale, [0.0, 1.0], ρ_scale, repeat(reshape(collect(ρ_scale), 1, :), 2, 1);
             colormap=Reverse(:grays), colorrange=(0, 1), interpolate=false)
    xlims!(ax_scale, 0, 1)
    ylims!(ax_scale, 0, 1)

    row = nrows + 1
    ax_c = Axis(fig[row, 1]; xlabel="Evaluación del objetivo", ylabel="Norma de compliance ‖C‖₈")
    ax_v = Axis(fig[row, 2]; xlabel="Evaluación del objetivo", ylabel="Fracción de volumen")
    lines!(ax_c, result.c; color=:royalblue, linewidth=2)
    lines!(ax_v, result.V; color=:seagreen, linewidth=2)
    hlines!(ax_v, [Vₜ]; color=:gray50, linestyle=:dash)
    for ax in (ax_c, ax_v)
        vlines!(ax, result.stage_ends[1:end-1] .+ 0.5; color=:gray65, linestyle=:dash)
    end

    limit = 1000*d_max
    D = limit .* sqrt.(max.(0, 1 .+ result.g_d))
    ax_d = Axis(fig[row+1, 1:2]; xlabel="Evaluación de la restricción", ylabel="D (mm)",
                title="Desplazamiento del anillo (límite = $limit mm)")
    for (Dᵢ, θ) in zip(eachcol(D), θ_deg)
        lines!(ax_d, Dᵢ; linewidth=2, label="$(θ)°")
    end
    Legend(fig[row+1, 3], ax_d)
    colsize!(fig.layout, 3, Fixed(90))
    D_min, D_max = extrema(D)
    margin = 0.1 * (D_max - D_min)
    if D_min - margin <= limit <= D_max + margin
        hlines!(ax_d, [limit]; color=:gray50, linestyle=:dash)
    end

    save(pdf, fig)
end
