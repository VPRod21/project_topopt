const COLUMN_WIDTH = 241    # 85 mm at 72 pt/in
const PAGE_WIDTH = 496      # 175 mm at 72 pt/in
const DENSITY_COLORS = Reverse(:grays)
const LIMIT_STYLE = (color=:gray35, linestyle=:dot, linewidth=1.05)
const CASE_COLORS = ("#0072B2", "#D55E00", "#009E73", "#8B4B9D",
                     "#916300", "#397B88", "#AF3764", "#4C573A")
const CASE_STYLES = (:solid, :dash, :dot, :dashdot, :dash,
                     :solid, :dashdot, :dot)
const CASE_LABELS = ("θ₁", "θ₂", "θ₃", "θ₄", "θ₅", "θ₆", "θ₇", "θ₈")

decimal(x) = replace(string(x), '.' => ',')
compact_decimal(x) = decimal(isinteger(x) ? Int(x) : x)
decimal_ticks(values) = compact_decimal.(round.(values; sigdigits=8))

function article_figure(width, height; fontsize=8, padding=5)
    return Figure(; size=(width, height), fontsize, figure_padding=padding,
                  fonts=(; regular="TeX Gyre Heros Makie"))
end

function history_axis(slot; xlabel, ylabel)
    return Axis(slot; xlabel, ylabel, xlabelpadding=3, ylabelpadding=3,
                xticklabelsize=7, yticklabelsize=7, ytickformat=decimal_ticks,
                xticksize=3, yticksize=3, xgridvisible=false,
                ygridcolor=(:black, 0.11), topspinevisible=false, rightspinevisible=false)
end

function stage_boundaries!(ax, stage_ends)
    vlines!(ax, stage_ends[1:end-1] .+ 0.5; color=:gray68, linestyle=:dot, linewidth=0.75)
end

function plot_cases!(ax, values)
    for (i, series) in enumerate(eachcol(values))
        lines!(ax, series; color=CASE_COLORS[i], linestyle=CASE_STYLES[i], linewidth=1.2,
               label=CASE_LABELS[i])
    end
end

function article_legend(slot, args...; nbanks=2)
    return Legend(slot, args...; nbanks, tellwidth=false, tellheight=true,
                  framevisible=false, labelsize=7, patchsize=(13, 7),
                  padding=(0, 0, 0, 0), rowgap=3, colgap=8, patchlabelgap=3)
end

function plot_mesh(mesh, nx, ny)
    fig = article_figure(COLUMN_WIDTH, 135; padding=(5, 10, 5, 5))
    points = [Point2f(1000x, 1000y) for (x, y) in mesh.nodes]
    xs = [pt[1] for pt in points[1:nx+1]]
    ys = [pt[2] for pt in points[1:nx+1:end]]
    ax = Axis(fig[1, 1]; aspect=DataAspect(), xlabel="x (mm)", ylabel="y (mm)",
              xticks=range(first(xs), last(xs); length=4), yticks=[first(ys), last(ys)],
              xticklabelsize=7, yticklabelsize=7, xlabelpadding=3, ylabelpadding=3,
              xticksize=3, yticksize=3, xgridvisible=false, ygridvisible=false,
              backgroundcolor=:white)
    solid_color, void_color = "#9C80BA", "#F3D48A"
    for (ids, color) in ((mesh.full_elements, solid_color), (mesh.void_elements, void_color))
        cells = [CairoMakie.Makie.GeometryBasics.Polygon(points[mesh.elements[i]]) for i in ids]
        poly!(ax, cells; color, strokewidth=0)
    end
    segments = Point2f[]
    for x in xs
        append!(segments, (Point2f(x, first(ys)), Point2f(x, last(ys))))
    end
    for y in ys
        append!(segments, (Point2f(first(xs), y), Point2f(last(xs), y)))
    end
    linesegments!(ax, segments; color=(:gray40, 0.4), linewidth=0.18)
    load_color, support_color = "#005989", "#9F3800"
    scatter!(ax, points[mesh.loaded]; color=load_color, markersize=0.8, strokewidth=0)
    scatter!(ax, points[mesh.supported]; color=support_color, markersize=0.8, strokewidth=0)
    limits!(ax, first(xs), last(xs), first(ys), last(ys))
    handles = [MarkerElement(color=load_color, marker=:circle, markersize=5),
               MarkerElement(color=support_color, marker=:circle, markersize=5),
               PolyElement(color=solid_color, strokecolor=:gray45, strokewidth=0.5),
               PolyElement(color=void_color, strokecolor=:gray45, strokewidth=0.5)]
    article_legend(fig[2, 1], handles, ["Nós com carga", "Nós fixos", "Elementos sólidos", "Elementos vazios"])
    rowgap!(fig.layout, 5)
    return fig
end

function plot_topology(result, nx, ny, p, β)
    nrows = cld(length(p), 2)
    shared_row = isodd(length(p))
    fig = article_figure(PAGE_WIDTH, 92*nrows + (shared_row ? 0 : 29); fontsize=9)
    for (i, ρ̄) in enumerate(result.ρ̄)
        row, col = divrem(i-1, 2) .+ 1
        title = "Etapa $i   ·   p = $(compact_decimal(p[i]))   ·   β = $(compact_decimal(β[i]))"
        ax = Axis(fig[row, col]; title, titlesize=9, titlefont=:regular, titlegap=4,
                  aspect=DataAspect(), spinewidth=0.55,
                  leftspinecolor=:gray48, rightspinecolor=:gray48,
                  bottomspinecolor=:gray48, topspinecolor=:gray48)
        hidedecorations!(ax)
        heatmap!(ax, 1:nx, 1:ny, reshape(ρ̄, nx, ny);
                 colormap=DENSITY_COLORS, colorrange=(0, 1), interpolate=false)
        limits!(ax, 0.5, nx+0.5, 0.5, ny+0.5)
    end
    scale_slot = shared_row ? fig[nrows, 2] : fig[nrows+1, 1:2]
    Colorbar(scale_slot; vertical=false, width=144, size=6,
             tellwidth=false, tellheight=!shared_row,
             limits=(0, 1), colormap=DENSITY_COLORS, ticks=([0, 0.5, 1], ["0", "0,5", "1"]),
             ticklabelsize=8, ticksize=3, spinewidth=0.5, label="Densidade física ρ̄", labelpadding=2)
    colgap!(fig.layout, 13)
    rowgap!(fig.layout, 6)
    return fig
end

function plot_compliance(result, compliance_scale)
    fig = article_figure(COLUMN_WIDTH, 177)
    ax = history_axis(fig[1, 1]; xlabel="Avaliação do objetivo", ylabel="Flexibilidade (N·m)")
    stage_boundaries!(ax, result.stage_ends)
    if size(result.c_cases, 2) == 1
        lines!(ax, compliance_scale .* result.c; color=CASE_COLORS[1], linewidth=1.35,
               label="q-norm · $(first(CASE_LABELS))")
    else
        lines!(ax, compliance_scale .* result.c; color=:gray15, linewidth=1.6, label="q-norm")
        plot_cases!(ax, compliance_scale .* result.c_cases)
    end
    article_legend(fig[2, 1], ax; nbanks=5)
    rowgap!(fig.layout, 5)
    return fig
end

function plot_discreteness(result)
    fig = article_figure(COLUMN_WIDTH, 147)
    ax = history_axis(fig[1, 1]; xlabel="Avaliação do objetivo", ylabel="Discretização (%)")
    stage_boundaries!(ax, result.stage_ends)
    lines!(ax, 100 .* (1 .- result.Mnd); color=CASE_COLORS[1], linewidth=1.35)
    ax.yticks = 0:25:100
    ylims!(ax, 0, 100)
    return fig
end

function plot_constraints(result, Vₜ, d_max)
    fig = article_figure(COLUMN_WIDTH, 261)
    ax_v = history_axis(fig[1, 1]; xlabel="Avaliação do objetivo", ylabel="Fração volumétrica")
    ax_d = history_axis(fig[2, 1]; xlabel="Avaliação da restrição", ylabel="Deslocamento (mm)")
    stage_boundaries!(ax_v, result.stage_ends)
    lines!(ax_v, result.V; color=:gray15, linewidth=1.35)
    hlines!(ax_v, [Vₜ]; LIMIT_STYLE...)
    limit = 1000d_max
    D = limit .* sqrt.(max.(0, 1 .+ result.g_d))
    plot_cases!(ax_d, D)
    hlines!(ax_d, [limit]; LIMIT_STYLE...)
    article_legend(fig[3, 1], ax_d; nbanks=4)
    rowgap!(fig.layout, 7)
    return fig
end

function plot_results(outdir, result, mesh, nx, ny, p, β, Vₜ, d_max; compliance_scale)
    mkpath(outdir)
    topology_pdf = joinpath(outdir, "topologia.pdf")
    save(topology_pdf, plot_topology(result, nx, ny, p, β); pt_per_unit=1)
    save(joinpath(outdir, "flexibilidade.pdf"), plot_compliance(result, compliance_scale); pt_per_unit=1)
    save(joinpath(outdir, "discretizacao.pdf"), plot_discreteness(result); pt_per_unit=1)
    save(joinpath(outdir, "restricoes.pdf"), plot_constraints(result, Vₜ, d_max); pt_per_unit=1)
    save(joinpath(outdir, "malha.pdf"), plot_mesh(mesh, nx, ny); pt_per_unit=1)
    return topology_pdf
end

function plot_path_2links(path=path_2links())
    (; l, targets, y_floor, q₁_deg, q₂_deg, θ₁_deg, θ₂_deg, residuals, residual_tol) = path
    poses = eachindex(targets)
    failed = .!(residuals .<= residual_tol)
    target_color, arm_color, failure_color = "#848B98", "#267383", "#B75531"

    fig = article_figure(PAGE_WIDTH, 225; padding=8)
    ax_path = history_axis(fig[1, 1]; xlabel="x (m)", ylabel="y (m)")
    ax_path.title = "(a) Configurações do braço"
    ax_path.titlesize = 8
    ax_path.aspect = DataAspect()
    ax_path.xtickformat = decimal_ticks
    ax_angles = history_axis(fig[1, 2]; xlabel="Configuração", ylabel="Ângulo (°)")
    ax_angles.title = "(b) Ângulos por configuração"
    ax_angles.titlesize = 8
    ax_angles.xticks = poses
    ax_angles.yticks = -200:100:200
    ylims!(ax_angles, -210, 210)

    hlines!(ax_path, [y_floor]; color=(:gray, 0.4), linewidth=0.8)
    lines!(ax_path, first.(targets), last.(targets);
           color=target_color, linestyle=:dash, linewidth=1.3)
    scatter!(ax_path, first.(targets), last.(targets); color=target_color, markersize=4.5)
    for i in poses
        q₁ = q₁_deg[i]
        ϕ = q₁ + q₂_deg[i]
        x_elbow, y_elbow = l * cosd(q₁), l * sind(q₁)
        x_tip, y_tip = x_elbow + l * cosd(ϕ), y_elbow + l * sind(ϕ)
        color = failed[i] ? failure_color : arm_color
        lines!(ax_path, [0, x_elbow, x_tip], [0, y_elbow, y_tip]; color=(color, 0.7), linewidth=1.2)
        scatter!(ax_path, [x_elbow, x_tip], [y_elbow, y_tip]; color, markersize=[3.75, 6])
        if failed[i]
            x_target, y_target = targets[i]
            lines!(ax_path, [x_tip, x_target], [y_tip, y_target];
                   color, linestyle=:dot, linewidth=1.1)
            vspan!(ax_angles, i - 0.4, i + 0.4; color=(color, 0.11))
        end
    end
    scatter!(ax_path, [0.0], [0.0]; color=:slategray, markersize=7.5)

    entries = [LineElement(color=target_color, linestyle=:dash), LineElement(color=arm_color)]
    labels = ["Trajetória-alvo", "Configuração convergida"]
    if any(failed)
        push!(entries, LineElement(color=failure_color))
        push!(labels, "Fora da tolerância")
    end
    article_legend(fig[2, 1], entries, labels; nbanks=2)

    for (values, color, linestyle, marker, label) in (
        (q₁_deg, arm_color, :solid, :circle, "q₁"),
        (q₂_deg, "#62539A", :solid, :rect, "q₂"),
        (θ₁_deg, "#0072B2", :dash, :utriangle, "θ₁"),
        (θ₂_deg, "#A45628", :dash, :diamond, "θ₂"),
    )
        scatterlines!(ax_angles, poses, values; color, linestyle, marker, label,
                      markersize=6, linewidth=1.2)
    end
    article_legend(fig[2, 2], ax_angles; nbanks=4)
    rowgap!(fig.layout, 5)
    return fig
end

function save_path_2links(path=path_2links(); outdir=joinpath(@__DIR__, "output"))
    fig = plot_path_2links(path)
    pdfdir = mkpath(joinpath(outdir, "pdf"))
    pdf = joinpath(pdfdir, "trajetoria_2r.pdf")
    save(pdf, fig; pt_per_unit=1)
    save(joinpath(outdir, "trayectoria_2r.png"), fig; px_per_unit=3)
    return pdf
end
