function plane_stress_Ke(ν)
    kk = SVector(  1/2 - ν/6,  1/8 + ν/8,  -1/4 - ν/12, -1/8 + 3ν/8,
                 -1/4 + ν/12, -1/8 - ν/8,          ν/6,  1/8 - 3ν/8)
    KE_IDX = @SMatrix [
        1 2 3 4 5 6 7 8;
        2 1 8 7 6 5 4 3;
        3 8 1 6 7 4 5 2;
        4 7 6 1 8 3 2 5;
        5 6 7 8 1 2 3 4;
        6 5 4 3 2 1 8 7;
        7 4 5 2 3 8 1 6;
        8 3 2 5 4 7 6 1]
    return 1 / (1 - ν^2) * kk[KE_IDX]
end

struct AssemblyPattern
    K::SparseMatrixCSC{Float64,Int}
    slot::Matrix{Int}
end

function AssemblyPattern(element_dofs, ndofs, free)
    M = 8
    nf = length(free)
    free_index = zeros(Int, ndofs)
    free_index[free] = 1:nf
    rows, cols = Int[], Int[]
    for dofs in element_dofs, j in 1:M, i in 1:M
        ri, cj = free_index[dofs[i]], free_index[dofs[j]]
        if ri > 0 && cj > 0
            push!(rows, ri); push!(cols, cj)
        end
    end
    K = sparse(rows, cols, ones(length(rows)), nf, nf)

    slot = zeros(Int, M*M, length(element_dofs))
    for (e, dofs) in enumerate(element_dofs), j in 1:M, i in 1:M
        ri, cj = free_index[dofs[i]], free_index[dofs[j]]
        (ri == 0 || cj == 0) && continue
        col = K.colptr[cj]:K.colptr[cj+1]-1          # rows of column cj, sorted
        p = col[searchsortedfirst(view(K.rowval, col), ri)]
        slot[i + (j-1)*M, e] = p
    end
    return AssemblyPattern(K, slot)
end

function assemble!(pat::AssemblyPattern, Kₑ::StaticMatrix, scale::AbstractVector)
    nz = pat.K.nzval
    fill!(nz, 0.0)
    for e in axes(pat.slot, 2)
        for k in axes(pat.slot, 1)
            s = pat.slot[k, e]
            s == 0 && continue                        # fixed DOF
            nz[s] += scale[e] * Kₑ[k]
        end
    end
    return pat.K
end

struct HelmholtzFilter
    KF_chol::SparseArrays.CHOLMOD.Factor{Float64,Int}
    Tₑₙ::SparseMatrixCSC{Float64,Int}
end

function HelmholtzFilter(mesh, rmin)
    # rmin is the physical filter radius divided by the element size.
    R = rmin/2/sqrt(3)
    # Q4 element matrices for unit element: Laplacian (K_lap/6) and mass (M/36)
    K1 = [ 4 -1 -2 -1;
          -1  4 -1 -2;
          -2 -1  4 -1;
          -1 -2 -1  4]*R^2/6
    K2 = [ 4  2  1  2;
           2  4  2  1;
           1  2  4  2;
           2  1  2  4]/36
    KFe = K1 + K2
    # Assemble KF ∈ ℝ^{nn×nn} from the element matrices (COO → CSC)
    coo_length = length(mesh.elements)*4^2
    vals = Vector{Float64}(undef,coo_length)
    rows = Vector{Int}(undef,coo_length)
    cols = Vector{Int}(undef,coo_length)
    n=1
    for elmt in mesh.elements
        for i=1:4,j=1:4
            rows[n]=elmt[i]
            cols[n]=elmt[j]
            vals[n]=KFe[i,j]
            n+=1
        end
    end
    KF=sparse(rows, cols, vals, length(mesh.nodes), length(mesh.nodes))
    KF_chol = cholesky(Symmetric(KF))
    # Tₑₙ ∈ ℝ^{ne×nn}: each element row has four entries 0.25
    coo_length = length(mesh.elements)*4
    vals = 0.25*ones(coo_length)
    rows = Vector{Int}(undef,coo_length)
    cols = Vector{Int}(undef,coo_length)
    n=1
    for (e,elmt) in enumerate(mesh.elements)
        for node in elmt
            rows[n]=e
            cols[n]=node
            n+=1
        end
    end
    Tₑₙ=sparse(rows, cols, vals, length(mesh.elements),length(mesh.nodes))
    return HelmholtzFilter(KF_chol,Tₑₙ)
end

(H::HelmholtzFilter)(ρ) = H.Tₑₙ*(H.KF_chol\(H.Tₑₙ'*ρ))

function project_density(ρ̃, β)
    tβ = tanh(β/2)
    return @. (tβ + tanh(β*(ρ̃-0.5))) / (2tβ)
end

function ∂ρ̄∂ρ̃(ρ̃, β)
    tβ = tanh(β/2)
    return @. β*(1 - tanh(β*(ρ̃-0.5))^2) / (2tβ)
end

E_SIMP(ρ̄ₑ, p, Eₘᵢₙ) = Eₘᵢₙ + ρ̄ₑ^p * (1 - Eₘᵢₙ)

∂E_SIMP∂ρ̄(ρ̄ₑ, p, Eₘᵢₙ) = p * ρ̄ₑ^(p-1) * (1 - Eₘᵢₙ)

function factor_stiffness(prob, ρ̄, p)
    scale = E_SIMP.(ρ̄, p, prob.Eₘᵢₙ) # ρ^p with soft kill
    K_ff = assemble!(prob.pat, prob.Kₑ, scale)
    return cholesky(Symmetric(K_ff))
end

function solve_state(prob, ρ̄, p, F; K=factor_stiffness(prob, ρ̄, p))
    dofs = (prob.free_dofs, axes(F)[2:end]...)
    U = zeros(size(F))
    U[dofs...] = K \ F[dofs...]
    return U
end

function impose_passive!(ρ, prob; full_value=1.0)
    ρ[prob.void_elements] .= 0.0
    ρ[prob.full_elements] .= full_value
    return ρ
end

# Prescribed input densities also participate in the filter.
function physical_density(prob, ρ, β)
    ρ̃ = prob.H(impose_passive!(copy(ρ), prob))
    ρ̄ = impose_passive!(project_density(ρ̃, β), prob)
    return ρ̃, ρ̄
end
