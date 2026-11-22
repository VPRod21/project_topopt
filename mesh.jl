function ring_mesh(lx, ly, nx, ny)
    dx, dy = lx/nx, ly/ny
    ids = LinearIndices((nx+1, ny+1))
    nodes = [SVector(i*dx, j*dy) for j in 0:ny for i in 0:nx]
    elements = [SVector(ids[i,j], ids[i+1,j], ids[i+1,j+1], ids[i,j+1])
                for j in 1:ny for i in 1:nx]

    rᵢ², rₒ² = (ly/5)^2, (1.2ly/5)^2
    left_radius²(x) = (x[1]-lx/6)^2 + (x[2]-ly/2)^2
    right_radius²(x) = (x[1]-5lx/6)^2 + (x[2]-ly/2)^2
    loaded = findall(x -> rᵢ² <= left_radius²(x) <= rₒ², nodes)
    supported = findall(x -> rᵢ² <= right_radius²(x) <= rₒ², nodes)

    centers = [sum(nodes[n] for n in e)/4 for e in elements]
    void_elements = union(findall(x -> left_radius²(x) < rᵢ², centers),
                          findall(x -> right_radius²(x) < rᵢ², centers))
    full_elements = union(findall(x -> rᵢ² <= left_radius²(x) <= rₒ², centers),
                          findall(x -> rᵢ² <= right_radius²(x) <= rₒ², centers))

    element_dofs = [SVector{8,Int}(Tuple(d for n in e for d in (2n-1, 2n))) for e in elements]
    ndofs = 2length(nodes)
    free_dofs = setdiff(1:ndofs, [d for n in supported for d in (2n-1, 2n)])
    dofs_y = 2 .* loaded
    dofs_x = dofs_y .- 1
    return (; nodes, elements, element_dofs, ndofs, free_dofs, dofs_x, dofs_y,
            void_elements, full_elements)
end
