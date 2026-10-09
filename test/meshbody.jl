# Needs WaterLilyMeshBodies in the environment, so it is not part of runtests.jl yet.
include("utils.jl")
using WaterLilyMeshBodies

@testset "diagonal closed MeshBody (3D)" begin
    T,dims = Float32,(40,40,40)
    file = joinpath(pkgdir(WaterLilyMeshBodies),"test","meshes","sphere.stl")
    # a bare closed MeshBody flood fills from scratch every measurement: the exact reference
    a,b = pair(dims,()->MeshBody(file;scale=16f0,map=(x,t)->x.-16,boundary=true))
    for (n,Δx) in enumerate(diff([[zero(SVector{3,T})]; diagonal(3)]))
        shift(tri) = tri .+ Δx
        a.body = update!(a.body,shift.(a.body.mesh),1f0)
        b.body = update!(b.body,shift.(b.body.body.mesh),1f0)
        @test same(a,b,T(n);x₀=SA{T}[16,16,16]) == (true,0)
    end
end
