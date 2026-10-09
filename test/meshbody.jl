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

# a NarrowBand of a SetBody: compare to the union of its parts, each measured on its own
@testset "diagonal SetBody with closed MeshBodies (3D)" begin
    T,dims,fd² = Float32,(56,40,40),6.25f0
    file = joinpath(pkgdir(WaterLilyMeshBodies),"test","meshes","sphere.stl")
    x₀ = SA{T}[16,16,16]; Δ = SA{T}[22,0,0] # second part offset
    sphere(Δ) = (b = MeshBody(file;scale=16f0,map=RigidMap(x₀,zero(x₀)),boundary=true); update!(b,(tri->tri.+Δ).(b.mesh)))
    ball() = AutoBody((x,t)->√sum(abs2,x-Δ)-6,RigidMap(x₀,zero(x₀)))
    σ(body) = (d = zeros(T,dims.+2); measure_sdf!(d,body,0f0;fastd²=fd²); d)
    for parts in (()->(sphere(zero(Δ)),sphere(Δ)), ()->(sphere(zero(Δ)),ball()))
        a,b = parts(); band = NarrowBand(+(parts()...),dims)
        for (n,Δx) in enumerate([[zero(x₀)]; diagonal(3)])
            a,b,band = setmap.((a,b,band);x₀=x₀+Δx)
            σu,σr = zeros(T,dims.+2),min.(σ(a),σ(b))
            @test fulls(()->measure_sdf!(σu,band,0f0;fastd²=fd²)) == (n==1) # only the first is full
            # same sign everywhere, same distance in the band (≈: an AutoBody's sdf isn't normalized, a SetBody's is)
            @test all(I->signbit(σu[I])==signbit(σr[I]) && (σr[I]^2≥fd² || σu[I]≈σr[I]),inside(σu))
        end
    end
end
