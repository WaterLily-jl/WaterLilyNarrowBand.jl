include("utils.jl")

# kirigami-like SetBody of wavy rings
function rings(map;R=12f0,n=4,thk=1.2f0)
    δR = R/n
    ring(i) = AutoBody(map) do (x,y,z),t
        r = hypot(y,z); δx = 2tanh(r/δR)*cos(4atan(z,y)+π*(i%2))
        hypot(x-δx,r-clamp(r,δR*(i-1)+thk,δR*i-thk))-thk
    end
    sum(ring,1:n)
end

@testset "rotating SetBody (3D)" begin
    T = Float32; x₀ = SA{T}[16,16,16]
    a,b = pair((40,32,32),()->rings(RigidMap(x₀,zero(x₀)));ϵ=0.5f0)
    for n in 1:20 # ½ cell steps at the tip radius ≈12
        a.body,b.body = setmap.((a.body,b.body);θ=SA{T}[0,0,0.04f0n])
        @test same(a,b,T(n);x₀) == (true,0)
    end
    @test count(b.body.active) < length(b.body.active)÷3
end

@testset "diagonal SetBody (3D)" begin
    T = Float32; x₀ = SA{T}[16,16,16]
    a,b = pair((48,40,40),()->rings(RigidMap(x₀,zero(x₀)));ϵ=0.5f0)
    for (n,Δx) in enumerate(diagonal(3))
        a.body,b.body = setmap.((a.body,b.body);x₀=x₀+Δx)
        @test same(a,b,T(n);x₀=x₀+Δx) == (true,0)
    end
    a.body,b.body = setmap.((a.body,b.body);x₀=x₀.+6) # beyond the safety factor
    @test same(a,b,13f0;x₀=x₀.+6) == (true,1)
end

@testset "diagonal AutoBody (2D)" begin
    T = Float32; x₀ = SA{T}[16,16]
    a,b = pair((48,40),()->AutoBody((x,t)->√sum(abs2,x)-6,RigidMap(x₀,0f0)))
    for (n,Δx) in enumerate(diagonal(2))
        a.body,b.body = setmap.((a.body,b.body);x₀=x₀+Δx)
        @test same(a,b,T(n)) == (true,0)
    end
end
