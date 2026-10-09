module WaterLilyNarrowBand

using WaterLily, Adapt
import WaterLily: AbstractBody, SetBody, measure, measure_sdf!, sdf, setmap, update!, @loop, inside, loc, δ

"""
    NarrowBand(body::AbstractBody, dims; k=4, mem=Array)

Wrap `body` to only measure it in an active region: the cells with `d²<fastd²` at the last
`measure_sdf!`, dilated by `k` cells. Elsewhere `d=±√fastd²`, with the sign from a warm-started
flood fill, and `measure` returns immediately.

The body must move less than a cell per measurement (the CFL condition). The active region
contains moves of up to `(k-1)/√3` cells in any direction, so the default `k=4` gives a safety
factor of 1.7. The first measurement, and any that moves beyond this, measure the whole domain
and flood fill from scratch, seeded by `outside!`.

  - `dims`: WaterLily domain size used to create the `Bool` cache arrays.
  - `k`: number of cells to dilate the active region.
  - `mem`: memory location of the cache. `Array`, `CuArray`, etc.
"""
struct NarrowBand{B<:AbstractBody,A<:AbstractArray{Bool}} <: AbstractBody
    body::B
    near::A
    reached::A
    farinside::A
    active::A
    k::Int
end
function NarrowBand(body::AbstractBody,dims;k=4,mem=Array)
    near,reached,farinside,active = ntuple(_->fill!(mem{Bool}(undef,dims.+2),true),4) # all active ⇒ full measure
    NarrowBand(body,near,reached,farinside,active,k)
end
Adapt.@adapt_structure NarrowBand
rewrap(a::NarrowBand,body) = NarrowBand(body,a.near,a.reached,a.farinside,a.active,a.k)
setmap(a::NarrowBand;kwargs...) = rewrap(a,setmap(a.body;kwargs...))
update!(a::NarrowBand,args...) = rewrap(a,update!(a.body,args...))

# measure d,n,V: skip inactive cells (x is assumed to be loc(i,I) for some i,I)
@inline function measure(a::NarrowBand,x::AbstractVector{T},t;fastd²=Inf) where T
    I = CartesianIndex(ntuple(i->floor(Int,x[i]+2),length(x)))
    if checkbounds(Bool,a.active,I) && !a.active[I]
        D = max(one(T),sqrt(T(fastd²)))
        return (a.reached[I] ? D : -D,zero(x),zero(x))
    end
    measure(a.body,x,t;fastd²)
end

"""
    measure_sdf!(d::AbstractArray, a::NarrowBand, t=0; fastd²=1)

Fill `d` with the signed distance from `a.body`, exact within `d² < max(fastd²,1)` and `±√fastd²` outside.
"""
function measure_sdf!(d::AbstractArray{T},a::NarrowBand,t=zero(T);fastd²=one(T)) where T
    fd² = max(T(fastd²),one(T))
    cold = all(a.active) || !narrow!(d,a,t,fd²)
    cold && full!(d,a,t,fd²)
    signs!(d,a,t,fd²,cold); activate!(a,d,fd²)
    return d
end

# measure the active region, false if the body moved out of it
function narrow!(d,a,t,fd²)
    (;body,active,farinside) = a
    @inside d[I] = active_sdf(body,active,I,t,fd²)
    touching!(farinside,d,active,fd²) && return false
    @inside farinside[I] = d[I]^2<fd²
    return any(farinside) # false if the body left the active region entirely
end
@inline active_sdf(body,active,I,t,fd²::T) where T = active[I] ? sdf(body,loc(0,I,T),t;fastd²=fd²) : √fd²

# any exact cell next to an inactive one?
touching!(s,d,active,fd²) = (fill!(s,false); @inside s[I] = d[I]^2<fd² && touches_inactive(active,I); any(s))
@inline function touches_inactive(active,I::CartesianIndex{N}) where N # ignoring ghosts
    for i in 1:N, s in (-1,1)
        J = I+s*δ(i,I)
        2 ≤ J[i] ≤ size(active,i)-1 && !active[J] && return true
    end; return false
end

# measure everywhere
full!(d,a,t,fd²) = (@debug "NarrowBand: full measurement"; fill!(a.active,true); narrow!(d,a,t,fd²))

# clear the flood fill seeds (d≥1) that `body` can't vouch for. Extend this for bodies without a reliable sign
# (e.g. MeshBody). Clearing is always safe: the flood fill reaches the outside cells from the band.
outside!(reached,body,t) = nothing
outside!(reached,body::SetBody,t) = (outside!(reached,body.a,t); outside!(reached,body.b,t))

# fill the signs: seed the flood fill where `body` vouches that d≥1 is outside (everywhere when cold, else
# only where d is exact), flood fill, then the reached cells are outside and the far-inside cells are inside
function signs!(d,a,t,fd²,cold)
    (;body,near,reached,farinside) = a
    @inside farinside[I] = d[I]≥1
    outside!(farinside,body,t)
    cold && fill!(reached,false)
    @inside reached[I] = cold || d[I]^2<fd² ? farinside[I] : reached[I]
    flood_fill!(near,reached,farinside,d)
    @inside d[I] = farinside[I] ? -abs(d[I]) : reached[I] ? abs(d[I]) : d[I]
end

# active region = exact cells dilated k times
function activate!(a,d,fd²)
    (;active,near) = a; fill!(active,false)
    @inside active[I] = d[I]^2<fd²
    for _ in 1:a.k
        @inside near[I] = flood_update(I,active)
        @inside active[I] = near[I]
    end
end

# Flood-fill to classify points as inside or outside a closed boundary
function flood_fill!(near, reached, scratch, sdf; cutoff=1f0)
    @. near = abs(sdf) < cutoff
    copyto!(scratch,reached)
    for _ in 1:min(size(near)...)÷2
        @inside scratch[I] = flood_update(I, reached, near)
        @inside reached[I] = flood_update(I, scratch, near)
        reached == scratch && break # converged
    end
    @. scratch = !near && !reached # far-inside
end
@inline function flood_update(I::CartesianIndex{d}, r, blocked=nothing) where d
    !isnothing(blocked) && blocked[I] && return false
    f = r[I]
    for i in 1:d
        f = f || r[I+δ(i,I)] || r[I-δ(i,I)]
    end; return f
end

export NarrowBand

end # module
