module WaterLilyNarrowBand

using WaterLily, Adapt
import WaterLily: AbstractBody, SetBody, measure, measure_sdf!, sdf, setmap, update!, @loop, inside, loc, δ

"""
    NarrowBand(body::AbstractBody, dims; k=4, interior=hasinterior(body), mem=Array)

Only measure `body` within `k` cells of its surface, using a warm-started flood fill for the sign elsewhere.

  - `dims`: WaterLily domain size.
  - `k`: dilation of the band. Covers moves of `(k-1)/√3` cells per measurement.
  - `interior`: set `false` to skip the flood fill for bodies thinner than the band.
  - `mem`: memory location of the cache. `Array`, `CuArray`, etc.
"""
struct NarrowBand{B<:AbstractBody,A<:AbstractArray{Bool},P} <: AbstractBody
    body::B
    near::A
    reached::A
    farinside::A
    active::A
    proj::P
    k::Int
    interior::Bool
end
function NarrowBand(body::AbstractBody,dims;k=4,interior=hasinterior(body),mem=Array)
    near,reached,farinside,active = ntuple(_->fill!(mem{Bool}(undef,dims.+2),true),4) # full first measure, all outside if !interior
    proj = map(n->fill!(mem{Bool}(undef,n+2),false),dims) # active cells projected on each axis, empty ⇒ full first measure
    NarrowBand(body,near,reached,farinside,active,proj,k,interior)
end
Adapt.@adapt_structure NarrowBand
rewrap(a::NarrowBand,body) = NarrowBand(body,a.near,a.reached,a.farinside,a.active,a.proj,a.k,a.interior)
setmap(a::NarrowBand;kwargs...) = rewrap(a,setmap(a.body;kwargs...))
update!(a::NarrowBand,args...) = rewrap(a,update!(a.body,args...))

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
    R = box(a)
    cold = isempty(R) || !narrow!(d,a,R,t,fd²)
    cold && (R = inside(d); full!(d,a,R,t,fd²))
    signs!(d,a,t,fd²,cold); activate!(a,d,R,fd²)
    return d
end

function narrow!(d,a,R,t,fd²)
    (;body,active,farinside) = a; D = √fd²
    @inside d[I] = D
    @loop d[I] = active_sdf(body,active,I,t,fd²) over I ∈ R
    touching!(farinside,d,active,R,fd²) && return false
    @loop farinside[I] = d[I]^2<fd² over I ∈ R
    return any(farinside)
end
@inline active_sdf(body,active,I,t,fd²::T) where T = active[I] ? sdf(body,loc(0,I,T),t;fastd²=fd²) : √fd²

touching!(s,d,active,R,fd²) = (fill!(s,false); @loop s[I] = d[I]^2<fd² && touches_inactive(active,I) over I ∈ R; any(s))
@inline function touches_inactive(active,I::CartesianIndex{N}) where N
    for i in 1:N, s in (-1,1)
        J = I+s*δ(i,I)
        2 ≤ J[i] ≤ size(active,i)-1 && !active[J] && return true
    end; return false
end

full!(d,a,R,t,fd²) = (@debug "NarrowBand: full measurement"; fill!(a.active,true); narrow!(d,a,R,t,fd²))
box(a) = CartesianIndices(map(span,a.proj))
@static isdefined(WaterLily,:bbox) && (WaterLily.bbox(a::NarrowBand,R) = box(a) ∩ R)
span(p) = (i = findall(Array(p)); isempty(i) ? (1:0) : first(i):last(i))
grow(R,k,d) = (δ = k*oneunit(first(R)); (first(R)-δ:last(R)+δ) ∩ inside(d))

# extend these for bodies without a reliable sign or interior (e.g. MeshBody)
outside!(reached,body,t) = nothing
outside!(reached,body::SetBody,t) = (outside!(reached,body.a,t); outside!(reached,body.b,t))
hasinterior(body) = true
hasinterior(body::SetBody) = hasinterior(body.a) || hasinterior(body.b)

function signs!(d,a,t,fd²,cold)
    a.interior || return
    (;body,near,reached,farinside) = a
    @inside farinside[I] = d[I]≥1
    outside!(farinside,body,t)
    cold && fill!(reached,false)
    @inside reached[I] = cold || d[I]^2<fd² ? farinside[I] : reached[I]
    flood_fill!(near,reached,farinside,d)
    @inside d[I] = farinside[I] ? -abs(d[I]) : reached[I] ? abs(d[I]) : d[I]
end

function activate!(a,d,R,fd²)
    (;active,near,proj) = a; fill!(active,false); foreach(p->fill!(p,false),proj)
    @loop active[I] = d[I]^2<fd² over I ∈ R
    R = grow(R,a.k,d)
    for _ in 1:a.k
        @loop near[I] = flood_update(I,active) over I ∈ R
        @loop set_active!(active,near,proj,I) over I ∈ R
    end
end
@inline set_active!(active,near,proj,I) = (active[I] = near[I]; near[I] && foreach((p,i)->(p[i] = true),proj,Tuple(I)); nothing)

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
