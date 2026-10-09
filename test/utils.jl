using WaterLilyNarrowBand, WaterLily, StaticArrays, Test, Logging
using WaterLily: pressure_force, viscous_force, pressure_moment

# simulations of the plain and wrapped body (separate bodies: a MeshBody is mutated by update!)
function pair(dims,body::Function;T=Float32,ϵ=1f0)
    sim(body) = Simulation(dims,ntuple(i->i==1 ? 1 : 0,length(dims)),16;body,T,ϵ)
    a,b = sim(body()),sim(NarrowBand(body(),dims))
    b.flow.u .= a.flow.u .= rand(T,size(a.flow.u))
    return a,b
end

# count the full measurements logged by f()
struct FullCounter <: AbstractLogger; n::Base.RefValue{Int}; end
Logging.min_enabled_level(::FullCounter) = Logging.Debug
Logging.shouldlog(::FullCounter,args...) = true
Logging.catch_exceptions(::FullCounter) = false
Logging.handle_message(c::FullCounter,level,msg,args...;kw...) = (msg=="NarrowBand: full measurement" && (c.n[] += 1); nothing)
fulls(f) = (c = FullCounter(Ref(0)); with_logger(f,c); c.n[])

# measure both: are all BDIM fields and integrated loads identical, and how many full measurements?
function same(a,b,t;x₀=nothing)
    measure!(a,t); n = fulls(()->measure!(b,t))
    b.flow.p .= a.flow.p .= rand(eltype(a.flow.p),size(a.flow.p))
    a.flow.μ₀==b.flow.μ₀ && a.flow.μ₁==b.flow.μ₁ && a.flow.V==b.flow.V &&
        pressure_force(a)==pressure_force(b) && viscous_force(a)==viscous_force(b) &&
        (isnothing(x₀) || pressure_moment(x₀,a)==pressure_moment(x₀,b)), n
end

# diagonal displacements: realistic ½ cell steps, then 1 and 1.5 cells (the safety factor)
diagonal(D) = cumsum([fill(0.5f0,10); 1f0; 1.5f0]) .* Ref(fill(1/√Float32(D),SVector{D}))
