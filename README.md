# WaterLilyNarrowBand

A companion package to [WaterLily.jl](https://github.com/WaterLily-jl/WaterLily.jl) that speeds up measuring moving bodies. `NarrowBand(body, dims)` wraps any `AbstractBody` so it is only evaluated in a thin band of cells around its surface, giving identical BDIM fields and forces.

```julia
] add https://github.com/weymouth/WaterLilyNarrowBand.jl
```

### Usage

```julia
using WaterLily, WaterLilyNarrowBand, StaticArrays

L = 32; dims = (4L,2L,2L)
body = AutoBody((x,t)->√sum(abs2,x)-L/2, RigidMap(SA{Float32}[L,L,L],SA{Float32}[0,0,0])) # any AbstractBody
sim = Simulation(dims, (1,0,0), L; body=NarrowBand(body,dims)) # NarrowBand(body,dims;mem=CuArray) on a GPU

sim.body = setmap(sim.body; θ=SA[0,0,0.1f0]) # setmap and update! work through the wrapper
sim_step!(sim; remeasure=true)
```

Outside the band, the distance sign comes from the warm-started flood fill of [WaterLilyMeshBodies](https://github.com/WaterLily-jl/WaterLilyMeshBodies.jl), and `measure` returns immediately, so the force and moment integrals skip the body there too. The band is dilated by `k=4` cells, enough for the body to move `(k-1)/√3≈1.7` cells per step in any direction: a safety factor of 1.7 on the CFL condition. Larger moves fall back to measuring the whole domain, which is correct, just not faster.

For a 16-ring kirigami `SetBody` rotating on a $384×256×96$ grid (median seconds per call, `Float32`):

| | `measure!` | `sim_step!(remeasure=true)` |
|---|---|---|
| CPU (4 threads) | 2.74 → 0.18 | 4.03 → 1.52 |
| GPU (RTX A2000 Laptop) | 0.383 → 0.044 | 0.597 → 0.258 |

The BDIM fields match the unwrapped body exactly on both backends. The cache costs 4 bytes per cell.
