"""
Struct for user supbroblem callbacks.

**Fields**

* `name::String` : User's name for the callback. Defaults to `string(h)`.
* `h::Function` : Callback function. See notes below for calling signature.
* `ext::Dict{Symbol,Any}` : Dictionary to store data between callback calls or needed parameters.

The callback function `h` must have the signature
    `h(ext::Dict{Symbol,Any}, sp::T, niter::Int, scenario_id::ScenarioID)  where T <: AbstractSubproblem`
where `ext` is the same dictionary given to the `Callback` constructor, `sp` is a concrete type of `AbstractSubproblem` (see `AbstractSubproblem`), `niter` is the current iteration and `scenario_id` is a scenario identifier (see `ScenarioID`).
"""
struct SubproblemCallback
    name::String
    h::Function
    ext::Dict{Symbol,Any}
end

function Base.show(io::IO, spcb::SubproblemCallback)
    print(io, spcb.name)
    return
end

"""
    SubproblemCallback(f::Function)

Creates a `SubproblemCallback` structure for function `f`.
"""
function SubproblemCallback(f::Function)
    return SubproblemCallback(string(f), f, Dict{Symbol,Any}())
end

"""
    SubproblemCallback(f::Function, ext::Dict{Symbol,Any})

Creates a `SubproblemCallback` structure for function `f` with the external data dictionary `ext`.
"""
function SubproblemCallback(f::Function, ext::Dict{Symbol,Any})
    return SubproblemCallback(string(f), f, ext)
end

"""
    spcb(f::Function)
    spcb(f::Function, ext::Dict{Symbol,Any})
    spcb(name::String, f::Function, ext::Dict{Symbol,Any})

Shorthand for `SubproblemCallback`.
"""
spcb(f::Function) = SubproblemCallback(f)
spcb(f::Function, ext::Dict{Symbol,Any}) = SubproblemCallback(f, ext)
spcb(name::String, f::Function, ext::Dict{Symbol,Any}) = SubproblemCallback(name, f, ext)
