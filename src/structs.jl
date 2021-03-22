
#### Exceptions ####

"""
Exception indicating that the specified method is not implemented for an interface.
"""
struct UnimplementedError <: Exception
    msg::String
end

#### Internal Types and Methods ####

struct Indexer
    next_index::Dict{NodeID, Index}
    indices::Dict{NodeID, Dict{String, Index}}
end

function Indexer()::Indexer
    idxr = Indexer(Dict{NodeID, Index}(),
                   Dict{NodeID, Dict{String,Index}}())
    return idxr
end

function _retrieve_and_advance_index(idxr::Indexer, nid::NodeID)::Index
    if !haskey(idxr.next_index, nid)
        idxr.next_index[nid] = Index(zero(INDEX))
    end
    idx = idxr.next_index[nid]
    idxr.next_index[nid] = _increment(idx)
    return idx
end

function index(idxr::Indexer, nid::NodeID, name::String)::Index

    if !haskey(idxr.indices, nid)
        idxr.indices[nid] = Dict{String, Index}()
    end
    node_vars = idxr.indices[nid]

    if haskey(node_vars, name)

        idx = node_vars[name]

    else

        idx = _retrieve_and_advance_index(idxr, nid)
        node_vars[name] = idx

    end

    return idx
end

mutable struct VariableInfo
    name::String
    xhat_id::XhatID
end

function VariableInfo(name::String,
                      nid::NodeID
                      )::VariableInfo
    return VariableInfo(name, nid, 0.0)
end

mutable struct HatVariable
    value::Float64 # Current value of variable
    vars::Set{VariableID} # All nonhat variable ids that contribute to this variable
end

HatVariable()::HatVariable = HatVariable(0.0, Set{VariableID}())
HatVariable(val::Float64, vid::VariableID) = HatVariable(val, Set{VariableID}([vid]))

function value(a::HatVariable)::Float64
    return a.value
end

function set_value(a::HatVariable, v::Float64)::Nothing
    a.value = v
    return
end

function add_variable(a::HatVariable, vid::VariableID)
    push!(a.vars, vid)
    return
end

function variables(a::HatVariable)::Set{VariableID}
    return a.vars
end

mutable struct ProblemData
    obj::Float64
    sts::MOI.TerminationStatusCode
    time::Float64
end

ProblemData() = ProblemData(0.0, MOI.OPTIMIZE_NOT_CALLED, 0.0)

struct ScenarioInfo
    pid::Int
    prob::Float64
    branch_vars::Dict{VariableID, Float64}
    leaf_vars::Dict{VariableID, Float64}
    w_vars::Dict{VariableID, Float64}
    xhat_vars::Dict{VariableID, Float64}
    problem_data::ProblemData
end

function ScenarioInfo(pid::Int,
                      prob::Float64,
                      branch_ids::Set{VariableID},
                      leaf_ids::Set{VariableID}
                      )::ScenarioInfo

    branch_map = Dict{VariableID, Float64}()
    w_dict = Dict{VariableID, Float64}()
    x_dict = Dict{VariableID, Float64}()
    for vid in branch_ids
        branch_map[vid] = 0.0
        w_dict[vid] = 0.0
        x_dict[vid] = 0.0
    end

    leaf_map = Dict{VariableID, Float64}(vid => 0.0 for vid in leaf_ids)

    return ScenarioInfo(pid,
                        prob,
                        branch_map,
                        leaf_map,
                        w_dict,
                        x_dict,
                        ProblemData())
end

function create_ph_dicts(sinfo::ScenarioInfo)
    return (sinfo.w_vars, sinfo.xhat_vars)
end

function objective_value(sinfo::ScenarioInfo)::Float64
    return sinfo.problem_data.obj
end

function retrieve_variable_value(sinfo::ScenarioInfo, vid::VariableID)::Float64
    if haskey(sinfo.branch_vars, vid)
        vi = sinfo.branch_vars[vid]
    else
        vi = sinfo.leaf_vars[vid]
    end
    return vi
end

struct PHIterate
    xhat::Dict{XhatID,Float64}
    x::Dict{VariableID,Float64}
    w::Dict{VariableID,Float64}
end

struct PHIterateHistory
    iterates::Dict{Int,PHIterate}
end

function PHIterateHistory()
    return PHIterateHistory(Dict{Int,PHIterate}())
end

function _save_iterate(phih::PHIterateHistory,
                       iter::Int,
                       phi::PHIterate,
                       )::Nothing
    phih.iterates[iter] = phi
    return
end

struct PHResidual
    abs_res::Float64
    rel_res::Float64
    xhat_sq::Float64
    x_sq::Float64
end

struct PHResidualHistory
    residuals::Dict{Int,PHResidual}
end

function PHResidualHistory()::PHResidualHistory
    return PHResidualHistory(Dict{Int,PHResidual}())
end

function residual_vector(phrh::PHResidualHistory)::Vector{Float64}
    if length(phrh.residuals) > 0
        max_iter = maximum(keys(phrh.residuals))
        return [phrh.residuals[k].abs_res for k in sort!(collect(keys(phrh.residuals)))]
    else
        return Vector{Float64}()
    end
end

function relative_residual_vector(phrh::PHResidualHistory)::Vector{Float64}
    if length(phrh.residuals) > 0
        max_iter = maximum(keys(phrh.residuals))
        return [phrh.residuals[k].rel_res for k in sort!(collect(keys(phrh.residuals)))]
    else
        return Vector{Float64}()
    end
end

function residual_components(phrh::PHResidualHistory)::NTuple{2,Vector{Float64}}
    if length(phrh.residuals) > 0
        max_iter = maximum(keys(phrh.residuals))
        sorted_keys = sort!(collect(keys(phrh.residuals)))
        xhat_sq = [phrh.residuals[k].xhat_sq for k in sorted_keys]
        x_sq = [phrh.residuals[k].x_sq for k in sorted_keys]
        return (xhat_sq, x_sq)
    else
        return (Vector{Float64}(), Vector{Float64}())
    end
end

function _save_residual(phrh::PHResidualHistory, iter::Int, res::PHResidual)::Nothing
    # @assert(!(iter in keys(phrh.residuals)))
    phrh.residuals[iter] = res
    return
end

#### Primary PH Data Structure ####

struct PHData{R <: AbstractPenaltyParameter}
    r::R
    scenario_tree::ScenarioTree
    scenario_map::Dict{ScenarioID, ScenarioInfo}
    xhat::Dict{XhatID, HatVariable}
    variable_data::Dict{VariableID, VariableInfo}
    indexer::Indexer
    iterate_history::PHIterateHistory
    residual_history::PHResidualHistory
    time_info::TimerOutputs.TimerOutput
end

function PHData(r::AbstractPenaltyParameter,
                tree::ScenarioTree,
                scen_proc_map::Dict{Int, Set{ScenarioID}},
                var_map::Dict{ScenarioID, Dict{VariableID, String}},
                time_out::TimerOutputs.TimerOutput
                )::PHData

    var_data = Dict{VariableID,VariableInfo}()
    xhat_dict = Dict{XhatID, HatVariable}()
    idxr = Indexer()

    scenario_map = Dict{ScenarioID, ScenarioInfo}()
    for (pid, scenarios) in pairs(scen_proc_map)
        for scen in scenarios

            branch_ids = Set{VariableID}()
            leaf_ids = Set{VariableID}()

            for (vid, vname) in pairs(var_map[scen])

                vnode = node(tree, vid.scenario, vid.stage)

                if vnode == nothing
                    error("Unable to locate node for variable id $vid.")
                end

                idx = index(idxr, vnode.id, vname)
                xhid = XhatID(vnode.id, idx)
                var_info = VariableInfo(vname, xhid)
                var_data[vid] = var_info

                if is_leaf(vnode)

                    push!(leaf_ids, vid)

                else

                    push!(branch_ids, vid)

                    if !haskey(xhat_dict, xhid)
                        xhat_dict[xhid] = HatVariable()
                    end
                    add_variable(xhat_dict[xhid], vid)

                end
            end

            scenario_map[scen] = ScenarioInfo(pid,
                                              tree.prob_map[scen],
                                              branch_ids,
                                              leaf_ids,
                                              )
        end
    end

    return PHData(r,
                  tree,
                  scenario_map,
                  xhat_dict,
                  var_data,
                  idxr,
                  PHIterateHistory(),
                  PHResidualHistory(),
                  time_out,
                  )
end

function Base.show(io::IO, phd::PHData)
    print(io, "ProgressiveHedging.jl data structure")
    return
end

function convert_to_variable_ids(phd::PHData, xid::XhatID)
    return variables(phd.xhat[xid])
end

function convert_to_xhat_id(phd::PHData, vid::VariableID)::XhatID
    return phd.variable_data[vid].xhat_id
end

function is_leaf(phd::PHData, xhid::XhatID)::Bool
    return is_leaf(phd.scenario_tree, xhid.node)
end

function name(phd::PHData, xid::XhatID)
    return name(phd, first(convert_to_variable_ids(phd, xid)))
end

function ph_variables(phd::PHData)::Dict{XhatID,HatVariable}
    return phd.xhat
end

function probability(phd::PHData, scenario::ScenarioID)::Float64
    return phd.scenario_map[scenario].prob
end

function residuals(phd::PHData)::Vector{Float64}
    return residual_vector(phd.residual_history)
end

function residual_components(phd::PHData)::NTuple{2,Vector{Float64}}
    return residual_components(phd.residual_history)
end

function relative_residuals(phd::PHData)::Vector{Float64}
    return relative_residual_vector(phd.residual_history)
end

function stage_id(phd::PHData, xid::XhatID)::StageID
    return phd.scenario_tree.tree_map[xid.node].stage
end

function scenario_bundle(phd::PHData, xid::XhatID)::Set{ScenarioID}
    return scenario_bundle(phd.scenario_tree, xid.node)
end

function scenarios(phd::PHData)::Set{ScenarioID}
    return scenarios(phd.scenario_tree)
end
