
const NODE_ID = Int
const SCENARIO_ID = Int
const STAGE_ID = Int
const INDEX = Int

struct NodeID
    value::NODE_ID
end
_value(nid::NodeID)::NODE_ID = nid.value
Base.isless(a::NodeID, b::NodeID) = _value(a) < _value(b)

struct ScenarioID
    value::SCENARIO_ID
end
_value(scid::ScenarioID)::SCENARIO_ID = scid.value
Base.isless(a::ScenarioID, b::ScenarioID) = _value(a) < _value(b)

const DUMMY_SCENARIO_ID = ScenarioID(-1)

struct StageID
    value::STAGE_ID
end
_value(sid::StageID)::STAGE_ID = sid.value
_increment(sid::StageID)::StageID = StageID(_value(sid) + one(STAGE_ID))
Base.isless(a::StageID, b::StageID) = _value(a) < _value(b)

const DUMMY_STAGE_ID = StageID(-1)

struct Index
    value::INDEX
end
_value(idx::Index)::INDEX = idx.value
_increment(index::Index)::Index = Index(_value(index) + one(INDEX))
Base.isless(a::Index, b::Index) = _value(a) < _value(b)

mutable struct Generator
    next_node_id::NODE_ID
    # next_stage_id::STAGE_ID
    next_scenario_id::SCENARIO_ID
end

function Generator()
    return Generator(0,0)
end

function _generate_node_id(gen::Generator)::NodeID
    id = NodeID(gen.next_node_id)
    gen.next_node_id += 1
    return id
end

function _generate_scenario_id(gen::Generator)::ScenarioID
    id = ScenarioID(gen.next_scenario_id)
    gen.next_scenario_id += 1
    return id
end

struct VariableID
    # scenario::ScenarioID # scenario to which this variable belongs
    stage::StageID # stage to which this variable belongs
    index::Index # coordinate in vector
end

function Base.isless(a::VariableID, b::VariableID)
    return (a.stage < b.stage ||
            (a.stage == b.stage && a.index < b.index))
end

struct UniqueVariableID
    scenario::ScenarioID # scenario to which this variable belongs
    stage::StageID # stage to which this variable belongs
    index::Index # coordinate in vector
end

function Base.isless(a::UniqueVariableID, b::UniqueVariableID)
    return (a.stage < b.stage ||
            (a.stage == b.stage &&
             (a.scenario < b.scenario ||
              (a.scenario == b.scenario && a.index < b.index))))
end

function UniqueVariableID(s::ScenarioID, vid::VariableID)
    return UniqueVariableID(s, vid.stage, vid.index)
end

struct XhatID
    node::NodeID
    index::Index
end

function Base.isless(a::XhatID, b::XhatID)
    return (a.node < b.node ||
            (a.node == b.node && a.index < b.index))
end

# struct Translator{A,B}
#     a_to_b::Dict{A,B}
#     b_to_a::Dict{B,A}
# end

# function Translator{A,B}() where {A,B}
#     return Translator(Dict{A,B}(), Dict{B,A}())
# end

# function add_pair(t::Translator{A,B}, a::A, b::B) where{A,B}
#     if !(a in keys(t.a_to_b)) && !(b in keys(t.b_to_a))
#         t.a_to_b[a] = b
#         t.b_to_a[b] = a
#     else
#         @error("One of $a or $b already maps to something")
#     end
#     return
# end

# function translate(t::Translator{A,B}, a::A) where {A,B}
#     return t.a_to_b[a]
# end

# function translate(t::Translator{A,B}, b::B) where {A,B}
#     return t.b_to_a[b]
# end

"""
Struct representing a node in a scenario tree.
"""
struct ScenarioNode
    id::NodeID # id of this node
    stage::StageID # stage of this node
    scenario_bundle::Set{ScenarioID} # scenarios that are indistiguishable
    variable_indices::Set{Index} # var indices
    parent::Union{Nothing, ScenarioNode}
    children::Set{ScenarioNode}
end

function _add_child(parent::ScenarioNode, child::ScenarioNode)
    push!(parent.children, child)
    return
end

function _create_node(gen::Generator, parent::Union{Nothing, ScenarioNode})
    nid = _generate_node_id(gen)
    stage = (parent==nothing ? StageID(1) : _increment(parent.stage))
    sn = ScenarioNode(nid, stage,
                      Set{ScenarioID}(), Set{Index}(),
                      parent, Set{ScenarioNode}())
    if parent != nothing
        _add_child(parent, sn)
    end
    return sn
end

function next_index(node::ScenarioNode)
    v_set = node.variable_indices
    idx = isempty(v_set) ? Index(1) : (_increment(maximum(v_set)))
    push!(node.variable_indices, idx)
    return idx
end

"""
Struct representing the scenario structure of a stochastic program.

Can be built up by the user using the functions `add_node` and `add_leaf`.

**Constructor**

ScenarioTree()

Default constructor generates the root node of the tree. Can get the root node with `root`.
"""
struct ScenarioTree
    root::ScenarioNode
    tree_map::Dict{NodeID, ScenarioNode} # map from NodeID to tree node
    stage_map::Dict{StageID, Set{NodeID}} # nodes in each stage
    prob_map::Dict{ScenarioID, Float64}
    id_gen::Generator
end

function ScenarioTree(root_node::ScenarioNode, gen::Generator)

    tree_map = Dict{NodeID, ScenarioNode}()
    stage_map = Dict{StageID, Set{NodeID}}()
    prob_map = Dict{ScenarioID, Float64}()

    st = ScenarioTree(root_node,
                      tree_map,
                      stage_map,
                      prob_map,
                      gen)

    _add_node(st, root_node)
    return st
end

function ScenarioTree()
    gen = Generator()
    rn = _create_node(gen, nothing)
    st = ScenarioTree(rn, gen)
    return st
end

"""
    root(tree::ScenarioTree)

Return the root node of the given ScenarioTree
"""
function root(tree::ScenarioTree)
    return tree.root
end

function last_stage(tree::ScenarioTree)
    return maximum(keys(tree.stage_map))
end

function is_leaf(node::ScenarioNode)
    return length(node.children) == 0 ||
        (length(node.children) == 1 && is_leaf(first(node.children)))
end

function is_leaf(tree::ScenarioTree, nid::NodeID)
    return is_leaf(tree.tree_map[nid])
end

function _add_node(tree::ScenarioTree, node::ScenarioNode)
    tree.tree_map[node.id] = node
    if node.stage in keys(tree.stage_map)
        push!(tree.stage_map[node.stage], node.id)
    else
        tree.stage_map[node.stage] = Set{NodeID}([node.id])
    end
    return
end

"""
    add_node(tree::ScenarioTree, parent::ScenarioNode)

Add a node to the ScenarioTree `tree` with parent node `parent`. Return the added node. If the node to add is a leaf, use `add_leaf` instead.
"""
function add_node(tree::ScenarioTree, parent::ScenarioNode)
    new_node = _create_node(tree.id_gen, parent)
    _add_node(tree, new_node)
    return new_node
end

function _add_scenario_to_bundle(tree::ScenarioTree, nid::NodeID, scid::ScenarioID)
    push!(tree.tree_map[nid].scenario_bundle, scid)
    return
end

"""
    add_leaf(tree::ScenarioTree, parent::ScenarioNode, probability<:Real)

Add a leaf to the ScenarioTree `tree` with parent node `parent`. The probability of this scenario occuring is given by `probability`. Returns the ScenarioID representing the scenario.
"""
function add_leaf(tree::ScenarioTree, parent::ScenarioNode, probability::R
                  ) where R <: Real
    leaf = add_node(tree, parent)
    scid = _assign_scenario_id(tree)
    tree.prob_map[scid] = probability

    node = leaf
    while node != nothing
        id = node.id
        _add_scenario_to_bundle(tree, id, scid)
        node = node.parent
    end
    return scid
end

_assign_scenario_id(tree::ScenarioTree)::ScenarioID = _generate_scenario_id(tree.id_gen)

function scenario_bundle(node::ScenarioNode)::Set{ScenarioID}
    return node.scenario_bundle
end

function scenario_bundle(tree::ScenarioTree, nid::NodeID)::Set{ScenarioID}
    return scenario_bundle(tree.tree_map[nid])
end

scenarios(tree::ScenarioTree) = tree.root.scenario_bundle

mutable struct VariableInfo
    ref::Union{Future,JuMP.VariableRef}
    name::String
    node_id::NodeID
    value::Float64
end

function VariableInfo(ref::Union{Future,JuMP.VariableRef},
                      name::String,
                      nid::NodeID)
    return VariableInfo(ref, name, nid, 0.0)
end

function value(vi::VariableInfo)::Float64
    return vi.value
end

mutable struct PHHatVariable
    value::Float64
end

# PHHatVariable() = PHHatVariable(Set{ScenarioID}(), 0.0)
PHHatVariable() = PHHatVariable(0.0)

function value(a::PHHatVariable)
    return a.value
end

mutable struct RefValuePair
    ref::Union{Future,Nothing}
    value::Float64
end

function RefValuePair()
    return RefValuePair(nothing, 0.0)
end

function ref(a::RefValuePair)::Union{Future,Nothing}
    return a.ref
end

function value(a::RefValuePair)::Float64
    return a.value
end

struct ScenarioVariableRecord
    p::Float64
    x::VariableInfo
    w::RefValuePair
end

function ScenarioVariableRecord(p::Float64, x::VariableInfo)::ScenarioVariableRecord
    return ScenarioVariableRecord(p, x, RefValuePair())
end

struct PHVariableRecord
    scen_bundle::Dict{ScenarioID,ScenarioVariableRecord}
    xhat::PHHatVariable
end

function PHVariableRecord()::PHVariableRecord
    return PHVariableRecord(Dict{ScenarioID,ScenarioVariableRecord}(), PHHatVariable())
end

function add_scenario_record(phvr::PHVariableRecord,
                             s::ScenarioID,
                             srec::ScenarioVariableRecord)
    phvr.scen_bundle[s] = srec
    return
end

function is_leaf_record(phvr::PHVariableRecord)
    return length(phvr.scen_bundle) == 1
end

function xhat_value(phvr::PHVariableRecord)::Float64
    return value(phvr.xhat)
end

function set_xhat_value(phvr::PHVariableRecord, value::Float64)::Nothing
    phvr.xhat.value = value
    return
end

struct ScenarioInfo
    proc::Int
    prob::Float64
    model::Future
    branch_vars::Dict{VariableID, VariableInfo}
    leaf_vars::Dict{VariableID, VariableInfo}
    w_vars::Dict{VariableID, RefValuePair}
    xhat_vars::Dict{XhatID, RefValuePair}
end

function retrieve_variable(sinfo::ScenarioInfo, vid::VariableID)::VariableInfo
    if haskey(sinfo.branch_vars, vid)
        vi = sinfo.branch_vars[vid]
    else
        vi = sinfo.leaf_vars[vid]
    end
    return vi
end

struct PHResidualHistory
    residuals::Dict{Int,Float64}
end

function PHResidualHistory()
    return PHResidualHistory(Dict{Int,Float64}())
end

function residual_vector(phrh::PHResidualHistory)::Vector{Float64}
    if length(phrh.residuals) > 0
        max_iter = maximum(keys(phrh.residuals))
        return [phrh.residuals[k] for k in sort!(collect(keys(phrh.residuals)))]
    else
        return Vector{Float64}()
    end
end

function save_residual(phrh::PHResidualHistory, iter::Int, res::Float64)::Nothing
    @assert(!(iter in keys(phrh.residuals)))
    phrh.residuals[iter] = res
    return
end

struct PHData
    r::Float64
    scenario_tree::ScenarioTree
    scenario_view::Dict{ScenarioID, ScenarioInfo}
    ph_view::Dict{XhatID, PHVariableRecord}
    time_info::TimerOutputs.TimerOutput
    residual_info::PHResidualHistory
end

function PHData(r::N, tree::ScenarioTree,
                scen_proc_map::Dict{ScenarioID, Int},
                probs::Dict{ScenarioID, Float64},
                submodels::Dict{ScenarioID, Future},
                var_map::Dict{ScenarioID, Dict{VariableID, VariableInfo}},
                time_out::TimerOutputs.TimerOutput
                ) where {N <: Number}

    scenario_view = Dict{ScenarioID, ScenarioInfo}()
    ph_view = Dict{XhatID, PHVariableRecord}()
    for (scid, model) in pairs(submodels)

        leaf_vars = Dict{VariableID, VariableInfo}()
        branch_vars = Dict{VariableID, VariableInfo}()
        w_vars = Dict{VariableID, RefValuePair}()
        xhat_vars = Dict{XhatID, RefValuePair}()

        p = probs[scid]

        for (vid, vinfo) in var_map[scid]

            xid = XhatID(vinfo.node_id, vid.index)

            if is_leaf(tree, vinfo.node_id)
                leaf_vars[vid] = vinfo
            else
                branch_vars[vid] = vinfo

                if !haskey(ph_view, xid)
                    ph_view[xid] = PHVariableRecord()
                end

                sv_rec = ScenarioVariableRecord(p, vinfo)
                add_scenario_record(ph_view[xid], scid, sv_rec)

                w_vars[vid] = sv_rec.w
                xhat_vars[xid] = RefValuePair()
            end

        end

        scenario_view[scid] = ScenarioInfo(scen_proc_map[scid],
                                           p,
                                           model,
                                           branch_vars,
                                           leaf_vars,
                                           w_vars,
                                           xhat_vars
                                           )

    end

    return PHData(float(r),
                  tree,
                  scenario_view,
                  ph_view,
                  time_out,
                  PHResidualHistory(),
                  )
end

function residuals(phd::PHData)::Vector{Float64}
    return residual_vector(phd.residual_info)
end

function save_residual(phd::PHData, iter::Int, res::Float64)::Nothing
    save_residual(phd.residual_info, iter, res)
    return
end

function scenario_bundle(phd::PHData, xid::XhatID)::Set{ScenarioID}
    return scenario_bundle(phd.scenario_tree, xid.node)
end

function stage_id(phd::PHData, xid::XhatID)::StageID
    return phd.scenario_tree.tree_map[xid.node].stage
end

function num_shared_variables(phd::PHData)::Int
    count = 0
    for (xhid, phv_rec) in pairs(phd.ph_view)
        if !is_leaf_record(phv_rec)
            count += 1
        end
    end
    return count
end

function convert_to_variable_id(phd::PHData, xid::XhatID)
    idx = xid.index
    stage = stage_id(phd, xid)
    scen = first(scenario_bundle(phd, xid))
    return (scen, VariableID(stage, idx))
end

function convert_to_xhat_id(phd::PHData, scid::ScenarioID, vid::VariableID)::XhatID
    vinfo = retrieve_variable(phd.scenario_view[scid], vid)
    return XhatID(vinfo.node_id, vid.index)
end
