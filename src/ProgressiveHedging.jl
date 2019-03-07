module ProgressiveHedging

import JuMP
import StructJuMP
import DataFrames

import MathOptInterface
const MOI = MathOptInterface

#export solve

# TODO: Find a way to return the objective value for each scenario

# TODO: Pass in actual error function for adding things in case something
# goes wrong. Currently just a function stub.

include("structs.jl")
include("utils.jl")

include("algorithm.jl")
include("setup.jl")

function solve(root_model::StructJuMP.StructuredModel,
               optimizer_factory::JuMP.OptimizerFactory,
               r::T; model_type::Type{M}=JuMP.Model, max_iter=500, atol=1e-8
               ) where {T <: Number, M <: JuMP.AbstractModel}
    # Initialization
    ph_data = initialize(root_model, r, optimizer_factory, M)
    
    # Solution
    niter = 0
    residual = atol + 1.0e10
    
    while niter < max_iter && residual > atol
        # Update Xhat values
        compute_and_save_xhat(ph_data)
        
        # Set initial values, fix cross model values (W and Xhat) and
        # solve the subproblems
        set_start_values(ph_data)
        fix_ph_variables(ph_data)
        solve_subproblems(ph_data)

        # Update W values
        compute_and_save_w(ph_data)

        # Update stopping criteria
        residual = compute_residual(ph_data)
        niter += 1
    end

    if niter >= max_iter
        @warn("Performed $niter iterations without convergence. " *
              "Consider increasing max_iter from $max_iter.")
    end

    # Post Processing
    #(soln_df, cost_dict) = retrieve_soln(ph_data)
    soln_df = retrieve_soln(ph_data)
    
    # return (niter, residual, soln_df, cost_dict, ph_data)
    return (niter, residual, soln_df, ph_data)
end

function build_extensive_form(root_model::StructJuMP.StructuredModel,
                              model::M) where {M <: JuMP.AbstractModel}
    (scen_tree, probs) = build_scenario_tree(root_model)

    sj_models = [root_model]
    probs = [1.0]
    obj = JuMP.GenericAffExpr{Float64, JuMP.variable_type(model)}()

    last_used = Index(0)
    var_map = Dict{VariableID, VariableInfo{JuMP.variable_type(model)}}()
    name_map = Dict{VariableID, String}()
    var_translator = Dict{NodeID, Dict{Int,Index}}()

    while !isempty(sj_models)
        sjm = popfirst!(sj_models)
        p = popfirst!(probs)

        for (id, cmod) in pairs(sjm.children)
            push!(sj_models, cmod)
            # Here's that Markov assumption again
            push!(probs, p * sjm.probability[id])
        end

        # Add variables
        node = translate(scen_tree, sjm)
        var_translator[node.id] = Dict{Int,Index}()
        scid = add_variables_extensive(model, var_map, name_map,
                                       scen_tree, var_translator,
                                       sjm, last_used)
        last_used = maximum(values(var_translator[node.id]))

        # Add constraints
        copy_constraints(model, sjm, scen_tree, var_translator,
                         scid, var_map)

        # Add to objective function
        obj += p * convert_expression(JuMP.variable_type(model),
                                      sjm.objective_function,
                                      scen_tree, var_translator,
                                      scid, var_map)
    end
    
    # Add objective function
    JuMP.set_objective(model, root_model.objective_sense, obj)

    return model
end

function solve_extensive_form(root_model::StructJuMP.StructuredModel,
                              optimizer_factory::JuMP.OptimizerFactory;
                              model::M=JuMP.Model()
                              ) where {M <: JuMP.AbstractModel}
    build_extensive_form(root_model, model)
    JuMP.optimize!(model, optimizer_factory)
    return model
end

end # module
