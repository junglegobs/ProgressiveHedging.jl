
function _report_values(var_dict::Dict{VariableID, VariableInfo}
                        )::Dict{VariableID, Float64}
    val_dict = Dict{VariableID, Float64}()
    for (vid, vinfo) in pairs(var_dict)
        val_dict[vid] = JuMP.value(fetch(vinfo.ref))
    end
    return val_dict
end

function retrieve_values(phd::PHData, leaf_mode::Bool)::Nothing

    map_sym = leaf_mode ? :leaf_vars : :branch_vars

    val_dict = Dict{ScenarioID, Future}()
    for (scid,sinfo) in phd.scenario_view
        vmap = getfield(sinfo, map_sym)
        val_dict[scid] = @spawnat(sinfo.proc, _report_values(vmap))
    end

    for (scid,fv) in pairs(val_dict)
        sinfo = phd.scenario_view[scid]
        var_values = fetch(fv)
        for (vid, value) in pairs(var_values)
            vmap = getfield(sinfo, map_sym)
            vmap[vid].value = value
        end
    end

    return
end

function compute_and_save(phd::PHData)::Tuple{Float64,Float64}

    x_res_sq = zeros(Threads.nthreads())
    w_res_sq = zeros(Threads.nthreads())

    Threads.@threads for xhid in collect(keys(phd.ph_view))

        phv_rec = phd.ph_view[xhid]
        @assert(length(phv_rec.scen_bundle) > 1)

        # Update xhat
        xhat = 0.0
        norm = 0.0

        for (s, sv_rec) in pairs(phv_rec.scen_bundle)
            xhat += sv_rec.p * value(sv_rec.x)
            norm += sv_rec.p
        end

        xhat = xhat / norm
        x_res_sq[Threads.threadid()] += (xhat - xhat_value(phv_rec))^2
        set_xhat_value(phv_rec, xhat)

        # Update w
        exp = 0.0

        for (s,sv_rec) in pairs(phv_rec.scen_bundle)
            kx = value(sv_rec.x) - xhat
            sv_rec.w.value += phd.r * kx

            w_res_sq[Threads.threadid()] += sv_rec.p * kx^2

            exp += sv_rec.p * sv_rec.w.value
        end

        if abs(exp) > 1e-6
            @warn("Conditional expectation of " *
                  "W[$(node.scenario_bundle),$(node.stage),$i] " *
                  "is non-zero: " * string(exp/norm))
        end

    end

    return (sum(x_res_sq), sum(w_res_sq))
end

function update_ph_variables(phd::PHData)::Tuple{Float64,Float64}
    retrieve_values(phd, false)
    return compute_and_save(phd)
end

function update_ph_leaf_variables(phd::PHData)::Nothing

    retrieve_values(phd, true)

    for (scid, sinfo) in pairs(phd.scenario_view)
        for (vid, vinfo) in pairs(sinfo.leaf_vars)
            xhid = convert_to_xhat_id(phd, scid, vid)
            @assert(!haskey(phd.ph_view, xhid))

            phvr = PHVariableRecord()
            set_xhat_value(phvr, vinfo.value)
            phd.ph_view[xhid] = phvr
        end
    end

    return
end

# Some MOI interfaces do not support setting of start values
function _set_start_values(model::JuMP.Model)::Nothing
    for var in JuMP.all_variables(model)
        if !JuMP.is_fixed(var)
            JuMP.set_start_value(var, JuMP.value(var))
        end
    end
    return
end

function set_start_values(phd::PHData)::Nothing

    @sync for (scid, sinfo) in pairs(phd.scenario_view)
        model = sinfo.model
        @spawnat(sinfo.proc, _set_start_values(fetch(model)))
    end

    return

end

function _fix_values(ph_vars::Vector{RefValuePair})::Nothing

    for phv in ph_vars
        JuMP.fix(fetch(phv.ref), phv.value, force=true)
    end

    return
end

function update_si_xhat(phd::PHData)::Nothing
    for (scid, sinfo) in pairs(phd.scenario_view)
        for (xhat_id, xhat_var) in pairs(sinfo.xhat_vars)
            xhat_var.value = value(phd.ph_view[xhat_id].xhat)
        end
    end
    return
end

function fix_ph_variables(phd::PHData)::Nothing

    update_si_xhat(phd)

    @sync for (scid, sinfo) in pairs(phd.scenario_view)
        w_array = collect(values(sinfo.w_vars))
        xhat_array = collect(values(sinfo.xhat_vars))

        @spawnat(sinfo.proc, _fix_values(w_array))
        @spawnat(sinfo.proc, _fix_values(xhat_array))
    end

    return
end

function solve_subproblems(phd::PHData)::Nothing

    # Find subproblem solutions--in parallel if we have the workers for it.
    # @sync will wait for all processes to complete
    @sync for (scen, sinfo) in pairs(phd.scenario_view)
        model = sinfo.model
        @spawnat(sinfo.proc, JuMP.optimize!(fetch(model)))
    end

    for (scen, sinfo) in pairs(phd.scenario_view)
        # MOI refers to the MathOptInterface package. Apparently this is made
        # accessible by JuMP since it is not imported here
        model = sinfo.model
        sts = fetch(@spawnat(sinfo.proc, JuMP.termination_status(fetch(model))))
        if sts != MOI.OPTIMAL && sts != MOI.LOCALLY_SOLVED &&
            sts != MOI.ALMOST_LOCALLY_SOLVED
            @error("Scenario $scen subproblem returned $sts.")
        end
    end

    return
end

function hedge(ph_data::PHData,
               max_iter::Int,
               atol::Float64,
               rtol::Float64,
               report::Int,
               save_res::Bool,
               warm_start::Bool,
               )::Tuple{Int,Float64}
    niter = 0
    report_flag = (report > 0)

    (xhat_res_sq, x_res_sq) = @timeit(ph_data.time_info, "Update PH Vars",
                                      update_ph_variables(ph_data))

    nsqrt = sqrt(length(ph_data.ph_view))
    xmax = max(maximum(abs.(xhat_value.(values(ph_data.ph_view)))), 1e-12)
    residual = sqrt(xhat_res_sq + x_res_sq) / nsqrt

    if report_flag
        @printf("Iter: %4d   AbsR: %12.6e   RelR: %12.6e   Xhat: %12.6e   X: %12.6e\n",
                niter, residual, residual/xmax,
                sqrt(xhat_res_sq)/nsqrt, sqrt(x_res_sq)/nsqrt
                )
        flush(stdout)
    end

    if save_res
        save_residual(ph_data, 0, residual)
    end
    
    while niter < max_iter && residual > atol && residual > rtol * xmax
        
        # Set initial values, fix cross model values (W and Xhat) and
        # solve the subproblems

        # Setting start values causes issues with some solvers
        if warm_start
            @timeit(ph_data.time_info, "Set start values",
                    set_start_values(ph_data))
        end
        @timeit(ph_data.time_info, "Fix PH variables",
                fix_ph_variables(ph_data))

        @timeit(ph_data.time_info, "Solve subproblems",
                solve_subproblems(ph_data))

        # Update xhat and w
        (xhat_res_sq, x_res_sq) = @timeit(ph_data.time_info, "Update PH Vars",
                                          update_ph_variables(ph_data))

        # Update stopping criteria -- xhat_res_sq measures the movement of
        # xhat values from k^th iteration to the (k+1)^th iteration while
        # x_res_sq measures the disagreement between the x variables and
        # its corresponding xhat variable (so lack of consensus amongst the
        # subproblems or violation of the nonanticipativity constraint)
        residual= sqrt(xhat_res_sq + x_res_sq) / nsqrt
        xmax = max(maximum(abs.(xhat_value.(values(ph_data.ph_view)))), 1e-12)
        
        niter += 1

        if report_flag && niter % report == 0
            @printf("Iter: %4d   AbsR: %12.6e   RelR: %12.6e   Xhat: %12.6e   X: %12.6e\n",
                    niter, residual, residual/xmax,
                    sqrt(xhat_res_sq)/nsqrt, sqrt(x_res_sq)/nsqrt
                    )
            flush(stdout)
        end

        if save_res
            save_residual(ph_data, niter, residual)
        end

    end

    @timeit(ph_data.time_info, "Update PH leaf variables",
            update_ph_leaf_variables(ph_data))

    if niter >= max_iter && residual > atol
        @warn("Performed $niter iterations without convergence. " *
              "Consider increasing max_iter from $max_iter.")
    end

    return (niter, residual)
end
