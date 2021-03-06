
@testset "Penalty Parameter Interface" begin
    struct Fake <: PH.AbstractPenaltyParameter end
    fake = Fake()
    xhid = PH.XhatID(PH.NodeID(0), PH.Index(0))
    phd = fake_phdata(3, 4; add_leaves=false)
    @test_throws PH.UnimplementedError PH.get_penalty_value(fake)
    @test_throws PH.UnimplementedError PH.get_penalty_value(fake, xhid)
    @test_throws PH.UnimplementedError PH.is_initial_value_dependent(Fake)
    @test_throws PH.UnimplementedError PH.is_subproblem_dependent(Fake)
    @test_throws PH.UnimplementedError PH.is_variable_dependent(Fake)
    @test_throws PH.UnimplementedError PH.penalty_map(fake)
    @test_throws PH.UnimplementedError PH.process_penalty_initial_value(fake, phd)
    @test_throws PH.UnimplementedError PH.process_penalty_subproblem(fake,
                                                                     phd,
                                                                     PH.scid(0),
                                                                     Dict{PH.VariableID,
                                                                          Float64}()
                                                                  )
end

@testset "Scalar Penalty" begin
    @test PH.is_initial_value_dependent(PH.ScalarPenaltyParameter) == false
    @test PH.is_subproblem_dependent(PH.ScalarPenaltyParameter) == false
    @test PH.is_variable_dependent(PH.ScalarPenaltyParameter) == false

    rflt = 10.0*rand()
    r = PH.ScalarPenaltyParameter(rflt)
    @test PH.get_penalty_value(r) == rflt
    @test PH.get_penalty_value(r, PH.XhatID(PH.NodeID(0), PH.Index(0))) == rflt
end

@testset "Proportional Penalty" begin
    @test PH.is_initial_value_dependent(PH.ProportionalPenaltyParameter) == false
    @test PH.is_subproblem_dependent(PH.ProportionalPenaltyParameter) == true
    @test PH.is_variable_dependent(PH.ProportionalPenaltyParameter) == true

    N = 2
    st = PH.two_stage_tree(N)
    var_map = Dict{PH.ScenarioID,Dict{PH.VariableID,PH.VariableInfo}}()
    coef_map = Dict{Int,Dict{PH.VariableID,Float64}}()
    for s in 0:(N-1)
        js = two_stage_model(PH.scid(s))
        svm = PH.report_variable_info(js, st)
        var_map[PH.scid(s)] = svm
        (br_vars, lf_vars) = PH._split_variables(st, collect(keys(svm)))
        coef_map[s] = PH.report_penalty_info(js,
                                             br_vars,
                                             PH.ProportionalPenaltyParameter
                                             )
    end

    c_rho = 10.0 * (rand() + 1e-12) # Avoid the possibility of c_rho == 0.0
    phd = PH.PHData(PH.ProportionalPenaltyParameter(c_rho),
                    st,
                    Dict(1 => Set([PH.scid(0), PH.scid(1)])),
                    var_map,
                    TimerOutputs.TimerOutput()
                    )

    for s in 0:(N-1)
        PH.process_penalty_subproblem(phd.r, phd, PH.scid(s), coef_map[s])
    end

    for (xhid, rho) in pairs(phd.r.penalties)
        @test st.tree_map[xhid.node].stage == PH.stid(1)
        @test isapprox(rho, c_rho) || isapprox(rho, phd.r.constant)
    end
end

@testset "SEP Penalty" begin
    @test PH.is_initial_value_dependent(PH.SEPPenaltyParameter) == true
    @test PH.is_subproblem_dependent(PH.SEPPenaltyParameter) == true
    @test PH.is_variable_dependent(PH.SEPPenaltyParameter) == true

    N = 2
    st = PH.two_stage_tree(N)
    var_map = Dict{PH.ScenarioID,Dict{PH.VariableID,PH.VariableInfo}}()
    var_vals = Dict{Int,Dict{PH.VariableID,Float64}}()
    coef_map = Dict{Int,Dict{PH.VariableID,Float64}}()
    for s in 0:(N-1)
        js = two_stage_model(PH.scid(s))
        svm = PH.report_variable_info(js, st)
        var_map[PH.scid(s)] = svm
        (br_vars, lf_vars) = PH._split_variables(st, collect(keys(svm)))
        PH.solve_subproblem(js)
        var_vals[s] = PH.report_values(js, collect(keys(svm)))
        coef_map[s] = PH.report_penalty_info(js,
                                             br_vars,
                                             PH.SEPPenaltyParameter
                                             )
    end

    phd = PH.PHData(PH.SEPPenaltyParameter(),
                    st,
                    Dict(1 => Set([PH.scid(0), PH.scid(1)])),
                    var_map,
                    TimerOutputs.TimerOutput()
                    )

    for s in 0:(N-1)
        PH._copy_values(phd.scenario_map[PH.scid(s)].branch_vars, var_vals[s])
        PH.process_penalty_subproblem(phd.r, phd, PH.scid(s), coef_map[s])
    end

    PH.compute_and_save_xhat(phd)
    PH.process_penalty_initial_value(phd.r, phd)

    value = 0.5 * (10.0 + 3.0)
    value = 0.5 * (10.0 - value) + 0.5 * (value - 3.0)
    rho_val = 1.0 / value
    for (xhid, rho) in pairs(phd.r.penalties)
        @test st.tree_map[xhid.node].stage == PH.stid(1)
        @test isapprox(rho, rho_val) || isapprox(rho, phd.r.default)
    end

    if isdefined(Main, :Xpress)
        N = 8
        st = PH.two_stage_tree(N)
        var_map = Dict{PH.ScenarioID,Dict{PH.VariableID,PH.VariableInfo}}()
        var_vals = Dict{Int,Dict{PH.VariableID,Float64}}()
        coef_map = Dict{Int,Dict{PH.VariableID,Float64}}()
        for s in 0:(N-1)
            js = two_stage_int(PH.scid(s))
            svm = PH.report_variable_info(js, st)
            var_map[PH.scid(s)] = svm
            (br_vars, lf_vars) = PH._split_variables(st, collect(keys(svm)))
            PH.solve_subproblem(js)
            var_vals[s] = PH.report_values(js, collect(keys(svm)))
            coef_map[s] = PH.report_penalty_info(js,
                                                 br_vars,
                                                 PH.SEPPenaltyParameter
                                                 )
        end

        phd = PH.PHData(PH.SEPPenaltyParameter(),
                        st,
                        Dict(1 => Set([PH.scid(s) for s in 0:(N-1)])),
                        var_map,
                        TimerOutputs.TimerOutput()
                        )

        for s in 0:(N-1)
            PH._copy_values(phd.scenario_map[PH.scid(s)].branch_vars, var_vals[s])
            PH.process_penalty_subproblem(phd.r, phd, PH.scid(s), coef_map[s])
        end

        PH.compute_and_save_xhat(phd)
        PH.process_penalty_initial_value(phd.r, phd)

        value = (10 - 4 + 1)
        rho_val = 2.0 / value
        for (xhid, rho) in pairs(phd.r.penalties)
            @test st.tree_map[xhid.node].stage == PH.stid(1)
            @test isapprox(rho, rho_val)
        end
    else
        println("Skipping SEP integer setting tests.")
    end
end

@testset "Reset penalty parameters" begin
    for r in (PH.ProportionalPenaltyParameter(1.0), PH.SEPPenaltyParameter(1.0))
        @test isempty(r.penalties) == true
        r.penalties[PH.XhatID(PH.NodeID(0), PH.Index(0))] = 1.0
        @test isempty(r.penalties) == false
        PH.reset_penalty_map(r)
        @test isempty(r.penalties) == true
    end
end
