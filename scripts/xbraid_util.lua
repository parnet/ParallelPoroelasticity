---
--- Generated by EmmyLua(https://github.com/EmmyLua)
--- Created by parnet.
--- DateTime: 28.05.21 15:44
---
---ug_load_script("tools/factory.lua")
xbraid_util = xbraid_util or {}

p_mgrit_cycle_type = "V" or "F" -- todo
p_mgrit_relax_type = "FCF" or "F-FCF" or "F" -- todo

function xbraid_util.set_cycle_type(braid, p_mgrit_cycle_type)
    if p_mgrit_cycle_type == "F" then
        print("MGRIT uses F-Cycle")
        braid:set_cycle_fmg()
    elseif p_mgrit_cycle_type == "V" then
        print("MGRIT uses V-Cycle")
    else
        print("Invalid MGRIT-cycle parameter - using V cycle")
    end
end

function xbraid_util.create_integratorFactory(name, domainDiscT, lsolver, nlsolver, biotErrorEst, endTime, orderOrTheta)
    if name == "T" then
        integrator_theta = ThetaIntegratorFactory()
        integrator_theta:set_domain(domainDiscT)
        integrator_theta:set_solver(lsolver)
        integrator_theta:set_level_theta(0, orderOrTheta)
        integrator_theta:set_level_theta(1, orderOrTheta)
        integrator_theta:set_level_theta(2, orderOrTheta)
        integrator_theta:set_level_theta(3, orderOrTheta)
        integrator_theta:set_level_theta(4, orderOrTheta)
        integrator_theta:set_level_theta(5, orderOrTheta)
        integrator_theta:set_level_theta(6, orderOrTheta)
        integrator_theta:set_level_theta(7, orderOrTheta)
        integrator_theta:set_level_theta(8, orderOrTheta)
        integrator_theta:set_level_theta(9, orderOrTheta)
        integrator_theta:set_level_theta(10, orderOrTheta)
        integrator = integrator_theta
        print("XBRAID fine integrator: using theta time step")
    elseif name == "CT" then
        integrator_theta = CachedThetaIntegratorFactory()
        integrator_theta:set_domain(domainDiscT)
        integrator_theta:set_solver(lsolver)
        integrator_theta:set_level_theta(0, orderOrTheta)
        integrator_theta:set_level_theta(1, orderOrTheta)
        integrator_theta:set_level_theta(2, orderOrTheta)
        integrator_theta:set_level_theta(3, orderOrTheta)
        integrator_theta:set_level_theta(4, orderOrTheta)
        integrator_theta:set_level_theta(5, orderOrTheta)
        integrator_theta:set_level_theta(6, orderOrTheta)
        integrator_theta:set_level_theta(7, orderOrTheta)
        integrator_theta:set_level_theta(8, orderOrTheta)
        integrator_theta:set_level_theta(9, orderOrTheta)
        integrator_theta:set_level_theta(10, orderOrTheta)
        integrator = integrator_theta
        print("XBRAID fine integrator: using cached theta time step")
    elseif name == "FT" then
        print("THETA " .. orderOrTheta)
        integrator_theta = FixedStepThetaIntegratorFactory()
        integrator_theta:set_domain(domainDiscT)
        integrator_theta:set_solver(lsolver)
        integrator_theta:set_level_num_steps(0, orderOrTheta)
        integrator_theta:set_level_num_steps(1, orderOrTheta)
        integrator_theta:set_level_num_steps(2, orderOrTheta)
        integrator_theta:set_level_num_steps(3, orderOrTheta)
        integrator_theta:set_level_num_steps(4, orderOrTheta)
        integrator_theta:set_level_num_steps(5, orderOrTheta)
        integrator_theta:set_level_num_steps(6, orderOrTheta)
        integrator_theta:set_level_num_steps(7, orderOrTheta)
        integrator_theta:set_level_num_steps(8, orderOrTheta)
        integrator_theta:set_level_num_steps(9, orderOrTheta)
        integrator_theta:set_level_num_steps(10, orderOrTheta)
        integrator = integrator_theta
        print("XBRAID fine integrator: using fixed step theta time step")
    elseif name == "B" then
        integrator_bdf = BDF_IntegratorFactory()
        integrator_bdf:set_domain(domainDiscT)
        integrator_bdf:set_solver(lsolver)
        integrator_bdf:set_level_order(10, orderOrTheta)
        integrator_bdf:set_level_order(9, orderOrTheta)
        integrator_bdf:set_level_order(8, orderOrTheta)
        integrator_bdf:set_level_order(7, orderOrTheta)
        integrator_bdf:set_level_order(6, orderOrTheta)
        integrator_bdf:set_level_order(5, orderOrTheta)
        integrator_bdf:set_level_order(4, orderOrTheta)
        integrator_bdf:set_level_order(3, orderOrTheta)
        integrator_bdf:set_level_order(2, orderOrTheta)
        integrator_bdf:set_level_order(1, orderOrTheta)
        integrator_bdf:set_level_order(0, orderOrTheta)
        integrator = integrator_bdf
        print("XBRAID fine integrator: using theta time step")
    elseif name == "X" then
        integrator_limex = LimexFactory()
        integrator_limex:set_domain_disc(domainDiscT)
        integrator_limex:set_solver(nlsolver)
        integrator_limex:set_error_estimator(biotErrorEst)
        integrator_limex:set_tol(1e-3)
        integrator_limex:set_dt_min(1e-20)
        integrator = integrator_limex
        print("XBRAID fine integrator: using limex ")
    elseif name == "C" then
        integrator_const = ConstStepLinearTimeIntegratorFactory()
        integrator_const:set_time_disc(ThetaTimeStep(domainDiscT))
        integrator_const:set_num_steps(1)
        integrator_const:set_solver(lsolver)
        integrator = integrator_const
        print("XBRAID fine integrator: using const step linear")
    elseif name == "S" then
        --time_stepper  =  LinearImplicitEuler(domainDiscT)
        integrator_simple = SimpleIntegratorFactory()
        --integrator_simple:set_time_stepper(time_stepper)
        integrator_simple:set_domain_disc(domainDiscT)
        integrator_simple:set_solver(nlsolver)
        integrator_simple:set_dt_min(endTime / 131072)
        integrator_simple:set_dt_max(endTime)
        integrator_simple:set_reduction_factor(0.2)
        integrator = integrator_simple
        print("XBRAID fine integrator: using simple integrator")
    elseif name == "L" then
        integrator_linear = LinearTimeIntegratorFactory()
        integrator_linear:set_time_disc(ThetaTimeStep(domainDiscT))
        integrator_linear:set_solver(lsolver)
        integrator = integrator_linear
        print("XBRAID fine integrator: using linear time integrator")
    elseif name == "A" then
        integrator_adaptive = TimeIntegratorLinearAdaptiveFactory()
        integrator_adaptive:set_time_stepper_1(ThetaTimeStep(domainDiscT))
        integrator_adaptive:set_time_stepper_2(ThetaTimeStep(domainDiscT))
        integrator_adaptive:set_time_step_min(endTime / 131072)
        integrator_adaptive:set_time_step_max(endTime)
        integrator_adaptive:set_tol(1e-3)
        integrator_adaptive:set_level_factor(orderOrTheta)
        integrator = integrator_adaptive
    elseif name == "D" then
        print("XBRAID coarse integrator: Discontinuity not implemented yet.")
        exit()
    else
        print("XBRAID coarse integrator: ERROR T or C or L (future or D)")
        exit()
    end
    return integrator
end

function xbraid_util.createFSTheta(domain, lsolver, theta, num_steps, threshold)
    integrator = FixedStepThetaIntegrator()
    integrator:set_theta(theta)
    integrator:set_num_steps(num_steps)
    integrator:set_reassemble_threshold(threshold)
    integrator:set_domain(domain)
    integrator:set_solver(lsolver)
    return integrator
end

function xbraid_util.createThetaStepper(domain, lsolver, theta, threshold)
    integrator = ThetaStepper()
    integrator:set_theta(theta)
    -- integrator:set_num_steps(num_steps)
    integrator:set_reassemble_threshold(threshold)
    integrator:set_domain(domain)
    integrator:set_solver(lsolver)
    return integrator
end

function xbraid_util.createNLTheta(domain, nlsolver, theta, num_steps, threshold)
    if num_steps == 1 then
        integrator = NLThetaIntegrator()
        integrator:set_domain(domain)
        integrator:set_solver(nlsolver)
        integrator:set_theta(theta)
       -- integrator:set_reassemble_threshold(threshold)
        return integrator
    else
        integrator = NLFixedStepThetaIntegrator()
        integrator:set_domain(domain)
        integrator:set_solver(nlsolver)
        integrator:set_theta(theta)
        integrator:set_num_steps(num_steps)
        -- integrator:set_reassemble_threshold(threshold)
        return integrator
    end

end

function xbraid_util.CreateNLLevel(app, domain, nlsolver, theta, num_steps, threshold)
    app:set_default_integrator(xbraid_util.createNLTheta(domain, nlsolver, theta, num_steps, threshold))

    print("Set Integrator Methods - Leveldependend")
    app:set_integrator(0, xbraid_util.createNLTheta(domain, nlsolver, theta, num_steps, threshold))

    app:set_integrator(1, xbraid_util.createNLTheta(domain, nlsolver, theta, num_steps, threshold))

    app:set_integrator(2, xbraid_util.createNLTheta(domain, nlsolver, theta, num_steps, threshold))

    app:set_integrator(3, xbraid_util.createNLTheta(domain, nlsolver, theta, num_steps, threshold))

    app:set_integrator(4, xbraid_util.createNLTheta(domain, nlsolver, theta, num_steps, threshold))

    app:set_integrator(5, xbraid_util.createNLTheta(domain, nlsolver, theta, num_steps, threshold))

    app:set_integrator(6, xbraid_util.createNLTheta(domain, nlsolver, theta, num_steps, threshold))

    app:set_integrator(7, xbraid_util.createNLTheta(domain, nlsolver, theta, num_steps, threshold))

    app:set_integrator(8, xbraid_util.createNLTheta(domain, nlsolver, theta, num_steps, threshold))

    app:set_integrator(9, xbraid_util.createNLTheta(domain, nlsolver, theta, num_steps, threshold))

    app:set_integrator(10, xbraid_util.createNLTheta(domain, nlsolver, theta, num_steps, threshold))

    app:set_integrator(11, xbraid_util.createNLTheta(domain, nlsolver, theta, num_steps, threshold))

    app:set_integrator(12, xbraid_util.createNLTheta(domain, nlsolver, theta, num_steps, threshold))

    app:set_integrator(13, xbraid_util.createNLTheta(domain, nlsolver, theta, num_steps, threshold))

    app:set_integrator(14, xbraid_util.createNLTheta(domain, nlsolver, theta, num_steps, threshold))

    app:set_integrator(15, xbraid_util.createNLTheta(domain, nlsolver, theta, num_steps, threshold))
end


function xbraid_util.CreateNLLevelFC_old(app, domain, nlsolver_fine,nlsolver_coarse, theta, num_steps, threshold)
    app:set_default_integrator(xbraid_util.createNLTheta(domain, nlsolver_coarse, theta, num_steps, threshold))

    print("Set Integrator Methods - Leveldependend")
    app:set_integrator(0, xbraid_util.createNLTheta(domain, nlsolver_fine, theta, num_steps, threshold))

    app:set_integrator(1, xbraid_util.createNLTheta(domain, nlsolver_coarse, theta, num_steps, threshold))

    app:set_integrator(2, xbraid_util.createNLTheta(domain, nlsolver_coarse, theta, num_steps, threshold))

    app:set_integrator(3, xbraid_util.createNLTheta(domain, nlsolver_coarse, theta, num_steps, threshold))

    app:set_integrator(4, xbraid_util.createNLTheta(domain, nlsolver_coarse, theta, num_steps, threshold))

    app:set_integrator(5, xbraid_util.createNLTheta(domain, nlsolver_coarse, theta, num_steps, threshold))

    app:set_integrator(6, xbraid_util.createNLTheta(domain, nlsolver_coarse, theta, num_steps, threshold))

    app:set_integrator(7, xbraid_util.createNLTheta(domain, nlsolver_coarse, theta, num_steps, threshold))

    app:set_integrator(8, xbraid_util.createNLTheta(domain, nlsolver_coarse, theta, num_steps, threshold))

    app:set_integrator(9, xbraid_util.createNLTheta(domain, nlsolver_coarse, theta, num_steps, threshold))

    app:set_integrator(10, xbraid_util.createNLTheta(domain, nlsolver_coarse, theta, num_steps, threshold))
end

function xbraid_util.CreateNLLevelFC(app, domain, lsolver_fine,lsolver_coarse, theta, num_steps, threshold)
    print("herexe")
    nlsolver_fine = factory.create_nlsolver(lsolver_fine)
    nlsolver_coarse_a = factory.create_nlsolver_coarse(lsolver_coarse)
    nlsolver_coarse_b = factory.create_nlsolver_coarse(lsolver_coarse)
    nlsolver_coarse_c = factory.create_nlsolver_coarse(lsolver_coarse)
    nlsolver_coarse_d = factory.create_nlsolver_coarse(lsolver_coarse)
    nlsolver_coarse_e = factory.create_nlsolver_coarse(lsolver_coarse)
    nlsolver_coarse_f = factory.create_nlsolver_coarse(lsolver_coarse)
    nlsolver_coarse_g = factory.create_nlsolver_coarse(lsolver_coarse)
    nlsolver_coarse_h = factory.create_nlsolver_coarse(lsolver_coarse)
    nlsolver_coarse_i = factory.create_nlsolver_coarse(lsolver_coarse)
    nlsolver_coarse_j = factory.create_nlsolver_coarse(lsolver_coarse)
    nlsolver_coarse_k = factory.create_nlsolver_coarse(lsolver_coarse)
    print("asdads")
    app:set_default_integrator(xbraid_util.createNLTheta(domain, nlsolver_coarse_a, theta, num_steps, threshold))

    print("Set Integrator Methods - Leveldependend")
    app:set_integrator(0, xbraid_util.createNLTheta(domain, nlsolver_fine, theta, num_steps, threshold))

    app:set_integrator(1, xbraid_util.createNLTheta(domain, nlsolver_coarse_b, theta, num_steps, threshold))

    app:set_integrator(2, xbraid_util.createNLTheta(domain, nlsolver_coarse_c, theta, num_steps, threshold))

    app:set_integrator(3, xbraid_util.createNLTheta(domain, nlsolver_coarse_d, theta, num_steps, threshold))

    app:set_integrator(4, xbraid_util.createNLTheta(domain, nlsolver_coarse_e, theta, num_steps, threshold))

    app:set_integrator(5, xbraid_util.createNLTheta(domain, nlsolver_coarse_f, theta, num_steps, threshold))

    app:set_integrator(6, xbraid_util.createNLTheta(domain, nlsolver_coarse_g, theta, num_steps, threshold))

    app:set_integrator(7, xbraid_util.createNLTheta(domain, nlsolver_coarse_h, theta, num_steps, threshold))

    app:set_integrator(8, xbraid_util.createNLTheta(domain, nlsolver_coarse_i, theta, num_steps, threshold))

    app:set_integrator(9, xbraid_util.createNLTheta(domain, nlsolver_coarse_j, theta, num_steps, threshold))

    app:set_integrator(10, xbraid_util.createNLTheta(domain, nlsolver_coarse_k, theta, num_steps, threshold))
end

function xbraid_util.CreateFSLevel(app, domain, lsolver, theta, num_steps, threshold)
    app:set_default_integrator(xbraid_util.createFSTheta(domain, lsolver, theta, num_steps, threshold))
    print("Set Integrator Methods - Leveldependend")
    app:set_integrator(0, xbraid_util.createFSTheta(domain, lsolver, theta, num_steps, threshold))

    app:set_integrator(1, xbraid_util.createFSTheta(domain, lsolver, theta, num_steps, threshold))

    app:set_integrator(2, xbraid_util.createFSTheta(domain, lsolver, theta, num_steps, threshold))

    app:set_integrator(3, xbraid_util.createFSTheta(domain, lsolver, theta, num_steps, threshold))

    app:set_integrator(4, xbraid_util.createFSTheta(domain, lsolver, theta, num_steps, threshold))

    app:set_integrator(5, xbraid_util.createFSTheta(domain, lsolver, theta, num_steps, threshold))

    app:set_integrator(6, xbraid_util.createFSTheta(domain, lsolver, theta, num_steps, threshold))

    app:set_integrator(7, xbraid_util.createFSTheta(domain, lsolver, theta, num_steps, threshold))

    app:set_integrator(8, xbraid_util.createFSTheta(domain, lsolver, theta, num_steps, threshold))

    app:set_integrator(9, xbraid_util.createFSTheta(domain, lsolver, theta, num_steps, threshold))

    app:set_integrator(10, xbraid_util.createFSTheta(domain, lsolver, theta, num_steps, threshold))
end

function xbraid_util.createBDF(domain, lsolver, order, threshold)
    integrator = BDF_Integrator()
    integrator:set_order(order)
    integrator:set_reassemble_threshold(threshold)
    integrator:set_domain(domain)
    integrator:set_solver(lsolver)
    return integrator
end

function xbraid_util.createBDFLevel(app, domain, lsolver, order, threshold)
    app:set_default_integrator(xbraid_util.createBDF(domain,
            lsolver, 4, threshold))
    print("Set Integrator Methods - Leveldependend")
    app:set_integrator(0, xbraid_util.createBDF(domain,
            lsolver, order, threshold))

    app:set_integrator(1, xbraid_util.createBDF(domain,
            lsolver, order, threshold))

    app:set_integrator(2, xbraid_util.createBDF(domain,
            lsolver, order, threshold))

    app:set_integrator(3, xbraid_util.createBDF(domain,
            lsolver, order, threshold))

    app:set_integrator(4, xbraid_util.createBDF(domain,
            lsolver, order, threshold))

    app:set_integrator(5, xbraid_util.createBDF(domain,
            lsolver, order, threshold))

    app:set_integrator(6, xbraid_util.createBDF(domain,
            lsolver, order, threshold))

    app:set_integrator(7, xbraid_util.createBDF(domain,
            lsolver, order, threshold))

    app:set_integrator(8, xbraid_util.createBDF(domain,
            lsolver, order, threshold))

    app:set_integrator(9, xbraid_util.createBDF(domain,
            lsolver, order, threshold))

    app:set_integrator(10, xbraid_util.createBDF(domain,
            lsolver, order, threshold))
end

function xbraid_util.set_relax_type(braid, p_mgrit_relax_type)
    if p_mgrit_relax_type == "F" then
        braid:set_n_relax(-1, 0)
        braid:set_n_relax(0, 0)
        print("MGRIT uses F-Relaxation")
    elseif p_mgrit_relax_type == "FCF" then
        braid:set_n_relax(-1, 1)
        braid:set_n_relax(0, 1)
        print("MGRIT uses FCF-Relaxation on all level")
    elseif p_mgrit_relax_type == "FFCF" then
        braid:set_n_relax(-1, 1)
        braid:set_n_relax(0, 0)
        print("MGRIT uses FCF-Relaxation on coarse level and F-Relaxation on fines level ")
    elseif p_mgrit_relax_type == "FCFF" then
        braid:set_n_relax(-1, 0)
        braid:set_n_relax(0, 1)
        print("MGRIT uses FCF-Relaxation on coarse level and F-Relaxation on fines level ")
    else
        print("Invalid MGRIT-relax parameter - using FCF relaxation")
    end
end

function coarsening_strategie(numTemporalProcs, initial_time_steps, cfactor_communication, cfactor_innerproc)

end

braid_adaptive_conv = {
    exactSolver = p_useResidual,
    forceConvergence = false,
    adaptiveSolver = p_adaptiveSolver,
    looseTol = 5e-6,
    tightTol = 5e-8,
    strongfirst = p_strongfirst,
}

function xbraid_util.SetConvCheck(braid, desc)
    braid:set_max_iterations(desc.conv_check.max_iter)
    if (desc.conv_check.absolute ~= nil) then
        braid:set_absolute_tol(desc.conv_check.absolute)
    else
        braid:set_relative_tol(desc.conv_check.reduction)
    end
end

function xbraid_util.CreateTimeHierarchy(braid, desc)
    if desc.default_cfactor ~= nil then
        braid:set_c_factor(-1, desc.default_cfactor)
    end

    if type(desc.cfactor) == "table" then
        for key, value in pairs(desc.cfactor) do
            braid:set_c_factor(key - 1, value)
        end
    else
        braid:set_c_factor(-1, desc.cfactor)
    end
end

function xbraid_util.SetBaseValues(app, desc, domain, scriptor)
    app:set_verbose(desc.verbose)
    app:set_start_time(desc.time.t_0)
    app:set_end_time(desc.time.t_end)
    app:set_number_of_timesteps(desc.time.n)
    app:set_time_values(desc.time.t_0, desc.time.t_end, desc.time.n)
    app:set_domain(domain)
    app:set_max_levels(desc.max_level)
    app:set_scriptor(scriptor)
end

function xbraid_util.CreateIntegratorFactory(desc, domain, scriptor)
    app = BraidIntegratorFactory()
    xbraid_util.SetBaseValues(app, desc, domain, scriptor)
    -- app:set_fine_time_integrator(fintegrator)
    -- app:set_coarse_time_integrator(cintegrator)
    return app
end

function xbraid_util.CreateIntegrator(desc, domain, scriptor)
    print("create app")
    app = BraidIntegrator()
    print("app created")
    xbraid_util.SetBaseValues(app, desc, domain, scriptor)
    print("base values set")
    return app
end


function xbraid_util.CreateNLIntegrator(desc, domain, scriptor)
    print("create app")
    app = BraidNLIntegrator()
    print("app created")
    xbraid_util.SetBaseValues(app, desc, domain, scriptor)
    print("base values set")
    return app
end



function xbraid_util.CreateTimeStepper(desc, domain, scriptor, solver)
    app = BraidTimeStepper()
    xbraid_util.SetBaseValues(app, desc, domain, scriptor)
    return app
end

function xbraid_util.CreateStepperLevel(app, domainDiscT, lsolver, theta, threshold)
    app:set_integrator(0, xbraid_util.createThetaStepper(domainDiscT,
            lsolver, theta, threshold))

    app:set_integrator(1, xbraid_util.createThetaStepper(domainDiscT,
            lsolver, theta, threshold))

    app:set_integrator(2, xbraid_util.createThetaStepper(domainDiscT,
            lsolver, theta, threshold))

    app:set_integrator(3, xbraid_util.createThetaStepper(domainDiscT,
            lsolver, theta, threshold))

    app:set_integrator(4, xbraid_util.createThetaStepper(domainDiscT,
            lsolver, theta, threshold))

    app:set_integrator(5, xbraid_util.createThetaStepper(domainDiscT,
            lsolver, theta, threshold))

    app:set_integrator(6, xbraid_util.createThetaStepper(domainDiscT,
            lsolver, theta, threshold))

    app:set_integrator(7, xbraid_util.createThetaStepper(domainDiscT,
            lsolver, theta, threshold))

    app:set_integrator(8, xbraid_util.createThetaStepper(domainDiscT,
            lsolver, theta, threshold))

    app:set_integrator(9, xbraid_util.createThetaStepper(domainDiscT,
            lsolver, theta, threshold))

    app:set_integrator(10, xbraid_util.createThetaStepper(domainDiscT,
            lsolver, theta, threshold))

end

function xbraid_util.CreateResidualStepper(desc, domain, scriptor, solver)
    app = BraidResidualStepper()
    xbraid_util.SetBaseValues(app, desc, domain, scriptor)
    app:set_solver(solver)
    return app
end

function xbraid_util.CreateExecutor(desc, communicator, app, logging)
    braid = BraidExecutor(communicator, app)
    braid:set_residual(desc.use_residual)
    braid:set_temporal_norm(desc.temporal_norm)

    xbraid_util.set_relax_type(braid, desc.mgrit_relax_type)
    -- braid:set_n_relax()
    xbraid_util.set_cycle_type(braid, desc.mgrit_cycle_type)
    -- braid:set_cycle_fmg()
    -- braid:set_cycle_nfmg()
    -- braid:set_cycle_nfmgv()
    xbraid_util.CreateTimeHierarchy(braid, desc)
    xbraid_util.SetConvCheck(braid, desc)
    -- braid:set_max_iterations(desc.conv_check.max_iter)
    braid:set_access_level(desc.access_level)
    braid:set_print_level(desc.print_level)
    braid:set_store_values(desc.store_values)
    braid:set_skip_downcycle_work(desc.skip_downcycle_work)
    braid:set_max_levels(desc.max_level)
    braid:set_min_coarse(desc.min_coarsening)
    braid:set_sequential(desc.sequential)
    braid:set_spatial_coarsen_and_refine(desc.spatial_coarsen_and_refine)
    braid:set_refine(desc.time_refinement)
    braid:set_max_refinements(desc.max_refinement)
    braid:set_print_file(desc.printfile)
    --todo braid:set_default_print_file()
    -- braid:set_output()
    braid:set_filename(desc.outputfile)
    braid:set_paralog(logging)

    if braid_desc.sync then
        braid:set_sync()
    end
    --braid:set_increase_max_levels()
    --braid:set_relax_only_cg()
    --braid:set_agg_c_factor()
    --braid:set_periodic()
    --braid:set_final_fc_relax()
    --braid:set_reverted_ranks()
    --braid:set_file_io_level()
    --braid:set_c_relax_weight()
    --braid:set_t_points_cutoff()
    --braid:set_full_residual_norm()
    --braid:set_time_grid()
    --braid:get_num_iteration()
    --braid:get_c_factor()
    --braid:get_residual_norms()
    --braid:get_num_level()
    --braid:get_warm_restart()
    --braid:get_distribution_lower()
    --braid:get_distribution_upper()
    --braid:get_distribution()
    --braid:get_id()
    --braid:get_app()
    --braid:set_settings()
    --braid:get_settings()
    --braid:apply_settings()
    --braid:create_access()
    --braid:set_filename()
    --braid:print_settings()
    --braid:print_summary()
    --braid:test()
    --braid:set_output()CreateIntegrator
    --braid:set_initializer()
    --braid:set_norm_provider()
    --braid:set_paralog()
    --braid:set_paralog_script()
    braid:set_richardson_estimation(desc.richardson_estimation, desc.richardson_extrapolation, desc.richardson_local_order)

    app:init()
    braid:print_settings()
    app:print_settings()
    print("XBraid Integrator created")
    return braid
end

function xbraid_util.SetFullVTKOutput(app)
    app:set_vtk_scriptor(scriptor)

    app:set_vtk_ustart_before(VTKScriptor(VTKOutput(), "u_step_before"))
    app:set_vtk_ustart_after(VTKScriptor(VTKOutput(), "u_step_after"))

    app:set_vtk_uend_before(VTKScriptor(VTKOutput(), "uend_before"))
    app:set_vtk_uend_after(VTKScriptor(VTKOutput(), "uend_after"))

    app:set_vtk_rhs(VTKScriptor(VTKOutput(), "rhs"))
    app:set_vtk_fstop(VTKScriptor(VTKOutput(), "fstop_p" .. communicator:get_temporal_rank()))

    app:set_vtk_resr_before(VTKScriptor(VTKOutput(), "resr_before"))
    app:set_vtk_resr_after(VTKScriptor(VTKOutput(), "resr_after"))

    app:set_vtk_resu_before(VTKScriptor(VTKOutput(), "resu_before"))
    app:set_vtk_resu_after(VTKScriptor(VTKOutput(), "resu_after"))

    app:set_vtk_norm(VTKScriptor(VTKOutput(), "norm_p" .. communicator:get_temporal_rank()))
end

--    app:set_approx_space(approxspace)

function xbraid_util.get_cfactor(str_cfactor)
    -- zero index is finest level
    cfactor = {}
    i = 1
    for s in str_cfactor:gmatch("[^_]+") do
        cfactor[i] = tonumber(s)
        i = i + 1
    end
    return cfactor
end
