---
--- Generated by EmmyLua(https://github.com/EmmyLua)
--- Created by parnet.
--- DateTime: 28.05.21 15:44
---
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


function xbraid_util.create_integrator(name, domainDiscT, lsolver, nlsolver, biotErrorEst,endTime, orderOrTheta)
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
        integrator_bdf:set_level_order(10,orderOrTheta)
        integrator_bdf:set_level_order(9,orderOrTheta)
        integrator_bdf:set_level_order(8,orderOrTheta)
        integrator_bdf:set_level_order(7,orderOrTheta)
        integrator_bdf:set_level_order(6,orderOrTheta)
        integrator_bdf:set_level_order(5,orderOrTheta)
        integrator_bdf:set_level_order(4,orderOrTheta)
        integrator_bdf:set_level_order(3,orderOrTheta)
        integrator_bdf:set_level_order(2,orderOrTheta)
        integrator_bdf:set_level_order(1,orderOrTheta)
        integrator_bdf:set_level_order(0,orderOrTheta)
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
        integrator_simple:set_dt_min(endTime/131072)
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
        integrator_adaptive:set_time_step_min(endTime/131072)
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

braid_integration_desc = {
    time = { t_0 = 0, t_end = 1, n = 100 },
    max_level = 15,
    integrator = limex, -- or table

    coarsening_factor = 2,
    mgrit_cycle_type = "V",
    mgrit_relax_type = "FCF",
    store_values = 0,
    print_level = 3,
    access_level = 1,
    verbose = true,

    richardson_estimation = true, --set_richardson_estimation
    richardson_extrapolation = true,
    richardson_local_order = 4,

    richardson_order = 4,

    temporalNorm = 3, -- {1,2,3}
    conv_check = {
        iterations = 10,
        -- reduction = 1e-9
        absolute = 5e-7
    }
}

function xbraid_util.SetConvCheck(braid, desc)
    braid:set_max_iterations(desc.conv_check.max_iter)
    if (desc.conv_check.absolute ~= nil) then
        braid:set_absolute_tol(desc.conv_check.absolute)
    else
        braid:set_reduction(desc.conv_check.reduction)
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

function xbraid_util.CreateBraidIntegrator(desc, communicator, logging, fintegrator, cintegrator, scriptor, domain)
    -- creating app
    app = BraidIntegratorFactory()
    -- set app base values
    app:set_verbose(desc.verbose)
    app:set_start_time(desc.time.t_0)
    app:set_end_time(desc.time.t_end)
    app:set_number_of_timesteps(desc.time.n)
    app:set_time_values(desc.time.t_0, desc.time.t_end, desc.time.n)
    app:set_domain(domain)
    -- app:set_start_vector()
    -- app:set_scriptor()
    app:set_max_levels(desc.max_level)
    app:set_fine_time_integrator(fintegrator)
    app:set_coarse_time_integrator(cintegrator)
    app:set_scriptor(scriptor)

    --app:set_vtk_scriptor(scriptor)
    --app:set_vtk_ustart_before(VTKScriptor(VTKOutput(),"u_step_before"))
    --app:set_vtk_ustart_after(VTKScriptor(VTKOutput(),"u_step_after"))
    --app:set_vtk_uend_before(VTKScriptor(VTKOutput(),"uend_before"))
    --app:set_vtk_uend_after(VTKScriptor(VTKOutput(),"uend_after"))

    app:set_vtk_residuum(VTKScriptor(VTKOutput(),"cresiduum_p" .. communicator:get_temporal_rank()))
    app:set_vtk_norm(VTKScriptor(VTKOutput(),"norm_p" .. communicator:get_temporal_rank()))
    -- set app specific values
    -- app:set_adapt_convergence()
    -- todo set time integration method

    -- creating executor
    braid = BraidExecutor(communicator, app)
    -- braid:apply()
    -- braid:test()
    -- braid:set_app()
    -- braid:get_app()

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
    braid:set_max_levels(desc.max_level) -- todo dublicate?
    braid:set_min_coarse(desc.min_coarsening)
    braid:set_sequential(desc.sequential)
    braid:set_spatial_coarsen_and_refine(desc.spatial_coarsen_and_refine)
    braid:set_refine(desc.time_refinement)
    braid:set_max_refinements(desc.max_refinement)
    braid:set_print_file(desc.printfile)
    -- braid:set_output()
    braid:set_filename(desc.outputfile)
    braid:set_paralog(logging)

    --braid:set_richardson_estimation(desc.richardson_estimation, desc.richardson_extrapolation, desc.richardson_local_order)

    app:init()
    braid:print_settings()
    app:print_settings()
    print("XBraid Integrator created")

    return braid
end

function xbraid_util.CreateBraidStepper(
        desc, communicator, logging, solver, scriptor, domain, approxspace)
    --cmin,cmax, fmin,fmax, iter)
    -- creating app
    app = BraidTimeStepper()
    -- set app base values

    -- cmin, cmax, fmin, fmax, iteration
    --app:set_iteration_param(cmin,cmax, fmin,fmax, iter)
    app:set_verbose(desc.verbose)

    app:set_start_time(desc.time.t_0)
    app:set_end_time(desc.time.t_end)
    app:set_number_of_timesteps(desc.time.n)
    app:set_time_values(desc.time.t_0, desc.time.t_end, desc.time.n)

    app:set_domain(domain)
    app:set_approx_space(approxspace)
    -- app:set_start_vector()
    -- app:set_scriptor()
    app:set_max_levels(desc.max_level)
    --app:set_fine_time_integrator(fintegrator)
    --app:set_coarse_time_integrator(cintegrator)
    app:set_scriptor(scriptor)
    app:set_solver(solver)

    --app:set_vtk_scriptor(scriptor)
    --app:set_vtk_ustart_before(VTKScriptor(VTKOutput(),"u_step_before"))
    --app:set_vtk_ustart_after(VTKScriptor(VTKOutput(),"u_step_after"))
    --app:set_vtk_uend_before(VTKScriptor(VTKOutput(),"uend_before"))
    --app:set_vtk_uend_after(VTKScriptor(VTKOutput(),"uend_after"))

    --app:set_vtk_resr_before(VTKScriptor(VTKOutput(),"resr_before"))
    --app:set_vtk_resr_after(VTKScriptor(VTKOutput(),"resr_after"))
    --app:set_vtk_resu_before(VTKScriptor(VTKOutput(),"resu_before"))
    --app:set_vtk_resu_after(VTKScriptor(VTKOutput(),"resu_after"))

    --app:set_vtk_rhs(VTKScriptor(VTKOutput(),"rhs"))
    --app:set_vtk_rhs_res(VTKScriptor(VTKOutput(),"rhs_res"))

    app:set_vtk_residuum(VTKScriptor(VTKOutput(),"residuum_p" .. communicator:get_temporal_rank()))
    app:set_vtk_norm(VTKScriptor(VTKOutput(),"norm_p" .. communicator:get_temporal_rank()))
    app:set_vtk_fstop(VTKScriptor(VTKOutput(),"fstop_p" .. communicator:get_temporal_rank()))

    -- set app specific values
    -- app:set_adapt_convergence()
    -- todo set time integration method

    -- creating executor
    braid = BraidExecutor(communicator, app)
    -- braid:apply()
    -- braid:test()
    -- braid:set_app()
    -- braid:get_app()


    -- braid:set_residual(desc.use_residual)
    braid:set_residual(true)
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
    braid:set_max_levels(desc.max_level) -- todo dublicate?
    braid:set_min_coarse(desc.min_coarsening)
    braid:set_sequential(desc.sequential)
    braid:set_spatial_coarsen_and_refine(desc.spatial_coarsen_and_refine)
    braid:set_refine(desc.time_refinement)
    braid:set_max_refinements(desc.max_refinement)
    braid:set_print_file(desc.printfile)
    -- braid:set_output()
    braid:set_filename(desc.outputfile)
    braid:set_paralog(logging)

    -- braid:set_richardson_estimation(desc.richardson_efstimation, desc.richardson_extrapolation, desc.richardson_local_order)

    app:init()
    braid:print_settings()
    app:print_settings()
    print("XBraid Integrator created")

    return braid
end

function xbraid_util.CreateBraidResidualStepper(
        desc, communicator, logging, solver, scriptor, domain, approxspace)
    --,scaling, adapt, full)
    -- creating app
    app = BraidResidualStepper()
    -- set app base values

    -- cmin, cmax, fmin, fmax, iteration
    --app:set_scaling(scaling,adapt,full)
    app:set_verbose(desc.verbose)

    app:set_start_time(desc.time.t_0)
    app:set_end_time(desc.time.t_end)
    app:set_number_of_timesteps(desc.time.n)
    app:set_time_values(desc.time.t_0, desc.time.t_end, desc.time.n)

    app:set_domain(domain)
    app:set_approx_space(approxspace)
    -- app:set_start_vector()
    -- app:set_scriptor()
    app:set_max_levels(desc.max_level)
    --app:set_fine_time_integrator(fintegrator)
    --app:set_coarse_time_integrator(cintegrator)
    app:set_scriptor(scriptor)
    app:set_solver(solver)

    --app:set_vtk_scriptor(scriptor)
    --app:set_vtk_ustart_before(VTKScriptor(VTKOutput(),"u_step_before"))
    --app:set_vtk_ustart_after(VTKScriptor(VTKOutput(),"u_step_after"))
    --app:set_vtk_uend_before(VTKScriptor(VTKOutput(),"uend_before"))
    --app:set_vtk_uend_after(VTKScriptor(VTKOutput(),"uend_after"))

    --app:set_vtk_resr_before(VTKScriptor(VTKOutput(),"resr_before"))
    --app:set_vtk_resr_after(VTKScriptor(VTKOutput(),"resr_after"))
    --app:set_vtk_resu_before(VTKScriptor(VTKOutput(),"resu_before"))
    --app:set_vtk_resu_after(VTKScriptor(VTKOutput(),"resu_after"))

    --app:set_vtk_rhs(VTKScriptor(VTKOutput(),"rhs"))
    --app:set_vtk_rhs_res(VTKScriptor(VTKOutput(),"rhs_res"))

    --app:set_vtk_norm(VTKScriptor(VTKOutput(),"residuum_p" .. communicator:get_temporal_rank()))

    -- set app specific values
    -- app:set_adapt_convergence()
    -- todo set time integration method

    -- creating executor
    braid = BraidExecutor(communicator, app)
    -- braid:apply()
    -- braid:test()
    -- braid:set_app()
    -- braid:get_app()


    braid:set_residual(true)
    -- braid:set_residual(true)
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
    braid:set_max_levels(desc.max_level) -- todo dublicate?
    braid:set_min_coarse(desc.min_coarsening)
    braid:set_sequential(desc.sequential)
    braid:set_spatial_coarsen_and_refine(desc.spatial_coarsen_and_refine)
    braid:set_refine(desc.time_refinement)
    braid:set_max_refinements(desc.max_refinement)
    braid:set_print_file(desc.printfile)
    -- braid:set_output()
    braid:set_filename(desc.outputfile)
    braid:set_paralog(logging)

    -- braid:set_richardson_estimation(desc.richardson_estimation, desc.richardson_extrapolation, desc.richardson_local_order)

    app:init()
    braid:print_settings()
    app:print_settings()
    print("XBraid Integrator created")

    return braid
end
