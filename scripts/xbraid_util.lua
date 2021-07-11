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
    access_level = 3,
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

function xbraid_util.SetConvCheck(braid,desc)
    braid:set_max_iterations(desc.conv_check.max_iter)
    if (desc.conv_check.absolute ~= nil) then
        braid:set_absolute_tol(desc.conv_check.absolute)
    else
        braid:set_reduction(desc.conv_check.reduction)
    end
end

function xbraid_util.CreateTimeHierarchy(braid,desc)
    if desc.default_cfactor ~= nil then
        braid:set_c_factor(-1,desc.default_cfactor)
    end

    if type(desc.cfactor) == "table" then
        for key,value in pairs(desc.cfactor) do
            braid:set_c_factor(key-1,value)
        end
    else
        braid:set_c_factor(-1,desc.cfactor)
    end
end

function xbraid_util.CreateBraidIntegrator(desc, communicator, logging, fintegrator,cintegrator,scriptor)
    -- creating app
    app = BraidIntegratorFactory()
    -- set app base values


    app:set_verbose(desc.verbose)

    app:set_start_time(desc.time.t_0)
    app:set_end_time(desc.time.t_end)
    app:set_number_of_timesteps(desc.time.n)
    app:set_time_values(desc.time.t_0,desc.time.t_end,desc.time.n)
    app:set_vtk_scriptor(scriptor)
    -- app:set_start_vector()
    -- app:set_scriptor()
    app:set_max_levels(desc.max_level)
    app:set_fine_time_integrator(fintegrator)
    app:set_coarse_time_integrator(cintegrator)
    app:set_scriptor(scriptor)
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

    xbraid_util.set_relax_type(braid,desc.mgrit_relax_type)
    -- braid:set_n_relax()
    xbraid_util.set_cycle_type(braid,desc.mgrit_cycle_type)
    -- braid:set_cycle_fmg()
    -- braid:set_cycle_nfmg()
    -- braid:set_cycle_nfmgv()
    xbraid_util.CreateTimeHierarchy(braid,desc)
    xbraid_util.SetConvCheck(braid,desc)
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

    braid:set_richardson_estimation(desc.richardson_estimation,desc.richardson_extrapolation,desc.richardson_local_order)

    app:init()
    braid:print_settings()
    app:print_settings()
    print("XBraid Integrator created")






    return braid
end


braid_timestepper_desc = {
    time = { t_0 = 0, t_end = 1, n = 100 },
    max_level = 15,
    integrator = impliciteuler, -- or table

    coarsening_factor = 2,
    mgrit_cycle_type = "V",
    mgrit_relax_type = "FCF",
    store_values = 0,
    print_level = 3,
    access_level = 3,
    verbose = true,

    temporalNorm = 3, -- {1,2,3}
    conv_check = {
        iterations = 100,
        -- reduction = 1e-9
        absolute = 5e-7
    }
}