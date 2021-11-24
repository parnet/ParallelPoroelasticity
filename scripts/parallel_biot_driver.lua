walltime = BraidTimer() -- get Time of whole execution
walltime:start()

local math_pi = 3.14159265359
ug_load_script("ug_util.lua")

num_spatial_procs = util.GetParamNumber("--npx", 1, "number of spatial procs (must divide totalproc number)") -- numSpatialProcs * numTimeProcs = numGlobalProcs
num_world_ranks = NumProcs()
space_time_communicator = SpaceTimeCommunicator()

if num_world_ranks % num_spatial_procs == 0 then
    space_time_communicator:split(num_spatial_procs)
    local num_temporal_procs = num_world_ranks / num_spatial_procs;
    local num_spatial_procs = num_spatial_procs
    print("temporal x spatial = world")
    print(num_temporal_procs.." x " .. num_spatial_procs .. " = " .. num_world_ranks)
else
    space_time_communicator:split(1)
    print("temporal x spatial = world")
    print(num_world_ranks.." x " .. 1 .. " = " .. num_world_ranks)
end

ug_load_script("util/load_balancing_util_2.lua")
ug_load_script("util/profiler_util.lua")
ug_load_script("plugins/Limex/limex_util.lua")
ug_load_script("xbraid_util.lua")

environment = util.GetParam("--env", "hawk", "local, hawk, gcsc")
print("ENVIRONMENT: " .. environment)

XARGS = {
    p_method = util.GetParam("--mode", "MGRIT", ""), -- SEQ CHK R NL
    p_redirect = util.GetParamNumber("--redirect", 1, "") == 1, -- 0,1

    p_num_time = util.GetParamNumber("--numtime", 32, " maximum number of levels"),
    p_max_iter = util.GetParamNumber("--maxiter", 100, " maximum number of iterations"),
    p_max_level = util.GetParamNumber("--maxlevel", 15, " maximum number of levels"),
    p_c_factor = util.GetParam("--cfactor", "2_2_2", "relaxation type FCF, FFCF or F-relaxation"),
    p_cycle = util.GetParam("--cycle", "V", " cycletype V-Cycle or F-Cycle "),
    p_relaxation = util.GetParam("--relax", "FCF", "relaxation type FCF, FFCF or F-relaxation"),

    p_driver = util.GetParam("--driver", "Integrator", "relaxation type FCF, FFCF or F-relaxation"),
    p_boolskipdown = util.GetParamNumber("--boolskipdown", 0, "relaxation type FCF, FFCF or F-relaxation") == 1,
    p_useResidual = util.GetParamNumber("--use-residual", 0, " 0 xbraid residual, 1 use residual") == 1,

    p_accesslevel = util.GetParamNumber("--accesslevel", 1, ""),
    p_printlevel = util.GetParamNumber("--printlevel", 1, ""),
    p_store_values = util.GetParamNumber("--store-level", 0, ""),
-- 1e-6, 1e-6 for rel
-- 1e-14, 1e-14 for cmpConvCheck
    p_tol_red_p = util.GetParamNumber("--tol-red-p", 1e-6, " 0 use residual, 1 xbraid residual"),
    p_tol_red_u = util.GetParamNumber("--tol-red-u", 1e-6, " 0 use residual, 1 xbraid residual"),
    p_tol_abs_p = util.GetParamNumber("--tol-abs-p", 1e-14, " 0 use residual, 1 xbraid residual"),
    p_tol_abs_u = util.GetParamNumber("--tol-abs-u", 1e-14, " 0 use residual, 1 xbraid residual"),

    p_tol_abs_braid = util.GetParamNumber("--tol-abs-braid", 1e-50, " 0 use residual, 1 xbraid residual"),
    p_tol_norm_braid = util.GetParamNumber("--tol-norm-braid", 2, " 0 use residual, 1 xbraid residual"),
    p_tol_rel_braid = util.GetParamNumber("--tol-rel-braid", 1e-50, " 0 use residual, 1 xbraid residual"),
}
PARGS = {
    p_napprox = util.GetParamNumber("--napprox", 512, "relaxation type FCF, FFCF or F-relaxation"),
}
IARGS = {
    method = util.GetParam("--integrator", 1, "relaxation type FCF, FFCF or F-relaxation"),
    theta = util.GetParamNumber("--theta", 1, "relaxation type FCF, FFCF or F-relaxation"),
    order = util.GetParamNumber("--order", 2, "relaxation type FCF, FFCF or F-relaxation"),
    num_step = util.GetParamNumber("--gridstep", 2, "relaxation type FCF, FFCF or F-relaxation"),
}
RARGS = {
    rich_est = util.GetParamNumber("--rich-est", 0, "relaxation type FCF, FFCF or F-relaxation") == 1,
    rich_ext = util.GetParamNumber("--rich-ext", 0, "relaxation type FCF, FFCF or F-relaxation") == 1,
    rich_order = util.GetParamNumber("--rich-order", 2, "relaxation type FCF, FFCF or F-relaxation"),
    time_refine = util.GetParamNumber("--trefine", 0, "relaxation type FCF, FFCF or F-relaxation") == 1,

    rich_refine = util.GetParamNumber("--rich-refine", 2, "relaxation type FCF, FFCF or F-relaxation"),
    rich_bound = util.GetParamNumber("--rich-bound", 1.1, "relaxation type FCF, FFCF or F-relaxation"),
    coarse = util.GetParamNumber("--coarse", 0, "relaxation type FCF, FFCF or F-relaxation"),
}
if XARGS.p_redirect then
    repl = ReplaceStandardStream()
    repl:set_space_time_comm(space_time_communicator)
    repl:apply()
end
ug_load_script("test.lua")
util.biot.CheckAssertions()


local dtFrac = util.GetParamNumber("--dtFrac", 1e-5, "time step size")
local dtMinFrac = util.GetParamNumber("--dtminFrac", 1e-2, "minimal admissible time step size")
local dtRed = util.GetParamNumber("--dtred", 0.5, "time step size reduction factor on divergence")
local numRefs = util.GetParamNumber("--num-refs", 3, "total number of refinements (incl. pre-Refinements)")
local paraStab = util.GetParamNumber("--stab", 4, "total number of refinements (incl. pre-Refinements)")
local endTimeFactor = util.GetParamNumber("--endtime", 0, "total number of refinements (incl. pre-Refinements)")

local paraPOrder = util.GetParamNumber("--porder", 1, "total number of refinements (incl. pre-Refinements)")
local paraUOrder = util.GetParamNumber("--uorder", 2, "total number of refinements (incl. pre-Refinements)")

local ARGS = {
    solverID = util.GetParam("--solver-id", "GMG"), --  "FixedStressEX", "UzawaMG", "UzawaSmoother","UzawaMGKrylov"
    smootherID = util.GetParam("--smoother-id", "uzawa"), --  "FixedStressEX", "UzawaMG", "UzawaSmoother","UzawaMGKrylov"

    bSteadyStateMechanics = not util.HasParamOption("--with-transient-mechanics"), -- OPTIONAL: transient mechanics
    MGCycleType = util.GetParam("--mg-cycle-type", "W", "V,F,W"),
    MGBaseLevel = util.GetParamNumber("--mg-base-level", 0, "some non-negative integer"),
    MGNumSmooth = util.GetParamNumber("--mg-num-smooth", 2, "some positive integer"),
    LimexTOL = util.GetParamNumber("--limex-tol", 1e-3, "TOL"),
    LimexNStages = util.GetParamNumber("--limex-num-stages", 4, "number of LIMEX stages q"),
}
print("MGNumSmooth=" .. ARGS.MGNumSmooth)
print("MGCycleType=" .. ARGS.MGCycleType)
print("MGBaseLevel=" .. ARGS.MGBaseLevel)
local problem = BarryMercerProblem2dCPU1("ux,uy", "p")
if num_spatial_procs < 8 then -- todo eliminate this quick fix when point source is fixed
    problem:set_source_scaling(2.0)
elseif num_spatial_procs == 8 then
    problem:set_source_scaling(1.0)
else
    problem:set_source_scaling(0.5)
end

problem:set_stab(paraStab)
problem:set_order(paraUOrder, paraPOrder)
local charTime = problem:get_char_time()  -- implemented by C++ object
print("characteristic time is " .. charTime)

startTime = 0.0 * charTime * math_pi
endTime =   2.0 * charTime * math_pi
if endTimeFactor > 0 then
    endTime = endTimeFactor * math_pi * charTime
end



print("@Integrate from " .. startTime .. " to " .. endTime)
local dt = dtFrac * charTime
local dtMin = dtMinFrac
local dtMax = endTime

local porder = problem:get_porder() or 1
local uorder = problem:get_uorder() or (porder + 1)
print("porder is " .. porder)
print("uorder is " .. uorder)

local dim = 2
local cpu = 1
InitUG(dim, AlgebraType("CPU", cpu));

local balancerDesc = {
    hierarchy = {
        type = "standard",
        minElemsPerProcPerLevel = 32,
        maxRedistProcs = 120,
        qualityRedistLevelOffset = 2,
        intermediateRedistributions = true,
        {
            upperLvl = 0,
            maxRedistProcs = 40
        },
        {
            upperLvl = 2,
            maxRedistProcs = 120
        },
    }
}
if XARGS.p_redirect then
    repl:undo()
    repl:apply()
end
local gridName = problem:get_gridname()
local mandatorySubsets = nil
local dom = util.CreateDomain(gridName, 0, mandatorySubsets)
util.refinement.CreateRegularHierarchy(dom, numRefs, true, balancerDesc)
local approxSpace = util.biot.CreateApproxSpace(dom, dim, uorder, porder)

print("FE discretization...")
local bSteadyStateMechanics = ARGS.bSteadyStateMechanics -- true
local domainDisc0 = DomainDiscretization(approxSpace)
problem:add_elem_discs(domainDisc0, bSteadyStateMechanics)
problem:add_boundary_conditions(domainDisc0, bSteadyStateMechanics)

local domainDiscT = DomainDiscretization(approxSpace)
problem:add_elem_discs(domainDiscT, bSteadyStateMechanics)
problem:add_boundary_conditions(domainDiscT, bSteadyStateMechanics)

local uzawaSchurUpdateDisc = DomainDiscretization(approxSpace)
problem:add_uzawa_discs(uzawaSchurUpdateDisc, bSteadyStateMechanics)
print("Discretization done!")
local u_start = GridFunction(approxSpace)
local gs = GaussSeidel()
local bgs = BackwardGaussSeidel()
bgs:enable_overlap(true)
gs:enable_overlap(true)
local uzawaSchurUpdateOp = AssembledLinearOperator()
uzawaSchurUpdateOp:set_discretization(uzawaSchurUpdateDisc)

function createUzawaIteration(sSchurCmp, aiForward, aiSchur, aiBackward, uzawaSchurUpdateOp, uzawaSchurWeight)
    local uzawa = UzawaBase(sSchurCmp)
    local weight = uzawaSchurWeight or 1.0
    if (aiForward) then
        uzawa:set_forward_iter(aiForward)
    end
    if (aiSchur) then
        uzawa:set_schur_iter(aiSchur)
    end
    if (aiBackward) then
        uzawa:set_backward_iter(aiBackward)
    end
    uzawa:set_schur_operator_update(uzawaSchurUpdateOp, weight)
    return uzawa
end

local uzawaWeight = 1.0

local preSmoother = createUzawaIteration("p", gs, Jacobi(0.66), nil, uzawaSchurUpdateOp, uzawaWeight)
local postSmoother = createUzawaIteration("p", nil, Jacobi(0.66), bgs, uzawaSchurUpdateOp, uzawaWeight)

if ARGS.smootherID == "uzawa3" then
    print("using uzawa 3")
    preSmoother = createUzawaIteration("p", SymmetricGaussSeidel(), SymmetricGaussSeidel(), nil, uzawaSchurUpdateOp, uzawaWeight)
    postSmoother = createUzawaIteration("p", nil, SymmetricGaussSeidel(), SymmetricGaussSeidel(), uzawaSchurUpdateOp, uzawaWeight)
end
local superLU = SuperLU() --LU()
local gmg = GeometricMultiGrid(approxSpace)
gmg:set_discretization(domainDiscT)
gmg:set_base_level(ARGS.MGBaseLevel)  -- was 1 in Cincy
gmg:set_base_solver(superLU)  -- was baseLU in Cincy
gmg:set_presmoother(preSmoother) --(jac)
gmg:set_postsmoother(postSmoother)
gmg:set_cycle_type(ARGS.MGCycleType) -- 1:V, 2:W -- "F"
gmg:set_num_presmooth(ARGS.MGNumSmooth)
gmg:set_num_postsmooth(ARGS.MGNumSmooth)
gmg:set_rap(true)  -- mandatory, if set_stationary
local transfer = StdTransfer()
transfer:enable_p1_lagrange_optimization(true)
gmg:set_transfer(transfer)
local p0 = 1.0
local solver = {}
tol_reduction_p = XARGS.p_tol_red_p
tol_absolute_p = XARGS.p_tol_abs_p
tol_reduction_u = XARGS.p_tol_red_u
tol_absolute_u = XARGS.p_tol_abs_u
-- local cmpConvCheck = CompositeConvCheck(approxSpace)
-- cmpConvCheck:set_component_check("ux", p0 * tol_absolute_u, tol_reduction_u)
-- cmpConvCheck:set_component_check("uy", p0 * tol_absolute_u, tol_reduction_u)
-- if (dim == 3) then
--     cmpConvCheck:set_component_check("uz", p0 * tol_absolute_u, tol_reduction_u)
-- end
-- cmpConvCheck:set_component_check("p", p0 * tol_absolute_p, tol_reduction_p)
-- cmpConvCheck:set_maximum_steps(100)
-- cmpConvCheck:set_verbose(true)
-- cmpConvCheck:set_supress_unsuccessful(false)

local convCheck = ConvCheck()
convCheck:set_maximum_steps(100)
convCheck:set_reduction(1e-8)
convCheck:set_minimum_defect(1e-14)
convCheck:set_verbose(true)
convCheck:set_supress_unsuccessful(false)

local convCheckCoarse = ConvCheck()
convCheckCoarse:set_maximum_steps(RARGS.coarse)
convCheckCoarse:set_reduction(1e-8)
convCheckCoarse:set_minimum_defect(1e-14)
convCheckCoarse:set_verbose(true)
convCheckCoarse:set_supress_unsuccessful(true)

local newtonCheck = ConvCheck()
newtonCheck:set_maximum_steps(1)
newtonCheck:set_reduction(5e-6)
newtonCheck:set_minimum_defect(1e-14)
newtonCheck:set_verbose(true)
newtonCheck:set_supress_unsuccessful(true)

local newtonCheckCoarse = ConvCheck()
newtonCheckCoarse:set_maximum_steps(1)
newtonCheckCoarse:set_reduction(5e-6)
newtonCheckCoarse:set_minimum_defect(1e-14)
newtonCheckCoarse:set_verbose(true)
newtonCheckCoarse:set_supress_unsuccessful(true)

solver["GMG"] = LinearSolver()
solver["GMG"]:set_preconditioner(gmg) -- gmg, dbgIter
solver["GMG"]:set_convergence_check(convCheck)

solver["GMGKrylov"] = BiCGStab()
solver["GMGKrylov"]:set_preconditioner(gmg) -- gmg, dbgIter
solver["GMGKrylov"]:set_convergence_check(convCheck) -- convCheck

solver["LU"] = LinearSolver()
solver["LU"]:set_preconditioner(LU())
solver["LU"]:set_convergence_check(convCheck)
local lsolver = solver[ARGS.solverID]
print("using "..ARGS.solverID)


local uzawaWeight = 1.0
-- coarse solver
local gmgCoarse = GeometricMultiGrid(approxSpace)
gmgCoarse:set_discretization(domainDiscT)
gmgCoarse:set_base_level(ARGS.MGBaseLevel)  -- was 1 in Cincy
gmgCoarse:set_base_solver(superLU)  -- was baseLU in Cincy
gmgCoarse:set_presmoother(preSmoother) --(jac)
gmgCoarse:set_postsmoother(postSmoother)
gmgCoarse:set_cycle_type(ARGS.MGCycleType) -- 1:V, 2:W -- "F"
gmgCoarse:set_num_presmooth(ARGS.MGNumSmooth)
gmgCoarse:set_num_postsmooth(ARGS.MGNumSmooth)
gmgCoarse:set_rap(true)  -- mandatory, if set_stationary
gmgCoarse:set_transfer(transfer)

solverCoarse = {}
-- local cmpConvCheckCoarse = CompositeConvCheck(approxSpace)
-- cmpConvCheckCoarse:set_component_check("ux", p0 * tol_absolute_u, tol_reduction_u)
-- cmpConvCheckCoarse:set_component_check("uy", p0 * tol_absolute_u, tol_reduction_u)
-- if (dim == 3) then
--     cmpConvCheckCoarse:set_component_check("uz", p0 * tol_absolute_u, tol_reduction_u)
-- end
-- cmpConvCheckCoarse:set_component_check("p", p0 * tol_absolute_p, tol_reduction_p)
-- cmpConvCheckCoarse:set_maximum_steps(RARGS.coarse)
-- cmpConvCheckCoarse:set_verbose(true)
-- cmpConvCheckCoarse:set_supress_unsuccessful(false)


solverCoarse["GMG"] = LinearSolver()
solverCoarse["GMG"]:set_preconditioner(gmgCoarse) -- gmg, dbgIter
solverCoarse["GMG"]:set_convergence_check(convCheckCoarse)

solverCoarse["GMGKrylov"] = BiCGStab()
solverCoarse["GMGKrylov"]:set_preconditioner(gmgCoarse) -- gmg, dbgIter
solverCoarse["GMGKrylov"]:set_convergence_check(convCheckCoarse)

local lsolverCoarse = solverCoarse[ARGS.solverID]



local vtk = VTKOutput()

print("initialization done.\n\n\n")
desc_conv_control = {
    type = "static",
    looseTol = 5e-6,
    tightTol = 5e-8,
    force_convergence = false
}
braid_desc = {
    driver = XARGS.p_driver, -- XARGS.p_driver == "IntegratorFactory"
    integrator = "FS",
    time = { t_0 = startTime,
             t_end = endTime,
             n = XARGS.p_num_time },
    cfactor = xbraid_util.get_cfactor(XARGS.p_c_factor),
    default_cfactor = XARGS.p_c_factor_default,
    max_level = XARGS.p_max_level,
    mgrit_cycle_type = XARGS.p_cycle,
    mgrit_relax_type = XARGS.p_relaxation,
    store_values = XARGS.p_store_values,
    print_level = XARGS.p_printlevel,
    access_level = XARGS.p_accesslevel,
    sequential = XARGS.p_method == "SEQMGRIT",
    temporal_norm = XARGS.p_tol_norm_braid, -- {1,2,3}
    conv_check = {
        max_iter = XARGS.p_max_iter,
        -- reduction = XARGS.p_tol_rel_braid
        absolute = XARGS.p_tol_abs_braid
    },
    skip_downcycle_work = XARGS.p_boolskipdown,
    max_refinement = 10,
    spatial_coarsen_and_refine = false,
    min_coarsening = 2,
    printfile = "000 " .. XARGS.p_driver .. "_" .. XARGS.p_num_time .. "_" .. XARGS.p_max_level .. "_" .. XARGS.p_cycle .. "_" .. XARGS.p_relaxation .. ".mgrit",
    outputfile = "integrator_out",
    use_residual = XARGS.p_useResidual,
    sync = RARGS.time_refine and RARGS.rich_est,
    time_refinement = RARGS.time_refine,
    richardson_estimation = RARGS.rich_est, --set_richardson_estimation
    richardson_extrapolation = RARGS.rich_ext,
    richardson_local_order = RARGS.rich_order,
    verbose = true,
}

local newtonSolver = NewtonSolver()
newtonSolver:set_linear_solver(lsolver)
newtonSolver:set_convergence_check(newtonCheck)



local newtonSolverCoarse = NewtonSolver()
newtonSolverCoarse:set_linear_solver(lsolverCoarse)
newtonSolverCoarse:set_convergence_check(newtonCheckCoarse)


local newtonSolverSec = NewtonSolver()
newtonSolverSec:set_linear_solver(lsolver)
newtonSolverSec:set_convergence_check(newtonCheck)

local nlsolver = newtonSolver
local nlsolver_coarse
print("RARGS.coarse " , RARGS.coarse)
if RARGS.coarse == 0 then
    print("##### using fine")
    nlsolver_coarse = nlsolver
else
    print("##### using coarse")
    nlsolver_coarse = newtonSolverCoarse
end

print(lsolver:config_string())

function myStepCallback0(u, step, time)
    -- problem:post_processing(u, step, time)
    -- io = PIOGridFunction()
    -- io:write(u,"solution_t"..step)
    vtk:print("PoroElasticityInitial.vtu", u, step, time)
end
print("Interpolation start values")
problem:interpolate_start_values(u_start, startTime)
print("Integrating from "..startTime.." to " .. endTime)


if ((ARGS.LimexNStages ~= 0)) then
    local dt0 = charTime * 1e-50
    print("Computing consistent initial value w/ dt0=" .. dt0)
    util.SolveNonlinearTimeProblem(u_start, domainDisc0, newtonSolverSec, myStepCallback0, "PoroElasticityInitial", "ImplEuler", 1, startTime, dt0, dt0, dt0, dtRed);
    print("initial value calculation done. \n\n\n\n\n")
end

log_job = Paralog() -- todo move to desc
log_job:set_comm(space_time_communicator)
log_job:set_file_name("joba")
log_job:init()

paralog_script = Paralog() -- todo move to desc
paralog_script:set_comm(space_time_communicator)
paralog_script:set_file_name("script")
paralog_script:init()

scr_cmp = BraidBiotCheckPrecomputed()
scr_cmp:set_log(log_job)
scr_cmp:set_solution_name(vtk, "sequential")
scr_cmp:set_diff_name(vtk, "error")
scr_cmp:set_vtk_write_mode(false,false)
scr_cmp:set_io_write_mode(false,false)
scr_cmp:set_num_ref(numRefs)
scr_cmp:set_max_index(128, braid_desc.time.n)

if environment == "hawk" then
    scr_cmp:set_base_path("/lustre/hpe/ws10/ws10.1/ws/igcmparn-mgrit/analyticsolution")
elseif environment == "gcsc" then
    scr_cmp:set_base_path("/home/mparnet/analyticsolution")
elseif environment == "local" then
    scr_cmp:set_base_path("/home/maro/hawk/analyticsolution")
    --scr_cmp:set_base_path("/home/maro/hawk/analyticsolution")
end

-- scr_biot = BraidBiotCheck()
-- scr_biot:set_problem(problem)
-- scr_biot:set_napprox(PARGS.p_napprox)

scr_vtk = VTKScriptor(vtk, "method")

if (XARGS.p_method == "SEQ") then
    timespan = braid_desc.time.t_end - braid_desc.time.t_0
    dt = timespan / braid_desc.time.n

    uapprox_tstart = u_start:clone()
    uapprox_tstop = u_start:clone()
    local tstop = braid_desc.time.t_0
    local tstart = braid_desc.time.t_0
    print("X\t\t", tstart, " \t ", tstop, " \t ", dt)

    --integrator = xbraid_util.createBDF(domainDiscT,                lsolver, IARGS.order, 1e-8)
    -- integrator = NLThetaIntegrator()
    -- integrator:set_domain(domainDiscT)
    -- integrator:set_solver(nlsolver)
    -- integrator:set_theta(1)
    -- integrator:set_num_steps(IARGS.num_step)
    integrator = xbraid_util.createNLTheta(domainDiscT,nlsolver,IARGS.theta,IARGS.num_step,1e-8)

    print("setup done ")
    time = BraidTimer()
    time:start()


    for i = 1, braid_desc.time.n do
        tstart = tstop
        tstop = tstop + dt
        uapprox_tstart = uapprox_tstop:clone()
        uapprox_tstop = uapprox_tstart:clone()
        integrator:init(uapprox_tstart)
        print("\nSeqStep: ", i, "\t\t from ", tstart, " to ", tstop, "  with dt=", dt)
        integrator:apply(uapprox_tstop, tstop, uapprox_tstart, tstart)
        outputval = uapprox_tstop:clone()

        -- scr_cmp:lua_write(outputval, i, tstop)
        scr_vtk:lua_write(outputval,i,tstop,0,1)
    end
    time:stop()
    integration_time = time:get()
    print("\n"..integration_time, "finished sequential timestepping with integrator")
elseif (XARGS.p_method == "NL") then

    dt = endTime / XARGS.p_num_time
    dtMin = dt
    dtMax = dt


    time = BraidTimer()
    time:start()
    test.SolveNonlinearTimeProblem(u_start, domainDiscT, nlsolver, myStepCallback0, "PoroElasticityTransient",
            "ImplEuler", 1, startTime, endTime, dt, dtMin, 0.5);
    time:stop()
    integration_time = time:get()
    print("\n"..integration_time, "finished sequential timestepping with integrator")

elseif (XARGS.p_method == "CHK") then


    base_path_1024 = "/home/maro/hawk/analyticsolution/num_ref_4/BarryMercer2D_"
    base_path_1448 = "/home/maro/hawk/analyticsolution_check/num_ref_4/BarryMercer2D_"

    time = BraidTimer()
    time:start()
    iogf = IOGridFunction()

    for i = 1, 128 do
        print("================" .. i .."================")
        u = u_start:clone()
        v = u_start:clone()
        print(base_path_1024 ..i ..".gridfunction")
        iogf:read(u, base_path_1024 ..i ..".gridfunction")
        print(base_path_1448 ..i ..".gridfunction")
        iogf:read(v, base_path_1448 ..i ..".gridfunction")
        -- scr_cmp:lua_compare(u,v,i,i/128*charTime*2*math_pi,0,0)
        scr_vtk:lua_write(outputval,i,tstop,0,1)
    end
    time:stop()
    integration_time = time:get()
    print("\n"..integration_time, "finished sequential timestepping with integrator")
elseif (XARGS.p_method == "R") then

    scriptor = BraidBiotCheck()
    scriptor:set_problem(problem)
    scriptor:set_napprox(PARGS.p_napprox)

    local tstop = braid_desc.time.t_end
    local tstart = braid_desc.time.t_0
    offset = 1
    dt_total = tstop - tstart
    t_N = braid_desc.time.n
    dt_fine = dt_total / t_N

    t_rank = space_time_communicator:get_temporal_rank()
    t_rank_total = space_time_communicator:get_temporal_size()
    t_proc = t_N / t_rank_total

    proc_offset = offset + t_proc * t_rank

    print(tstart, "\n", tstop, "\n", t_N, "\n\n", offset, "\n", proc_offset, "\n", dt_total, "\n", dt_fine, "\n\n", t_rank, "\n", t_proc, "\n", proc_offset, "\n\n\n")
    time = BraidTimer()
    for i = proc_offset, proc_offset + t_proc - 1 do
        -- print("get_physical_memory_consumed " .. get_physical_memory_consumed())
        -- print("get_world_memory_consumed " .. get_world_memory_consumed())
        ctime = tstart + i * dt_fine
        print(t_rank, "\t", i, "\t", ctime)
        outputval = u_start:clone()
        scriptor:lua_write(outputval, i, ctime, 0, 0)
    end
    time:stop()
    integration_time = time:get()
    print("\ntimer " .. integration_time)
else
    print("Residual ", braid_desc.use_residual)


    for i = 1, #braid_desc.cfactor do
        scr_cmp:set_c_factor(i - 1, braid_desc.cfactor[i])
    end
    bscriptor = scr_vtk
    if braid_desc.driver == "IntegratorFactory" then
        app = xbraid_util.CreateIntegratorFactory(braid_desc,
                domainDiscT,
                bscriptor,
                fine_integrator,
                coarse_integrator
        )

        braid = xbraid_util.CreateExecutor(braid_desc,
                space_time_communicator,
                app,
                log_job
        )
    elseif braid_desc.driver == "Integrator" then
        app = xbraid_util.CreateIntegrator(braid_desc,
                domainDiscT,
                bscriptor
        )

        xbraid_util.CreateNLLevelFC(app,
                    domainDiscT,
                    nlsolver,
                    nlsolver_coarse,
                    IARGS.theta,
                    IARGS.num_step,
                    1e-8)

        app:set_ref_factor(RARGS.rich_refine)
        app:set_threshold(RARGS.rich_bound)
        braid = xbraid_util.CreateExecutor(braid_desc,
                space_time_communicator,
                app,
                log_job
        )
    elseif braid_desc.driver == "NLIntegrator" then
        print("NLDRIVER")
        app = xbraid_util.CreateNLIntegrator(braid_desc,
                domainDiscT,
                bscriptor
        )
        -- app:set_iter(start_iter, target_iter, fulliter) -- not implemented
        app:set_conv_check(convCheck)
        app:set_tol(1e-3, 1e-14)


        xbraid_util.CreateNLLevelFC(app,
                domainDiscT,
                nlsolver,
                nlsolver_coarse,
                IARGS.theta,
                IARGS.num_step,
                1e-8)

        app:set_ref_factor(RARGS.rich_refine)
        app:set_threshold(RARGS.rich_bound)
        braid = xbraid_util.CreateExecutor(braid_desc,
                space_time_communicator,
                app,
                log_job
        )
    elseif braid_desc.driver == "TimeStepper" then
        app = xbraid_util.CreateTimeStepper(braid_desc,
                domainDiscT,
                bscriptor
        )

        print("Set Stepper Methods - Leveldependend")
        xbraid_util.CreateStepperLevel(app,
                domainDiscT,
                lsolver,
                IARGS.theta, 1e-8)

        braid = xbraid_util.CreateExecutor(braid_desc,
                space_time_communicator,
                app,
                log_job
        )

    elseif XARGS.p_driver == "ResidualStepper" then
        app = xbraid_util.CreateBraidResidualStepper(braid_desc,
                domainDiscT,
                bscriptor,
                lsolver
        )

        braid = xbraid_util.CreateExecutor(braid_desc,
                space_time_communicator,
                app,
                log_job
        )
    else
        print("integrator type not supported " .. XARGS.p_driver)
    end

    braid:set_paralog_script(paralog_script)
    v = u_start:clone()
    sv_init = StartValueInitializer()
    sv_init:set_start_vector(u_start)
    braid:set_initializer(sv_init)
    -- if braid_desc.use_residual then
    --    print("Using euclidian norm")
    --    l2norm = BraidEuclidianNorm()
    --    braid:set_norm_provider(l2norm)
    -- else
    --    print("Using biot norm")
    --    bio_norm = BiotBraidSpatialNorm() --BraidEuclidianNorm()
    --    bio_norm:set_order(4, 2)
    --    bio_norm:set_parameter(1.0, 1000000/7, 250000/7)
    --    braid:set_norm_provider(bio_norm)
    -- end
    norm_displ  = BiotBraidDisplacementNorm()
    norm_displ:set_log(log_job)
    braid:set_norm_provider(norm_displ)
    time = BraidTimer()
    time:start()
    braid:apply(u_start, endTime, u_start, startTime)
    time:stop()

    braid:print_settings()
    braid:print_summary()
    log_job:release()
    integration_time = time:get()
    end
if XARGS.p_redirect then
    repl:undo() -- give back the std::cout to terminal
end

walltime:stop()

if space_time_communicator:get_temporal_rank() == 0 then
    print("\n" .. integration_time .. " seconds for time integration")
    print("\n" .. walltime:get() .. " seconds wall time  (rank=" .. space_time_communicator:get_temporal_rank() .. ")")
end

