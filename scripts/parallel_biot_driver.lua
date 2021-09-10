---
--- Generated by EmmyLua(https://github.com/EmmyLua)
--- Created by parnet.
--- DateTime: 28.05.21 15:38
---
walltime = BraidTimer() -- get Time of whole execution
walltime:start()

-- constant definition
local math_pi = 3.14159265359

-- load lua scripts - global
ug_load_script("ug_util.lua")
ug_load_script("util/load_balancing_util_2.lua")
ug_load_script("util/profiler_util.lua")

-- load lua scripts - plugins
ug_load_script("plugins/Limex/limex_util.lua")

-- load lua scripts - local
ug_load_script("xbraid_util.lua") -- load neccessary XBraid lua interfaces


-- PARALLEL [[
-- XBraid Arguments
XARGS = {
    num_spatial_procs = util.GetParamNumber("--npx", 1, "number of spatial procs (must divide totalproc number)"), -- numSpatialProcs * numTimeProcs = numGlobalProcs
    p_num_time = util.GetParamNumber("--numtime", 32, " maximum number of levels"),
    p_max_iter = util.GetParamNumber("--maxiter", 100, " maximum number of iterations"),

    p_max_level = util.GetParamNumber("--maxlevel", 15, " maximum number of levels"),
    p_c_factor = util.GetParam("--cfactor", "2_2_2", "relaxation type FCF, FFCF or F-relaxation"),
    p_cycle = util.GetParam("--cycle", "V", " cycletype V-Cycle or F-Cycle "),
    p_relaxation = util.GetParam("--relax", "FCF", "relaxation type FCF, FFCF or F-relaxation"),

    p_driver = util.GetParam("--driver", "IntegratorFactory", "relaxation type FCF, FFCF or F-relaxation"),

    pp_skip_downcylce = util.GetParam("--skip", "", "relaxation type FCF, FFCF or F-relaxation"),

    p_useResidual = util.GetParamNumber("--use-residual", 0, " 0 xbraid residual, 1 use residual"),

    p_sequential_exec = util.GetParam("--sequential", "", ""),

    p_tol_reduction = util.GetParamNumber("--tol-reduction", 1e-16, " 0 use residual, 1 xbraid residual"),
    p_tol_absolute = util.GetParamNumber("--tol-absolute", 1e-20, " 0 use residual, 1 xbraid residual"),
}

PARGS = {
    p_napprox = util.GetParamNumber("--napprox", 512, "relaxation type FCF, FFCF or F-relaxation"),
}

IARGS = {
    orderOrTheta = util.GetParamNumber("--orderOrTheta", 1, "relaxation type FCF, FFCF or F-relaxation"),
}

num_world_ranks = NumProcs()
space_time_communicator = SpaceTimeCommunicator()

if num_world_ranks % XARGS.num_spatial_procs == 0 then
    space_time_communicator:split(XARGS.num_spatial_procs)
    num_temporal_procs = num_world_ranks / XARGS.num_spatial_procs;
    num_spatial_procs = XARGS.num_spatial_procs
    print("Using: " .. num_spatial_procs .. " of " .. num_world_ranks .. " for spatial")
else
    space_time_communicator:split(1)
    print("Using: " .. 1 .. " of " .. num_world_ranks .. " for spatial") -- todo exit?
end

repl = ReplaceStandardStream()
repl:set_space_time_comm(space_time_communicator)
repl:apply()

-- PARALLEL ]]

util.biot.CheckAssertions()

-- TIMES AND TIME-STEPPING
-- local startTime = util.GetParamNumber("--start", 0.0, "end time")
-- local endTime = util.GetParamNumber("--end", 1e+5, "end time") - override by characteristic time
local dtFrac = util.GetParamNumber("--dtFrac", 1e-5, "time step size")
local dtMinFrac = util.GetParamNumber("--dtminFrac", 1e-2, "minimal admissible time step size")
-- local dtMaxFrac = util.GetParamNumber("--dtmaxFrac", 0.1, "minimal admissible time step size (as fraction of tend)")
local dtRed = util.GetParamNumber("--dtred", 0.5, "time step size reduction factor on divergence")
-- REFINEMENT
-- local numPreRefs = util.GetParamNumber("--numPreRefs", 0, "number of pre-Refinements (before distributing grid)")
local numRefs = util.GetParamNumber("--num-refs", 3, "total number of refinements (incl. pre-Refinements)") --4 --
local paraStab = util.GetParamNumber("--stab", 4, "total number of refinements (incl. pre-Refinements)") --4 --
local paraPOrder = util.GetParamNumber("--porder", 1, "total number of refinements (incl. pre-Refinements)") --4 --
local paraUOrder = util.GetParamNumber("--uorder", 2, "total number of refinements (incl. pre-Refinements)") --4 --

local ARGS = {
    solverID = util.GetParam("--solver-id", "GMG"), --  "FixedStressEX", "UzawaMG", "UzawaSmoother","UzawaMGKrylov"
    useVTK = util.HasParamOption("--with-vtk", "Plot VTK"),
    useDebugIter = util.HasParamOption("--with-debug-iter", "Activate debug solver."),
    bSteadyStateMechanics = not util.HasParamOption("--with-transient-mechanics"), -- OPTIONAL: transient mechanics
    MGCycleType = util.GetParam("--mg-cycle-type", "W", "V,F,W"),
    MGBaseLevel = util.GetParamNumber("--mg-base-level", 0, "some non-negative integer"),
    MGNumSmooth = util.GetParamNumber("--mg-num-smooth", 2, "some positive integer"),
    MGSmootherType = util.GetParam("--mg-smoother-type", "uzawa3", "uzawa,cgs"),
    MGDebugLevel = util.GetParam("--mg-debug-level", 0, "some non-negative integer"),
    LimexTOL = util.GetParamNumber("--limex-tol", 1e-3, "TOL"),
    LimexNStages = util.GetParamNumber("--limex-num-stages", 4, "number of LIMEX stages q"),
}

print("MGSmootherType=" .. ARGS.MGSmootherType)
print("MGNumSmooth=" .. ARGS.MGNumSmooth)
print("MGCycleType=" .. ARGS.MGCycleType)
print("MGBaseLevel=" .. ARGS.MGBaseLevel)
print("MGDebugLevel=" .. ARGS.MGDebugLevel)

GetLogAssistant():set_debug_level("LIB_DISC_MULTIGRID", ARGS.MGDebugLevel);


local problem = BarryMercerProblem2dCPU1("ux,uy", "p")
problem:set_stab(paraStab)
problem:set_order(paraUOrder, paraPOrder)
if (not problem) then
    print("ERROR: Problem '" .. ARGS.problemID .. "' not found")
    quit()
end

local charTime = problem:get_char_time()  -- implemented by C++ object
print("characteristic time is " .. charTime)

startTime = 0.0
endTime = 2.0 * charTime * math_pi

print("Integrate from " .. startTime .. " to " .. endTime)

local dt = dtFrac * charTime
local dtMin = dtMinFrac
local dtMax = endTime
--local doSteadyState = false
local doTransient = true

----------------------------------
--  Settings
----------------------------------

local dim = 2
local cpu = 1    -- default: block dim+1

-- Order for Ansatz functions.
local porder = problem:get_porder() or 1
local uorder = problem:get_uorder() or (porder + 1)

print("porder is " .. porder)
print("uorder is " .. uorder)

InitUG(dim, AlgebraType("CPU", cpu));

--------------------------------------------------------------------------------
-- Domain / ApproximationSpace setup
--------------------------------------------------------------------------------

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
} -- balancerDesc
repl:undo()
repl:apply() -- reapply
-- Create, Load, Refine and Distribute Domain

local gridName = problem:get_gridname()
local mandatorySubsets = nil
local dom = util.CreateDomain(gridName, 0, mandatorySubsets)
util.refinement.CreateRegularHierarchy(dom, numRefs, true, balancerDesc)

-----------------------------------------------------------------
--  Approximation Space
-----------------------------------------------------------------
local approxSpace = util.biot.CreateApproxSpace(dom, dim, uorder, porder)

--------------------------------------------------------------------------------
-- Problem Setup
--------------------------------------------------------------------------------
print("FE discretization...")
local bSteadyStateMechanics = ARGS.bSteadyStateMechanics -- true

-- For computing consistent initial values.
local domainDisc0 = DomainDiscretization(approxSpace)
problem:add_elem_discs(domainDisc0, bSteadyStateMechanics)  -- implemented by C++ object
problem:add_boundary_conditions(domainDisc0, bSteadyStateMechanics)  -- implemented by C++ object

-- For time-dependent problem.
local domainDiscT = DomainDiscretization(approxSpace)
problem:add_elem_discs(domainDiscT, bSteadyStateMechanics)  -- implemented by C++ object
problem:add_boundary_conditions(domainDiscT, bSteadyStateMechanics)  -- implemented by C++ object

-- For Uzawa fixed-stress smoother.
local uzawaSchurUpdateDisc = DomainDiscretization(approxSpace)
problem:add_uzawa_discs(uzawaSchurUpdateDisc, bSteadyStateMechanics)  -- implemented by C++ object
print("Discretization done!")


--------------------------------------------------------------------------------
--  Algebra
--------------------------------------------------------------------------------
local u_start = GridFunction(approxSpace)
local dbgVector = u_start:clone()
--------------------------------------------------------------------------------


----------------------------------------
-- create algebraic Preconditioner
----------------------------------------

--local jac = Jacobi()
--jac:set_damp(0.66)
local gs = GaussSeidel()
--local sgs = SymmetricGaussSeidel()
local bgs = BackwardGaussSeidel()
bgs:enable_overlap(true)
gs:enable_overlap(true)
--sgs:enable_overlap(true)
--local ilu = ILU()
--ilu:set_beta(-0.5);
--local ilut = ILUT()
--ilut:set_threshold(1e-3)

--local egs_weights = u:clone();
--egs_weights:set(1.0);
--Interpolate(0.1, egs_weights, "p")

--local egs = ElementGaussSeidel() -- patches per node
--egs:select_schur_cmp({ "p" }, 4.0)
--egs:set_relax(0.125)

--local cgs = ComponentGaussSeidel(1.0, { "p" }) -- patches per node
--cgs:set_alpha(1.0)
--cgs:set_beta(1.0)
--- 0 > 0.25  (beta=0.0: no pressure change) -- 1.0: works
--cgs:set_weights(true)

--local ssc_vanka_space
--if (dim == 2) then
--    ssc_vanka_space = VertexCenteredVankaSubspace2dCPU1({ "p" }, { "ux", "uy" })
--else
--    ssc_vanka_space = VertexCenteredVankaSubspace3dCPU1({ "p" }, { "ux", "uy", "uz" })
--end

--local ssc = SequentialSubspaceCorrection(1.0)
--ssc:set_vertex_subspace(ssc_vanka_space)

local dbgWriter = GridFunctionDebugWriter(approxSpace)
local uzawaSchurUpdateOp = AssembledLinearOperator()
uzawaSchurUpdateOp:set_discretization(uzawaSchurUpdateDisc)



--- Factory for Uzawa iteration.
-- @function createUzawaIteration
-- @param #string sSchurCmp  Schur complement will be built for this unknown.
-- @param aiForward Approximate Inverse (forward problem)
-- @param aiSchur Approximate Inverse (Schur complement)
-- @param aiBackward Approximate Inverse (backward problem)
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
    -- uzawa:set_debug(dbgWriter)
    return uzawa
end

--local uzawaForward = {}
--local uzawaBackward = {}

local uzawaWeight = 1.0
--local uzawaForward1 = createUzawaIteration("p", gs, Jacobi(0.66), nil, uzawaSchurUpdateOp, uzawaWeight)
--local uzawaBackward1 = createUzawaIteration("p", nil, Jacobi(0.66), bgs, uzawaSchurUpdateOp, uzawaWeight)

--local uzawaForward2 = createUzawaIteration("p", SymmetricGaussSeidel(), Jacobi(0.66), nil, uzawaSchurUpdateOp, uzawaWeight)
--local uzawaBackward2 = createUzawaIteration("p", nil, Jacobi(0.66), SymmetricGaussSeidel(), uzawaSchurUpdateOp, uzawaWeight)

local uzawaForward_3 = createUzawaIteration("p", SymmetricGaussSeidel(), SymmetricGaussSeidel(), nil, uzawaSchurUpdateOp, uzawaWeight)
local uzawaBackward_3 = createUzawaIteration("p", nil, SymmetricGaussSeidel(), SymmetricGaussSeidel(), uzawaSchurUpdateOp, uzawaWeight)

--uzawaForward[0] = createUzawaIteration("p", Jacobi(0.5), Jacobi(0.5), nil, uzawaSchurUpdateOp, uzawaWeight)
--uzawaBackward[0] = createUzawaIteration("p", nil, Jacobi(0.5), Jacobi(0.5), uzawaSchurUpdateOp, uzawaWeight)

--uzawaForward[4] = createUzawaIteration("p", CG(Jacobi()), CG(Jacobi()), nil, uzawaSchurUpdateOp, uzawaWeight)
--uzawaBackward[4] = createUzawaIteration("p", nil, CG(Jacobi()), CG(Jacobi()), uzawaSchurUpdateOp, uzawaWeight)
--[[
local pi=LinearIteratorProduct()
pi:add_iterator(SymmetricGaussSeidel())
pi:add_iterator(SymmetricGaussSeidel())

local uzawaForward2 = createUzawaIteration("p", pi, SymmetricGaussSeidel(), nil, uzawaSchurUpdateOp, uzawaWeight)
local uzawaBackward2 = createUzawaIteration("p", nil, SymmetricGaussSeidel(), pi, uzawaSchurUpdateOp, uzawaWeight)
--]]

--local uzawaForward = createUzawaIteration("p", Jacobi(0.66), Jacobi(0.66), nil, uzawaSchurUpdateOp, uzawaWeight)
--local uzawaBackward = createUzawaIteration("p", nil, Jacobi(0.66), Jacobi(0.66), uzawaSchurUpdateOp, uzawaWeight)
--local uzawaSym = createUzawaIteration("p", gs, sgs, bgs, uzawaSchurUpdateOp, uzawaWeight)
--local uzawaBackward = createUzawaIteration("p", nil, Jacobi(0.5), Jacobi(0.66), uzawaSchurOp, uzawaWeight)
--local uzawa = uzawaForward

local preSmoother
local postSmoother

if (ARGS.MGSmootherType == "uzawa3") then
    preSmoother = uzawaForward_3
    postSmoother = uzawaBackward_3
    --elseif (ARGS.MGSmootherType == "uzawa0") then
    --    preSmoother = uzawaForward[0]
    --    postSmoother = uzawaBackward[0]
    --elseif (ARGS.MGSmootherType == "uzawa2") then
    --    preSmoother = uzawaForward2
    --    postSmoother = uzawaBackward2
    --elseif (ARGS.MGSmootherType == "uzawa") then
    -- preSmoother = uzawaForward1
    -- postSmoother = uzawaBackward1
    --elseif (ARGS.MGSmootherType == "uzawa4") then
    --    preSmoother = uzawaForward[4]
    --    postSmoother = uzawaBackward[4]
    --elseif (ARGS.MGSmootherType == "cgs") then
    --    preSmoother = cgs
    --    postSmoother = cgs
    --elseif (ARGS.MGSmootherType == "vanka-ssc") then
    --    preSmoother = ssc
    --    postSmoother = ssc
    --elseif (ARGS.MGSmootherType == "sgs") then
    --    preSmoother = sgs
    --    postSmoother = sgs
else
    quit()
end

-------------------------
-- create GMG
-------------------------

-- Base Solver
local baseConvCheck = ConvCheck()
baseConvCheck:set_maximum_steps(5000)
baseConvCheck:set_reduction(1e-12)
baseConvCheck:set_verbose(false)

--local base = BiCGStab()
--base:set_preconditioner(jac)
--base:set_convergence_check(baseConvCheck)

--local baseCG = CG()
--baseCG:set_preconditioner(jac)
--baseCG:set_convergence_check(baseConvCheck)

-- exact base solver
--local baseLU = LU()
local superLU = SuperLU()


-- Geometric Multi Grid
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


-- gmg:set_debug(dbgWriter)
-- if (problem.bRAP) then print ("gmg:bRAP=true")
--else print ("gmg:bRAP=false") end
-- gmg:set_debug(dbgWriter)



local transfer = StdTransfer()
transfer:enable_p1_lagrange_optimization(true)
--transfer:set_debug(dbgWriter)
gmg:set_transfer(transfer)

--local gmgP = GeometricMultiGrid(approxSpace)
-- gmgP:set_discretization(domainDiscP)
--gmgP:set_base_level(numPreRefs)  -- was 1 in Cincyj
--gmgP:set_base_solver(baseLU)  -- was baseLU in Cincy
--gmgP:set_presmoother(sgs)
--gmgP:set_postsmoother(sgs)
---gmgP:set_cycle_type("V") -- 1:V, 2:W -- "F"
--gmgP:set_num_presmooth(3)
----gmgP:set_num_postsmooth(3)
--gmgP:set_rap(true)  -- mandatory, if set_stationary


--local gmgU = GeometricMultiGrid(approxSpace)
-- gmgU:set_discretization(domainDiscU)
--gmgU:set_base_level(numPreRefs)  -- was 1 in Cincyj
--gmgU:set_base_solver(baseLU)  -- was baseLU in Cincy
--gmgU:set_presmoother(sgs)
--gmgU:set_postsmoother(sgs)
--gmgU:set_cycle_type("V") -- 1:V, 2:W -- "F"
--gmgU:set_num_presmooth(3)
--gmgU:set_num_postsmooth(3)
--gmgU:set_rap(true)  -- mandatory, if add_ionary
--gmgU:set_debug(dbgWriter)

--local uzawaTotal = createUzawaIteration("p", ILUT(1e-8), ILUT(1e-8), nil, uzawaSchurUpdateOp, 1.0)      -- ???
--local fixedStressLU = createUzawaIteration("p", nil, ILUT(1e-12), ILUT(1e-12), uzawaSchurUpdateOp, 1.0)
--local fixedStressSuperLU = createUzawaIteration("p", nil, SuperLU(), SuperLU(), uzawaSchurUpdateOp, 1.0)
-- local fixedStressLU    = createUzawaIteration("p", nil, SuperLU(), SuperLU(), uzawaSchurUpdateOp, 1.0)
-- local fixedStressMG = createUzawaIteration("p", nil, gmgU, gmgP, uzawaSchurUpdateOp, 1.0)


--local transfer = StdTransfer()
--transfer:enable_p1_lagrange_optimization(false)
--gmg:set_transfer(transfer)

--------------------------------
-- debug solver /iter
--------------------------------


local p0 = 1.0



--------------------------------
-- create and choose a Solver
--------------------------------
local solver = {}


-- GMG

-- tol_reduction = 1e-16
-- tol_absolute = 1e-22
tol_reduction = XARGS.p_tol_reduction
tol_absolute = XARGS.p_tol_absolute

local cmpConvCheck = CompositeConvCheck(approxSpace)
cmpConvCheck:set_component_check("ux", p0 * tol_absolute, tol_reduction)
cmpConvCheck:set_component_check("uy", p0 * tol_absolute, tol_reduction)
if (dim == 3) then
    cmpConvCheck:set_component_check("uz", p0 * tol_absolute, tol_reduction)
end
cmpConvCheck:set_component_check("p", p0 * tol_absolute, tol_reduction)
cmpConvCheck:set_maximum_steps(100)
cmpConvCheck:set_verbose(false)

solver["GMG"] = LinearSolver()
solver["GMG"]:set_preconditioner(gmg) -- gmg, dbgIter
solver["GMG"]:set_convergence_check(cmpConvCheck)


-- LU
local convCheck = ConvCheck()
convCheck:set_maximum_steps(100)
convCheck:set_reduction(1e-8)
convCheck:set_minimum_defect(1e-14)
convCheck:set_verbose(false)

solver["LU"] = LinearSolver()
solver["LU"]:set_preconditioner(LU())
solver["LU"]:set_convergence_check(convCheck)

local lsolver = solver[ARGS.solverID]
local vtk = VTKOutput()

-- Init error estimator.
local biotErrorEst
biotErrorEst = util.biot.CreateDefaultErrorEst(dim)

if (problem.error_estimator) then
    biotErrorEst = problem:error_estimator()
else
    biotErrorEst = util.biot.CreateDefaultErrorEst(dim)
end

--------------------------------------------------------------------------------
--  Solve transient (linear) problem
--------------------------------------------------------------------------------
print("initialization done.\n\n\n")

-- adaptive convergence
desc_conv_control = {
    type = "static",
    looseTol = 5e-6,
    tightTol = 5e-8,
    force_convergence = false
}

local boolskipdown = true
if XARGS.pp_skip_downcylce == "NO" then
    boolskipdown = false
end

-- PARALLEL [[
braid_desc = {
    type = "integrator",
    time = { t_0 = startTime, t_end = endTime, n = XARGS.p_num_time }, --math.ceil((endTime-startTime)/dt) },
    cfactor = xbraid_util.get_cfactor(XARGS.p_c_factor), --{XARGS.p_c_factor,2,2,2,2}, -- 0 finest level,
    --cfactor = 2,
    default_cfactor = XARGS.p_c_factor_default,
    max_level = XARGS.p_max_level,

    integrator = limex, -- todo or table

    mgrit_cycle_type = XARGS.p_cycle,
    mgrit_relax_type = XARGS.p_relaxation,
    store_values = 0,
    print_level = 3,
    access_level = 1,

    sequential = false, -- todo change for parallel

    temporal_norm = 3, -- {1,2,3}
    conv_check = {
        max_iter = XARGS.p_max_iter,
        -- reduction = 1e-9
        absolute = 5e-7
    },

    skip_downcycle_work = boolskipdown,
    time_refinement = false,
    max_refinement = 10,
    spatial_coarsen_and_refine = false,
    min_coarsening = 2,

    printfile = "000 " .. XARGS.p_driver .. "_" .. XARGS.p_num_time .. "_" .. XARGS.p_max_level .. "_" .. XARGS.p_cycle .. "_" .. XARGS.p_relaxation .. ".mgrit",
    outputfile = "integrator_out",
    -- output = Scriptor or multiscriptor if table
    -- store_operator
    verbose = true,
    use_residual = false,

    richardson_estimation = false, --set_richardson_estimation
    richardson_extrapolation = false,
    richardson_local_order = 2,
}


vtk_scriptor = VTKScriptor(vtk, "access")
-- PARALLEL ]]

local newtonCheck = ConvCheck()
newtonCheck:set_maximum_steps(10)
newtonCheck:set_minimum_defect(1e-14)
newtonCheck:set_reduction(5e-6)
newtonCheck:set_verbose(true)

local newtonCheck2 = CompositeConvCheck(approxSpace)
newtonCheck2:set_component_check("ux", p0 * 1e-7, 5e-8)
newtonCheck2:set_component_check("uy", p0 * 1e-7, 5e-8)
if (dim == 3) then
    newtonCheck2:set_component_check("uz", p0 * 1e-7, 5e-8)
end
newtonCheck2:set_component_check("p", p0 * 1e-9, 5e-8)
newtonCheck2:set_maximum_steps(2)

local newtonSolver = NewtonSolver()
newtonSolver:set_linear_solver(lsolver)
newtonSolver:set_convergence_check(newtonCheck)

local nlsolver = newtonSolver
print(lsolver:config_string())

if (doTransient) then
    print("Interpolation start values")
    problem:interpolate_start_values(u_start, startTime)

    print("Integrating from 0.0 to " .. endTime)

    local charTime = problem:get_char_time()
    dt = 1e-2 * charTime
    dtMin = 1e-2 * 1e-2 * charTime

    --
    local myclock = CuckooClock()
    local stepClock = CuckooClock()

    if ((ARGS.LimexNStages > 0)) then
        local dt0 = charTime * 1e-50
        print("Computing consistent initial value w/ dt0=" .. dt0)
        util.SolveNonlinearTimeProblem(u_start, domainDisc0, nlsolver, myStepCallback0, "PoroElasticityInitial", "ImplEuler", 1, startTime, dt0, dt0, dt0, dtRed);
        print("initial value calculation done. \n\n\n\n\n")
    end

    if (XARGS.p_sequential_exec == "X") then
        vtk_scriptor = VTKScriptor(vtk, "output")
        timespan = braid_desc.time.t_end - braid_desc.time.t_0
        dt = timespan / braid_desc.time.n

        uapprox_tstart = u_start:clone()
        uapprox_tstop = u_start:clone()
        local tstop = braid_desc.time.t_0
        local tstart = braid_desc.time.t_0
        print("X\t\t", tstart, " \t ", tstop, " \t ", dt)
        integrator = xbraid_util.creadFSTheta(domainDiscT,
                lsolver, 1, 2, 1e-8)

        time = BraidTimer()
        time:start()

        -- outputval = uapprox_tstop:clone()
        -- vtk_scriptor:lua_write(outputval, 0, tstop, 0, 0)

        for i = 1, braid_desc.time.n do
            tstart = tstop
            tstop = tstop + dt
            uapprox_tstart = uapprox_tstop:clone()
            uapprox_tstop = uapprox_tstart:clone()
            integrator:init(uapprox_tstart)
            integrator:prepare(uapprox_tstart)
            print("SeqStep: ", i, "\t\t from ", tstart, " to ", tstop, "  with dt=", dt)
            integrator:apply(uapprox_tstop, tstop, uapprox_tstart, tstart)

            -- outputval = uapprox_tstop:clone()
            -- scriptor:lua_write(outputval,i,tstop,0,0)

            -- outputval = uapprox_tstop:clone()
            -- vtk_scriptor:lua_write(outputval, i, tstop, 0, 0)
        end
        time:stop()
        integration_time = time:get()
        print(integration_time, "finished sequential timestepping with integrator")
    elseif (XARGS.p_sequential_exec == "R") then

        scriptor = BraidBiotCheck()
        scriptor:set_problem(problem)
        scriptor:set_napprox(PARGS.p_napprox)

        local tstop = braid_desc.time.t_end
        local tstart = braid_desc.time.t_0
        offset = 1
        dt_total = tstop - tstart
        t_N = braid_desc.time.n
        dt_fine = dt_total / t_N

        print(space_time_communicator:get_temporal_rank())

        t_rank = space_time_communicator:get_temporal_rank()
        t_rank_total = space_time_communicator:get_temporal_size()
        t_proc = t_N / t_rank_total

        proc_offset = offset + t_proc * t_rank

        print(tstart, "\n", tstop, "\n", t_N, "\n\n", offset, "\n", proc_offset, "\n", dt_total, "\n", dt_fine, "\n\n", t_rank, "\n", t_proc, "\n", proc_offset, "\n\n\n")
        time = BraidTimer()
        -- for i = proc_offset, proc_offset+t_proc-1 do
        i = 128
        ctime = tstart + i * dt_fine
        print(t_rank, "\t", i, "\t", ctime)
        outputval = u_start:clone()
        scriptor:lua_write(outputval, i, ctime, 0, 0)
        -- end
        time:stop()
        integration_time = time:get()
        print("timer " .. integration_time)
    elseif (ARGS.LimexNStages == 0) then
        print("Solving predefined step sizes (TESTING) ")
        -- Execute linear solver test suite.
        convCheck:set_reduction(1e-10)
        convCheck:set_maximum_steps(100)

        local dtTestSet = { 1.0, 0.1, 0.01, 1e-3, 1e-4, 1e-6, 1e-8, 0.0 }
        for index, dtvalue in ipairs(dtTestSet) do
            dt = dtvalue * charTime
            endTime = dt
            print("%DTFACTOR=\t" .. dtvalue .. "\tindex=\t" .. index)
            problem:interpolate_start_values(u_start, startTime)
            myclock:tic()
            util.SolveNonlinearTimeProblem(u_start, domainDiscT, nlsolver, myStepCallback0, "SolverTest" .. index,
                    "ImplEuler", 1, startTime, endTime, dt, dt, dtRed);
            print("MYCLOCK=" .. myclock:cuckoo() .. "; " .. myclock:toc())
        end


    elseif (ARGS.LimexNStages == 1) then
        print("Linear time Stepping")
        function myStepCallback0(u, step, time)
            -- problem:post_processing(u, step, time)
            vtk:print("PoroElasticity.vtu", u, step, time)
        end
        time = BraidTimer()
        time:start()
        time:stop()
        -- STANDARD (implicit Euler) time-stepping.
        local bCheckpointing = false

        local tstop = braid_desc.time.t_end
        local tstart = braid_desc.time.t_0
        dt_total = tstop - tstart
        t_N = braid_desc.time.n
        dt = dt_total / t_N
        dtMin = dt / 2

        util.SolveLinearTimeProblem(u_start, domainDiscT, lsolver, nil, "PoroElasticityTransient",
                "ImplEuler", 1, startTime, endTime, dt, dtMin, 0.5,
                bCheckpointing, myStepCallback0);

        --util.SolveNonlinearTimeProblem(u, domainDiscT, nlsolver, myStepCallback0, "PoroElasticityTransient",
        --             "ImplEuler", 1, startTime, endTime, dt, dtMin, 0.5);
        integration_time = time:get()
    else
        print("Using Limex for Time Stepping")
        newtonCheck:set_maximum_steps(1)
        newtonCheck:set_supress_unsuccessful(true)
        -- Create & configure LIMEX descriptor.
        local limexDesc = {
            nstages = ARGS.LimexNStages,
            steps = { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 },

            domainDisc = domainDiscT,
            nonlinSolver = nlsolver,

            tol = ARGS.LimexTOL,
            dt = dt,
            dtmin = dtMin,
            dtmax = dtMax,
        }

        ARGS.useVTK = true -- todo adapt
        if (ARGS.useVTK) then
            local vtkFull = VTKOutput()
            local vtkobserver = VTKOutputObserver("PoroElasticityLimex.vtk", vtkFull)
            --limex:attach_observer(vtkobserver)
        end

        local luaobserver = LuaCallbackObserver()

        function myLuaLimexPostProcess(step, time, currdt)
            print("Time per step :" .. stepClock:toc()) -- get time for last step
            local usol = luaobserver:get_current_solution()
            problem:post_processing(usol, step, time) -- todo uncommment to use
            stepClock:tic() -- reset timing
            return 0;
        end

        luaobserver:set_callback("myLuaLimexPostProcess")
        --limex:attach_observer(luaobserver)
        -- Solve problem using LIMEX.
        -- myclock:tic()
        -- limex:apply(u, endTime, u, startTime)
        -- print("CDELTA=" .. myclock:toc())
        print("Residual ", braid_desc.use_residual)

        logging = Paralog() -- todo move to desc
        logging:set_comm(space_time_communicator)
        logging:set_file_name("joba")

        vtk_scriptor = VTKScriptor(vtk, "access")
        if XARGS.p_driver == "IntegratorFactory" then
            app = xbraid_util.CreateIntegratorFactory(braid_desc,
                    domainDiscT,
                    vtk_scriptor,
                    fine_integrator,
                    coarse_integrator
            )

            braid = xbraid_util.CreateExecutor(braid_desc,
                    space_time_communicator,
                    app,
                    logging
            )

        elseif XARGS.p_driver == "Integrator" then
            print("Create Integrator")
            app = xbraid_util.CreateIntegrator(braid_desc,
                    domainDiscT,
                    vtk_scriptor
            )
            print("Set Integrator Methods - Default")
            app:set_default_integrator(xbraid_util.creadFSTheta(domainDiscT,
                    lsolver, 1, 1, 1e-8))
            print("Set Integrator Methods - Leveldependend")
            app:set_integrator(0, xbraid_util.creadFSTheta(domainDiscT,
                    lsolver, 1, 1, 1e-8))

            app:set_integrator(1, xbraid_util.creadFSTheta(domainDiscT,
                    lsolver, 1, 1, 1e-8))

            app:set_integrator(2, xbraid_util.creadFSTheta(domainDiscT,
                    lsolver, 1, 1, 1e-8))

            app:set_integrator(3, xbraid_util.creadFSTheta(domainDiscT,
                    lsolver, 1, 1, 1e-8))

            app:set_integrator(4, xbraid_util.creadFSTheta(domainDiscT,
                    lsolver, 1, 1, 1e-8))

            app:set_integrator(5, xbraid_util.creadFSTheta(domainDiscT,
                    lsolver, 1, 1, 1e-8))

            app:set_integrator(6, xbraid_util.creadFSTheta(domainDiscT,
                    lsolver, 1,1, 1e-8))

            app:set_integrator(7, xbraid_util.creadFSTheta(domainDiscT,
                    lsolver, 1, 1, 1e-8))

            app:set_integrator(8, xbraid_util.creadFSTheta(domainDiscT,
                    lsolver, 1, 1, 1e-8))

            app:set_integrator(9, xbraid_util.creadFSTheta(domainDiscT,
                    lsolver, 1, 1, 1e-8))

            app:set_integrator(10, xbraid_util.creadFSTheta(domainDiscT,
                    lsolver, 1, 1, 1e-8))

            print("Create Braid Object")
            braid = xbraid_util.CreateExecutor(braid_desc,
                    space_time_communicator,
                    app,
                    logging
            )
            print("Finished")
        elseif XARGS.p_driver == "TimeStepper" then
            app = xbraid_util.CreateTimeStepper(braid_desc,
                    domainDiscT,
                    vtk_scriptor,
                    lsolver
            )

            braid = xbraid_util.CreateExecutor(braid_desc,
                    space_time_communicator,
                    app,
                    logging
            )

        elseif XARGS.p_driver == "ResidualStepper" then
            app = xbraid_util.CreateBraidResidualStepper(braid_desc,
                    domainDiscT,
                    vtk_scriptor,
                    lsolver
            )

            braid = xbraid_util.CreateExecutor(braid_desc,
                    space_time_communicator,
                    app,
                    logging
            )
        else
            print("integrator type not supported " .. XARGS.p_driver)
        end

        script_logging = Paralog() -- todo move to desc
        script_logging:set_comm(space_time_communicator)
        script_logging:set_file_name("script")
        script_logging:init()
        braid:set_paralog_script(script_logging)

        v = u_start:clone()

        sv_init = StartValueInitializer()
        sv_init:set_start_vector(u_start)
        braid:set_initializer(sv_init)

        if braid_desc.p_useResidual then
            print("Using euclidian norm")
            l2norm = BraidEuclidianNorm()

            braid:set_norm_provider(l2norm)
        else
            print("Using biot norm")
            bio_norm = BiotBraidSpatialNorm() --BraidEuclidianNorm()
            bio_norm:set_order(4,2)
            bio_norm:set_parameter(1.0, 142857,35714.3)

            braid:set_norm_provider(bio_norm)
        end




        time = BraidTimer()
        time:start()
        --print("starttttt ")
        --space_time_communicator:sleep(100000000)


        braid:apply(u_start, endTime, u_start, startTime)

        time:stop()
        braid:print_settings()
        braid:print_summary()
        logging:release()

        integration_time = time:get()
    end
end -- doTransient
repl:undo()
-- PARALLEL [[
walltime:stop()

if space_time_communicator:get_temporal_rank() == 0 then
    print(integration_time .. " seconds for time integration")
    print(walltime:get() .. " seconds wall time  (rank=" .. space_time_communicator:get_temporal_rank() .. ")")
end

-- PARALLEL ]]

