walltime = BraidTimer() -- get Time of whole execution
walltime:start()

local math_pi = 3.14159265359

ug_load_script("ug_util.lua")
num_world_ranks = NumProcs()
space_time_communicator = SpaceTimeCommunicator()
num_spatial_procs = util.GetParamNumber("--npx", 1, "number of spatial procs (must divide totalproc number)") -- numSpatialProcs * numTimeProcs = numGlobalProcs
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

    p_tol_red_p = util.GetParamNumber("--tol-red-p", 1e-6, " 0 use residual, 1 xbraid residual"),
    p_tol_red_u = util.GetParamNumber("--tol-red-u", 1e-6, " 0 use residual, 1 xbraid residual"),
    p_tol_abs_p = util.GetParamNumber("--tol-abs-p", 1e-14, " 0 use residual, 1 xbraid residual"),
    p_tol_abs_u = util.GetParamNumber("--tol-abs-u", 1e-14, " 0 use residual, 1 xbraid residual"),
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
}
if XARGS.p_redirect then
    repl = ReplaceStandardStream()
    repl:set_space_time_comm(space_time_communicator)
    repl:apply()
end
ug_load_script("test.lua")
util.biot.CheckAssertions()


-- TIMES AND TIME-STEPPING
local startTime  = util.GetParamNumber("--start", 0.0, "end time")
local endTime    = util.GetParamNumber("--end", 1e+5, "end time")
local dtFrac     = util.GetParamNumber("--dtFrac", 1e-5, "time step size")
local dtMinFrac  = util.GetParamNumber("--dtminFrac", 1e-2, "minimal admissible time step size")
local dtMaxFrac  = util.GetParamNumber("--dtmaxFrac", 0.1, "minimal admissible time step size (as fraction of tend)")
local dtRed      = util.GetParamNumber("--dtred", 0.5, "time step size reduction factor on divergence")


-- REFINEMENT
local numPreRefs   = util.GetParamNumber("--numPreRefs", 0, "number of pre-Refinements (before distributing grid)")
local numRefs      = util.GetParamNumber("--num-refs", 3, "total number of refinements (incl. pre-Refinements)") --4 --



local ARGS = {
    problemID = util.GetParam("--problem-id", "deleeuw2d"), -- cryer3d‚
    solverID =  util.GetParam("--solver-id", "GMG"),  --  "FixedStressEX", "UzawaMG", "UzawaSmoother","UzawaMGKrylov"

    useVTK =  util.HasParamOption("--with-vtk", "Plot VTK"),
    useDebugIter =  util.HasParamOption("--with-debug-iter", "Activate debug solver."),
    -- doCheck =  util.HasParamOption("--with-check", ""),

    bSteadyStateMechanics = not util.HasParamOption("--with-transient-mechanics"), -- OPTIONAL: transient mechanics

    MGCycleType = util.GetParam("--mg-cycle-type", "W", "V,F,W"),
    MGBaseLevel = util.GetParamNumber("--mg-base-level", 0, "some non-negative integer"),
    MGNumSmooth = util.GetParamNumber("--mg-num-smooth", 2, "some positive integer"),
    MGSmootherType =  util.GetParam("--mg-smoother-type", "uzawa", "uzawa,cgs"),
    MGDebugLevel =  util.GetParam("--mg-debug-level", 0, "some non-negative integer"),

    -- LIMEX
    LimexTOL     = util.GetParamNumber("--limex-tol", 1e-3, "TOL"),
    LimexNStages = util.GetParamNumber("--limex-num-stages", 4, "number of LIMEX stages q"),

}



print ("MGSmootherType="..ARGS.MGSmootherType)
print ("MGNumSmooth="..ARGS.MGNumSmooth)
print ("MGCycleType="..ARGS.MGCycleType)
print ("MGBaseLevel="..ARGS.MGBaseLevel)
print ("MGDebugLevel="..ARGS.MGDebugLevel)

GetLogAssistant():set_debug_level("LIB_DISC_MULTIGRID", ARGS.MGDebugLevel);
--SetDebugLevel("LIB_DISC_MULTIGRID", 0)

-- Set parameters
local kperm   = 1e-0 -- m/s 1e-3
local poro    = 0.2
local nu      = 0.25

local EYoung  = 2.0 * 1e+2                  -- kPa 2.0 * 1e+4
local Kmedium = EYoung/(3.0*(1.0-2.0*nu))
local Kfluid  = 2.2 * 1e+6                  -- kPa -- 2.2 * 1e+6 --

print ("Kmedium = "..Kmedium)
print ("Kfluid  = "..Kfluid)

--deleeuw2d--deleeuw2d -- cryer3d --cryer2d -- mandel3d --, mandel--, cryer3d
local problemList = {
    -- legacy style
    ["deleeuw2d"] = deleeuw2d,
    ["deleeuw3d"] = deleeuw3d,
    ["deleeuw3dTet"] = deleeuw3dTet,
    ["cryer3d"] = cryer3d,
    ["cryer3dTet"] = cryer3dTet,
    ["footing2D"] = footing2D,
    ["footing2D_tri"] = footing2D_tri,
    ["footing3D"] = footing3D,

    ["bm2D_tri"] = barrymercer2D_tri,

    -- new style
    ["bm2D_new"] = BarryMercerProblem2dCPU1("ux,uy", "p"),
}

local problem = problemList[ARGS.problemID]
if (not problem) then
    print ("ERROR: Problem '".. ARGS.problemID.. "' not found")
    quit()
end

if (problem.parse_cmd_args) then
    problem:parse_cmd_args()
end

--problem:init(kperm, poro, nu, 1.0/Kmedium, 1.0/Kfluid, 0.0)
if(problem.init) then
    problem:init(kperm, poro, nu, 1.0/Kmedium, 0.0, 0.0)
end

print("Get charTime")
local charTime = problem:get_char_time()  -- implemented by C++ object
print("charTime="..charTime)

startTime = 0.0
endTime   = 2.0*charTime

local dt  = dtFrac*charTime
local dtMin = dtMinFrac
local dtMax = endTime


if (problem == mandel) then
    local time = 1e-4
    while (time<=10.0) do
        problem:create_test_data(time*charTime)
        time = time *10.0;
    end
end


local doSteadyState = false
local doTransient = true

----------------------------------
----------------------------------
--  Settings
----------------------------------
----------------------------------

local dim = problem.dim or 2
local cpu = problem.cpu or 1    -- default: block dim+1

-- Order for Ansatz functions.
local porder = problem.porder or problem:get_porder() or 1
local uorder = problem.uorder or problem:get_uorder() or (porder+1)

InitUG(dim, AlgebraType("CPU", cpu));



-- OUTPUT-ASSISTANT FOR SEVERAL PROCESSES
GetLogAssistant():enable_file_output(true, "output_p_"..ProcRank()..
        "_Lev"..numRefs..".txt")
GetLogAssistant():set_debug_level("SchurDebug", 7);
--------------------------------------------------------------------------------
--------------------------------------------------------------------------------
-- Domain / ApproximationSpace setup
--------------------------------------------------------------------------------
--------------------------------------------------------------------------------


local balancerDesc = {
    hierarchy = {
        --type = "standard",
        -- maxRedistProcs = ARGS.redistProcs,

        -- minElemsPerProcPerLevel = ARGS.minElemsPerProcPerLevel,
        -- qualityRedistLevelOffset = ARGS.qualityRedistLevelOffset,
        -- intermediateRedistributions = ARGS.intermediateRedistributions,

        type            = "standard",
        minElemsPerProcPerLevel   = 32,
        maxRedistProcs        = 120,
        qualityRedistLevelOffset  = 2,
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


-- Create, Load, Refine and Distribute Domain
local gridName = problem.gridName or problem:get_gridname()
-- local dom = problem:create_domain(numRefs, numPreRefs)

local mandatorySubsets = problem.mandatorySubsets or nil
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
print("done!")


--------------------------------------------------------------------------------
--------------------------------------------------------------------------------
--  Algebra
--------------------------------------------------------------------------------
local u = GridFunction(approxSpace)
local dbgVector=u:clone()
--------------------------------------------------------------------------------


----------------------------------------
-- create algebraic Preconditioner
----------------------------------------

local jac = Jacobi()
jac:set_damp(0.66)
local gs = GaussSeidel()
local sgs = SymmetricGaussSeidel()
local bgs = BackwardGaussSeidel()
bgs:enable_overlap(true)
gs:enable_overlap(true)
sgs:enable_overlap(true)
local ilu = ILU()
--ilu:set_beta(-0.5);
local ilut = ILUT()
ilut:set_threshold(1e-3)

--local egs_weights = u:clone();
--egs_weights:set(1.0);
--Interpolate(0.1, egs_weights, "p")

local egs = ElementGaussSeidel() -- patches per node
egs:select_schur_cmp({"p"}, 4.0)
egs:set_relax(0.125)

local cgs = ComponentGaussSeidel(1.0, {"p"}) -- patches per node
cgs:set_alpha(1.0)
cgs:set_beta(1.0) --- 0 > 0.25  (beta=0.0: no pressure change) -- 1.0: works
cgs:set_weights(true)


local ssc_vanka_space
if (dim == 2) then
    ssc_vanka_space = VertexCenteredVankaSubspace2dCPU1({"p"}, {"ux", "uy"})
else
    ssc_vanka_space = VertexCenteredVankaSubspace3dCPU1({"p"}, {"ux", "uy", "uz"})
end


local ssc = SequentialSubspaceCorrection(1.0)
ssc:set_vertex_subspace(ssc_vanka_space)

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
    if (aiForward) then uzawa:set_forward_iter(aiForward)  end
    if (aiSchur) then uzawa:set_schur_iter(aiSchur) end
    if (aiBackward) then uzawa:set_backward_iter(aiBackward)  end

    uzawa:set_schur_operator_update(uzawaSchurUpdateOp, weight)
    -- uzawa:set_debug(dbgWriter)
    return uzawa
end

local uzawaForward = {}
local uzawaBackward = {}

local uzawaWeight = 1.0
local uzawaForward1 = createUzawaIteration("p", gs, Jacobi(0.66), nil, uzawaSchurUpdateOp, uzawaWeight)
local uzawaBackward1 = createUzawaIteration("p", nil, Jacobi(0.66), bgs, uzawaSchurUpdateOp, uzawaWeight)

local uzawaForward2 = createUzawaIteration("p", SymmetricGaussSeidel(), Jacobi(0.66), nil, uzawaSchurUpdateOp, uzawaWeight)
local uzawaBackward2 = createUzawaIteration("p", nil, Jacobi(0.66), SymmetricGaussSeidel(), uzawaSchurUpdateOp, uzawaWeight)

uzawaForward[3] = createUzawaIteration("p", SymmetricGaussSeidel(), SymmetricGaussSeidel(), nil, uzawaSchurUpdateOp, uzawaWeight)
uzawaBackward[3] = createUzawaIteration("p", nil, SymmetricGaussSeidel(), SymmetricGaussSeidel(), uzawaSchurUpdateOp, uzawaWeight)

uzawaForward[0] = createUzawaIteration("p", Jacobi(0.5), Jacobi(0.5), nil, uzawaSchurUpdateOp, uzawaWeight)
uzawaBackward[0] = createUzawaIteration("p", nil, Jacobi(0.5), Jacobi(0.5), uzawaSchurUpdateOp, uzawaWeight)

uzawaForward[4] = createUzawaIteration("p", CG(Jacobi()), CG(Jacobi()), nil, uzawaSchurUpdateOp, uzawaWeight)
uzawaBackward[4] = createUzawaIteration("p", nil, CG(Jacobi()), CG(Jacobi()), uzawaSchurUpdateOp, uzawaWeight)
--[[
local pi=LinearIteratorProduct()
pi:add_iterator(SymmetricGaussSeidel())
pi:add_iterator(SymmetricGaussSeidel())

local uzawaForward2 = createUzawaIteration("p", pi, SymmetricGaussSeidel(), nil, uzawaSchurUpdateOp, uzawaWeight)
local uzawaBackward2 = createUzawaIteration("p", nil, SymmetricGaussSeidel(), pi, uzawaSchurUpdateOp, uzawaWeight)
--]]

--local uzawaForward = createUzawaIteration("p", Jacobi(0.66), Jacobi(0.66), nil, uzawaSchurUpdateOp, uzawaWeight)
--local uzawaBackward = createUzawaIteration("p", nil, Jacobi(0.66), Jacobi(0.66), uzawaSchurUpdateOp, uzawaWeight)
local uzawaSym = createUzawaIteration("p", gs, sgs, bgs, uzawaSchurUpdateOp, uzawaWeight)
--local uzawaBackward = createUzawaIteration("p", nil, Jacobi(0.5), Jacobi(0.66), uzawaSchurOp, uzawaWeight)
local uzawa = uzawaForward



local preSmoother
local postSmoother

if (ARGS.MGSmootherType == "uzawa") then
    preSmoother  = uzawaForward1
    postSmoother = uzawaBackward1
elseif (ARGS.MGSmootherType == "uzawa0") then
    preSmoother  = uzawaForward[0]
    postSmoother = uzawaBackward[0]
elseif (ARGS.MGSmootherType == "uzawa2") then
    preSmoother  = uzawaForward2
    postSmoother = uzawaBackward2
elseif (ARGS.MGSmootherType == "uzawa3") then
    preSmoother  = uzawaForward[3]
    postSmoother = uzawaBackward[3]
elseif (ARGS.MGSmootherType == "uzawa4") then
    preSmoother  = uzawaForward[4]
    postSmoother = uzawaBackward[4]
elseif (ARGS.MGSmootherType == "cgs") then
    preSmoother  = cgs
    postSmoother = cgs
elseif (ARGS.MGSmootherType == "vanka-ssc") then
    preSmoother  = ssc
    postSmoother = ssc
elseif (ARGS.MGSmootherType == "sgs") then
    preSmoother  = sgs
    postSmoother = sgs
else
    quit()
end

-------------------------
-- create GMG
-- TODO: Move to
-------------------------

-- Base Solver
local	baseConvCheck = ConvCheck()
baseConvCheck:set_maximum_steps(5000)
baseConvCheck:set_reduction(1e-12)
baseConvCheck:set_verbose(false)

local	base = BiCGStab()
base:set_preconditioner(jac)
base:set_convergence_check(baseConvCheck)

local	baseCG = CG()
baseCG:set_preconditioner(jac)
baseCG:set_convergence_check(baseConvCheck)

-- exact base solver
local	baseLU = LU()
local	superLU = SuperLU()


-- Geometric Multi Grid
local	gmg = GeometricMultiGrid(approxSpace)
gmg:set_discretization(domainDiscT)
gmg:set_base_level(ARGS.MGBaseLevel)  -- was 1 in Cincy
gmg:set_base_solver(superLU)  -- was baseLU in Cincy
gmg:set_presmoother(preSmoother) --(jac)
gmg:set_postsmoother(postSmoother)
gmg:set_cycle_type(ARGS.MGCycleType) -- 1:V, 2:W -- "F"
gmg:set_num_presmooth(ARGS.MGNumSmooth)
gmg:set_num_postsmooth(ARGS.MGNumSmooth)
gmg:set_rap(problem.bRAP or true)  -- mandatory, if set_stationary


-- gmg:set_debug(dbgWriter)
if (problem.bRAP) then print ("gmg:bRAP=true")
else print ("gmg:bRAP=false") end
-- gmg:set_debug(dbgWriter)



local transfer = StdTransfer()
transfer:enable_p1_lagrange_optimization(true)
--transfer:set_debug(dbgWriter)
gmg:set_transfer(transfer)

local gmgP = GeometricMultiGrid(approxSpace)
-- gmgP:set_discretization(domainDiscP)
gmgP:set_base_level(numPreRefs)  -- was 1 in Cincyj
gmgP:set_base_solver(baseLU)  -- was baseLU in Cincy
gmgP:set_presmoother(sgs)
gmgP:set_postsmoother(sgs)
gmgP:set_cycle_type("V") -- 1:V, 2:W -- "F"
gmgP:set_num_presmooth(3)
gmgP:set_num_postsmooth(3)
gmgP:set_rap(true)  -- mandatory, if set_stationary


local gmgU = GeometricMultiGrid(approxSpace)
-- gmgU:set_discretization(domainDiscU)
gmgU:set_base_level(numPreRefs)  -- was 1 in Cincyj
gmgU:set_base_solver(baseLU)  -- was baseLU in Cincy
gmgU:set_presmoother(sgs)
gmgU:set_postsmoother(sgs)
gmgU:set_cycle_type("V") -- 1:V, 2:W -- "F"
gmgU:set_num_presmooth(3)
gmgU:set_num_postsmooth(3)
gmgU:set_rap(true)  -- mandatory, if add_ionary
gmgU:set_debug(dbgWriter)

local uzawaTotal       = createUzawaIteration("p", ILUT(1e-8), ILUT(1e-8), nil, uzawaSchurUpdateOp, 1.0)      -- ???
local fixedStressLU    = createUzawaIteration("p", nil, ILUT(1e-12), ILUT(1e-12), uzawaSchurUpdateOp, 1.0)
local fixedStressSuperLU    = createUzawaIteration("p", nil, SuperLU(), SuperLU(), uzawaSchurUpdateOp, 1.0)
-- local fixedStressLU    = createUzawaIteration("p", nil, SuperLU(), SuperLU(), uzawaSchurUpdateOp, 1.0)
local fixedStressMG    = createUzawaIteration("p", nil, gmgU, gmgP, uzawaSchurUpdateOp, 1.0)


--local transfer = StdTransfer()
--transfer:enable_p1_lagrange_optimization(false)
--gmg:set_transfer(transfer)

--------------------------------
-- debug solver /iter
--------------------------------

local p0 = 1.0
if (problem.modelParameter and problem.modelParameter.p0) then
    p0 = problem.modelParameter.p0
end


local cmpConvCheck = CompositeConvCheck(approxSpace)
cmpConvCheck:set_component_check("ux", p0*1e-14, 1e-6)
cmpConvCheck:set_component_check("uy", p0*1e-14, 1e-6)
if (dim==3) then
    cmpConvCheck:set_component_check("uz", p0*1e-14, 1e-6)
end
cmpConvCheck:set_component_check("p", p0*1e-14, 1e-6)
cmpConvCheck:set_maximum_steps(100)
cmpConvCheck:set_supress_unsuccessful(true)

local cmpConvCheck2 = CompositeConvCheck(approxSpace)
cmpConvCheck2:set_component_check("ux", p0*1e-12, 1e-6)
cmpConvCheck2:set_component_check("uy", p0*1e-12, 1e-6)
if (dim==3) then
    cmpConvCheck2:set_component_check("uz", p0*1e-12, 1e-6)
end
cmpConvCheck2:set_component_check("p", p0*1e-12, 1e-6)
cmpConvCheck2:set_maximum_steps(50)

cmpConvCheck2 = ConvCheck(200, 1e-25, 1e-20)

local dbgSolver = LinearSolver()
dbgSolver:set_preconditioner(gmg) -- cgs, gmg, uzawa
dbgSolver:set_convergence_check(cmpConvCheck2)
dbgSolver:set_convergence_check(cmpConvCheck)

local dbgIter= DebugIterator()
dbgIter:set_preconditioner(gmg)  -- gmg is the 'real' preconditioner
dbgIter:set_solver(dbgSolver)
dbgIter:set_solution(dbgVector)
dbgIter:set_random_bounds(-5e-6, 5e-6)
dbgIter:set_debug(dbgWriter)  -- print t_0 anf t_N

--------------------------------
-- create and choose a Solver
--------------------------------

local solver = {}

local convCheck = ConvCheck()
convCheck:set_maximum_steps(50)
convCheck:set_reduction(1e-8)
convCheck:set_minimum_defect(1e-14)
-- convCheck = cmpConvCheck  -- for DEBUGGING purposes

local iluSolver = LinearSolver()
iluSolver:set_preconditioner(ilut)
iluSolver:set_convergence_check(convCheck)

local jacSolver = LinearSolver()
jacSolver:set_preconditioner(jac)
jacSolver:set_convergence_check(convCheck)

solver["UzawaSmoother"] = LinearSolver()
solver["UzawaSmoother"]:set_preconditioner(uzawaForward2)
solver["UzawaSmoother"]:set_convergence_check(convCheck)

solver["GMG"] = LinearSolver()
solver["GMG"]:set_preconditioner(dbgIter) -- gmg, dbgIter
solver["GMG"]:set_convergence_check(convCheck) -- cmpConvCheck
solver["GMGKrylov"] = BiCGStab()
solver["GMGKrylov"]:set_preconditioner(gmg) -- gmg, dbgIter
solver["GMGKrylov"]:set_convergence_check(convCheck) -- cmpConvCheck


solver["FixedStressEX"] = LinearSolver()
solver["FixedStressEX"]:set_preconditioner(fixedStressSuperLU)
solver["FixedStressEX"]:set_convergence_check(convCheck)

solver["FixedStressEXKrylov"] = BiCGStab()
solver["FixedStressEXKrylov"]:set_preconditioner(fixedStressSuperLU)
solver["FixedStressEXKrylov"]:set_convergence_check(convCheck)

solver["FixedStressMG"] = LinearSolver() -- BiCGStab()
solver["FixedStressMG"]:set_preconditioner(fixedStressMG)
solver["FixedStressMG"]:set_convergence_check(convCheck)

solver["SuperLU"] = SuperLU() -- SuperLU

solver["LU"] = LinearSolver()
solver["LU"]:set_preconditioner(LU())
solver["LU"]:set_convergence_check(convCheck)

local myIter = gmg
if (ARGS.useDebugIter) then myIter =  dbgIter end

local bicgstabSolver = BiCGStab()
bicgstabSolver:set_preconditioner(myIter) --(gmg)
bicgstabSolver:set_convergence_check(convCheck)

local cgSolver = CG()
cgSolver:set_preconditioner(myIter) --(gmg)
cgSolver:set_convergence_check(convCheck)


local gmresSolver = GMRES(3)
gmresSolver:set_preconditioner(myIter) -- gmg, dbgIter
gmresSolver:set_convergence_check(convCheck)

local sluSolver = SuperLU()

-- local luSolver = LinearSolver()
-- luSolver:set_preconditioner(LU())
-- luSolver:set_convergence_check(convCheck)

-- Select solver.
local lsolver = solver[ARGS.solverID]
--solver = jacSolver
--lsolver = iluSolver
--lsolver = gmgSolver
--lsolver = cgSolver
-- lsolver = bicgstabSolver
--lsolver = gmresSolver
--lsolver:set_compute_fresh_defect_when_finished(true)
-- lsolver = sluSolver

local vtk=VTKOutput()
vtk:select_nodal("p", "PNodal")
if (dim == 2) then vtk:select({"ux", "uy"}, "uNodal") end
if (dim == 3) then vtk:select({"ux", "uy", "uz"}, "uNodal") end

--vtk:select_element( displacementEqDisc:displacement(), "DispElem")
--vtk:select_element( displacementEqDisc:divergence(), "DivElem")
--vtk:select_element( flowEqDisc:gradient()dbgSolver, "GradP")
--vtk:select(massLinker, "Mass")

-- Init error estimator.
local biotErrorEst
if (problem.error_estimator) then
    biotErrorEst = problem:error_estimator()
else
    biotErrorEst = util.biot.CreateDefaultErrorEst(dim)
end



--------------------------------------------------------------------------------
--------------------------------------------------------------------------------
--  Solve transient (linear) problem
--------------------------------------------------------------------------------
--------------------------------------------------------------------------------
if (doTransient) then



    local newtonCheck = ConvCheck()
    newtonCheck:set_maximum_steps(10)
    newtonCheck:set_minimum_defect(1e-14)
    newtonCheck:set_reduction(5e-6)
    newtonCheck:set_verbose(true)

    local newtonCheck2 = CompositeConvCheck(approxSpace)
    newtonCheck2:set_component_check("ux", p0*1e-7, 5e-6)
    newtonCheck2:set_component_check("uy", p0*1e-7, 5e-6)
    if (dim==3) then
        newtonCheck2:set_component_check("uz", p0*1e-7, 5e-6)
    end
    newtonCheck2:set_component_check("p", p0*1e-9, 5e-6)
    newtonCheck2:set_maximum_steps(2)

    local newtonSolver = NewtonSolver()
    newtonSolver:set_linear_solver(lsolver)
    newtonSolver:set_convergence_check(newtonCheck)

    --[[
    local lineSearch = StandardLineSearch();
    lineSearch:set_maximum_steps(6)
    lineSearch:set_accept_best(true)
    newtonSolver:set_line_search(lineSearch)
    --]]

    --newtonSolver:set_debug(dbgWriter)

    local nlsolver = newtonSolver

    print(lsolver:config_string())

    if (problem.check) then problem:check(u) end

    print("Interpolation start values")
    problem:interpolate_start_values(u, startTime)

    -- Create callback.
    function myStepCallback0(u, step, time)
        -- problem:post_processing(u, step, time)
        vtk:print("PoroElasticity.vtu", u, step, time)
    end

    print ("Integrating from 0.0 to "..endTime)


    --dt =dt*1e-4*problem:get_char_time() -- smaller => more complicated

    local charTime = problem:get_char_time()
    dt = 1e-2*problem:get_char_time()
    dtMin = 1e-2*dt

    --
    local myclock = CuckooClock()
    local stepClock = CuckooClock()

    if (( ARGS.LimexNStages > 0)) then
        local dt0 = charTime*1e-50
        print("Computing consistent initial value w/ dt0="..dt0)
        util.SolveNonlinearTimeProblem(u, domainDisc0, nlsolver, myStepCallback0, "PoroElasticityInitial",
                "ImplEuler", 1, startTime, dt0, dt0, dt0, dtRed);


    end



    if ( ARGS.LimexNStages==0) then

        -- Execute linear solver test suite.
        convCheck:set_reduction(1e-10)
        convCheck:set_maximum_steps(100)

        local dtTestSet = {1.0, 0.1, 0.01, 1e-3, 1e-4, 1e-6, 1e-8, 0.0}
        for index,dtvalue in ipairs(dtTestSet) do
            dt = dtvalue*charTime
            endTime = dt
            print("%DTFACTOR=\t"..dtvalue.."\tindex=\t"..index)
            problem:interpolate_start_values(u, startTime)
            myclock:tic()
            test.SolveNonlinearTimeProblem(u, domainDiscT, nlsolver, myStepCallback0, "SolverTest"..index,
                    "ImplEuler", 1, startTime, endTime, dt, dt, dtRed);
            print("MYCLOCK="..myclock:cuckoo().."; "..myclock:toc())
        end


    elseif ( ARGS.LimexNStages==1) then
        -- STANDARD (implicit Euler) time-stepping.

        util.SolveNonlinearTimeProblem(u, domainDiscT, nlsolver, myStepCallback0, "PoroElasticityTransient",
                "ImplEuler", 1, startTime, endTime, dt, dtMin, 0.5,false, myStepCallback0);

        --util.SolveNonlinearTimeProblem(u, domainDiscT, nlsolver, myStepCallback0, "PoroElasticityTransient",
        --             "ImplEuler", 1, startTime, endTime, dt, dtMin, 0.5);

    else
        --  ARGS.LimexNStages > 1

        -- LIMEX time-stepping.

        -- Adjust NEWTON for LIMEX.
        newtonCheck:set_maximum_steps(1)
        newtonCheck:set_supress_unsuccessful(true)

        -- Create & configure LIMEX descriptor.
        local limexDesc = {
            nstages = ARGS.LimexNStages,
            steps = {1,2,3,4,5,6,7,8,9,10},

            domainDisc=domainDiscT,
            nonlinSolver = nlsolver,

            tol = ARGS.LimexTOL,
            dt = dt,
            dtmin = dtMin,
            dtmax = dtMax,

            -- gammaDiscOPT= gammaTensorDisc,  -- no gamma for linear problem
        }


        -- Call factory.
        local limex = util.limex.CreateIntegrator(limexDesc)
        limex:add_error_estimator(biotErrorEst)
        -- limex:set_tolerance(0.001)
        -- limex:set_time_step(dt)
        -- limex:set_dt_min(dtMin)
        -- limex:set_dt_max(dtMax)
        limex:set_stepsize_safety_factor(0.25)
        limex:set_increase_factor(10)
        limex:set_stepsize_greedy_order_factor(0.0)

        limex:disable_matrix_cache()        -- This problem is linear
        --limex:set_time_derivative(udot)   --


        -- Create (& attach) observers.
        if (ARGS.useVTK) then
            local vtkFull = VTKOutput()
            local vtkobserver = VTKOutputObserver("PoroElasticityLimex.vtk", vtkFull)
            limex:attach_observer(vtkobserver)
        end

        local luaobserver = LuaCallbackObserver()

        function myLuaLimexPostProcess(step, time, currdt)
            print ("Time per step :"..stepClock:toc()) -- get time for last step
            local usol=luaobserver:get_current_solution()
            problem:post_processing(usol, step, time)
            stepClock:tic() -- reset timing
            return 0;
        end

        luaobserver:set_callback("myLuaLimexPostProcess")
        limex:attach_observer(luaobserver)



        -- Solve problem using LIMEX.
        myclock:tic()
        limex:apply(u, endTime, u, startTime)
        print("CDELTA="..myclock:toc())
    end

end -- doTransient

--------------------------------------------------------------------------------
--------------------------------------------------------------------------------
--  Solve linear, steady state problem (easy!)
--------------------------------------------------------------------------------
--------------------------------------------------------------------------------

if (doSteadyState) then

    local A = MatrixOperator()
    u = GridFunction(approxSpace)
    local b = GridFunction(approxSpace)
    u:set(0.0)
    b:set(0.0)

    -- 1. assemble matrix and rhs
    domainDiscT:assemble_linear(A, b)

    -- 2. set dirichlet values in start iterate
    u:set(0.0)
    domainDiscT:adjust_solution(u)

    -- 3. init solver for linear Operator
    lsolver:init(A)

    SaveMatrixForConnectionViewer(u, A, "Stiffness.mat")

    -- 4. apply solver
    u:set_random(0.0, 1.0)
    lsolver:apply_return_defect(u,b)
    vtk:print("PoroElasticitySteadyState", u, 0, 0.0)


end  --  doSteadyState




walltime:stop()

if space_time_communicator:get_temporal_rank() == 0 then
    print(integration_time .. " seconds for time integration")
    print(walltime:get() .. " seconds wall time  (rank=" .. space_time_communicator:get_temporal_rank() .. ")")
end
