factory = factory or {}
function factory.create_uzawa_iteration(sSchurCmp, aiForward, aiSchur, aiBackward, uzawaSchurUpdateOp, uzawaSchurWeight)
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
function factory.trash_uzawa()
    if ARGS.smootherID == "uzawa1" then
        print("using uzawa 1")
        preSmoother = create_uzawa_iteration("p", gs, SymmetricGaussSeidel(), nil, uzawaSchurUpdateOp, uzawaWeight)
        postSmoother = create_uzawa_iteration("p", nil, SymmetricGaussSeidel(), bgs, uzawaSchurUpdateOp, uzawaWeight)
    end
    --if ARGS.smootherID == "uzawa2" then
    --    print("using uzawa 2")
    --    preSmoother = create_uzawa_iteration("p", SymmetricGaussSeidel(), Jacobi(0.66), nil, uzawaSchurUpdateOp, uzawaWeight)
    --    postSmoother = create_uzawa_iteration("p", nil, Jacobi(0.66), SymmetricGaussSeidel(), uzawaSchurUpdateOp, uzawaWeight)
    --end
    if ARGS.smootherID == "uzawa3" then
        print("using uzawa 3")
        preSmoother = create_uzawa_iteration("p", SymmetricGaussSeidel(), SymmetricGaussSeidel(), nil, uzawaSchurUpdateOp, uzawaWeight)
        postSmoother = create_uzawa_iteration("p", nil, SymmetricGaussSeidel(), SymmetricGaussSeidel(), uzawaSchurUpdateOp, uzawaWeight)
    end
    if ARGS.smootherID == "uzawa5" then
        print("using uzawa 1")
        preSmoother = createUzawaIteration("p", SymmetricGaussSeidel(), SymmetricGaussSeidel(), nil, uzawaSchurUpdateOp, uzawaWeight)
        postSmoother = createUzawaIteration("p", nil, SymmetricGaussSeidel(), bgs, uzawaSchurUpdateOp, uzawaWeight)
    end

    if ARGS.smootherID == "uzawa6" then
        print("using uzawa 1")
        preSmoother = createUzawaIteration("p", gs, SymmetricGaussSeidel(), nil, uzawaSchurUpdateOp, uzawaWeight)
        postSmoother = createUzawaIteration("p", nil, SymmetricGaussSeidel(), SymmetricGaussSeidel(), uzawaSchurUpdateOp, uzawaWeight)
    end
end

function factory.create_presmoother()

end


function factory.create_lsolver(selector, approxSpace,domainDiscT,preSmoother,postSmoother)
    local p0 = 0.0
    local transfer = StdTransfer()
    transfer:enable_p1_lagrange_optimization(true)

    local superLU = SuperLU() --LU()
    local gmg = GeometricMultiGrid(approxSpace)
    gmg:set_discretization(domainDiscT)
    gmg:set_base_level(0)  -- was 1 in Cincy
    gmg:set_base_solver(superLU)  -- was baseLU in Cincy
    gmg:set_presmoother(preSmoother) --(jac)
    gmg:set_postsmoother(postSmoother)
    gmg:set_cycle_type("W") -- 1:V, 2:W -- "F"
    gmg:set_num_presmooth(2)
    gmg:set_num_postsmooth(2)
    gmg:set_rap(true)  -- mandatory, if set_stationary

    gmg:set_transfer(transfer)

    local cmpConvCheck = CompositeConvCheck(approxSpace)
    cmpConvCheck:set_component_check("ux", p0 * 1e-14, 1e-6)
    cmpConvCheck:set_component_check("uy", p0 * 1e-14, 1e-6)
    if (dim == 3) then
        cmpConvCheck:set_component_check("uz", p0 * 1e-14, 1e-6)
    end
    cmpConvCheck:set_component_check("p", p0 * 1e-14, 1e-6)
    cmpConvCheck:set_maximum_steps(200)
    cmpConvCheck:set_verbose(true)
    cmpConvCheck:set_supress_unsuccessful(true)

    local convCheck = ConvCheck()
    convCheck:set_maximum_steps(200)
    convCheck:set_reduction(1e-8)
    convCheck:set_minimum_defect(1e-14)
    convCheck:set_verbose(true)
    convCheck:set_supress_unsuccessful(false)



    local solver = {}

    solver["GMG"] = LinearSolver()
    solver["GMG"]:set_preconditioner(gmg) -- gmg, dbgIter
    solver["GMG"]:set_convergence_check(convCheck)

    solver["GMGKrylov"] = BiCGStab()
    solver["GMGKrylov"]:set_preconditioner(gmg) -- gmg, dbgIter
    solver["GMGKrylov"]:set_convergence_check(convCheck) -- convCheck

    solver["SuperLU"] = LinearSolver()
    solver["SuperLU"]:set_preconditioner(SuperLU())
    solver["SuperLU"]:set_convergence_check(convCheck)

    solver["LU"] = LinearSolver()
    solver["LU"]:set_preconditioner(LU())
    solver["LU"]:set_convergence_check(convCheck)
    return solver[selector]
end
function factory.create_nlsolver(lsolver)

    local newtonCheck = ConvCheck()
    newtonCheck:set_maximum_steps(1)
    newtonCheck:set_reduction(9.999999999999e-1)
    newtonCheck:set_minimum_defect(1e-14)
    newtonCheck:set_verbose(true)
    newtonCheck:set_supress_unsuccessful(true)

    local newtonSolver = NewtonSolver()
    newtonSolver:set_linear_solver(lsolver)
    newtonSolver:set_convergence_check(newtonCheck)
    return newtonSolver
end
function factory.create_lsolver_coarse(selector,approxSpace,domainDiscT,coarse,preSmoother,postSmoother)
    if coarse == 0 then
        coarse = 100
    end

    local p0 = 1.0

    local cmpConvCheckCoarse = CompositeConvCheck(approxSpace)
    cmpConvCheckCoarse:set_component_check("ux", p0 * 1e-14, 1e-6)
    cmpConvCheckCoarse:set_component_check("uy", p0 * 1e-14, 1e-6)
    if (dim == 3) then
        cmpConvCheckCoarse:set_component_check("uz", p0 * 1e-14, 1e-6)
    end
    cmpConvCheckCoarse:set_component_check("p", p0 * 1e-14, 1e-6)
    cmpConvCheckCoarse:set_maximum_steps(coarse)
    cmpConvCheckCoarse:set_verbose(true)
    cmpConvCheckCoarse:set_supress_unsuccessful(true)

    local convCheckCoarse = ConvCheck()
    convCheckCoarse:set_maximum_steps(coarse)
    convCheckCoarse:set_reduction(1e-8)
    convCheckCoarse:set_minimum_defect(1e-14)
    convCheckCoarse:set_verbose(true)
    convCheckCoarse:set_supress_unsuccessful(true)



    local transfer = StdTransfer()
    transfer:enable_p1_lagrange_optimization(true)
    local superLU = SuperLU() --LU()

    -- coarse solver
    local gmgCoarse = GeometricMultiGrid(approxSpace)
    gmgCoarse:set_discretization(domainDiscT)
    gmgCoarse:set_base_level(0)  -- was 1 in Cincy
    gmgCoarse:set_base_solver(superLU)  -- was baseLU in Cincy
    gmgCoarse:set_presmoother(preSmoother) --(jac)
    gmgCoarse:set_postsmoother(postSmoother)
    gmgCoarse:set_cycle_type("W") -- 1:V, 2:W -- "F"
    gmgCoarse:set_num_presmooth(2)
    gmgCoarse:set_num_postsmooth(2)
    gmgCoarse:set_rap(true)  -- mandatory, if set_stationary
    gmgCoarse:set_transfer(transfer)


    local solverCoarse = {}
    solverCoarse["GMG"] = LinearSolver()
    solverCoarse["GMG"]:set_preconditioner(gmgCoarse) -- gmg, dbgIter
    solverCoarse["GMG"]:set_convergence_check(convCheckCoarse)

    solverCoarse["GMGKrylov"] = BiCGStab()
    solverCoarse["GMGKrylov"]:set_preconditioner(gmgCoarse) -- gmg,
    --solverCoarse["GMGKrylov"]:set_convergence_check(convCheckCoarse)
    solverCoarse["GMGKrylov"]:set_convergence_check(convCheckCoarse)

    solverCoarse["LU"] = LinearSolver()
    solverCoarse["LU"]:set_preconditioner(LU())
    solverCoarse["LU"]:set_convergence_check(convCheckCoarse)

    return solverCoarse[selector]
end



function factory.create_nlsolver_coarse(lsolver)
    local newtonCheckCoarse = ConvCheck()
    newtonCheckCoarse:set_maximum_steps(1)
    newtonCheckCoarse:set_reduction(5e-6)
    newtonCheckCoarse:set_minimum_defect(1e-14)
    newtonCheckCoarse:set_verbose(true)
    newtonCheckCoarse:set_supress_unsuccessful(true)

    local newtonSolverCoarse = NewtonSolver()
    newtonSolverCoarse:set_linear_solver(lsolver)
    newtonSolverCoarse:set_convergence_check(newtonCheckCoarse)
    return newtonSolverCoarse
end