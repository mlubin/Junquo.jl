module MIQCQP

    # Standard interface
    require(joinpath(Pkg.dir("MathProgBase"),"src","MathProgSolverInterface.jl"))
    importall MathProgSolverInterface

    type MIQCQPModel    
        A
        collb
        colub
        obj

        rowlb
        rowub
        sense

        qobj_rowidx
        qobj_colidx
        qobj_val

        qcons

        sol
        objval
    end
    MIQCQPModel() = MIQCQPModel(nothing,nothing,nothing,nothing,
                                nothing,nothing,nothing,
                                nothing,nothing,nothing,
                                Any[],nothing,0.0)


    type QuadCon
        linearidx
        linearval
        quadrowidx
        quadcolidx
        quadval
        sense
        rhs
    end

    include("solver.jl")


    include("MIQCQPSolverInterface.jl")
end