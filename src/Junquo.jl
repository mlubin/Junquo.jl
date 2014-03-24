#############################################################################
# Junquo
# (Ju)lia
# (n)onconvex 
# (qu)adratically constrained adratic program 
# (o)ptimizer
# 
# A (mixed-integer) nonconvex quadratically constrained quadratic 
# program (QCQP) solver. http://github.com/IainNZ/Junquo.jl
#############################################################################

module Junquo

    # Mathprogbase interface machinery
    require(joinpath(Pkg.dir("MathProgBase"),"src","MathProgSolverInterface.jl"))
    importall MathProgSolverInterface


    type JunquoModel    
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
    JunquoModel() = JunquoModel(nothing,nothing,nothing,nothing,
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


    include("JunquoSolverInterface.jl")
end