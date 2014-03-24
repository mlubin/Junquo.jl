#############################################################################
# Junquo
# (Ju)lia
# (n)onconvex 
# (qu)adratically constrained quadratic program 
# (o)ptimizer
# 
# A (mixed-integer) nonconvex quadratically constrained quadratic 
# program (QCQP) solver. http://github.com/IainNZ/Junquo.jl
#############################################################################

# The model type used internally by JuMP
type JunquoMathProgModel <: AbstractMathProgModel
    inner::JunquoModel
end
JunquoMathProgModel(;options...) = return JunquoMathProgModel(JunquoModel())

# The solver the user can select when building a model in JuMP
export JunquoSolver
immutable JunquoSolver <: AbstractMathProgSolver
    options
end
JunquoSolver(;kwargs...) = JunquoSolver(kwargs)


#############################################################################
# BEGIN MATHPROGBASE INTERFACE
#############################################################################

model(s::JunquoSolver) = JunquoMathProgModel(;s.options...)

function loadproblem!(m::JunquoMathProgModel, A, collb, colub, obj, rowlb, rowub, sense)
    m.inner.A     = A
    m.inner.collb = collb
    m.inner.colub = colub
    m.inner.obj   = obj
    m.inner.rowlb = rowlb
    m.inner.rowub = rowub
    m.inner.sense = sense
end

function setquadobjterms!(m::JunquoMathProgModel, rowidx, colidx, quadval)
    m.inner.qobj_rowidx   = rowidx
    m.inner.qobj_colidx   = colidx
    m.inner.qobj_val      = quadval
end

function addquadconstr!(m::JunquoMathProgModel, linearidx, linearval,
                                                quadrowidx, quadcolidx, quadval,
                                                sense, rhs)
    push!(m.inner.qcons, QuadCon(linearidx, linearval,
                                        quadrowidx, quadcolidx, quadval,
                                        sense, rhs))
end

optimize!(m::JunquoMathProgModel) = solveJunquo(m.inner)
status(m::JunquoMathProgModel) = :Optimal
getobjval(m::JunquoMathProgModel) = m.inner.objval
getsolution(m::JunquoMathProgModel) = m.inner.sol
    