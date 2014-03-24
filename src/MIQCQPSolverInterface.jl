export MIQCQPSolver

type MIQCQPMathProgModel <: AbstractMathProgModel
    inner::MIQCQPModel
end
function MIQCQPMathProgModel(;options...)
    return MIQCQPMathProgModel(MIQCQPModel())
end

immutable MIQCQPSolver <: AbstractMathProgSolver
    options
end
MIQCQPSolver(;kwargs...) = MIQCQPSolver(kwargs)


model(s::MIQCQPSolver) = MIQCQPMathProgModel(;s.options...)

function loadproblem!(m::MIQCQPMathProgModel, A, collb, colub, obj, rowlb, rowub, sense)
    m.inner.A     = A
    m.inner.collb = collb
    m.inner.colub = colub
    m.inner.obj   = obj
    m.inner.rowlb = rowlb
    m.inner.rowub = rowub
    m.inner.sense = sense
end

function setquadobjterms!(m::MIQCQPMathProgModel, rowidx, colidx, quadval)
    m.inner.qobj_rowidx   = rowidx
    m.inner.qobj_colidx   = colidx
    m.inner.qobj_val      = quadval
end

function addquadconstr!(m::MIQCQPMathProgModel, linearidx, linearval,
                                                quadrowidx, quadcolidx, quadval,
                                                sense, rhs)
    push!(m.inner.qcons, QuadCon(linearidx, linearval,
                                        quadrowidx, quadcolidx, quadval,
                                        sense, rhs))
end

function optimize!(m::MIQCQPMathProgModel)
    solveMIQCQP(m.inner)
end

function status(m::MIQCQPMathProgModel)
    return :Optimal
end

getobjval(m::MIQCQPMathProgModel) = m.inner.objval

function getsolution(m::MIQCQPMathProgModel)
    return m.inner.sol
end
    