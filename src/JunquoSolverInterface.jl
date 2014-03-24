export JunquoSolver

type JunquoMathProgModel <: AbstractMathProgModel
    inner::JunquoModel
end
function JunquoMathProgModel(;options...)
    return JunquoMathProgModel(JunquoModel())
end

immutable JunquoSolver <: AbstractMathProgSolver
    options
end
JunquoSolver(;kwargs...) = JunquoSolver(kwargs)


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

function optimize!(m::JunquoMathProgModel)
    solveJunquo(m.inner)
end

function status(m::JunquoMathProgModel)
    return :Optimal
end

getobjval(m::JunquoMathProgModel) = m.inner.objval

function getsolution(m::JunquoMathProgModel)
    return m.inner.sol
end
    