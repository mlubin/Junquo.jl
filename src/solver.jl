using JuMP
using Gurobi

type Linearization
    product::Variable
    left::Variable
    right::Variable
end
function Linearization(m::Model, xleft::Variable, xright::Variable)
    left_lb     = getLower(xleft)
    left_ub     = getUpper(xleft)
    right_lb    = getLower(xright)
    right_ub    = getUpper(xright)
    xproduct    = Variable(m, min(left_lb*right_lb, left_lb*right_ub,
                                  left_ub*right_lb, left_ub*right_ub),
                              max(left_lb*right_lb, left_lb*right_ub,
                                  left_ub*right_lb, left_ub*right_ub), 0)
    return Linearization(xproduct, xleft, xright)
end
function addBounds!(m::Model, lin::Linearization)
    x = lin.left;  y = lin.right
    xL = getLower(x);   xU = getUpper(x)
    yL = getLower(y);   yU = getUpper(y)
    @addConstraint(m, lin.product >= yL*x + xL*y - xL*yL)
    @addConstraint(m, lin.product >= yU*x + xU*y - xU*yU)
    @addConstraint(m, lin.product <= yL*x + xU*y - xU*yL)
    @addConstraint(m, lin.product <= yU*x + xL*y - xL*yU)
end
function printError(lin::Linearization)
    p_val = getValue(lin.product)
    L_val = getValue(lin.left)
    R_val = getValue(lin.right)
    #println(lin.product.col,    " ", p_val, " ",
    #        lin.left.col,       " ", L_val, " ",
    #        lin.right.col,      " ", R_val, " ",
    #        abs(p_val - L_val*R_val))
    println(p_val, " ~ ", L_val, " * ", R_val, " = ", L_val*R_val, 
            ", err=", abs(p_val - L_val*R_val), 
            " (", lin.product.col, ",", lin.left.col, ",", lin.right.col, ")")
end
function calcError(lin::Linearization)
    p_val = getValue(lin.product)
    L_val = getValue(lin.left)
    R_val = getValue(lin.right)
    return p_val != 0 ? abs(p_val - L_val*R_val)/abs(p_val) : 0
end



function solveRelax(data::MIQCQPModel, collb, colub)
    num_var = length(collb)
    m = Model(solver=GurobiSolver(OutputFlag=0))
    @defVar(m, collb[i] <= x[i=1:num_var] <= colub[i])

    


    # Analyze quad constraints to find bilinear terms
    lins = Dict{(Int,Int),Linearization}()
    for c in data.qcons
        num_quad = length(c.quadrowidx)
        for j in 1:num_quad
            left_ind  = int64(c.quadrowidx[j])
            right_ind = int64(c.quadcolidx[j])
            key = (left_ind,right_ind)
            if !(key in lins)
                lins[key] = Linearization(m, x[left_ind], x[right_ind])
            end
        end
    end
    # and objective
        num_quad = length(data.qobj_rowidx)
        for j in 1:num_quad
            left_ind  = int64(data.qobj_rowidx[j])
            right_ind = int64(data.qobj_colidx[j])
            key = (left_ind,right_ind)
            if !(key in lins)
                lins[key] = Linearization(m, x[left_ind], x[right_ind])
            end
        end

    # Build linearized model
    if data.qobj_val != nothing
        setObjective(m, data.sense, dot(x,data.obj) + 
            sum([data.qobj_val[j] * lins[(data.qobj_rowidx[j],data.qobj_colidx[j])].product
                    for j in 1:length(data.qobj_val)])
        )
    else
        setObjective(m, data.sense, dot(x,data.obj))
    end

    # Linear constraints
    for lin_row = 1:data.A.m
        lhs = AffExpr()
        row = data.A[lin_row,:]
        for j in 1:data.A.n
            if row[j] != 0.0
                push!(lhs, row[j], x[j])
            end
        end
        if data.rowlb[lin_row] == -Inf
            # <= constraint
            addConstraint(m, lhs <= data.rowub[lin_row])
        elseif data.rowub[lin_row] == +Inf
            # >= constraint
            addConstraint(m, lhs >= data.rowlb[lin_row])
        elseif data.rowlb[lin_row] == data.rowub[lin_row]
            # == constraint
            addConstraint(m, lhs == data.rowlb[lin_row])
        end
    end

    # Quadratic constraints
    for quad_row = 1:length(data.qcons)
        con = data.qcons[quad_row]
        lhs = AffExpr()
        for lin_ind = 1:length(con.linearidx)
            push!(lhs, con.linearval[lin_ind], x[con.linearidx[lin_ind]])
        end
        for quad_ind = 1:length(con.quadrowidx)
            L_ind = con.quadrowidx[quad_ind]
            R_ind = con.quadcolidx[quad_ind]
            push!(lhs, con.quadval[quad_ind], lins[(L_ind,R_ind)].product)
        end
        con.sense == '<'  &&  addConstraint(m, lhs <= con.rhs)
        con.sense == '>'  &&  addConstraint(m, lhs >= con.rhs)
    end

    #println(m)

    # Linearization
    for leftright in keys(lins)
        addBounds!(m, lins[leftright])
    end

    status = solve(m,suppress_warnings=true)
    if status == :Infeasible
        return -Inf, nothing, nothing, nothing
    end


    for leftright in keys(lins)
        #printError(lins[leftright])
    end

    errors = zeros(num_var)
    for leftright in keys(lins)
        L,R=leftright
        err=calcError(lins[leftright])
        errors[L] += err
        errors[R] += err
    end

    for j = 1:num_var
        #println(j, "  ", errors[j])
    end
    max_error, max_error_j = findmax(errors)

    return getObjectiveValue(m), getValue(x)[:], max_error, max_error_j
end


function solveMIQCQP(data::MIQCQPModel)
    srand(10)
    const ERR_TOL = 1e-5
    const GAP_TOL = 0.001

    # Always make it a maximization problem
    flipped_obj = false
    if data.sense == :Min
        data.obj        = -1 .* data.obj
        if data.qobj_val != nothing
            data.qobj_val   = -1 .* data.qobj_val
        end
        data.sense      = :Max
        flipped_obj     = true
    end

    # Normalize indices so left/row is <= right/col
    for c in data.qcons
        num_quad = length(c.quadrowidx)
        new_rowidx = zeros(Int64, num_quad)
        new_colidx = zeros(Int64, num_quad)
        for j in 1:num_quad
            left_ind  = int64(c.quadrowidx[j])
            right_ind = int64(c.quadcolidx[j])
            if left_ind > right_ind
                left_ind, right_ind = right_ind, left_ind
            end
            new_rowidx[j] = left_ind
            new_colidx[j] = right_ind
        end
        c.quadrowidx = new_rowidx
        c.quadcolidx = new_colidx
    end
        num_quad = data.qobj_rowidx == nothing ? 0 : length(data.qobj_rowidx)
        new_rowidx = zeros(Int64, num_quad)
        new_colidx = zeros(Int64, num_quad)
        for j in 1:num_quad
            left_ind  = int64(data.qobj_rowidx[j])
            right_ind = int64(data.qobj_colidx[j])
            if left_ind > right_ind
                left_ind, right_ind = right_ind, left_ind
            end
            new_rowidx[j] = left_ind
            new_colidx[j] = right_ind
        end
        data.qobj_rowidx = new_rowidx
        data.qobj_colidx = new_colidx

    # Prepare the nonlinear solver
    nlp_model = buildNLP(data)

    collb_queue = Any[]
    colub_queue = Any[]
    relax_queue = Any[]
    error_queue = Any[]
    push!(collb_queue, data.collb)
    push!(colub_queue, data.colub)
    push!(relax_queue, +Inf)
    push!(error_queue, +Inf)

    obj_upper_bound = +Inf
    obj_lower_bound = -Inf
    incumbent = nothing

    # While still problems to solve
    iter = 0
    old_gap     = nothing
    last_action = nothing
    ratio = 1
    while length(collb_queue) > 0  
        if iter % 20 == 0
            print(lpad("ITER",       5," ")," |")       # 7
            print(lpad("LBOUND",    10," ")," |")       # 12  19
            print(lpad("UBOUND",    10," ")," |")       # 12  31
            print(lpad("CURRENT",   10," ")," |")       # 12  43
            print(lpad("GAP%",      10," ")," |")       # 12  55
            println(" ACTION")
        end
        iter += 1

        # Update upper bound based on remaining unsolved nodes
        obj_upper_bound = maximum(relax_queue)
        gap = abs(obj_upper_bound - obj_lower_bound)/abs(obj_lower_bound)
        if gap < GAP_TOL
            println("Gap closed at Iter $iter, gap = $(gap*100)% < $(GAP_TOL*100)%, final bound $obj_lower_bound")
            break
        end

        # Node selection strategy:
        # "Pick node with smallest error" (depth-first, raises lower bound)
        # "Pick node with biggest relaxation" (breadth-first, lowers upper bound)
        selected_node_index = 0
        if incumbent == nothing #|| (incumbent != nothing && rand() < 0.1)
            selected_node_index = indmin(error_queue)
        else
            selected_node_index = indmax(relax_queue)
        end
        collb = splice!(collb_queue, selected_node_index)
        colub = splice!(colub_queue, selected_node_index)
        relax = splice!(relax_queue, selected_node_index)
        _     = splice!(error_queue, selected_node_index)

        # Check if already bad
        if relax <= obj_lower_bound
            print(lpad(string(iter),5," ")," |")
            println(" FATHOM BASED ON PARENT $relax")
            continue
        end


        z, x, max_err, max_err_j = solveRelax(data, collb, colub)

        print(lpad(string(iter),5," ")," |")                            # 7
        print(lpad(string(round(obj_lower_bound,4)),10," ")," |")       # 12  19
        print(lpad(string(round(obj_upper_bound,4)),10," ")," |")       # 12  31
        print(lpad(string(round(z,              4)),10," ")," |")       # 12  43
        print(lpad(string(round(gap*100,3)),         9," "),"% |")      # 12  55

        # Check if infeasible
        if z == -Inf && x == nothing
            println(" INFEASIBLE")
            continue
        end

        # Check if we can fathom
        if z <= obj_lower_bound
            println(" FATHOM (Z<=LBOUND)")
            continue
        end

        # Check if error below tolerance
        if max_err <= ERR_TOL
            # Try to improve
            new_z = solveNLP(nlp_model, data, x)
            #new_z = -Inf #solveNLP(nothing, data, x)
            if new_z > z
                print("   IMP [$(round(z,4)),$(round(new_z,4))]")
                z = new_z
            else
                print(" NOIMP")
            end

            if z > obj_lower_bound
                obj_lower_bound = z
                incumbent = x[:]
                println(" NEW LBOUND")
            else
                println(" FEASIBLE")
            end
            continue
        end

        # Branch
        k = max_err_j
        mid_point = (collb[k] + colub[k])/2
        branch_point = 0.15*mid_point + 0.85*x[k]
        println(" BRANCH ON $k [$(round(collb[k],4)),$(round(x[k],4)),$(round(colub[k],4))]")
        down_branch_collb       = collb[:]
        down_branch_colub       = colub[:]
        down_branch_colub[k]    = branch_point    
        #println("    Down branch: was [$(collb[k]), $(colub[k])], now [$(down_branch_collb[k]), $(down_branch_colub[k])]")
        push!(collb_queue, down_branch_collb)
        push!(colub_queue, down_branch_colub)
        push!(relax_queue, z)
        push!(error_queue, max_err)

        up_branch_collb         = collb[:]
        up_branch_colub         = colub[:]
        up_branch_collb[k]      = branch_point
        #println("    Up   branch: was [$(collb[k]), $(colub[k])], now [$(up_branch_collb[k]), $(up_branch_colub[k])]")
        push!(collb_queue, up_branch_collb)
        push!(colub_queue, up_branch_colub)
        push!(relax_queue, z)
        push!(error_queue, max_err)

        #iter >= 10 && break
    end  # while

    data.sol    = incumbent
    data.objval = obj_lower_bound * (flipped_obj ? -1 : +1)
end


function buildNLP(data::MIQCQPModel)
    num_var = length(data.collb)
    nlp = Model()
    @defVar(nlp, data.collb[i] <= x[i=1:num_var] <= data.colub[i])

    # OBJECTIVE
    if data.qobj_val != nothing
        num_quad = length(data.qobj_val)
        linobj = data.obj
        qobj   = data.qobj_val
        rowx   = [x[data.qobj_rowidx[j]] for j=1:num_quad]
        colx   = [x[data.qobj_rowidx[j]] for j=1:num_quad]
        @setNLObjective(nlp, data.sense, 
            sum{linobj[j] * x[j], j=1:num_var} + 
            sum{qobj[j] * rowx[j] * colx[j], j=1:num_quad}
        )
    else
        setObjective(nlp, data.sense, dot(x,data.obj))
    end

    # LINEAR CONSTRAINTS
    for lin_row = 1:data.A.m
        lhs = AffExpr()
        row = data.A[lin_row,:]
        for j in 1:data.A.n
            if row[j] != 0.0
                push!(lhs, row[j], x[j])
            end
        end
        lb = data.rowlb[lin_row]
        ub = data.rowub[lin_row]
        if lb == -Inf
            # <= constraint
            addConstraint(nlp, lhs <= ub)
        elseif ub == +Inf
            # >= constraint
            addConstraint(nlp, lhs >= lb)
        elseif lb == ub
            # == constraint
            addConstraint(nlp, lhs == lb)
        end
    end

    # QUADRATIC CONSTRAINTS
    for quad_row = 1:length(data.qcons)
        con = data.qcons[quad_row]::QuadCon
        num_lin     = length(con.linearval)
        num_quad    = length(con.quadval)
        rhs = con.rhs
        linval = con.linearval
        linx = [x[con.linearidx[lin_ind]] for lin_ind = 1:num_lin]
        qval = con.quadval
        qvar1 = [x[con.quadrowidx[q_ind]] for q_ind = 1:num_quad]
        qvar2 = [x[con.quadcolidx[q_ind]] for q_ind = 1:num_quad]
        if con.sense == '<'
            @addNLConstraint(nlp,
                sum{ linval[lin_ind]*linx[lin_ind], lin_ind = 1:num_lin} +
                sum{ qval[q_ind]*qvar1[q_ind]*qvar2[q_ind], q_ind = 1:num_quad}
                <= rhs
            )
        elseif con.sense == '>'
            @addNLConstraint(nlp,
                sum{ linval[lin_ind]*linx[lin_ind], lin_ind = 1:num_lin} +
                sum{ qval[q_ind]*qvar1[q_ind]*qvar2[q_ind], q_ind = 1:num_quad}
                >= rhs
            )
        end
    end

    return nlp
end

function solveNLP(nlp::Model, data::MIQCQPModel, start_x::Vector{Float64})
#function solveNLP(m, data::MIQCQPModel, start_x::Vector{Float64})


    for j = 1:length(start_x)
        nlp.colVal[j] = start_x[j]
    end

    options=Dict()
    options["print_level"] = 0
    solve(nlp, IPoptOptions=options)  # with IPOPT
    return getObjectiveValue(nlp)
end