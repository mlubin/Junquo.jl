using JuMP
using Junquo

m = Model(solver=JunquoSolver())

L = zeros(15)
L[1] =  8.0
L[2] = 43.0
L[3] =  3.0

U = zeros(15)
U[1] = 21.0
U[2] = 57.0
U[3] = 16.0
for k in 1:4
    U[3*k+1] =  90.0
    U[3*k+2] = 120.0
    U[3*k+3] =  60.0
end

@defVar(m, L[i] <= x[i=1:15] <= U[i])

setObjective(m, :Min, 
  sum([(2.3     * x[3*k+1]   +
        0.0001  * x[3*k+1]^2 + 
        1.7     * x[3*k+2]   +
        0.0001  * x[3*k+2]^2 +
        2.2     * x[3*k+3]   +
        0.00015 * x[3*k+3]^2) for k in 0:4]))

# constr1
for j in 1:4
    @addConstraint(m, x[3*j+1] - x[3*j-2] + 7 <= 13)
    @addConstraint(m, x[3*j+1] - x[3*j-2] + 7 >=  0)
end

# constr2
for j in 1:4
    @addConstraint(m, x[3*j+2] - x[3*j-1] + 7 <= 14)
    @addConstraint(m, x[3*j+2] - x[3*j-1] + 7 >=  0)
end

# constr3
for j in 1:4
    @addConstraint(m, x[3*j+3] - x[3*j  ] + 7 <= 13)
    @addConstraint(m, x[3*j+3] - x[3*j  ] + 7 >=  0)
end

@addConstraint(m, x[1] + x[2] + x[3]    >= 60)
@addConstraint(m, x[4] + x[5] + x[6]    >= 50)
@addConstraint(m, x[7] + x[8] + x[9]    >= 70)
@addConstraint(m, x[10] + x[11] + x[12] >= 85)
@addConstraint(m, x[13] + x[14] + x[15] >= 100)

solve(m)

println(getObjectiveValue(m), " z")