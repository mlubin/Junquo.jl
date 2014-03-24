using JuMP
using MIQCQP

m = Model(solver=MIQCQPSolver())
@defVar(m, 0 <= x1 <= 1)
@defVar(m, 0 <= x2 <= 1)
@defVar(m, 0 <= x3 <= 1)
@defVar(m, 0 <= x4 <= 1)
@defVar(m, 0 <= x5 <= 1)
@defVar(m, 0 <= x6 <= 1)
@defVar(m, 0 <= x7 <= 1)
@defVar(m, 0 <= x8 <= 1)
@defVar(m, 0 <= x9 <= 1)

addConstraint(m, 
    1.41*x2*x6 + 1.72*x2*x7 + 1.57*x2*x9 - 0.88*x3*x7 + 1.78*x4*x6 + (-1.08*x5*
    x6) - 1.43*x5*x7 + 1.75*x5*x9 - 0.18*x2 + 0.77*x3 - 0.5*x4 - 0.55*x5 + 0.58
    *x6 - 0.61*x7 + 0.21*x9 - 0.17*x8 <= 48.14)

addConstraint(m, 
    0.56*x3*x9 - 1.83*x3*x8 - 1.25*x2*x8 + (-1.4*x4*x7) - 0.31*x4*x9 + 0.93*x5*
    x6 + 0.24*x5*x9 + 0.88*x2 + 0.67*x3 - 0.37*x4 - 0.26*x5 + 0.39*x6 - 0.34*x7
     - 0.13*x8 - 0.3*x9 <= 22.3)

addConstraint(m, 
    1.96*x2*x9 - 0.97*x2*x7 + (-1.53*x3*x6) - 0.23*x3*x8 + 1.34*x3*x9 + 0.26*x4
    *x6 - 0.16*x4*x9 + 0.81*x5*x9 + x2 + 0.01*x3 + 0.43*x4 + 0.14*x5 + 0.77*x6
     + 0.57*x7 - 0.54*x8 - 0.09*x9 <= 52.73)

addConstraint(m, 
    1.88*x2*x6 + 1.05*x2*x9 + (-1.95*x3*x6) - 1.06*x3*x8 + 1.24*x4*x8 + 0.79*x5
    *x8 - 0.29*x5*x6 + 0.77*x5*x9 - 0.75*x2 + 0.6*x3 + 0.84*x4 + 0.99*x5 - 0.71
    *x6 - 0.73*x8 - 0.86*x9 - 0.05*x7 <= 75.04)

addConstraint(m, 
    (-0.75*x2*x7) - 1.13*x2*x9 + 0.12*x3*x6 - 1.47*x3*x7 + 1.3*x4*x7 - 1.22*x4*
    x8 + 1.51*x5*x6 - 1.83*x5*x9 - 0.84*x2 + 0.9*x3 + 0.25*x4 - 0.25*x5 - 0.15*
    x6 + 0.09*x7 - 0.4*x8 - 0.2*x9 <= 76.59)

addConstraint(m, 
    0.09*x2*x6 + 0.26*x2*x7 + 1.37*x2*x8 + 1.52*x3*x7 - 1.63*x3*x8 + 0.88*x4*x7
     - 1.79*x4*x8 + 0.28*x5*x7 - 0.71*x2 - 0.83*x3 + 0.63*x4 - 0.21*x5 + 0.58*
    x6 + 0.41*x7 - 0.27*x8 + 0.23*x9 <= 73.3)

addConstraint(m, 
    (-1.7*x2*x7) - 0.18*x2*x8 + (-1.09*x3*x8) - 1.46*x3*x9 - 1.72*x4*x6 + (-
    1.09*x5*x7) - 1.67*x5*x8 - 1.91*x5*x9 + 0.36*x2 + 0.67*x3 + 0.92*x4 - 0.33*
    x5 - 0.51*x6 - 0.14*x7 - 0.24*x8 - 0.43*x9 <= 43.03)

addConstraint(m, 
    (-1.06*x3*x8) - 1.19*x3*x9 - 1.96*x2*x9 + 1.93*x4*x6 - 0.32*x4*x7 + (-0.24*
    x5*x7) - 0.85*x5*x8 - 0.21*x5*x9 - 0.88*x2 + 0.79*x3 - 0.01*x4 - 0.08*x5 + 
    0.55*x6 + 0.66*x7 + 0.62*x8 + 0.78*x9 <= 63.61)

setObjective(m, :Min,
              (-1.04*x2*x7) - 0.93*x2*x9 + (-1.11*x3*x7) - 1.52*x3*x8 + 1.26*x4
              *x8 - 1.6*x4*x6 - 1.62*x4*x9 + 1.43*x5*x6 + 0.21*x2 - 0.75*x3 + 
              0.93*x4 + 0.88*x5 + 0.35*x6 + 0.67*x7 + 0.74*x8 + 0.92*x9)

solve(m)

println(getValue(x2), " x2")
println(getValue(x3), " x3")
println(getValue(x4), " x4")
println(getValue(x5), " x5")
println(getValue(x6), " x6")
println(getValue(x7), " x7")
println(getValue(x8), " x8")
println(getValue(x9), " x9")
println(getObjectiveValue(m), " z")