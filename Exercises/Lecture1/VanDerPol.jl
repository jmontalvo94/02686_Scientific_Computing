using DifferentialEquations

## Driver

μ = 10
x₀ = [2.0 0.0]
t = 5.0*μ
prob = ODEProblem(VanDerPolModel,t,x₀,μ)
options = odeset("Jacobian',@JacVanDerPol,'RelTol',1.0e-6,'AbsTol",1.0e-6)
[T,X]=ode15s(@VanDerPolModel,[0 5*mu],x0,options,mu)

plot(T,X[:,1])
title("\mu="+mu)
xlabel("time")
ylabel("x_[1](t)")

plot(T,X[:,2])
title("\mu="+mu)
xlabel("time")
ylabel("x_[2](t)")

## Model

function VanDerPolModel(t,x,μ)
# VANDERPOL Implementation of the Van der Pol model
#
# Syntax: xdot = VanDerPol[t,x,μ]
    ẋ=zeros(2,1)
    ẋ[1] = x[2]
    ẋ[2] = μ*(1-x[1]*x[1])*x[2]-x[1]
    return ẋ
end

function JacVanDerPol(t,x,μ)
# JACVANDERPOL Jacobian for the Van der Pol Equation
#
# Syntax: Jac = JacVanDerPol(t,x,μ)

    Jac = zeros(2,2)
    Jac[2,1] = -2*μ*x[1]*x[2]-1.0
    Jac[1,2] = 1.0
    Jac[2,2] = μ*(1-x[1]*x[1])
    return Jac
end
