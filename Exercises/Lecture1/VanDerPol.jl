using DifferentialEquations
using Plots

## Initialize variables

μ = 10.0
x₀ = [2.0 0.0]
tspan = (0.0, 5.0*μ)

## Model

function VanDerPol(ẋ,x₀,tspan,μ)
# VANDERPOL Implementation of the Van der Pol model
#
# Syntax: ẋ = VanDerPol(x,t,μ)
    ẋ = zeros(2)
    x₁, x₂ = x₀
    ẋ[1] = x₂
    ẋ[2] = μ*(1-x₁*x₁)*x₂-x₁
    return ẋ
end

function JacVanDerPol(Jac,x₀,tspan,μ)
# JACVANDERPOL Jacobian for the Van der Pol Equation
#
# Syntax: Jac = JacVanDerPol(x,t,μ)
    Jac = zeros(2, 2)
    x₁, x₂ = x₀
    Jac[2,1] = -2*μ*x₁*x₂-1.0
    Jac[1,2] = 1.0
    Jac[2,2] = μ*(1-x₁*x₁)
    return Jac
end

ẋ = ODEFunction(VanDerPol;jac=JacVanDerPol)

## Define the problem and solve

prob = ODEProblem(ẋ,x₀,tspan,μ)
sol = solve(prob, QNDF(), reltol=1e-6, abstol=1e-6)

## Data Visualization

plot(sol)
