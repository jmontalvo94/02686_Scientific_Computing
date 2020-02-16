%% Driver

% VanDerPol System

% Initialize variables
tf = 50;
t0 = 0;
N = 10000;
x0 = [2;0];
mu=5;
dt =(tf-t0)/N;

% EE
[T,X] = EulerExplicitMethodModel(@VanDerPolModel,t0,tf,N,x0,mu);

%% Data Visualization

% VanDerPol System
tiledlayout(2,1)
for i=1:size(x0,1)
    nexttile
    plot(T,X(:,i), 'LineWidth',2)
    title("\mu="+mu)
    xlabel('time')
    ylabel("x_{"+i+"}(t)")
end
tiledlayout(1,1)
nexttile
plot(X(:,1),X(:,2), 'LineWidth',2)
title("\mu="+mu)
xlabel('x_{1}(t)')
ylabel('x_{2}(t)')

%% Model

function [T,X] = EulerExplicitMethodModel(fun,t0,tf,N,x0,varargin)

% Compute step size and allocate memory
dt = (tf-t0)/N;
nx = size(x0,1);
X = zeros(nx,N+1);
T = zeros(1,N+1);

% Eulers Explicit Method
T(:,1) = t0;
X(:,1) = x0;
for k=1:N
    f = feval(fun,T(k),X(:,k),varargin{:});
    T(:,k+1) = T(:,k) + dt;
    X(:,k+1) = X(:,k) + dt*f;
end

% Form a nice table for the result
T = T';
X = X';
end

function xdot = VanDerPolModel(t,x,mu)
% VANDERPOL Implementation of the Van der Pol model
%
% Syntax: xdot = VanDerPol(t,x,mu)
xdot=zeros(2,1);
xdot(1) = x(2);
xdot(2) = mu*(1-x(1)*x(1))*x(2)-x(1);
end