%% Driver

% Initialize variables
tf = 10;
t0 = 0;
N = 100;
x0 = [1;1;1];
lambda = [1;0;-1];
dt =(tf-t0)/N;

% EE
[T,X] = EulerExplicitMethodModel(@LinearSystemModel,t0,tf,N,x0,lambda);

%% Data Visualization

%Linear System
tiledlayout(1,1)
plot(T,X,'LineWidth',2)
axis([t0 tf -1 10])
title("Step size: \Deltat = "+dt)
legend({"\lambda= "+lambda(1),"\lambda= "+lambda(2),"\lambda= "+lambda(3)},'Location','northeast')
xlabel('t')
ylabel('x(t)')

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

function xdot = LinearSystemModel(t,x,lambda)
% LINEAR Implementation of a simple scalar linear system
%
% Syntax: xdot = LinearSystem(t,x,lambda)
xdot=zeros(3,1);
xdot(1,1) = lambda(1)*x(1);
xdot(2,1) = lambda(2)*x(2);
xdot(3,1) = lambda(3)*x(3);
end