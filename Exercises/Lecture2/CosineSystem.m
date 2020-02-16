%% Driver

% Cosine System

% Initialize variables
tf = 2*pi*10;
t0 = 0;
N = 1000;
x0 = 1;
dt =(tf-t0)/N;

% EE
[T,X] = EulerExplicitMethodModel(@CosineSystemModel,t0,tf,N,x0);

%% Data Visualization

% Cosine System
tiledlayout(1,1)
plot(T,X,'LineWidth',2)
axis([t0 tf -inf inf])
title("Step size: \Deltat = "+dt)
legend({'cos(t)x(t)'},'Location','northeast')
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

function xdot = CosineSystemModel(t,x)
% COSINE Implementation of a simple cosine system
%
% Syntax: xdot = CosineSystem(t,x)
xdot=zeros(1);
xdot(1) = cos(t)*x(1);
end