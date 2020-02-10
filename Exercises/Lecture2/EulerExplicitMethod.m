%% Driver

% Initialize variables
tb=100;
ta=0;
N=10;
x0=[0;0];
t0=0;

EulerExplicitMethodModel

%% Model

function [T,X] = EulerExplicitMethodModel(fun,ta,tb,N,x0,varargin)

% Compute step size and allocate memory
dt = (tb-ta)/N;
nx = size(x0,1);
X = zeros(nx,N+1);
T = zeros(1,N+1);

% Eulers Explicit Method
T(:,1) = t0;
X(:,1) = x0;
for k=1:N
    f = feval(fun,T(k),X(:,k),vararginf:g);
    T(:,k+1) = T(:,k) + dt;
    X(:,k+1) = X(:,k) + dt*f;
end

% Form a nice table for the result
T = T';
X = X';
end