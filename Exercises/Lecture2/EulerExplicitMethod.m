%% Driver

% Linear System

% Initialize variables
% tf = 10;
% t0 = 0;
% N = 100;
% x0 = [1;1;1];
% lambda = [1;0;-1];
% dt =(tf-t0)/N;

% EE
%[T,X] = EulerExplicitMethodModel(@LinearSystem,t0,tf,N,x0,lambda);

% Cosine System

% Initialize variables
% tf = 2*pi*10;
% t0 = 0;
% N = 1000;
% x0 = 1;
% dt =(tf-t0)/N;

% EE
%[T,X] = EulerExplicitMethodModel(@CosineSystem,t0,tf,N,x0);

% VanDerPol System

% Initialize variables
% tf = 50;
% t0 = 0;
% N = 10000;
% x0 = [2;0];
% mu=5;
% dt =(tf-t0)/N;
% 
% % EE
% [T,X] = EulerExplicitMethodModel(@VanDerPol,t0,tf,N,x0,mu);

% PreyPredator System

% Initialize variables
tf = 50;
t0 = 0;
N = 1000;
x0 = [2;2];
a=1;
b=1;
dt =(tf-t0)/N;

% EE
[T,X] = EulerExplicitMethodModel(@PreyPredator,t0,tf,N,x0,a,b);

%% Data Visualization

%Linear System
% tiledlayout(1,1)
% plot(T,X,'LineWidth',2)
% axis([t0 tf -1 10])
% title("Step size: \Deltat = "+dt)
% legend({"\lambda= "+lambda(1),"\lambda= "+lambda(2),"\lambda= "+lambda(3)},'Location','northeast')
% xlabel('t')
% ylabel('x(t)')

% Cosine System
% tiledlayout(1,1)
% plot(T,X,'LineWidth',2)
% axis([t0 tf -inf inf])
% title("Step size: \Deltat = "+dt)
% legend({'cos(t)x(t)'},'Location','northeast')
% xlabel('t')
% ylabel('x(t)')

% VanDerPol System
% tiledlayout(2,1)
% for i=1:size(x0,1)
%     nexttile
%     plot(T,X(:,i), 'LineWidth',2)
%     title("\mu="+mu)
%     xlabel('time')
%     ylabel("x_{"+i+"}(t)")
% end
% tiledlayout(1,1)
% nexttile
% plot(X(:,1),X(:,2), 'LineWidth',2)
% title("\mu="+mu)
% xlabel('x_{1}(t)')
% ylabel('x_{2}(t)')

% PreyPredator System
tiledlayout(2,1)
for i=1:size(x0,1)
    nexttile
    plot(T,X(:,i), 'LineWidth',2)
    title("a = "+a+", b = "+b+", \Deltat = "+dt)
    xlabel('time')
    ylabel("x_{"+i+"}(t)")
end
tiledlayout(1,1)
nexttile
plot(X(:,1),X(:,2), 'LineWidth',2)
title("a = "+a+", b = "+b+", \Deltat = "+dt)
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

function xdot = LinearSystem(t,x,lambda)
% LINEAR Implementation of a simple scalar linear system
%
% Syntax: xdot = LinearSystem(t,x,lambda)
xdot=zeros(3,1);
xdot(1,1) = lambda(1)*x(1);
xdot(2,1) = lambda(2)*x(2);
xdot(3,1) = lambda(3)*x(3);
end

function xdot = CosineSystem(t,x)
% COSINE Implementation of a simple cosine system
%
% Syntax: xdot = CosineSystem(t,x)
xdot=zeros(1);
xdot(1) = cos(t)*x(1);
end

function xdot = VanDerPol(t,x,mu)
% VANDERPOL Implementation of the Van der Pol model
%
% Syntax: xdot = VanDerPol(t,x,mu)
xdot=zeros(2,1);
xdot(1) = x(2);
xdot(2) = mu*(1-x(1)*x(1))*x(2)-x(1);
end

function xdot = PreyPredator(t,x,a,b)
% PREYPREDATOR The Prey-Predator Model
%
% Syntax: xdot = PreyPredator(t,x,a,b)
xdot = zeros(2,1);
xdot(1) = a*(1-x(2))*x(1);
xdot(2) = -b*(1-x(1))*x(2);
end