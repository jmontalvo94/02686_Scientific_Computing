%% Driver

% Set parameters
alfa=0.04;
beta=1E4;
gamma=3E7;
p=[alfa;beta;gamma];

% Simulation
t0=0;
tf=100;
x0 = [1;0;0];

options = odeset('Jacobian',@JacRobertson,'RelTol',1.0e-6,'AbsTol',1.0e-6);
[T,X]=ode15s(@RobertsonModel,[t0 tf],x0,options,p);

%% Data Visualization
tiledlayout(3,1)
nexttile
plot(T,X(:,1))
nexttile
plot(T,X(:,2))
nexttile
plot(T,X(:,3))

%% Model

function xdot = RobertsonModel(t,x,p)
% ROBERTSON The Robertson's Chemical Reaction Problem Model
%
% Syntax: xdot = RobertsonModel(t,x,p)

% Differential Equations
xdot = zeros(3,1);
xdot(1) = -p(1)*x(1)+p(2)*x(2)*x(3);
xdot(2) = p(1)*x(1)-p(2)*x(2)*x(3)-p(3)*x(2)^2;
xdot(3) = p(3)*x(2)^2;
end

function Jac = JacRobertson(t,x,p)
Jac = zeros(3,3);
Jac(1,1) = -p(1);
Jac(1,2) = p(2)*x(3);
Jac(1,3) = p(2)*x(2);
Jac(2,1) = p(1);
Jac(2,2) = -p(2)*x(3)-2*p(3)*x(2);
Jac(2,3) = -p(2)*x(2);
Jac(3,1) = 0;
Jac(3,2) = 2*p(3)*x(2);
Jac(3,3) = 0;
end