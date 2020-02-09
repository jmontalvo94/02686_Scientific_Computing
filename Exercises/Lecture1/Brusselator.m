%% Driver
a = 1;
b = 3;
x0 = [0; 0]; % Try with different combinations of 0 and 4
t0=0;
tf=100;


options = odeset('Jacobian',@JacBrusselator,'RelTol',1.0e-6,'AbsTol',1.0e-6);
[T,X]=ode15s(@BrusselatorModel,[t0 tf],x0,options,a,b);

%% Data Visualization

tiledlayout(2,1)
nexttile
plot(T,X(:,1))
nexttile
plot(T,X(:,2))

tiledlayout(1,1)
plot(X(:,1),X(:,2))

%% Model

function xdot = BrusselatorModel(t,x,a,b)
% BRUSSELATOR The Brusselator Problem Model
%
% Syntax: xdot = BrusselatorPredator(t,x,a,b)
xdot = zeros(2,1);
xdot(1) = a+x(1)^2*x(2)-(b+1)*x(1);
xdot(2) = b*x(1)-x(1)^2*x(2);
end

function Jac = JacBrusselator(t,x,a,b)
Jac = zeros(2,2);
Jac(1,1) = 2*x(1)*x(2)-b+1;
Jac(2,1) = b-2*x(1)*x(2);
Jac(1,2) = x(1)^2;
Jac(2,2) = -x(1)^2;
end
