%% Driver

% --------------------------------------------------------------
% Parameters
% --------------------------------------------------------------
a1 = 1.2272; %[cm2] Area of outlet pipe 1
a2 = 1.2272; %[cm2] Area of outlet pipe 2
a3 = 1.2272; %[cm2] Area of outlet pipe 3
a4 = 1.2272; %[cm2] Area of outlet pipe 4

A1 = 380.1327; %[cm2] Cross sectional area of tank 1
A2 = 380.1327; %[cm2] Cross sectional area of tank 2
A3 = 380.1327; %[cm2] Cross sectional area of tank 3
A4 = 380.1327; %[cm2] Cross sectional area of tank 4

gamma1 = 0.45; % Flow distribution constant. Valve 1
gamma2 = 0.40; % Flow distribution constant. Valve 2

g = 981; %[cm/s2] The acceleration of gravity
rho = 1.00; %[g/cm3] Density of water

p = [a1; a2; a3; a4; A1; A2; A3; A4; gamma1; gamma2; g; rho];

% --------------------------------------------------------------
% Simulation scenario
% --------------------------------------------------------------

% Normal scenario
% t0 = 0.0; % [s] Initial time
% tf = 20*60; % [s] Final time
% 
% m10 = 0.0; % [g] Liquid mass in tank 1 at time t0
% m20 = 0.0; % [g] Liquid mass in tank 2 at time t0
% m30 = 0.0; % [g] Liquid mass in tank 3 at time t0
% m40 = 0.0; % [g] Liquid mass in tank 4 at time t0
% 
% F1 = 300; % [cm3/s] Flow rate from pump 1
% F2 = 300; % [cm3/s] Flow rate from pump 2
% 
% x0 = [m10; m20; m30; m40];
% u = [F1; F2];

% Discrete scenario with noise
t0 = 0.0; % [s] Initial time
tf = 20*60; % [s] Final time
Ts = 10; % [s] Sample Time
t = [t0:Ts:tf]'; % [s] Sample instants
N = length(t);

m10 = 0; % [g] Liquid mass in tank 1 at time t0
m20 = 0; % [g] Liquid mass in tank 2 at time t0
m30 = 0; % [g] Liquid mass in tank 3 at time t0
m40 = 0; % [g] Liquid mass in tank 4 at time t0

F1 = 300; % [cm3/s] Flow rate from pump 1
F2 = 300; % [cm3/s] Flow rate from pump 2

x0 = [m10; m20; m30; m40];
u = [repmat(F1,1,N); repmat(F2,1,N)];

% Process Noise
Q = [20^2 0;0 40^2];
Lq = chol(Q,'lower');
w = Lq*randn(2,N);

% Measurement Noise
R = eye(4);
Lr = chol(R,'lower');
v = Lr*randn(4,N);

% --------------------------------------------------------------
% Compute the solution / Simulate
% --------------------------------------------------------------

% Solve the system of differential equations [normal]
%[T,X] = ode45(@FourTankSystemModel,[t0 tf],x0,[],u,p);

% Solve the system of differential equations [discrete with noise]
nx = 4; nu = 2; ny = 4; nz = 2;
x = zeros(nx,N);
y = zeros(ny,N);
z = zeros(nz,N);
X = zeros(0,nx);
T = zeros(0,1);
x(:,1) = x0;

for k = 1:N-1
    y(:,k) = FourTankSystemSensor(x(:,k),p)+v(:,k); % Sensor function
    z(:,k) = FourTankSystemOutput(x(:,k),p); % Output function
    [Tk,Xk] = ode45(@FourTankSystemModel,[t(k) t(k+1)],x(:,k),[],...
        u(:,k)+w(:,k),p);
    x(:,k+1) = Xk(end,:)';
    T = [T; Tk];
    X = [X; Xk];
end

k = N;
y(:,k) = FourTankSystemSensor(x(:,k),p)+v(:,k); % Sensor function
z(:,k) = FourTankSystemOutput(x(:,k),p); % Output function

% --------------------------------------------------------------
% Additional variables for plotting
% --------------------------------------------------------------

% Define help variables
[nT,nX] = size(X);
a = p(1:4,1)';
A = p(5:8,1)';

% Compute the measured variables
H = zeros(nT,nX);
for i=1:nT
    H(i,:) = X(i,:)./(rho*A);
end

% Compute the flows out of each tank
Qout = zeros(nT,nX);
for i=1:nT
    Qout(i,:) = a.*sqrt(2*g*H(i,:));
end

%% Data Visualization
% --------------------------------------------------------------
% Plotting
% --------------------------------------------------------------

Figure2by2(T,X/1000,nX,'Masses','m[kg]') %Divided over 1000 to transform to kg
Figure2by2(T,H,nX,'Height','h[cm]')
Figure2by2(T,Qout,nX,'Outflow Rate','q[cm^3/s]')

%% Model

function xdot = FourTankSystemModel(t,x,u,p)
% FOURTANKSYSTEM Model dx/dt = f(t,x,u,p) for 4-tank System
%
% This function implements a differential equation model for the
% 4-tank system.
%
% Syntax: xdot = FourTankSystem(t,x,u,p)

% Unpack states, MVs, and parameters
m = x; % Mass of liquid in each tank [g]
F = u; % Flow rates in pumps [cm3/s]
a = p(1:4,1); % Pipe cross sectional areas [cm2]
A = p(5:8,1); % Tank cross sectional areas [cm2]
gamma = p(9:10,1); % Valve positions [-]
g = p(11,1); % Acceleration of gravity [cm/s2]
rho = p(12,1); % Density of water [g/cm3]

% Inflows
qin = zeros(4,1);
qin(1,1) = gamma(1)*F(1); % Inflow from valve 1 to tank 1 [cm3/s]
qin(2,1) = gamma(2)*F(2); % Inflow from valve 2 to tank 2 [cm3/s]
qin(3,1) = (1-gamma(2))*F(2); % Inflow from valve 2 to tank 3 [cm3/s]
qin(4,1) = (1-gamma(1))*F(1); % Inflow from valve 1 to tank 4 [cm3/s]

% Outflows
h = m./(rho*A); % Liquid level in each tank [cm]
qout = a.*sqrt(2*g*h); % Outflow from each tank [cm3/s]

% Differential equations
xdot = zeros(4,1);
xdot(1,1) = rho*(qin(1,1)+qout(3,1)-qout(1,1)); % Mass balance Tank 1
xdot(2,1) = rho*(qin(2,1)+qout(4,1)-qout(2,1)); % Mass balance Tank 2
xdot(3,1) = rho*(qin(3,1)-qout(3,1)); % Mass balance Tank 3
xdot(4,1) = rho*(qin(4,1)-qout(4,1)); % Mass balance Tank 4
end



function y = FourTankSystemSensor(x,p)
% FOURTANKSYSTEMSENSOR Level for each tank in the four tank system
%
% Syntax: y = FourTankSystemSensor(x,p)

% Extract states and parameters
m = x;
A = p(5:8,1);
rho = p(12,1);

% Compute level in each tank
rhoA = rho*A;
y = m./rhoA;
end

function z = FourTankSystemOutput(x,p)
% FOURTANKSYSTEMOUTPUT Level for the lower tanks in the four tank system
%
% Syntax: z = FourTankSystemOutput(x,p)

% Extract states and parameters
m = x(1:2,1);
A = p(5:6,1);
rho = p(12,1);

% Compute level in each tank
rhoA = rho*A;
z = m./rhoA;
end

function PlottedFigure = Figure2by2(a,b,n,name,yname)
PlottedFigure = figure;
t=tiledlayout(2,2);
title(t,name)
for i=1:n
    nexttile
    plot(a(:,1)/60,b(:,i)) % Divided over 60 to use minutes instead of second
    title("Tank "+i)
    xlabel('time[min]')
    ylabel(yname)
end
end