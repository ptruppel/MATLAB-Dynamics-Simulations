    function [] = PhillipTruppelli_ConvergenceAnalysis()
%EULERS_METHOD_PROJECTILE_MOTION_DEMO This script simulates the exact
%and numericaly simulated motion of a simple pendulum. It uses matlabs internal 
%ODE23 solver as well as a local function using Euler's Method to simulate 
%the motion. 
clear;
clc;
%% Define Properties
%theta_o = 10*(pi/180);    %uncomment change the theta_o value to switch cases
theta_o = pi/1.1;     %comment this out then
theta_o_dot= 0;

% projectile properties
g = 9.81;       % m/s - acceleration due to gravity
M = 1;          % kg - mass
L = 1;          % m - length of string

% control run time
N_points = 100000;                 % number of points to simulate in Eulers Method

%ODE23 options - used below in ODE23 solver
ode_options = odeset('RelTol',1e-6,...
    'AbsTol',1e-6);

%% Exact Solution
% initial velocity components
w_o=sqrt(g/L);
A = theta_o;
B = (theta_o_dot/w_o);

% angular displacement vs time
theta = @(t) A*cos(w_o*t)+B*sin(w_o*t);             %anonymous function for theta
theta_dot = @(t) -A*w_o*sin(w_o*t)+B*w_o*cos(w_o*t);%anonymous function for theta_dot

% Generate time vector to simulate exact motion
t = linspace(0, 10, N_points);
%% Setup ODE Function Handle
% Matlab's ODE solvers require you to give them function handles of only 2
% variables. They must have the time t first, and the state of the system
% as the second argument. 
odefun = @(t,State) my_ode(t, State, g, L);

% my_ode is a local function defined below. It takes 3 things as input
% arguments. It takes 1) t, 2) state, 3) g, and 4) L
%% Evaluate Motion - ODE23
% input quantities 
State_0 = [theta_o; theta_o_dot];       % initial state
tspan = [0, 10];             % initial / end time to simulate

% Run the ODE solver 
[t_1,State_1] = ode23(odefun, tspan, State_0, ode_options);

t_plot=t_1;
State_plot=State_1;

% Extract trajectory from state vector for easier use
theta_1 = State_1(:,1);
theta_dot = State_1(:,2);

%% Evaluate Motion - Euler's Method
% calculate step size - Euler's method uses a fixed step size, which we
% need to tell it before it can run.
delta_t_array = logspace(-5,-1,30);
delT0 = 10 / N_points;
% Run Euler's Method
energy = @(S_final) ((L.*(1-cos(S_final(1))))*M*g)+((1/2).*M.*(L.*S_final(2)).^2);
%[t_2, State_2] = eulers_method(odefun, tspan, State_0, delT0);
[del_S, del_E, run_time]=convergence_analysis(delta_t_array,odefun, tspan, State_0,energy);

%% Plot Convergence 
h=figure(1);
clf;

subplot(2,1,1);
grid on

yyaxis left
h(1)=loglog(delta_t_array(1:29),abs(del_S));  
h(1).LineWidth = 2;
h(1).Color='b';
h(1).LineStyle='-';
title('\theta Convergence');
xlabel('\Deltat (s)');
ylabel('\Delta{\theta}_f (rad)');

yyaxis right
h(2)=loglog(delta_t_array, run_time);           % results from ODE 23
h(2).LineWidth = 1.5;
h(2).Color='r';
h(2).LineStyle='-';
ylabel('Run Time (s)');
%
subplot(2,1,2);
grid on;

yyaxis left
h(3)=loglog(delta_t_array(1:29),abs(del_E));  
h(3).LineWidth = 2;
h(3).Color='b';
h(3).LineStyle='-';
title('Energy Convergence');
ylabel('\DeltaE_f (J)');
xlabel('\Deltat (s)');

yyaxis right
h(4)=loglog(delta_t_array, run_time);           % results from ODE 23
h(4).LineWidth = 1.5;
h(4).Color='r';
h(4).LineStyle='-';
ylabel('Run Time (s)');

%Print a png with 200 dots per inch resolutionn
print(gcf, 'Convergence_Plots','-dpng','-r200');
 end

%% Local Functions
function dSdt = my_ode(t, State, g, L)  
%MY_ODE calculate the general velocity vector for 2D projectile motion.
%   This function takes the second order ODE for a projectile and then
%   converts it into two first order differential equations.

% Rate of change of State at time t
dSdt(1,1) = State(2);       %velocities
dSdt(2,1) = -(g/L)*sin(State(1));
end

function [State,t] = eulers_method(odefun, tspan, State_0, delT0)
%EULERS_METHOD simulates the response of a system of ODES with a specified
%time step. 

% Determine relevant numbers
nVars = numel(State_0);                 %number of variables
t = (tspan(1) : delT0 : tspan(2))'; 
nT = numel(t);         %number of time steps

% Initialize storage
State = zeros(nT, nVars);   % State(:,1) = x(t) and State(:,2) = y(t), etc.
State(1,:) = State_0;       % set initial values

for i=1 : (nT-1)
    
    % Rate of change of State at given instant
     dSdt = odefun(t(i),State(i,:));
         
    % Store Updated State
     State(i+1,1) = State(i,1)+dSdt(1,1)*delT0;
     State(i+1,2) = State(i,2)+dSdt(2,1)*delT0;
end
end
%% Convergence Analysis
function [del_S, del_E, run_time]=convergence_analysis(delta_t_array,odefun, tspan, State_0,energy)
n = numel(delta_t_array);
del_S=zeros([1 29]);
theta_final=zeros([1 29]);
run_time=zeros([1 30]);
energy_final=zeros([1 30]);
for j=1:n
    tic;
    delta_t_array(j);
    S=eulers_method(odefun,tspan,State_0,delta_t_array(j));
    run_time(j)=toc;
    S_final=S(end,:);
    theta_final(j)=S_final(1);
    energy_final(j)=energy(S_final);
    
end
del_S=diff(theta_final);
del_E=diff(energy_final);

end
