    function [] = PhillipTruppelli_Simple_Pendulum()
%EULERS_METHOD_PROJECTILE_MOTION_DEMO This script simulates the exact
%and numericaly simulated motion of a simple pendulum. It uses matlabs internal 
%ODE23 solver as well as a local function using Euler's Method to simulate 
%the motion. 
clear;
clc;
%% Define Properties
%theta_o = 10*(pi/180);    %uncomment change the theta_o value to switch cases
theta_o = 80*(pi/180);     %comment this out then
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

%% Energy Functions 
% Define relevant energy expressions 
% U = @(t) (L.*(1-cos(theta(t))))*M*g;                    %potential
% T = @(t) (1/2).*M.*(L.*theta_dot(t)).^2;                %kinetic
% E = @(t) U(t) + T(t);                                   %total 

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
tspan = [0, 10];              % initial / end time to simulate

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
delT0 = 10 / N_points;

% Run Euler's Method
[t_2, State_2] = eulers_method(odefun, tspan, State_0, delT0);
plot(t_2, State_2(:,1));
% Extract trajectory from state vector for easier use
theta_2 = State_2(:,1);

%% Plot Motion
% Plot the trajectory that was simulated using each method to compare the
% results. 
figure(1)
clf
hold on
title('Theta Initial is 80 Degrees') %uncomment this to do 80 degrees
%title('Theta Initial is 10 Degrees'); %comment this out
%Plot trajectories
 h(1) = plot(t, theta(t));            % plot exact equations 
 h(1).LineWidth = 2;
 h(1).Color='r';
 h(1).LineStyle='-';
 
 h(2) = plot(t_1, theta_1);           % results from ODE 23
 h(2).LineWidth = 1.5;
 h(2).Color='b';
 h(2).LineStyle='--';
 
 h(3) = plot(t_2, theta_2);           % results from Euler's Method
 h(3).LineWidth = 1;
 h(3).Color='m';
 h(3).LineStyle='--';
 
 h(1).DisplayName = 'Exact Solution';
 h(2).DisplayName = 'Eulers Method';
 h(3).DisplayName = 'ode23';
 legend();
 
%Axis appearannce
 grid on
 xlabel('Time (s)');
 ylabel('Theta(t) (rad)');

 hold off 

% Print a png with 200 dots per inch resolutionn
%print(gcf, 'theta_vs_time_state_2','-dpng','-r200');

%% Plot Energy
% Plot the total energy of each particle that was simulated using to 
% compare the results. 
% save_animation = true;

% You should be able to do this entirely on your own 


%% Animation
save_animation = true;
file_name = 'pendulum motion';
movie_dur = 7;                  % s - time for movie to last
n = 2e2;                        % number of time points to simulate
frame_dur = movie_dur / n;      % s - frame duration

% Control Appearance
line_width = 1.5;       % linewidth 
marker_size = 6;        % size of projectile

X = sin(theta_2(1:1:end));
Y = -cos(theta_2(1:1:end))+(L*1.2);

Xmax = max(X);
Ymax = max(Y);

fig = figure(2);    % create 
clf  
subplot(2,1,1);
ax = gca;           % store handle for axis 

hold on                         	%allows multiple lines to be plotted
h1 = plot(0,(1.2*L),'ro','MarkerFaceColor','r');
h2 = plot3(X(1), Y(1), 0,'ro','MarkerFaceColor','r');    
h3 = plot([0,X(1)], [(1.2*L),Y(1)]);
title('Simple Pendulum');
xlabel('X (m)');
ylabel('Y (m)');

%Add a grid 
grid on

%Axis limits - set these BEFORE animating to make sure the axis is fixed
%while animating the motion. If this isn't done, then the axis will be
%constantly updating while you animatae the motion and it looks terrible 
ax.XLim = [-1.1, 1.1];
ax.YLim = [0, 1.2];
pbaspect([2,1,1]);
% Update Appearance
h3.LineWidth = 2;
h2.MarkerSize = marker_size;

ax.FontWeight = 'bold';
ax.FontSize = 12;

U = (L.*(1-cos(State_2(:,1))))*M*g;                    %potential
T = (1/2).*M.*(L.*State_2(:,2)).^2;                %kinetic
E = U + T;      

subplot(2,1,2);
hold on
h4 = plot(U,'m');
h5 = plot(T,'b');
h6 = plot(E,'r');
ax2=gca;
ax2.XLim = [0, 10];
ax2.YLim = [0, E(1)*1.2];
title('Animation of Energy');
xlabel('Time (s)');
ylabel('Energy (J)');
 h4.DisplayName = 'Potential';
 h5.DisplayName = 'Kinetic';
 h6.DisplayName = 'Total';
 leg = legend;
 leg.Location='bestoutside';
frames{n} = getframe(gcf);  %initial storage for movie frames
for i = 1:100:numel(t)
    
    %Update Trajectories   
    h2.XData = X(i);        %ball at time t = current
    h2.YData = Y(i);
    h3.XData = [0,X(i)];
    h3.YData = [(1.2*L),Y(i)];
    h4.XData=t(1:i);
    h4.YData=U(1:i);
    h5.XData=t(1:i);
    h5.YData=T(1:i);
    h6.XData=t(1:i);
    h6.YData=E(1:i);
    drawnow                     %update the plot
    
    % Store the frame. This will be written to a video file in the next cell
    if save_animation
        frames{i} = getframe(gcf);
    end
end

if save_animation
        
    v = VideoWriter(['animated_simple_pendulum','.mp4']);    %Create a video file object
    
    v.FrameRate = 1 / frame_dur;            %Set the frame rate
    
    open(v)                                 %Open video to start writing
    for i = 1:100:numel(t)
        writeVideo(v,frames{i})             %Write frames
    end
    close(v)                                %Close video / create final file
end

%%



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

function [t, State] = eulers_method(odefun, tspan, State_0, delT0)
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
