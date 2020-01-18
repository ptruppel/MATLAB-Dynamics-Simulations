%EULERS_METHOD_PROJECTILE_MOTION_DEMO This script simulates the exact
%and numericaly simulated motion of a simple pendulum. It uses matlabs internal 
%ODE23 solver as well as a local function using Euler's Method to simulate 
%the motion. 
clear;
clc;
%% Define Properties
theta_o = pi/2;
theta_o_dot= 0;

% projectile properties
g = 9.81;       % m/s - acceleration due to gravity
Mu = 6;         % mass ratio (M/m)
m = 1;          % swinging mass
M = 6;          % rising and falling mass
R_o = 1.5;      % m - length of string
r_o = 1.5;      % initial length of swinging arm
r_o_dot = 0;    % The swinging mass starts at rest
% control run time
N_points = 120000;                 % number of points to simulate in Eulers Method
total_frames = 1500;
step_frames = N_points/total_frames;
%ODE23 options - used below in ODE23 solver
ode_options = odeset('RelTol',1e-6,...
    'AbsTol',1e-6);

t = linspace(0, 10, N_points);
%% Setup ODE Function Handle
% Matlab's ODE solvers require you to give them function handles of only 2
% variables. They must have the time t first, and the state of the system
% as the second argument. 
odefun = @(t,State) my_ode(t, State, g, m, M);

% input quantities 
State_0 = [theta_o; theta_o_dot; r_o; r_o_dot];       % initial state
tspan1 = [0,1];
tspan2 = [0, 50];             % initial / end time to simulate


%% Evaluate Motion - Euler's Method
% calculate step size - Euler's method uses a fixed step size, which we
% need to tell it before it can run.
delta_t_array = logspace(-7,-2,30);
delT0 = 10 / N_points;
% Run Euler's Method
energy = @(S_final) (.5*M*S_final(4)^2)+(.5*m*((S_final(4)^2)+(S_final(3)^2)*(S_final(2)^2)))+...
    (M*g*S_final(3))-(m*g*S_final(3)*cos(S_final(1)));
[State, t] = eulers_method(odefun, tspan2, State_0, delT0);
Eulers_State = State;
[del_r_div_r_0, del_E_div_r_0, del_S_div_r_0, run_time]=convergence_analysis(delta_t_array,odefun, tspan1, State_0,energy);

% Plot Convergence 
h=figure(1);
clf;

subplot(3,1,1);
grid on;

yyaxis left
h(3)=loglog(delta_t_array(1:29),abs(del_r_div_r_0));  
h(3).LineWidth = 2;
h(3).Color='b';
h(3).LineStyle='-';
title('Convergence');
ylabel('\Deltar / r_0');
yticks([10^-8 10^-6 10^-4 10^-2]);

yyaxis right
h(4)=loglog(delta_t_array, run_time);           % results from ODE 23
h(4).LineWidth = 1.5;
h(4).Color='r';
h(4).LineStyle='-';
ylabel('Run Time (s)');
yticks([10^-4 10^-2 10^0 10^2]);

subplot(3,1,2);
grid on

yyaxis left
h(1)=loglog(delta_t_array(1:29),abs(del_E_div_r_0));  %E and Theta data are switched
h(1).LineWidth = 2;
h(1).Color='b';
h(1).LineStyle='-';
ylabel('\Delta{\theta}/{\theta}_0');
yticks([10^-8 10^-6 10^-4 10^-2]);

yyaxis right
h(2)=loglog(delta_t_array, run_time);           % results from ODE 23
h(2).LineWidth = 1.5;
h(2).Color='r';
h(2).LineStyle='-';
ylabel('Run Time (s)');
yticks([10^-4 10^-2 10^0 10^2]);

subplot(3,1,3);
grid on;

yyaxis left
h(3)=loglog(delta_t_array(1:29),abs(del_S_div_r_0));  %E and Theta data are switched
h(3).LineWidth = 2;
h(3).Color='b';
h(3).LineStyle='-';
ylabel('\DeltaE/E_0');
xlabel('\Deltat (s)');
yticks([10^-8 10^-6 10^-4 10^-2]);

yyaxis right
h(4)=loglog(delta_t_array, run_time);
h(4).LineWidth = 1.5;
h(4).Color='r';
h(4).LineStyle='-';
ylabel('Run Time (s)');
yticks([10^-4 10^-2 10^0 10^2]);

%Print a png with 200 dots per inch resolution
print(gcf, 'convergence_analysis','-dpng','-r200');

%% Trajectory
x_traj = State(:,3).*sin(State(:,1));
y_traj = -State(:,3).*cos(State(:,1));

m=figure(2);
clf;
ax=gca;

m(3)=plot(x_traj,y_traj);  
m(3).LineWidth = .5;
m(3).Color='b';
m(3).LineStyle='-';
title('\theta_{0}=90^{o} -- \mu=6 -- t_{f}=50 sec');
ax.XLim = [-1.5, 1.5];
ax.YLim = [-1.5, 1];
ylabel('y(m)');
xlabel('x(m)');
daspect([1, 1, 1]);
grid on;

%Print a png with 200 dots per inch resolution
print(gcf, 'trajectory_vs_time','-dpng','-r200');

%% Animation
save_animation = true;
file_name = 'Swinging Atwood Machine';
movie_dur = 7;                  % s - time for movie to last
n = 2e2;                        % number of time points to simulate
frame_dur = movie_dur / n;      % s - frame duration

% Control Appearance
line_width = 1.5;       % linewidth 
marker_size = 6;        % size of projectile

Y_M = State(:,3)-3;
X_M = -1*ones(size(Y_M));

Y_m = y_traj;%-State(:,3).*cos(State(:,1));
X_m = x_traj+1;%State(:,3).*sin(State(:,1));

fig = figure(3);    % create 
clf  
ax = gca;           % store handle for axis 

hold on                         	%allows multiple lines to be plotted
small_traj = plot(X_m(1),Y_m(1),'k');
small_mass = plot(X_m(1),Y_m(1),'ro','MarkerFaceColor','r');
large_mass = plot(X_M(1), Y_M(1),'ro','MarkerFaceColor','r', 'MarkerSize', 15);    
small_line = plot([1,X_m(1)], [0,Y_m(1)], 'b');
large_line = plot([-1,1], [0,0], 'b');
Horiz_line = plot([-1,1], [0,0], 'b');
black_dot1 = plot(-1, 0, 'ko', 'MarkerFaceColor','k');
black_dot2 = plot(1, 0, 'ko', 'MarkerFaceColor','k');
title('Animated Trajectory - \mu=M/m=6');
xlabel('X(m)');
ylabel('Y(m)');

%Add a grid 
grid on

%Axis limits - set these BEFORE animating to make sure the axis is fixed
%while animating the motion. If this isn't done, then the axis will be
%constantly updating while you animatae the motion and it looks terrible 
ax.XLim = [-3, 3];
ax.YLim = [-3, 3];
daspect([1 ,1 ,1]);
% Update Appearance
h3.LineWidth = 2;
h2.MarkerSize = marker_size;

ax.FontWeight = 'bold';
ax.FontSize = 12;     

frames{n} = getframe(gcf);  %initial storage for movie frames
for i = 1:step_frames:N_points
    
    %Update Trajectories   
    
    small_traj.XData = X_m(1:i);
    small_traj.YData = Y_m(1:i);
    small_mass.XData = X_m(i);        %ball at time t = current
    small_mass.YData = Y_m(i);
    large_mass.XData = X_M(i);
    large_mass.YData = Y_M(i);
    small_line.XData = [1,X_m(i)];
    small_line.YData = [0,Y_m(i)];
    large_line.XData = [-1,-1];
    large_line.YData = [0,Y_M(i)];
    horiz_line.XData = [-1,1];
    horiz_line.YData = [0,0];
    black_dot1.XData = -1;
    black_dot1.YData = 0;
    black_dot2.XData = 1;
    black_dot2.YData = 0;
    
    
    drawnow                     %update the plot
    
    % Store the frame. This will be written to a video file in the next cell
    if save_animation
        frames{i} = getframe(gcf);
    end
end

if save_animation
        
    v = VideoWriter(['animated_swinging_atwood','.mp4']);    %Create a video file object
    
    v.FrameRate = 1 / frame_dur;            %Set the frame rate
    
    open(v)                                 %Open video to start writing
    for i = 1:step_frames:N_points
        writeVideo(v,frames{i})             %Write frames
    end
    close(v)                                %Close video / create final file
end


%%
% Local Functions
function dSdt = my_ode(t, State, g, m, M)  
% MY_ODE calculate the general velocity vector for 2D projectile motion.
%   This function takes the second order ODE for a projectile and then
%   converts it into two first order differential equations.

% Rate of change of State at time t
dSdt(1,1) = State(2);   %theta dot
dSdt(2,1) = ((-2*State(4)*State(2))-(g*sin(State(1))))/State(3);%theta dot dot
dSdt(3,1) = State(4); %r dot
dSdt(4,1) = ((m*State(3)*State(2)^2)-(M*g)+(m*g*cos(State(1))))/(M+m); %r dot dot
end

function [State,t] = eulers_method(odefun, tspan, State_0, delT0)
% EULERS_METHOD simulates the response of a system of ODES with a specified
% time step. 

% Determine relevant numbers
nVars = numel(State_0);                 %number of variables
t = (tspan(1) : delT0 : tspan(2))'; 
nT = numel(t);         %number of time steps

% Initialize storage
State = zeros(nT, nVars);   % State(:,1) = x(t) and State(:,2) = y(t), etc.
State(1,:) = State_0;       % set initial values
for i=1 : (nT-1)
    
%     Rate of change of State at given instant
     dSdt = odefun(t(i),State(i,:));
         
%     Store Updated State
     State(i+1,1) = State(i,1)+dSdt(1,1)*delT0; %theta
     State(i+1,2) = State(i,2)+dSdt(2,1)*delT0; %theta dot
     State(i+1,3) = State(i,3)+dSdt(3,1)*delT0; %r
     State(i+1,4) = State(i,4)+dSdt(4,1)*delT0; %r dot
end
end
% Convergence Analysis
function [del_r_div_r_0, del_S_div_r_0, del_E_div_r_0, run_time]=convergence_analysis(delta_t_array,odefun, tspan, State_0,energy)
n = numel(delta_t_array);
del_S=zeros([1 n]);
theta_final=zeros([1 n]);
r_final=zeros([1 n]);
run_time=zeros([1 n]);
energy_final=zeros([1 n]);
for j=1:n
    disp(j);
    tic;
    S=eulers_method(odefun,tspan,State_0,delta_t_array(j));
    run_time(j)=toc;
    S_final=S(end,:);
    theta_final(j)=S_final(1);
    r_final(j)=S_final(3);
    energy_final(j)=energy(S_final);
end
 
del_r_div_r_0=diff(r_final)/State_0(3);
del_S_div_r_0=diff(theta_final)/State_0(1);
del_E_div_r_0=diff(energy_final)/((6*9.81*State_0(3))-(1*9.81*State_0(3)*cos(State_0(1))));
end


