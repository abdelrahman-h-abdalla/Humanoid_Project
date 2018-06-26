clear; clc; close all;

% Initialization and definition
% ---------------------------------------

delta = 0.02;
pre_time = 0.6;
h_com = 0.26;
g = 9.8;
omega = sqrt(g/h_com);
ss_time = 0.2;
ds_time = 0.1;
w = 0.04;
first_support_foot = 'left'; %doesn't work

parameters1 = Parameters(0.26, 0.01, 0.2, 0.1, first_support_foot, 0.6, 60, 0.04, 0.125, 0.05, 0.025, 'periodic')
hCom 0.26
delta 0.01
singleSupportDuration 0.2
doubleSupportDuration 0.1
firstSupportFoot first_support_foot
predictionTime 0.6
nSamples 60
footSize 0.04
feasibilityL 0.125
feasibilityX 0.05
feasibilityY 0.025
tailType 'periodic'

sim_time = 30;
N_sim = round(sim_time/delta);

%x, xd, xdd, y, yd, ydd, theta, time
initial_state = [0, 0, 0, 0, 0, 0, 0, 0];

solv1 = Solver(initial_state, parameters1);

solv1.set_footsteps([0.0,-0.1;...
                     0.05,0.1;...
                     0.05,-0.1;...
                     0.05,0.1;...
                     0.05,-0.1;...
                     0.05,0.1;...
                     0.05,-0.1]);

plot_options = PlotOptions();
plot_options.plot_com = 1;
plot_options.plot_zmp = 1;
plot_options.plot_pred_zmp = 1;
plot_options.plot_footsteps = 1;
plot_options.plot_pred_footsteps = 1;
plot_options.plot_pred_zmp_constraints = 1;
plot_options.plot_pred_footstep_constraints = 1;
plot_options.plot_orientation = 0;
solv1.set_plot_options(plot_options);

% Set Gains: QZd, QJ, QVx, QVy, QZx, QZy, Qfsx, Qfsy, Q_omega, Q_obstacle
solv1.set_gains(1,0,10,10,0,0,0,0);

solv1.set_plot_limits(-0.1,0.8,-0.2,0.2);

% Simulation cycle
% ---------------------------------------

disp('Simulation cycle')

for i = 1 : N_sim 

    exit_var = solv1.cycle(i);
    if exit_var
        break
    end
    
    % plot
    clf
    
    hold on
    solv1.plot();
    hold off
    grid on
    
    drawnow
end
