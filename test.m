%% Initialisering
clear; clc; close all; 

addpath(genpath('Resources'))
addpath(genpath('Library'))

%% 
load('./Resources/sample_data/0001.mat');

%% Check EKG has right side up:
x_adjusted = change_EKG_direction(x);

%% Process EKG into surrogate signals: HRV, EDR
tic 
tm = RRTachogram.wrapper(x_adjusted, fs) ;
toc
%% plotting
ts_up = tm.getTM('EKG up');
EKG_up = tm.getX('EKG up');

% plot
fig = figure;
fig.Position = [10 100 1500 600];

ax(1) = subplot(211);
plot(ts_up/60, EKG_up); hold on;
plot(ts_up(tm.getTM('qrs up'))/60, EKG_up(tm.getTM('qrs up')), '*')
xlabel('min');
legend('EKG', 'qrs');
grid minor;

ax(2) = subplot(212); 
plot(tm.getTM('HRV')/60, tm.getX('HRV')); hold on;
plot(tm.getTM('HRV')/60, tm.getDis('HRV'));
xlabel('min');
legend('HRV', 'dis');
grid minor;

linkaxes(ax, 'x');
xlim([ts_up(1)/60 ts_up(end)/60 ])
