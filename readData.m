%% Script to read data and plot the T-S Diagram
clc; clear; close all
filenames = {"isobaric23.txt","isobaric45.txt","isobaric51.txt"};

%% 2-> 3
filename = filenames{1};
data = readtable(filename);
T_23 = flip(data.Temperature_K_);
p_23 = flip(data.Pressure_MPa_);
s_23 = flip(data.Entropy_J_g_K_);
h_23 = flip(data.Enthalpy_kJ_kg_);
s2 = s_23(1);
s3 = s_23(end);
T2 = T_23(1);
T3 = T_23(end);

%%  4 -> 5
filename = filenames{2};
data = readtable(filename);
T_45 = data.Temperature_K_;
p_45 = data.Pressure_MPa_;
s_45 = data.Entropy_J_g_K_;
h_45 = data.Enthalpy_kJ_kg_;
s4 = s_45(1);
s5 = s_45(end);
T4 = T_45(1);
T5 = T_45(end);


%% 5 -> 1
filename = filenames{3};
data = readtable(filename);
T_51 = data.Temperature_K_;
p_51 = data.Pressure_MPa_;
s_51 = data.Entropy_J_g_K_;
h_51 = data.Enthalpy_kJ_kg_;
T1 = T_51(end);
s1 = s_51(end);

%% Estimate shit by cubic spline (or use polynomial if you want)
% 1->2 and 2->3 not isobaric nor isothermal nor isochoric, no database

T_12 = linspace(T1, T2, 50); 
s_12 = spline([T1 T2], [s1 s2], T_12);
% s_12 = polyval(polyfit([T1 T2], [s1 s2], 3), T_12);

T_34 = linspace(T3, T4, 50);
s_34 = spline([T3 T4], [s3 s4], T_34);


%% Plot
figure;
hold on
plot(s_12, T_12, 'k')
plot(s_23, T_23, 'k')
plot(s_34, T_34, 'k')
plot(s_45, T_45, 'k')
plot(s_51,T_51, 'k')

% Plot and label key state points
plot(s1, T1, 'o', 'MarkerSize', 6, 'MarkerFaceColor', 'k')
text(s1, T1, ' 1', 'FontSize', 10, 'VerticalAlignment', 'bottom')

plot(s2, T2, 'o', 'MarkerSize', 6, 'MarkerFaceColor', 'k')
text(s2, T2, ' 2', 'FontSize', 10, 'VerticalAlignment', 'bottom')

plot(s3, T3, 'o', 'MarkerSize', 6, 'MarkerFaceColor', 'k')
text(s3, T3, ' 3', 'FontSize', 10, 'VerticalAlignment', 'bottom')

plot(s4, T4, 'o', 'MarkerSize', 6, 'MarkerFaceColor', 'k')
text(s4, T4, ' 4', 'FontSize', 10, 'VerticalAlignment', 'bottom')

plot(s5, T5, 'o', 'MarkerSize', 6, 'MarkerFaceColor', 'k')
text(s5, T5, ' 5', 'FontSize', 10, 'VerticalAlignment', 'bottom')

% Labels and grid
xlabel('Entropy [J/g·K]')
ylabel('Temperature [K]')
grid on
title('T–S Diagram')