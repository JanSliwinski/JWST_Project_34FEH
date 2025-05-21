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
v_23 = flip(data.Volume_m3_kg_);

s2 = s_23(1);
s3 = s_23(end);
T2 = T_23(1);
T3 = T_23(end);
v2 = v_23(1);
v3 = v_23(end);
p2 = p_23(1);
p3 = p_23(end);

%%  4 -> 5
filename = filenames{2};
data = readtable(filename);
T_45 = data.Temperature_K_;
p_45 = data.Pressure_MPa_;
s_45 = data.Entropy_J_g_K_;
h_45 = data.Enthalpy_kJ_kg_;
v_45 = data.Volume_m3_kg_;

s4 = s_45(1);
s5 = s_45(end);
T4 = T_45(1);
T5 = T_45(end);
v4 = v_45(1);
v5 = v_45(end);
p4 = p_45(1);
p5 = p_45(end);

%% 5 -> 1
filename = filenames{3};
data = readtable(filename);
T_51 = data.Temperature_K_;
p_51 = data.Pressure_MPa_;
s_51 = data.Entropy_J_g_K_;
h_51 = data.Enthalpy_kJ_kg_;
v_51 = data.Volume_m3_kg_;

v1 = v_51(end);
T1 = T_51(end);
s1 = s_51(end);
p1 = p_51(end);

%% Estimate shit by cubic spline (or use polynomial if you want)
% 1->2 and 2->3 not isobaric nor isothermal nor isochoric, no database

T_12 = linspace(T1, T2, 50); 
s_12 = spline([T1 305 T2], [s1 24.5 s2], T_12);
% s_12 = polyval(polyfit([T1 T2], [s1 s2], 3), T_12);

p_12 = linspace(p1, p2, 50);
v_12 = spline([p1 0.8 p2], [v1 0.8 v2], p_12);

T_34 = linspace(T3, T4, 50);
s_34 = spline([T3 11 T4], [s3 5.4 s4], T_34);
p_34 = linspace(p3, p4, 50);
v_34 = spline([p3 p4], [v3 v4], p_34);


%% Plot T-s
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

%% Plot p-v
figure;
hold on
plot(v_12, p_12, 'k')
plot(v_23, p_23, 'k')
plot(v_34, p_34, 'k')
plot(v_45, p_45, 'k')
plot(v_51, p_51, 'k')

% Plot and label key state points
plot(v1, p1, 'o', 'MarkerSize', 6, 'MarkerFaceColor', 'k')
text(v1, p1, ' 1', 'FontSize', 10, 'VerticalAlignment', 'bottom')

plot(v2, p2, 'o', 'MarkerSize', 6, 'MarkerFaceColor', 'k')
text(v2, p2, ' 2', 'FontSize', 10, 'VerticalAlignment', 'bottom')

plot(v3, p3, 'o', 'MarkerSize', 6, 'MarkerFaceColor', 'k')
text(v3, p3, ' 3', 'FontSize', 10, 'VerticalAlignment', 'bottom')

plot(v4, p4, 'o', 'MarkerSize', 6, 'MarkerFaceColor', 'k')
text(v4, p4, ' 4', 'FontSize', 10, 'VerticalAlignment', 'bottom')

plot(v5, p5, 'o', 'MarkerSize', 6, 'MarkerFaceColor', 'k')
text(v5, p5, ' 5', 'FontSize', 10, 'VerticalAlignment', 'bottom')

% Labels and grid
xlabel('Specific Volume [m^3/kg]')
ylabel('Pressure [MPa]')
grid on
title('p-v Diagram')

