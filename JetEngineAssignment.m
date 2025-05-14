clear ;close all;clc;
warning off
%% To make sure that matlab will find the functions. You must change it to your situation 
relativepath_to_generalfolder='General'; % relative reference to General folder (assumes the folder is in you working folder)
addpath(relativepath_to_generalfolder); 
%% Load Nasadatabase
TdataBase=fullfile('General','NasaThermalDatabase');
load(TdataBase);
%% Nasa polynomials are loaded and globals are set. 
%% values should not be changed. These are used by all Nasa Functions. 
global Runiv Pref
Runiv=8.314472;
Pref=1.01235e5; % Reference pressure, 1 atm!
Tref=298.15;    % Reference Temperature
%% Some convenient units
kJ=1e3;kmol=1e3;dm=0.1;bara=1e5;kPa = 1000;kN=1000;kg=1;s=1;
%% Given conditions. 
%  For the final assignment take the ones from the specific case you are supposed to do.                  
v1=200; Tamb=300; P3overP2=8; Pamb=100*kPa; mfurate=0.58*kg/s; AF=170.35; % Group 184 settings
cFuel='H2'; % Pick H2 as the fuel
%% Select species for the case at hand
iSp = myfind({Sp.Name},{cFuel,'O2','CO2','H2O','N2'});                      % Find indexes of these species
SpS=Sp(iSp);                                                                % Subselection of the database in the order according to {'Gasoline','O2','CO2','H2O','N2'}
NSp = length(SpS);
Mi = [SpS.Mass];
%% Air composition
Xair = [0 0.21 0 0 0.79];                                                   % Order is important. Note that these are molefractions
MAir = Xair*Mi';                                                            % Row times Column = inner product 
Yair = Xair.*Mi/MAir;                                                       % Vector. times vector is Matlab's way of making an elementwise multiplication
%% Fuel composition
Yfuel = [1 0 0 0 0];                                                        % Only fuel
%% Range of enthalpies/thermal part of entropy of species
TR = [200:1:3000];NTR=length(TR);
for i=1:NSp                                                                 % Compute properties for all species for temperature range TR 
    hia(:,i) = HNasa(TR,SpS(i));                                            % hia is a NTR by 5 matrix
    sia(:,i) = SNasa(TR,SpS(i));                                            % sia is a NTR by 5 matrix
end
hair_a= Yair*hia';                                                          % Matlab 'inner product': 1x5 times 5xNTR matrix muliplication, 1xNTR resulT -> enthalpy of air for range of T 
sair_a= Yair*sia';                                                          % same but this thermal part of entropy of air for range of T
% whos hia sia hair_a sair_a                                                  % Shows dimensions of arrays on commandline
%% Two methods are presented to 'solve' the conservation equations for the Diffusor
%-------------------------------------------------------------------------
% ----> This part shows the interpolation method
% Bisection is in the next 'cell'
%-------------------------------------------------------------------------
% [1-2] Diffusor :: Example approach using INTERPOLATION
cMethod = 'Interpolation Method';
sPart = 'Diffusor';
T1 = Tamb;
P1 = Pamb;
Rg = Runiv/MAir;
for i=1:NSp
    hi(i)    = HNasa(T1,SpS(i));
end
h1 = Yair*hi';
v2 = 0;
h2 = h1+0.5*v1^2-0.5*v2^2;                                                  % Enhalpy at stage: h2 > h1 due to kinetic energy
T2 = interp1(hair_a,TR,h2);                                                 % Interpolate h2 on h2air_a to approximate T2. Pretty accurate
for i=1:NSp
    hi2(i)    = HNasa(T2,SpS(i));
    si1(i)    = SNasa(T1,SpS(i));
    si2(i)    = SNasa(T2,SpS(i));
end
h2check = Yair*hi2';                                                        % Single value (1x5 times 5x1). Why do I do compute this h2check value? Any ideas?
s1thermal = Yair*si1';
s2thermal = Yair*si2';
lnPr = (s2thermal-s1thermal)/Rg;                                            % ln(P2/P1) = (s2-s1)/Rg , see lecture (s2 are only the temperature integral part of th eentropy)
Pr = exp(lnPr);
P2 = P1*Pr;
S1  = s1thermal - Rg*log(P1/Pref);                                          % Total specific entropy
S2  = s2thermal - Rg*log(P2/Pref);
% Print to screen
fprintf('\n%14s\n',cMethod);
fprintf('Stage  ||%14s        [unit]\n      NR|%9i %9i\n',sPart,1,2);
fprintf('-------------------------------------\n');
fprintf('%8s| %9.2f %9.2f  [K]\n','Temp',T1,T2);
fprintf('%8s| %9.2f %9.2f  [kPa]\n','Press',P1/kPa,P2/kPa);
fprintf('%8s| %9.2f %9.2f  [m/s]\n','v',v1,v2);
fprintf('---  H/S    -------------------------\n');
fprintf('%8s| %9.2f %9.2f  [kJ/kg]\n','h',h1/kJ,h2/kJ);
fprintf('%8s| %9.2f %9.2f  [kJ/kg/K]\n','Total S',S1/kJ,S2/kJ);
T2int = T2;
%% Here starts your part (compressor,combustor,turbine and nozzle). ...
%% Compressor (2 -> 3)
s2Part = 'Compressor';
% Assume isentropic compression & assume the air velocity is negligible.
P3 = P2 * (P3overP2);
v3 = v2;
S3 = S2;

% Obtain T3 by interpolation
s3thermal = S3 + Rg*log(P3 / Pref);
T3 = interp1(sair_a, TR, s3thermal);

% Obtain the Enthalpy
for i = 1:NSp
    hi3(i) = HNasa(T3,SpS(i));
end
h3 = Yair * hi3';

% Print to screen
fprintf('\n%14s\n',cMethod);
fprintf('Stage  ||%14s        [unit]\n      NR|%9i %9i\n',s2Part,2,3);
fprintf('-------------------------------------\n');
fprintf('%8s| %9.2f %9.2f  [K]\n','Temp',T2,T3);
fprintf('%8s| %9.2f %9.2f  [kPa]\n','Press',P2/kPa,P3/kPa);
fprintf('%8s| %9.2f %9.2f  [m/s]\n','v',v2,v3);
fprintf('---  H/S    -------------------------\n');
fprintf('%8s| %9.2f %9.2f  [kJ/kg]\n','h',h2/kJ,h3/kJ);
fprintf('%8s| %9.2f %9.2f  [kJ/kg/K]\n','Total S',S2/kJ,S3/kJ);
%% Combustion (3 -> 4)
s3Part = 'Combustor';
% Assume isobaric combustion & constant volume explosion. 
% Assume the air velocity is negligible.
v4 = v3;
P4 = P3; 

% Mass flow rates in the combustor (Stage 3)
mdot_air_3 = AF * mfurate;                 
mdot_total_3 = mfurate + mdot_air_3;      
mdot_N2_3 = Yair(5) * mdot_air_3;          
mdot_O2_3 = Yair(2) * mdot_air_3;         
    
% Mass fractions of the fuel-air mixture before combustion at Stage 3
Y_mixture_3 = [mfurate/mdot_total_3, mdot_O2_3/mdot_total_3, 0, 0, mdot_N2_3/mdot_total_3];
M_mixture_3 = 1 /sum(Y_mixture_3 ./ Mi); % Molar mass of the mixture       
R_g3 = Runiv / M_mixture_3; 

% Molar flows before combustion
ndot_fuel_3 = mfurate / Mi(1);  
ndot_air_3 = mdot_air_3 / MAir; 
ndot_O2_3 = ndot_air_3 * Xair(2); 
ndot_N2_3 = ndot_air_3 * Xair(5); 

% Stoichiometric reaction: H2 + 0.5 O2 -> H2O
ndot_O2_consumed = 0.5 * ndot_fuel_3;
ndot_H2O_4 = ndot_fuel_3;             
ndot_O2_4 = ndot_O2_3 - ndot_O2_consumed;  
ndot_N2_4 = ndot_N2_3;  % N2 remains unchanged (inert gas)

% Stoichiometric air-fuel ratio calculation
AFR_stoichiometric = (0.5 * Mi(2) + (0.5 * Xair(5) / Xair(2)) * Mi(5)) / Mi(1);

% Equivalence ratio
EquivalenceRatio = AFR_stoichiometric / AF;

% Total molar flow after combustion 
ndot_total_4 = ndot_H2O_4 + ndot_O2_4 + ndot_N2_4;

% Mass flow rates after combustion
mdot_H2O_4 = ndot_H2O_4 * Mi(4);       
mdot_O2_4 = ndot_O2_4 * Mi(2);         
mdot_N2_4 = ndot_N2_4 * Mi(5);         
mdot_total_4 = mdot_H2O_4 + mdot_O2_4 + mdot_N2_4; 

% Mass fractions after combustion
Y_mixture_4 = [0, mdot_O2_4/mdot_total_4, 0, mdot_H2O_4/mdot_total_4, mdot_N2_4/mdot_total_4];

% Calculate the new Specific Gas Constant
M_mixture_4 = 1 /sum(Y_mixture_4 ./ Mi); % Molar mass of the mixture       
R_g4 = Runiv / M_mixture_4; 

% Calculating Internal Energies using NASA Polynomials
for i = 1:NSp  
    ui3(i) = UNasa(T3, SpS(i)); 
end

u3 = ui3 * Y_mixture_3'; 
u4 = u3; % Constant volume explosion

for i = 1:NSp 
    uia(:,i) = UNasa(TR, SpS(i)); 
end
umixture_a = Y_mixture_4 * uia'; 

% Interpolate to find T4
T4 = interp1(umixture_a, TR, u4);

% Find Enthalpy and Entropy post-combustion
for i = 1 : NSp
    hi4(i) = HNasa(T4, SpS(i));
    si4(i) = SNasa(T4, SpS(i));
end
h4 =  Y_mixture_4 * hi4';
S4 =  Y_mixture_4 * si4';

% Print to screen
fprintf('\n%14s\n',cMethod);
fprintf('Stage  ||%14s        [unit]\n      NR|%9i %9i\n',s3Part,2,3);
fprintf('-------------------------------------\n');
fprintf('%8s| %9.2f %9.2f  [K]\n','Temp',T3,T4);
fprintf('%8s| %9.2f %9.2f  [kPa]\n','Press',P3/kPa,P4/kPa);
fprintf('%8s| %9.2f %9.2f  [m/s]\n','v',v3,v4);
fprintf('---  H/S    -------------------------\n');
fprintf('%8s| %9.2f %9.2f  [kJ/kg]\n','h',h3/kJ,h4/kJ);
fprintf('%8s| %9.2f %9.2f  [kJ/kg/K]\n','Total S',S3/kJ,S4/kJ);

%% Turbine (4 -> 5)
s4Part = 'Turbine';
% Assume isentropic process & negligible air velocity
S5 = S4;
v5 = v4;

% Calculating enthalpy based on the first law
h5 = h4 - (mdot_total_3/mdot_total_4) * (h3 - h2);

% Enthalpy and Entropy ranges for the mixture
hmixture_a = Y_mixture_4 * hia';
smixture_a = Y_mixture_4 * sia';

% Interpolation to find T5
T5 = interp1(hmixture_a, TR, h5);

% Find the pressure 
for i = 1 : NSp
    si5(i) = SNasa(T5,SpS(i));
end
s5thermal = Y_mixture_4 * si5';
s4thermal = Y_mixture_4 * si4';
P5 = P4 * exp((s5thermal-s4thermal)/R_g4);

fprintf('\n%14s\n',cMethod);
fprintf('Stage  ||%14s        [unit]\n      NR|%9i %9i\n',s4Part,4,5);
fprintf('-------------------------------------\n');
fprintf('%8s| %9.2f %9.2f  [K]\n','Temp',T4,T5);
fprintf('%8s| %9.2f %9.2f  [kPa]\n','Press',P4/kPa,P5/kPa);
fprintf('%8s| %9.2f %9.2f  [m/s]\n','v',v4,v5);
fprintf('---  H/S    -------------------------\n');
fprintf('%8s| %9.2f %9.2f  [kJ/kg]\n','h',h4/kJ,h5/kJ);
fprintf('%8s| %9.2f %9.2f  [kJ/kg/K]\n','Total S',S4/kJ,S5/kJ);
%% Nozzle (5 -> 6)
s5part = 'Nozzle';
% Assume isentropic process
S6 = S5;
P6 = Pamb;

% Obtain T6 using interpolation
s6thermal = s5thermal + R_g4 * log(P6/P5);
T6 = interp1(smixture_a, TR, s6thermal);

% Calculate enthalpy and velocity
for i = 1:NSp
    hi6(i) = HNasa(T6, SpS(i));
end
h6 = Y_mixture_4 * hi6';
v6 = sqrt(2 * (h5 - h6));

fprintf('\n%14s\n',cMethod);
fprintf('Stage  ||%14s        [unit]\n      NR|%9i %9i\n',s5part,5,6);
fprintf('-------------------------------------\n');
fprintf('%8s| %9.2f %9.2f  [K]\n','Temp',T5,T6);
fprintf('%8s| %9.2f %9.2f  [kPa]\n','Press',P5/kPa,P6/kPa);
fprintf('%8s| %9.2f %9.2f  [m/s]\n','v',v5,v6);
fprintf('---  H/S    -------------------------\n');
fprintf('%8s| %9.2f %9.2f  [kJ/kg]\n','h',h5/kJ,h6/kJ);
fprintf('%8s| %9.2f %9.2f  [kJ/kg/K]\n','Total S',S5/kJ,S6/kJ);

%% Final summary for all stages
fprintf('\n%67s\n', sprintf('Final Summary (Used the %s)', cMethod)); 
fprintf(' Stage ||   Temp [K]  |  Press [kPa] |  Velocity [m/s] |  h [kJ/kg]   |  S [kJ/kg/K]\n');
fprintf('------------------------------------------------------------------------------------\n');
fprintf('   1   ||  %9.2f  |  %9.2f   |   %9.2f     |  %9.2f   |   %9.2f \n', T1, P1/kPa, v1, h1/kJ, S1/kJ);
fprintf('   2   ||  %9.2f  |  %9.2f   |   %9.2f     |  %9.2f   |   %9.2f \n', T2, P2/kPa, v2, h2/kJ, S2/kJ);
fprintf('   3   ||  %9.2f  |  %9.2f   |   %9.2f     |  %9.2f   |   %9.2f \n', T3, P3/kPa, v3, h3/kJ, S3/kJ);
fprintf('   4   ||  %9.2f  |  %9.2f   |   %9.2f     |  %9.2f   |   %9.2f \n', T4, P4/kPa, v4, h4/kJ, S4/kJ);
fprintf('   5   ||  %9.2f  |  %9.2f   |   %9.2f     |  %9.2f   |   %9.2f \n', T5, P5/kPa, v5, h5/kJ, S5/kJ);
fprintf('   6   ||  %9.2f  |  %9.2f   |   %9.2f     |  %9.2f   |   %9.2f \n', T6, P6/kPa, v6, h6/kJ, S6/kJ);
