%% AER E 344 Lab 4 Calculations
% Section 1
% Group 1
clear, clc
%% Constants
P_atm = 99000; %pascals
R = 287; %temp constant J/kg*K
rho_water = 997; %density of water kg/m^3
T = 22 + 273.15; %Kelvin
rho = P_atm/(R*T); %density of air kg/m^3
K = 1.14393; % may need to change this
mu = 1.825*10^(-5); % dynamic viscosity of air
dia = 0.7*0.0254; % diameter (2R) of the cyllinder

%% Parsing the file values and averaging them
% I moved all of the files into the working directory for easier access
% NOTE: matlab does NOT read "pwd" as the directory in which this file is
% saved, instead it reads as the CURRENT FOLDER THAT IS OPEN IN THE CURRENT
% FOLDER TAB!!!! This is so dumb... who the fuck creates a programming
% language that does this???
CSVFiles = dir(fullfile(pwd, '*.csv'));

% get the number of files
[numFiles, junk] = size(CSVFiles);

% initialize an array to hold all of the averages of the pressures
pressureAvg = [];
P_A_P_EAvg = [];
for i = 1:numFiles
    % import the table
    importedCSV = readtable(CSVFiles(i).name);
    % grab the "first" 18 columns (there is a single column that is the ms)
    pressureTable = importedCSV(:, 2:19); 
    % then change them from a table to an array... I fucking hate matlab...
    pressure = pressureTable{:,:};
    % get the average of each
    averagedP = mean(pressure);
    % then transpose it to make it a column
    pressureAvg(:,i) = averagedP';
    % Get PA and PE and do the same thing to them
    P_A_P_ETable = importedCSV(:,20:21);
    P_A_P_E = P_A_P_ETable{:,:};
    averagedP_A_P_E = mean(P_A_P_E);
    P_A_P_EAvg(:,i) = averagedP_A_P_E';
end

% Write the values to a csv file to import into LaTeX doc
writetable(array2table(pressureAvg), 'outputFiles/avgP.csv', 'Delimiter', ';');
writetable(array2table(P_A_P_EAvg), 'outputFiles/avgP_A_P_E.csv',...
    'Delimiter', ';');

%% Get Dynamic Pressure and Velocity

% init q array
q = [];
% v = sqrt(2(P_A - P_E)/rho*k)
v = [];
for i = 1:numFiles
    P_AmP_E = P_A_P_EAvg(1, i) - P_A_P_EAvg(2, i);
    q(i) = K*P_AmP_E;
    v(i) = sqrt((2*P_AmP_E)/(rho * K));
end

%% Get Cp
% C_P = average pressure of each point/dynamic pressure at that air speed
% C_P_theoretical = 1-4*sin(theta)*sin(theta)
[rows, cols] = size(pressureAvg);
C_P = zeros(rows, cols);
for i = 1:cols
    for j = 1:rows
        C_P(j, i) = pressureAvg(j,i)/q(i);
    end
end
C_P_theoretical = 1-4*sin((-180:5:180)*pi/180).*sin((-180:5:180)*pi/180);

%% Get D 

%integral -pi to pi, P*cos(theta)*R*dtheta
% = -deltaTheta/2 * sum_{i = 1}^{n} C_pi*cos(theta)
% start theta from 360 and subtract by 20 until theta = 0
% 
    
C_D = [];

for i = 1:cols
    theta = 360;
    for j = 1:rows
        sum = C_P(j,i)* cos(theta*(pi/180));
        theta = theta - 20;
    end
    C_D(i) = sum*(-(20*pi/180)/2);
end
    
 
%% Get C_D as function of Re
% Re = (rho*V_inf*D)/mu
% D = C_D*rho*V_inf^2*2*R*0.5
D = [];
Re = [];
for i = 1:cols
    D(i) = C_D(i)*rho*(v(i)^2)*dia*0.5;
    Re(i) = (rho*v(i)*D(i))/mu;
end


%% Plots

for plotNum = 1:cols
    figure(plotNum);
    plot(C_P(:,plotNum));    
end

figure(plotNum);
plot(C_D, Re)  
xlabel("Coefficient of Drag");
ylabel("Reynolds Number");
title("Coefficient of Drag as a function of Reynolds Number");









