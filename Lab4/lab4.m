%% AER E 344 Lab 4 Calculations
% Section 1
% Group 1
clear, clc
%% Constants
P_atm = 99000; %pascals
R = 287; %temp constant J/kg*K
rho_water = 997; %density of water kg/m^3
T = 22 + 273.15; %Kelvin
rho = P_atm*R*T; %
K = 1.1; % may need to change this

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
    % grab the "first" 16 columns (there is a single column that is the ms)
    first = importedCSV(:, 2:17); 
    % grab the "second" 2 columns (there are a lot of values for
    % temperature)
    second = importedCSV(:, 35:36);
    % concatenate them...
    combinedTable = [first, second];
    % then change them from a table to an array... I fucking hate matlab...
    combined = combinedTable{:,:};
    % get the average of each
    averagedP = mean(combined);
    % then transpose it to make it a column
    pressureAvg(:,i) = averagedP';
    
    % Get PA and PE and do the same thing to them
    P_A_P_ETable = importedCSV(:,37:38);
    P_A_P_E = P_A_P_ETable{:,:};
    averagedP_A_P_E = mean(P_A_P_E);
    P_A_P_EAvg(:,i) = averagedP_A_P_E';
end

% Write the values to a csv file to import into LaTeX doc
writetable(array2table(pressureAvg), 'outputFiles/avgP.csv', 'Delimiter', ';');
writetable(array2table(P_A_P_EAvg), 'outputFiles/avgP_A_P_E.csv', 'Delimiter', ';');

%% Get Dynamic Pressure and Velocity

% init q array
q = [];
% v = sqrt(2(P_A - P_E)/rho*k)
v = [];
for i = 1:numFiles
    P_AmP_E = P_A_P_EAvg(1, i) - P_A_P_EAvg(2, i);
    
    q(i) = P_AmP_E;
    v(i) = sqrt((2*P_AmP_E)/(rho * K));
    
end

%% Get Cp

[rows, cols] = size(pressureAvg);
Cp = zeros(rows, cols);
for i = 1:cols
    for j = 1:rows
        Cp(j, i) = pressureAvg(j,i)/q(i);
    end
end
    
    
    
    
    
    








