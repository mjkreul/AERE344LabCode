%% AER E 344 Lab 3 Calculations
% Section 1
% Group 1
clear, clc
pause on
%% Constants
P_atm = 99000; %pascals
R = 287; %temp constant J/kg*K
rho_water = 997; %density of water kg/m^3
Setra1C = 746.52; %Pa/volt
Setra2C = 248.84; %Pa/volt
inH2OtoPa = 248.84; %Pa
inHgtoPa = 3386.39; %Pa
T = 22 + 273.15; %Kelvin
v = 10; % m/s^2


%% Parsing the file values
% I moved all of the files into the working directory for easier access
calibrationfiles = dir(fullfile('.', '*iw.txt'));
testingfiles = dir(fullfile('.', '*in.txt'));

% Display the names
calibrationfiles.name;
testingfiles.name;

% get the number of files for the calibration and other test
[lengthCal, junk] = size(calibrationfiles);
[lengthTest, junk] = size(testingfiles);

% Get the values from the given files
[calvolt, calvolttime, calH2OHeight] = parseTestFiles(lengthCal, calibrationfiles);

[testvolt, testvolttime, testDist] = parseTestFiles(lengthTest, testingfiles);


%% Converting all of the values
calPressures = calH2OHeight * inH2OtoPa;
% rho = 

% q = dynamicP(v, rho);





