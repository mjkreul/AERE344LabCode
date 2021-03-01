%% AER E 344 Lab 4 Calculations
% Section 1
% Group 1
clear, clc
%% Constants
P_atm = 99000; %pascals
R = 287; %temp constant J/kg*K
rho_water = 997; %density of water kg/m^3
T = 22 + 273.15; %Kelvin


%% Parsing the file values
% I moved all of the files into the working directory for easier access
CSVFiles = dir(fullfile('.', '*.csv'));
% testingfiles = dir(fullfile('.', '*in.txt'));


% Display the names
CSVFiles.name;


% get the number of files for the calibration and other test
[numFiles, junk] = size(CSVFiles);

pressureAvg = [];

for i = 1:numFiles
    importedCSV = readtable(CSVFiles(i).name);
    first = importedCSV(:, 2:17);
    second = importedCSV(:, 35:36);
    combined = [first, second];
    [rows, cols] = size(combined);
    for j = 1:cols
        pressureAvg(:,i) = mean(combined(:,j));
    end
end

% Get the values from the given files
% [calvolt, calvolttime, calH2OHeight] = parseTestFilesLab4(lengthCal, files);




%% Converting all of the values
% calPressures = calH2OHeight * inH2OtoPa;
% rho = P_atm/(R*T);
% 
% %% Calibration Data
% 
% %extracting average values of voltage
% calvoltavg = zeros(10,1);
% for i = 1:10
%     calvoltavg(i) = mean(calvolt(:,i));
% end
% 
% %linear regression for calculating test data
% p = polyfit(calvoltavg,calPressures',1);
% C = p(1);
% v0 = p(2);
% 
% %plotting voltage vs pressure
% figure(1)
% plot(calvoltavg,calPressures)
% xlabel('Sensor Voltage')
% ylabel('Pressure Applied to Sensor (Pa)')
% title('Voltage vs Applied Pressure')
% 
% %% Test Data
% 
% %assign empty matrix
% avgtestvals = zeros(25,2);
% 
% %calculating average values for voltage
% for i = 1:13
%     avgtestvals(i,1) = (i*0.5)-0.5;
%     avgtestvals(i,2) = mean(testvolt(:,i));
% end
% 
% %mirroring average values for 2nd half of section
% for i = 14:25
%     avgtestvals(i,1) = (i*0.5)-0.5;
%     revref = 26 - i;
%     avgtestvals(i,2) = avgtestvals(revref,2);
% end
% 
% %calculating voltage -> dynamic pressure -> velocity
% testpressures = (avgtestvals(:,2)*C)+v0;
% testvelocities = sqrt(2*rho*testpressures);
% 
% %plotting dynamic pressure cross section
% figure(2)
% plot(testpressures,avgtestvals(:,1))
% xlabel('Dynamic Pressure (Pa)')
% ylabel('Distance From Wall (in)')
% title('Tunnel Dynamic Pressure Cross-Section')
% 
% %plotting velocity cross section
% figure(3)
% plot(testvelocities,avgtestvals(:,1))
% xlabel('Velocity (m/s)')
% ylabel('Distance From Wall (in)')
% title('Tunnel Velocity Cross-Section')