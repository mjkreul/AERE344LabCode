%% AER E 344 Lab 4 Calculations
% Section 1
% Group 1
clear, clc, clf
pause on
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

CSVFiles.name;

names = [];
% get the number of files
[numFiles, junk] = size(CSVFiles);

% initialize an array to hold all of the averages of the pressures
pressureAvg = [];
P_A_P_EAvg = [];
for i = 1:numFiles
    % import the table
    % get the name from the input struct
    stringname = CSVFiles(i).name;
    % set the name as the file name minus the "_iw.txt" or "_in.txt"
    names(1,i) = str2double(strtok(stringname,"."));
    % open the file
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



%% Get Dynamic Pressure and Velocity
% init q array
q = [];
% v = sqrt(2(P_A - P_E)/rho*k)
v = [];
for i = 1:numFiles
    P_AmP_E = P_A_P_EAvg(1, i) - P_A_P_EAvg(2, i);
    q(i) = K*P_AmP_E;
    v(i) = sqrt((2*q(i))/(rho));
end

%% Get Cp
% C_P = average pressure of each point - P_E for that speed/dynamic 
% pressure at that air speed
% C_P_theoretical = 1-4*sin(theta)*sin(theta)
[rows, cols] = size(pressureAvg);
C_P = zeros(rows, cols);
for i = 1:cols
    for j = 1:rows
        C_P(j, i) = (pressureAvg(j,i) - P_A_P_EAvg(2,i))/q(i);
    end
end
C_P_theoretical = 1-4*sind((-180:0.5:180)).*sind((-180:0.5:180));

%% Get D and C_D
% integral -pi to pi, P*cos(theta)*R*dtheta
% = -deltaTheta/2 * sum_{i = 1}^{n} C_pi*cos(theta)
% start theta from 360 and subtract by 20 until theta = 0
% 
C_D = [];
for i = 1:cols
    theta = 0;
    sum = 0;
    for j = 1:rows
        sum = sum + C_P(j,i)* cosd(theta);
        theta = theta + 20;
    end
    C_D(i) = sum*(-(pi/18));
end
 
%% Get Re
% Re = (rho*V_inf*D)/mu
% D = C_D*rho*V_inf^2*2*R*0.5
Re = [];
for i = 1:cols
    Re(i) = (rho*v(i)*dia)/mu;
end

%% Tables 
% Write the values to a csv file to import into LaTeX doc
writetable(array2table(round(pressureAvg,4)), 'outputFiles/avgP.csv', 'Delimiter', ',');
writetable(array2table(round(P_A_P_EAvg,4)), 'outputFiles/avgP_A_P_E.csv',...
    'Delimiter', ',');
writetable(array2table(round(C_P, 4)), 'outputFiles/C_P.csv', 'Delimiter', ',');
writetable(array2table(round(C_D,4)), 'outputFiles/C_D.csv', 'Delimiter', ',');
writetable(array2table(round(Re,4)), 'outputFiles/Re.csv', 'Delimiter', ',');

%% Plots

% Plotting all of the C_P for each intake 
for plotNum = 1:cols
    figure(plotNum);
    plot(-180:20:160, C_P(:,plotNum),'-o', 'LineWidth', 2.0);  
    grid on
    xlabel("Air Inlet Angles (degrees)");
    ylabel("Coefficient of Pressure");
    title(sprintf("Coefficient of pressure at %dHz",names(1,plotNum)) );
    set(gca,'FontSize',14)
    % Output figure to a png file (Note this only works in matlab2020a and
    % newer...
    exportgraphics(gcf,sprintf('outputfiles/f%d.png',plotNum),'Resolution',300)
%     pause(1);
    
end

% Plotting the theoretical C_P
figure(plotNum+1)
plot(-180:0.5:180, C_P_theoretical, 'LineWidth', 2.0);
grid on
xlabel("Angle (Degrees)")
ylabel("Theoretical Coefficient of Pressure");
title("Theoretical Coefficient of Pressure");
set(gca,'FontSize',14)
exportgraphics(gcf,sprintf('outputfiles/f%d.png',plotNum+1),'Resolution',300)

% Plotting C_D and Re
figure(plotNum+2);
plot(Re, C_D, '-o', 'LineWidth', 2.0)  
grid on
ylabel("Coefficient of Drag");
xlabel("Reynolds Number");
title("C_D as a function of Reynolds Number");
set(gca,'FontSize',14)
exportgraphics(gcf,sprintf('outputfiles/f%d.png',plotNum+2),'Resolution',300)








