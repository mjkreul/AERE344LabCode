%% AER E 344 Lab 6 Calculations
% Section 1
% Group 1
clear, clc, clf
pause on
clear, clc, clf
close all

%% Constants
P_atm = 99000; %pascals
R = 287; %temp constant J/kg*K
rho_water = 997; %density of water kg/m^3
T = 21 + 273.15; %Kelvin
rho = P_atm/(R*T); %density of air kg/m^3
K = 1.14393; % may need to change this
mu = 1.825*10^(-5); % dynamic viscosity of air
c = 0.101; %chord length (m)

%% Opening and parsing data

ai1 = load('n4.csv');
ai2 = load('0.csv');
ai3 = load('4.csv');
ai4 = load('8.csv');
ai5 = load('10.csv');
ai6 = load('12.csv');

% coords = load('coordinates.txt');
% x = coords(:,1);
% y = coords(:,2);

for i = 2:44
    a1(i-1) = mean(ai1(:,i));
    a2(i-1) = mean(ai2(:,i));
    a3(i-1) = mean(ai3(:,i));
    a4(i-1) = mean(ai4(:,i));
    a5(i-1) = mean(ai5(:,i));
    a6(i-1) = mean(ai6(:,i));
end


P_A_P_ETable(1,1) = a1(:,42);
P_A_P_ETable(1,2) = a1(:,43);
P_A_P_ETable(2,1) = a2(:,42);
P_A_P_ETable(2,2) = a2(:,43);
P_A_P_ETable(3,1) = a3(:,42);
P_A_P_ETable(3,2) = a3(:,43);
P_A_P_ETable(4,1) = a4(:,42);
P_A_P_ETable(4,2) = a4(:,43);
P_A_P_ETable(5,1) = a5(:,42);
P_A_P_ETable(5,2) = a5(:,43);
P_A_P_ETable(6,1) = a6(:,42);
P_A_P_ETable(6,2) = a6(:,43);


PA = P_A_P_ETable(:,1);
PE = P_A_P_ETable(:,2);

p(:,1) = a1(1:41)';
p(:,2) = a2(1:41)';
p(:,3) = a3(1:41)';
p(:,4) = a4(1:41)';
p(:,5) = a5(1:41)';
p(:,6) = a6(1:41)';


numFiles = 6;

%% Get Dynamic Pressure and Velocity

% init q array
q = [];
% v = sqrt(2(P_A - P_E)/rho*k)
v = [];
for i = 1:numFiles
    P_AmP_E = PA(i) - PE(i);
    q(i) = K*P_AmP_E;
    v(i) = sqrt((2*P_AmP_E)/(rho * K));
end

%% Get Cp and velocity of each place
% C_P = average pressure of each point/dynamic pressure at that air speed
% C_P_theoretical = 1-4*sin(theta)*sin(theta)
C_P = zeros(41, numFiles);
u = zeros(41, numFiles);
for i = 1:numFiles
    for j = 1:41
        C_P(j, i) = p(j,i)/q(i);
        u(j,i) = sqrt(2*(p(j,i)/q(i) - PE(i))/rho);
    end
end

%% Plotting Cp values

angle = [-4, 0, 4, 8, 10, 12,];
for i = 1:numFiles
    figure(i)
    plot(0:2:80,C_P(:,i))
    xticks(0:4:80)
    xlabel('Rake Positions (mm)')
    ylabel('Coefficient of Pressure, C_P')
    title(['Angle of Attack ', num2str(angle(i)), char(176)])
%     exportgraphics(gcf,sprintf('outputfiles/f%d.png',i+1),'Resolution',300)
end

%% Calculating Drag

%preallocating variables

for i = 1:numFiles %calculating coefficients
    %coefficient of drag
    dragSum = 0;
    for j = 1:41
        dragSum = dragSum + ((u(j,i)/v(i))*(1 - (u(j,i)/v(i))))*0.002; 
    end
    Cd(i) = dragSum*(2/c);
end
%% Plotting Drag

% 
figure(10)
plot(angle,Cd)
title('Coefficient of Drag vs Angle of Attack')
xlabel(['Angle of Attack (', char(176), ')'])
ylabel('C_d')


%% Hotwire Calibration



%% Sending calculated data to csv files
% writetable(array2table(round(p,4)), 'outputFiles/avgP.csv', 'Delimiter', ',');
% writetable(array2table(round([PA(:), PE(:)],4)), 'outputFiles/avgP_A_P_E.csv',...
%     'Delimiter', ',');
% 
% writetable(array2table(round(C_P, 4)), 'outputFiles/C_P.csv', 'Delimiter', ',');
% writetable(array2table(round(Cd', 4)), 'outputFiles/C_D.csv', 'Delimiter', ',');
% 
% writetable(array2table(round(Cl', 4)), 'outputFiles/C_l.csv', 'Delimiter', ',');
%  
% 
% writetable(array2table(round(Cm', 4)), 'outputFiles/C_m.csv', 'Delimiter', ',');

