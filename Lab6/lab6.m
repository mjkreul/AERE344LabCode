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
inH2OtoPa = 248.84; %Pa

%% Opening and parsing data

ai1 = load('n4.csv');
ai2 = load('0.csv');
ai3 = load('4.csv');
ai4 = load('6.csv');
ai5 = load('8.csv');
ai6 = load('10.csv');
ai7 = load('12.csv');

for i = 2:44
    a1(i-1) = mean(ai1(:,i));
    a2(i-1) = mean(ai2(:,i));
    a3(i-1) = mean(ai3(:,i));
    a4(i-1) = mean(ai4(:,i));
    a5(i-1) = mean(ai5(:,i));
    a6(i-1) = mean(ai6(:,i));
    a7(i-1) = mean(ai7(:,i));
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
P_A_P_ETable(7,1) = a7(:,42);
P_A_P_ETable(7,2) = a7(:,43);

PA = P_A_P_ETable(:,1);
PE = P_A_P_ETable(:,2);

p(:,1) = a1(1:41)';
p(:,2) = a2(1:41)';
p(:,3) = a3(1:41)';
p(:,4) = a4(1:41)';
p(:,5) = a5(1:41)';
p(:,6) = a6(1:41)';
p(:,7) = a7(1:41)';

numFiles = 7;

%% Get Dynamic Pressure and Velocity

% init q array
q = [];
% v = sqrt(2(P_A - P_E)/rho*k)
v = [];
for i = 1:numFiles
    % get pa - pe
    P_AmP_E = PA(i) - PE(i);
    % dynamic pressure = k*pa - pe
    q(i) = K*P_AmP_E;
    % velocity of airflow = sqrt(2*q/rho)
    v(i) = sqrt((2*K*P_AmP_E)/(rho));
end

%% Get Cp and velocity of each place
% C_P = average pressure of each point/dynamic pressure at that air speed
% C_P_theoretical = 1-4*sin(theta)*sin(theta)
C_P = zeros(41, numFiles);
u = zeros(41, numFiles);
for i = 1:numFiles
    for j = 1:41
        % cp = pressure/dynamic pressure
        C_P(j, i) = p(j,i)/q(i);
        % v = sqrt(2*(dynamic pressure at that point)/density)
        % dynamic pressure = ptotal-pstatic
        % pstatic = pe
        % ptotal = pressure at that point
        u(j,i) = sqrt(2*(p(j,i) - PE(i))/rho);
    end
end

%% Plotting Cp values

angle = [-4, 0, 4, 6, 8, 10, 12];
for i = 1:numFiles
    figure(i)
    % plot figures against
    plot(C_P(:,i), 0:2:80)
    yticks(0:4:80)
    grid on
    ylabel('Rake Positions (mm)')
    xlabel('Coefficient of Pressure, C_P')
    title(['Angle of Attack ', num2str(angle(i)), char(176)])
    exportgraphics(gcf,sprintf('outputfiles/f%d.png',i),'Resolution',300)
end

% superimposed figures
figure(numFiles+1)
hold on
for i = 1:numFiles
    plot(C_P(:,i), 0:2:80)
end
grid on
hold off
yticks(0:4:80)
ylabel('Rake Positions (mm)')
xlabel('Coefficient of Pressure, C_P')
title('All AoA C_P Superimposed');
superimpcplgd =legend(['-4',char(176), ' AoA'],['0',char(176), ' AoA'],...
    ['4',char(176), ' AoA'], ['6',char(176), ' AoA'],['8',char(176), ' AoA'],...
    ['10',char(176), ' AoA'], ['12',char(176), ' AoA']);
superimpcplgd.Location = 'Northwest';
exportgraphics(gcf,sprintf('outputfiles/f10.png'),'Resolution',300)

%% Calculating Drag

Cd =[];

%  cd = (2/c)\integral ((u(y)/u_inf)*(1 - u(y)/u_inf))dy
% represented by trapezoidal rule of sum i = 2 to n
% h*(f(x)_{i-1} + f(x)_{i})/2 
% where f(x) = ((u(y)/u_inf)*(1 - u(y)/u_inf))*(2/c)
for i = 1:numFiles %calculating coefficients
    %coefficient of drag
    dragSum = 0;
    trapTerm = 0;
    %trapezoidal rule for getting the coefficient of drag
    for j = 2:40
        trapTerm = trapTerm + ((u(j,i)/v(i))*(1 - u(j,i)/v(i)))*(2/c);
        %function of the left term f(x)_{i-1}
%         trapTerm1 = ((u(j-1,i)/v(i))*(1 - u(j-1,i)/v(i)))*(2/c);
%         %function of the right term f(x)_{i}
%         trapTerm2 = ((u(j,i)/v(i))*(1 - u(j,i)/v(i)))*(2/c);
        %apply trapezoidal rule
%         tempIntegral = ((trapTerm1 + trapTerm2)/2)*0.002;
        %add value to sum
%         dragSum = dragSum + tempIntegral;
    end
    trapTerm = 0.002*(trapTerm + (((u(41,i)/v(i))*(1 - u(41,i)/v(i))) + ((u(1,i)/v(i))*(1 - u(1,i)/v(i))))/2);
    Cd(i) = trapTerm;
end

%% Plotting Drag
figure(numFiles+2)
plot(angle,Cd)
grid on
title('Coefficient of Drag vs Angle of Attack')
xlabel(['Angle of Attack (', char(176), ')'])
ylabel('C_d')
exportgraphics(gcf,sprintf('outputfiles/f%d.png',i+1),'Resolution',300)

%% Getting Re
% for i =1:numFiles
%     Re(i) = (rho*v(i)*c*abs(sind(angle(i))))/mu;
% end


%% Hotwire Calibration
calibrationfiles = dir(fullfile('.', '*iw.txt'));

calibrationfiles.name;
% get the number of files for the calibration and other test
[lengthCal, junk] = size(calibrationfiles);

% Get the values from the given files
[calvolt, calvolttime, calH2OHeight] = parseTestFiles(lengthCal, calibrationfiles);

%% Changing H2O values to pascals and getting velocity

qCal = zeros(1, lengthCal);
vCal = zeros(1, lengthCal);
for i = 1:lengthCal
    qCal(i) = calH2OHeight(i)*inH2OtoPa;
    vCal(i) = sqrt((qCal(i)*2)/rho);
end

%% Getting average of voltage for the hotwire
voltAvg = zeros(1, lengthCal);

for i = 1:lengthCal
    voltAvg(i) = mean(calvolt(:,i));
end
%% Getting polyfit line

coeff = polyfit(vCal, voltAvg, 4); 

xVals = 0:(14/100):14;

yVals = polyval(coeff, xVals);

figure(numFiles+3)

plot(voltAvg, vCal, 'o');
hold on
plot(yVals, xVals)% I know this looks wrong but trust me its correct
grid on
hold off
ylabel('Velocity (m/s)');
xlabel('Voltage (V)');
title('Velocity vs. Voltage for Hot Wire');
exportgraphics(gcf,sprintf('outputfiles/f%d.png',numFiles+2),'Resolution',300)


%% Sending calculated data to csv files
% writetable(array2table(round(p,4)), 'outputFiles/avgP.csv', 'Delimiter', ',');
% writetable(array2table(round([PA(:), PE(:)],4)), 'outputFiles/avgP_A_P_E.csv',...
%     'Delimiter', ',');
% 
% writetable(array2table(round(C_P, 4)), 'outputFiles/C_P.csv', 'Delimiter', ',');
writetable(array2table(round(Cd', 4)), 'outputFiles/C_D.csv', 'Delimiter', ',');
% 
% writetable(array2table(round(u, 4)), 'outputFiles/u.csv', 'Delimiter', ',');
% 
% writetable(array2table(round(Cd', 4)), 'outputFiles/C_D.csv', 'Delimiter', ',');

