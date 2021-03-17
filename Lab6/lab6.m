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

for i = 2:46
    a1(i-1) = mean(ai1(:,i));
    a2(i-1) = mean(ai2(:,i));
    a3(i-1) = mean(ai3(:,i));
    a4(i-1) = mean(ai4(:,i));
    a5(i-1) = mean(ai5(:,i));
    a6(i-1) = mean(ai6(:,i));
end


P_A_P_ETable(1,1) = a1(:,44);
P_A_P_ETable(1,2) = a1(:,45);
P_A_P_ETable(2,1) = a2(:,44);
P_A_P_ETable(2,2) = a2(:,45);
P_A_P_ETable(3,1) = a3(:,44);
P_A_P_ETable(3,2) = a3(:,45);
P_A_P_ETable(4,1) = a4(:,44);
P_A_P_ETable(4,2) = a4(:,45);
P_A_P_ETable(5,1) = a5(:,44);
P_A_P_ETable(5,2) = a5(:,45);
P_A_P_ETable(6,1) = a6(:,44);
P_A_P_ETable(6,2) = a6(:,45);


PA = P_A_P_ETable(:,1);
PE = P_A_P_ETable(:,2);

p(:,1) = a1(1:43)';
p(:,2) = a2(1:43)';
p(:,3) = a3(1:43)';
p(:,4) = a4(1:43)';
p(:,5) = a5(1:43)';
p(:,6) = a6(1:43)';


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

%% Get Cp
% C_P = average pressure of each point/dynamic pressure at that air speed
% C_P_theoretical = 1-4*sin(theta)*sin(theta)
C_P = zeros(43, numFiles);
for i = 1:numFiles
    for j = 1:43
        C_P(j, i) = p(j,i)/q(i);
    end
end

% C_P_theoretical = 1-4*sin((-180:5:180)*pi/180).*sin((-180:5:180)*pi/180);

%% Plotting Cp values

angle = [-4, 0, 4, 8, 10, 12,];
for i = 1:numFiles
    figure(i)
    plot(1:43,C_P(:,i))
%     hold on
%     plot(21:43,C_P(21:43,i),'r')
%     set(gca,'Ydir','reverse')
%     hold off
%     legend('Lower Surface','Upper Surface')
    xlabel('x/c')
    ylabel('Coefficient of Pressure, C_P')
    title(['Angle of Attack ', num2str(angle(i)), char(176)])
%     exportgraphics(gcf,sprintf('outputfiles/f%d.png',i+1),'Resolution',300)
end

%% Calculating Lift, Drag, and Moment

%preallocating variables
Pseg = zeros(43:numFiles);
dx = zeros(numFiles);
dy = zeros(numFiles);
xseg = zeros(numFiles);
yseg = zeros(numFiles);
N = zeros(numFiles);
A = zeros(numFiles);
M = zeros(numFiles);

%calculating p_i+1/2, delta X, delta Y, x_i+1/2, and y_i+1/2
for j=1:numFiles
    for i=1:43
        if i == 43 %edge value as described in lab instructions
            Pseg(i,j) = (p(1,j)+p(length(p),j)/2);
            dx(i) = (x(1)-x(43))/2;
            dy(i) = (y(1)-y(43))/2;
            xseg(i) = (x(1)+x(43))/2;
            yseg(i) = (y(1)+y(43))/2;
        else
            Pseg(i,j) = (p(i,j)+p(i+1,j))/2;
            dx(i) = (-x(i)+x(i+1));
            dy(i) = (-y(i)+y(i+1));
            xseg(i) = (x(i)+x(i+1))/2;
            yseg(i) = (y(i)+y(i+1))/2;
        end
        
    end
end

for i = 1:numFiles
    
    %summing N A and M for each aoa
    N(i) = 0;
    A(i) = 0;
    M(i) = 0;
    for j = 1:43
        N(i) = Pseg(j,i)*dx(j)+N(i);
        A(i) = -Pseg(j,i)*dy(j)+A(i);
        M(i) = -(Pseg(j,i)*dx(j))*xseg(j)-Pseg(j,i)*dy(j)*yseg(j)+M(i);
    end
    
    %calculating lift and drag forces at each aoa
    L(i) = N(i)*cosd(angle(i))-A(i)*sind(angle(i)); %lift
    D(i) = N(i)*sind(angle(i))+A(i)*cosd(angle(i)); %drag
end

for i = 1:numFiles %calculating coefficients
%     Cl(i) = L(i)/q(i); %coefficient of lift
    Cd(i) = D(i)/q(i); %coefficient of drag
%     Cm(i) = M(i)/q(i); %coefficient of moment about LE
end
%% Plotting Lift, Drag, Moment

% figure(9)
% plot(angle,Cl)
% title('Coefficient of Lift vs Angle of Attack')
% xlabel(['Angle of Attack (', char(176), ')'])
% ylabel('C_l')
% exportgraphics(gcf,sprintf('outputfiles/f9.png'),'Resolution',300)
% 
% figure(10)
% plot(angle,Cd)
% title('Coefficient of Drag vs Angle of Attack')
% xlabel(['Angle of Attack (', char(176), ')'])
% ylabel('C_d')
% exportgraphics(gcf,sprintf('outputfiles/f10.png'),'Resolution',300)
% 
% figure(11)
% plot(angle,Cm)
% title('Coefficient of Moment vs Angle of Attack')
% xlabel(['Angle of Attack (', char(176), ')'])
% ylabel('C_m_,_  _L_E')
% exportgraphics(gcf,sprintf('outputfiles/f11.png'),'Resolution',300)

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

