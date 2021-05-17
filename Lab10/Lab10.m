
%% AER E 344 Lab 10
% Section 1
% Group 1
% @author: Matt Kreul, Andrew Fung
clear,  clc
%% Constants
P_atm = 1010*0.0145038; %psi
R = 287; %temp constant J/kg*K
rho_water = 997; %density of water kg/m^3
T = 25 + 273.15; %Kelvin
rho = P_atm/(R*T); %density of air kg/m^3
K = 1.14393; % may need to change this
mu = 1.825*10^(-5); % dynamic viscosity of air
c = 0.101; %chord length (m)
inH2OtoPa = 248.84; %Pa
gamma = 1.4;
throatPos = 5;

%% Inlets

inlPos = [-4.00, -1.50 , -0.30, -0.18, 0.00, 0.15, 0.30, 0.45, 0.60, ...
    0.75, 0.90, 1.05, 1.20, 1.35, 1.45]; %inches
inlArea = [0.800, 0.529, 0.480, 0.478, 0.476, 0.497, 0.518, ...
    0.539, 0.560, 0.581, 0.599, 0.616, 0.627, 0.632, 0.634]; %inches squared
shockLocation = [inf, inf, inf, 15, 14, 6];

%% Experimental

% 1) pretty much anytime ~2-54s is under expanded
% 2) ~55s looks like 3rd crit
% 3) ~1:56 -> 116s
% 4) ~2:26 -> 146s 
% 5) ~2:55 -> 175s 
% 6) ~3:57 -> 237s 
caseTimes = [30, 55, 116, 146, 175, 237];

data = load('PressureData.csv');

case1 = data(caseTimes(1),:)*0.000145038 + P_atm;
case2 = data(caseTimes(2),:)*0.000145038 + P_atm;
case3 = data(caseTimes(3),:)*0.000145038 + P_atm;
case4 = data(caseTimes(4),:)*0.000145038 + P_atm;
case5 = data(caseTimes(5),:)*0.000145038 + P_atm;
case6 = data(caseTimes(6),:)*0.000145038 + P_atm;


%% Case 1
%P=P=P 12e
% M1 = M2 = Me
ptotal1 = case1(throatPos)*((1 + ((gamma - 1)*0.5))^(gamma/(gamma - 1)));
Mach1 = real(sqrt(((((ptotal1./case1).^((gamma - 1)/gamma)) - 1)*2)/(gamma -1)));


%% Case 2

ptotal2 = case2(throatPos)*((1 + ((gamma - 1)*0.5))^(gamma/(gamma - 1)));
Mach2 = real(sqrt(((((ptotal2./case2).^((gamma - 1)/gamma)) - 1)*2)/(gamma -1)));

%% Case 3
ptotal3 = case3(throatPos)*((1 + ((gamma - 1)*0.5))^(gamma/(gamma - 1)));
Mach3 = real(sqrt(((((ptotal3./case3).^((gamma - 1)/gamma)) - 1)*2)/(gamma -1)));


%% Case 4
ptotal4 = case4(throatPos)*((1 + ((gamma - 1)*0.5))^(gamma/(gamma - 1)));
Mach4 = real(sqrt(((((ptotal4./case4).^((gamma - 1)/gamma)) - 1)*2)/(gamma -1)));

Mach4_2 = getM2(Mach4(15), gamma);
ptotal4_2 = case4(15)*((1 + ((gamma - 1)*0.5)*(Mach4_2^2))^(gamma/(gamma - 1)));
Mach4(16) = 0;
Mach4(16)  = real(sqrt(((((ptotal4_2./case4(end)).^((gamma - 1)/gamma)) - 1)*2)/(gamma -1)));


%% Case 5
ptotal5 = case5(throatPos)*((1 + ((gamma - 1)*0.5))^(gamma/(gamma - 1)));
Mach5 = real(sqrt(((((ptotal5./case5).^((gamma - 1)/gamma)) - 1)*2)/(gamma -1)));

Mach5_2 = getM2(Mach5(10), gamma);
ptotal5_2 = case5(10)*((1 + ((gamma - 1)*0.5)*(Mach5_2^2))^(gamma/(gamma - 1)));
Mach5(16) = 0;
Mach5(11:end)  = real(sqrt(((((ptotal5_2./case5(10:end)).^((gamma - 1)/gamma)) - 1)*2)/(gamma -1)));

%% Case 6

ptotal6 = case6(throatPos)*((1 + ((gamma - 1)*0.5))^(gamma/(gamma - 1)));
Mach6 = real(sqrt(((((ptotal6./case6).^((gamma - 1)/gamma)) - 1)*2)/(gamma -1)));

Mach6_2 = getM2(Mach6(6),gamma);
ptotal6_2 = case6(6)*((1 + ((gamma - 1)*0.5)*(Mach6_2^2))^(gamma/(gamma - 1)));
Mach6(16) = 0;
Mach6(7:end) = real(sqrt(((((ptotal6_2./case5(6:end)).^((gamma - 1)/gamma)) - 1)*2)/(gamma -1)));

%% END EXPERIMENTAL CALCULATIONS

%% BEGIN THEORETICAL CALCULATIONS
% Andrew Fung

%% Area Ratios
%Area
a = [0.8,0.529,0.48,0.478,0.476,0.497,0.518,0.539,0.56,0.581,0.599,0.616,0.627,0.632,0.634];
%a/astar
aastar = a/a(5);
%% Theoretical Formulas

%A/A*, points up to shock
AAstar = @(M) sqrt((1/(M^2))*(((2/2.4)*(1+(0.2*(M^2))))^(2.4/0.4)));

%mach after shock
M2 = @(M1) sqrt((1+(0.2*(M1^2)))/(1.4*(M1^2)-0.2));

%Virtual throat cross section area
A2star = @(M2,As) sqrt((M2^2)*(As^2)*(((2/2.4)*(1+(0.2*(M2^2))))^(2.4/0.4))); 

%Ptotal/Pstatic
PtPs = @(M) (1+(0.2*(M^2)))^(1.4/0.4);

%Total pressure 2, in or outside nozzle
Pt2 = @(Me) 14.7*((1+(0.2*(Me^2)))^(1.4/0.4));

%Total pressure 1, shock inside nozzle
Pt1in = @(M2,M1) ((1+1.4*(M2^2))/(1+1.4*(M1^2)))*((1+(0.2*(M1^2))^(1.4/0.4)))*((1+(0.2*(M2^2))^(1.4/0.4)));

%% Theoretical Calculations

for i = 1:15
    Msub(i) = flowisentropic(1.4,aastar(i),'sub');
    Msup(i) = flowisentropic(1.4,aastar(i),'sup');
end

MT = zeros(15,4); %theoretical mach
%3rd crit
MT(1:5,1) = Msub(1:5);
MT(6:15,1) = Msup(6:15);
%2nd crit
MT(1:5,2) = Msub(1:5);
MT(6:15,2) = Msup(6:15);
%In nozzle
MT(1:5,3) = Msub(1:5);
MT(6:14,3) = Msup(6:14);
M2in = M2(Msup(14));
VirAAstarin = a(15)/A2star(M2in,a(15));
MT(15,3) = flowisentropic(1.4,VirAAstarin,'sub');
%1st crit
MT(1:5,4) = Msub(1:5);
MT(6,4) = Msup(6);
M21st = M2(Msup(6));
for i = 7:15
    curVirAAstar1st = a(i)/A2star(M21st,a(6));
    MT(i,4) = flowisentropic(1.4,curVirAAstar1st,'sub');
end

PtT = zeros(15,4);
PsT = zeros(15,4);
%3rd crit
for i = 1:15
    PtT(i,1) = Pt2(MT(15,1));
end
%2nd crit
for i = 1:15
    PtT(i,2) = Pt2(MT(15,2));
end
%In nozzle
for i = 1:14
    PtT(i,3) = 2*(Pt2(MT(15,3)));
end
PtT(15,3) = 2*Pt1in(MT(15,3),MT(14,3));
%1st crit
for i = 1:15
    PtT(i,4) = Pt2(MT(15,4));
end

%all static pressures
for i = 1:4
    for j = 1:15
        curPtPs = PtPs(MT(j,i));
        PsT(j,i) = PtT(j,i)/curPtPs;
    end
end

%% Making total experimental pressures a matrix
PtE = zeros(16,4);
PtE(:,1) = ptotal2;
PtE(:,2) = ptotal4;
for i = 1:14
    PtE(i,3) = ptotal5;
end
PtE(15,3) = ptotal5_2;
for i = 1:6
    PtE(i,4) = ptotal6;
end
for i = 7:15
    PtE(i,4) = ptotal6_2;
end

%% graphing pressure vs inlet position
%
figure(1)
plot(inlPos, case2)
grid on
pt2(1:15) = ptotal2;
hold on
plot(inlPos, pt2)
xlabel('Inlet Positions (inches)');
ylabel('Pressure (psi)');
title('Case 2, 3rd Critical: Pressure vs. Location');
legend('Static Pressure','Total Pressure','Location','southwest')
set(gca,'FontSize',14)
exportgraphics(gcf,'f1.png','Resolution',300)
hold off

figure(2)
inlPos4 = inlPos;
inlPos4(16) = inlPos(15)
plot(inlPos, case4)
grid on
pt4(1:15) = ptotal4;
pt4(16) = ptotal4_2;
hold on
plot(inlPos4, pt4)
xlabel('Inlet Positions (inches)');
ylabel('Pressure (psi)');
title('Case 4, 2nd Critical: Pressure vs. Location');
legend('Static Pressure','Total Pressure','Location','southwest')
set(gca,'FontSize',14)
exportgraphics(gcf,'f2.png','Resolution',300)
hold off

figure(3)
inlPos5(1:9) = inlPos(1:9);
inlPos5(10:11) = inlPos(10);
inlPos5(12:16) = inlPos(11:end)
plot(inlPos, case5)
grid on
pt5(1:10) = ptotal5;
pt5(11:16) = ptotal5_2;
hold on
plot(inlPos5, pt5)
xlabel('Inlet Positions (inches)');
ylabel('Pressure (psi)');
title('Case 5, Normal Shock in Nozzle: Pressure vs. Location');
legend('Static Pressure','Total Pressure','Location','southwest')
set(gca,'FontSize',14)
exportgraphics(gcf,'f3.png','Resolution',300)
hold off

figure(4)
inlPos6(1:5) = inlPos(1:5);
inlPos6(6:7) = inlPos(6);
inlPos6(8:16) = inlPos(7:end)
plot(inlPos, case6)
grid on
pt6(1:6) = ptotal6;
pt6(7:16) = ptotal6_2;
hold on
plot(inlPos6, pt6)
xlabel('Inlet Positions (inches)');
ylabel('Pressure (psi)');
title('Case 6, 1st Critical: Pressure vs. Location');
legend('Static Pressure','Total Pressure','Location','southwest')
set(gca,'FontSize',14)
exportgraphics(gcf,'f4.png','Resolution',300)
hold off

figure(5)
plot(inlPos, PsT(:,1))
grid on
hold on
plot(inlPos, PtT(:,1))
xlabel('Inlet Positions (inches)');
ylabel('Pressure (psi)');
title('Theoretical Case 2, 3rd Critical: Pressure vs. Location');
legend('Static Pressure','Total Pressure','Location','southwest')
set(gca,'FontSize',14)
exportgraphics(gcf,'f5.png','Resolution',300)
hold off

figure(6)
plot(inlPos, PsT(:,2))
grid on
hold on
plot(inlPos, PtT(:,2))
xlabel('Inlet Positions (inches)');
ylabel('Pressure (psi)');
title('Theoretical Case 4, 2nd Critical: Pressure vs. Location');
legend('Static Pressure','Total Pressure','Location','southwest')
set(gca,'FontSize',14)
exportgraphics(gcf,'f6.png','Resolution',300)
hold off

figure(7)
plot(inlPos, PsT(:,3))
grid on
hold on
plot(inlPos, PtT(:,3))
xlabel('Inlet Positions (inches)');
ylabel('Pressure (psi)');
title('Theoretical Case 5, Shock in Nozzle: Pressure vs. Location');
legend('Static Pressure','Total Pressure','Location','southwest')
set(gca,'FontSize',14)
exportgraphics(gcf,'f7.png','Resolution',300)
hold off

figure(8)
plot(inlPos, PsT(:,4))
grid on
hold on
plot(inlPos, PtT(:,4))
xlabel('Inlet Positions (inches)');
ylabel('Pressure (psi)');
title('Theoretical Case 6, 1st Critical: Pressure vs. Location');
legend('Static Pressure','Total Pressure','Location','southwest')
set(gca,'FontSize',14)
exportgraphics(gcf,'f8.png','Resolution',300)
hold off

figure(9)
plot(inlPos, Mach2)
grid on
hold on
plot(inlPos, MT(:,1))
xlabel('Inlet Position (inches)')
ylabel('Mach Number')
title('Case 2, 3rd Critical: Theoretical vs. Experimental Mach Number')
legend('Experimental','Theoretical')
set(gca,'FontSize',14)
exportgraphics(gcf,'f9.png','Resolution',300)
hold off

figure(10)
plot(inlPos, Mach4(2:16))
grid on
hold on
plot(inlPos, MT(:,2))
xlabel('Inlet Position (inches)')
ylabel('Mach Number')
title('Case 4, 2nd Critical: Theoretical vs. Experimental Mach Number')
legend('Experimental','Theoretical')
set(gca,'FontSize',14)
exportgraphics(gcf,'f10.png','Resolution',300)
hold off

figure(11)
plot(inlPos, Mach5(2:16))
grid on
hold on
plot(inlPos, MT(:,3))
xlabel('Inlet Position (inches)')
ylabel('Mach Number')
title('Case 5, Shock in Nozzle: Theoretical vs. Experimental Mach Number')
legend('Experimental','Theoretical')
set(gca,'FontSize',14)
exportgraphics(gcf,'f11.png','Resolution',300)
hold off

figure(12)
plot(inlPos, Mach6(2:16))
grid on
hold on
plot(inlPos, MT(:,4))
xlabel('Inlet Position (inches)')
ylabel('Mach Number')
title('Case 6, 1st Critical: Theoretical vs. Experimental Mach Number')
legend('Experimental','Theoretical')
set(gca,'FontSize',14)
exportgraphics(gcf,'f12.png','Resolution',300)
hold off

figure(13)
plot(inlPos, case2)
grid on
hold on
plot(inlPos, PsT(:,1))
plot(inlPos, PtT(:,1))
plot(inlPos, pt2)
xlabel('Inlet Positions (inches)');
ylabel('Pressure (psi)');
title('Case 2, 3rd Critical: Pressure vs. Location, Theoretical and Measured');
legend('Measured Static Pressure','Measured Total Pressure',...
    'Theoretical Static Pressure','Theoretical Total Pressure', 'Location','southwest')
set(gca,'FontSize',14)
exportgraphics(gcf,'f13.png','Resolution',300)
hold off

figure(14)
plot(inlPos, case4)
grid on
hold on
plot(inlPos, PsT(:,2))
plot(inlPos, PtT(:,2))
plot(inlPos4, pt4)
xlabel('Inlet Positions (inches)');
ylabel('Pressure (psi)');
title('Case 4, 2nd Critical: Pressure vs. Location, Theoretical and Measured');
legend('Measured Static Pressure','Measured Total Pressure',...
    'Theoretical Static Pressure','Theoretical Total Pressure', 'Location','southwest')
set(gca,'FontSize',14)
exportgraphics(gcf,'f14.png','Resolution',300)
hold off

figure(15)
plot(inlPos, case5)
grid on
hold on
plot(inlPos, PsT(:,3))
plot(inlPos, PtT(:,3))
plot(inlPos5, pt5)
xlabel('Inlet Positions (inches)');
ylabel('Pressure (psi)');
title('Case 5, Normal Shock in Nozzle: Pressure vs. Location,  Theoretical and Measured');
legend('Measured Static Pressure','Measured Total Pressure',...
    'Theoretical Static Pressure','Theoretical Total Pressure', 'Location','southwest')
set(gca,'FontSize',14)
exportgraphics(gcf,'f15.png','Resolution',300)
hold off

figure(16)
plot(inlPos, case6)
grid on
hold on
plot(inlPos, PsT(:,4))
plot(inlPos, PtT(:,4))
plot(inlPos6, pt6)
xlabel('Inlet Positions (inches)');
ylabel('Pressure (psi)');
title('Case 6, 1st Critical: Pressure vs. Location,  Theoretical and Measured');
legend('Measured Static Pressure','Measured Total Pressure',...
    'Theoretical Static Pressure','Theoretical Total Pressure', 'Location','southwest')
set(gca,'FontSize',14)
exportgraphics(gcf,'f16.png','Resolution',300)
hold off

%}

