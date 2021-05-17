%% AerE 344 Lab 8 Calculations
% Section 1
% Group 1
% Author: Andrew Fung

clear, clc, clf

%% Constants
P_atm = 99000; %pascals
R = 287; %temp constant J/kg*K
rho_water = 997; %density of water kg/m^3
T = 25 + 273.15; %Kelvin
rho = P_atm/(R*T); %density of air kg/m^3
K = 1.14393; % may need to change this
mu = 1.825*10^(-5); % dynamic viscosity of air
c = 0.101; %chord length (m)
inH2OtoPa = 248.84; %Pa

%% Importing and Parsing Data

%defining empty variables
distance(1:11) = [0:1:10];
distance(12:22) = [15:5:65];
y = [4:4:120];
pitot = zeros(2,22);
data = zeros(30,22);

for i = 1:22
    %loading files
    curDistance = distance(i);
    if ispc == 1
        fileName = strcat('Data\',num2str(curDistance),'in.csv');
    else
        % changed to add this since I have a Mac and the way files are
        % organized are a little different
        % https://www.mathworks.com/help/matlab/ref/ispc.html
        fileName = strcat('Data/',num2str(curDistance),'in.csv');
    end
    curDataSet = load(fileName);
    %calculating mean and eliminating unneeded values
    %Don't remember why I did it this exact way but it works ig
%     curMean = curDataSet(1,:);
    % Get the pitot tube values
    pitot(:,i) = mean(curDataSet(:,32:33))';
%     curMean(1) = [];
%     curMean(17:33) = [];
%     curMean(31:48) = [];
    curMean = mean(curDataSet(:,2:31));
    curMean(12) = curMean(11);
    for j = 17:22 %This just fills a funky area according to Ram's suggestion
        curMean(j) = curMean(16);
    end
    data(:,i) = curMean';
    % place the rest of the values in data 
%     data(:,i) = mean(curDataSet(:,2:31))';
end





% ----------------------------------------------------------------
% -------------CHANGE THIS WHEN IT WORKS---------------
% ----------------------------------------------------------------

% Changing the fucked values
% These are just baseless guesses so I can make the rest of the code 
% without having imaginary velocities

% ----------------------------------------------------------------
% ----------------------------------------------------------------
% ----------------------------------------------------------------
% data(1:3,2) = [-50;-24;-10];
% data(1:4,3) = [-70;-60;-50;-40];
% data(1:4,4) = [-70;-60;-50;-40];
% -----------------------------------------------------------------
% ----------------------------------------------------------------
% ----------------------------------------------------------------







%% Calculating Velocity
Uinf = zeros(22,1);
U = zeros(30,22);
for i = 1:22
    dP = pitot(2,i) - pitot(1,i); %Ptot-Pstat
    curUinf = sqrt(2*dP/rho);
    Uinf(i,1) = curUinf; %Saves UInfinity
    for j = 1:30 %Gets velocities for all points
        U(j,i) = sqrt(2*( abs(data(j,i)-pitot(1,i)))/rho);
    end
end
pppoopoo = Uinf
Uinf = mean(Uinf);
blCutoff = 0.99*Uinf;

%% Plot 1 data

BLT = zeros(22,1);
BLTindex = zeros(22,1);
BLT(1:2) = 0;
BLTindex(1:2) = 1;

%finding BL thickness
for i = 3:22
    blFound = 0;
    j = 1; %index searches through velocities
    while blFound == 0 && j <= 30
        curU = U(j,i);
        if curU >= mean(blCutoff) %|| j == 29 %Change this to change boundary layer definition
            blFound = 1;
            BLT(i) = y(j);
            BLTindex(i) = j;
        end
        j = j+1;
    end
end

%y/delta
yBLT = zeros(30,22);
%U/Uinf
uuinf = zeros(30,22);
for i = 1:22
    delta = BLT(i); %BL Thickness at x distance
    for j = 1:30
        yBLT(j,i) = y(j)/delta; %constant y positions of rake/current BL thickness
%         uuinf(j,i) = U(j,i)/Uinf(i);
        uuinf(j,i) = U(j,i)/Uinf;
    end
end

%% 

momthicc = zeros(1,22); %momentum thickness
Cd = zeros(1,22);
Cf = zeros(1,22);
for i = 1:22 %calculating drag coefficients
    trapTerm = 0;
    CurBLT = BLTindex(i);
%     j = 1;
    %trapezoidal rule for getting the coefficient of drag
    for j = 2:CurBLT-1
%         trapTerm = trapTerm + ((U(j,i)/Uinf(i))*(1 - U(j,i)/Uinf(i)));
%         j = j+1;
        trapTerm = trapTerm + ((U(j,i)/Uinf)*(1 - U(j,i)/Uinf));
    end
%     trapTerm = 0.004*(trapTerm + (((U(CurBLT,i)/Uinf(i))*(1 - U(CurBLT,i)/Uinf(i)))...
%         + ((U(1,i)/Uinf(i))*(1 - U(1,i)/Uinf(i))))/2);
    trapTerm = 0.004*(trapTerm + (((U(CurBLT,i)/Uinf)*(1 - U(CurBLT,i)/Uinf))...
        + ((U(1,i)/Uinf)*(1 - U(1,i)/Uinf)))/2);
    momthicc(i) = trapTerm;
    Cd(i) = momthicc(i)*2/distance(i);
end
for i = 1:22
    dtdx = dthetadx(momthicc,distance);
    Cf(i) = 2*dtdx(i);
end

%% Theoretical Values

% theoretical thicknesses
Re = 10^5;
Tlam = 25.4*distance(1:11)*(5/sqrt(Re)); %Thickness Laminar
Ttur(1) = Tlam(end);
Ttur(2:12) = 25.4*distance(12:22)*(0.37/((Re)^0.2));%Thickness Turbulent

TRe = 1.651*rho*mean(Uinf)/mu; %Theoretical Reynolds
TCd = 0.074/(TRe^0.2); %Theoretical Cd




%% Plotting

%y/delta vs u/uinf
figure(1)
plot(uuinf(:,1),yBLT(:,1))
hold on
for i = 2:22 %Lab instructions only ask for 10, but keeping this as all for now
    plot(uuinf(:,i),yBLT(:,i))
end
hold off
xlabel('U/U_\infty')
ylabel('y/\delta')

%boundary layer thicknesses
figure(2)
plot(distance,BLT,'-o')
hold on
plot(distance(1:11),Tlam)
plot(distance(11:22),Ttur)
hold off
title('Boundary Layer Thickness')
xlabel('Distance From Leading Edge (in)')
ylabel('Bounday Layer Thickness (mm)')
legend('Experimental Thickness','Theoretical Laminar Thickness','Theoretial Turbulent Thickness')

%
figure(3)
plot(distance,Cd,'-o')
title('Coefficient of Drag vs X')
xlabel('Distance From Leading Edge (in)')
ylabel('Coefficient of Drag (C_d)')

figure(4)
plot(distance,Cf,'-o')
title('Local Shear Stress Coefficient vs X')
xlabel('Distance From Leading Edge (in)')
ylabel('Local Shear Stress Coefficient (C_f)')


%% Functions

function ddx = dthetadx(theta,x)
%Finitie difference using different nonlinear matrix for x values

if length(theta) == length(x)
    ddx = zeros(1,length(x));
    ndif = 0; %numerator difference
    ddif = 0; %denominator difference
    for i = 1:length(x)
        if i == 1
            ndif = theta(i+1) - theta(i);
            ddif = x(i+1) - x(i);
            ddx(i) = ndif/ddif;
        elseif i ~= length(x)
            ndif = theta(i+1) - theta(i-1);
            ddif = x(i+1) - x(i-1);
            ddx(i) = ndif/ddif;
        else
            ndif = theta(i) - theta(i-1);
            ddif = x(i) - x(i-1);
            ddx(i) = ndif/ddif;
        end
    end
else
    ddx = 0;
    disp('Momentum Thickness, X Values, dimension mismatch')
end
end

