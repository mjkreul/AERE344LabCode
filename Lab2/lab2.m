%% AER E 344 Lab 2 Calculations
% Section 1 
% Group 1 
clear, clc
pause on
%% Constants
P_atm = 99000; %pascals
R = 287; %temp constant J/kg*K
rho_water = 997; %density of water kg/m^3

%% Beginning operations

% load the data
load lab2data.csv; 

% takes the tempuratures from the data file and puts them into vector T 
% while converting it to kelvin 
T = lab2data(:, end) + 273.15;

% Get the motor frequency
mot_freq = lab2data(:,1); % Hz

% get all of the height data and convert it to meters
h_atm = lab2data(:,2) * 0.0254; %meters
h_A = lab2data(:,3) * 0.0254; %meters
h_E = lab2data(:,4) * 0.0254; %meters
h_total = lab2data(:,5) * 0.0254; %meters
h_static = lab2data(:,6) * 0.0254; %meters


%% Writing conversion to .dat file

metarr = [mot_freq, h_atm, h_A, h_E, h_total, h_static, T];
writematrix(metarr, "lab2datamet.csv")

%% Initialize arrays
% temperature is constant so we can do this but make it so it's more
% modular later
rho_air = P_atm/(R * T(1)); % kg/m^3 

% Every vector is the same size so this will be the size we can use to
% calculate out delta p in a loop.
[row, col] = size(h_A); 

% initialize vectors 
dp1 = zeros(row, col);
dp2 = zeros(row, col);
dp3 = zeros(row, col);
dp4 = zeros(row, col);

% delta p
pa_pe = zeros(row, col);

% dynamic pressure q_T
pt_ps = zeros(row, col);

% velocity 
v = zeros(row, col);

%% Calculations 
for i = 1:row
    % finding all dp for every pressure gauge 
    dp1(i) = deltaP(rho_water, h_A(i), h_A(1), h_atm(i), h_atm(1));
    dp2(i) = deltaP(rho_water, h_E(i), h_E(1), h_atm(i), h_atm(1));
    dp3(i) = deltaP(rho_water, h_total(i), h_total(1), h_atm(i), h_atm(1));
    dp4(i) = deltaP(rho_water, h_static(i), h_static(1), h_atm(i), h_atm(1));
    
    % finding pa - pe and q_t
    pa_pe(i) = dp1(i) - dp2(i);
    pt_ps(i) = dp3(i) - dp4(i);
    
    % getting the velocity at each point
    v(i) = sqrt((2 * pt_ps(i))/ rho_air); %m/s
end


%% Plotting
% For the P_A - P_E vs. q_T graph

% getting the best fit line for it
coefficients = polyfit(pa_pe, pt_ps, 1);
xFit = linspace(min(pa_pe), max(pa_pe), 100);
yFit = polyval(coefficients , xFit);

%plot
figure(1)
scatter(pa_pe, pt_ps, 'ro')
hold on
grid on
plot(xFit, yFit, 'b-')
title("P_A - P_E vs. q_T")
xlabel("Change in Pressure P_A - P_E (Pa)")
ylabel("Dynamic Pressure q_T (Pa)")
legend("Experimental Data", "Linear Fit Curve")
hold off 

% print out the best fit line for this
fprintf("The best fit line for P_A - P_E vs. q_T is:\n\t %fx + %f\n", coefficients(1), coefficients(2));

pause(1)
% For the V vs motor hz graph

% getting the best fit line for it
coefficients = polyfit(mot_freq, v, 1);
xFit = linspace(min(mot_freq), max(mot_freq), 100);
yFit = polyval(coefficients , xFit);

%plot
figure(2)
scatter(mot_freq, v, 'ro')
hold on
grid on
plot(xFit, yFit, 'b-')
title("Motor Frequency vs. Flow Velocity")
xlabel("Motor Frequency (Hz)")
ylabel("Flow Velocity (m/s)")
legend("Experimental Data", "Linear Fit Curve")
hold off

% print out the best fit line for this
fprintf("The best fit line for Motor Frequency vs. Flow Velocity is:\n\t %fx + %f\n", coefficients(1), coefficients(2));
