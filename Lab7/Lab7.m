%% AER E 344 Lab 7 Calculations
% Section 1
% Group 1
clear, clc, clf
pause on

%% Constants

c = 0.101; %chord length (m)
xdist = [2.5,2.7,2.9,3.1,3.3,3.5,3.7,3.9,4.1,4.3,4.5,4.7,4.9,5.1,5.3,5.5]'; %x distances of sensor
y = xdist - 4;
aoavals = [0,4,8,12]; %angle of attack values
C = [0.2884,-0.4925,3.4814,-13.672,13.192]; %calibration coefficients
Uinf = 25.83; %free stream flow velocity at 20 Hz
t = [0:5/9999:5];

%% Opening files parsed by parseLab7Files.m

aoa0 = load('aoa0.csv');
aoa4 = load('aoa4.csv');
aoa8 = load('aoa8.csv');
aoa12 = load('aoa12.csv');

%% Averaging data and getting velocity
avgvolt = zeros(16,4);
v = zeros(16,4);

for i = 1:16 %averaging values
    avgvolt(i,1) = mean(aoa0(:,i));
    avgvolt(i,2) = mean(aoa4(:,i));
    avgvolt(i,3) = mean(aoa8(:,i));
    avgvolt(i,4) = mean(aoa12(:,i));
end

for i = 1:16
    for j = 1:4
        v(i,j) = hotwireconv(avgvolt(i,j));
    end
end

%% Calculating turbulence intensity
TI = zeros(16,4);

for i = 1:16 %Angle of Attack 0
    velsum = 0;
    curVel = v(i,1);
    
    for j = 1:10000 %each point of data
        itVel = hotwireconv(aoa0(j,i));
        velsum = velsum + ((itVel-curVel)^2);
    end
    
    curTI = (sqrt(velsum/10000))/curVel;
    TI(i,1) = curTI;
end

for i = 1:16 %Angle of Attack 4
    velsum = 0;
    curVel = v(i,2);
    
    for j = 1:10000 %each point of data
        itVel = hotwireconv(aoa4(j,i));
        velsum = velsum + ((itVel-curVel)^2);
    end
    curTI = (sqrt(velsum/10000))/curVel;
    TI(i,2) = curTI;
end

for i = 1:16 %Angle of Attack 8
    velsum = 0;
    curVel = v(i,3);
    
    for j = 1:10000 %each point of data
        itVel = hotwireconv(aoa8(j,i));
        velsum = velsum + ((itVel-curVel)^2);
    end
    
    curTI = (sqrt(velsum/10000))/curVel;
    TI(i,3) = curTI;
end

for i = 1:16 %Angle of Attack 12
    velsum = 0;
    curVel = v(i,4);
    
    for j = 1:10000 %each point of data
        itVel = hotwireconv(aoa12(j,i));
        velsum = velsum + ((itVel-curVel)^2);
    end
    
    curTI = (sqrt(velsum/10000))/curVel;
    TI(i,4) = curTI;
end

%% Calculating y/delta & U/Ue

uue = v/Uinf;


%% Drag

for i = 1:4 %calculating drag coefficients
    %coefficient of drag
    dragSum = 0;
    trapTerm = 0;
    
    %trapezoidal rule for getting the coefficient of drag
    for j = 2:15
        trapTerm = trapTerm + ((v(j,i)/Uinf)*(1 - v(j,i)/Uinf))*(2/c);
    end
    
    trapTerm = (0.00508)*(trapTerm + (((v(16,i)/Uinf)*(1 - v(16,i)/Uinf)) + ((v(1,i)/Uinf)*(1 - v(1,i)/Uinf)))/2);
    Cd(i) = trapTerm;
end

%% FFT Transform

%getting all velocities
vaoa0 = zeros(10000,16);
vaoa4 = zeros(10000,16);
vaoa8 = zeros(10000,16);
vaoa12 = zeros(10000,16);
for i = 1:10000
    for j = 1:16
        vaoa0(i,j) = hotwireconv(aoa0(i,j));
        vaoa4(i,j) = hotwireconv(aoa4(i,j));
        vaoa8(i,j) = hotwireconv(aoa8(i,j));
        vaoa12(i,j) = hotwireconv(aoa12(i,j));
    end
end

%aoa0
j = 1; %outwake interation
k = 1; %inwake iteration
for i = 1:16
    curVel = v(i,1);
    if curVel > 24.77
        for m = 1:10000
            outWake0(m,j) = vaoa0(m,j);
        end
        j = j+1;
    else
        for m = 1:10000
            inWake0(m,k) = vaoa0(m,k);
        end
        k = k+1;
    end
end

%aoa4
j = 1; %outwake interation
k = 1; %inwake iteration
for i = 1:16
    curVel = v(i,2);
    if curVel > 23.9
        for m = 1:10000
            outWake4(m,j) = vaoa0(m,j);
        end
        j = j+1;
    else
        for m = 1:10000
            inWake4(m,k) = vaoa0(m,k);
        end
        k = k+1;
    end
end

%aoa8
j = 1; %outwake interation
k = 1; %inwake iteration
for i = 1:16
    curVel = v(i,3);
    if curVel > 24.55
        for m = 1:10000
            outWake8(m,j) = vaoa0(m,j);
        end
        j = j+1;
    else
        for m = 1:10000
            inWake8(m,k) = vaoa0(m,k);
        end
        k = k+1;
    end
end

%aoa12
j = 1; %outwake interation
k = 1; %inwake iteration
for i = 1:16
    curVel = v(i,4);
    if curVel > 24.77
        for m = 1:10000
            outWake12(m,j) = vaoa0(m,j);
        end
        j = j+1;
    else
        for m = 1:10000
            inWake12(m,k) = vaoa0(m,k);
        end
        k = k+1;
    end
end

Fs = 2000;
f = Fs*(0:(10000/2))/10000;

%fft of angle of attack 0
Yout0 = fft(outWake0-mean(outWake0));
P2out0 = 2*abs(Yout0/10000);
P1out0 = P2out0(1:10000/2+1);

Yin0 = fft(inWake0-mean(inWake0));
P2in0 = 2*abs(Yin0/10000);
P1in0 = P2in0(1:10000/2+1);

%fft of angle of attack 4
Yout4 = fft(outWake4-mean(outWake4));
P2out4 = 2*abs(Yout4/10000);
P1out4 = P2out4(1:10000/2+1);

Yin4 = fft(inWake4-mean(inWake4));
P2in4 = 2*abs(Yin4/10000);
P1in4 = P2in4(1:10000/2+1);

%fft of angle of attack 8
Yout8 = fft(outWake8-mean(outWake8));
P2out8 = 2*abs(Yout8/10000);
P1out8 = P2out8(1:10000/2+1);

Yin8 = fft(inWake8-mean(inWake8));
P2in8 = 2*abs(Yin8/10000);
P1in8 = P2in8(1:10000/2+1);

%fft of angle of attack 12
Yout12 = fft(outWake12-mean(outWake12));
P2out12 = 2*abs(Yout12/10000);
P1out12 = P2out12(1:10000/2+1);

Yin12 = fft(inWake12-mean(inWake12));
P2in12 = 2*abs(Yin12/10000);
P1in12 = P2in12(1:10000/2+1);


%% Plotting

runPlots = 1; %for turning off plotting if debugging something unrelated
if runPlots == 1
    figure(1)
    plot(aoavals,Cd,'-*')
    title('Coefficient of Drag vs Angle of Attack')
    xlabel('Angle of Attack')
    ylabel('Coefficient of Drag (C_d)')
    
    figure(2) %AOA 0
    subplot(3,2,1)
    plot(y,v(:,1),'-*')
    title('Mean Velocity Profile')
    xlabel('Y/\theta')
    ylabel('U_m_e_a_n(m/s)')
    subplot(3,2,2)
    plot(y,TI(:,1),'-*')
    title('Turbulence Intensity')
    xlabel('Y/\theta')
    ylabel('TI')
    subplot(3,2,3)
    plot(t,vaoa0(:,1))
    title('U Velocity Plot Outside Wake')
    xlabel('Time (s)')
    ylabel('U (m/s)')
    subplot(3,2,4)
    plot(t,vaoa0(:,10))
    title('U Velocity Plot Inside Wake')
    xlabel('Time (s)')
    ylabel('U (m/s)')
    subplot(3,2,5)
    plot(f,P1out0)
    title('Single-Sided Spectrum Outside Wake')
    xlabel('Frequency (Hz)')
    subplot(3,2,6)
    plot(f,P1in0)
    title('Single-Sided Spectrum Inside Wake')
    xlabel('Frequency (Hz)')
    sgtitle(['0' char(176) ' Angle of Attack'])
    
    figure(3) %AOA 4
    subplot(3,2,1)
    plot(y,v(:,2),'-*')
    title('Mean Velocity Profile')
    xlabel('Y/\theta')
    ylabel('U_m_e_a_n(m/s)')
    subplot(3,2,2)
    plot(y,TI(:,2),'-*')
    title('Turbulence Intensity')
    xlabel('Y/\theta')
    ylabel('TI')
    subplot(3,2,3)
    plot(t,vaoa4(:,1))
    title('U Velocity Plot Outside Wake')
    xlabel('Time (s)')
    ylabel('U (m/s)')
    subplot(3,2,4)
    plot(t,vaoa4(:,8))
    title('U Velocity Plot Inside Wake')
    xlabel('Time (s)')
    ylabel('U (m/s)')
    subplot(3,2,5)
    plot(f,P1out4)
    title('Single-Sided Spectrum Outside Wake')
    xlabel('Frequency (Hz)')
    subplot(3,2,6)
    plot(f,P1in4)
    title('Single-Sided Spectrum Inside Wake')
    xlabel('Frequency (Hz)')
    sgtitle(['4' char(176) ' Angle of Attack'])
    
    figure(4) %AOA 8
    subplot(3,2,1)
    plot(y,v(:,3),'-*')
    title('Mean Velocity Profile')
    xlabel('Y/\theta')
    ylabel('U_m_e_a_n(m/s)')
    subplot(3,2,2)
    plot(y,TI(:,3),'-*')
    title('Turbulence Intensity')
    xlabel('Y/\theta')
    ylabel('TI')
    subplot(3,2,3)
    plot(t,vaoa8(:,1))
    title('U Velocity Plot Outside Wake')
    xlabel('Time (s)')
    ylabel('U (m/s)')
    subplot(3,2,4)
    plot(t,vaoa8(:,8))
    title('U Velocity Plot Inside Wake')
    xlabel('Time (s)')
    ylabel('U (m/s)')
    subplot(3,2,5)
    plot(f,P1out8)
    title('Single-Sided Spectrum Outside Wake')
    xlabel('Frequency (Hz)')
    subplot(3,2,6)
    plot(f,P1in8)
    title('Single-Sided Spectrum Inside Wake')
    xlabel('Frequency (Hz)')
    sgtitle(['8' char(176) ' Angle of Attack'])
    
    figure(5) %AOA 12
    subplot(3,2,1)
    plot(y,v(:,4),'-*')
    title('Mean Velocity Profile')
    xlabel('Y/\theta')
    ylabel('U_m_e_a_n(m/s)')
    subplot(3,2,2)
    plot(y,TI(:,4),'-*')
    title('Turbulence Intensity')
    xlabel('Y/\theta')
    ylabel('TI')
    subplot(3,2,3)
    plot(t,vaoa12(:,1))
    title('U Velocity Plot Outside Wake')
    xlabel('Time (s)')
    ylabel('U (m/s)')
    subplot(3,2,4)
    plot(t,vaoa12(:,8))
    title('U Velocity Plot Inside Wake')
    xlabel('Time (s)')
    ylabel('U (m/s)')
    subplot(3,2,5)
    plot(f,P1out12)
    title('Single-Sided Spectrum Outside Wake')
    xlabel('Frequency (Hz)')
    subplot(3,2,6)
    plot(f,P1in12)
    title('Single-Sided Spectrum Inside Wake')
    xlabel('Frequency (Hz)')
    sgtitle(['12' char(176) ' Angle of Attack'])
end

%% Functions

function v = hotwireconv(volt)
%conversion from voltage to velocity using the given equation
v = (0.2889*(volt^4))-(0.4295*(volt^3))+(3.4814*(volt^2))-(13.672*volt)+13.192;

end