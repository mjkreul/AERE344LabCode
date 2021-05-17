%% AerE 344 Lab 11
% Section 1
% Group 1

clear, clc

%% Constants

P_atm = 99000; %pascals
R = 287; %temp constant J/kg*K
T = 25 + 273.15; %Kelvin
rho = P_atm/(R*T); %density of air kg/m^3
c = 0.101; %chord length, m
Uinf = 10; %Uinfinity, m/s
Rec = 68000; %chord reynolds number

%% Importing and parsing data
% reformatting CSVs from parselab11files

%initializing variables
iAOA4 = zeros(3750,4,100);
iAOA8 = zeros(3750,4,100);
iAOA12 = zeros(3750,4,100);
iAOA16 = zeros(3750,4,100);

%loading outputted files
raw4 = load('AOA4.csv');
raw8 = load('AOA8.csv');
raw12 = load('AOA12.csv');
raw16 = load('AOA16.csv');

kill = 0;
i = 1; %index for position in data
j = 1; %index for 3d vector
while kill == 0
    if i < 401
        iAOA4(:,:,j) = raw4(:,i:i+3);
        iAOA8(:,:,j) = raw8(:,i:i+3);
        iAOA12(:,:,j) = raw12(:,i:i+3);
        iAOA16(:,:,j) = raw16(:,i:i+3);
        
        i = i+4;
        j = j+1;
    else
        kill = 1;
    end
end

%% Average Velocity

%initialize average variables
avgAOA4 = zeros(3750,4);
avgAOA8 = zeros(3750,4);
avgAOA12 = zeros(3750,4);
avgAOA16 = zeros(3750,4);

%adding every value to average
for i = 1:100
    for j = 1:3750
        for k = 1:4
            avgAOA4(j,k) = avgAOA4(j,k)+iAOA4(j,k,i);
            avgAOA8(j,k) = avgAOA8(j,k)+iAOA8(j,k,i);
            avgAOA12(j,k) = avgAOA12(j,k)+iAOA12(j,k,i);
            avgAOA16(j,k) = avgAOA16(j,k)+iAOA16(j,k,i);
        end
    end
end

%dividing entire sums by n
avgAOA4 = avgAOA4/100;
avgAOA8 = avgAOA8/100;
avgAOA12 = avgAOA12/100;
avgAOA16 = avgAOA16/100;

%% Turbulent Velocity Fluctuations

tvfX = zeros(3750,4);
tvfY = zeros(3750,4);
n = 100;

%AOA 4
for i = 1:3750
    sumX = 0;
    sumY = 0;
    for j = 1:100
        sumX = sumX + (iAOA4(i,3,j)-avgAOA4(i,3));
        sumY = sumY + (iAOA4(i,4,j)-avgAOA4(i,4));
    end
    tvfX(i,1) = sqrt((sumX^2)/n);
    tvfY(i,1) = sqrt((sumY^2)/n);
end

%AOA 8
for i = 1:3750
    sumX = 0;
    sumY = 0;
    for j = 1:100
        sumX = sumX + (iAOA8(i,3,j)-avgAOA8(i,3));
        sumY = sumY + (iAOA8(i,4,j)-avgAOA8(i,4));
    end
    tvfX(i,2) = sqrt((sumX^2)/n);
    tvfY(i,2) = sqrt((sumY^2)/n);
end

%AOA 12
for i = 1:3750
    sumX = 0;
    sumY = 0;
    for j = 1:100
        sumX = sumX + (iAOA12(i,3,j)-avgAOA12(i,3));
        sumY = sumY + (iAOA12(i,4,j)-avgAOA12(i,4));
    end
    tvfX(i,3) = sqrt((sumX^2)/n);
    tvfY(i,3) = sqrt((sumY^2)/n);
end

%AOA 16
for i = 1:3750
    sumX = 0;
    sumY = 0;
    for j = 1:100
        sumX = sumX + (iAOA16(i,3,j)-avgAOA16(i,3));
        sumY = sumY + (iAOA16(i,4,j)-avgAOA16(i,4));
    end
    tvfX(i,4) = sqrt((sumX^2)/n);
    tvfY(i,4) = sqrt((sumY^2)/n);
end

%% Turbulent Kinetic Energy Distributions
ufTKE = zeros(3750,4); %unformatted TKE

for i = 1:3750
    for j = 1:4
        ufTKE(i,j) = 0.5*rho*((tvfX(i,j)^2)+(tvfY(i,j)^2));
    end
end


%% Reformatting data
% Placing data in a spatially correct matrix
% place in matrix = coordinates

% creating positions
minx = min(avgAOA4(:,1));
maxx = max(avgAOA4(:,1));
xstep = (maxx-minx)/74;
xrange = (minx:xstep:maxx);

miny = min(avgAOA4(:,2));
maxy = max(avgAOA4(:,2));
ystep = (maxy-miny)/49;
yrange = (miny:ystep:maxy);
yrange = flip(yrange);

[x,y] = meshgrid(xrange,yrange);

% initializing matricies
Vx = zeros(50,75,4);
Vy = zeros(50,75,4);

% reassigning velocity values
for i = 1:3750
    %finding matching indicies
    xind = find(abs(x(1,:)-avgAOA4(i,1))<=0.01);
    yind = find(abs(y(:,1)-avgAOA4(i,2))<=0.01);
    
    %creating velocity plot vectors
    Vx(yind,xind,1) = avgAOA4(i,3);
    Vy(yind,xind,1) = avgAOA4(i,4);
    Vx(yind,xind,2) = avgAOA8(i,3);
    Vy(yind,xind,2) = avgAOA8(i,4);
    Vx(yind,xind,3) = avgAOA12(i,3);
    Vy(yind,xind,3) = avgAOA12(i,4);
    Vx(yind,xind,4) = avgAOA16(i,3);
    Vy(yind,xind,4) = avgAOA16(i,4);
    
    %creating TKE vectors
    TKE(yind,xind,1) = ufTKE(i,1);
    TKE(yind,xind,2) = ufTKE(i,2);
    TKE(yind,xind,3) = ufTKE(i,3);
    TKE(yind,xind,4) = ufTKE(i,4);
end

%% Calculating vorticity
vort = zeros(50,75,4);

for i = 3:48
    for j = 3:73
        for k = 1:4
%             T1 = 0.5*xstep*(Vx(i-1,j-1,k)+2*Vx(i,j-1,k)+Vx(i+1,j-1,k));
%             T2 = 0.5*ystep*(Vy(i+1,j-1,k)+2*Vy(i+1,j,k)+Vy(i+1,j+1,k));
%             T3 = 0.5*xstep*(Vx(i+1,j+1,k)+2*Vx(i,j+1,k)+Vx(i-1,j+1));
%             T4 = 0.5*ystep*(Vy(i-1,j+1,k)+2*Vy(i-1,j,k)+Vy(i-1,j-1,k));
%             Gamma = T1+T2-T3-T4;
            dvdx = (-Vy(i + 2, j, k) + 8*Vy(i + 1, j, k) - ...
                8*Vy(i - 1, j, k) + Vy(i - 2, j, k))/(12 * xstep);
            dudy = (-Vx(i, j + 2, k) + 8*Vx(i, j + 1, k) - ...
                8*Vx(i, j - 1, k) + Vx(i, j - 2, k))/(12 * ystep);
            Gamma = dvdx - dudy;
            vort(i,j,k) = Gamma;
        end
    end
end

vort = vort/(4*xstep*ystep);


%% Calculating wake profile
% 1/2 chord = 50.5mm, which is past the edge of recording, so the furthest
% data points are used instead

WP(:,1) = sqrt(Vx(:,1,1).^2+Vy(:,1,1).^2);
WP(:,2) = sqrt(Vx(:,1,2).^2+Vy(:,1,2).^2);
WP(:,3) = sqrt(Vx(:,1,3).^2+Vy(:,1,3).^2);
WP(:,4) = sqrt(Vx(:,1,4).^2+Vy(:,1,4).^2);

%% Plotting

figure(1)
contourf(x,y,vort(:,:,1))
hold on
quiver(x,y,Vx(:,:,1),Vy(:,:,1), 'm')
hold off
title(['Mean Flow Field, AOA 4',char(176)])
xlabel('X Position (mm)')
ylabel('Y Position (mm)')
exportgraphics(gcf,'f1.png','Resolution',300)

figure(2)
contourf(x,y,vort(:,:,2))
hold on
quiver(x,y,Vx(:,:,2),Vy(:,:,2), 'm')
hold off
title(['Mean Flow Field, AOA 8',char(176)])
xlabel('X Position (mm)')
ylabel('Y Position (mm)')
exportgraphics(gcf,'f2.png','Resolution',300)

figure(3)
contourf(x,y,vort(:,:,3))
hold on
quiver(x,y,Vx(:,:,3),Vy(:,:,3), 'm')
hold off
title(['Mean Flow Field, AOA 12',char(176)])
xlabel('X Position (mm)')
ylabel('Y Position (mm)')
exportgraphics(gcf,'f3.png','Resolution',300)

figure(4)
contourf(x,y,vort(:,:,4))
hold on
quiver(x,y,Vx(:,:,4),Vy(:,:,4), 'm')
hold off
title(['Mean Flow Field, AOA 16',char(176)])
xlabel('X Position (mm)')
ylabel('Y Position (mm)')
exportgraphics(gcf,'f4.png','Resolution',300)

figure(5)
contourf(x,y,TKE(:,:,1))
title('Turbulent kinetic energy distribution')
xlabel('X Position (mm)')
ylabel('Y Position (mm)')
exportgraphics(gcf,'f5.png','Resolution',300)

figure(6)
contourf(x,y,TKE(:,:,2))
title('Turbulent kinetic energy distribution')
xlabel('X Position (mm)')
ylabel('Y Position (mm)')
exportgraphics(gcf,'f6.png','Resolution',300)

figure(7)
contourf(x,y,TKE(:,:,3))
title('Turbulent kinetic energy distribution')
xlabel('X Position (mm)')
ylabel('Y Position (mm)')
exportgraphics(gcf,'f7.png','Resolution',300)

figure(8)
contourf(x,y,TKE(:,:,4))
title('Turbulent kinetic energy distribution')
xlabel('X Position (mm)')
ylabel('Y Position (mm)')
exportgraphics(gcf,'f8.png','Resolution',300)

figure(9)
plot(WP(:,1), y)
title(['Wake Profile X/C = 0.5C, AOA 4',char(176)])
xlabel('U/U_{inf}')
ylabel('Y Position (mm)')
exportgraphics(gcf,'f9.png','Resolution',300)

figure(10)
plot(WP(:,2), y)
title(['Wake Profile X/C = 0.5C, AOA 8',char(176)])
xlabel('U/U_{inf}')
ylabel('Y Position (mm)')
exportgraphics(gcf,'f10.png','Resolution',300)

figure(11)
plot(WP(:,3), y)
title(['Wake Profile X/C = 0.5C, AOA 12',char(176)])
xlabel('U/U_{inf}')
ylabel('Y Position (mm)')
exportgraphics(gcf,'f11.png','Resolution',300)

figure(12)
plot(WP(:,4), y)
title(['Wake Profile X/C = 0.5C, AOA 16',char(176)])
xlabel('U/U_{inf}')
ylabel('Y Position (mm)')
exportgraphics(gcf,'f12.png','Resolution',300)







