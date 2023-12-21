%% Assignment 2 main.m
%Andreas Wenger & Vinzenz Goetz
%Units = [kg,m,°C,s]

close all
clear all
clc

%Linestyles
lStyle = ["-g", "--b", ":m", "-r.", "-.k", "--m", ":r"];

%this defines mesh sized to be investigated
meshSize = 10:3:25;


%this section defines the constants used in the problem
k = 100;
h = 10;
P = 0.4;
L = 1;
cp = 890;
rho = 2770;
T_p0 = 20;
dt = 0.01;
timesteps = 1000000;
A_t = 0.01;
T_w = 300;
T_ar = 20;
%m is m^2 in this case, shorthand because a longer variable is hard to
%write
m = h*P/(k*A_t);
n = h*P/(cp*A_t);

%this is defining some constants which help in the analytical solution of
%the problem
mprime=sqrt(m);
a = h*exp(mprime*L)+mprime*k*exp(mprime*L);
bprime = h*exp(-mprime*L)-mprime*k*exp(-mprime*L);
gam = (T_w-T_ar)/((1-a/bprime));


%this is for no convection on the face
%gam = (T_w-T_ar)/(1+exp(-2*sqrt(m)*L));

%this initializes the error vectors with the right length
errRel = zeros(1,length(meshSize));

for i = 1:length(meshSize)
    T0 = ones(meshSize(i),1)*20;
    %here, the system matrix as well as source terms are initialized
    dx = L/meshSize(i);
    % calculate element mass
    M_p0 = rho*dx*A_t;
    %system matrix is created
    A = createA(meshSize(i), cp, h, P, L, A_t, dt, k, M_p0);
    %source terms are created
    b = sourceTerms(meshSize(i), L, P, k, h, A_t, T_ar, T_w, cp, M_p0, T0, dt);

    %vector for x coordinate is initialized
    x = dx/2:dx:L-dx/2;

    %initialize a residual vector over time and an array for T in each time
    %step
    res = ones(1,10000);
    T_last = T_p0;
    T_time = (zeros(timesteps,meshSize(i)));


    for j = 1:timesteps
        T = A\b;
        T_time(j,:) = T;
        b = sourceTerms(meshSize(i), L, P, k, h, A_t, T_ar, T_w, cp, M_p0, T, dt);
        res(j) = T_last-T(end);
        T_last = T(end);
    end

    %initialize x and y coordinate arrays for surface plotting
    o = ones(j,1).*x;
    p = ones(1,meshSize(i)).*(1:1:timesteps)';

    % plot residual
    figure(4)
    grid
    semilogy(abs(res))
    hold on

    % shift the Temperature and x vectors around to reflect first and last
    % point on the face of the fin
    T_p = T(end);
    T_f = ((2*k/dx)*T_p+h*T_ar)/((2*k/dx)+h);
    T(2:end+1) = T;
    T(1) = T_w;
    T(end+1) = T_f;
    x(2:end+1) = x;
    x(1) = 0;
    x(end+1) = L;
    
    %plot T over x
    figure(1)
    plot(x,T,lStyle(i),LineWidth=1.5);
    hold on

    
    %calculate analytical solution for given mesh size
    anaSol = gam.*exp(mprime.*x)-gam.*a./bprime.*exp(-mprime*x)+T_ar;
    
    %fill error vector
    errRel(i) = abs(anaSol(end)-T(end))/abs(anaSol(end));

    errLen = abs(anaSol-T')./anaSol;
    figure(2)
    plot(x, errLen,lStyle(i),LineWidth=1.5)
    hold on
    
    ltext(i) = strcat("Mesh elements = ", num2str(meshSize(i)));

end

% plot surface of T,t and x
figure(5)
surf(o,p,T_time);
ylim([0 200]);
ylabel("Time [s]");
xlabel("x [m]");
zlabel("T [°C]");

%output of temperature at the end of the fin
fprintf('The temperature at the tip of the fin T_f = %.3f °C for a meshsize of %d elements \n', T_f, meshSize(end))

%plot errors and the analytical solution
figure(2)
title("Relative error over length")
xlabel("Length [m]");
ylabel("Relative error");
legend(ltext,Location="northwest");
grid
fontsize(13,"points")
saveas(2, "errLen.png")


figure(1)
grid
plot(x,anaSol,LineWidth=1.5)
xlabel("x [m]");
ylabel("T [°C]");
ltext(end+1) = "analytical solution";
legend(ltext);
title("Temperature over length");
%the fontsize command requires version R2022a or later
fontsize(13,"points")
saveas(1, "T.png")

figure(3)
loglog(meshSize, errRel)
title("Relative error on the last point")
xlabel("Mesh elements");
ylabel("Relative error");
grid
fontsize(13,"points")
saveas(3, "errElems.png")

hold off