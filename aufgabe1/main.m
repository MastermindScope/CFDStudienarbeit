%% Assignment 1 main.m
%Andreas Wenger & Vinzenz Goetz
%Units = [kg,m,°C]
%
%functions:
%[output1, output2, ...] = fuction(input1, input2, ...)
%matrix = createA(meshSize, k, h, P, L, At, T_w, T_ar)
%vector = createFields(meshSize)
%vector = initFields(vector, Phi0)
%vector = sourceTerms(meshSize, L, P, k, h, At, T_ar, T_w)
%[boolean, scalar] = residual(eps, vectorOld, vectorNew, norm)
%graph = visualize(vectorx,vectorT)
close all
clear all
clc

%Linestyles
lStyle = ["-g", "--b", ":m", "-r.", "-.k", "--m", ":r"];

%this defines mesh sized to be investigated
meshSize = 5:3:20;


%this section defines the constants used in the problem
k = 100;
h = 10;
P = 0.4;
L = 1;
A_t = 0.01;
T_w = 300;
T_ar = 20;
T0 = 20;
%m is m^2 in this case, shorthand because a longer variable is hard to
%write
m = h*P/(k*A_t);

%this is defining some constants which help in the analytical solution of
%the problem
mprime=sqrt(m);
a = h*exp(mprime*L)+mprime*k*exp(mprime*L);
bprime = h*exp(-mprime*L)-mprime*k*exp(-mprime*L);
gam = (T_w-T_ar)/((1-a/bprime));


%this is for no convection on the face
%gam = (T_w-T_ar)/(1+exp(-2*sqrt(m)*L));

%this initializes the error vectors with the right length
errAbs = zeros(1,length(meshSize));
errRel = zeros(1,length(meshSize));

for i = 1:length(meshSize)
    %here, the system matrix as well as source terms are initialized
    dx = L/meshSize(i);
    %system matrix is created
    A = createA(meshSize(i), k, h, P, L, A_t, T_w, T_ar);
    %source terms are created
    b = sourceTerms(meshSize(i), L, P, k, h, A_t, T_ar, T_w);
    x = dx/2:dx:L-dx/2;
    
    %system is solved and T obtained
    T = A\b;
    T_p = T(end);
    T_f = ((2*k/dx)*T_p+h*T_ar)/((2*k/dx)+h);
    T(2:end+1) = T;
    T(1) = T_w;
    T(end+1) = T_f;
    x(2:end+1) = x;
    x(1) = 0;
    x(end+1) = L;
    
    figure(1)
    plot(x,T,lStyle(i),LineWidth=1.5);
    hold on

    
  
    
    %calculate analytical solution for given mesh size
    anaSol = gam.*exp(mprime.*x)-gam.*a./bprime.*exp(-mprime*x)+T_ar;
    
    %fill error vector
    errRel(i) = norm(anaSol-T)/norm(anaSol)/length(T);

    errLen = abs(anaSol-T')./anaSol;
    figure(2)
    plot(x, errLen,lStyle(i),LineWidth=1.5)
    hold on
    
    ltext(i) = strcat("Mesh elements = ", num2str(meshSize(i)));

end

%output of temperature at the end of the fin
fprintf('The temperature at the tip of the fin = %.3f °C for a meshsize of %d elements \n', T_f, meshSize(end))

%plot errors and the analytical solution
figure(2)
title("Relative error over length")
xlabel("Length [m]");
ylabel("Relative error");
legend(ltext);
grid
fontsize(13,"points")
saveas(2, "errLen.png")



figure(1)
plot(x,anaSol,LineWidth=1.5)
xlabel("x [m]");
ylabel("T [°C]");
ltext(end+1) = "analytical solution";
legend(ltext);
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
