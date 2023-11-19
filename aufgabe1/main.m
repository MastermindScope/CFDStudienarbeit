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

%this defines mesh sized to be investigated
meshSize = 5:5:1000;


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
    x = linspace(0,L,meshSize(i));
    
    %system is solved and T obtained
    T = A\b;
    plot(x,T);
    hold on
    
    %calculate analytical solution for given mesh size
    anaSol = gam.*exp(mprime.*x)-gam.*a./bprime.*exp(-mprime*x)+T_ar;
    
    %fill error vectors
    errAbs(i) = abs(anaSol(end)-T(end));
    errRel(i) = norm(anaSol(end)-T(end))/norm(anaSol(end));

end

%plot errors and the analytical solution
plot(x,anaSol)
xlabel("x [m]");
ylabel("T [°C]");
figure
semilogx(meshSize, errAbs)
title("Absolute error on the last point")
xlabel("Mesh elements");
ylabel("Absolute error");
grid
figure
semilogx(meshSize, errRel)
title("Relative error on the last point")
xlabel("Mesh elements");
ylabel("Relative error");
grid

hold off