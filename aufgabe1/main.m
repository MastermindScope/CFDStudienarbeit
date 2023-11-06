%% Assignment 1 main.m
%Andreas Wenger & Vinzenz Goetz
%Units = [kg,m,°C]
%
%A*T+b = T', do this very often
%functions:
%[output1, output2, ...] = fuction(input1, input2, ...)
%matrix = createA(meshSize, k, h, P, L, At, T_w, T_ar)
%vector = createFields(meshSize)
%vector = initFields(vector, Phi0)
%vector = sourceTerms(meshSize, L, P, k, h, At, T_ar, T_w)
%[boolean, scalar] = residual(eps, vectorOld, vectorNew, norm)
%graph = visualize(vectorx,vectorT)

meshSize = input("Input Mesh Size >>>  ");
eps = input("Input Max Residual >>>   ");
if meshSize == 0
    meshSize = 20;
end

if eps == 0
    eps = 1e-6;
end

k = 100;
h = 10;
P = 0.4;
L = 1;
A_t = 0.01;
T_w = 100;
T_ar = 10;
T0 = 3;
dx = L/meshSize;
m = h*P/(k*A_t);


A = createA(meshSize, k, h, P, L, A_t, T_w, T_ar);
T = createFields(meshSize);
T = initFields(T, T0);
b = sourceTerms(meshSize, L, P, k, h, A_t, T_ar, T_w);
x = linspace(0,L,meshSize);

res = inf;
count = 0;

while (res ~= 1)
    disp(count);
    count = count + 1;
    Tprime = T;
    T = (A*T+b);
    [res, beepBoop] = residual(eps, Tprime, T, 2);
    if count > 100
        res = 1;
    end
end


plot(x,T)










