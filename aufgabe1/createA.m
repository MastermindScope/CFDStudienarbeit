function A = createA(meshSize, k, h, P, L, A_t, T_w, T_ar)
    dx = L/meshSize;
    m = h*P/(k*A_t);
    A = zeros(meshSize);
    A = A+diag(-ones(meshSize-1,1)/dx,1);
    A = A+diag(-ones(meshSize-1,1)/dx,-1);
    A = A+diag(ones(meshSize,1)*(2/dx+m*dx),0);
    A(1,1) = 3/dx+m*dx;
    A(1,2) = -1/dx;
    A(meshSize,meshSize) = h/k/(1+h*dx/(2*k));
    A(meshSize,meshSize-1) = -1/dx;
end