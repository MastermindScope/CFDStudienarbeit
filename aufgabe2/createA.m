function A = createA(meshSize, cp, h, P, L, A_t, dt, k, Mp)
    %define dx
    dx = L/meshSize;
    %define m
    n = h*P/(cp*A_t);
    %initialize A as sparse matrix of the right size
    A = sparse(meshSize,meshSize);
    %create a vector of ones
    ons = ones(meshSize-1,1);
    %create a vector with ascending integers
    rowCol = linspace(1,meshSize-1,meshSize-1)';
    %set values of the sparse matrix
    A = A+sparse(rowCol,rowCol+ons,-ons/dx*k/cp,meshSize,meshSize);
    A = A+sparse(rowCol+ons,rowCol,-ons/dx*k/cp,meshSize,meshSize);
    %another vector of ones
    ons = ones(meshSize,1);
    %another vector of ascending integers
    rowCol = linspace(1,meshSize,meshSize)';
    %set further values in the matrix
    A = A+sparse(rowCol,rowCol,ons*(2*k/dx/cp+n*dx+Mp/dt));
    %replace wrong entries
    A(1,1) = 3*k/(dx*cp)+n*dx+Mp/dt;
    A(1,2) = -k/(dx*cp);
    A(meshSize,meshSize) = Mp/dt+2*k*h/(cp*(2*k+h*dx))+k/cp/dx+n*dx;
    A(meshSize,meshSize-1) = -k/dx/cp;
end