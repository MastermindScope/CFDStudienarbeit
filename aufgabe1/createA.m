function A = createA(meshSize, k, h, P, L, A_t)
    %define dx
    dx = L/meshSize;
    %define m
    m = h*P/(k*A_t);
    %initialize A as sparse matrix of the right size
    A = sparse(meshSize,meshSize);
    %create a vector of ones
    ons = ones(meshSize-1,1);
    %create a vector with ascending integers
    rowCol = linspace(1,meshSize-1,meshSize-1)';
    %set values of the sparse matrix
    A = A+sparse(rowCol,rowCol+ons,-ons/dx,meshSize,meshSize);
    A = A+sparse(rowCol+ons,rowCol,-ons/dx,meshSize,meshSize);
    %another vector of ones
    ons = ones(meshSize,1);
    %another vector of ascending integers
    rowCol = linspace(1,meshSize,meshSize)';
    %set further values in the matrix
    A = A+sparse(rowCol,rowCol,ons*(2/dx+m*dx));
    %replace wrong entries
    A(1,1) = 3/dx+m*dx;
    A(1,2) = -1/dx;
    A(meshSize,meshSize) = h/k/(1+h*dx/(2*k))+1/dx+m*dx;
    A(meshSize,meshSize-1) = -1/dx;
end