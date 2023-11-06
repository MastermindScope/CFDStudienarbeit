function b = sourceTerms(meshSize, L, P, k, h, A_t, T_ar, T_w)
    b = ones(meshSize,1);
    dx = L/meshSize;
    m = h*P/(k*A_t);
    b = b*m*T_ar*dx;
    b(1) = m*T_ar*dx+2*T_w/dx;
    if meshSize > 1
        b(meshSize) = m*T_ar*dx+h/k*T_ar/(1+h*dx/(2*k));
    end
end