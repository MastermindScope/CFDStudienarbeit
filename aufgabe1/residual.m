function [endIt, res] = residual(eps, vectorOld, vectorNew, normP)
    if size(vectorOld) ~= size(vectorNew)
        error("Please give same size vectors");
    end
    res = norm(vectorOld-vectorNew,normP);
    if eps > res
        endIt = 1;
    else
        endIt = 0;
    end
end