function [U,P] = initial(n)
    U = zeros(2*n*(n-1),1);
    P = zeros(n*n,1);