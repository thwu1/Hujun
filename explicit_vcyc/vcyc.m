function [U_out,P_out] = vcyc(A,B,LU,LP,I,U_in,P_in,F_in)
% v1 = 2;
% v2 = 0;
% a = 1;
% method = @uzawa;
v1 = 10;
v2 = 0;
a = 0.2;
method = @inexact_uzawa;
n = floor(sqrt(length(F_in)/2))+1;
L = log2(n);
% Initialization
F = cell(1,L);
U = cell(1,L);
P = cell(1,L);
F{1} = F_in;
U{1} = U_in;
P{1} = P_in;
    [ U{1},P{1} ] = method( A{1},B{1},U{1},P{1},F{1},v1,a );
    for i = 2:L
        F_res = F{i-1} - A{i-1}*U{i-1} - B{i-1}*P{i-1};
        F{i} = I{i-1}*F_res;
        [U{i},P{i}] = initial(n/2^(i-1));
        [U{i},P{i}] = method(A{i},B{i},U{i},P{i},F{i},v1,a);
    end
    for i = L-1:-1:1
        U{i} = U{i} + LU{i}*U{i+1};
        P{i} = P{i} + LP{i}*P{i+1};
        [U{i},P{i}] = method(A{i},B{i},U{i},P{i},F{i},v2,0.001);
    end

U_out = U{1};
P_out = P{1};
