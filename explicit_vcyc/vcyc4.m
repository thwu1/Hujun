function [U_out] = vcyc4(U_in,F_in)
v1 = 2;
v2 = 0;
a = 1;
method = @GS;
% v1 = 10;
% v2 = 0;
% a = 0.2;
% method = @inexact_uzawa;
n = floor(sqrt(length(F_in)/2))+1;
L = log2(n);
% Initialization
F = cell(1,L);
U = cell(1,L);
I = cell(1,L-1);
LU = cell(1,L-1);
A = cell(1,L);
B = cell(1,L);
F{1} = F_in;
for i = 1:L
    [A{i},B{i}] = mac(n/2^(i-1));
    [U{i},~] = initial(n/2^(i-1));
end
for i = 1:L-1
    [LU{i},~] = prol(n/(2^i));
end
for i = 1:L-1
    I{i} = rest(n/2^i);
end
U{1} = U_in;
    U{1} = method( A{1},B{1},U{1},P{1},F{1} );
    for i = 2:L
        F_res = F{i-1} - A{i-1}*U{i-1};
        F{i} = I{i-1}*F_res;
        U{i} = initial(n/2^(i-1));
        U{i} = method(A{i},B{i},U{i},P{i},F{i});
    end
    for i = L-1:-1:1
        U{i} = U{i} + LU{i}*U{i+1};
        U{i} = method(A{i},B{i},U{i},P{i},F{i});
    end
U_out = U{1};
