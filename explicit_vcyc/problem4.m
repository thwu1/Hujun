clear all;
for k = 6:11
global I LU LP A B F n p;
n = 2^k;
L = log2(n);
a = 1.1;
v = 2;
p = 5;

% Initialization
I = cell(1,L-1);
LU = cell(1,L-1);
LP = cell(1,L-1);
A = cell(1,L);
B = cell(1,L);
F = get_F(n);
for i = 1:L
    [A{i},B{i}] = mac(n/2^(i-1));
end
for i = 1:L-1
    [LU{i},LP{i}] = prol(n/(2^i));
end
for i = 1:L-1
    I{i} = rest(n/2^i);
end
[U,P] = initial(n);



tic
itt = 0;
res_norm = 1;
while res_norm > 1e-8
    itt = itt + 1;
    for i = 1:v
U = vcyc(U,F-B{1}*P);
    end
res_norm = norm(F - A{1}*U - B{1}*P,2)/n^2;
if res_norm < 1e-8
    break;
end
% tic
% [U,P] = method(A{1},B{1},U,P,F,v,a);
% toc
P = P + a*(B{1}'*U);
end
time = toc;
err = get_error(U,P);
fprintf("N=%d Time=%f err=%e Iter=%d a:%f v:%d p:%d\n",n,time,err,itt,a,v,p);

end



function U_out = vcyc(U_in,F_in)
global A LU I n p;
v1 = 1;
v2 = 0;
L = log2(n);
% Initialization
F = cell(1,L);
U = cell(1,L);
F{1} = F_in;
U{1} = U_in;
    U{1} = GS4(A{1},U{1},F{1},v1);
    for i = 2:p
        F_res = F{i-1} - A{i-1}*U{i-1};
        F{i} = I{i-1}*F_res;
        U{i} = zeros(2*(n/2^(i-1))*(n/2^(i-1)-1),1);
        if i < p
        U{i} = GS4(A{i},U{i},F{i},v1);
        else 
            U{i} = A{i}\F{i};
        end
    end
    for i = p-1:-1:1
        U{i} = U{i} + LU{i}*U{i+1};
        U{i} = GS4(A{i},U{i},F{i},v2);
    end
U_out = U{1};
end

function U = GS4(A,U,F,v)
for i = 1:v
    U = U + tril(A)\(F - A*U);
end
end
