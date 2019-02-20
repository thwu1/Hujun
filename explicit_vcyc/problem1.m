addpath('multigrid_fun');
addpath('error_fun');
addpath('algorithm');
clear all;
v1 = 1;
v2 = 0;

method = @dgs2;
for k = 6:11
n = 2^k;
L = log2(n)-1;
% Initialization
F = cell(1,L);
G = cell(1,L);
U = cell(1,L);
P = cell(1,L);
I = cell(1,L-1);
LU = cell(1,L-1);
LP = cell(1,L-1);
A = cell(1,L);
B = cell(1,L);
Ap = cell(1,L);
F{1} = get_F(n);
for i = 1:L
    [A{i},B{i}] = mac(n/2^(i-1));
    [U{i},P{i}] = initial(n/2^(i-1));
    Ap{i} = B{i}'*B{i};
end
for i = 1:L-1
    [LU{i},LP{i}] = prol(n/(2^i));
end
for i = 1:L-1
    I{i} = rest(n/2^i);
    Ig{i} = rest_G(n/2^i);
end

tic
itt = 0;
res_norm = 1;
while res_norm > 1e-8
    itt = itt + 1;
    G{1} = zeros(size(P{1}));
    [ U{1},P{1} ] = method( A{1},B{1},Ap{1},U{1},P{1},F{1},G{1},v1 );
    G{1} = G{1} - B{1}'*U{1};
    for i = 2:L
        F_res = F{i-1} - A{i-1}*U{i-1} - B{i-1}*P{i-1};
        if i == 2
        res_norm = norm(F_res,2)/norm(F{1},2);
%         fprintf("%e\n",res_norm);
        if res_norm < 1e-8
            break;
        end
        end
        F{i} = I{i-1}*F_res;
        G{i} = Ig{i-1}*G{i-1};
        [U{i},P{i}] = initial(n/2^(i-1));
        [U{i},P{i}] = method(A{i},B{i},Ap{i},U{i},P{i},F{i},G{i},v1);
        G{i} = G{i} - B{i}'*U{i};
    end
    for i = L-1:-1:1
        U{i} = U{i} + LU{i}*U{i+1};
        P{i} = P{i} + LP{i}*P{i+1};
        [U{i},P{i}] = method(A{i},B{i},Ap{i},U{i},P{i},F{i},G{i},v2);
    end
end
time = toc;
err = get_error(U{1},P{1});
fprintf("N=%d Time=%f err=%e Iter=%d\n",n,time,err,itt);
end
