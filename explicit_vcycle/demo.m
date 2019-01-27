clear all;
v1 = 5;
v2 = 1;
L = 7; % layers of Multigrid method
n = 128;
% Initialization
F = cell(1,L);
U = cell(1,L);
P = cell(1,L);
I = cell(1,L-1);
LU = cell(1,L-1);
LP = cell(1,L-1);
A = cell(1,L);
B = cell(1,L);
F{1} = get_F(n);
for i = 1:L
    [A{i},B{i}] = mac(n/2^(i-1));
    [U{i},P{i}] = initial(n/2^(i-1));
end
for i = 1:L-1
    [LU{i},LP{i}] = prol(n/(2^i));
end
for i = 1:L-1
    I{i} = rest(n/2^i);
end

tic
itt = 0;
res_norm = 1;
while res_norm > 1e-8
    itt = itt + 1;
    fprintf("Start iteration : %d  res_norm = %f\n",itt,res_norm);
    [ U{1},P{1} ] = uzawa( A{1},B{1},U{1},P{1},F{1},v1,0.1 );
    for i = 2:L
        F_res = F{i-1} - A{i-1}*U{i-1} - B{i-1}*P{i-1};
        if i == 2
        res_norm = norm(F_res,2)/n^2;
        end
        F{i} = I{i-1}*F_res;
        [U{i},P{i}] = initial(n/2^(i-1));
        [U{i},P{i}] = uzawa(A{i},B{i},U{i},P{i},F{i},v1,0.1);
    end
    for i = L-1:-1:1
        U1 = LU{i}*U{i+1};
        P1 = LP{i}*P{i+1};
        U{i} = U{i} + U1;
        P{i} = P{i} + P1;
        [U{i},P{i}] = uzawa(A{i},B{i},U{i},P{i},F{i},v1,0.01);
    end
end
time = toc;
err = get_error(U{1},P{1});
fprintf("Parameter:\nL = %d\nv1 = %d\nv2 = %d\nerror = %f\niteration: %d \nTime = %f\n",L,v1,v2,err,itt,time);

