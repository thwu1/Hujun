clear all;
v1 = 2;
v2 = 0;
a = 1;

% v1 = 5;
% v2 = 1;
% a = 0.1;

method = @uzawa;

n = 128;
L = log2(n);
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

itt = 0;
res_norm = 1;
while res_norm > 1e-8
    itt = itt + 1;
%     fprintf("Start iteration : %d  res_norm = %f\n",itt,res_norm);
tic
    [ U{1},P{1} ] = method( A{1},B{1},U{1},P{1},F{1},v1,a );
    toc
    for i = 2:L
        F_res = F{i-1} - A{i-1}*U{i-1} - B{i-1}*P{i-1};
        if i == 2
        res_norm = norm(F_res,2)/n^2;
        fprintf("%e\n",res_norm);
        if res_norm < 1e-8
            break;
        end
        end
        F{i} = I{i-1}*F_res;
        [U{i},P{i}] = initial(n/2^(i-1));
        [U{i},P{i}] = method(A{i},B{i},U{i},P{i},F{i},v1,a);
    end
    for i = L-1:-1:1
        U{i} = U{i} + LU{i}*U{i+1};
        P{i} = P{i} + LP{i}*P{i+1};
        [U{i},P{i}] = method(A{i},B{i},U{i},P{i},F{i},v2,0.001);
    end
end
% time = toc;
err = get_error(U{1},P{1});
% fprintf("Parameter:\nL = %d\nv1 = %d\nv2 = %d\nerror = %f\niteration: %d \nTime = %f\n",L,v1,v2,err,itt,time);
fprintf("N=%d Time=f err=%e Iter=%d\n",n,err,itt);

