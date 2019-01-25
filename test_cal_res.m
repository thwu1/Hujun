% This file is used to test cal_res
% It returns cal_res(u_exact,v_exact,p_exact,f,g)
n = 64;
[f,g] = get_const(n);

u_exact = zeros(n,n+1);
v_exact = zeros(n+1,n);
p_exact = zeros(n,n);

for i = 1:n
    for j = 1:n+1
u_exact(i,j) = fun_u( (j-1)/n,(2*i-1)/(2*n) );
    end
end

for i = 1:n+1
    for j = 1:n
v_exact(i,j) = fun_v( (2*j-1)/(2*n),(i-1)/n );
    end
end

for i = 1:n
    for j = 1:n
p_exact(i,j) = fun_p( (2*j-1)/(2*n),(2*i-1)/(2*n) );
    end
end

[r1,r2] = cal_res(u_exact,v_exact,p_exact,f,g);
norm([r1,r2'],'fro')/n

