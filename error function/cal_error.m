function error = cal_error( u,v,p )
% calculate the error between u,v,p and the exact solution
fun_u = @(x,y)( (1-cos(2*pi*x))*sin(2*pi*y) );
fun_v = @(x,y)(-(1-cos(2*pi*y))*sin(2*pi*x) );
fun_p = @(x,y)( (x^3)/3 - 1/12 );

n = size(p,1);

u_exact = zeros(n,n+1);
v_exact = zeros(n+1,n);
% p_exact = zeros(n,n);

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
% for i = 1:n
%     for j = 1:n
% p_exact(i,j) = fun_p( (2*j-1)/(2*n),(2*i-1)/(2*n) );
%     end
% end

% error = norm(u_exact-u,'fro')^2 + norm(v_exact-v,'fro')^2 + norm(p_exact-p,'fro');
error = norm(u_exact-u,'fro')^2 + norm(v_exact-v,'fro')^2;
% max(max(u_exact - u))
error = sqrt(error)/n;