function err = get_error(U,P)
n = sqrt(length(P));
u = U(1:n^2-n);
v = U(n^2-n+1:end);
u = [zeros(n,1),reshape(u,n,n-1),zeros(n,1)];
v = [zeros(1,n);reshape(v,n-1,n);zeros(1,n)];
p = reshape(P,n,n);
err = cal_error(u,v,p);

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
for i = 1:n
    for j = 1:n
p_exact(i,j) = fun_p( (2*j-1)/(2*n),(2*i-1)/(2*n) );
    end
end
% imshow(abs(p-p_exact)/max(max(abs(p-p_exact))));
error = norm(u_exact-u,'fro')^2 + norm(v_exact-v,'fro')^2;
error = sqrt(error)/n;