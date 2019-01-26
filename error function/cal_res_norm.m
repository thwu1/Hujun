function res_norm = cal_res_norm(u,v,p)
% calculate norm(rh,'fro')
% [u0,v0,p0] is exact solution
n = size(p,1);
fun_b = @(x)( 2*n*sin(pi/n)*( 1-cos(2*pi*x) ) );
fun_l = @(y)(-2*n*sin(pi/n)*( 1-cos(2*pi*y) ) );
fun_f = @(x,y)(-4*pi^2*(2*cos(2*pi*x)-1) * sin(2*pi*y) + x^2 );
fun_g = @(x,y)( 4*pi^2*(2*cos(2*pi*y)-1) * sin(2*pi*x) );
fun_u = @(x,y)( (1-cos(2*pi*x))*sin(2*pi*y) );
fun_v = @(x,y)(-(1-cos(2*pi*y))*sin(2*pi*x) );
fun_p = @(x,y)( (1/3)*x^3 - 1/12 );
f0 = zeros( n,n+1 );
g0 = zeros( n+1,n );
p0 = zeros(n,n);
u0 = zeros(n,n+1);
v0 = zeros(n+1,n);

for i = 1:n
    for j = 1:n+1
u0(i,j) = fun_u( (j-1)/n,(2*i-1)/(2*n) );
    end
end

for i = 1:n+1
    for j = 1:n
v0(i,j) = fun_v( (2*j-1)/(2*n),(i-1)/n );
    end
end
for i = 1:n
    for j = 1:n
p0(i,j) = fun_p( (2*j-1)/(2*n),(2*i-1)/(2*n) );
    end
end

% calculate b,t
for i = 1:n+1
b(i) = fun_b( (i-1)/n );
end
t = b;
% calculate l,r
for i = 1:n+1
l(i) = fun_l( (i-1)/n );
end
r = l;
% calculate f
for i = 1:n
    for j = 1:n+1
f0(i,j) = fun_f( (j-1)/n,(2*i-1)/(2*n) );
    end
end
% add boundary condition to f
for j = 1:n+1
f0(1,j) = f0(1,j) - n*b(j);
f0(n,j) = f0(n,j) + n*t(j);
end
% calculate g
for i = 1:n+1
    for j = 1:n
g0(i,j) = fun_g( (2*j-1)/(2*n),(i-1)/n );
    end
end
% add boundary condition to g
for i = 1:n+1
g0(i,1) = g0(i,1) - n*l(i);
g0(i,n) = g0(i,n) + n*r(i);
end


[r1,r2] = cal_res(u,v,p,f0,g0);
% [r3,r4] = cal_res(u0,v0,p0,f0,g0);
% fprintf("rh   : %f\nexact: %f\n",norm([r1,r2'],'fro'),norm([r3,r4'],'fro'));
res_norm = norm([r1,r2'],'fro');