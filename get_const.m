function [ f,g ] = get_const( n )
%
fun_f = @(x,y)(-4*pi^2*(2*cos(2*pi*x)-1) * sin(2*pi*y) + x^2 );
fun_g = @(x,y)( 4*pi^2*(2*cos(2*pi*y)-1) * sin(2*pi*x) );
fun_b = @(x)( 2*pi*( 1-cos(2*pi*x) ) );
fun_l = @(y)(-2*pi*( 1-cos(2*pi*y) ) );

f = zeros( n,n+1 );
g = zeros( n+1,n );
b = zeros( n+1,1 );
t = zeros( n+1,1 );
l = zeros( n+1,1 );
r = zeros( n+1,1 );

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
f(i,j) = fun_f( (j-1)/n,(2*i-1)/(2*n) );
    end
end
% add boundary condition to f
for j = 1:n+1
f(1,j) = f(1,j) - n*b(j);
f(n,j) = f(n,j) + n*t(j);
end
% calculate g
for i = 1:n+1
    for j = 1:n
g(i,j) = fun_g( (2*j-1)/(2*n),(i-1)/n );
    end
end
% add boundary condition to g
for i = 1:n+1
g(i,1) = g(i,1) - n*l(i);
g(i,n) = g(i,n) + n*r(i);
end

