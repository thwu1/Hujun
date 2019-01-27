function [ u,v ] = implicit_gs( u,v,p,f,g )
% first step of DGS
n = size(p,1);
h = 1/n;
% update u
for i = 2:n-1
    for j = 2:n
            u(i,j) = (1/4)*( h^2 * f(i,j) - h*( p(i,j) - p(i,j-1) ) + u(i,j+1) + u(i,j-1) + u(i-1,j) + u(i+1,j) );
    end
end

i = 1;
for j = 2:n
            u(i,j) = (1/3)*( h^2 * f(i,j) - h*( p(i,j) - p(i,j-1) ) + u(i,j+1) + u(i,j-1) + u(i+1,j) );
end

i = n;
for j = 2:n
            u(i,j) = (1/3)*( h^2 * f(i,j) - h*( p(i,j) - p(i,j-1) ) + u(i,j+1) + u(i,j-1) + u(i-1,j) );
end

i = n;
for j = n:-1:2
            u(i,j) = (1/3)*( h^2 * f(i,j) - h*( p(i,j) - p(i,j-1) ) + u(i,j+1) + u(i,j-1) + u(i-1,j) );
end
i = 1;
for j = n:-1:2
            u(i,j) = (1/3)*( h^2 * f(i,j) - h*( p(i,j) - p(i,j-1) ) + u(i,j+1) + u(i,j-1) + u(i+1,j) );
end
for i = n-1:-1:2
    for j = n:-1:2
            u(i,j) = (1/4)*( h^2 * f(i,j) - h*( p(i,j) - p(i,j-1) ) + u(i,j+1) + u(i,j-1) + u(i-1,j) + u(i+1,j) );
    end
end



% update v
for i = 2:n
    for j = 2:n-1
            v(i,j) = (1/4)*( h^2 * g(i,j) - h*( p(i,j) - p(i-1,j) ) + v(i,j+1) + v(i,j-1) + v(i-1,j) + v(i+1,j) );
    end 
end

j = 1;
for i = 2:n
            v(i,j) = (1/3)*( h^2 * g(i,j) - h*( p(i,j) - p(i-1,j) ) + v(i,j+1) + v(i-1,j) + v(i+1,j) );
end
j = n;
for i = 2:n
            v(i,j) = (1/3)*( h^2 * g(i,j) - h*( p(i,j) - p(i-1,j) ) + v(i,j-1) + v(i-1,j) + v(i+1,j) );
end
for i = n:-1:2
            v(i,j) = (1/3)*( h^2 * g(i,j) - h*( p(i,j) - p(i-1,j) ) + v(i,j-1) + v(i-1,j) + v(i+1,j) );
end
j = 1;
for i = n:-1:2
            v(i,j) = (1/3)*( h^2 * g(i,j) - h*( p(i,j) - p(i-1,j) ) + v(i,j+1) + v(i-1,j) + v(i+1,j) );
end 

for i = 2:n
    for j = 2:n-1
            v(i,j) = (1/4)*( h^2 * g(i,j) - h*( p(i,j) - p(i-1,j) ) + v(i,j+1) + v(i,j-1) + v(i-1,j) + v(i+1,j) );
    end 
end
