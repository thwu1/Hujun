function [ f_out,g_out ] = cal_res( u,v,p,f,g )
%
n = size(p,1);
% update f
f_out(:,1) = zeros(n,1);
f_out(:,n+1) = zeros(n,1);
        i = 1;
        for j = 2:n
            f_out(i,j) = f(i,j) + n^2 * ( u(i+1,j) + u(i,j+1) + u(i,j-1) - 3*u(i,j) ) - n * ( p(i,j) - p(i,j-1) );
        end 
        i = n;
        for j = 2:n
            f_out(i,j) = f(i,j) + n^2 * ( u(i-1,j) + u(i,j+1) + u(i,j-1) - 3*u(i,j) ) - n * ( p(i,j) - p(i,j-1) );
        end
        for i = 2:n-1
            for j = 2:n
            f_out(i,j) = f(i,j) + n^2 * ( u(i+1,j) + u(i-1,j) + u(i,j+1) + u(i,j-1) - 4*u(i,j) ) - n * ( p(i,j) - p(i,j-1) );
            end
        end

% update g
g_out(1,:) = zeros(1,n);
g_out(n+1,:) = zeros(1,n);
        j = 1;
        for i = 2:n
            g_out(i,j) = g(i,j) + n^2 * ( v(i,j+1) + v(i-1,j) + v(i+1,j) - 3*v(i,j) ) - n * ( p(i,j) - p(i-1,j) );
        end
        j = n;
        for i = 2:n
            g_out(i,j) = g(i,j) + n^2 * ( v(i,j-1) + v(i-1,j) + v(i+1,j) - 3*v(i,j) ) - n * ( p(i,j) - p(i-1,j) );
        end 
        for i = 2:n
            for j = 2:n-1
            g_out(i,j) = g(i,j) + n^2 * ( v(i,j+1) + v(i,j-1) + v(i-1,j) + v(i+1,j) - 4*v(i,j) ) - n * ( p(i,j) - p(i-1,j) );
            end 
        end