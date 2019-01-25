function [ f_out,g_out ] = cal_res( u,v,p,f,g )
%
n = size(p,1);

% % update b
% b = n*( u(2,:) - u(1,:) );
% % update t
% t = n*( u(n,:) - u(n-1,:) );
% % update l
% l = n*( v(:,2) -v(:,1) );
% % update r
% r = n*( v(:,n) - v(:,n-1) );

% update f
f_out(:,1) = zeros(n,1);
f_out(:,n+1) = zeros(n,1);
for i = 1:n
    for j = 2:n
        
        if i == 1
            f_out(i,j) = f(i,j) + n^2 * ( u(i+1,j) + u(i,j+1) + u(i,j-1) - 3*u(i,j) ) - n * ( p(i,j) - p(i,j-1) );
        end % end if
        if i == n
            f_out(i,j) = f(i,j) + n^2 * ( u(i-1,j) + u(i,j+1) + u(i,j-1) - 3*u(i,j) ) - n * ( p(i,j) - p(i,j-1) );
        end % end if
        if i ~= 1 && i ~= n
            f_out(i,j) = f(i,j) + n^2 * ( u(i+1,j) + u(i-1,j) + u(i,j+1) + u(i,j-1) - 4*u(i,j) ) - n * ( p(i,j) - p(i,j-1) );
        end % end if
        
    end % end for
end % end for

% update g
g_out(1,:) = zeros(1,n);
g_out(n+1,:) = zeros(1,n);
for i = 2:n
    for j = 1:n
        
        if j == 1
            g_out(i,j) = g(i,j) + n^2 * ( v(i,j+1) + v(i-1,j) + v(i+1,j) - 3*v(i,j) ) - n * ( p(i,j) - p(i-1,j) );
        end % end if
        if j == n
            g_out(i,j) = g(i,j) + n^2 * ( v(i,j-1) + v(i-1,j) + v(i+1,j) - 3*v(i,j) ) - n * ( p(i,j) - p(i-1,j) );
        end % end if
        if j ~= 1 && j ~= n
            g_out(i,j) = g(i,j) + n^2 * ( v(i,j+1) + v(i,j-1) + v(i-1,j) + v(i+1,j) - 4*v(i,j) ) - n * ( p(i,j) - p(i-1,j) );
        end % end if
        
    end % end for
end % end for