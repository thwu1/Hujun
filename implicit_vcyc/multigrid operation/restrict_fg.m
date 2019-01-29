function [ f_out,g_out ] = restrict_fg( f,g )
% restrict operator for the multigrid method

n = size( f,1 );

f_out = zeros(n/2,n/2+1);
g_out = zeros(n/2+1,n/2);

f_out(:,2:end-1) = (1/4)*( f(1:2:end,3:2:end-2) + f(2:2:end,3:2:end-2) ) + (1/8)*( f(1:2:end,2:2:end-3) + f(1:2:end,4:2:end-1) + f(2:2:end,2:2:end-3) + f(2:2:end,4:2:end-1) );
g_out(2:end-1,:) = (1/4)*( g(3:2:end-2,1:2:end) + g(3:2:end-2,2:2:end) ) + (1/8)*( g(2:2:end-3,1:2:end) + g(4:2:end-1,1:2:end) + g(2:2:end-3,2:2:end) + g(4:2:end-1,2:2:end) );

end