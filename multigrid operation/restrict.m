function [ u_out,v_out,p_out ] = restrict( u,v,p )
% restrict operator for the multigrid method

n = size( u,1 );

u_out = zeros(n/2,n/2+1);
v_out = zeros(n/2+1,n/2);
p_out = zeros(n/2,n/2);

u_out(:,2:end-1) = (1/4)*( u(1:2:end,3:2:end-2) + u(2:2:end,3:2:end-2) ) + (1/8)*( u(1:2:end,2:2:end-3) + u(1:2:end,4:2:end-1) + u(2:2:end,2:2:end-3) + u(2:2:end,4:2:end-1) );
v_out(2:end-1,:) = (1/4)*( v(3:2:end-2,1:2:end) + v(3:2:end-2,2:2:end) ) + (1/8)*( v(2:2:end-3,1:2:end) + v(4:2:end-1,1:2:end) + v(2:2:end-3,2:2:end) + v(4:2:end-1,2:2:end) );
p_out = (1/4)*( p(1:2:end,1:2:end) + p(1:2:end,2:2:end) + p(2:2:end,1:2:end) + p(2:2:end,2:2:end) );

end