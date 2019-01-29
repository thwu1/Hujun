function [ u_out,v_out,p_out ] = lifting( u,v,p )
% lifting operator for the multigrid method
n = size(p,1);
p_out = zeros(2*n,2*n);
u_out = zeros(2*n,2*n+1);
v_out = zeros(2*n+1,2*n);

p_out(1:2:end,1:2:end) = p;
p_out(2:2:end,1:2:end) = p;
p_out(1:2:end,2:2:end) = p;
p_out(2:2:end,2:2:end) = p;

u_out(1,1:2:end) = (1/2)*u(1,:);
u_out(end,1:2:end) = (1/2)*u(end,:);
u_out(2:2:end-1,1:2:end) = (3/4)*u(1:end-1,:) + (1/4)*u(2:end,:);
u_out(3:2:end,1:2:end) = (1/4)*u(1:end-1,:) + (3/4)*u(2:end,:);
u_out(:,2:2:end-1) = (1/2)*( u_out(:,1:2:end-1) + u_out(:,3:2:end) );

v_out(1:2:end,1) = (1/2)*v(:,1);
v_out(1:2:end,end) = (1/2)*v(:,end);
v_out(1:2:end,2:2:end-1) = (3/4)*v(:,1:end-1) + (1/4)*v(:,2:end);
v_out(1:2:end,3:2:end) = (1/4)*v(:,1:end-1) + (3/4)*v(:,2:end);
v_out(2:2:end-1,:) = (1/2)*( v_out(1:2:end-1,:) + v_out(3:2:end,:) );
