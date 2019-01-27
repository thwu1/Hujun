function [u0,v0,p0] = vcycle_uzawa2(u0,v0,p0,v1,v2,L,a,n)
addpath('algorithm');
addpath('error function');
addpath('multigrid operation');
% Parameter
% v1 = 1;
% v2 = 100; 
% v3 = 100;
% L = 2   ; % layers of Multigrid method
% a = 0.02; % Parameter of uzawa smoothor
% n = 512;
% Initialization
f = cell(1,L);
g = cell(1,L);
u = cell(1,L);
v = cell(1,L);
p = cell(1,L);

[f{1},g{1}] = get_const(n);
% start main program
tic
err0 = 1;
err1 = 0;
err2 = 1;
iteration = 0;

% [u{L},v{L},p{L}] = initialize(n/2^(L-1));
% for i = L:-1:2
%    [ u{i},v{i},p{i} ] = uzawa( u{i},v{i},p{i},f{i},g{i},v2,a );
%    [ u{i-1},v{i-1},p{i-1} ] = lifting( u{i},v{i},p{i} );
% end % end for
% [ u{1},v{1},p{1} ] = uzawa( u{1},v{1},p{1},f{1},g{1},v2,a );

u{1} = u0;
v{1} = v0;
p{1} = p0;

while abs(err0 - err1)/n > 1e-5
err0 = err1;
iteration = iteration + 1;
fprintf("Start iteration : %d\n",iteration);
[ u{1},v{1},p{1} ] = uzawa( u{1},v{1},p{1},f{1},g{1},v1,a );
for i = 2:L
    [u{i},v{i},p{i}] = initialize(n/(2^(i-1)));
    [ f_out,g_out ] = cal_res(u{i-1},v{i-1},p{i-1},f{i-1},g{i-1});
    [ f{i},g{i} ] = restrict_fg( f_out,g_out );
    [ u{i},v{i},p{i} ] = uzawa( u{i},v{i},p{i},f{i},g{i},v1,a );
end % end for

for i = L:-1:2
   [ u{i},v{i},p{i} ] = uzawa( u{i},v{i},p{i},f{i},g{i},v2,a );
   [u_out,v_out,p_out] = lifting( u{i},v{i},p{i} );
   u{i-1} = u{i-1} + u_out;
   v{i-1} = v{i-1} + v_out;
   p{i-1} = p{i-1} + p_out;
end % end for

[ u{1},v{1},p{1} ] = uzawa( u{1},v{1},p{1},f{1},g{1},v2,a );

err2 = cal_error(u{1},v{1},p{1});
err1 = cal_res_norm(u{1},v{1},p{1});
fprintf("error:%f\nrh:%f\n",err2,err1);
end % end while
toc
fprintf("Parameter:\nN = %d\nL = %d\nv1 = %d\nv2 = %d\na = %f\nerror = %f\niteration: %d \n",n,L,v1,v2,a,err2,iteration);
u0 = u{1};
v0 = v{1};
p0 = p{1};
