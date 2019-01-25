addpath('algorithm');
addpath('error function');
addpath('multigrid operation');
% Parameter
v1 = 20;
v2 = 20; 
L = 5; % layers of Multigrid method
a = -0.01; % Parameter of uzawa smoothor
n = 256;
% Initialization
f = cell(1,L);
g = cell(1,L);
u = cell(1,L);
v = cell(1,L);
p = cell(1,L);

[u{1},v{1},p{1},~,~] = initialize(n);
for i = 1:L
    n0 = n/( 2^(i-1) );
    [f{i},g{i}] = get_const(n0);
end

tic
err0 = 0;
err1 = 1;
iteration = 0;
[u{L},v{L},p{L},~,~] = initialize(n/2^(L-1));
for i = L:-1:2
   [ u{i},v{i},p{i} ] = uzawa( u{i},v{i},p{i},f{i},g{i},v2,a );
   [ u{i-1},v{i-1},p{i-1} ] = lifting( u{i},v{i},p{i} );
end % end for
[ u{1},v{1},p{1} ] = uzawa( u{1},v{1},p{1},f{1},g{1},v2,a );

while abs(err0 - err1) > 1e-8
    err0 = err1;
%      [r1,r2] = cal_res( u{1},v{1},p{1},f{1},g{1} );
%     if norm([r1;r2']/n^2,'fro') < 1e-2
%         break;
%     end
%  fprintf("r_h:%f\n",norm([r1;r2']/n^2,'fro'));
iteration = iteration + 1;
fprintf("Start iteration : %d\n",iteration);
[ u{1},v{1},p{1} ] = uzawa( u{1},v{1},p{1},f{1},g{1},v1,a );
for i = 2:L
    [ u{i},v{i},p{i} ] = restrict( u{i-1},v{i-1},p{i-1} );
    [ u{i},v{i},p{i} ] = uzawa( u{i},v{i},p{i},f{i},g{i},v1,a );
%     [ f0,g0 ] = cal_res( u{i},v{i},p{i},f{i},g{i} );
end % end for

for i = L:-1:2
   [ u{i},v{i},p{i} ] = uzawa( u{i},v{i},p{i},f{i},g{i},v2,a );
   [ u{i-1},v{i-1},p{i-1} ] = lifting( u{i},v{i},p{i} );
end % end for
[ u{1},v{1},p{1} ] = uzawa( u{1},v{1},p{1},f{1},g{1},v2,a );
err1 = cal_error(u{1},v{1},p{1});
end % end while
toc
fprintf("Parameter:\nL = %d\nv1 = %d\nv2 = %d\na = %f\nerror = %f\niteration: %d \n",L,v1,v2,a,err1,iteration);

