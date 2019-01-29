function  [u_out,v_out,p_out] = vcycle_uzawa(v1,v2,L,a,n)
addpath('algorithm');
addpath('error function');
addpath('multigrid operation');
% Parameter
% v1 = 1;
% v2 = 300; 
% L = 5   ; % layers of Multigrid method
% a = 0.01; % Parameter of uzawa smoothor
% n = 128;
% Initialization
f = cell(1,L);
g = cell(1,L);
u = cell(1,L);
v = cell(1,L);
p = cell(1,L);

[u{1},v{1},p{1}] = initialize(n);
for i = 1:L
    n0 = n/( 2^(i-1) );
    [f{i},g{i}] = get_const(n0);
end
% start main program
tic
err0 = 1;
err1 = 0;
err2 = 1;
iteration = 0;

[u{L},v{L},p{L}] = initialize(n/2^(L-1));
for i = L:-1:2
   [ u{i},v{i},p{i} ] = uzawa( u{i},v{i},p{i},f{i},g{i},v2,a );
   [ u{i-1},v{i-1},p{i-1} ] = lifting( u{i},v{i},p{i} );
end % end for
[ u{1},v{1},p{1} ] = uzawa( u{1},v{1},p{1},f{1},g{1},v2,a );


while abs(err2)/n > 1e-5
% err0 = err1;
iteration = iteration + 1;
fprintf("Start iteration : %d\n",iteration);
[ u{1},v{1},p{1} ] = implicit_uzawa( u{1},v{1},p{1},f{1},g{1},v1,a );
for i = 2:L
    [ u{i},v{i},p{i} ] = restrict( u{i-1},v{i-1},p{i-1} );
    [ u{i},v{i},p{i} ] = implicit_uzawa( u{i},v{i},p{i},f{i},g{i},v1,a );
end % end for

for i = L:-1:2
   [ u{i},v{i},p{i} ] = implicit_uzawa( u{i},v{i},p{i},f{i},g{i},v2,a );
   [ u{i-1},v{i-1},p{i-1} ] = lifting( u{i},v{i},p{i} );
end % end for

[ u{1},v{1},p{1} ] = implicit_uzawa( u{1},v{1},p{1},f{1},g{1},v2,a );

err2 = cal_error(u{1},v{1},p{1});
% err1 = cal_res_norm(u{1},v{1},p{1});
fprintf("error:%f\n",err2);
end % end while
toc
fprintf("Parameter:\nN = %d\nL = %d\nv1 = %d\nv2 = %d\na = %f\nerror = %f\niteration: %d \n",n,L,v1,v2,a,err2,iteration);
u_out = u{1};
v_out = v{1};
p_out = p{1};
