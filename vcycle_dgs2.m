addpath('algorithm');
addpath('error function');
addpath('multigrid operation');
% Parameter
v1 = 1;
v2 = 200; 
L = 1   ; % layers of Multigrid method
n = 128;
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
   [ u{i},v{i},p{i} ] = implicit_dgs( u{i},v{i},p{i},f{i},g{i},v2 );
   [ u{i-1},v{i-1},p{i-1} ] = lifting( u{i},v{i},p{i} );
end % end for
[ u{1},v{1},p{1} ] = implicit_dgs( u{1},v{1},p{1},f{1},g{1},v2 );


while abs(err0 - err1)/n > 1e-8
err0 = err1;
iteration = iteration + 1;
fprintf("Start iteration : %d\n",iteration);
[ u{1},v{1},p{1} ] = implicit_dgs( u{1},v{1},p{1},f{1},g{1},v1 );
for i = 2:L
    [u{i},v{i},p{i}] = initialize(n/(2^(i-1)));
    [ f_out,g_out ] = cal_res(u{i-1},v{i-1},p{i-1},f{i-1},g{i-1});
    [ f{i},g{i} ] = restrict_fg( f_out,g_out );
    [ u{i},v{i},p{i} ] = implicit_dgs( u{i},v{i},p{i},f{i},g{i},v1 );
end % end for

for i = L:-1:2
   [ u{i},v{i},p{i} ] = implicit_dgs( u{i},v{i},p{i},f{i},g{i},v2 );
   [u_out,v_out,p_out] = lifting( u{i},v{i},p{i} );
   u{i-1} = u{i-1} + u_out;
   v{i-1} = v{i-1} + v_out;
   p{i-1} = p{i-1} + p_out;
end % end for

[ u{1},v{1},p{1} ] = implicit_dgs( u{1},v{1},p{1},f{1},g{1},v2 );

err2 = cal_error(u{1},v{1},p{1});
err1 = cal_res_norm(u{1},v{1},p{1});
fprintf("error:%f\nrh:%f\n",err2,err1);
end % end while
toc
fprintf("Parameter:\nN = %d\nL = %d\nv1 = %d\nv2 = %d\na = %f\nerror = %f\niteration: %d \n",n,L,v1,v2,a,err2,iteration);

