v1 = 5;
v2 = 1;
L = 1;
n = 64;
max_iteration = 100;
f = cell(1,L+1);
g = cell(1,L+1);
u = cell(1,L+1);
v = cell(1,L+1);
p = cell(1,L+1);

[ f{1},g{1} ] = get_const(n);

for i = 1:L+1
    n0 = n/( 2^(i-1) );
    [ u{i},v{i},p{i},~,~ ] = initialize( n0 );
end

while max_iteration
    
    [ r1,r2 ] = cal_res( u{1},v{1},p{1},f{1},g{1} );
    
%     if norm([r1;r2'],'fro') < 1e-2
%         break;
%     end
fprintf("r_h:%f\n",norm([r1;r2'],'fro'));
    max_iteration = max_iteration - 1;

for i = 1:L
    [ u{i},v{i},p{i} ] = DGS( u{i},v{i},p{i},f{i},g{i},v1 );
    [ f0,g0 ] = cal_res( u{i},v{i},p{i},f{i},g{i} );
    [ f{i+1},g{i+1} ] = restrict( f0,g0 );
end % end for

for i = L:-1:1
   [ u_out,v_out,p_out ] = lifting( u{i+1},v{i+1},p{i+1} );
   u{i} = u{i} + u_out;
   v{i} = v{i} + v_out;
   p{i} = p{i} + p_out;
   [ u{i},v{i},p{i} ] = DGS( u{i},v{i},p{i},f{i},g{i},v2 );
end % end for
fprintf("error:%f\n",cal_error(u{1},v{1},p{1}));
end % end while


