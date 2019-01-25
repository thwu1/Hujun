v1 = 15;
v2 = 1;
L = 1;
n = 64;
max_iteration = 100;

[ f,g ] = get_const(n);
[ u,v,p,~,~ ] = initialize(n);

while max_iteration
    max_iteration = max_iteration - 1;
    [ u,v,p ] = DGS( u,v,p,f,g,v1 );
fprintf("error:%f\n",cal_error(u,v,p));
end % end while


