% This file is used to test uzawa method
v1 = 10;
n = 64;
max_iteration = 500;

[f,g] = get_const(n);
[u,v,p,~,~] = initialize(n);
tic
while max_iteration
    max_iteration = max_iteration - 1;
    [u,v,p] = uzawa(u,v,p,f,g,v1,-0.002);
  fprintf("error:%f\n",cal_error(u,v,p));
end % end while
toc

