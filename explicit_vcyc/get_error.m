function err = get_error(U,P)
n = sqrt(length(P));
u = U(1:n^2-n);
v = U(n^2-n+1:end);
u = [zeros(n,1),reshape(u,n,n-1),zeros(n,1)];
v = [zeros(1,n);reshape(v,n-1,n);zeros(1,n)];
p = reshape(P,n,n);
err = cal_error(u,v,p);