function [LU,LP] = prol(n)
% n->2n
J = sparse(2*n,n);
for j = 2:n-1
    J(2*j-1,j) = 3/4;
    J(2*j,j) = 3/4;
    J(2*j-2,j) = 1/4;
    J(2*j+1,j) = 1/4;
end
J(1,1) = 1/2;
J(2,1) = 3/4;
J(3,1) = 1/4;
J(end,end) = 1/2;
J(end-1,end) = 3/4;
J(end-2,end) = 1/4;
K = sparse(2*n-1,n-1);
for j = 1:n-1
    K(2*j,j) = 1;
    K(2*j-1,j) = 1/2;
    K(2*j+1,j) = 1/2;
end
L1 = kron(K,J);
L2 = kron(J,K);
LU = blkdiag(L1,L2);
L = sparse(n,n);
for i = 1:n
    L(2*i-1,i) = 1;
    L(2*i,i) = 1;
end
LP = kron(L,L);


