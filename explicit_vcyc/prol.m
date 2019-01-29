function [LU,LP] = prol(n)
% n->2n
A = sparse(2*n,n);
for j = 2:n-1
    A(2*j-1,j) = 3/4;
    A(2*j,j) = 3/4;
    A(2*j-2,j) = 1/4;
    A(2*j+1,j) = 1/4;
end
A(1,1) = 1/2;
A(2,1) = 3/4;
A(3,1) = 1/4;
A(end,end) = 1/2;
A(end-1,end) = 3/4;
A(end-2,end) = 1/4;
B = sparse(2*n-1,n-1);
for j = 1:n-1
    B(2*j,j) = 1;
    B(2*j-1,j) = 1/2;
    B(2*j+1,j) = 1/2;
end
L1 = kron(B,A);
L2 = kron(A,B);
LU = blkdiag(L1,L2);
B = sparse(n,n);
for i = 1:n
    B(2*i-1,i) = 1;
    B(2*i,i) = 1;
end
LP = kron(B,B);


