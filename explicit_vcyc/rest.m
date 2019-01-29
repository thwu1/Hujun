function R = rest(n)
% 2n->n
A1 = sparse(n-1,2*n-1);
A2 = sparse(n-1,2*n-1);
B = sparse(n,2*n);
for i = 1:n-1
    A1(i,2*i-1) = 1;
    A1(i,2*i+1) = 1;
    A2(i,2*i) = 1;
end
for i = 1:n
    B(i,2*i-1) = 1;
    B(i,2*i) = 1;
end
I1 = 0.25*kron(A2,B) + 0.125*kron(A1,B);
I2 = 0.25*kron(B,A2) + 0.125*kron(B,A1);
R = blkdiag(I1,I2);

