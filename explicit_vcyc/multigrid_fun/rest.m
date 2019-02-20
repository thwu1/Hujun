function R = rest(n)
% 2n->n
H = sparse(n-1,2*n-1);
B = sparse(n,2*n);
for i = 1:n-1
    H(i,2*i-1) = 1/8;
    H(i,2*i+1) = 1/8;
    H(i,2*i) = 1/4;
end
for i = 1:n
    B(i,2*i-1) = 1;
    B(i,2*i) = 1;
end
I1 = kron(H,B);
I2 = kron(B,H);
R = blkdiag(I1,I2);

