function Rg = rest_G(n)
% 2n->n
A = sparse(n,2*n);
for i = 1:n
    A(i,2*i-1) = 1/2;
    A(i,2*i) = 1/2;
end
Rg = kron(A,A);