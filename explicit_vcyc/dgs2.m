function [U,P] = dgs2(A,B,U,P,F,v)
for i = 1:v
G = B'*U;
U = U + tril(A)\(F - A*U - B*P);
Ap = B'*B;
dq = tril(Ap)\(G - B'*U);
U = U + B*dq;
P = P - B'*(B*dq);
end