function [U,P] = uzawa(A,B,U,P,F,v,a)
n = sqrt(length(P));
for i = 1:v
U = A\(F - B*P);
P = P + a*(B'*U);
end