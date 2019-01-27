function [U,P] = uzawa(A,B,U,P,F,v,a)
n = sqrt(length(P));
for i = 1:v
U = GS(A,B,U,P,F);
P = P + a*(B'*U);
end