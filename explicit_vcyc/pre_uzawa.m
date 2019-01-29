function [U,P] = pre_uzawa(A,B,U,P,F,v,a)
for i = 1:v
    [U,P] = vcyc(U,P,F);
    P = P + a*(B'*U);
end