function U = GS4(A,U,F,v)
for i = 1:v
    U = U + tril(A)\(F - A*U);
end
