m = 512;
Error=0;
for i=1:10
    A = randn(m);
    B = hess(A'*A)*10^20;
    [V1,D1] = eig(B);
    [V2,D2] = DivideAndConquer(B);
    Error = (norm(D1-diag(D2)));
end
Error=Error/10
