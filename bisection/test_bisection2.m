m = 64; %4,8,16,32,64,128,256,512
Error=0;
a = 1;
b = 512; %2, 4,8,16,32,64,128,256
Eigenvalues_found =0;
for i=1:10
    A = randn(m);
    B = hess(A'*A);      % normal numbers
    E = eig(B);
    ANS = bisection(B,a,b,10e-14,-1,-1);% normal numbers
    
    E_in_range=zeros(0,1);
    for j=1:m
        if (E(j)>=a && E(j)<=b)
            E_in_range(end+1,1)=E(j);
        end
    end
    
    if size(ANS,1)~=size(E_in_range,1)
        fprintf('BAD! %d %d \n',size(ANS,1),size(E_in_range,1));
        ANS
        E_in_range
    end
    Eigenvalues_found = Eigenvalues_found + size(ANS,1);
    Error= Error + norm(E_in_range-ANS);
end
Error=Error/10
Eigenvalues_found