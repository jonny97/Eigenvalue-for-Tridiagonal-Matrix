m = 512; %4,8,16,32,64,128,256,512
Error=0;
for i=1:10
    A = randn(m);
    B = hess(A'*A)*10^40; % huge numbers
    %B = hess(A'*A);      % normal numbers

    E = eig(B);
    %ANS = bisection(B,0,10000,10e-14,-1,-1);% normal numbers
    ANS = bisection(B,0,10000*10^40,10e-14,-1,-1);% huge numbers
    if size(ANS,1)~=m
        fprintf('BAD!\n');
    end
    Error= Error + norm(E(1:size(ANS,1))-ANS);
end
Error=Error/10