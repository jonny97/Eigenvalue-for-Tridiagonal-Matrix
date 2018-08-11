function eig=QuickBisection(A,a,b,tol)
    %%
    % Using Sturm sequences and Bisection to solve for all eigenvalues of
    % A in [a,b].
    % This is a naive version which has no detection of bad edge cases.
    % precondition: the eigenvalues of A must be unique in context of
    % tolerence. A is tridiagonal and PSD
    
    %%
    [num_b,~] = Sturm(A,b,tol);
    [num_a,~]           = Sturm(A,a,tol);
    eig                 = zeros(num_b-num_a,1);
    
    % if there is no eigvalue in [a,b), no need to ocontinue
    if (num_b>num_a)
        if (b-a<tol && num_b-num_a==1)
            eig(1)=(b+a)/2;
        else
            middle = (b+a)/2;
            num_middle = Sturm(A,middle,tol);
            if (num_middle-num_a>=1)
                eig(1:num_middle-num_a) = QuickBisection(A,a,middle,tol);
            end
            if (num_b>=num_middle+1)
                eig(num_middle-num_a+1:num_b-num_a) = QuickBisection(A,middle,b,tol);
            end
        end
    end
end