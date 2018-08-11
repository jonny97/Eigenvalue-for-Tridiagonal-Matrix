function eig=bisection(A,a,b,tol,num_a,num_b)
    %%
    % Using Sturm sequences and Bisection to solve for all eigenvalues of
    % A in [a,b]
    % precondition: the eigenvalues of A must be unique in context of
    % tolerance. A is tridiagonal and PSD. Therefore, all eigenvalues of A
    % are distinct.
    
    %%
    %fprintf('calling bisect %d %d\n',a,b);
    b_is_eigval=0;
    if (num_a==-1 && num_b==-1)
        [num_b,b_is_eigval] = Sturm(A,b,tol);
        [num_a,~]           = Sturm(A,a,tol);
    end
    eig                 = zeros(num_b-num_a,1);
    % deal with edge case: b being eigenvalue
    if (b_is_eigval==1)
        eig(end+1)=b;
    end
    
    % if there is no eigvalue in [a,b), no need to ocontinue
    if (num_b>num_a)
        if ((b-a)/a<tol && num_b-num_a==1)
            eig(1)=(b+a)/2;
        else
            middle = (b+a)/2;

            num_middle = Sturm(A,middle,tol);
            %pause;
            if (num_middle-num_a>=1)
                % temp =bisection(A,a,middle,tol);
                % disp(num_middle-num_a)
                % disp(size(temp));
                % eig(1:num_middle-num_a) = temp;
                eig(1:num_middle-num_a) = bisection(A,a,middle,tol,num_a,num_middle);
            end
            if (num_b>=num_middle+1)
                eig(num_middle-num_a+1:num_b-num_a) = bisection(A,middle,b,tol,num_middle,num_b);
            end
        end
    end
    %disp(eig);
    %fprintf('ending bisect %d %d\n',a,b);
end