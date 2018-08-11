function [num,is_eigval]=Sturm(A,lambda,tol)
    % Returns the number of eigenvalues of A between (-inf,lambda)
    % Using standard Sturm method
    % fprintf('check sturm %0.30f\n',lambda);
    is_eigval=0;
    m = size(A,1);
    %A= A-lambda*speye(m,m);
    %for i=1:m
    %    A(i,i)=A(i,i)-lambda;
    %end
    
    
    det = zeros(m,1);
    det(1) = (A(1,1)-lambda);
    det(2) = (A(1,1)-lambda)*(A(2,2)-lambda)-A(1,2)*A(2,1);
    for i = 3:m
        det(i) = (A(i,i)-lambda)* det(i-1) - A(i,i-1) * A(i-1,i) * det(i-2);
        % shift numbers so that the calculation will not below up
        % Note that, my code does not implement cases for underflow,which
        % can be really tricky as the machine cannot tell if it is a zero
        % or it is a really small number. The best way is to scale the
        % input matrix if it has really small eigenvalues.
        if (abs(det(i-1)/det(1))>10^20 )
            det(i)=det(i)/abs(det(i-1));
            det(i-1)=det(i-1)/abs(det(i-1));
        end
    end

    num = 0;
    if (det(1)<0)
        num=1;
    end

    for i = 2:m
        if (det(i) < 0 && det(i-1)>=0)
            num=num+1;
        end
        if (det(i) > 0 && det(i-1)<=0)
            num=num+1;
        end
    end
    
    % If the input lambda is close enough to the eigenvalue, we will 
    % shift the input by little to aviod dealing with boundary cases
    % Thus, this function finds the number of eigenvalues in (-inf,lambda)
    % Note: This is unlikely to happen
    if (abs(det(end))<tol)
        num = Sturm(A, lambda*(1-tol)-tol,tol);
        is_eigval=1;
    end
end