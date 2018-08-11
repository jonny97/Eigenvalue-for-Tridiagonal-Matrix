function [Q,eigvalue]=DivideAndConquer(A)
    %%
    % Using the idea of divide and conquer to solve all eigenvalues of A.
    % Return: All eigenvalues of A; eigvector matrix Q
    % precondition: A is m by m (m>=2), symmetric, tridiagonal.
    % 1. If A is not irreducable, reduce A
    % 2. Divide: find eigenvalues/eigenvector of upperhalf/ lowerhalf
    % 3. MergeSort eigenvalues to be eig_sorted
    % 4. Conquer: Find eigenvalue of D +/- w*w'
    % 5. Conquer: Given eigenvalue, find eigenvector in O(m^2)
    
    m=size(A,1);
    if m<=2
        [Q,D]   = eig(A);
        eigvalue= diag(D);
        return
    end
    
    % Ruduce
    for i=1:m-1
        if (A(i,i+1)==0)
            [Q1,eig_top] = DivideAndConquer(A(1:i,1:i));
            [Q2,eig_bot] = DivideAndConquer(A(i+1:m,i+1:m));
            eigvalue     = [eig_top;eig_bot];
            Q            = [Q1 zeros(i,m-i);zeros(m-i,i) Q2];      
            return;
        end
    end
    
    % Divide
    n       = floor(m/2);
    beta    = A(n+1,n);
    topmat  = A(1:n,1:n);topmat(end,end)=topmat(end,end)-beta;
    botmat  = A(n+1:m,n+1:m);botmat(1,1)=botmat(1,1)-beta;
    [Q1,eig_top] = DivideAndConquer(topmat);
    [Q2,eig_bot] = DivideAndConquer(botmat);
    q1           = Q1(end,:)';
    q2           = Q2(1,:)';
    
    %MergeSort eigenvalues
    i=0;j=0;
    eig_sorted=zeros(m,1);
    while i+j<m
       if (j>=m-n) 
           eig_sorted(i+j+1)=eig_top(i+1);
           i=i+1;
       elseif (i>=n)
           eig_sorted(i+j+1)=eig_bot(j+1);
           j=j+1; 
       elseif eig_top(i+1)<eig_bot(j+1)
           eig_sorted(i+j+1)=eig_top(i+1);
           i=i+1;
       else
           eig_sorted(i+j+1)=eig_bot(j+1);
           j=j+1; 
       end
    end
    
    %Conquer
    if beta>0
        w            = sqrt(beta)*[q1;q2];
        eigvalue     = SolveSecularEq([eig_top;eig_bot],w,eig_sorted);
    else
        w            = sqrt(-beta)*[q1;q2];
        eigvalue     = SolveSecularEqNeg([eig_top;eig_bot],w,eig_sorted);
    end
    Q   = zeros(m,m);
    for i=1:m
        Q(:,i)  = (diag([eig_top;eig_bot])-eigvalue(i)*eye(m))\w;
        Q(:,i)  = Q(:,i)/norm(Q(:,i));
    end
    Q   = [Q1 zeros(n,m-n);zeros(m-n,n) Q2]*Q;
end