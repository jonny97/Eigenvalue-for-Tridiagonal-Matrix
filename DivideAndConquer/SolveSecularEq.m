function eig=SolveSecularEq(D,w,Dsort)
    %%
    % SolveSecularEq in Treften Page 231, (30.12)
    % 
    m=size(w,1);
    eig=zeros(m,1);
    for i=1:m-1
        eig(i)=(Dsort(i)+Dsort(i+1))/2;
        left=Dsort(i);
        right=Dsort(i+1);
        while true
            if (abs(left-right)<10e-10*(abs(left)+abs(right)))
            % Given the interval is so small, no need to dig more
                break        
            else
                %% Bisection only
                %temp = 1+w'*((diag(D)-eye(m)*eig(i))\w); 
                %% Stable Newton:
                shifted_diag = diag(D)-eye(m)*eig(i);
                temp = 1+w'*(shifted_diag\w); 
                dirivative = w'*(shifted_diag\(shifted_diag\w)); % always positive! ( Given w!= vec 0)
                % The tangent line will be y = dirivative * (x - eig(i)) + temp
                intercept  = eig(i) - temp / dirivative;
                
                if abs(temp/dirivative)<10e-10*(abs(left)+abs(right))
                % Not much to improve in this case: Newton updates will be
                % ~10e-10
                    break
                
                elseif (intercept>0.99*left+0.01*right && intercept<right*0.99+0.01*left )
                    % This is a slight correction to avoid Newton goes to
                    % other curves or going to bad lambda that give a
                    % overflow/underflow of f(lambda).
                    eig(i) = intercept;
                 
                elseif temp>0
                    % Bisection
                    right=eig(i);
                    eig(i)=(eig(i)+left)/2;
                else
                    % Bisection
                    left=eig(i);
                    eig(i)=(eig(i)+right)/2;              
                end
            end
        end
        %fprintf('end\n');

    end
    % Edge case:    
    eig(m)=2*Dsort(m);
    left=Dsort(m);right=inf;
    while true
        temp = 1+w'*((diag(D)-eye(m)*eig(m))\w);
        if (abs(left-right)<10e-10*(abs(left)+abs(right)))
            break
        elseif temp>0
            right=eig(m);
            eig(m)=(eig(m)+left)/2;
        else
            left=eig(m);
            if right==inf
                eig(m)=2*left;              
            else
                eig(m)=(eig(m)+right)/2;   
            end
        end
    end  
end