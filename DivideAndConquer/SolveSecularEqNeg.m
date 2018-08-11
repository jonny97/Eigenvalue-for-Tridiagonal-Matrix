function eig=SolveSecularEqNeg(D,w,Dsort)
    %%
    % SolveSecularEq in Treften Page 231, (30.12) With negative case
    % 
    m=size(w,1);
    eig=zeros(m,1);
    for i=2:m
        eig(i)=(Dsort(i-1)+Dsort(i))/2;
        left=Dsort(i-1);
        right=Dsort(i);
        while true
            if (abs(left-right)<10e-10*(abs(left)+abs(right)))
                break
            else
                %% Bisection only
                %temp = 1-w'*((diag(D)-eye(m)*eig(i))\w); 
                %% Stable Newton:
                shifted_diag = diag(D)-eye(m)*eig(i);
                temp = 1 - w'*(shifted_diag\w); 
                dirivative = -w'*(shifted_diag\(shifted_diag\w)); % always positive! ( Given w!= vec 0)
                % The tangent line will be y = dirivative * (x - eig(i)) + temp
                intercept  = eig(i) - temp / dirivative;
 
                if abs(temp/dirivative)<10e-10*(abs(left)+abs(right))
                % Not much to improve in this case: Newton updates will be
                % ~10e-10                    
                    break
                
                elseif (intercept>0.98*left+0.02*right && intercept<right*0.98+0.02*left )
                    % This is a slight correction to avoid Newton goes to
                    % other curves or going to bad lambda that give a
                    % overflow/underflow of f(lambda).
                    eig(i) = intercept;                            


                elseif temp<0
                    right=eig(i);
                    eig(i)=(eig(i)+left)/2;
                else
                    left=eig(i);
                    eig(i)=(eig(i)+right)/2;              
                end
            end
        end
    end
    % Edge case:
    eig(1)=-2*abs(Dsort(1));
    left=-inf;right=Dsort(1);
    while true
        temp = 1-w'*((diag(D)-eye(m)*eig(1))\w);
        if (abs(left-right)<10e-10*(abs(left)+abs(right)))
            break
        elseif temp>0
            left=eig(1);
            eig(1)=(eig(1)+right)/2;
        else
            right=eig(1);
            if left==-inf
                eig(1)=-2*abs(right);              
            else
                eig(1)=(eig(1)+left)/2;   
            end
        end
    end
end