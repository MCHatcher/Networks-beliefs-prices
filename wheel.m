%----------------
%Wheel network
%----------------

A = eye(n);  %Agents listen to themselves

for i=1:n
    for j=1:n
    
        if i>1
            A(i,i-1) = 1;
        end
    
        if i<n
            A(i,i+1) = 1;
        end
    
        if i==n
            A(i,1) = 1;
        end

        if i==1
            A(i,n) = 1;
        end
    
    end 
end
