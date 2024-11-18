%-----------------------------------------------------------------------
%Star network with two stars and a wheel-like structure
%-----------------------------------------------------------------------
A = zeros(n,n);

for i=1:n
    
    %if i>1
        A(i,1) = 1;
    %end  
    %if i<n
        A(i,n) = 1;
    %end  
    if i==n
        A(i,1) = 1;
    end
    if i==1
        A(i,n) = 1;
    end
    
A(i,i) = 1;  %Agents listen to themselves

end

