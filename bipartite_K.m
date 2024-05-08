%----------------------
%Bipartite network
%----------------------
for i=1:K
    for j=1:K
        A(i,j) = 0;     
    end
    
    for j=K+1:n
        A(i,j) = 1;    
    end
        
end  

for i=K+1:n
    for j=1:K
        A(i,j) = 1;    
    end

    for j=K+1:n
        A(i,j) = 0;  
    end   
end   

for i=1:n
        A(i,i) = 1;  %Agents listen to themselves
end