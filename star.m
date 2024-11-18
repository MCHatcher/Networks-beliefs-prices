%---------------
%Star network
%---------------

A = zeros(n,n);

for i=1:n
    A(i,1) = 1;
    A(1,i) = 1;

    A(i,i) = 1;   %Agents listen to themselves
end