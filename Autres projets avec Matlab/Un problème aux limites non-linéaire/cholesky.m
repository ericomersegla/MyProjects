function L=cholesky(A)
[m,n]=size(A);
L=zeros(m,n);
if m~=n 
    disp('matrice n`est pas carrée')
    L=[];
    return
end
for i=1:m
    L(i,i)=sqrt(A(i,i)-sum(abs(L(i,1:i-1)).^2));
    if conj(L(i,i))~=L(i,i) %L(i,i) n'est pas réel
      disp('matrice n`est pas HDP')
      L=[];
        return
    end  
    if abs(L(i,i))<1e-16  %L(i,i) zéro
        disp('matrice n`est pas HDP')
        L=[];
        return
    end
    for j=i+1:m
            L(j,i)=conj((A(i,j)-sum(L(i,1:i-1).*conj(L(j,1:i-1))))/L(i,i));
    end
end