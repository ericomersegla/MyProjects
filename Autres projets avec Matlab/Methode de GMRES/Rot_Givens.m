function [H,c_given,s_given] = Rot_Givens(H,c,s,n)

   for j = 1:n-1
     H_sortie=c(j)*H(j)+s(j)*H(j+1);
     H(j+1)=c(j)*H(j+1)-s(j)*H(j);
     H(j)=H_sortie;
            
    end
    [c_given,s_given] = rotation(H(n), H(n+1));

    H(n)=c_given*H(n)+s_given*H(n+1);
    H(n+1)=0.0;
end


function [c,s] = rotation(a,b)

      r = sqrt(a^2 + b^2);
      c = a / r;
      s =b / r;
end
