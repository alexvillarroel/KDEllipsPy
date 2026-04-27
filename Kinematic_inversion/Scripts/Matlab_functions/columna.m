function c=columna(a)
%      function c=columna(a)
% coloca el vector como columna la direccion mayor.
% 6-11-01 rbk modifica para numero complejo no da conjugado

[m,n]=size(a);

if m < n
c=a.';
else
c=a;
end


return
