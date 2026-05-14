function t=time(a,Fs)
% t=time(a,Fs)
% genera vector de tiempo.
% t : vector de tiempo de longitud length(a)
%     Si a es escalar genera vector de longitud a
% Fs: frecuencia de muestreo


if length(a)==1
   t=(0:a-1)/Fs;
else
   t=(0:(length(a)-1))/Fs;
end
t=t';
