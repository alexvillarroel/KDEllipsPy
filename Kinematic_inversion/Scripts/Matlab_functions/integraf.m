function [ix]= integraf(d,Fs);
%
% [ix]= integraf(d,Fs)
%       Integra la funcion d muestrada a una frecuencia de Fs
%       utilizando el  metodo del filtrado con coeficientes de
%       trapecio
%

[m,n]=size(d);
d=columna(d);

xn=length(d);
dt=1/Fs;
b=[1 1]/2;
a=[1 -1];

ix=filter(b,a,d)*dt;
ix=ix-ix(1).*ones(xn,1);


