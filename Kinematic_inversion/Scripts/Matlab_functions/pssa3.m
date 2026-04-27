function [ac,v,d]=pssa3(ac,Fs,N,f1,f2,nophase,izero,graf,titulo);
%  [ac,v,d]=pssa(acc,Fs,N,f1,f2,filfilt(0 o 1),zeros(0/1),graf(0/1))
% 
% ac, v, d    Aceleracion, velocidad y desplazamiento corregidos
% Fs          Frecuencia de muestreo de acc
% N = 4        Orden del filtro butterworth
% f1,f2       Frecuencia de corte del filtro
% f1 = 0.167
% f2 = 25
% nophase (0 causal, 1 acausal)    Procesa en ambas direcciones todas las se�ales (se
%             duplica orden del filtro. Si 1 se realiza filtfilt 
% zeros       Agreaga zeros (1) al inicio y final del registro de acc
% graf        0 no grafica 1 graf , [def 1]
%                                            

% figure
% fi = fopen('desp.dat', 'w')
if nargin <8
    graf=1;
end

%aceleracion original
acc=ac;

[m,n]=size(ac);
if m < n
    mm=m;
    m=n;
    n=mm;
end

%% dcorrea 00 de correcion
if     ((f1=='NO')&(f2=='NO'))
    filtrar=0;
else
    filtrar=1;
end;

% --------------------------------------------------------------------------------
% PROCESO

if filtrar 
    
    % ---------------------------------------------------------
    % LA SE�AL ES FILTRADA DE ACUERDO A LOS VALORES DE F1 Y F2
    % ---------------------------------------------------------
    
    if     ((f1=='NO')&(f2~='NO')) % FILTRO PASA-BAJOS
        w=f2/(Fs/2);
        [b,a]=butter(N,w);
    elseif ((f1~='NO')&(f2=='NO')) % FILTRO PASA-ALTOS
        w=f1/(Fs/2);
        [b,a]=butter(N,w,'high');
    elseif ((f1~='NO')&(f2~='NO')) % FILTRO PASA-BANDA
        w=[ f1 f2 ]/(Fs/2);
        [b,a]=butter(N,w);
    end;
    
    %%%%%%%%%% dcorrea 00
    
    %w=[f1 f2]/Fs*2;
    %[b,a]=butter(N,w);
    
%     
%     if graf
%         b;
%         a;
%     end
     nfi=length(b)*10;
    
    if izero==1   
        ac=[zeros(1,nfi) detrend(ac)' zeros(1,nfi)]';
    else 
        ac=detrend(ac);
    end
    
    if nophase
        ac=filtfilt(b,a,ac);
    else
        ac=filter(b,a,detrend(ac));
    end
    v=integraf(ac,Fs);
    v=detrend(v);
   if nophase
       v=filtfilt(b,a,detrend(v));
   else
       v=filter(b,a,detrend(v));
   end
    d=integraf(v,Fs);
  if nophase
      d=filtfilt(b,a,detrend(d));
  else
      d=filter(b,a,detrend(d));
%       fprintf(fi,'%f\n',d)
%       fclose(fi)
  end
%  elimina zeros que se agregaron
    if izero
        ac=ac(nfi+1:nfi+m);
        v=v(nfi+1:nfi+m);
        d=d(nfi+1:nfi+m);
%         fprintf(fi,'%f\n',d)
%         fclose(fi)
    end
    
else
    
    % ---------------------------------------------------------
    % LA SE�AL NO SER� FILTRADA
    % ---------------------------------------------------------
    
    ac=detrend(ac);
    v=integraf(ac,Fs);
    d=integraf(v,Fs);
    fprintf(fi,'%f\n',d)
    fclose(fi)
 end;  
%ac=ac/980;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%
% OJO
%%%%%%%%%%%%%
%v=v.^2; d=d.^2; ac=ac.^2;
%v=v/max(abs(v)); d=d/max(abs(d)); ac=ac/max(abs(ac));



t=time(ac,Fs);

if graf
%    hold off
%    clf
    
    iflag=0;
%    if iflag==1
%        subplot(211), plot(t,ac,t,ac), grid, ylabel('Aceleracion [g]'), xlabel('Tiempo (seg)')
%        subplot(212), plot(t,v), grid, ylabel('Velocidad [cm/seg]'), xlabel('Tiempo (seg)'),
%        pause;clf;
%        subplot(211), plot(t,d), grid, ylabel('Desplazamiento [cm]'), xlabel('Tiempo (seg)')
%    end
    
    if iflag==0
        subplot(311), plot(t,ac,'k'), grid, ylabel('Aceleracion [cm/s2]'), xlabel('Tiempo (seg)'),title(titulo)
        subplot(312), plot(t,v,'k'), grid, ylabel('Velocidad [cm/seg]'), xlabel('Tiempo (seg)'),title('Velocidad')
        subplot(313),plot(t,d,'k'), grid, ylabel('Desplazamiento [cm]'), xlabel('Tiempo (seg)') ,title('Desplazamiento')
%        subplot(514), plot(at,d,t,v), grid, ylabel('Velocidad [m/seg]'), xlabel('Tiempo (seg)') %,title('Velocidad')
%        subplot(515), plot(t,d,t,ac), grid, ylabel('Desplazamiento [m]'), xlabel('Tiempo (seg)') %,title('Desplazamiento')
    end
end

% d=v;
% Max_A(acc/980)
% Max_A(v)
% Max_A(d)