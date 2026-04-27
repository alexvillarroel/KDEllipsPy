clear all
close all
load models.dat
load misfits.dat

misfits_cte=0.5;
tt_cte=30;
tt_cte2=0.5;
aa=0;
ll=0;
chi=100000;
% 
% fin=length(misfits)/3;
% for i=1:fin
%     j=(i)*3;
%    if and((models(j)>5),misfits(i)<0.55)
%        aa=models(j);
%        ll=i;
%        chi=misfits(i);  
%        dm=models(j-1);
%        vr=models(j);
%    end
%    
% end
% chi
% aa
% ll
% dm
% vr
%  misfits=misfits(1:100);
% %  misfits(20000)=0;
%  models=models(1:300);
 for i=1:length(models)/3
    models_vr(i)=models((i-1)*3+3);
    models_dm(i)=models((i-1)*3+2);
%     models_alfa(i)=models((i-1)*7+3);
    models_a(i)=1000*models((i-1)*3+1);
%     models_b(i)=1000*models((i-1)*7+2);
%     models_c(i)=pi*models_a(i)*models_b(i);
    tt(i)=tt_cte;
        if misfits(i)>misfits_cte
        misfits(i)=misfits_cte;
        tt(i)=tt_cte2;
        end
end
% misfits(length(misfits))=0;
mu=62.5*10^9;

for k=1:length(models)/3
aa=0; 
ccc=0;
    for i=1:round(models_a(k)/200)
        for j=1:round(0.53*models_a(k)/200)
              ll=(i-1)*round(0.53*models_a(k)/200)+j;
              ii=(i-1)*500;
              jj=(j-1)*500;
              ddd=(ii^2/models_a(k)^2)+(jj^2/(0.53*models_a(k))^2);
              dd=exp(-ddd);
              if ddd<1
              aa=dd*500*500+aa;
              end
        end
    end

    Mo(k)=4*mu*aa*models_dm(k);
    if Mo(k)>3E19
        Mo(k)=3E19;
    end
    
end
% Mo=(2/3)*(log10(Mo)-9.1);



% plot(models_a,models_b,'*')

% figure
% scatter(models_a,models_b,tt,models_alfa,'filled');
% colorbar
% misfits(length(misfits))=0;
figure
scatter(models_vr,models_dm,tt,misfits,'filled');

colorbar
xlabel('Vr')
ylabel('Dm')


% figure
% 
% scatter(models_dm.^3,Mo,tt,misfits,'filled');
% 
% colorbar
% xlabel('Dm3')
% ylabel('Momento')

figure
scatter3(models_dm,models_vr,Mo,tt,misfits,'filled');

colorbar
xlabel('Dm')
ylabel('vr')
zlabel('Momento')

figure

scatter(models_vr,Mo,tt,misfits,'filled');

colorbar
xlabel('Vr')
ylabel('Momento')




% figure
% scatter(((models_a+models_b)/2).^3,models_dm,tt,misfits,'filled');
% 
% colorbar
% xlabel('L caracteristico')
% ylabel('Dm')

figure
scatter(((models_vr)),models_a/1000,tt,misfits,'filled');

colorbar
xlabel('Vr')
ylabel('Largo a [km]')

% figure
% scatter(((models_b)),models_dm,tt,misfits,'filled');
% 
% colorbar
% xlabel('L b')
% ylabel('Dm')
