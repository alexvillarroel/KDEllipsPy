clear all
close all

models = load('models.dat','r');
misfits = load('misfits.dat','r');

n=length(models);

a = models((1:7:n),:);
b = models((2:7:n),:);
or = models((3:7:n),:);
dn = models((4:7:n),:);
ang = models((5:7:n),:);
maxslip = models((6:7:n),:);
vr = models((7:7:n),:);


plot(vr,maxslip,'*')
xlabel('Velocidad de Ruptura [km/s]')
ylabel('Slip maximo [m]')
xlim([1;4])
ylim([0;2])


mu = 62.5*10^9;

M0 = (62.5*10^9)*(pi.*(a*1000/2).*(b*1000/2)).*(0.63212055*maxslip);

figure

plot(vr,M0,'*')
xlabel('Velocidad de Ruptura [km/s]')
ylabel('Momento sismico [N]')
