function [lam, gamma, nu, AICc] = MLE_telegraph(xdata, N)
global  xdata tdata tt
S=length(xdata); 
tt=[]; tdata=[];
tt=N.*xdata; tdata=0:1:S-1;
for i=1:10
    lain(i)=0.1+4.9*rand();
    gain(i)=0.1+4.9*rand();
    vin(i)=5+25*rand();
end

la_e = []; ga_e = []; v_e = []; laga_e = []; Smin = []; 
for i=1:10
tic
b=[log(lain(i))  log(gain(i))  log(vin(i))];
[bmin1, Smin1] = fminsearch(@S1fun_pm1,b); 
toc
la_e = [la_e exp(bmin1(1))];
ga_e = [ga_e exp(bmin1(2))];
v_e = [v_e exp(bmin1(3))];
Smin = [Smin Smin1];
end

Smin_min = min(Smin); 
Smin_min_loction1 = find(Smin==Smin_min,1,'first');
lam = la_e(Smin_min_loction1); 
gamma = ga_e(Smin_min_loction1); 
nu = v_e(Smin_min_loction1);
k=3;
AICc=calculate_AICc(Smin_min, k, N);
end