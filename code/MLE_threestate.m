function [lam1, lam2, gamma, nu, AICc] = MLE_threestate(xdata, N)
global  xdata tdata tt
S=length(xdata); 
tt=[]; tdata=[];
tt=N.*xdata; tdata=0:1:S-1;

for i=1:10
    lain(i)=0.1+4.9*rand();
    gain(i)=0.1+4.9*rand();
    vin(i)=5+25*rand();
end

la1_e = []; la2_e=[]; ga_e = []; v_e = []; Smin = []; 
for i = 1:10
la1in(i)= 2*lain(i); la2in(i)=2*lain(i);
tic
 b3=[log(la1in(i))  log(la2in(i))  log(gain(i))  log(vin(i))];
[bmin3, Smin3] = fminsearch(@S3fun_pm1,b3); 
toc
la1_e = [la1_e exp(bmin3(1))];   
la2_e = [la2_e exp(bmin3(2))]; 
ga_e = [ga_e exp(bmin3(3))]; 
v_e = [v_e exp(bmin3(4))];
Smin = [Smin Smin3]; 
end

Smin_min = min(Smin3);
Smin_min_loction3 = find(Smin3==Smin_min,1,'first'); 
lam1 = la1_e(Smin_min_loction3);  
lam2 = la2_e(Smin_min_loction3);   
gamma = ga_e(Smin_min_loction3); 
nu = v_e(Smin_min_loction3);
k=4;
AICc=calculate_AICc(Smin_min, k, N);
end