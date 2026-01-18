function [lam1, lam2, q1, gamma, nu, AICc] = MLE_crosstalk(xdata, N)
global  xdata tdata tt
S=length(xdata); 
tt=[]; tdata=[];
tt=N.*xdata; tdata=0:1:S-1;
for i=1:10
    lain(i)=0.1+4.9*rand();
    gain(i)=0.1+4.9*rand();
    vin(i)=5+25*rand();
end
la1_e = []; la2_e=[]; q1_e=[]; ga_e = []; v_e = []; Smin = [];
for i=1:10
 la1in(i)= lain(i) ;  la2in(i)= lain(i)  ; q1in=0.5 ;   
b2=[log(la1in(i))  log(la2in(i))  log(1/q1in-1)  log(gain(i))  log(vin(i))];
[bmin2, Smin4] = fminsearch(@S2fun_pm1,b2);
la1_e = [la1_e exp(bmin2(1))]; 
la2_e=[la2_e exp(bmin2(2))]; 
q1_e=[q1_e 1/(1+exp(bmin2(3)))]; 
ga_e = [ga_e exp(bmin2(4))]; 
v_e = [v_e exp(bmin2(5))];
Smin = [Smin Smin4];
end

Smin_min = min(Smin);
Smin_min_loction2 = find(Smin==Smin_min,1,'first'); 
lam1 = la1_e(Smin_min_loction2); 
lam2=la2_e(Smin_min_loction2); 
q1=q1_e(Smin_min_loction2); 
gamma = ga_e(Smin_min_loction2); 
nu = v_e(Smin_min_loction2); 
k=4;
AICc=calculate_AICc(Smin_min, k, N);
end