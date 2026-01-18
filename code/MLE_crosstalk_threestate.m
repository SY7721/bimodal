function [lam, kappa1, kappa2, q1, gamma, nu, AICc] = MLE_crosstalk_threestate(xdata, N)
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
tic
b4=[log(la1in(i))  log(la2in(i))  log(1/q1in-1)  log(gain(i))  log(vin(i))];
[bmin4, Smin4] = fminsearch(@S4fun_pm1,b4);
toc
la1_e = [la1_e exp(bmin4(1))]; 
la2_e=[la2_e exp(bmin4(2))]; 
q1_e=[q1_e 1/(1+exp(bmin4(3)))]; 
ga_e = [ga_e exp(bmin4(4))]; 
v_e = [v_e exp(bmin4(5))];
Smin = [Smin Smin4]; 
end

Smin_min = min(Smin); 
Smin_min_loction = find(Smin==Smin_min,1,'first'); 
lam = la1_e(Smin_min_loction);  
kappa1=la2_e(Smin_min_loction); 
kappa2=kappa1; 
q1=q1_e(Smin_min_loction); 
gamma = ga_e(Smin_min_loction); 
nu = v_e(Smin_min_loction);
k=5;
AICc=calculate_AICc(Smin_min, k, N);
end