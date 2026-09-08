function parameter=MLE_p(model)

    %%% initial values
    for i=1:10
        lain(i)=0.1+4.9*rand();
        gain(i)=0.1+4.9*rand();
        vin(i)=5+45*rand();
        betain(i)=rand;   
    end
    
    if model==1  %%% Fitted by beta+telegraph
       la_e = []; ga_e = []; v_e = [];  Smin = []; beta_e=[];
       for i=1:10  
          
           b1=[log(lain(i))  log(gain(i))  log(vin(i))   log(1/betain(i)-1)];
           [bmin1, Smin1] = fminsearch(@S1fun_BM_pm,b1); % 调用搜索函数
          
           la_e = [la_e exp(bmin1(1))];
           ga_e = [ga_e exp(bmin1(2))];
           v_e = [v_e exp(bmin1(3))];
           beta_e=[beta_e 1/(1+exp(bmin1(4)))];
           Smin = [Smin Smin1];
       end
       Smin_min = min(Smin);
       Smin_min_loction = find(Smin==Smin_min,1,'first');
       latr = la_e(Smin_min_loction);
       gatr = ga_e(Smin_min_loction);
       beta = beta_e(Smin_min_loction);
       vtr = v_e(Smin_min_loction);
       parameter=[latr gatr vtr beta  Smin_min];

    elseif model==2
        lar1_e = []; lar2_e=[]; qr1_e=[]; gar_e = []; vr_e = []; Smin = []; beta_e1=[];
        for i=1:10
            la1in(i)= lain(i) ;  la2in(i)= lain(i)  ; q1in=0.5 ;  
          
            b2=[log(la1in(i))  log(la2in(i))  log(1/q1in-1)  log(gain(i))  log(vin(i)) log(1/betain(i)-1)];
            [bmin2, Smin2] = fminsearch(@S2fun_BM_pm,b2);
           
            lar1_e = [lar1_e exp(bmin2(1))]; 
            lar2_e=[lar2_e exp(bmin2(2))]; 
            qr1_e=[qr1_e 1/(1+exp(bmin2(3)))];
            gar_e = [gar_e exp(bmin2(4))];
            vr_e = [vr_e exp(bmin2(5))];
            beta_e1 = [beta_e1 1/(1+exp(bmin2(6)))];
            Smin = [Smin Smin2];
         end
         Smin_min1 = min(Smin);
         Smin_min_loction2 = find(Smin==Smin_min1,1,'first');
         la1 = lar1_e(Smin_min_loction2);
         la2=lar2_e(Smin_min_loction2);
         q1=qr1_e(Smin_min_loction2);
         ga = gar_e(Smin_min_loction2);
         v = vr_e(Smin_min_loction2);
         beta=beta_e1(Smin_min_loction2);
         parameter=[la1 la2 q1 ga v beta Smin_min1];

    elseif model==3  %%% Fitted by p+three-state
        la1_e = []; la2_e=[]; ga_e = []; v_e = []; Smin = []; beta_e=[];
 
        for i = 1:10
            la1in(i)= 2*lain(i); la2in(i)=2*lain(i);
          
            b3=[log(la1in(i))  log(la2in(i))  log(gain(i))  log(vin(i)) log(1/betain(i)-1)];
            [bmin3, Smin3] = fminsearch(@S3fun_BM_pm,b3); % 调用搜索函数

            la1_e = [la1_e exp(bmin3(1))];  
            la2_e = [la2_e exp(bmin3(2))];
            ga_e = [ga_e exp(bmin3(3))];
            v_e = [v_e exp(bmin3(4))];
            beta_e=[beta_e 1/(1+exp(bmin3(5)))];
            Smin = [Smin Smin3];
        end
        Smin_min = min(Smin);
        Smin_min_loction = find(Smin==Smin_min,1,'first');
        latr1 = la1_e(Smin_min_loction);
        latr2 = la2_e(Smin_min_loction);
        gatr = ga_e(Smin_min_loction);
        beta = beta_e(Smin_min_loction);
        vtr = v_e(Smin_min_loction);
        parameter=[latr1 latr2 gatr vtr beta  Smin_min];

    elseif model==4  %%Fitted by p+three-state cross-talk
         la_e = []; kappa_e=[]; q1_e=[]; ga_e = []; v_e = []; Smin = [];  beta_e=[];

         for i=1:10
             la1in(i)= lain(i) ;  la2in(i)= lain(i); q1in=0.5 ;   

             b4=[log(la1in(i))  log(la2in(i))  log(1/q1in-1)  log(gain(i))  log(vin(i))  log(1/betain(i)-1)];
             [bmin4, Smin4] = fminsearch(@S4fun_BM_pm,b4);

             la_e = [la_e exp(bmin4(1))]; 
             kappa_e=[kappa_e exp(bmin4(2))]; 
             q1_e=[q1_e 1/(1+exp(bmin4(3)))];
             ga_e = [ga_e exp(bmin4(4))];
             v_e = [v_e exp(bmin4(5))];
             beta_e=[beta_e 1/(1+exp(bmin4(6)))];
             Smin = [Smin Smin4];
         end
         Smin_min1 = min(Smin);
         Smin_min_loction2 = find(Smin==Smin_min1,1,'first');
         la1 = la_e(Smin_min_loction2);
         kappa1=kappa_e(Smin_min_loction2);
         q1=q1_e(Smin_min_loction2);
         ga = ga_e(Smin_min_loction2);
         v = v_e(Smin_min_loction2);
         beta=beta_e(Smin_min_loction2);
         parameter=[la1 kappa1  kappa1 q1 ga v beta Smin_min1];

    elseif model==5
        la_e = []; mu_e=[];   ga_e = [];   v_e = [];   Smin = []; beta_e=[];
        for i = 1:10
            muin = 0.1+4.9*rand();

            b5=[log(lain(i))  log(gain(i))  log(vin(i)) log(muin) log(1/betain(i)-1)];
            [bmin5, Smin5] = fminsearch(@S5fun_BM_pm,b5);

            la_e = [la_e exp(bmin5(1))]; 
            ga_e = [ga_e exp(bmin5(2))];
            v_e = [v_e exp(bmin5(3))];
            mu_e=[mu_e exp(bmin5(4))];
            beta_e=[beta_e 1/(1+exp(bmin5(5)))];
            Smin = [Smin Smin5];
        end
        Smin_min1 = min(Smin);
        Smin_min_loction = find(Smin==Smin_min1,1,'first');
        la = la_e(Smin_min_loction);
        ga = ga_e(Smin_min_loction);
        v = v_e(Smin_min_loction);
        mu=mu_e(Smin_min_loction);
        beta=beta_e(Smin_min_loction);
        parameter=[la ga v mu beta Smin_min1];

    elseif model==6
        la_e = []; nu_e=[];  ga_e = [];  v_e = []; Smin = []; beta_e=[];
 
        for i = 1:10
            nuin=0.1+4.9*rand();

            b6=[log(lain(i))  log(gain(i))  log(vin(i)) log(nuin) log(1/betain(i)-1)];
            [bmin6, Smin6] = fminsearch(@S6fun_BM_pm,b6); % 调用搜索函数

             la_e = [la_e exp(bmin6(1))];  
             ga_e = [ga_e exp(bmin6(2))];
             v_e = [v_e exp(bmin6(3))];
             nu_e = [nu_e exp(bmin6(4))];
             beta_e=[beta_e 1/(1+exp(bmin6(5)))];
             Smin = [Smin Smin6];
         end

             Smin_min1 = min(Smin);
             Smin_min_loction = find(Smin == Smin_min1,1,'first');
             lar = la_e(Smin_min_loction);    
             ga = ga_e(Smin_min_loction);
             v = v_e(Smin_min_loction);
             nu = nu_e(Smin_min_loction);
             beta=beta_e(Smin_min_loction);
             parameter=[lar ga v nu beta Smin_min1];

    end
end