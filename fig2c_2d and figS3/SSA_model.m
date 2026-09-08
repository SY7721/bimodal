function ydata = SSA_model(parameter, model, N1)
if model == 1
    lam=parameter(1);  gamma=parameter(2); nu=parameter(3); 
delta = 1;
%%steady time
TT=200;  H=200;  eps=0.01;
TP(1)=TT/H;
for i=2:1:H
   TP(i)=TT/H+TP(i-1) ;
end
odefun=@(t,y)[             
gamma*y(2)-lam*y(1);
lam*y(1)-gamma*y(2);
nu*y(2)-delta*y(3);
-(gamma+delta)*y(4)+lam*y(5)+nu*y(2);
gamma*y(4)-(lam+delta)*y(5);
-2*delta*y(6)+delta*y(3)+nu*y(2)+2*nu*y(4);
];
y00=[1,0,0,0,0,0];
options=odeset('reltol',1e-8,'abstol',1e-8);
[t,y]=ode45(odefun,[0,TT],y00,options);
 for j=1:1:H
   for i=1:1:length(t)
   if t(i)>=TP(j)
   break
   end
   k(j)=i+1;
   end
end
for j=1:1:H
meandata(j)=y(k(j),3);
seconddata(j)=y(k(j),6);
end
k=H; 
for j=1:1:H-1
   if  abs(meandata(H)-meandata(H-j))/meandata(H)<=eps
      k=k-1;
   else
       break
   end
end
meantime=k*(TT/H);
k1=H; 
for j=1:1:H-1
   if  abs(seconddata(H)-seconddata(H-j))/seconddata(H)<eps
      k1=k1-1;
   else
       break
   end
end
secondtime=k1*(TT/H);
tend=max(meantime,secondtime);

T = tend;               %% steady time
M=3*ceil(nu)+1;    %%the largest mRNA count
X(1) = 0; t(1) = 0;   % % the initial number of mRNA, the initial time
t1= T/5;     
t2 =2*T/5; 
t3=3*T/5; 
t4=4*T/5; 
t5=T; 

 S1 = []; S2 = []; S3=[]; S4=[]; 
S5=[];
        
for i = 1:N1
 s(1)=1;
    n = 1; 
    Xall = [X(1)]; Sall = [s(1)]; tall = [t(1)];
 
  while t(n) <= T ;
        h1 = 1; c1 = nu; a1 = h1*c1; % generate
        h2 = 1; c2 = gamma; a2 = h2*c2; % ON --> I1
        h3 = 1; c3 = lam; a3 = h3*c3;  % I1 --> ON 
        h4 = X(n); c4 = delta; a4 = h4*c4;  % decay
        n = n+1;
        
        if s(n-1) == 0; % ON
            a0 = a1+a2+a4;
            r1 = rand; r2=rand;
            tau = -log(r1)/a0;  % time interval in which nothing occurs
            if a0*r2 <= a1; % generate occurs
                X(n)=X(n-1)+1;  s(n)=0;
            elseif a0*r2 <=a1+a2;  % transition occurs
                X(n)=X(n-1);  s(n)=1;
            else % decay
                X(n)=X(n-1)-1; s(n)=0;
            end
            t(n) = t(n-1) + tau;
        else s(n-1) == 1; % I1
            a0 = a3+a4;
            r1=rand;  r2=rand;
            tau=-log(r1)/a0;
            if a0*r2 <= a3; % transition occurs
                X(n)=X(n-1); s(n)=0;
            else % decy occurs
                X(n)=X(n-1)-1; s(n)=1;
            end
            t(n) = t(n-1) + tau;
        end
        Xall = [Xall X(n)]; Sall = [Sall s(n)]; tall = [tall t(n)];
    end
    if t(n)<T;
        Xall=[Xall X(n)]; Sall=[Sall s(n)]; tall=[tall T];
    end
    for j=2:n
        if t(j-1) <=t1 && t(j) >t1;
            S1 = [S1 X(j-1)];
        end
        if t(j-1) <=t2 && t(j) >t2;
            S2 = [S2 X(j-1)];
        end
        if t(j-1) <=t3 && t(j) >t3;
            S3 = [S3 X(j-1)];
        end
        if t(j-1) <=t4 && t(j) >t4;
            S4 = [S4 X(j-1)];
        end
        if t(j-1) <=t5 && t(j) >t5;
            S5 = [S5 X(j-1)];
        end
    end
   
end
ydata = S5;

elseif model ==2
    lam1=parameter(1); lam2=parameter(2); q1=parameter(3); gamma=parameter(4); nu=parameter(4); q2=1-q1; delta=1;
vT=nu*delta;
gaT=gamma*delta;
la1=lam1*delta;
la2=lam2*delta;
ta=[]; yy=[]; 
TT=200;  H=200;  eps=0.01;
TP(1)=TT/H;
for i=2:1:H
   TP(i)=TT/H+TP(i-1) ;
end 
odefun=@(ta,yy)[
q1*gaT*yy(3)-la1*yy(1);
q2*gaT*yy(3)-la2*yy(2);
la1*yy(1)+la2*yy(2)-gaT*yy(3);
vT*yy(3)-delta*yy(4);
-(gaT+delta)*yy(5)+la1*yy(6)+la2*yy(7)+vT*yy(3);
q1*gaT*yy(5)-(la1+delta)*yy(6);
q2*gaT*yy(5)-(la2+delta)*yy(7);
-2*delta*yy(8)+delta*yy(4)+vT*yy(3)+2*vT*yy(5)
];
y00=[q1,q2,0,0,0,0,0,0];
options=odeset('reltol',1e-6,'abstol',1e-8);
[ta,yy]=ode45(odefun,[0,TT],y00,options);
for j=1:1:H
   for i=1:1:length(ta)
   if ta(i)>=TP(j)
   break
   end
   k(j)=i+1;
   end
end
for j=1:1:H
meandata(j)=yy(k(j),4);
seconddata(j)=yy(k(j),8);
end
k=H; 
for j=1:1:H-1
   if  abs(meandata(H)-meandata(H-j))/meandata(H)<=eps
      k=k-1;
   else
       break
   end
end
meantime=k*(TT/H);
k1=H; 
for j=1:1:H-1
   if  abs(seconddata(H)-seconddata(H-j))/seconddata(H)<eps
      k1=k1-1;
   else
       break
   end
end
secondtime=k1*(TT/H);

T=6*max(meantime,secondtime);  %% steady time
 M=3*ceil(nu)+1;
X(1) = 0; t(1) = 0;  
t1= 0; t2 =T/10; t3=2*T/10; t4=3*T/10; t5=4*T/10;t6=5*T/10; 
t7=6*T/10; t8=7*T/10; t9=8*T/10; t10=9*T/10; 
t11=T; 
S1 = []; S2 = []; S3=[]; S4=[]; S5=[];S6 = []; 
S7 = []; S8=[]; S9=[]; S10=[];S11 = [];

for i = 1:N1
    if rand<=q1 %  初始状态
        s(1)=1; 
    else
        s(1)=2;
    end
    n = 1; 
    Xall = [X(1)]; Sall = [s(1)]; tall = [t(1)];
    
    while t(n) <= T
        h1 = 1; c1 = nu; a1 = h1*c1; % generate
        h2 = 1; c2 = gamma; a2 = h2*c2; % ON --> 
        h3 = 1; c3 = lam1; a3 = h3*c3;  % I1 --> ON
        h4 = 1; c4 = lam2; a4 = h4*c4;  % I2 --> ON 
        h5 = X(n); c5 = delta; a5 = h5*c5;  % decay
        n = n+1;
        if s(n-1) == 0 % ON
            a0 = a1+a2+a5;
            r1 = rand; r2=rand;
            tau = -log(r1)/a0;  % time interval in which nothing occurs
            if a0*r2 <= a1 % generate occurs
                X(n)=X(n-1)+1;  s(n)=0;
            elseif a0*r2 <=a1+a2  % transition occurs
                X(n)=X(n-1);
                if rand <=q1
                    s(n)=1;
                else
                    s(n)=2;
                end
            else % decay
                X(n)=X(n-1)-1; s(n)=0;
            end
            t(n) = t(n-1) + tau;
        elseif s(n-1) == 1 % I1
            a0 = a3+a5;
            r1=rand;  r2=rand;
            tau=-log(r1)/a0;
            if a0*r2 <= a3 % transition occurs
                X(n)=X(n-1); s(n)=0;
            else % decy occurs
                X(n)=X(n-1)-1; s(n)=1;
            end
            t(n) = t(n-1) + tau;
        else  % I2
            a0 = a4+a5;
            r1=rand;  r2=rand;
            tau=-log(r1)/a0;
            if a0*r2 <= a4  % transition occurs
                X(n)=X(n-1); s(n)=0;
            else  % decay occurs
                X(n)=X(n-1)-1; s(n)=2;
            end
            t(n) = t(n-1) + tau;
        end
        Xall = [Xall X(n)]; Sall = [Sall s(n)]; tall = [tall t(n)];
    end
    if t(n)<T
        Xall=[Xall X(n)]; Sall=[Sall s(n)]; tall=[tall T];
    end
    for j=2:n
        if t(j-1) <=t1 && t(j) >t1
            S1 = [S1 X(j-1)];
        end
        if t(j-1) <=t2 && t(j) >t2
            S2 = [S2 X(j-1)];
        end
        if t(j-1) <=t3 && t(j) >t3
            S3 = [S3 X(j-1)];
        end
        if t(j-1) <=t4 && t(j) >t4
            S4 = [S4 X(j-1)];
        end
        if t(j-1) <=t5 && t(j) >t5
            S5 = [S5 X(j-1)];
        end
        if t(j-1) <=t6 && t(j) >t6
            S6 = [S6 X(j-1)];
        end
        if t(j-1) <=t7 && t(j) >t7
            S7 = [S7 X(j-1)];
        end
        if t(j-1) <=t8 && t(j) >t8
            S8 = [S8 X(j-1)];
        end
        if t(j-1) <=t9 && t(j) >t9
            S9 = [S9 X(j-1)];
        end
        if t(j-1) <=t10 && t(j) >t10
            S10 = [S10 X(j-1)];
        end
        if t(j-1) <=t11 && t(j) >t11
            S11 = [S11 X(j-1)];
        end
    end   
end
ydata = S11;

elseif model == 3
    lam1=parameter(1); lam2=parameter(2);  gamma=parameter(3); nu=parameter(4);  delta=1;
nu1=nu;
T = 200;    %% time at steady state
M=3*ceil(nu)+1;       
X(1) = 0; t(1) = 0;  % initial number and time
t1= 10;     
t2 =20; 
t3=30; 
t4=50; 
t5=T;  %% steady state
S1 = []; S2 = []; S3=[]; S4=[]; S5=[];
        
for i = 1:N1
 s(1)=1;
    n = 1; 
    Xall = [X(1)]; Sall = [s(1)]; tall = [t(1)];
  while t(n) <= T
        h1 = 1; c1 = nu1; a1 = h1*c1; % generate
        h2 = 1; c2 = gamma; a2 = h2*c2; % ON --> I1
        h3 = 1; c3 = lam1; a3 = h3*c3;  % I1 --> I2
        h4 = 1; c4 = lam2; a4 = h4*c4;  % I2 --> ON 
        h5 = X(n); c5 = delta; a5 = h5*c5;  % decay
        n = n+1;
            
        if s(n-1) == 0 % ON
            a0 = a1+a2+a5;
            r1 = rand; r2=rand;
            tau = -log(r1)/a0;  % time interval in which nothing occurs
            if a0*r2 <= a1 % generate occurs
                X(n)=X(n-1)+1;  s(n)=0;
            elseif a0*r2 <=a1+a2  % transition occurs
                X(n)=X(n-1);  s(n)=1;
            else % decay
                X(n)=X(n-1)-1; s(n)=0;
            end
            t(n) = t(n-1) + tau;
        elseif s(n-1) == 1 % I1
            a0 = a3+a5;
            r1=rand;  r2=rand;
            tau=-log(r1)/a0;
            if a0*r2 <= a3 % transition occurs
                X(n)=X(n-1); s(n)=2;
            else % decy occurs
                X(n)=X(n-1)-1; s(n)=1;
            end
            t(n) = t(n-1) + tau;
        else  % I2
            a0 = a4+a5;
            r1=rand;  r2=rand;
            tau=-log(r1)/a0;
            if a0*r2 <= a4  % transition occurs
                X(n)=X(n-1); s(n)=0;
            else  % decay occurs
                X(n)=X(n-1)-1; s(n)=2;
            end
            t(n) = t(n-1) + tau;
        end
        Xall = [Xall X(n)]; Sall = [Sall s(n)]; tall = [tall t(n)];
    end
    if t(n)<T
        Xall=[Xall X(n)]; Sall=[Sall s(n)]; tall=[tall T];
    end
    for j=2:n
        if t(j-1) <=t5 && t(j) >t5
            S5 = [S5 X(j-1)];
        end
    end
end
ydata = S5;

elseif model ==4
    lam=parameter(1); kappa1=parameter(2); kappa2=parameter(3); q1 = parameter(4); gamma=parameter(5); nu=parameter(6); q2=1-q1; delta=1;

TT=200;  H=200;  eps=0.01;
TP(1)=TT/H;
for i=2:1:H
   TP(i)=TT/H+TP(i-1) ;
end
odefun=@(t,y)[
q1*gamma*y(4) - kappa1*y(1);
q2*gamma*y(4) - lam*y(2);
kappa1*y(1) - kappa2*y(3)
lam*y(2) + kappa2*y(3) - gamma*y(4);
nu*y(4) - delta*y(5);
];
y00=[q1,q2,0,0,0];
options=odeset('reltol',1e-8,'abstol',1e-8);
[t,y]=ode45(odefun,[0,TT],y00,options);
 y(end,5)
 for j=1:1:H
   for i=1:1:length(t)
   if t(i)>=TP(j)
   break
   end
   k(j)=i+1;
   end
end
for j=1:1:H
meandata(j)=y(k(j),5);
end
k2=H; 
for j=1:1:H 
   if  abs(meandata(H)-meandata(H-j))/meandata(H)<=eps
      k2=k2-1;
   else
       break
   end
end
meantime=k2*(TT/H);
tend = 2*max(meantime)

T = tend; 
M=3*ceil(nu)+1;
t8=tend;
X(1) = 0; t(1) = 0;  
S8=[]; 
for i = 1:N1
    if rand<=q1 
        s(1)=1; 
    else
        s(1)=2;
    end
    n = 1; 
    Xall = [X(1)]; Sall = [s(1)]; tall = [t(1)];
    while t(n) <= T
        h1 = 1; c1 = nu; a1 = h1*c1; % generate
        h2 = 1; c2 = gamma; a2 = h2*c2; % ON --> 
        h3 = 1; c3 = kappa1; a3 = h3*c3;  % I1 --> I3
        h4 = 1; c4 = kappa2; a4 = h4*c4;  % I3 --> ON
        h5 = 1; c5 = lam; a5 = h5*c5;  % I2 --> ON 
        h6 = X(n); c6 = delta; a6 = h6*c6;  % decay
        n = n+1;
        if s(n-1) == 0 % ON
            a0 = a1+a2+a6;
            r1 = rand; r2=rand;
            tau = -log(r1)/a0;  % time interval in which nothing occurs
            if a0*r2 <= a1 % generate occurs
                X(n)=X(n-1)+1;  s(n)=0;
            elseif a0*r2 <=a1+a2  % transition occurs
                X(n)=X(n-1);
                if rand <=q1  % transition occurs -->I1
                    s(n)=1;
                else
                    s(n)=2;  % transition occurs -->I2
                end
            else % decay
                X(n)=X(n-1)-1; s(n)=0;
            end
            t(n) = t(n-1) + tau;
        elseif s(n-1) == 1 % I1
            a0 = a3+a6;
            r1=rand;  r2=rand;
            tau=-log(r1)/a0;
            if a0*r2 <= a3 % transition occurs -->I3
                X(n)=X(n-1); s(n)=3;
            else % decy occurs
                X(n)=X(n-1)-1; s(n)=1;
            end
            t(n) = t(n-1) + tau;
        elseif  s(n-1) == 2 % I2
            a0 = a5+a6;
            r1=rand;  r2=rand;
            tau=-log(r1)/a0;
            if a0*r2 <= a5  % transition occurs ON
                X(n)=X(n-1); s(n)=0;
            else  % decay occurs
                X(n)=X(n-1)-1; s(n)=2;
            end
            t(n) = t(n-1) + tau;
        else % I3
            a0 = a4+a6;
            r1=rand;  r2=rand;
            tau=-log(r1)/a0;
            if a0*r2 <= a4 % transition occurs -->ON
                X(n)=X(n-1); s(n)=0;
            else % decy occurs
                X(n)=X(n-1)-1; s(n)=3;
            end
            t(n) = t(n-1) + tau;
        end
        Xall = [Xall X(n)]; Sall = [Sall s(n)]; tall = [tall t(n)];
    end
    if t(n)<T
        Xall=[Xall X(n)]; Sall=[Sall s(n)]; tall=[tall T];
    end
    for j=2:n
        if t(j-1) <=t8 && t(j) >t8
            S8 = [S8 X(j-1)];
        end
    end
    
end
ydata = S8;

end
end
