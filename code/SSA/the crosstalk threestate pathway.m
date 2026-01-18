function the crosstalk threestate pathway

clc
close all
clear all
N=1000;  %% sample size
lam = lam2_cross(zushu);  kappa1 = 0.2;kappa2 = 0.2; q1 = q1_cross(zushu); gamma = ga_cross(zushu); nu = v_cross(zushu); delta = 1; 

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
for i = 1:N
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
    
for i = 1:M
    y8(i) = sum(S8==(i-1))/N;
end

xdata=y8;

end

