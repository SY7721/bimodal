clc
clear all
N1=100; %% sample size
 lam1=0.1+4.9*rand(); lam2=0.1+4.9*rand();  gamma=0.1+4.9*rand(); nu=5+45*rand();  delta=1;
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
for i = 1:M
    y(i) = sum(S5==(i-1))/N1;
end
xdata=y;


