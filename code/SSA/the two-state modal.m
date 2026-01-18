
N1=100;   %% sample size
 lam=0.1+4.9*rand(); gamma=0.1+4.9*rand(); nu=5+45*rand(); delta = 1;
 
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

for i = 1:M
    yy(i) = sum(S5==(i-1))./N1;
end
xdata=[]; tdata=[]; tt=[];
xdata=yy;
tt=N1.*xdata;
S=numel(xdata);
tdata=[];
for i=1:1:S;
tdata(i)=i-1;
end
xdata;
