function S2 = S2fun_pm1(b2)
global tdata xdata  tt

rho0 = exp(b2(5)); rho1=0; rho2=0; d = 1; 
% transition matrix with all diagnoal elements being zero
K=[0 exp(b2(4))/(1+exp(b2(3))) exp(b2(3))*exp(b2(4))/(1+exp(b2(3)));
   exp(b2(1)) 0 0;
   exp(b2(2)) 0 0];

N =2*(max(tdata))+10; num = 3*N;
Q = zeros(num);
one = ones(num,1);
dist = zeros(1,N);
for i = 1:N-1
    Q(i+1,i) = i*d;
    Q(i+1+N,i+N) = i*d;
    Q(i+1+2*N,i+2*N) = i*d;
    
    Q(i,i+1) = rho0;
    Q(i+N,i+1+N) = rho1;
    Q(i+2*N,i+1+2*N) = rho2;
end
for i = 1:N
    for j=1:3
        for k=1:3
            Q(i+(j-1)*N,i+(k-1)*N) = K(j,k);
        end
    end
end
temp = Q*one;
for i = 1:num
    Q(i,i) = -temp(i);
end
Qmod = Q;
for i = 1:num
    Qmod(i,num) = 1;
end
vec = zeros(1,num);
vec(num) = 1;
mu = vec/Qmod;
for i = 1:N
    dist(i) = dist(i)+mu(i)+mu(i+N)+mu(i+2*N);
end
for i=1:N
    if dist(i)<=1e-15
        dist(i)=0;
    end
end
for t = 1:1:max(tdata)+1;
      x(t)=dist(t);
end

S2 = 0;
for i = 1:length(tdata)
       S2 = S2 - tt(i).*log(x(i)+(1e-17));
end
end




