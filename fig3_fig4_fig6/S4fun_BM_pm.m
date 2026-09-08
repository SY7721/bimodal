function S4 = S4fun_BM_pm(b4)
    global  xdata tdata tt     
    lam=exp(b4(1)); kappa1=exp(b4(2)); kappa2=kappa1; q1=1/(1+exp(b4(3))); gamma=exp(b4(4)); nu=exp(b4(5));
     if any([lam, kappa1, gamma, nu] >= 1e10) || any([lam, kappa1, gamma, nu] <= 1e-10)
        S4 = 1e30;
        return;
    end
q2=1-q1;
delta = 1; d=delta;
K=[0 q1*gamma q2*gamma 0;
   0 0 0 kappa1;
   lam 0 0 0;
   kappa2 0 0 0];
M=max(tdata);
N =2*M+10; num = 4*N;
Q = zeros(num);
one = ones(num,1);
dist = zeros(1,N);
for i = 1:N-1
    Q(i+1,i) = i*d;
    Q(i+1+N,i+N) = i*d;
    Q(i+1+2*N,i+2*N) = i*d;
    Q(i+1+3*N,i+3*N) = i*d;
    Q(i,i+1) = nu; 
end
for i = 1:N
    for j=1:4
        for k=1:4
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
    dist(i) = dist(i)+mu(i)+mu(i+N)+mu(i+2*N)+mu(i+3*N);
end
    for i=1:N
        if dist(i)<1e-15
           dist(i)=0;
        end
    end

    beta_b = 1/(1+ exp(b4(end)));
    x = zeros(1, M+1);
    for t = 1:M+1
        i_range = t:N;
        x(t) = x(t) + sum(binopdf(t-1, i_range-1, beta_b).*dist(i_range));
    end
    x(x<=0) = 1e-17;
    S0 = 0;
    for i = 1:M+1
        S0 = S0 - tt(i).*log(x(i));
    end
   S4 = S0;
end
