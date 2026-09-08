function dist = fsp_model(parameter, model, M)

if model == 1 
    rho =parameter(3); d = 1; sigma1 = parameter(1); sigma0 = parameter(2);
N = 2*M+10;
num = 2*N;
Q = zeros(num);
one = ones(num,1);
for i = 1:N-1
    Q(i+1,i) = i*d;
    Q(i+N+1,i+N) = i*d;
    Q(i+N,i+N+1) = rho;
end
for i = 1:N
    Q(i,i+N) = sigma1;
    Q(i+N,i) = sigma0;
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
ssd = vec/Qmod;
dist = ssd(1:N)+ssd(N+1:2*N);
dist(dist <= 1e-15) = 0;

elseif model ==2
       latr1=parameter(1); latr2=parameter(2); qr1=parameter(3); gar=parameter(4); vr=parameter(5); 
       qr2=1-qr1;
       rho1=0; rho2=0; d = 1; 
K=[0 gar*qr1 (1-qr1)*gar;
   latr1 0 0;
   latr2 0 0];
N =2*M+10; num = 3*N;
Q = zeros(num);
one = ones(num,1);
dist = zeros(1,N);
for i = 1:N-1
    Q(i+1,i) = i*d;
    Q(i+1+N,i+N) = i*d;
    Q(i+1+2*N,i+2*N) = i*d;
    
    Q(i,i+1) = vr;
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
dist(dist <= 1e-15) = 0;

elseif model == 3
    la1=parameter(1);  la2=parameter(2);   ga=parameter(3);  v=parameter(4);
    rho0 = v; rho1=0; rho2=0;
K=[0 0 ga;
   la1 0 0;
   0 la2 0];
N = 2*M+10; d=1;
num = 3*N;
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
dist(dist <= 1e-15) = 0;

elseif model ==4
    lam=parameter(1); kappa1=parameter(2); kappa2=kappa1; q1=parameter(4); gamma=parameter(5); nu=parameter(6);
    q2=1-q1; delta = 1; d=delta;
K=[0 q1*gamma q2*gamma 0;
   0 0 0 kappa1;
   lam 0 0 0;
   kappa2 0 0 0];
N = 2*M+10; num = 4*N;
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
dist(dist <= 1e-15) = 0;

end
end