function S1 = S1fun_pm1(b)
%%%%¼«´óËÆÈ»

global tdata xdata  tt

rho =exp(b(3)); d = 1; sigma1 = exp(b(1)); sigma0 = exp(b(2));
N = 2*(max(tdata))+10;
% FSP
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
for i=1:N
    if dist(i)<=1e-15
        dist(i)=0;
    end
end 
for t = 1:1:max(tdata)+1;
      x(t)=dist(t);
end

S1 = 0;
for i = 1:length(tdata)
       S1 = S1 - tt(i).*log(x(i)+(1e-17));
end
end