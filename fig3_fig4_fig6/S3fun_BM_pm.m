function S3 = S3fun_BM_pm(b3)
    global  xdata tdata tt     
    rho0 = exp(b3(4)); rho1=0; rho2=0;
     if any([exp(b3(1)), exp(b3(2)), exp(b3(3)), exp(b3(4))]>=1e+10) || any([exp(b3(1)), exp(b3(2)), exp(b3(3)), exp(b3(4))]<=1e-10) 
        S3=1e30;
        return;
    end
    d = 1; 
     % transition matrix with all diagnoal elements being zero
    K=[0 0 exp(b3(3));
       exp(b3(2)) 0 0;
       0 exp(b3(1)) 0];
    M=max(tdata);
    N =2*M+10; num = 3*N;
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
        if dist(i)<1e-15
           dist(i)=0;
        end
    end

    beta_b = 1/(1+ exp(b3(end)));
    x = zeros(1, M+1);
    for t = 1:M+1
        i_range = t:N;
        x(t) = x(t) + sum(binopdf(t-1, i_range-1, beta_b).*dist(i_range));
    end
    x(x<=0) = 1e-17;
%     figure;
%     bar(xdata)
%     hold on
%     plot(x)
    S0 = 0;
    for i = 1:M+1
        S0 = S0 - tt(i).*log(x(i));
    end
    S3 = S0;
end
