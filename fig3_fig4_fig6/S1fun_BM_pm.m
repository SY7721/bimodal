function S1 = S1fun_BM_pm(b1)
    global  xdata tdata tt     
    M=max(tdata); %b1=log(parameter(1:4));
    rho = exp(b1(3)) ; d = 1; sigma1 =  exp(b1(1)); sigma0 =  exp(b1(2)); 
    if any([sigma0, sigma1, rho] >= 1e10) || any([sigma0, sigma1, rho] <= 1e-10)
        S1 = 1e30;
        return;
    end
    % FSP
    N=2*M+10;
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

    dist(dist<1e-15)=0;

    beta_b =1/(1+ exp(b1(4)));
    x = zeros(1, M+1);
   for t = 1:M+1
       i_range = t:N;
       x(t) = x(t) + sum(binopdf(t-1, i_range-1, beta_b).*dist(i_range));
   end
   x(x<=0) = 1e-17;
    % figure;
    % bar(xdata)
    % hold on
    % plot(x)
    S0 = 0;
    for i = 1:M+1
        S0 = S0 - tt(i).*log(x(i));
    end
    S1=S0;
end
