function S2 = S2fun_BM_pm(b2)
       global  xdata tdata tt 
       M=max(tdata);
       d = 1; 
       rho0 = exp(b2(5)); rho1=0; rho2=0; 
       if any([rho0, exp(b2(1)), exp(b2(2)), exp(b2(4))] >= 1e+10) || ...
          any([rho0, exp(b2(1)), exp(b2(2)), exp(b2(4))] <= 1e-10)
          S2 = 1e30;
          return;
       end
       K=[0 exp(b2(4))/(1+exp(b2(3))) exp(b2(3))*exp(b2(4))/(1+exp(b2(3)));
          exp(b2(1)) 0 0;
          exp(b2(2)) 0 0];
       N=2*M+10;
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
        dist(dist<1e-15)=0;

        beta_b = 1/(1+ exp(b2(6)));
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
S2=S0;

end





