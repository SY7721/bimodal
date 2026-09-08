function S5 = S5fun_BM_pm(b5)
        global tdata xdata  tt

        rho =exp(b5(3)); d = 1; sigma1 = exp(b5(1)); sigma0 = exp(b5(2)); mu=exp(b5(4));
        M=max(tdata);
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
            Q(i,i+N) = sigma1+(i-1)*mu;
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

        beta_b =1/(1+ exp(b5(5)));
        x = zeros(1, M+1);
        for t = 1:M+1
            i_range = t:N;
            x(t) = x(t) + sum(binopdf(t-1, i_range-1, beta_b).*dist(i_range));
        end

         S0 = 0;
         for i = 1:M+1
             S0 = S0 - tt(i).*log(x(i)+1e-17);
         end
         if rho>=1e+10 || sigma1>=1e+10 || sigma0>=1e+10 || rho<=1e-10 || sigma1<=1e-10 || sigma0<=1e-10 || mu>=1e+10 || mu<=1e-10 || sum(x)<0.99
             S5=S0+1e+30;
         else
             S5=S0;
         end




%% 数值计算

%% plot result of the integration
% figure(1)
% plot(tdata,xdata,'r.','MarkerSize',20);
% hold on
% plot(tdata,x,'k','linewidth',2);
%  hold off
% drawnow



end