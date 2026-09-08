function [mean,  fano,   skewness, kurtosis, HD, y]=caculate_moment(parameter,model)
        global tdata xdata  tt 
      
        M=length(xdata)-1;
     if model==1   %% telegraph
        la=parameter(1); ga=parameter(2); v=parameter(3);
        N = 2*M+10; d=1; num = 2*N;
        Q = zeros(num);
        one = ones(num,1);
        for i = 1:N-1
            Q(i+1,i) = i*d;
            Q(i+N+1,i+N) = i*d;
            Q(i+N,i+N+1) = v;
        end
        for i = 1:N
            Q(i,i+N) = la;
            Q(i+N,i) = ga;
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

     elseif model==2
         la1=parameter(1);  la2=parameter(2);  q1=parameter(3);  ga=parameter(4);  v=parameter(5);
         rho1=0; rho2=0; d = 1; 
         K=[0 ga*q1 (1-q1)*ga;
            la1 0 0;
            la2 0 0];
         N =2*M+10; num = 3*N;
         Q = zeros(num);
         one = ones(num,1);
         dist = zeros(1,N);
         for i = 1:N-1
             Q(i+1,i) = i*d;
             Q(i+1+N,i+N) = i*d;
             Q(i+1+2*N,i+2*N) = i*d;  
             Q(i,i+1) = v;
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
     elseif model==3
         la1=parameter(1);  la2=parameter(2);   ga=parameter(3);  v=parameter(4);
         rho0 = v; rho1=0; rho2=0;
         d = 1; 
         K=[0 0 ga;
            la1 0 0;
            0 la2 0];
        N = 2*M+10;  num = 3*N;
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

     elseif model==4
         lam=parameter(1); kappa1=parameter(2); kappa2=parameter(3); q1=parameter(4); gamma=parameter(5); nu=parameter(6);
         q2=1-q1;
         d=1;
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
     end

     dist(dist<=1e-15) = 0;
     S = M+1;
     x = dist(1:S);    
     mean=0; twom1=0; threem=0; fourm=0;
    for t=1:1:S
        mean=mean+(t-1).*x(t);
        twom1=twom1+((t-1)^2).*x(t);
        threem=threem+((t-1)^3).*x(t);
        fourm=fourm+((t-1)^4).*x(t);
    end
    fano=twom1/mean-mean;
    Var=0; threec=0;fourc=0;
    for i=1:1:S
        Var=Var+((i-1-mean)^2).*x(i);
        threec=threec+((i-1-mean)^3).*x(i);
        fourc=fourc+((i-1-mean)^4).*x(i);
    end
    skewness=threec/(Var^(3/2));    
    kurtosis=fourc/(Var^2)-3;
    
    %%%% HD
     H1=0;
     for t=1:S
         H1=H1+(sqrt(xdata(t))-sqrt(x(t)))^2;
     end
     HD=sqrt(H1/2);


     %%% type of distribution
     flag=0; nn=N;
     for i=3:N
         if dist(i-1)>=dist(i) 
            flag=flag+0;
         else  
            flag=flag+1;
            nn=i;
         end
     end
 
     flag1=0;
     if nn==N
        flag1=0;
     elseif nn<N
        for i=nn:N-1
            if dist(i)>dist(i+1)
               flag1=flag1+1;
            end
        end
     end
 

if dist(1)<dist(2) %%单
    y=1;
end
if dist(1)==dist(2) && flag>0
    y=1;
end
if dist(1)>dist(2) && flag>0 && flag1>0 %%双峰
    y=3;
end
if dist(1)>dist(2) && flag>0 && flag1==0 %%双峰
    y=2;
end
if dist(1)>=dist(2) && flag==0 %%递减
    y=2;
end

    
end

