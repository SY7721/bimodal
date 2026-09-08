function [mean,  fano,   skewness, kurtosis,  y]=caculate_moment(parameter,model)
        global tdata xdata  tt 
      
        M=length(xdata)-1;
     if model==1   %% telegraph
        la=parameter(1); ga=parameter(2); v=parameter(3);
        N = 3*M+10; d=1; num = 2*N;
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
         N =3*M+10; num = 3*N;
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
     end
     dist(dist<=1e-15) = 0;
     S=N;
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

