
clc
clear all

[data,name]=xlsread("Fibroblasts_c57_fitmethod_nofixp.csv");
[data1,name1]=xlsread("Fibroblasts_c57_distribution.csv");
AICc=[data(:,10) data(:,24) ];
AICc1=min(AICc,[],2);
type=[data(:,16) data(:,30) ];
row=find( (AICc(:,1)==AICc1 & type(:,1)==3) | (AICc(:,2)==AICc1 & type(:,2)==3) );
data_b=data(row ,:); name_b=name(row+1,1); AICc_b=AICc(row,:); AICc_b1=AICc1(row); 
data2=data_b; name2=name_b;

is_member = ismember(name1, name2); % 找name1中哪些在name2中
data3 = data1(:,is_member); 
clearvars -except data3 data2 name2

% %%% Ecoli
% [data2,name2]=xlsread("Ecoli fitmethod.xlsx");
% [data,name]=xlsread("data_bimodal.xlsx");
% data3=data(26,1:31);
% 
% %%% HIV
% % HIV
% [data2,name2]=xlsread("HIV.xlsx", 'fit'); 
% [data,name]=xlsread("HIV.xlsx", 'counts');
% data3=data(:,end);
SS = length(name2);
kb1 = nan(SS,1); kb2 = nan(SS,1);
for zushu=1: SS  
    
%     data1=clearnan(data3);
%     for j=1:max(data1)+1
%         xdata(j)=sum(data1==(j-1))/length(data1);
%     end
    clear xdata
   xdata=clearnan(data3(5:end,zushu));
   
   S=length(xdata);
   S1=S;
%    figure(zushu)
%    bar(0:S1-1, xdata(1:S1))
%    hold on
%%% two-state
   parameter=data2(zushu,6:9);
la1=parameter(1); ga1=parameter(2); v1=parameter(3); %mu=parameter(4);
M=S-1;
N = 2*M+10; d=1;

% FSP
num = 2*N;
Q = zeros(num);
one = ones(num,1);
for i = 1:N-1
    Q(i+1,i) = i*d;
    Q(i+N+1,i+N) = i*d;
    Q(i+N,i+N+1) = v1;
end
for i = 1:N
    Q(i,i+N) = la1;
    Q(i+N,i) = ga1;
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
 
for t=1:N
    if dist(t)<=1e-15
        dist(t)=0;
    end
end
x=dist;
%  beta_b =parameter(end);
%     x = zeros(1, M+1);
%    for t = 1:M+1
%        i_range = t:N;
%        x(t) = x(t) + sum(binopdf(t-1, i_range-1, beta_b).*dist(i_range));
%    end
%    flag=0;
%  for i=3:M+1
%      if x(i-1)>=x(i) 
%         flag=flag+0;
%      else  
%         flag=flag+1;        
%      end
%  end
 hhigh =NaN; hlow =NaN; valley =NaN;
 hhigh=x(1);
for i=3:M
    if x(i-1)<x(i) && x(i)>=x(i+1)       
        hlow = min([hhigh,x(i)]);
        hhigh = max([hhigh,x(i)]);
    end 
end
for i=2:M
    if x(i)<x(i-1) && x(i)<x(i+1)
        valley = x(i);
    end
end
if all(~isnan([hlow valley hhigh]))
   kb1(zushu)=(hlow-valley)/hhigh;
end

% plot(0:S1-1,x(1:S1 ),'b')
% hold on
clearvars -except data2 name2 zushu data3 S S1 M kb1 kb2 SS

%%%两路径

parameter=data2(zushu,18:23);
latr1=parameter(1); latr2=parameter(2); qr1=parameter(3); gar=parameter(4); vr=parameter(5); 
qr2=1-qr1;
rho1=0; rho2=0;

d = 1; 
% transition matrix with all diagnoal elements being zero
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

for t=1:N
    if dist(t)<=1e-15
        dist(t)=0;
    end
end
x=dist;
% beta_b =parameter(end);
%     x = zeros(1, M+1);
%    for t = 1:M+1
%        i_range = t:N;
%        x(t) = x(t) + sum(binopdf(t-1, i_range-1, beta_b).*dist(i_range));
%    end
% 
%    flag=0;
%  for i=3:M+1
%      if x(i-1)>=x(i) 
%         flag=flag+0;
%      else  
%         flag=flag+1;        
%      end
%  end
  hhigh =NaN; hlow =NaN; valley =NaN;
 hhigh=x(1);
for i=3:M
    if x(i-1)<x(i) && x(i)>=x(i+1)       
        hlow = min([hhigh,x(i)]);
        hhigh = max([hhigh,x(i)]);
    end 
end
for i=2:M
    if x(i)<x(i-1) && x(i)<x(i+1)
        valley = x(i);
    end
end
if all(~isnan([hlow valley hhigh]))
   kb2(zushu)=(hlow-valley)/hhigh;
end
% plot(0:S1-1,x(1:S1 ),'black')
% hold on
clearvars -except data2 name2 zushu data3 S S1 M kb1 kb2 SS

% %%%threestate
% 
% parameter = data2(zushu,32:36);
% 
% la1=parameter(1);  la2=parameter(2);   ga=parameter(3);  v=parameter(4);
% rho0 = v; rho1=0; rho2=0;
% 
% % transition matrix with all diagnoal elements being zero
% K=[0 0 ga;
%    la1 0 0;
%    0 la2 0];
% 
% N = 2*M+10; d=1;
% num = 3*N;
% Q = zeros(num);
% one = ones(num,1);
% dist = zeros(1,N);
% for i = 1:N-1
%     Q(i+1,i) = i*d;
%     Q(i+1+N,i+N) = i*d;
%     Q(i+1+2*N,i+2*N) = i*d;
%     
%     Q(i,i+1) = rho0;
%     Q(i+N,i+1+N) = rho1;
%     Q(i+2*N,i+1+2*N) = rho2;
% end
% for i = 1:N
%     for j=1:3
%         for k=1:3
%             Q(i+(j-1)*N,i+(k-1)*N) = K(j,k);
%         end
%     end
% end
% 
% temp = Q*one;
% for i = 1:num
%     Q(i,i) = -temp(i);
% end
% 
% Qmod = Q;
% for i = 1:num
%     Qmod(i,num) = 1;
% end
% vec = zeros(1,num);
% vec(num) = 1;
% mu = vec/Qmod;
% 
% for i = 1:N
%     dist(i) = dist(i)+mu(i)+mu(i+N)+mu(i+2*N);
% end
% 
% for t=1:N
%     if dist(t)<=1e-15
%         dist(t)=0;
%     end
% end
% %x=dist;
% beta_b =parameter(end);
%     x = zeros(1, M+1);
%    for t = 1:M+1
%        i_range = t:N;
%        x(t) = x(t) + sum(binopdf(t-1, i_range-1, beta_b).*dist(i_range));
%    end
% 
% plot(0:S1-1,x(1:S1 ),'green')
% hold on
% clearvars -except data2 name2 zushu data3 S S1 M
% 
% %%% cross-talk three-state
% parameter=data2(zushu, 45:51);
% lam=parameter(1); kappa1=parameter(2); kappa2=kappa1; q1=parameter(4); gamma=parameter(5); nu=parameter(6);
% 
% q2=1-q1;
% delta = 1; d=delta;
% K=[0 q1*gamma q2*gamma 0;
%    0 0 0 kappa1;
%    lam 0 0 0;
%    kappa2 0 0 0];
% N = 2*M+10; num = 4*N;
% Q = zeros(num);
% one = ones(num,1);
% dist = zeros(1,N);
% for i = 1:N-1
%     Q(i+1,i) = i*d;
%     Q(i+1+N,i+N) = i*d;
%     Q(i+1+2*N,i+2*N) = i*d;
%     Q(i+1+3*N,i+3*N) = i*d;
%     Q(i,i+1) = nu;
%    
% end
% for i = 1:N
%     for j=1:4
%         for k=1:4
%             Q(i+(j-1)*N,i+(k-1)*N) = K(j,k);
%         end
%     end
% end
% temp = Q*one;
% for i = 1:num
%     Q(i,i) = -temp(i);
% end
% Qmod = Q;
% for i = 1:num
%     Qmod(i,num) = 1;
% end
% vec = zeros(1,num);
% vec(num) = 1;
% mu = vec/Qmod;
% for i = 1:N
%     dist(i) = dist(i)+mu(i)+mu(i+N)+mu(i+2*N)+mu(i+3*N);
% end
% for t=1:N
%     if dist(t)<=1e-15
%         dist(t)=0;
%     end
% end
% %x=dist;
% beta_b =parameter(end);
%     x = zeros(1, M+1);
%    for t = 1:M+1
%        i_range = t:N;
%        x(t) = x(t) + sum(binopdf(t-1, i_range-1, beta_b).*dist(i_range));
%    end
% 
% plot(0:S1-1,x(1:S1 ),'red')
% title(['筛为crosstalk',name2(zushu)])
% legend('样本数据',  'twostate', 'crosstalk', 'threestate', 'crosstalk-threestate')
% clearvars -except data2 name2 zushu data3  kb1
%       
end


