clc
clear all

[data, name] = xlsread("CTT1+0.4 Nacl.xlsx");
time = [0 1  10 15 35 40  50 55];
max_data = max(data(:));

la = 0.0426945; ga=	0.371866737	; v= 17.35678991; 
delta=0.063;
M = max(2 * max_data+10, max_data + 1);
    I = speye(M);
    Z = sparse(M, M);
    n = (0:M-1).';
    Adeg = spdiags(-delta * n, 0, M, M);
    if M > 1
        Adeg = Adeg + sparse((1:M-1).', (2:M).', ...
                            delta * (1:M-1).', M, M);
    end
    Aprod = -v * I;
    if M > 1
        Aprod = Aprod + sparse((2:M).', (1:M-1).', ...
                              v * ones(M-1, 1), M, M);
    end
    A00 = Adeg - la * I;
    A01 = ga * I;
    A10 = la * I;
    A11 = Adeg + Aprod - ga * I;
    A = [A00, A01; A10, A11];
    Pinit = zeros(2 * M, 1);
    Pinit(1) = 1;
    ode_options = odeset( ...
        'RelTol', 1e-5, ...
        'AbsTol', 1e-10, ...
        'Jacobian', A, ...
        'JPattern', spones(A));

    try
        [~, P] = ode15s(@(t, y) A * y, time(:), Pinit, ode_options);
    catch
        S1 = 1e30;
        return
    end
    if any(~isfinite(P(:)))
        S1 = 1e30;
        return
    end
    P_total = P(:, 1:M) + P(:, M+1:2*M);
clearvars -except P_total time max_data data
  
parameter = [0.000174926	0.290584908	0.066562045	1.301519093	30.72509239];
kon1 =parameter(1); kon2 = parameter(2); q1   =parameter(3);  koff =parameter(4); kb =parameter(5); 

q2 = 1 - q1; delta = 0.063;
    M = max(5 * max_data+10, max_data + 1);

    I = speye(M);
    Z = sparse(M, M);
    n = (0:M-1).';
    Adeg = spdiags(-delta * n, 0, M, M);
    if M > 1
        Adeg = Adeg + sparse((1:M-1).', (2:M).', ...
                            delta * (1:M-1).', M, M);
    end
    Aprod = -kb * I;
    if M > 1
        Aprod = Aprod + sparse((2:M).', (1:M-1).', ...
                              kb * ones(M-1, 1), M, M);
    end
    A11 = Adeg - kon1 * I;
    A22 = Adeg - kon2 * I;
    A33 = Adeg + Aprod - koff * I;

    A = [ ...
        A11,        Z, q1 * koff * I; ...
           Z,      A22, q2 * koff * I; ...
        kon1 * I, kon2 * I, A33];

    Pinit = zeros(3 * M, 1);
    Pinit(1) = q1;
    Pinit(M + 1) = q2;

    ode_options = odeset( ...
        'RelTol', 1e-5, ...
        'AbsTol', 1e-10, ...
        'Jacobian', A, ...
        'JPattern', spones(A));

    try
        [~, P] = ode15s(@(t, y) A * y, time(:), Pinit, ode_options);
    catch
        S2 = 1e30;
        return
    end

    if any(~isfinite(P(:)))
        S2 = 1e30;
        return
    end

    P_total1 = P(:, 1:M) + P(:, M+1:2*M) + P(:, 2*M+1:3*M);

    %%% figS1(a)
for i=1:length(time)
    a = clearnan(data(:,i));
    CN = length(a);
    for j = 1:max(a)+1
        xdata(j)  = sum(a == (j-1))./CN;
    end
   
    L = length(xdata);
    dist1 = P_total(i, 1:L);
    dist2 = P_total1(i, 1:L);
    LL = [14 14 104 204 218 214 188 130];
    yy = [1 1 0.14 0.1 0.09 0.09 0.16 0.18];
    SS =LL(i);
    figure(i)
    bar(0: SS-1, xdata(1:SS))
    hold on
    plot(0: SS-1,dist1(1:SS), 'r')
    hold on 
    plot(0: SS-1,dist2(1:SS), 'b')
    hold off
    xlim([-3 SS-4])
    ylim([0 yy(i)])
    clear xdata tt L dist1 dist2  a  CN
legend("CTT1", "Two-state", "cross-talk pathway")
end


%%% figS1(b)
clc 
clear all
%%% Ecoli
data=xlsread('Ecoli fitmethod.xlsx', 'counts');
xdata = data(2,1:31);
clear data
CN  = 100; % cell number
tt= round(CN.*xdata);
S=numel(xdata);
for i=1:1:S;
    tdata(i)=i-1;
end
M=S-1;
figure(2)
S1=M;
bar(0:S1-1,xdata(1:S1))
hold on

global  xdata tt tdata 
model = 1;
parameter = MLE(model);
Smin_min=parameter(end);
AICc = caculate_AICc(Smin_min,3,CN);
parameter1 = parameter(1:end-1);
dist = fsp_model(parameter1, model, M);
plot(0:S1-1, dist(1:S1), 'r')
hold on

clear parameter parameter1
model = 2;
parameter = MLE(model);
Smin_min=parameter(end);
AICc1 = caculate_AICc(Smin_min,4,CN);
parameter1 = parameter(1:end-1);
dist1 = fsp_model(parameter1, model, M);
plot(0:S1-1, dist1(1:S1), 'black')
hold on
legend('E.coli', 'two-state', 'cross-talk pathway')

clc
clear all
%%% HIV
xdata =[]; tt=[]; tdata = [];
[data2,name2]=xlsread("HIV.xlsx", 'counts');
data1=clearnan(data2); 
for j=1:max(data1)+1
    xdata(j)=sum(data1==(j-1))/length(data1);
end
CN = length(data1);
tt= round(CN.*xdata);
S=numel(xdata);
for i=1:1:S;
    tdata(i)=i-1;
end
M=S-1;

figure(3)
S1=M;
bar(0:S1-1,xdata(1:S1))
hold on

global  xdata tt tdata 
model = 1;
parameter = MLE(model);
Smin_min=parameter(end);
AICc = caculate_AICc(Smin_min,3,CN);
parameter1 = parameter(1:end-1);
dist = fsp_model(parameter1, model, M);
plot(0:S1-1, dist(1:S1), 'r')
hold on

clear parameter parameter1
model = 2;
parameter = MLE(model);
Smin_min=parameter(end);
AICc1 = caculate_AICc(Smin_min,4,CN);
parameter1 = parameter(1:end-1);
dist1 = fsp_model(parameter1, model, M);
plot(0:S1-1, dist1(1:S1), 'black')
hold on
legend('HIV', 'two-state', 'cross-talk pathway')



