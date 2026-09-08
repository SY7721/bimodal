clc
clear all
close all

nn = 3; %% 1: F c57; 2: F cast; 3: E c57; 4: E cast
 if nn == 1    
    [data,name]=xlsread("Fibroblasts_c57_fitmethod_nofixp.csv");    
    [~,name1] = xlsread("bimodal_gene.xlsx", 'F c57');  
    nname = 'F c57';
elseif nn == 2   
    [data,name]=xlsread("Fibroblasts_cast_fitmethod_nofixp.csv");    
    [~,name1] = xlsread("bimodal_gene.xlsx", 'F cast');   
    nname = 'F cast'; 
elseif nn == 3    
    [data,name]=xlsread("Embryonic_c57_fitmethod_nofixp.csv");    
    [~,name1] = xlsread("bimodal_gene.xlsx", 'E c57');    
    nname = 'E c57';  
elseif nn == 4    
    [data,name]=xlsread("Embryonic_cast_fitmethod_nofixp.csv");  
    [~,name1] = xlsread("bimodal_gene.xlsx", 'E cast');  
    nname = 'E cast'; 
 end 
 [~, row, ~] =intersect(name(2:end,1), name1(2:end,1), 'stable'); 
 data_F1=data(row,:);  name_F1=name(row+1,1);
 AICc_F1 = [data_F1(:,10) data_F1(:,24) ]; %  data(:,65) data(:,78)
 AICc_F11 = min(AICc_F1,[],2);

 data_tw = data_F1(AICc_F1(:,1) == AICc_F11, :);   name_tw = name_F1(AICc_F1(:,1) == AICc_F11, :); 
 data_c =  data_F1(AICc_F1(:,2) == AICc_F11, :);   name_c = name_F1(AICc_F1(:,2) == AICc_F11, :);  
 clear data name name1 row AICc_F11 AICc_F1
  if nn == 1
       [data,name]=xlsread("Fibroblasts_c57_distribution.csv"); 
       yy = [0 5];
       yy1 = [2 23];
  elseif nn ==2
       [data,name]=xlsread("Fibroblasts_cast_distribution.csv");
        yy = [0 4];
        yy1 = [2 18.7];
  elseif nn == 3
       [data,name]=xlsread("Embryonic_c57_distribution.csv"); 
        yy = [0 5];
        yy1 = [1 16.5];
  elseif  nn == 4
       [data,name]=xlsread("Embryonic_cast_distribution.csv"); 
        yy = [0 6.8];
        yy1 = [1 14.5];
  end
[coom1,~, idx1] = intersect(name_tw, name, 'stable');
data1 = data(:, idx1);
[coom2, ~, idx2] = intersect(name_c, name, 'stable');
data2 = data(:, idx2);

%%% fig6(c)
 for alpha = 1 : length(data_tw(:,1))
    global xdata
    xdata=[];
    xdata = clearnan(data1(5:end, alpha));
    parameter = data_tw(alpha, 6:8);
    la = parameter(1); ga=parameter(2); v=parameter(3);  delta=1;
    model = 1;
    [mean_steady1, fano_steady1, ~, ~, ~  ] = caculate_moment(parameter, model);
    fano1(alpha) =fano_steady1;
    mean1(alpha) = mean_steady1;
end
for alpha = 1 : length(data_c(:,1))
    global xdata
    xdata=[];
    xdata = clearnan(data2(5:end, alpha));
    parameter = data_c(alpha, 18:22);
     kon1=parameter(1); kon2=parameter(2); q1=parameter(3); koff=parameter(4); kb=parameter(5); q2=1-q1;  delta=1;
     model = 2; 
     tend = 300;  
     [mean_steady2,fano_steady2, ~, ~, ~  ] = caculate_moment(parameter, model);
     fano2(alpha) =fano_steady2;
     mean2(alpha) = mean_steady2;
end
group1 = [
    repmat({'twostate'}, size(clearnan(fano1'),1), 1);   
    repmat({'crosstalk'}, size(clearnan(fano2'),1), 1)];
figure(1)
boxplot([clearnan(fano1'); clearnan(fano2')], group1)
ylim(yy1)
ylabel('fano')
title(nname)
 clearvars -except data1 data2 data_c data_tw mean1 mean2 nname yy

 t_mean1 = nan(length(data_tw(:,1)), 1);  
for alpha = 1 : length(data_tw(:,1))
    global xdata
    xdata=[];
    xdata = clearnan(data1(5:end, alpha));
    parameter = data_tw(alpha, 6:8);
    la = parameter(1); ga=parameter(2); v=parameter(3);  delta=1;
    M = 3*(length(xdata) + 10); 
    dt = 0.0001;
    model =1 ;
    tend = 200;%steady_time(parameter, model)+30;
    N = round(tend/dt) + 1;   
    P0 = zeros(N, M+1);
    P1 = zeros(N, M+1);
    P0(1,1) = 1;  % 初始状态：OFF, mRNA=0

    m_values = 0:M;
    time = (0:N-1) * dt;
    for t = 2:N
        % i=1 (mRNA=0)
        P0(t,1) = P0(t-1,1) + dt .* (ga.*P1(t-1,1) - la.*P0(t-1,1) + delta.*P0(t-1,2));
        P1(t,1) = P1(t-1,1) + dt .* (la.*P0(t-1,1) - (v+ga).*P1(t-1,1) + delta.*P1(t-1,2));
    
        i = 2:M;
        P0(t,i) = P0(t-1,i) + dt .* (ga.*P1(t-1,i) - ((i-1).*delta + la).*P0(t-1,i) + i.*delta.*P0(t-1,i+1));
        P1(t,i) = P1(t-1,i) + dt .* (la.*P0(t-1,i) - (v + (i-1).*delta + ga).*P1(t-1,i) + i.*delta.*P1(t-1,i+1) + v.*P1(t-1,i-1));

        P_total = P0(t,:) + P1(t,:);
        mean_steady1 = mean1(alpha);
        target_mean = 0.5 * mean_steady1;
        mean_mRNA1 = sum(m_values .* P_total);
        if  mean_mRNA1 >= target_mean && isnan(t_mean1(alpha))
            t_mean1(alpha) = time(t);
            break
        end
      
    end
    clear N  P0 P1 time m_values
end

t_mean2 = nan(length(data_c(:,1)), 1); 
for alpha = 1 : length(data_c(:,1))
     clearvars -except data1 data2 data_c data_tw mean1 mean2 nname t_mean2 t_mean1 yy alpha P_total1
     global xdata
     xdata=[];
     xdata = clearnan(data2(5:end, alpha)); 
     parameter = data_c(alpha, 18:22);
     kon1=parameter(1); kon2=parameter(2); q1=parameter(3); koff=parameter(4); kb=parameter(5); q2=1-q1;  delta=1;
     tend = 300;  
     mean_steady2 = mean2(alpha);
     target_mean = 0.5 * mean_steady2;
   
     M = 5*(length(xdata) + 20); 
     I = speye(M);
     Z = sparse(M, M);
     n = (0:M-1).';
     Adeg = spdiags(-delta * n, 0, M, M);
     if M > 1
         Adeg = Adeg + sparse((1:M-1).', (2:M).', ...
                            delta * (1:M-1).', M, M);
     end
     production_diagonal = [  -kb * ones(M-1, 1);  0];
     Aprod = spdiags(production_diagonal, 0, M, M);
     if M > 1
        Aprod = Aprod + sparse((2:M).', (1:M-1).', ...
                              kb * ones(M-1, 1), M, M);
     end
     A11 = Adeg - kon1 * I;
     A22 = Adeg - kon2 * I;
     A33 = Adeg + Aprod - koff * I;
     A = [     A11,        Z, q1 * koff * I; ...
                 Z,      A22, q2 * koff * I; ...
          kon1 * I, kon2 * I, A33];

     P01 = zeros(3 * M, 1);
     P01(1) = q1;
     P01(M + 1) = q2;
     rate_reference = max([kon1, kon2, koff, kb, delta, 1]);
     t_min = max(1e-3 / rate_reference, 1e-14);
     t_log_end = min(1, tend);
     if t_log_end > t_min
         t_log = logspace( log10(t_min),  log10(t_log_end), 100000);
     else
         t_log = t_log_end;
     end
     if tend > t_log_end
         t_linear = linspace(t_log_end, tend, 100000);
     else
         t_linear = [];
     end
     t_check = unique([0, t_log, t_linear]);

    options = odeset( ...
    'RelTol', 1e-7, ...
    'AbsTol', 1e-12, ...
    'Jacobian', A, ...
    'NonNegative', 1:(3*M));

     p_current = P01;
     t_previous = 0;
     chunk_size = 50;
     stop_flag = false;
     solver_failed = false;
     max_tail_probability = 0;
   % mean_mRNA2 = zeros(length(t_check),1);  fano_mRNA2 = zeros(length(t_check),1); ttt = zeros(length(t_check),1);
     for chunk_start = 2:chunk_size: length(t_check)
        chunk_end = min(  chunk_start + chunk_size - 1,  length(t_check));
        current_indices = chunk_start:chunk_end;
        current_tspan = [  t_previous, t_check(current_indices)];
        try
            [t_segment, Y_segment] = ode15s(  @(t, p) A * p,  current_tspan,  p_current,  options);
        catch ME
            solver_failed = true;
            break; 
        end
      for time_index = 2:length(t_segment)
            current_time = t_segment(time_index);
            p = Y_segment(time_index, :).';
            p(p < 0) = 0;
        P1_current = p(1:M);
        P2_current = p(M+1:2*M);
        P3_current = p(2*M+1:3*M);
        P_total1 =  P1_current +  P2_current +  P3_current;
        number_tail_states = min(5, M);
        tail_probability = sum( P_total1(M-number_tail_states+1:M));
        max_tail_probability = max( max_tail_probability,  tail_probability);

        mean_mRNA2_current = n.' * P_total1;
        if isnan(t_mean2(alpha)) &&   mean_mRNA2_current >= target_mean
             t_mean2(alpha) = current_time;
             p_current = p;
             stop_flag = true;
             break;
        end
      end
     if solver_failed
            break;
     end
    if stop_flag
        break;
    end
     p_current = Y_segment(end, :).';
     p_current(p_current < 0) = 0;
     t_previous = t_segment(end);
 end
 if solver_failed
        continue;
 end

 alpha
end

group2 = [
    repmat({'twostate'}, size(clearnan(t_mean1),1), 1);   
    repmat({'crosstalk'}, size(clearnan(t_mean2),1), 1)];
figure(2)
boxplot([clearnan(t_mean1); clearnan(t_mean2)], group2)
ylim(yy)
ylabel('Response time')
title(nname)






