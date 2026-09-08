clc
clear all
close all

%%% fig6a
nn = 1; %% 1: F c57; 2: F cast; 3: E c57; 4: E c57
kk = 1; %% 1: bimodal; 2: unimodal
if nn == 1
    [data,name]=xlsread("Fibroblasts_c57_fitmethod_nofixp.csv");
    [~,name1] = xlsread("expression of bimodal genes.xlsx", 'F c57');
elseif nn == 2
    [data,name]=xlsread("Fibroblasts_cast_fitmethod_nofixp.csv");
    [~,name1] = xlsread("expression of bimodal genes.xlsx", 'F cast');
elseif nn == 3
    [data,name]=xlsread("Embryonic_c57_fitmethod_nofixp.csv");
    [~,name1] = xlsread("expression of bimodal genes.xlsx", 'E c57');
elseif nn == 4
    [data,name]=xlsread("Embryonic_cast_fitmethod_nofixp.csv");
    [~,name1] = xlsread("expression of bimodal genes.xlsx", 'E cast');
end
[~, row, ~] =intersect(name(2:end,1), name1(2:end,1), 'stable');
if kk == 1
    data_F1=data(row,:);  name_F1=name(row+1,1);
elseif kk == 2
    data(row,:) = []; name(row+1,:) =[];  
    data_F1 = data;  name_F1 = name;
end
AICc_F1 = [data_F1(:,10) data_F1(:,24) ]; %  data(:,65) data(:,78)
AICc_F11 = min(AICc_F1,[],2);
clear data name row  name1 

para_F1=[ data_F1(:, 6:8), data_F1(:, 18:22) ];
mean = [data_F1(AICc_F1(:,1) == AICc_F11,1);  data_F1(AICc_F1(:,2) == AICc_F11,1)];
bf=[ para_F1(AICc_F1(:,1) == AICc_F11 , 1 );  1./( para_F1(AICc_F1(:,2) == AICc_F11,6)./para_F1(AICc_F1(:,2) == AICc_F11,4) + (1-para_F1(AICc_F1(:,2) == AICc_F11,6))./para_F1(AICc_F1(:,2) == AICc_F11,5) ) ];
bz=[ para_F1(AICc_F1(:,1) == AICc_F11,3)./para_F1(AICc_F1(:,1) == AICc_F11,2);  para_F1(AICc_F1(:,2) == AICc_F11,8)./para_F1(AICc_F1(:,2) == AICc_F11,7) ];
k1=length(data_F1(AICc_F1(:,1) == AICc_F11,1));
k2=length(data_F1);
figure(1)
h1 = scatter(log10(mean(1:k1)), log10(bf(1:k1)),  30, 'red',  'Marker', 'o');
hold on
h2 = scatter(log10(mean(k1+1:k2)), log10(bf(k1+1:k2)),  30, 'blue',  'Marker', 'o');
hold off
ylim([-3 1.1])
xlabel('log10(mean)')
ylabel('log10(bf)')
legend('two-state', 'cross-talk pathway')
figure(2)
h1 = scatter(log10(mean(1:k1)), log10(bz(1:k1)),  30, 'red',  'Marker', 'o');
hold on
h2 = scatter(log10(mean(k1+1:k2)), log10(bz(k1+1:k2)),  30, 'blue',  'Marker', 'o');
hold off
ylim([-0.3 3.2])
xlabel('log10(mean)')
ylabel('log10(bz)')
legend('telegraph', 'cross-talk pathway')


%%% figs7
clc
clear all
[data,~]=xlsread("Fibroblasts_c57_fitmethod_nofixp.csv");
data_F1=data; 
clear data name row type AICc1 AICc

[data,~]=xlsread("Fibroblasts_cast_fitmethod_nofixp.csv");
data_F2=data; 
clear data name row type AICc1 AICc

[data,~]=xlsread("Embryonic_c57_fitmethod_nofixp.csv");
data_E1=data; 
clear data name row type AICc1 AICc

[data,~]=xlsread("Embryonic_cast_fitmethod_nofixp.csv");
data_E2=data; 
clear data name row type AICc1 AICc


for i=1:4 % 6
    F1(i)=length(find(AICc_F1(:,i)==AICc_F11))./length(AICc_F11);
    F2(i)=length(find(AICc_F2(:,i)==AICc_F21))./length(AICc_F21);
    E1(i)=length(find(AICc_E1(:,i)==AICc_E11))./length(AICc_E11);
    E2(i)=length(find(AICc_E2(:,i)==AICc_E21))./length(AICc_E21);
end


para_F1=[ data_F1(:, 6:8), data_F1(:, 18:22), data_F1(:, 32:35), data_F1(:, 45:50)   ]; 
para_F2=[ data_F2(:, 6:8), data_F2(:, 18:22), data_F2(:, 32:35), data_F2(:, 45:50)  ]; 
para_E1=[ data_E1(:, 6:8), data_E1(:, 18:22), data_E1(:, 32:35), data_E1(:, 45:50)  ];
para_E2=[ data_E2(:, 6:8), data_E2(:, 18:22), data_E2(:, 32:35), data_E2(:, 45:50)  ];
mean_F1 = data_F1(: , 1); mean_F2 = data_F2(: , 1); mean_E1 = data_E1(: , 1); mean_E2 = data_E2(: , 1);
bf_F1=[para_F1(:,1) , 1./( para_F1(:,6)./para_F1(:,4) + (1-para_F1(:,6))./para_F1(:,5) ), 1./( 1./para_F1(:,9)+ 1./para_F1(:,10)), 1./(2.*para_F1(:,16)./para_F1(:,14)+(1-para_F1(:,16))./para_F1(:,13)) ];
bz_F1=[ para_F1(:,3)./para_F1(:,2),  para_F1(:,8)./para_F1(:,7),  para_F1(:,12)./para_F1(:,11),  para_F1(:,18)./para_F1(:,17) ];
bf_F2=[para_F2(:,1) , 1./( para_F2(:,6)./para_F2(:,4) + (1-para_F2(:,6))./para_F2(:,5) ), 1./( 1./para_F2(:,9)+ 1./para_F2(:,10)), 1./(2.*para_F2(:,16)./para_F2(:,14)+(1-para_F2(:,16))./para_F2(:,13)) ];
bz_F2=[ para_F2(:,3)./para_F2(:,2),  para_F2(:,8)./para_F2(:,7),  para_F2(:,12)./para_F2(:,11),  para_F2(:,18)./para_F2(:,17) ];
bf_E1=[para_E1(:,1) , 1./( para_E1(:,6)./para_E1(:,4) + (1-para_E1(:,6))./para_E1(:,5) ), 1./( 1./para_E1(:,9)+ 1./para_E1(:,10)), 1./(2.*para_E1(:,16)./para_E1(:,14)+(1-para_E1(:,16))./para_E1(:,13)) ];
bz_E1=[ para_E1(:,3)./para_E1(:,2),  para_E1(:,8)./para_E1(:,7),  para_E1(:,12)./para_E1(:,11),  para_E1(:,18)./para_E1(:,17) ];
bf_E2=[para_E2(:,1) , 1./( para_E2(:,6)./para_E2(:,4) + (1-para_E2(:,6))./para_E2(:,5) ), 1./( 1./para_E2(:,9)+ 1./para_E2(:,10)), 1./(2.*para_E2(:,16)./para_E2(:,14)+(1-para_E2(:,16))./para_E2(:,13)) ];
bz_E2=[ para_E2(:,3)./para_E2(:,2),  para_E2(:,8)./para_E2(:,7),  para_E2(:,12)./para_E2(:,11),  para_E2(:,18)./para_E2(:,17) ];


n=1;
figure(3)
h1 = scatter(log10(mean_F1), log10(bf_F1(:,n)),  40, 'red',  'Marker', 'o');
hold on
h2 = scatter(log10(mean_F2), log10(bf_F2(:,n)), 40, 'blue',  'Marker', '^');
hold on
h3 = scatter(log10(mean_E1), log10(bf_E1(:,n)),  40, 'green',  'Marker', '*');
hold on
h4 = scatter(log10(mean_E2), log10(bf_E2(:,n)), 40, 'cyan',  'Marker', 's');
hold off
ylim([-3 1.2])
title('fitted by the telegraph model')
figure(4)
h1 = scatter(log10(mean_F1), log10(bz_F1(:,n)),  40, 'red',  'Marker', 'o');
hold on
h2 = scatter(log10(mean_F2), log10(bz_F2(:,n)), 40, 'blue',  'Marker', '^');
hold on
h3 = scatter(log10(mean_E1), log10(bz_E1(:,n)),  40, 'green',  'Marker', '*');
hold on
h4 = scatter(log10(mean_E2), log10(bz_E2(:,n)), 40, 'cyan',  'Marker', 's');
hold off
ylim([-0.3 4.2])
title('fitted by the telegraph model')

n=2;
figure(5)
h1 = scatter(log10(mean_F1), log10(bf_F1(:,n)),  40, 'red',  'Marker', 'o');
hold on
h2 = scatter(log10(mean_F2), log10(bf_F2(:,n)), 40, 'blue',  'Marker', '^');
hold on
h3 = scatter(log10(mean_E1), log10(bf_E1(:,n)),  40, 'green',  'Marker', '*');
hold on
h4 = scatter(log10(mean_E2), log10(bf_E2(:,n)), 40, 'cyan',  'Marker', 's');
hold off
%ylim([-3 1.5])
title('fitted by the cross-talk pathway model')
figure(6)
h1 = scatter(log10(mean_F1), log10(bz_F1(:,n)),  40, 'red',  'Marker', 'o');
hold on
h2 = scatter(log10(mean_F2), log10(bz_F2(:,n)), 40, 'blue',  'Marker', '^');
hold on
h3 = scatter(log10(mean_E1), log10(bz_E1(:,n)),  40, 'green',  'Marker', '*');
hold on
h4 = scatter(log10(mean_E2), log10(bz_E2(:,n)), 40, 'cyan',  'Marker', 's');
hold off
%ylim([-0.3 3])
title('fitted by the cross-talk pathway model')

n=3;
figure(7)
h1 = scatter(log10(mean_F1), log10(bf_F1(:,n)),  40, 'red',  'Marker', 'o');
hold on
h2 = scatter(log10(mean_F2), log10(bf_F2(:,n)), 40, 'blue',  'Marker', '^');
hold on
h3 = scatter(log10(mean_E1), log10(bf_E1(:,n)),  40, 'green',  'Marker', '*');
hold on
h4 = scatter(log10(mean_E2), log10(bf_E2(:,n)), 40, 'cyan',  'Marker', 's');
hold off
%ylim([-3 1.5])
title('fitted by the three-state model')
figure(8)
h1 = scatter(log10(mean_F1), log10(bz_F1(:,n)),  40, 'red',  'Marker', 'o');
hold on
h2 = scatter(log10(mean_F2), log10(bz_F2(:,n)), 40, 'blue',  'Marker', '^');
hold on
h3 = scatter(log10(mean_E1), log10(bz_E1(:,n)),  40, 'green',  'Marker', '*');
hold on
h4 = scatter(log10(mean_E2), log10(bz_E2(:,n)), 40, 'cyan',  'Marker', 's');
hold off
%ylim([-0.3 3])
title('fitted by the three-state model')

n=4;
figure(9)
h1 = scatter(log10(mean_F1), log10(bf_F1(:,n)),  40, 'red',  'Marker', 'o');
hold on
h2 = scatter(log10(mean_F2), log10(bf_F2(:,n)), 40, 'blue',  'Marker', '^');
hold on
h3 = scatter(log10(mean_E1), log10(bf_E1(:,n)),  40, 'green',  'Marker', '*');
hold on
h4 = scatter(log10(mean_E2), log10(bf_E2(:,n)), 40, 'cyan',  'Marker', 's');
hold off
%ylim([-3 1.5])
title('fitted by the cross-talk three-state model')
figure(10)
h1 = scatter(log10(mean_F1), log10(bz_F1(:,n)),  40, 'red',  'Marker', 'o');
hold on
h2 = scatter(log10(mean_F2), log10(bz_F2(:,n)), 40, 'blue',  'Marker', '^');
hold on
h3 = scatter(log10(mean_E1), log10(bz_E1(:,n)),  40, 'green',  'Marker', '*');
hold on
h4 = scatter(log10(mean_E2), log10(bz_E2(:,n)), 40, 'cyan',  'Marker', 's');
hold off
%ylim([-0.3 3])
title('fitted by the cross-talk three-state model')


