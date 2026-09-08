clc
clear all

[~, name1] = xlsread('bimodal_gene.xlsx','F c57');
name_F_c57 = name1(2:end,1);
clear name1
[~, name1] = xlsread('bimodal_gene.xlsx','F cast');
name_F_cast = name1(2:end,1);
clear name1
[~, name1] = xlsread('bimodal_gene.xlsx','E c57');
name_E_c57 = name1(2:end,1);
clear name1
[~, name1] = xlsread('bimodal_gene.xlsx','E cast');
name_E_cast = name1(2:end,1);
clear name1

union_F = union(name_F_c57, name_F_cast, 'stable');
union_E = union(name_E_c57, name_E_cast, 'stable');

unique_F = setdiff(union_F, union_E, 'stable'); %% unique_F 在 union_F 中但不在 union_E 中的元素
unique_E = setdiff(union_E, union_F, 'stable');

clc
clear all
[data, name] = xlsread("bimodal_gene.xlsx", 'E cast');
[data1, name1] = xlsread("bimodal_nop.xlsx");
name_nop = name(2:end,1);
name_p = name1(2:end,4);

unique_F = setdiff(name_nop, name_p, 'stable'); 
[common,~,~] = intersect(name_nop, name_p, 'stable');


% [~,name] = xlsread('gene name.xlsx','F c57');
% [~, name1] = xlsread('gene name.xlsx','E c57');
% 
% L1= sum(~cellfun(@isempty, name(:,1)));
% L2= sum(~cellfun(@isempty, name(:,3)));
% name_F = [name(2 : L1 ,1); name(2 : L2 ,3)];
% 
% L4= sum(~cellfun(@isempty, name1(:,1)));
% L5= sum(~cellfun(@isempty, name1(:,3)));
% name_E = [name1(2 : L4 ,1); name1(2 : L5 ,3)];
% 
% clearvars -except name_E name_F
% [common_gene, idx1, idx2] = intersect(name_F, name_E, 'stable');
% 
% name_F(idx1) = []; 
% name_E(idx2) = [];
