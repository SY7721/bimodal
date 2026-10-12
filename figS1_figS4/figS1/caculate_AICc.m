function A=caculate_AICc(Smin_min,k,CN)
         A=2*k+2*Smin_min+2*k*(k+1)./(CN-k-1);
end