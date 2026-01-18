function A=calculate_AICc(Smin_min, k, n)
A=2.*k+2.*Smin_min+2.*k.*(k+1)./(n-k-1);
end