
function [lam, gamma, nu]=Moment_of_method(xdata)
S=length(xdata); 
e1=0;e2=0;e3=0;
for i=1:1:S;
    e1=(i-1).*xdata(i)+e1;
    e2=((i-1).^2).*xdata(i)-(i-1).*xdata(i)+e2;
    e3=((i-1).^3).*xdata(i)-3*((i-1).^2).*xdata(i)+2*(i-1).*xdata(i)+e3;
end
    e1;  %the first moment 
    e2;  %the second moment
    e3;  %the third moment
    r1=e1;r2=e2/e1;r3=e3/e2;
    lam=2*r1*(r3-r2)/(r1*r2-2*r1*r3+r2*r3);
    gamma=2*(r2-r1)*(r1-r3)*(r3-r2)/((r1*r2-2*r1*r3+r2*r3)*(r1-2*r2+r3));
    nu=(-r1*r2+2*r1*r3-r2*r3)/(r1-2*r2+r3);
end