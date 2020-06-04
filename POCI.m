function test=POCI(y,x,m,k)
% Phillips-Ouliaris cointegration test
% y is a T*1 vector of observations
% x is a T*q matrix of observations
% m specifies the deterministic trend function
%   m=1: a constant
%   m=2: a constant and a linear time trend
% k is the number of lags of dy(t)

[T,q]=size(x);
switch m
    case 1
        z=ones(T,1);
        cv=[-3.90 -3.34 -3.05 % q=1
            -4.30 -3.74 -3.45 % q=2
            -4.65 -4.10 3.81]; % q=3
    case 2
        z=[ones(T,1),(1:T)'];
        cv=[-4.33 -3.78 -3.50 % q=1
            -4.67 -4.12 -3.83 % q=2
            -4.97 -4.43 -4.14]; % q=3
end

Z=[z,x];
u=y-Z*((Z'*Z)\(Z'*y));
Du=zeros(T-k-1,k); % since AR order is k+1
for idn=1:k
    Du(:,idn)=u(k+2-idn:T-idn,1)-u(k+1-idn:T-idn-1,1); % Dependent variable starts from K+2
end

if k==0
    uh=u(2:T,1);
    ul=u(1:T-1,1);
else
    uh=u(k+2:T)-Du*((Du'*Du)\(Du'*u(k+2:T)));
    ul=u(k+1:T-1)-Du*((Du'*Du)\(Du'*u(k+1:T-1)));
end

ahat=(ul'*ul)\(ul'*uh);
s2=((uh-ahat*ul)'*(uh-ahat*ul))/(T-2*k-2-m-q)
test=(ahat-1)/sqrt(s2/(ul'*ul));


disp('**********************************************')
disp('Reject the null of No CI if test < cv')
disp('critical values for size 1%, 5%, and 10% are')
disp(num2str(cv(q,:)))
disp('**********************************************')

end