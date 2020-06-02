% Unit Rate test
function test=ADFk(y,m,k)
% Augmented D-F Unit Root Test
% y is a column vector of observations
% m specifies the deterministic function(or trend), mu
% m=0: no deterministic function(or trend)
% m=1: a constant
% m=2: a constant and a linear trend
% k is the number of the lags of dy(t)
% AR order: k+1
T=size(y,1);
switch m
    case 1  % when m=1
        z=ones(T,1);
    case 2  % when m=2
        z=[ones(T,1),(1:T)'];
end
Dy=zeros(T-k-1,k); % since AR order is k+1
for idn=1:k
    Dy(:,idn)=y(k+2-idn:T-idn,1)-y(k+1-idn:T-idn-1,1); % Dependent variable starts from K+2
end

switch m
    case 0 % when m=0, no deterministic trend
        if k==0  % when k=0, just take 2nd to last regressor
            yh=y(2:T);
            yl=y(1:T-1);
        else
            Z=Dy;
            yh=y(k+2:T)-Z*((Z'*Z)\(Z'*y(k+2:T)));
            yl=y(k+1:T-1)-Z*((Z'*Z)\(Z'*y(k+1:T-1)));
        end
    otherwise
        if k==0
            Z=z(2:T,:);
        else
            Z=[z(k+2:T,:),Dy];
        end
        yh=y(k+2:T)-Z*((Z'*Z)\(Z'*y(k+2:T)));
        yl=y(k+1:T-1)-Z*((Z'*Z)\(Z'*y(k+1:T-1)));
end
ahat=(yl'*yl)\(yl'*yh);
s2=((yh-ahat*yl)'*(yh-ahat*yl))/(T-2*k-2-m);
test=(ahat-1)/sqrt(s2/(yl'*yl));

cv=[-2.58 -2.33 -1.95; % for m=0
    -3.43 -3.12 -2.86; % for m=1
    -3.96 -3.66 -3.41]; % for m=2

disp('**********************************************')
disp('Reject the null of UR if test < cv')
disp('critical values for size 1%, 5%, and 10% are')
disp(num2str(cv(m+1,:)))
disp('**********************************************')
end