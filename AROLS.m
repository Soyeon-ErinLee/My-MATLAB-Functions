function [beta,Sig2,Cov]=AROLS(y,p)
T=size(y,1);
Y=y(p+1:T,1);
if p>0
    x=zeros(T-p,p);
    for j=1:p
        x(:,j)=y(p+1-j:T-j,1);
%        x=[X,y(p+1-j,T-j,1)];
    end
    X=[ones(T-p,1),x];
elseif p==0
    X=ones(T-p,1);
end
    beta=(X'*X)\(X'*Y); %inv쓰면 속도 느려지니까 동일한 backlash 사용
    U=Y-X*beta;
    Sig2=(U'*U)/(T-p-1);
    Cov=Sig2*inv(X'*X);
end