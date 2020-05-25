function [Wald,pval]=GCtest(y1,y2,p)
% Testing for the null hypothesis that y2 does not Granger cause y1
% Wald is the value of the test stat
% pval is the p-value
% B'==[v,A1,...,Ap]
% v가 k x 1, A들이 k x k이기에 B'는 k x 1+kp matrix.
% once you vectorize B, B'의 each row가 하나의 column 아래 착착 쌓임(stack).
% 즉 K*(1+k*p) x 1 "by one vector"가 되는 것.
% Cmat shows the position which coefficients of are tested for zero in B' one by one
% 선형모형방법론에서 beta에 곱해지는 Cmat와 비슷한 역할인 듯.
% kxk matrix인 A 쪼개서
% | k1xk1   k1xk2 |  <- 0으로 만들고 싶은 부분
% | k2xk1   k2xk2 |
% s번째 row, A(r)에 있는 애 test하려면 
% (s-1)(1+K*p) + (r-1)*K+1 + k1 + q(반복) 번째 애를 를 1로 만들고 test

k1=size(y1,2);  
k2=size(y2,2);
K=k1+k2;
y=[y1,y2];
[B,SIG,~,ZZ]=olsvar(y,p);
alpha=reshape(B,K*(1+K*p),1); % alpha=vec(B)
SIGa = kron(SIG,inv(ZZ)); %covariance of beta_hat=vec(B)y
C=[];
Cmat=zeros(K,1+K*p);
for s=1:k1   % row info
    for r=1:p
        for q=1:k2   % column info
            Ctemp=zeros(1,(1+K*p)*K);
            Ctemp(1,(s-1)*(1+K*p)+(r-1)*K+1+k1+q)=1;
            Cmat=Cmat+reshape(Ctemp,1+K*p,K)';
            C=[C;Ctemp];
        end
    end
end

Wald=(C*alpha)'*inv(C*SIGa*C')*(C*alpha);
pval=1-cdf('chi2',Wald,size(C,1));

disp('H0: y2 does not Granger cause y1')
disp('')
disp(['The p-value is' num2str(pval)])
disp('')
pos=repmat('*',K,1+K*p); %처음엔 *로 이루어진 matrix 만들기
pos((Cmat==1))=['0']; % Cmat==1인 곳을 0으로 바꿈
disp(['Imposed zero restrictions are of the form'])
disp(['[v,A1,...,Ap] = '])
disp(pos)
