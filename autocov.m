function cov=autocov(X,taumax)
%compute the sample autocovariances up to order taumax
%X is a column vector of observations
[T,k]=size(X);
if k>1
    disp('Warning! X should be a column vector')
    return
else
   Xd=X-ones(T,1)*mean(X);
   cov=zeros(taumax+1,1);
   for n=0:taumax
       cov(n+1,1)=Xd(n+1:T,1)'*Xd(1:T-n,1)/T;
   end
end 
end