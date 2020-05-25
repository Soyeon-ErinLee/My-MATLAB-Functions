function corr=autocorr(X,taumax)
[~,k]=size(X);
if k>1
    disp('Warning! X should be a column vector')
    return
else
    acov=autocov(X,taumax);
    corr=acov(2:taumax+1,1)/acov(1,1);
end
end