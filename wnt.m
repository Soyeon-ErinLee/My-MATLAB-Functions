function [BoxPierce,pvalue_BP,LjungBox,pvalue_LB]=wnt(X,taumax)
[T,k]=size(X);
if k>1
    disp('Warning! X should be a column vector')
    return
else
    acorr=autocorr(X,taumax);
    BoxPierce=T*(acorr'*acorr);
    pvalue_BP=(1-cdf('Chi2',BoxPierce,taumax))*100;
    s=(T-1:-1:T-taumax)';
    LjungBox=T*(T+2)*((acorr./s)'*acorr); %./ = element by element division
    pvalue_LB=(1-cdf('Chi2',LjungBox,taumax))*100;
end
end