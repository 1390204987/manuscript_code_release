%2024/12/2 by zn calculate PC1 explained varience in crease with time(under mTDRanalysis)
function topEV = gettopEV(Xb_reshape,Bp,startlen,startpoint)
    [n,T] = size(Bp);
    Yp_initial = Bp(:,1:1+startlen).*Xb_reshape; %Yp_initial = (n,startlen,Nmax); 
    Yp_initial = reshape(Yp_initial,n,[]);
    [coeff,~,var] = pca(Yp_initial');

    nchunk = T-startlen+1-startpoint;
    topEV = (var(1))/sum(var);
    
    % Incremental fitting
    for j = 1:nchunk
        ibegin = 1;
        iend = ibegin+startlen+j-1;
        Yp_current = Bp(:,ibegin:iend).*Xb_reshape; %Yp_initial = (n,startlen,Nmax);
        Yp_current = reshape(Yp_current,n,[]);
        [coeff,~,var] = pca(Yp_current');
        topEV = [topEV;(var(1))/sum(var)];
    end

end