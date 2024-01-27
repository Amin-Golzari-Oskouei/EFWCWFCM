function j_fun = object_fun(N,d,k,Cluster_elem,landa,M,fuzzy_degree,W,z,p,X,gama)
    for j=1:k
        distance(j,:,:) = (1-exp((-1.*repmat(landa,N,1)).*((X-repmat(M(j,:),N,1)).^2)));
        WBETA = transpose(z(j,:));
        WBETA(WBETA==inf)=0;
        dNK(:,j) = reshape(distance(j,:,:),[N,d]) * WBETA * W(1,j)^p ;
    end
    j_fun = sum(sum(dNK .* transpose(Cluster_elem.^fuzzy_degree))) + (gama^-1 * sum(sum((z) .* log2(z))));
end

