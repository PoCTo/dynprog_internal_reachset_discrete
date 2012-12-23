function UV = A_svd(A)
    [U,~,V] = svd(A);
    UV = U*V;