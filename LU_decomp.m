function LU_decomp=LU_decomp(A)
    n=size(A,1);
    P=eye(n);
    L=eye(n);
    U=A;
    for i=1:n-1%第i行消元
        %寻找主元
        [~,k]=max(abs(U(i:n,i)));
        k=k+i-1;
        %交换行
        P=row_swap(P,i,k);
        L=row_swap(L,i,k);
        L=col_swap(L,i,k);
        U=row_swap(U,i,k);
        %消元
        for m=i+1:n
            L(m,i)=U(m,i)/U(i,i);
            U=row_addition(U,i,-L(m,i),m);
        end
    end
    LU_decomp.P=P;
    LU_decomp.L=L;
    LU_decomp.U=U;
end