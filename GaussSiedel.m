function [ sol, all_sol, N_iter ] = GaussSiedel( A, bi, x0, iter_max, epsilon )

    %term of stop
    treshold=epsilon;

    %number of equation
    n=length(bi);

    %iteration
    N=iter_max;

    %solution
    x=x0;
    k=2;
    stop=false;

    while stop~=true && k<=N
        x(k,:)=x(k-1,:);
        for i=1:n
            cal=0;
            for j=1:n
                if j~=i
                cal=cal+A(i,j)*x(k,j);
                end
            end
            x(k,i)=(-cal+bi(i))*1/A(i,i);
        end
        if abs(sum(x(k,:)-x(k-1,:)))< treshold
           stop=true;
        end
        k=k+1;
    end

    sol=x(end,:);
    all_sol=x;
    N_iter=k-1;
end

