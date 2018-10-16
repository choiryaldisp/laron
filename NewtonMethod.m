function [ sol, x, k ] = NewtonMethod( func, X, x, epsilon, max_iter)

    Nfun=length(func);
    Nvar=length(X);
    TOL=epsilon;
    k=1;
    stop=false;

    %Jacobian Matrix
    for i=1:Nfun
        for j=1:Nvar
            diffMatrix{i,j}=diff(func{i},X(j));
        end
    end
    
    
    while(stop~=true && k<max_iter)
        %subtitute the X* into Jacobian matrix and function
        for i=1:Nfun
            for j=1:Nvar
                jacobi(i,j)=double(subs(diffMatrix{i,j},X,x(k,:)));
            end
            Fx(i)=double(subs(func{i},X,x(k,:)));
        end
    
        %calculate the new X*
        y=(-1*(jacobi))\Fx';
        x(k+1,:)=x(k,:)+y';
        
        %stop condition
        if (sum(abs(x(k+1,:)-x(k,:)))<=TOL)
            stop=true;
        end
        k=k+1;
    end
    
    sol=x(end,:);
end

