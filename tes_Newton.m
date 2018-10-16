close 
clear
clc

%Choose one equation from equation's option and set parameters
equation=3;
epsilon=10^-5;
max_iter=30;

%Equation's Option
if(equation==1)   
   X = sym('x%d', [1 2]); 
   fun={X(1)^2+X(2)^2-4, X(1)*X(2)-1};
   lb=[-3 -3];ub=[3 3];
   axis=[-3 3];
elseif(equation==2)
   X = sym('x%d', [1 2]); 
   fun={exp(X(1)-X(2))-sin(X(1)+X(2)), (X(1)^2)*(X(2)^2)-cos(X(1)+X(2))};
   lb=[-10 -10];ub=[10 10];
   axis=[-10 10]; 
elseif(equation==3)
   X = sym('x%d', [1 2]); 
   fun={0.5*sin(X(1)*X(2)) - 0.25*X(2)/pi - 0.5*X(1), (1-0.25/pi)*(exp(2*X(1))-exp(1)) + exp(1)*X(2)/pi - 2*exp(1)*X(1)};
   lb=[-0.5 -20];ub=[3 5];
   axis=[-0.5 3 -20 5];
else
   X = sym('x%d', [1 6]); 
   fun={X(1)+(X(2)^2*X(4)*X(6))/4+0.75,... 
       X(2)+0.405*exp(1+X(1)*X(2))-1.405,...
       X(3)-X(4)*X(6)/2+1.5, ... 
       X(4)-0.605*exp(1-X(3)^2-0.395),...
       X(5)-X(2)*X(6)/2+1.5,...
       X(6)-X(1)*X(5)};
   lb=[-5 -5 -5 -5 -5 -5];ub=[5 5 5 5 5 5];
end

%Set number of point for spread point
Nvar=length(X);
spreadPoint=100;
r=rand(spreadPoint,Nvar);

%Spreads point to find initialization point that give root of equation
for i=1:spreadPoint
    for j=1:length(lb)
        x(i,j)=lb(j)+(ub(j)-lb(j))*r(i,j);
    end
end

sol=[];
k=[];
for i=1:spreadPoint 
    %Calculate the solution using Newton method
    [ temp, ~, n_iter] = NewtonMethod( fun, X, x(i,:), epsilon, max_iter);
    
    %Filter the solution if have the same or similar solution 
    stat=0;
    for j=1:Nvar
        if temp(j)>=lb(j) && temp(j)<=ub(j)
            stat=stat+0;
        else
            stat=stat+1;
        end
    end
    
    %stop condition
    if (stat==0 && n_iter<max_iter)
        sol=[sol;temp];
        k=[k;n_iter];
    end
    i
end

[row,col]=size(sol);

c=zeros(row,1);

%filter the solution if out of the domain of search 
for i=1:row-1
    for j=i+1:row
        if (sum(abs(sol(i,:)-sol(j,:)))<=10^-3)
            c(j)=c(j)+1;
        else 
            c(j)=c(j)+0;
        end
    end
end

sol=sol(c==0,:);

%Visualization of the equation
if equation ~=4
    ezplot(fun{1}, axis);
    hold on;
    ezplot(fun{2}, axis);
    hold on;
    plot(sol(:,1),sol(:,2),'*');
    title('Newton')
    xlabel('X1')
    ylabel('X2')
end