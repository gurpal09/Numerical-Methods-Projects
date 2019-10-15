% Newton Raphson Method
 clear all
 close all
 clc
 % Change here for different functions
 f=@(x) cos(x)-3*x+1
 %this is the derivative of the above function
 df=@(x) -sin(x)-3
 % Change lower limit 'a' and upper limit 'b'
 a=0; b=1;
 x=a;
 for i=1:1:100
 x1=x-(f(x)/df(x));
 x=x1;
 end
 sol=x;
 fprintf('Approximate Root is %.15f',sol)
 a=0;b=1;
 x=a;
 er(5)=0;
 for i=1:1:5
 x1=x-(f(x)/df(x));
 x=x1;
 er(i)=x1-sol;
 end
 plot(er)
 xlabel('Number of iterations')
 ylabel('Error')
 title('Error Vs. Number of iterations')