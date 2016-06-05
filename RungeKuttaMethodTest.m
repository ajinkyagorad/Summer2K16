xx%Implementation of Runge Kutta Method
%---------------------------------
% y(i+1) = y(i) +(Weighted Avg of slopes*h)
% K1*W1+K2*W2+K3*W3+...+Kv*Wv : v-stage method
% K1 = h*f(xi,yi)
% K2 = h*f(xi+c2*h,yi+a21*K1)
%For  2 stage Method
% y(i+1) = y(i)+ h f(xi+0.5h,yi+0.5*h*f(xi,yi))
%-----------------------------------
clear;
h = 0.01;
y = 2; % initial value
y_old = y; 
A=0;
for n = 0:1000;
    y = y_old+h*f(n*h+0.5*h,y_old*h*f(n*h,y_old));
    A=[A y];
    y_old=y;
end
plot(A);
    