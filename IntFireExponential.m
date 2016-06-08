% Init Neuron Type values
%C(pF) gL(ns) EL(mV) VT (mV) ?T (mV) a(nS) ?w(ms) b(pA) Vr(mV)
type = 'IB';
if type == 'RS'
   C =  200E-12;   gL =10E-9; EL=-70E-3; VT =-50E-3; deltaT = 2E-3;
   a = 2E-9;       Tw=30E-3;  b=0;      Vr = -58E-3;
elseif type == 'IB'
     C =  130E-12;   gL =18E-9; EL=-58E-3; VT =-50E-3; deltaT = 2E-3;
   a = 4E-9;       Tw=150E-3;  b=120E-12;      Vr = -50E-3;
elseif type == 'CH'
   C =  200E-12;   gL =10E-9; EL=-58E-3; VT =-50E-3; deltaT = 2E-3;
   a = 2E-9;       Tw=120E-3;  b=100E-12;      Vr = -46E-3;
end
Vpeak=0;

N = 4; % No of neurons
h = 0.1E-3; % time step
Tmax = 0.5; %max time
M = ceil(Tmax/h); % max time index i.e no. of columns (no. of delta_T time intervals that amount to T)


I = zeros(N,M);
for i = 1:N
    %I(i,1:M) = 400^-12 ; %(1 + alpha*i)*Ic;
     for k=1:M
         I(i,k) = (1-heaviside(k-M/2))*i*20*10^-12 ;
     end
    % I(i,1:M) = 600*10^-12 ;
end

y = zeros(N,M);
U = zeros(N,M);
% init condition
y(:,1) = Vr;
U(:,1) = 0;


F = @(t,V,U,I) (-gL*(V-EL)+gL*deltaT*exp((V-VT)/deltaT)-U+I)/C ; 
G = @(t,V,U) (a*(V-EL)-U)/Tw; 

%runge kutta implementation 4th order
for i = 1:M-1 
    k1 = F(i,       y(:,i),          U(:,i),          (I(:,i)+I(:,i+1))/2 );
    L1 = G(i,       y(:,i),          U(:,i));
    k2 = F(i+0.5*h, y(:,i)+0.5*h*k1, U(:,i)+0.5*h*L1, (I(:,i)+I(:,i+1))/2 );
    L2 = G(i+0.5*h, y(:,i)+0.5*h*k1, U(:,i)+0.5*h*L1);
    k3 = F(i+0.5*h, y(:,i)+0.5*h*k2, U(:,i)+0.5*h*L2, (I(:,i)+I(:,i+1))/2 );
    L3 = G(i+0.5*h, y(:,i)+0.5*h*k2, U(:,i)+0.5*h*L2);
    k4 = F(i+h,     y(:,i)+h*k3,     U(:,i)+h*L3,     (I(:,i)+I(:,i+1))/2 );
    L4 = G(i+h,     y(:,i)+h*k3,     U(:,i)+h*L3);
    
    y(:,i+1) = y(:,i) + (1/6)*(k1+2*(k2+k3)+k4)*h ;
    U(:,i+1) = U(:,i) + (1/6)*(L1+2*(L2+L3)+L4)*h ;
    
    y((y(:,i+1)>=0),i+1) = Vr;
    U((y(:,i+1)>=0),i+1) = U((y(:,i+1)>=0),i+1) + b;
end

t = h:h:Tmax;
% plot both on one plot only
figure(1)
subplot(2,1,1)
for i=1:N
    plot(t,I(i,:))
    hold on
end
title('current')
xlabel('Time (in s)')
ylabel('Current in A')
hold off
subplot(2,1,2)
for i=1:N
    figure(1)
    plot(t,y(i,:))
    hold on
end
hold off
title('Spiking Neuron')
xlabel('Time (in s)')
ylabel('Membrane Potential (in Volts)')