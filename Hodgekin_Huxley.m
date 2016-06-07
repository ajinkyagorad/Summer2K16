% constants
C = 1E-6;     % F/cm2
ENa = 50E-3;  % V
EK = -77E-3; 
Et = -55E-3;
gNa = 120E-3; % S/cm2
gK = 36E-3 ;  % S/cm2
gt = 0.3E-3;  % S/cm2
I0 = 15E-6;   % A/cm2

%Timings
h = 0.01E-3;  % 0.01 ms
Tmax = 150E-3;% 150 ms
size = int16(Tmax/h);
M = size;
T = 30E-3;    % 30ms

% free variables
V = zeros(1,size);
I = zeros(1,size);
m_conc = zeros(1,size);
n_conc = zeros(1,size);
h_conc = zeros(1,size);
w = zeros(4,M); % w = [ m; n; h; V]
       

% input current
I(int16((2*T)/h) : int16((3*T)/h)) = I0 ; 

% definations of functions

alpha_n = @(V) 0.01*(V*E3+55)/(1-exp(-(V*E3+55)/10));
alpha_m = @(V) 0.1*(V*E3+40)/(1-exp(-(V*E3+40)/10));
alpha_h = @(V) 0.07*(1-exp(-0.05*(V*E3+65))) ;
beta_n = @(V)  0.125*exp(-(V*E3+65)/80);
beta_m = @(V)  4*exp(-0.0556*(V*E3+65));
beta_h = @(V)  1/(1+exp(-0.1*(V*E3+35)));
iNa_t = @(m, h, V) gNa*m^3 * h* (V-ENa);
iK_t  = @(n,V) gK*n^4 * (V-EK);
il_t  = @(V) gl*(V-El);

% diffferntial equations governing the spiking of neuron
dmdt = @(t, m, n, h, V) alpha_m(V)*(1-m)-beta_m(V)*m;
dndt = @(t, m, n, h, V) alpha_n(V)*(1-n)-beta_n(V)*n;
dhdt = @(t, m, n, h, V) alpha_h(V)*(1-h)-beta_h(V)*h;
dVdt = @(t, m, n, h, V, Iext) ( -iNa_t(m, h, V) - iK_t(n, V) - il_t(V) + Iext)/C;

f = @(t,y,I) [dmdt(t, y(1), y(2), y(3), y(4)), dndt(t, y(1), y(2), y(3), y(4)), dhdt(t, y(1), y(2), y(3), y(4)), dVdt(t, y(1), y(2), y(3), y(4), I)] 

% init conditions
f_init = @(t,y) [dmdt(t, y(1), y(2), y(3), y(4)), dndt(t, y(1), y(2), y(3), y(4)), dhdt(t, y(1), y(2), y(3), y(4)), dVdt(t, y(1), y(2), y(3), y(4), 0)] 
w_init = fsolve(f_init,[0 0 0 0]);
w(:,1) = w_init(:);
% order of computation
%V_=>alpha_nmh, beta_nmh,=>m,n,l=>iNa,K,l=> V

%runge kutta implementation
for i = 1:M-1 
    k1 = h*f(i    , w(:,i)        , I(i));
    k2 = h*f(i+h/2, w(:,i)+0.5*k1 , I(i));
    k3 = h*f(i+h/2, w(:,i)+0.5*k2 , I(i)); 
    k4 = h*f(i+h  , w(:,i)+k3     , I(i));
    w(:,i+1) = w(:,i) + (k1 + 2*k2 + 2*k3 + k4)/6;
end
    