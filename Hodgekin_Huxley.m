% constants
C = 1E-6; %uF/cm2
ENa = 50E-3; % V
EK = -77E-3; 
Et = -55E-3;
gNa = 120E-3; % S/cm2
gK = 36E-3 ; % S/cm2
gt = 0.3E-3; % S/cm2

%Timings
h = 0.1E-3;;
Tmax = 1;
size = Tmax/h;

% free variables
V = zeros(1,size);
I =  zeros(1,size);
iNa =  zeros(1,size);
iK = zeros(1,size);
il =  zeros(1,size);
m = zeros(1,size);
n =  zeros(1,size);
h = zeros(1,size);

% definations of functions
I0 = 1E-3;
Iext = I0*ones(size);
dVdt = @(iNa,iK,il,Iext) (-iNa-iK-il+Iext)/C;
dxdt = @(x,a,b) a*(1-x)-b*x;

alpha_n = @(V) 0.01*(V*E3+55)/(1-exp(-(V*E3+55)/10));
alpha_m = @(V) 0.1*(V*E3+40)/(1-exp(-(V*E3+40)/10));
alpha_h = @(V) 0.07*(1-exp(-0.05(V*E3+65)));
beta_n = @(V)  0.125*exp(-(V*E3+65)/80);
beta_m = @(V)  4*exp(-0.0556(V*E3+65));
beta_h = @(V)  1/(1+exp(-0.1(V*E3+35)));
iNa_t = @(m,h,V) gNa*m^3 * h* (V-ENa);
iK_t  = @(n,V) gK*n^4 * (V-EK);
il_t  = @(V) gl*(V-El);

% order of computation
%V_=>alpha_nmh, beta_nmh,=>m,n,l=>iNa,K,l=> V
for i=1:h:Tmax
    K