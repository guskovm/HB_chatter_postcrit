function hb_chatter_postcrit
% HB_CHATTER_POSTCRIT.M - script accompanying the manuscript
% "Harmonic balance based simulation of chatter in turning including interrupted cut"
% (c) Mikhail GUSKOV, 2025, under LGPL-3.0 license
% PIMM laboratory, Arts & Métiers Institute of Technology, CNAM, CNRS UMR8006

global k z V nout n N Xh mu Ti Tdi 
N = 5; % retained harmonics number
n = 64; % time discretization
mu = 5; % delay terms number
t = (1:n)'*2*pi/n; % time variable array
Nh = 2*N+1; % number of harmonuc unknowns

z = .04; % damping ratio
V = sqrt(1+2*z);  % stability lobe low point frequency
D = 2*(pi-atan(1/V));  % stability lobe low point delay term
D00 = D; V00=V;

% Fourier transform matrix generation
T = ones(1,n); % X = T*x
Ti = 0.5*ones(n,1); % x = Ti*X
Tdi = zeros(n,1); % x' = Tdi*X
Tddi = zeros(n,1); % x'' = Tddi*X
for ii = 1:N
    T = [T; cos(ii*t'); sin(ii*t')];
    Ti = [Ti, cos(ii*t), sin(ii*t)];
    Tdi = [Tdi, -ii*sin(ii*t), ii*cos(ii*t)];
    Tddi = [Tddi, -ii^2*cos(ii*t), -ii^2*sin(ii*t)];
end
T = T*2/n;

% continuation loop variables initialization
SS = []; S = 0; nout = 20; io = 0; INC = 0;
ds = .003; dds = ds;
etol = min(1e-9,ds/10); ii = 1;
dsmax = ds; dsmin = etol*10; dsred = 2/3; dsinc = pi/exp(1);
itmax = 30; itmin = 3; itJ = 3;

% first point
V0 = V; k = ((V^2-1)^2+4*V^2*z^2)/(2*(V^2-1)); K = k; K0 = K;
X0 = zeros(2*N+1,1); X0(1) = -2*g(0); X = X0; X(2) = dds/2;

% MAIN LOOP 
CONT = 1; STOP = 0; Kmin=2*z*(1+z); Kmax=4*Kmin; 
XX = X0'; VV = V0; KK = K0; % results arrays
while CONT && ~STOP, CONT = V0 >0 & K < Kmax & norm(X) < 5 & ii < 1e4;     
    itit = 0; ee = 1;
    while (itit<itmax)&&(ee>etol)
        % current residual evaluation
        L = 1; L_ = L*0;
        for il = 1:N,
            L(2*il+(0:1),2*il+(0:1))=[-V^2*il^2+1, 2*V*z*il; -2*V*z*il, -V^2*il^2+1];
            L_(2*il+(0:1),2*il+(0:1))=[-2*V*il^2, 2*z*il; -2*z*il, -2*V*il^2];
        end
        D = D00*V/V00; 
        k = K; 
        ss = sx(X,D);
        g0 = g(Ti*X-ss);
        G = T*g0;
        R = [L*X + G; (X-X0).'*(X-X0) + (V-V0)^2 + (K-K0)^2 - dds^2; X(3)];

        % Newton=Raphson iteration step
        if mod(itit,itJ)==0, dd = 1e-10;
            dG = zeros(Nh);
            dXi = 1e-10;
            for ix = 1:Nh
                dX = X*0; dX(ix) = dXi; ssd = sx(X+dX,D);
                dG(:,ix) = T*(g(Ti*(X+dX)-ssd)-g0)/dXi;
            end
            J = [L + dG, L_*X + T*(g(Ti*X-ssd)-g0)/dd, G/k;
                 2*(X-X0)', 2*(V-V0), 2*(K-K0);
                 [0 0 1 zeros(1,2*N)]];
        end
        Y_ = [X; V; K]; dY = - J\R; if norm(dY)>(dds/2), dY = dY/10; end
        Y = Y_ + dY; ee = norm(R)+norm(Y-Y_)/norm(Y_)/10;
        X = Y(1:end-2); V = Y(end-1); K = Y(end);
        itit = itit+1;
        if ((itit >= itmax)||(dds>dsmax))&&~isempty(SS)&&(dds>dsmin)
            dds = dds*dsred; S = SS(end) + dds;
            if ii<=2, Xp = X; Vp = V; Kp = K+dds;
            else
                Xp = (XX(ii,:)+(XX(ii,:)-XX(ii-1,:))*dsred)';
                Vp = VV(ii)+(VV(ii)-VV(ii-1))*dsred;
                Kp = KK(ii)+(KK(ii)-KK(ii-1))*dsred;
            end
            itit = 0;
        end
    end
    if (itit <itmin)&&(dds<dsmax), dds = dds*dsinc; INC = 1; end
    if dds<dsmin, STOP = 1; end

    % time marching check
    if (ii>2&&ii<5)||mod(ii,nout)==0, io = io+1; iii(io) = ii+1; 
        Xh = X;
        ddeopt = ddeset('maxstep',diff(t(1:2))/2,'reltol',1e-8,'abstol',1e-5);
        sol = dde23(@ddefn, D*(1:mu), @ddehist,t',ddeopt);
        xt = interp1(sol.x,sol.y(1,:),t); xh = Ti*Xh;
        Eht(io) = rms(xh-xt)/rms(xh);
    end
    
    % increment, storing
    ii = ii+1;
    XX(ii,:) = X'; X0 = X;  VV = [VV V]; KK = [KK K]; SS = [SS S]; S=S+dds;
        
    % further prediction
    if ii<=2, Xp = X; Vp = V; Kp = K; Xp(2) = Xp(2)+dds;
    else
        Xp = (XX(ii,:)+(XX(ii,:)-XX(ii-1,:))*((1-INC)+dsinc*INC))';
        Vp = VV(ii)+(VV(ii)-VV(ii-1))*((1-INC)+dsinc*INC);
        Kp = KK(ii)+(KK(ii)-KK(ii-1))*((1-INC)+dsinc*INC);
    end; INC=0;
    X = Xp; V = Vp; K = Kp; V0 = VV(ii); K0 = KK(ii);
end

%% plots
figure(1); plot(KK,XX); 
legstr{1} = 'a_0' ; 
for il = 1:N
    legstr{2*il} = num2str(il, 'a_{%g}');
    legstr{2*il+1} =  num2str(il, 'b_{%g}');
end
grid on;  axis tight; xlabel('\kappa'); ylabel('x_i'), legend(legstr);
figure(2); plot(KK,VV); 
grid on; axis tight; xlabel('\kappa'); ylabel('\nu')
figure(3); pht = KK(iii);  xlab = '\kappa'; 
semilogy(pht,[Eht],'.-'); grid on; axis tight
xlabel(xlab); ylabel('\epsilon = RMS(x_h-x_t)/RMS(x_h)');

function ff = g(X)
% Cutting force [Dombovari'11 p332 eq(4)]
global k
H2 = -0.20015; H3 = 0.01946;
ff = k*((X + H2*X.^2 + H3*X.^3).*(X>-1)+(H2-H3-1).*(X<=-1));

function s = sx(X,D)
% Surface computation
global mu n Ti N 
t = (1:n)'*2*pi/n; md = 1; x = Ti*X; 
iic =1; s = 100*ones(size(x));
while ~isempty(iic) && md<=mu
    TimD = TiDf(md*D,t,N); xmD = TimD*X;
    iic = find((s+1)>(xmD+md));
    s(iic) = xmD(iic)+(md-1);
    md = md+1;
end

function TiD=TiDf(D,t,N)
% D-delayed inverse Fourier transform matrix
n=length(t);
TiD = zeros(n,2*N+1);
TiD(1:n,1)= 0.5*ones(n,1); 
for iv = 1:N
	TiD(1:n,2*iv+[0 1]) = [cos(iv*t-D), sin(iv*t-D)];
end

function ff = ddefn(t,y,yD) % y = [x, x'];
% DDE derivative function
global z V mu
s = min(yD(1,:) + (0:(mu-1)));
ff = [y(2);(- 2*z*V*y(2) - y(1) - g(y(1)-s))/V^2];

function y = ddehist(t) % y = [x, x'];
% DDE initial set
global Xh N
y = [Xh(1)/2; 0];
for ii = 1:N
    y(1) = y(1)+   [ cos(ii*t'), sin(ii*t')]*Xh(2*ii + (0:1));
    y(2) = y(2)+ii*[-sin(ii*t'), cos(ii*t')]*Xh(2*ii + (0:1));
end