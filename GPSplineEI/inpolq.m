clear;clc;
% ================================================= %
% Settings
% ================================================= %
% cost function
fun = @(x)3*x(1)^2+x(2)^2+55*x(3)^2+2*x(4)^2+x(5)^2;
gra = @(x)2*[3*x(1),x(2),55*x(3),2*x(4),x(5)]';
% ======================================================== %
% Just try data from bfgs (last step 4>5)
% ======================================================== %
pt4 = [0.8561,4.6657,-0.0115,-2.0839,3.6095]';
pt5 = -[0.1822,0.0984,0.0093,0.3879,0.0761]';
[storetest,d] = cursamples(12,pt4,pt5,fun);

funaa = @(alphar)fun(pt4 + alphar*d);
alpharo = fminbnd(funaa,0,10);

[slopeu4,xxu4,yyu4] = sl(fun(pt4), 0, gra(pt4), d);
[slopeu5,xxu5,yyu5] = sl(fun(pt5), alpharo, gra(pt5), d);


% add some noise to the objective function
obj4p = fun(pt4)+randn(1);
obj5p = fun(pt5)+randn(1);
gra4p = gra(pt4)+randn(1);
gra5p = gra(pt5)+randn(1);
[slopep4,xxp4,yyp4] = sl(obj4p, 0, ...
    gra4p, d);
[slopep5,xxp5,yyp5] = sl(obj5p, alpharo,...
    gra5p, d);

% just extra point
objp = fun(pt4+10*d)+randn(1);
grap = gra(pt4+10*d)+randn(1);
[slopep,xxp,yyp] = sl(objp, 10, ...
    grap, d);

objp2 = fun(pt4+2*d)+randn(1);
grap2 = gra(pt4+2*d)+randn(1);
[slopep2,xxp2,yyp2] = sl(objp2, 2, ...
    grap2, d);

% ======================================================== %
% interpolate with equation 7
% ======================================================== %
% Settings 
Y = [obj5p, objp2, objp, slopep5,slopep2, slopep]';
ts = [alpharo,2,10]';
% [Y,ts] = trainsample(fun,gra,3,d,12,pt4);
% Y = [objp2, objp, slopep2, slopep]';
% ts = [2,10]';

samplets = [0:0.001:12]';
samplegs = zeros(length(samplets),1);
sampleos = zeros(length(samplets),1);
for i = 1:length(samplets)
    tmobj = fun(pt4+samplets(i)*d);
    tmgra = gra(pt4+samplets(i)*d);
    [tsl,~,~] = sl(tmobj, samplets(i), tmgra, d);
    samplegs(i) = tsl;
    sampleos(i) = tmobj;
end
var = 1*ones(length(ts))';

% 1. posterior mean and variance
[mu,G] = posmean(Y,ts,samplets,var);
cov = poscov(G,ts,samplets);
f = [mu+2*sqrt(cov), flipdim(mu-2*sqrt(cov),1)];

% 2. expected improvement
EI = expectim(mu,cov);

% 3. Probabilistic Wolfe Conditions for Termination
[dmu,Gp] = dmean(Y,ts,samplets,var);
[Ks,Kds] = kvector(ts,samplets);
[mat,mbt,Caat,Cbbt,Cabt] = wolfcoef(Y,G,Gp,Ks,Kds,samplets);
[alim,blim,rhot] = coefcdf(mat,mbt,Caat,Cbbt,Cabt);

ind = blim == real(blim);
blim = blim(ind); alim = alim(ind); rhot = rhot(ind);
newsamplets = samplets(ind);
pt = zeros(1,length(blim));
for i = 1:length(blim)
    pt(i) = bvn(alim(i),inf,blim(i),inf,rhot(i));
end

% 4. marginal likelihood for optimization
L = marglik(ts,Y); % training data
% Ls= marglik(samplets,mu); % samples


% ======================================================== %
% Plot
% ======================================================== %
subplot(3,1,1)
fill([samplets; flipdim(samplets,1)]', f, [7 7 7]/8)
hold on;
plot(storetest(1,:),storetest(2,:));
plot(0,fun(pt4),'x')
plot(alpharo,fun(pt5),'o')
plot(samplets,mu)
plot(xxu4,yyu4,'LineWidth',2)
plot(xxu5,yyu5,'LineWidth',2)
plot(0,obj4p,'x');plot(xxp4,yyp4,'LineWidth',2);
plot(alpharo,obj5p,'x');plot(xxp5,yyp5,'LineWidth',2);
plot(10,objp,'o');plot(xxp,yyp,'LineWidth',2);
plot(2,objp2,'o');plot(xxp2,yyp2,'LineWidth',2);

xlabel('\alpha')
ylabel('f(\alpha)')
xlim([0,12])
ylim([-50,100])
title('GP with cubic spline')
grid on;

subplot(3,1,2)
plot(samplets,EI);
xlabel('\alpha')
ylabel('u_{EI}(\alpha)')
xlim([0,12])
ylim([-5,15])
title('Expected improvement')
grid on;

subplot(3,1,3)
plot(newsamplets,pt)
xlim([0,12])
xlabel('\alpha')
ylabel('p_{Wolfe}(\alpha)')
title('Wolfe Probability for Termination')
grid on;
% ======================================================== %
% function
% ======================================================== %
function [slope,xx,yy] = sl(obj, alpha, grad, dir)
    slope = grad'*dir;
    intercept = obj - slope*alpha;
    xx = [alpha-0.5,alpha+0.5];
    yy = [slope*xx(1)+intercept,slope*xx(2)+intercept];
end
function [alim,blim,rhot] = coefcdf(mat,mbt,Caat,Cbbt,Cabt)
alim=-mat./sqrt(Caat); blim=-mbt./sqrt(Cbbt);
rhot = Cabt./sqrt(Caat.*Cbbt);
end

function [mat,mbt,Caat,Cbbt,Cabt] = wolfcoef(Y,G,Gp,Ks,Kds,samplets)
ls = length(samplets); i0 = find(samplets==0);
mat = zeros(1,ls); mbt = zeros(1,ls);
Caat = zeros(1,ls); Cbbt = zeros(1,ls);
Cabt = zeros(1,ls);

c1 = 0.01; c2 = 0.5;

% Fundamental stuff 
mu0 = G(i0,:)*Y; mup0 = Gp(i0,:)*Y;
[tk00,tkd00,tdk00,tdkd00] = kernel(0,0);
% K00
k00 = tk00 - G(i0,:)*Ks(:,i0);
kd00 = tkd00 - G(i0,:)*Kds(:,i0);
dk00 = tdk00 - Gp(i0,:)*Ks(:,i0);
dkd00= tdkd00- Gp(i0,:)*Kds(:,i0);

for i = 1:ls
    it = find(samplets==samplets(i));
    % basic ingredients
    mut = G(it,:)*Y; mupt = G(it,:)*Y;
    % mat,mbt
    mat(i) = mu0-mut+c1*samplets(i)*mup0;
    mbt(i) = mupt-c2*mup0;
    
    % basic ingredients
    [tk0t,tkd0t,tdk0t,tdkd0t] = kernel(0,samplets(i));
    [tktt,tkdtt,tdktt,tdkdtt] = kernel(samplets(i),samplets(i));
    % K0t
    k0t = tk0t - G(i0,:)*Ks(:,it);
    kd0t = tkd0t - G(i0,:)*Kds(:,it);
    dk0t = tdk0t - Gp(i0,:)*Ks(:,it);
    dkd0t= tdkd0t- Gp(i0,:)*Kds(:,it);
    % Ktt
    ktt = tktt - G(it,:)*Ks(:,it);
    kdtt = tkdtt - G(it,:)*Kds(:,it);
    dktt = tdktt - Gp(it,:)*Ks(:,it);
    dkdtt= tdkdtt- Gp(it,:)*Kds(:,it);
    
    % C
    Caat(i)=k00+dkd00*(c1*samplets(i))^2+ktt+2*(c1*samplets(i)*(kd00-dk0t)-k0t);
    Cbbt(i)=dkd00*c2^2-2*c2*dkd0t+dkdtt;
    Cabt(i)=-c2*(kd00+c1*samplets(i)*dkd00)+dk0t*c2+kd0t+c1*samplets(i)*dkd0t-kdtt;
    
end

end

function [mu,G] = posmean(Y,ts,samplets,var)
mu = zeros(1,length(samplets));
lt = length(ts); G = zeros(length(samplets),2*lt);

K = zeros(2*lt,2*lt);
vI  = diag(var); [Xm,Ym] = meshgrid(ts);
ks = zeros(lt,lt);kds = zeros(lt,lt);
dks = zeros(lt,lt);dkds = zeros(lt,lt);

for i = 1:length(Xm)
    for j = 1:length(Xm)
        [ks(i,j),kds(i,j),dks(i,j),dkds(i,j)] = kernel(Xm(i,j),...
                                        Ym(i,j));
    end
end
K(1:lt,1:lt) = ks+vI;
K(1:lt,lt+1:2*lt) = kds;
K(lt+1:2*lt,1:lt) = dks; K(lt+1:2*lt,lt+1:2*lt) = dkds;
Kinv = inv(K);

for k = 1:length(samplets)
    tmkvec = zeros(1,2*lt);
    for l = 1:lt
        [tmk,~,tmdk,~]=kernel(samplets(k),ts(l));
        tmkvec(l) = tmk; tmkvec(l+lt) = tmdk;
    end
    mu(k) = tmkvec*Kinv*Y;
    G(k,:) = tmkvec*Kinv;
end

end

function [dmu,Gp] = dmean(Y,ts,samplets,var)
dmu = zeros(1,length(samplets));
lt = length(ts); Gp = zeros(length(samplets),2*lt);

K = zeros(2*lt,2*lt);
vI  = diag(var); [Xm,Ym] = meshgrid(ts);
ks = zeros(lt,lt);kds = zeros(lt,lt);
dks = zeros(lt,lt);dkds = zeros(lt,lt);

for i = 1:length(Xm)
    for j = 1:length(Xm)
        [ks(i,j),kds(i,j),dks(i,j),dkds(i,j)] = kernel(Xm(i,j),...
                                        Ym(i,j));
    end
end
K(1:lt,1:lt) = ks+vI; K(1:lt,lt+1:2*lt) = kds;
K(lt+1:2*lt,1:lt) = dks; K(lt+1:2*lt,lt+1:2*lt) = dkds;
Kinv = inv(K);

for k = 1:length(samplets)
    tmkvec = zeros(1,2*lt);
    for l = 1:lt
        [~,tmkd,~,tmdkd]=kernel(samplets(k),ts(l));
        tmkvec(l) = tmkd; tmkvec(l+lt) = tmdkd;
    end
    dmu(k) = tmkvec*Kinv*Y;
    Gp(k,:) = tmkvec*Kinv;
end
end

function [Ks,Kds] = kvector(ts,samplets)
lt = length(ts); ls = length(samplets);
Ks = zeros(2*lt,ls); Kds = zeros(2*lt,ls);

for k = 1:ls
    tmkvec1 = zeros(2*lt,1);
    tmkvec2 = zeros(2*lt,1);
    for l = 1:lt
        [tmk,tmkd,tmdk,tmdkd]=kernel(samplets(k),ts(l));
        tmkvec1(l) = tmk; tmkvec1(l+lt) = tmdk;
        tmkvec2(l) = tmkd; tmkvec2(l+lt) = tmdkd;
    end
    Ks(:,k) = tmkvec1;
    Kds(:,k) = tmkvec2;
end
end

function cov = poscov(G,ts,samplets)
lt = length(ts); cov = zeros(1,length(samplets));
for i = 1:length(G)
    tmkvec = zeros(2*lt,1);
    for l = 1:lt
        [stmk,~,~,~]=kernel(samplets(i),samplets(i));
        [tmk,~,tmdk,~]=kernel(samplets(i),ts(l));
        tmkvec(l) = tmk; tmkvec(l+lt) = tmdk;
    end
    cov(i) = stmk - G(i,:)*tmkvec;
end

end

function [k,kd,dk,dkd] = kernel(x,y) %cubic spline
% assume theta = 1;
theta = 1000;
tau = 10;
if x >= y % x = t, y = t', tau = 10
    k = (y+tau)^3/3+0.5*(x-y)*(y+tau)^2;
    kd = x*y-0.5*y^2;
    dk = 0.5*y^2;
    dkd = y;
elseif x < y
    k = (x+tau)^3/3+0.5*(y-x)*(x+tau)^2;
    kd = 0.5*x^2;
    dk = x*y-0.5*x^2;
    dkd = x;
end
k = theta*k;kd = theta*kd;dkd = theta*dkd;
dk = theta*dk;
end

function EI = expectim(mu,cov)
rho = min(mu);
tmp1 = 0.5*(rho-mu).*(1+erf((rho-mu)./real(sqrt(2*cov))));
tmp2 = sqrt(cov/2/pi).*exp(-(rho-mu).^2/2./cov);
EI = tmp1+tmp2;
end

function L = marglik(ts,Y)
lt = length(ts); m = 10; % m = tau
B = [m^3/3, m^2/2; m^2/2 m];
H = zeros(2,lt); H(1,:) = 1; H(2,:) = ts';

[Xm,Ym] = meshgrid(ts); 
K = zeros(lt,lt);
for i = 1:length(Xm)
    for j = 1:length(Ym)
        [tk,~,~,~] = kernel(Xm(i,j),Ym(i,j));
        K(i,j) = tk;
    end
end
y = Y(1:lt);

Kinv = inv(K); A = H*Kinv*H'; Ainv = inv(A);
DET = 0.5*(log(det(K))+log(det(A))+log(det(B)));
Z = Kinv-Kinv*H'*Ainv*H*Kinv;
L = -(0.5*y'*Z*y+DET+0.5*lt*log(2*pi));
end

function [storetest,d] = cursamples(xmax,x0,xf,fun)
d = (xf-x0)/norm(xf-x0);
storetest = zeros(2,length(0:0.001:xmax)); j = 0;
for alphar = 0:0.001:xmax
    j = j+1;
    tmpx = x0 + alphar*d;
    tmpobj = fun(tmpx);
    storetest(1,j) = alphar;
    storetest(2,j) = tmpobj;
end
% funaa = @(alphar)fun(x0 + alphar*d);
% alpharo = fminbnd(funaa,0,xmax);
end

function [Y,ts] = trainsample(fun,gra,num,dir,max,x0)
Y = zeros(2*num,1); ts = zeros(num,1);
randstep = max*rand(1,num); randnoise = randn(1,num);

for i = 1:num
    objp = fun(x0+randstep(i)*dir)+randnoise(i);
    grap = gra(x0+randstep(i)*dir);
    [slopep,~,~] = sl(objp, randstep(i), ...
    grap, dir);
    Y(i) = objp; Y(i+num) = slopep;
    ts(i) = randstep(i);
end
end

