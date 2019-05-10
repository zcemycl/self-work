clear all; clc;
% BFGS  
% ================================================= %
% Settings
% ================================================= %
% cost function
fun = @(x)3*x(1)^2+x(2)^2+55*x(3)^2+2*x(4)^2+x(5)^2;
gra = @(x)2*[3*x(1),x(2),55*x(3),2*x(4),x(5)]';

% initialization points
tol = 1e-4;
H0 = diag(ones(1,5));
% x0 = randi([100 100000],1,5)';

x0 = [10,0,24,31,11]';
% x0n= x0+10*randn(5,1);
dx = [5.3767,18.3389,-22.5885,8.6217,3.1877]';
x0n = x0+dx;

% memory
M = x0n;

% % try carl minimizer 
% [Xt fXt Ct] = minimize(x0n, 'testobj', 25);
% ================================================= %
% Line search
% ================================================= %
i = 0; xnew = x0n; Hnew = H0;
while norm(gra(xnew)) > tol
    i = i+1;
    Hold = Hnew; xold = xnew;
    [Hnew,xnew,bnew] = hesspt(Hold,xold,fun,gra);
    if not(bnew)
        break
    else
        M = [M,xnew];
    end    
end
% ================================================= %
% Plot
% ================================================= %
MF = zeros(1,5);
nmx = zeros(1,5);
for i = 1:length(MF)
    MF(i) = fun(M(:,i));
    nmx(i) = norm(M(:,i));
end
x1st = M(1,:);
plot(x1st,MF,'-o');
hold on
for ii = 1:length(M(1,:))
    text(x1st(ii),MF(ii),num2str(ii),'Color','r')
end
% just a polynomial fit along x1
p = polyfit(x1st,MF,2);
xm = linspace(-max(x1st),max(x1st));
ym = polyval(p,xm);
plot(xm,ym);
% cubic spline fit
xx = -max(x1st)*0.5:0.1:max(x1st);
yy = spline(x1st,MF,xx);
plot(xx,yy)

% plot(Xt(1),fXt)

% legend('iterative','polynomial','spline',...
%     'Location','north')
xlabel('x_1')
ylabel('f(x)')


% along only one direction
figure(2)
d = (M(:,5)-M(:,4))/norm(M(:,5)-M(:,4));
storetest = zeros(2,1001); j = 0;
for alphar = -20:0.01:30
    j = j+1;
    tmpx = M(:,4) + alphar*d;
    tmpobj = fun(tmpx);
    storetest(1,j) = alphar;
    storetest(2,j) = tmpobj;
end
plot(storetest(1,:),storetest(2,:));
funaa = @(alphar)fun(M(:,4) + alphar*d);
alpharo = fminbnd(funaa,0,10);
% ================================================= %
% GP
% ================================================= %
% meanfunc = [];
% covfunc = @covSEiso;
% likfunc = @likGauss;
% % test sample
% xs = linspace(-15,15,10000)'; % predictive mean
% 
% hyp = struct('mean',[], 'cov',[-0.5,0],'lik',0);
% hyp2 = minimize(hyp, @gp, -1000, @infGaussLik, ...
%     meanfunc,covfunc, likfunc, x1st,MF);
% [mu s2] = gp(hyp2, @infGaussLik,meanfunc,...
%     covfunc,likfunc, x1st,MF,xs);
% 
% f = [mu+2*sqrt(s2); flipdim(mu-2*sqrt(s2),1)];
% fill([xs; flipdim(xs,1)], f, [7 7 7]/8);
% plot(xs, mu); 
% ================================================= %
% Functions
% ================================================= %
function bool = wolfe(fun,gra,x0,x1,d,a)
c1 = 1e-4; c2 = 0.9;
bool1 = fun(x1)<=fun(x0)+c1*a*d'*gra(x0);
bool2 = norm(d'*gra(x1))<=c2*norm(d'*gra(x0));
if and(bool1,bool2)
    bool = true(1);
end
end

function [H1,x1,bool1] = hesspt(H0,x0,fun,gra)
d0 = -inv(H0)*gra(x0); 
funw = @(a,d,x)fun(x+a*d);
funa = @(a)funw(a,d0,x0);
alpha0 = fminbnd(funa,-max(x0),max(x0));
disp(alpha0)
try
    s0 = alpha0*d0; x1 = x0+s0;
catch
    bool1 = false(1); H1 = []; x1 = [];
    return;
end
q1 = gra(x1)-gra(x0);
H1 = H0+q1*q1'/(q1'*s0)-s0*s0'/(s0'*s0);
try
    bool1 = wolfe(fun,gra,x0,x1,d0,alpha0);
catch
    bool1 = false(1);
end

end
