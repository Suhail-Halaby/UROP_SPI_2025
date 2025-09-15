clear
clc
close all

y = [6;4;2;3;5;8;10;12;8;6;2;4;3;2;4;2;3;5;8;10;12;8;6;2;4;3;2];
x = linspace(0,length(y)-1,length(y));

xp = linspace(0,length(y),4*length(y));

hold on
plot(x,y,':.k',MarkerSize=15,MarkerEdgeColor='r')
yline(0)



%% Stats
mu = mean(y);
sigma = var(y);
noise = 1e-6;
lscale = 1.2;
%% Evaluating Kernel


RBF = @(x1, x2) sigma^2 * exp( - ( (x1 - x2).^2 ) / (2*lscale^2) );

% sample variance term

Kxx = zeros(length(x));

for i = 1:length(x)
    for j = 1:length(x)
        Kxx(i,j) = RBF(x(i),x(j));
    end
end

% covariance term

Kxxp = zeros(length(x),length(xp));

for i = 1:length(x)
    for j = 1:length(xp)
        Kxxp(i,j) = RBF(x(i),xp(j));
    end
end

% posterior variance term

Kxpxp = zeros(length(xp));

for i = 1:length(xp)
    for j = 1:length(xp)
        Kxpxp(i,j) = RBF(xp(i),xp(j));
    end
end


A = Kxx + noise*eye(size(Kxx));
B = Kxpxp;
C = Kxxp;

clear Kxx Kxxp Kxpxp

% Schol Complements

S1 = B - C'*inv(A)*C;
S2 = B - C'*(A\C);
%S3 = B - C'*chol(A)*C;

%% Posterior Predictive Mean

mu_star = mu + C'*inv(A)*(y-mu);
plot(xp,mu_star,'k',LineWidth=1.5);

%% Confidence Band

std_star = sqrt(diag(S2));

upper = mu_star + 2*std_star;
lower = mu_star - 2*std_star;


plot(xp,upper)
plot(xp,lower)


fill([xp; flipud(xp)], [upper'; flipud(lower')], [0.8 0.8 1], 'FaceAlpha', 0.3, 'EdgeColor', 'b');
