A = readmatrix('priceandrainfalldata.xlsx');
t = A(:,1);
Rt = A(:,2);
yt = A(:,3);

n = 90;
Q = diag([1 0 0]); %noise of xt
R = 1; %noise of yt
Ft = @(Rt)[0.6 0 Rt; 0 1 0; 0 0 1]; %state transition map
At = [0.5 1 0]; %state-to-observation map

xa = zeros(3,1,n+1); %analysis means
xf = zeros(3,1,n+1); %forecast means
Pa = zeros(3,3,n+1); %analysis cov mat
Pf = zeros(3,3,n+1); %forecast cov mat

xa(:,1,1) = [12;0;0]; %X0 = 12, prior mu, prior beta ~N(0,100)
Pa(:,:,1) = diag([0;100;100]);

%Kalman filter loop
for k = 2:(n+1) %starting at t=1
    %forecast
    F=Ft(Rt(k-1));
    xf(:,1,k) = F*xa(:,1,k-1);
    Pf(:,:,k) = F*Pa(:,:,k-1) * F' + Q;

    %analysis
    S = At*Pf(:,:,k)*At' + R;
    xa(:,1,k) = xf(:,1,k) + Pf(:,:,k)*(At')/S*(yt(k-1)-At*xf(:,1,k)); %update mean
    Pa(:,:,k) = Pf(:,:,k) - Pf(:,:,k)*(At')/S*At*Pf(:,:,k); %update cov

end

%question 4a 
%trajectory plot for all components
tgrid = 1:n;
xavec = squeeze(xa(:,1,2:end));
figure; 
plot(tgrid, xavec', 'LineWidth',1); 
title('Trajectory of x_t, \mu, and \beta');
legend('x_t','\mu','\beta');
grid on;

%individual trajectory plots with 95% confidence band
for i=1:3
    variable = {'x','\mu','\beta'};
    subscript = {'t',' ',' '};

    xa_vec = squeeze(xa(i,1,2:end));
    Pa_vec = squeeze(Pa(i,i,2:end));
    
    xa_vec = xa_vec(:).';
    Pa_vec = Pa_vec(:).';

    index = 1:n;
    xconf = [index index(end:-1:1)];
    yconf = [xa_vec+norminv(0.975)*sqrt(Pa_vec) ...
        xa_vec(end:-1:1)-norminv(0.975)*sqrt(Pa_vec(end:-1:1))];
    figure;
    p = fill(xconf,yconf,'green', 'DisplayName', '95% Confidence Band');
    p.FaceColor = [0.8 1 0.8];      
    p.EdgeColor = 'none'; 
    hold on;

    h = plot((1:n), xa_vec,'-*k', 'DisplayName', 'Analysis mean');

    % creating plot titles with proper letter formatting
    switch i
        case 1
            title('Trajectory of $x_t$', 'Interpreter','latex')
            ylabel('Values of $x_t$', 'Interpreter','latex')
        case 2
            title('Trajectory of $\mu$', 'Interpreter','latex')
            ylabel('Values of $\mu$', 'Interpreter','latex')
        case 3
            title('Trajectory of $\beta$', 'Interpreter','latex')
            ylabel('Values of $\beta$', 'Interpreter','latex')
    end

    xlabel('Time Steps, t')
    legend([h p], 'Location','northwest');    

    grid on;
    hold off;
end

%question 4b
%posterior distribution of mu and beta
mean_mu_beta = xa(2:3,1,n+1); mean_mu_beta % mean(mu,beta) at t=n 
cov_mu_beta = Pa(2:3,2:3,n+1); cov_mu_beta %cov(mu, beta) at t=n

sd_mu_beta = sqrt(diag(cov_mu_beta)); sd_mu_beta   %sd of mu and beta

z = norminv(0.975);
lower_bound = mean_mu_beta - z*(sd_mu_beta);
upper_bound = mean_mu_beta + z*(sd_mu_beta);
CI = [lower_bound upper_bound];
%displaying the CI in a table 
disp(array2table(CI, 'VariableNames', {'Lower','Upper'}, 'RowNames', {'mu','beta'}));

%question 4c
phi = (0:1:360)*(2*pi)/360;
Z = sqrt(chi2inv(0.95,2))*[cos(phi); sin(phi)]; % 95% boundary circle
ellipseconf = mean_mu_beta + sqrtm(cov_mu_beta)*Z;

%plotting 95% confidence region
figure;
p = fill(ellipseconf(1,:),ellipseconf(2,:),'green', 'DisplayName', '95% Confidence Region');
p.FaceColor = [0.8 1 0.8];
p.EdgeColor = 'none';
hold on;
h = plot(mean_mu_beta(1), mean_mu_beta(2),'*k', 'DisplayName', 'Posterior Mean Vector');
legend([h p], 'Location','southeast');
xlabel('\mu'); ylabel('\beta');
title('95% credible region for posterior parameters (\mu, \beta)');
grid on;
hold off;

%question 4d
a  = [1; -10]; %where R91 = 10 
mean_theta = (a')*mean_mu_beta; mean_theta %mean of theta 
var_theta = (a')*cov_mu_beta * a; var_theta %variance of theta
sd_theta = sqrt(var_theta); sd_theta %sd of theta

CI_lower = mean_theta - z*sd_theta;
CI_upper = mean_theta + z*sd_theta;
CI_theta = [CI_lower CI_upper]; CI_theta %95% credible interval 



