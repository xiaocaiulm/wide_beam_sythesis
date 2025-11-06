%% MM Epigraph Optimization for Antenna Array Weights
% Designs transmit weights for a ULA via an MM (Majorization–Minimization)
% epigraph formulation. Supports both continuous modulus and discrete-phase
% constraints, and can sweep the penalty parameter to study performance.
% Requirements: CVX toolbox must be installed and on the MATLAB path.
% Outputs: figures and a .mat file with selected weights.
N = 64;   %transmitting antennas
power =1;
width =60;
%% Radar parameters
delta=180;
theta=-pi/2:(pi/delta):pi/2;
target_DoA=0; 
w_index=0;
c=3e8;
fc=3.2e9;
lamda=c/fc; 
spacing=0.5*lamda;  % antenna spacing
spac = spacing/lamda;
beam_width=width*(delta/180);
loc = [spacing*[0:N-1]'];
A = kron(sin(theta'), loc(:)') ;
A = exp(2*pi*1i/lamda*A);
a = A.';
% for tt=1:N
%     for jj=1:length(theta) 
%         a(tt,jj)=exp(j*2*pi*spac*(tt-ceil((N)/2))*sin(theta(jj)));  %index (antennas / angle)   a(θ)=[1 ej2πΔsin(θ) ... ej2π(N −1)Δsin(θ)]  Δ=spacing/lamda=1/2
%     end
% end

theta_low = delta/2 +1 + target_DoA*(delta/180) -width/2 *(delta/180);
theta_up = delta/2 +1 + target_DoA*(delta/180) +width/2*(delta/180);

% y_cvx =a' * Rd;

% minimize y s.t. theta \in (theta_low,theta_up)
% y = y_cvx([theta_low:theta_up],:);   %21*400 size

a_target_complex= a(:,[theta_low:theta_up]);

a_target = getH(a_target_complex.',beam_width+1,N)';


% |a_target'* w *w'*a_target|



% sigma = 0.05;

zi = @(i,w) sum((getH(a_target_complex(:,i).',1,N) * w).^2); 
p = @(x) x'*x; %power
f = @(i,w) zi(i,w);

zi_cell = @(w) arrayfun(@(i) zi(i,w), 1:(beam_width+1),'UniformOutput', false);
min_z = @(w) min(cell2mat(zi_cell(w)));

% f = @(w) sigma * log(sum(exp(p(a_target'*w)/sigma)));

% ------2. A NEGATIVE SQUARE PENALTY METHOD--------------
% lambda = 0.1;
F =  @(i,w) f(i,w) ; %% '-'

% discrete set --> continuous
L = 4;
discrete_theta = [0:L-1]*2*pi/L + pi/L;
s = exp(1i*discrete_theta);


% -------------3. MM method ---------------------------
% object function: min F(w) s.t. w ~ K =convhull(s)
%initialize w0
% w0_complex = w_gsc;
w0_complex = s(randi(length(s), N,1)).'; % intin
% ialize w0
% for tt=1:N
%     w0_complex(tt)=exp(1i*2*pi*spac*(tt-ceil((N)/2))*sin((target_DoA)*pi/180));  %index (antennas / angle)   a(θ)=[1 ej2πΔsin(θ) ... ej2π(N −1)Δsin(θ)]  Δ=spacing/lamda=1/2
% end
w0 = getw(w0_complex);

% construct a upper bound for F(w)


% Compute gradient component for each column and sum them up
gradient_zi = @(i,w) getH(a_target_complex(:,i).',1,N)'*getH(a_target_complex(:,i).',1,N) * w;   % ai^T ai w
grad_f = @(i,w) (gradient_zi(i,w))';
gradientF = @(i,w) (grad_f(i,w));%

G = @(i, x, w) F(i,w) + gradientF(i,w)* (x - w);%

ai = @(i) getH(a_target_complex(:,i).',1,N);
G2 = @(i,x,w) F(i,w) + 2*(ai(i)*x-ai(i)*w)'*(ai(i)*w);

% penalty = @(x,w) lambda* w'*w + lambda * 2* (x-w)'*w;


% y0 = a' * w0_complex;
% plot(theta*180/pi, 20*log10(abs(y0)), 'LineWidth', 1.5,'color','r');
% hold on
% y = a' * x;
% figure(2)
% plot([-90:90], 20*log10(abs(y)),'LineWidth', 1.5,'color','b');
% hold on

%% Continuous phase optimization
% Iterative optimization
tolerance =1e-4; % A small threshold
maxIter =200; % Maximum number of iterations
w_prev = w0;
T_prev =0;
for i = 1:maxIter
    cvx_begin quiet
        variable x(N,1) complex
        variable combi(L, N) nonnegative
        variable T
        maximize T        
        subject to
        for be = 1: (beam_width+1)
%          T <= G(be,getw(x),w_prev);
         T <= G2(be,getw(x),w_prev);
        end
        for ii = 1: N
            abs(x(ii)) <= 1;
        end
    cvx_end
    
    % Update
    w_next = getw(x);
    w_next_compolex = x;
    T_next = T;
    % Check for convergence
    if mod(i,5) ==0 
        disp(i);
        disp(sum(abs(w_next-w_prev)));
%         y = a' * w_next_compolex;
%         plot([-90:90], 20*log10(abs(y)),'LineWidth', 1.5,'color','k');
%         hold on
    end
%     if abs(T_prev-T_next) < tolerance
     if sum(abs(w_next-w_prev)) < tolerance
        break;
    end
    w_prev = w_next;
    T_prev = T_next;

end
%%
figure ;
w_cont = w_next_compolex; 
y_con = a.' * w_cont;
plot(theta*180/pi, 20*log10(abs(y_con)),'LineWidth', 1.5,'color','k');
hold on


%% Discrete phase optimization
lambda = 0.1;
penalty = @(x,w) lambda* w'*w + lambda * 2* (x-w)'*w;
% penalty = @(x,w) lambda* w'*w + lambda * 2* (x-w)'*w;
% Iterative optimization
tolerance = 1e-4; % A small threshold
maxIter = 500; % Maximum number of iterations
w_prev = w0;
for i = 1:maxIter
    cvx_begin quiet
        variable x(N,1) complex
        variable combi(L, N) nonnegative
        variable T
        maximize T + penalty(getw(x),w_prev)
        subject to
        for be = 1: (beam_width+1)
         T <= G2(be,getw(x),w_prev);
        end
        for ii = 1: N
            x(ii) == s * combi(:,ii);
            sum(combi(:,ii)) == 1;
        end
    cvx_end
    
    % Update
    w_next = getw(x);
    w_next_compolex = x;
    T_next = T;
    % Check for convergence
    if mod(i,5) ==0 
        disp(i);
        disp(sum(abs(w_next-w_prev)));
    end
    if (abs(T_prev-T_next) < tolerance)
        break;
    end
 
    w_prev = w_next;
    T_prev = T_next;
end

w_results = projection(w_next_compolex, s);

figure;
y_con = a' * w_cont;
plot(theta*180/pi, 20*log10(abs(y_con)),'LineWidth', 1.5,'color','k');
hold on

y1 = a' * w_next_compolex;
y2 = a' * w_results;
plot([-90:90], 20*log10(abs(y1)),'LineWidth', 1.5,'color','g');
hold on
plot([-90:90], 20*log10(abs(y2)),'LineWidth', 1.5,'color','b');
title(['lambda = ',num2str(lambda)],'FontSize',16)
legend('continuous','before projection','after projection')

%% Discrete phase with lambda sweep

ll = 0.1:0.02:3;

min_power_before_proj = zeros(size(ll));
min_power_after_proj = zeros(size(ll));
avg_power_before_proj = zeros(size(ll));
avg_power_after_proj = zeros(size(ll));
w_before = zeros(N, 21);
w_after = zeros(N,21);
w_after_c = zeros(N,21);
kk =1;
T_prev = 0;

for lambda = ll
penalty = @(x,w) lambda* w'*w + lambda * 2* (x-w)'*w;
% penalty = @(x,w) lambda* norm(x,1);
% Iterative optimization
tolerance = 1e-5; % A small threshold
maxIter = 150; % Maximum number of iterations
%w_prev = w0;
w_prev = getw(w_results);
for i = 1:maxIter
    cvx_begin quiet
        variable x(N,1) complex
        variable combi(L, N) nonnegative
        variable T
        maximize T + penalty(getw(x),w_prev)
        subject to
        for be = 1: (beam_width+1)
         T <= G2(be,getw(x),w_prev);
        end
        for ii = 1: N
            x(ii) == s * combi(:,ii);
            sum(combi(:,ii)) == 1;
        end
    cvx_end
    
    % Update
    w_next = getw(x);
    w_next_compolex = x;
    T_next = T;
    % Check for convergence
%     if mod(i,5) ==0 
%         disp(i);
%         disp(sum(abs(w_next-w_prev)));
% %         y = a' * w_next_compolex;
% %         plot([-90:90], 20*log10(abs(y)),'LineWidth', 1.5,'color','k');
% %         hold on
%     end
    if abs(T_prev +lambda*norm(w_prev)^2 -T_next -lambda* norm(w_next)^2 ) < tolerance
        break;
    end
 
    w_prev = w_next;
    T_prev = T_next;
end

w_results = projection(w_next_compolex, s);
w_results_c = projection_con(w_next_compolex);

w_before(:,kk) = w_next_compolex;
w_after_c(:,kk) = w_results_c;
w_after(:,kk) = w_results;

y1 = a' * w_next_compolex; % before projection
y2 = a' * w_results;  %after project_L_quan
y2c = a'* w_results_c;


y1_target = y1([theta_low:theta_up],:);
y2_target = y2([theta_low:theta_up],:);
y2_target_c = y2c([theta_low:theta_up],:);

p1 = 20*log10(abs(y1_target));
p2 = 20*log10(abs(y2_target));

min_power_before_proj(kk) = min(p1);
min_power_after_proj(kk) = min(p2);
avg_power_before_proj(kk) = mean(p1);
avg_power_after_proj(kk) = mean(p2);
% disp(kk);

kk = kk +1;
end


%%
for kk = 1:length(ll)
y1 = a' * w_before(:,kk); % before projection
y2 = a' * w_after(:,kk);  %after project_L_quan
y2c = a'* w_after_c(:,kk);

y1_target = y1([theta_low:theta_up],:);
y2_target = y2([theta_low:theta_up],:);
y2_target_c = y2c([theta_low:theta_up],:);

p1 = 20*log10(abs(y1_target));
p2 = 20*log10(abs(y2_target));

min_power_before_proj(kk) = min(p1);
min_power_after_proj(kk) = min(p2);
avg_power_before_proj(kk) = mean(p1);
avg_power_after_proj(kk) = mean(p2);
end

dec =min_power_after_proj; % can be adjusted

[~,m] = max(dec);

y1 = a' * w_before(:,m); % before projection
y2 = a' * w_after(:,m);  %after project
y3 = a' * w_after_c(:,m);  %after project continuous phase
y_con = a' * w_cont;
figure;
hold on
plot(theta*180/pi, 20*log10(abs(y_con)),'LineWidth', 1.5,'color','k');
hold on
plot(theta*180/pi, 20*log10(abs(y1)),'LineWidth', 1.5,'color','g');
hold on
plot(theta*180/pi, 20*log10(abs(y3)),'LineWidth', 1.5,'color','b');
hold on
plot(theta*180/pi, 20*log10(abs(y2)),'LineWidth', 1.5,'color','r');
title(['L= ', num2str(L), ', \lambda = ', num2str(ll(m))], 'FontSize', 16);
legend('per antenna power constraint','convex hull','constant modulus','discrete phase')
xlabel('\theta [degree]');
ylabel('Beam pattern [dB]')
%%
figure;
plot(ll,min_power_before_proj,ll,min_power_after_proj,ll,avg_power_before_proj,ll,avg_power_after_proj);
legend('min power before projection','min power after projection','average power before projection','average power after projection');
xlabel('\lambda');
ylabel('Minimum power at target region[dB]')

%%
figure;
plot(ll,min_power_before_proj,ll,min_power_after_proj);
legend('constant modulus','discrete phase');
%%
filename = sprintf('lambda=%d.mat',m);
save(filename, 'w_cont', 'w_before', 'w_after','w_after_c');

%% plot power

N = [64,128,256,384,512];
tar = [19.740627202670563 23.603767842579966 29.333348882510062 32.690671327891515 35.115865533186046];
side = [10.870517576770226,13.850256086603839,16.179742138991934,17.875083047617764, 18.954895527434665];
% ratio = [7.709229286838878   9.448245636739712  20.670961321618165 30.308108131270018 41.313976731161979];
ratio = tar - side;
figure
yyaxis left
plot(N,tar,'LineWidth', 1.5)
hold on
plot(N,side,'LineWidth', 1.5)
xlabel('Number of elements $N$')
ylabel('Average power [dB]')
yyaxis right
plot(N,ratio,'LineWidth', 1.5)
legend('Average power at target area','Average power at side area','Target-to-side power ratio')
ylabel('Ratio')

