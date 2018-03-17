%% Lab 8 - Symbol Timing Synchronization for Binary Baseband PAM
% Nicholas McKibben
% ECEn 485
% 2018-03-17

clear;
close all;

% This is to get PLL functions
if ~exist('lab5','dir')
    fprintf('Adding lab5 to path...\n');
    addpath('../lab5','-end');
end

% Params
M = 2;
N = 8;
E = 1;
A = sqrt(E);

% Matched filter design
beta = 1/2;
span = 12;
b = rcosdesign(beta,span,N);

% Loop filter design
BT = .01;
zeta = 1/sqrt(2);
order = 2;
k0 = 1; kp = 1;
[ ~,~,K1,K2 ] = LF(order,zeta,BT,1,k0,kp);

%% Detector time
load('bbtrdata.mat');
n = bbtrdata(1,:);
r_t = bbtrdata(2,:);

% Now do the filterin'
x_t = conv(b,r_t);

% Downsample by N/2
xk = x_t(1:N/2:end);

%% Loop through x(k)
% Initial conditions
ff_sa_prev = 0;
ff_sc_prev = 0;
ff_sf_prev = 0;
lf_sd_prev = 0;
strobe = 0;
eta = 0;
CNT_next = 0;
mu(1) = 0;

ted_sa_prev = 0;
ted_sa_prev_prev = 0;
k = 1;

for ii = 2:numel(xk)
    
    CNT = CNT_next;
    
    % Farrow Filter - cubic interpolation
    ff_sa = xk(ii);
    ff_sb = ff_sa*(1/6);
    ff_sc = ff_sa_prev; % previous sa value
    ff_sd = ff_sc*(-1/2);
    ff_se = ff_sb + ff_sd;
    ff_sf = ff_sc_prev; % previous sc value
    ff_sg = ff_sf*(1/2);
    ff_sh = ff_se + ff_sg;
    ff_si = ff_sf_prev; % previous sf value
    ff_sj = ff_si*(1/6);
    ff_sk = ff_sh + ff_sj;
    ff_sl = mu(k);
    ff_sm = ff_sk*ff_sl;
    
    ff_sn = ff_sa_prev;
    ff_so = ff_sn*(-1/2);
    ff_sp = ff_sc_prev;
    ff_sq = ff_so - ff_sp;
    ff_sr = ff_sf_prev;
    ff_ss = ff_sr*(1/2);
    ff_st = ff_sq + ff_ss;
    ff_su = ff_st + ff_sm;
    
    ff_sv = ff_sl*ff_su;
    
    ff_sw = ff_sa*(-1/6);
    ff_sx = ff_sa_prev;
    ff_sy = ff_sx + ff_sw;
    ff_sz = ff_sc_prev;
    ff_saa = ff_sz*(-1/2);
    ff_sac = ff_sy + ff_saa;
    ff_sab = ff_sf_prev;
    ff_sad = ff_sab*(-1/3);
    ff_sae = ff_sac + ff_sad;
    ff_sal = ff_sae + ff_sv;
    ff_saf = ff_sl*ff_sal;
    
    ff_sag = ff_sa_prev;
    ff_sah = ff_sc_prev;
    farrow_output = ff_saf + ff_sah;
    
    % update previous values
    ff_sa_prev = ff_sa;
    ff_sc_prev = ff_sc;
    ff_sf_prev = ff_sf;
    
    % Decision
    if strobe
        a(k) = A*sign(farrow_output);
    end
    
    % Zero-crossing TED
    % If strobe is low, TED outputs 0
    if ~strobe
        e(ii) = 0;
    else
        ted_sa = farrow_output;
        ted_sb = ted_sa_prev;
        ted_sc = ted_sa_prev_prev;
        ted_sd = sign(ted_sc);
        ted_se = sign(ted_sa);
        ted_sf = ted_sd - ted_se;
        ted_sg = ted_sb*ted_sf;
    
        e(ii) = ted_sg;
        
        % update previous values
        ted_sa_prev_prev = ted_sa_prev;
        ted_sa_prev = ted_sa;
    end
   
    % Loop filter it
    lf_sa = e(ii)*K1;
    lf_sb = e(ii)*K2;
    lf_sc = lf_sd_prev;
    lf_sd = lf_sb + lf_sc;
    lf_output = lf_sa + lf_sd;

    % update previous values
    lf_sd_prev = lf_sd;
    
    % mod 1 counter
    W = lf_output + 1/2;
    eta1 = eta - W;
    
    % update previous values
    eta = eta1;
    
    % update strobe
    strobe = mod(ii,M);
    
    if strobe
        k = k + 1;
        mu(k) = mod(eta/W,1);
    end
end

% Book's way
% CNT_next = 0;
% mu_next = 0;
% underflow = 0;
% vi = 0;
% TEDBuff = [ 0 0 ];
% k = 1;
% ted_e = [];
% mus = [];
% xk = xk.';
% for n = 2:length(xk)-2
%     CNT = CNT_next;
%     mu = mu_next;
%     mus(end+1) = mu;
%     
%     idx = n+2:-1:n-1;
%     v2 = 1/2*[ 1 -1 -1 1 ]*xk(idx);
%     v1 = 1/2*[ -1 3 -1 -1 ]*xk(idx);
%     v0 = xk(n);
%     xI = (mu*v2 + v1)*mu + v0;    
%     if underflow == 1
%         e = TEDBuff(1)*(sign(TEDBuff(2)) - sign(xI));
%         xx(k) = xI;
%         k = k + 1;
%         
%         ted_e(end+1) = e;
%     else
%         e = 0;
%     end
%     
%     vp = K1*e;
%     vi = vi + K2*e;
%     v = vp + vi;
%     W = 1/2 + v;
%     
%     CNT_next = CNT - W;
%     if CNT_next < 0
%         CNT_next = 1 + CNT_next;
%         underflow = 1;
%         mu_next = CNT/W;
%     else
%         underflow = 0;
%         mu_next = mu;
%     end
%     TEDBuff = [ xI; TEDBuff(1) ];
% end
% a = sign(xx);

%% Decode Message
SYNC = [ 1 0 1 1 0 1 0 0 1 0 1 1 0 1 0 0 ];
DataL = 1288;
m = a;
m(m < 0) = 0;
mstr = strjoin(string(m),'');
idx = strfind(mstr,strjoin(string(SYNC),''));

% try each possibility
for ii = 1:numel(idx)
    start = idx(end-ii+1) + numel(SYNC);
    o = 0;
    try
        m2 = m(start+o:start+DataL-1+o);
        message = m2ascii(m2,M);
        fprintf('Message %d: %s\n',ii,message);
    catch
    end
end

%% Let's make some pretty plots
figure(1);
plot(e);
title('TED output');
xlabel('samples');
ylabel('e(t)');

figure(2);
plot(mu);
title('Fractional Interpolation Interval');
xlabel('k');
ylabel('\mu');