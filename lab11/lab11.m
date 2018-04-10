%% Lab 11 - QPSK: Over-the-Air With Pluto Radios
% Nicholas McKibben
% ECEn 485
% 2018-04-08

clear;
close all;

% This is to get PLL functions
if ~exist('lab5','dir')
    fprintf('Adding lab5 to path...\n');
    addpath('../lab5','-end');
end
% This is to get m2ascii function
if ~exist('lab0','dir')
    fprintf('Adding lab0 to path...\n');
    addpath('../lab0','-end');
end
% This is to get pam helper functions
if ~exist('pam','dir')
    fprintf('Adding PAM to path...\n');
    addpath('../pam','-end');
end

%% Params
M = 4;
N = 4;
E = 1;
A = sqrt((3*E)/(2*(M - 1)));
w0 = 2*pi*.3;

% Matched filter design
beta = 1/2;
span = 12;
b = rcosdesign(beta,span,N);

% Loop filter design
BT = .01;
zeta = 1/sqrt(2);
order = 2;
k0 = 1; kp = 2.7*2;
[ ~,~,K1,K2 ] = LF(order,zeta,BT,1,k0,kp);

% DDS Loop filter design
dds_BT = .01;
dds_zeta = 1/sqrt(2);
dds_order = 2;
dds_k0 = 1; dds_kp = kp;
[ ~,~,dds_K1,dds_K2 ] = LF(dds_order,dds_zeta,dds_BT,1,dds_k0,dds_kp);

%% Detector
% Get the datar
load('TestInput.mat');
% r_t = qpsktrdata(2,:);
t = 0:numel(r);

% Mix it down
LOx = @(x) sqrt(2)*cos(w0*x);
LOy = @(x) -sqrt(2)*sin(w0*x);
% Ir_t = r_t.*LOx(t);
% Qr_t = r_t.*LOy(t);
% G = .001;
G = sqrt(2)*1e-3;
Ir_t = G*real(r);
Qr_t = G*imag(r);

% Do the filterin'
x_t = conv(b,Ir_t);
y_t = conv(b,Qr_t);

% Downsample by N/2
xk = x_t(1:N/2:end);
yk = y_t(1:N/2:end);

%% Loop through x(k)
% Initial conditions
ff_sa_prev0 = 0;
ff_sc_prev0 = 0;
ff_sf_prev0 = 0;
lf_sd_prev0 = 0;

ff_sa_prev1 = 0;
ff_sc_prev1 = 0;
ff_sf_prev1 = 0;
lf_sd_prev1 = 0;

strobe = 0;
eta = 0;
CNT_next = 0;
mu(1) = 0;

ted_sa_prev0 = 0;
ted_sa_prev_prev0 = 0;
ted_sa_prev1 = 0;
ted_sa_prev_prev1 = 0;

k = 1;

% DDS Initial conditions
dds_sr = 1;
dds_skp = 0;
dds_sop = 0;

% These must be initialized because of the strobe signal
dds_sh = 0;
w0dds = 0;

for ii = 2:numel(xk)
    
    % CCW Rotation
    in = xk(ii) + 1j*yk(ii);
    dds_sa = in*dds_sr;
    xk(ii) = real(dds_sa);
    yk(ii) = imag(dds_sa);
    
    CNT = CNT_next;
    
    % Farrow Filter for a0 - cubic interpolation
    ff_sa0 = xk(ii);
    ff_sb0 = ff_sa0*(1/6);
    ff_sc0 = ff_sa_prev0; % previous sa value
    ff_sd0 = ff_sc0*(-1/2);
    ff_se0 = ff_sb0 + ff_sd0;
    ff_sf0 = ff_sc_prev0; % previous sc value
    ff_sg0 = ff_sf0*(1/2);
    ff_sh0 = ff_se0 + ff_sg0;
    ff_si0 = ff_sf_prev0; % previous sf value
    ff_sj0 = ff_si0*(1/6);
    ff_sk0 = ff_sh0 + ff_sj0;
    ff_sl0 = mu(k);
    ff_sm0 = ff_sk0*ff_sl0;
    
    ff_sn0 = ff_sa_prev0;
    ff_so0 = ff_sn0*(-1/2);
    ff_sp0 = ff_sc_prev0;
    ff_sq0 = ff_so0 - ff_sp0;
    ff_sr0 = ff_sf_prev0;
    ff_ss0 = ff_sr0*(1/2);
    ff_st0 = ff_sq0 + ff_ss0;
    ff_su0 = ff_st0 + ff_sm0;
    
    ff_sv0 = ff_sl0*ff_su0;
    
    ff_sw0 = ff_sa0*(-1/6);
    ff_sx0 = ff_sa_prev0;
    ff_sy0 = ff_sx0 + ff_sw0;
    ff_sz0 = ff_sc_prev0;
    ff_saa0 = ff_sz0*(-1/2);
    ff_sac0 = ff_sy0 + ff_saa0;
    ff_sab0 = ff_sf_prev0;
    ff_sad0 = ff_sab0*(-1/3);
    ff_sae0 = ff_sac0 + ff_sad0;
    ff_sal0 = ff_sae0 + ff_sv0;
    ff_saf0 = ff_sl0*ff_sal0;
    
    ff_sag0 = ff_sa_prev0;
    ff_sah0 = ff_sc_prev0;
    farrow_output0 = ff_saf0 + ff_sah0;
    
    % update previous values
    ff_sa_prev0 = ff_sa0;
    ff_sc_prev0 = ff_sc0;
    ff_sf_prev0 = ff_sf0;
    
    % Farrow Filter for a1 - cubic interpolation
    ff_sa1 = yk(ii);
    ff_sb1 = ff_sa1*(1/6);
    ff_sc1 = ff_sa_prev1; % previous sa value
    ff_sd1 = ff_sc1*(-1/2);
    ff_se1 = ff_sb1 + ff_sd1;
    ff_sf1 = ff_sc_prev1; % previous sc value
    ff_sg1 = ff_sf1*(1/2);
    ff_sh1 = ff_se1 + ff_sg1;
    ff_si1 = ff_sf_prev1; % previous sf value
    ff_sj1 = ff_si1*(1/6);
    ff_sk1 = ff_sh1 + ff_sj1;
    ff_sl1 = mu(k);
    ff_sm1 = ff_sk1*ff_sl1;
    
    ff_sn1 = ff_sa_prev1;
    ff_so1 = ff_sn1*(-1/2);
    ff_sp1 = ff_sc_prev1;
    ff_sq1 = ff_so1 - ff_sp1;
    ff_sr1 = ff_sf_prev1;
    ff_ss1 = ff_sr1*(1/2);
    ff_st1 = ff_sq1 + ff_ss1;
    ff_su1 = ff_st1 + ff_sm1;
    
    ff_sv1 = ff_sl1*ff_su1;
    
    ff_sw1 = ff_sa1*(-1/6);
    ff_sx1 = ff_sa_prev1;
    ff_sy1 = ff_sx1 + ff_sw1;
    ff_sz1 = ff_sc_prev1;
    ff_saa1 = ff_sz1*(-1/2);
    ff_sac1 = ff_sy1 + ff_saa1;
    ff_sab1 = ff_sf_prev1;
    ff_sad1 = ff_sab1*(-1/3);
    ff_sae1 = ff_sac1 + ff_sad1;
    ff_sal1 = ff_sae1 + ff_sv1;
    ff_saf1 = ff_sl1*ff_sal1;
    
    ff_sag1 = ff_sa_prev1;
    ff_sah1 = ff_sc_prev1;
    farrow_output1 = ff_saf1 + ff_sah1;
    
    % update previous values
    ff_sa_prev1 = ff_sa1;
    ff_sc_prev1 = ff_sc1;
    ff_sf_prev1 = ff_sf1;
    
    % Decision
    if strobe
        a0(k) = A*sign(farrow_output0);
        a1(k) = A*sign(farrow_output1);
        
        xp(ii) = farrow_output0;
        yp(ii) = farrow_output1;
        
        % PED
        dds_sf = a0(k)*imag(dds_sa);
        dds_sg = a1(k)*real(dds_sa);
        dds_sh = dds_sf - dds_sg;
        
        % save the error
        dds_e(k) = dds_sh;
    else
        dds_sh = 0;
    end
    
    % DDS Loop Filter
    dds_si = dds_sh*dds_K1;
    dds_sj = dds_sh*dds_K2;
    dds_sl = dds_skp; % previous sk
    dds_sk = dds_sj + dds_sl;
    dds_sm = dds_si + dds_sk;
    
    % DDS
    dds_sn = dds_sm*dds_k0;
    dds_sp = dds_sop; % previous so
    dds_so = dds_sn + w0dds + dds_sp;
    dds_sq = cos(dds_sp) + 1j*sin(dds_sp);
    dds_sr = conj(dds_sq);
    
    % Update previous DDS vals
    dds_skp = dds_sk;
    dds_sop = dds_so;
    
    
    ted_sa0 = farrow_output0;
    ted_sb0 = ted_sa_prev0;
    ted_sc0 = ted_sa_prev_prev0;
    ted_sd0 = sign(ted_sc0);
    ted_se0 = sign(ted_sa0);
    ted_sf0 = ted_sd0 - ted_se0;
    ted_sg0 = ted_sb0*ted_sf0;
    
    % update previous values
    ted_sa_prev_prev0 = ted_sa_prev0;
    ted_sa_prev0 = ted_sa0;
    
    ted_sa1 = farrow_output1;
    ted_sb1 = ted_sa_prev1;
    ted_sc1 = ted_sa_prev_prev1;
    ted_sd1 = sign(ted_sc1);
    ted_se1 = sign(ted_sa1);
    ted_sf1 = ted_sd1 - ted_se1;
    ted_sg1 = ted_sb1*ted_sf1;
    
    % update previous values
    ted_sa_prev_prev1 = ted_sa_prev1;
    ted_sa_prev1 = ted_sa1;
    
    if strobe
        e(ii) = ted_sg0 + ted_sg1;
    else
        e(ii) = 0;
    end
   
    % Loop filter it
    lf_sa = e(ii)*K1;
    lf_sb = e(ii)*K2;
    lf_sc = lf_sd_prev0;
    lf_sd = lf_sb + lf_sc;
    lf_output = lf_sa + lf_sd;

    % update previous values
    lf_sd_prev0 = lf_sd;
    
    % mod 1 counter
    W = lf_output + 1/2;
    eta1 = eta - W;
    
    % update previous values
    eta = eta1;
    
    % update strobe
    CNT_next = CNT - W;
    if CNT_next < 0
        CNT_next = CNT_next + 1;
        strobe = 1;
        k = k + 1;
        mu(k) = CNT/W;
    else
        strobe = 0;
    end
end

%% Differential Encoding
a0 = a0; % Choose sign to get the right answer
a1 = a1;  % Choose sign to get the right answer

d_hat0 = zeros(size(a0));
d_hat1 = zeros(size(a1));
d_hat0(a0 > 0) = 1;
d_hat1(a1 > 0) = 1;

d_hat = [ d_hat0; d_hat1 ];
d_hat = char(strjoin(string(d_hat(:).'),''));

% assume 00 at beginning
d_hat = strcat('00',d_hat);

idx = 0;
for ii = (1+log2(M)):log2(M):numel(d_hat)
    idx = idx + 1;
    cur = strcat(d_hat(ii),d_hat(ii+1));
    prev =  strcat(d_hat(ii-2),d_hat(ii-1));
    
    if strcmp(prev,'00') && strcmp(cur,'00')
        b(idx) = 0;
    elseif strcmp(prev,'00') && strcmp(cur,'01')
        b(idx) = 2;
    elseif strcmp(prev,'00') && strcmp(cur,'10')
        b(idx) = 1;
    elseif strcmp(prev,'00') && strcmp(cur,'11')
        b(idx) = 3;
        
    elseif strcmp(prev,'01') && strcmp(cur,'00')
        b(idx) = 1;
    elseif strcmp(prev,'01') && strcmp(cur,'01')
        b(idx) = 0;
    elseif strcmp(prev,'01') && strcmp(cur,'10')
        b(idx) = 3;
    elseif strcmp(prev,'01') && strcmp(cur,'11')
        b(idx) = 2;
        
        
    elseif strcmp(prev,'10') && strcmp(cur,'00')
        b(idx) = 2;
    elseif strcmp(prev,'10') && strcmp(cur,'01')
        b(idx) = 3;
    elseif strcmp(prev,'10') && strcmp(cur,'10')
        b(idx) = 0;
    elseif strcmp(prev,'10') && strcmp(cur,'11')
        b(idx) = 1;
        
    elseif strcmp(prev,'11') && strcmp(cur,'00')
        b(idx) = 3;
    elseif strcmp(prev,'11') && strcmp(cur,'01')
        b(idx) = 1;
    elseif strcmp(prev,'11') && strcmp(cur,'10')
        b(idx) = 2;
    elseif strcmp(prev,'11') && strcmp(cur,'11')
        b(idx) = 0;
    end
end

%% Find Unique Word
uw = ones(1,64);
idx = strfind(strjoin(string(b),''),strjoin(string(uw),''));

DataL = 3206/2;

% try each possibility
for ii = 1:numel(idx)
    start = idx(end-ii+1) + numel(uw);
    o = 0;
    try
        s_hat = b(start+o:start+DataL-1+o);
        fprintf('Message %d: %s\n',ii,m2ascii(s_hat,M));
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

%% Eye Diagram
eyediagram(x_t(end-DataL:end),N);
title('x(t)');
eyediagram(y_t(end-DataL:end),N);
title('y(t)');

%% Constellation
scatterplot(xp(end-DataL:end) + 1j*yp(end-DataL:end));
title('Constellation');

%% DDS Error
figure;
plot(dds_e);
title('DDS Output');