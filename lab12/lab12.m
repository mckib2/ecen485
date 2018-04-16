%% Lab 12
% Nicholas McKibben
% ECEn 485
% 2018-04-16

close all;
clear;

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

% (1,2) Separate I,Q components with appropriate gain
% G = sqrt(2)*1e-3;
G = 4e-3;
% G = 5e-3;
% load('TestInput.mat');

% Let's recieve the signal
RX = sdrrx('Pluto');
RX.RadioID = 'usb:0';
RX.CenterFrequency = 425e6;
RX.GainSource = 'AGC Slow Attack';
RX.ChannelMapping = 1;
RX.BasebandSampleRate = 4*40e3;
RX.OutputDataType = 'int16';
RX.SamplesPerFrame = 36600;

% r = [];
% for ii = 1:20
%     r = [ r; RX() ];
% end
load('data.mat');

Ir = G*real(r);
Qr = G*imag(r);

% Calibrating input to get avg radius of sqrt(2)...
p = rcosdesign(1/2,12,4);
% x = G*conv(r,p);
% plot(x(6*4+1:6*4+1000));
% grid on;
% axis square;

% (3) Choose samples/symbol
M = 4;
N = 4;
E = 2;
% E = sqrt(8);
A = sqrt((3*E)/(2*(M - 1)));

% (4) Match filter and downsample by 2
x_t = conv(p,Ir);
y_t = conv(p,Qr);
xk = x_t(1:N/2:end);
yk = y_t(1:N/2:end);

% (5) Design Loop filter for Carrier Phase sync
cp_BT = 0.01; % may need to increase
cp_zeta = 1/sqrt(2);
k0 = 1; ted_kp = 2.7*2;
[ ~,~,ted_K1,ted_K2 ] = LF(2,cp_zeta,cp_BT,1,k0,ted_kp);

% (6) Design the symbol timing Loop filter
st_BT = 0.01;
st_zeta = 1/sqrt(2);
ped_kp = 1;
[ ~,~,ped_K1,ped_K2 ] = LF(2,st_zeta,st_BT,1,k0,ped_kp);

% Make some input
in = xk + 1j*yk;

% DDS Initial Conditions
dds_sr = 1;
dds_sop = 0;
w0dds = 0;

% TED Initial Conditions
ted_sa_prev0 = 0;
ted_sa_prev_prev0 = 0;
ted_sa_prev1 = 0;
ted_sa_prev_prev1 = 0;
ted_lf_sd_prev0 = 0;

% PED Initial Conditions
ped_lf_skp = 0;

% Farrow Interpolator Initial Conditions
ff_prev0 = [ 0 0 0 ];
ff_prev1 = [ 0 0 0 ];

% Initialize Interpolation stuffs
k = 1;
mu(1) = 0;
strobe = 0;
eta = 0;
CNT_next = 0;

for ii = 2:numel(in)
    
    % Update Counter
    CNT = CNT_next;

    % CCW Rotation
    inp = in(ii)*dds_sr;
    
    % Get the rotated input ready for input to the farrow interpolators
    xp = real(inp);
    yp = imag(inp);
    
    % Run the interpolator for xp
    [ ff_out0,ff_prev0 ] = finterp(xp,mu(k),ff_prev0);
%     v2 = 1/2*[ 1 -1 -1 1]*real(in(ii:-1:ii-3));
%     v1 = 1/2*[ -1 3 -1 -1]*real(
    
    % Run the interpolator for yp
    [ ff_out1,ff_prev1 ] = finterp(yp,mu(k),ff_prev1);
    
    % These farrow outputs get passed to the TED,PED if the strobe is high
    if strobe
        
        % TED Stuff
        ted_sa0 = ff_out0;
        ted_sb0 = ted_sa_prev0;
        ted_sc0 = ted_sa_prev_prev0;
        ted_sd0 = sign(ted_sc0);
        ted_se0 = sign(ted_sa0);
        ted_sf0 = ted_sd0 - ted_se0;
        ted_sg0 = ted_sb0*ted_sf0;

        % update previous values
        ted_sa_prev_prev0 = ted_sa_prev0;
        ted_sa_prev0 = ted_sa0;

        ted_sa1 = ff_out1;
        ted_sb1 = ted_sa_prev1;
        ted_sc1 = ted_sa_prev_prev1;
        ted_sd1 = sign(ted_sc1);
        ted_se1 = sign(ted_sa1);
        ted_sf1 = ted_sd1 - ted_se1;
        ted_sg1 = ted_sb1*ted_sf1;

        % update previous values
        ted_sa_prev_prev1 = ted_sa_prev1;
        ted_sa_prev1 = ted_sa1;
        
        % TED Output
        ted_e(ii) = ted_sg0 + ted_sg1;
        
        % PED Stuff
        
        % Get our symbol estimates from the interpolation outputs
        a0(k) = A*sign(ff_out0);
        a1(k) = A*sign(ff_out1);
        
        % Save to make constellation
        xkk(k) = ff_out0;
        ykk(k) = ff_out1;
        
        if ~mod(k,5000)
            figure(1);
            plot(xkk,ykk,'kp');
            axis([ -3 3 -3 3 ]);
            title(sprintf('%d',k));
        end
        
        % ML Error
        ped_sf = a0(k)*ff_out1;
        ped_sg = a1(k)*ff_out0;
        ped_e(ii) = ped_sf - ped_sg;
        
    else
        % This is the hold part of the strobe business
        
        % TED,PED output 0
        ted_e(ii) = 0;
        ped_e(ii) = 0;
        
    end
    
    % PED Loop Filter
    ped_lf_si = ped_e(ii)*ped_K1;
    ped_lf_sj = ped_e(ii)*ped_K2;
    ped_lf_sl = ped_lf_skp; % previous sk
    ped_lf_sk = ped_lf_sj + ped_lf_sl;
    ped_lf_sm = ped_lf_si + ped_lf_sk;
    
    % Update
    ped_lf_skp = ped_lf_sk;
    
    % DDS Filter
    dds_sn = ped_lf_sm*k0;
    dds_sp = dds_sop; % previous so
    dds_so = dds_sn + w0dds + dds_sp;
    dds_sq = cos(dds_sp) + 1j*sin(dds_sp);
    dds_sr = conj(dds_sq);
    dds_out(ii) = dds_sq;
    
    % Update previous vals
    dds_sop = dds_so;
    
    % TED Loop Filter
    ted_lf_sa = ted_e(ii)*ted_K1;
    ted_lf_sb = ted_e(ii)*ted_K2;
    ted_lf_sc = ted_lf_sd_prev0;
    ted_lf_sd = ted_lf_sb + ted_lf_sc;
    ted_lf_output = ted_lf_sa + ted_lf_sd;
    
    % update previous values
    ted_lf_sd_prev0 = ted_lf_sd;
    
    % mod 1 counter
    W = ted_lf_output + 1/2;
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

%% Let's make some pretty plots
figure(1);
plot(ted_e);
title('TED Output');
xlabel('samples');
ylabel('TED e(t)');

figure(2);
plot(mu);
title('Fractional Interpolation Interval');
xlabel('k');
ylabel('\mu');

figure(3);
plot(ped_e);
title('PED Output');
xlabel('sample');
ylabel('PED  e(n)');

figure(4);
plot(real(dds_out(1:8e3)),'DisplayName','Real'); hold on;
plot(imag(dds_out(1:8e3)),'DisplayName','Imag');
title('DDS Output');
legend show;

%% Constellation
scatterplot(xkk(end-1e4:end) + 1j*ykk(end-1e4:end));
title('Constellation');

%% Differential Encoding
m0 = -a0; % Choose sign to get the right answer
m1 = -a1;  % Choose sign to get the right answer

d_hat0 = zeros(size(m0));
d_hat1 = zeros(size(m1));
d_hat0(m0 > 0) = 1;
d_hat1(m1 > 0) = 1;

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

% Find Unique Word
uw = ones(1,32)*3;
%uw = repmat([ 0 1 ],[ 1,16/2 ]);
idx = strfind(strjoin(string(b),''),strjoin(string(uw),''));
%disp(idx);
% idx = strfind(d_hat,strjoin(string(uw),''));

DataL = 3206/2;
%%
% try each possibility
for ii = 1:numel(idx)
    start = idx(end-ii+1) + numel(uw);
    o = 0;
    try
        s_hat = b(start+o:start+DataL-1+o);
        message{ii} = m2ascii(s_hat,M);
        fprintf('Message %d: %s\n',ii,message{ii});
    catch
    end
end