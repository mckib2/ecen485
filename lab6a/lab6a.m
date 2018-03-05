%% Lab 6a - BPSK Unique Word
% Nicholas McKibben
% ECEn 485
% 2018-02-27

clear;
close all;

% This is to get m2ascii function
if ~exist('lab0','dir')
    fprintf('Adding lab0 to path...\n');
    addpath('../lab0','-end');
end
% This is to get PLL functions
if ~exist('lab5','dir')
    fprintf('Adding lab5 to path...\n');
    addpath('../lab5','-end');
end

% Params
M = 2;
N = 8;
E = 1;
w0 = 2*pi*.2;
beta = .5;
span = 12;
uw = [ 0 0 0 1 0 1 1 0 ];
DataL = 224;
repeats = 4; % The packet is repeated to allow your PLL to acquire and lock
A = sqrt(E);
R = 1;
Fs = R*N;

% Grab the filter coefficients
b = rcosdesign(beta,span,N);

% PLL
order = 2;
k0 = 1; kp = 1;
BT = .01;
zeta = 1/sqrt(2);

% Loop filter
[ lf_b,lf_a ] = LF(order,zeta,BT,N,k0,kp);
lf_zi = zeros(1,max(length(lf_a),length(lf_b))-1);

% DDS
[ dds_b,dds_a ] = DDS(k0);
dds_zi = zeros(1,max(length(dds_a),length(dds_b))-1);

% Design the detector
load('bpskcruwdata.mat');
r_t  = bpskcruwdata(2,:);
to = bpskcruwdata(1,:);

% Make a local oscillator with which to mix
LOx = @(x) sqrt(2)*cos(w0*x);
LOy = @(x) -1*sqrt(2)*sin(w0*x);

% Now do the mixin'
Ir_t = r_t.*LOx(to);
Qr_t = r_t.*LOx(to);

% Now do the filterin'
fdelay = span/(2*R);
x_t = filter(fliplr(b),1,[ Ir_t zeros(1,fdelay*Fs) ]);
y_t = filter(fliplr(b),1,[ Qr_t zeros(1,fdelay*Fs) ]);

% Compensate for group delay
x_t = x_t(fdelay*Fs+1:end);
y_t = y_t(fdelay*Fs+1:end);

% Make sure everything's looking hunky dory
% eyediagram(x_t,N,1);
% title('x(t)');
% eyediagram(y_t,N,1);
% title('y(t)');
% eh, you can see it mostly

% Do the downsampling
xk = downsample(x_t,N);
yk = downsample(y_t,N);

% Let's try some visualization to get an idea of what we're up against
R = @(x,y,t) [ cos(t) -sin(t); sin(t) cos(t) ]*[ x; y ];
for k = 1:numel(xk)
    
    flipped = R(xk(k),yk(k),pi);
    len(k) = pdist([ xk(k) yk(k); flipped.' ]);
    
%     figure(1);
%     plot([ xk(k) flipped(1) ],[ yk(k) flipped(2) ],'k');
%     hold on; plot(xk(k),yk(k),'ko'); hold off;
%     axis([ -1.5 1.5 -1.5 1.5 ]);
%     drawnow;
end

zers = 17:52:numel(xk);
s = 1;
for k = 1:numel(xk)
    if ismember(k,zers)
        s = s*-1;
    end
    
    a(k) = s*sign(xk(k));
    
    len(k) = len(k)*s;
end

% Feed into PLL
e = zeros(1,numel(xk)); theta_hat = 0;
out = e; v = e; a = e; in = e; prev = 0; s = 1;
for k = 3:(numel(xk)-1)
    % CCW rotation on xk,yk to get xp,yp
%     tmp = R(xk(k),yk(k),angle(out(k)))*s;
    tmp = R(xk(k),yk(k),0)*s;
    xp = tmp(1);
    yp = tmp(2);
    
    % Now make a decision
    a(k) = sign(xp); % technically multiply by A, too, but A = 1
  
    % Now find the phase error
    % Heuristic
    %e(k) = atan2(yp,xp) - atan2(0,a(k));
    % ML
    e(k) = yp*a(k);

    % Flip signs where we see pi phase shift
    if (e(k-2) > e(k-1)) && (e(k-1) < e(k))
        a(k-1) = a(k-1)*-1;
        a(k) = a(k)*-1;
        s = s*-1;
    end

    % Can't seem to get this to work...
    [ v(k),lf_zi ] = filter(lf_b,lf_a,e(k),lf_zi); % Run through the loop filter
    [ theta_hat,dds_zi ] = filter(dds_b,dds_a,v(k) + w0,dds_zi); % Run through the DDS
    out(k+1) = exp(1j*(theta_hat)); % construct our output
end

% % For fun, let's look at the PLL plots
figure(4);
show(numel(xk)-1,numel(xk)-1,e,in,out,theta_hat);

%% Check for the unique word
m = a;
m(m < 0) = 0;
mstr = strjoin(string(m),'');
idx = strfind(mstr,strjoin(string(uw),''));

if isempty(idx)
    fprintf('Checking for pi phase..\n');
    m = -a;
    m(m < 0) = 0;
    mstr = strjoin(string(m),'');
    idx = strfind(mstr,strjoin(string(uw),''));
end

% try each possibility
for ii = 1:numel(idx)
    start = idx(end-ii+1) + numel(uw);
    o = 0;
    try
        m2 = m(start+o:start+DataL-1+o);
        fprintf('Message %d: %s\n',ii,m2ascii(m2,M));
    catch
    end
end
