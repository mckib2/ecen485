%% Lab 3b
% Nicholas McKibben
% ECEn 485
% 2018-02-10

clear;
close all;

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

% Params
R = 1000; % data rate - just for funsies
N = 8;
Fs = R*N;
beta = 0.5;
span = 12;
b = rcosdesign(beta,span,N);
b = b/max(abs(b)); % normalize filter coefficients
w0 = 2*pi*1/4;

% Find A for QPSK
E = 9;
M = 4;
A = sqrt((3*E)/(2*(M - 1)));

% Define our local oscillators
LOx = @(x) sqrt(2)*cos(w0*x);
LOy = @(x) -1*sqrt(2)*sin(w0*x);

%% Test signal time
i_t = [ A A -A -A ]/sqrt(2);
q_t = [ A -A A -A ]/sqrt(2);

% For the polyphase filter, append span/2 zeros at the end to flush all
% the useful samples out of the filter
i_t = [ i_t zeros(1,span/2) ];
q_t = [ q_t zeros(1,span/2) ];

% Now do the upsampling
i_t = upsample(i_t,N);
q_t = upsample(q_t,N);

% Now we can filter
I_t = filter(b,1,i_t);
Q_t = filter(b,1,q_t);

% Compensate for the raised cosine filter group delay by delaying the
% input signal
fdelay = span/(2*R);
I_t = I_t(fdelay*Fs+1:end);
Q_t = Q_t(fdelay*Fs+1:end);

% For our test signal, it was 4 long - mix in the LO
to = 1000*(0:numel([ 0 0 0 0 ])*N - 1)/Fs;
I_t = I_t.*LOx(to);
Q_t = Q_t.*LOy(to);

% Now add 'em together and you've got a modulated signal
s_t = I_t + Q_t;

%% Now for the detector part
% r_t = s_t; % test signal

% Load in the no-kidding real data
load('qpskdata.mat');
r_t = qpskdata(2,:);
DataL = 49; % the message is this long
to = qpskdata(1,:); % new time vector

% Mix it up
Ir_t = r_t.*LOx(to);
Qr_t = r_t.*LOy(to);

% Filter - make sure to flush out all the useful samples at the end
x_t = filter(fliplr(b),1, [ Ir_t zeros(1,fdelay*Fs) ]);
y_t = filter(fliplr(b),1, [ Qr_t zeros(1,fdelay*Fs) ]);

% Compensate for group delay
x_t = x_t(fdelay*Fs+1:end);
y_t = y_t(fdelay*Fs+1:end);

% Downsample
xk = downsample(x_t,N);
yk = downsample(y_t,N);

% And finally pick out the samples we were looking for
xk = xk(N:DataL+N-1);
yk = yk(N:DataL+N-1);

% Now do some deciding - all we need to do is check signs!
xyk = double(sign([ xk; yk ]));

keys = { '11', '-11', '1-1', '-1-1' };
vals = { 0,1,2,3 };
lut = LUT(keys,vals);

s_hat = zeros(1,size(xyk,2));
for nn = 1:size(xyk,2)
    s_hat(nn) = lut.forward(char(strjoin(string(xyk(:,nn).'),'')));
end

% Display the message
m = m2ascii(s_hat,M);
disp(m);

%% Get the eye diagram and constellation
eyediagram(x_t,N,1);
title('x(t)');
eyediagram(y_t,N,1);
title('y(t)');

figure(3);
plot(xk,yk,'.'); hold on;
D = 5.25;
plot([ D D -D -D ],[ D -D D -D ],'+');
title('Costellation, QPSK');
xlabel('\phi_0(t)');
ylabel('\phi_1(t)');