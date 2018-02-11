%% Lab 3a - Binary Phase Shift Keying
% Nicholas McKibben
% ECEn 485
% 2018-02-05

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

% s = [ 1 0 0 1 1 1 0 1 ];
N = 8;
M = 2;
E = 1; A = amp(E,M);
beta = 1;
span = 12;
sps = 8;
b = rcosdesign(beta,span,sps);
% b = ones(1,N)/sqrt(N);
Omega0 = 2*pi*1/4;

s = randi([ 0 (M-1) ],1,200);

keyss = {  0,1 };
vals = { -A,A };
LUT2 = LUT(keyss,vals);
n = 0:(N*length(s) - 1);

r = modulator(s,LUT2,b,N);

% Matched filtering, spit a zero out front to make everything work
x = filter(fliplr(b),1,[ zeros([ size(r,1) 1 ]) r ]);
xk = downsample(x(:,N:end).',N);

% Now make a decision
d = cell2mat(keys(LUT2.reverse)); % amplitudes
[ ~,idx ] = min((d - xk).^2,[],2);
a = d(idx);
s_hat = cell2mat(values(LUT2.reverse,num2cell(a)));

s_hat = s_hat(12:end);

% eyediagram(x,N,1,1);
% scatterplot(xk,N);

% Check to make sure everything is hunky dory
if ~isequal(s(1:numel(s_hat)),s_hat)
    fprintf('Something has gone terribly wrong...\n');
end

%% Exercise
load('bpskdata.mat');
r = bpskdata(2,:);
n = bpskdata(1,:);
% [ s,~,x,xk ] = demodulator(r.*cos(Omega0*n),b,LUT2,N);

% Matched filtering, spit a zero out front to make everything work
x = filter(fliplr(b),1,[ zeros([ size(r,1) 1 ]) r.*cos(Omega0*n) ]);
xk = downsample(x(:,N:end).',N);

% Now make a decision
d = cell2mat(keys(LUT2.reverse)); % amplitudes
[ ~,idx ] = min((d - xk).^2,[],2);
a = d(idx);
s = cell2mat(values(LUT2.reverse,num2cell(a)));

message = m2ascii(s(end-154:end-N),M);
disp(message);

eyediagram(x,N,1,1);
scatterplot(xk,N);