%% Lab 0
% Nicholas McKibben
% ECEn 485
% 2018-01-22

%% Test bin2ascii
test = '1000001'; % A
out = bin2ascii(test);
fprintf('bin2ascii: %s\n',out);

%% Test m2ascii
test = [ 6 0 7 0 5 4 3 ];
out = m2ascii(test,8);
fprintf('m2ascii: %s\n',out); % a b c