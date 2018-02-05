%% Demodulator
function [ s,a,x,xk ] = demodulator(r,b,lut,N)

    % Matched filtering, spit a zero out front to make everything work
    x = filter(fliplr(b),1,[ 0 r ]);

    % Downsample
    xk = downsample(x,N);
    xk = xk(2:end); % boot the leading zero we added for filtering
    
    % Make a decision
    d = cell2mat(keys(lut.reverse)); % amplitudes
    
    % Realize these are the indices where the Euclidean squared norm is
    % minimized
    [ ~,idx ] = min((d.' - xk).^2);
    
    % Now grab the coefficients a_hat using those indices
    a = d(idx);
    
    % Reverse look up table to get the recovered signal
    s = cell2mat(values(lut.reverse,num2cell(a)));
end