%% Demodulator
function [ s,a,x ] = demodulator(r,b,lut,N)

    % Matched filter
    x = filter(b,1,[0 r]);

    % Downsample
    xk = downsample(x,N);
    xk = xk(2:end);
    
    % Make a decision
    d = cell2mat(keys(lut.reverse));
    a = zeros(1,numel(xk));
    s = a;
    for ii = 1:numel(xk)
        test_sig = d - xk(ii);
        [ ~,k ] = min(abs(test_sig));
        a(ii) = d(k);
        s(ii) = lut.reverse(a(ii));
    end
end