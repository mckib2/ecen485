%% Modulator
function [ r ] = modulator(s,lut,b,N)
    % Upsame the mappings by N
    x = upsample(cell2mat(values(lut.forward,num2cell(s))),N);
    
    % Pulse-shape filtering
    r = filter(b,1,x);
end