%% Modulator
function [ r ] = modulator(s,lut,b,N)
    keys = cell(size(s));
    for ii = 1:numel(s)
        keys{ii} = s(ii);
    end
    x = upsample(cell2mat(values(lut.forward,keys)),N);
    r = filter(b,1,x);
end