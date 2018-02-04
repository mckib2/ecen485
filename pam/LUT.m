function [ lut ] = LUT(keys,vals)
    % The forward look up table
    flut = containers.Map(keys,vals);
    % The reverse look up table
    rlut = containers.Map(vals,keys);
    
    lut.forward = flut;
    lut.reverse = rlut;
end