function [ lut ] = LUT(keys,vals)
    % The LUT has two directions, unfortunately MATLAB's Map object only
    % goes one way, so we'll make a struct that contains both

    % The forward look up table
    flut = containers.Map(keys,vals);
    
    % The reverse look up table
    rlut = containers.Map(vals,keys);
    
    % Put it in a struct
    lut.forward = flut;
    lut.reverse = rlut;
end