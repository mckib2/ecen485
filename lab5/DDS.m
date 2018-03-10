function [ b,a ] = DDS(k0)
    % Give filter parameters for the DDS
    b = [ 0 k0 ];
    a = [ 1 -1 ];
end