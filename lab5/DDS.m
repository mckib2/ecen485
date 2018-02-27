function [ b,a ] = DDS(k0)
    % Give filter parameters for the DDS
%     b = [ k0 2*pi*.01 ];
%     a = [ 1 -1 ];
    b = [ 0 k0 ];
    a = [ 1 -1 ];
end