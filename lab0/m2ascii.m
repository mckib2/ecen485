function [ out ] = m2ascii(in,M)
%     % Calculate the number of bits
%     N = log2(M);
%     
%     % Run decimal to binary conversion
%     % Consider using dec2bin()
%     in = dec2bin(in,N);
%     
%     % Run bin2ascii
%     out = bin2ascii(in.');
    
    % But we can do it in a oneliner
    out = bin2ascii(dec2bin(in,log2(M)).');
end