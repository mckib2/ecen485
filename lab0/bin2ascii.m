function [ out ] = bin2ascii(in)
%     % chop input into 7 bit chunks for ASCII conversion
%     % consider the 'reshape' command
%     in = reshape(in,7,[]).';
% 
%     % Convert to String
%     % consider using num2str
%     in = num2str(in);
% 
%     % Run binary to decimal conversion
%     % consider using bin2dec
%     out = bin2dec(in);
% 
%     % Cast the result as a char
%     out = char(out).';

    % But we can do it in a oneliner:
    out = char(bin2dec(num2str(reshape(in,7,[]).'))).';
end