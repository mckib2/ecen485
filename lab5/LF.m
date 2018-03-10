function [ b,a,K1,K2 ] = LF(order,zeta,BT,N,k0,kp)
    % Generate loop filter coefficients.
    
    if order == 2
        den = 1 + (2*zeta/N)*(BT/(zeta + 1/(4*zeta))) + ...
            (BT/(N*(zeta + 1/(4*zeta))))^2;

        K1_num = (4*zeta/N)*(BT/(zeta + 1/(4*zeta)));
        K1 = K1_num/den;
        K1 = K1/(k0*kp);

        K2_num = (4/N^2)*(BT/(zeta + 1/(4*zeta)))^2;
        K2 = K2_num/den;
        K2 = K2/(k0*kp);

        b = [ (K1+K2) -K1 ];
        a = [ 1 -1 ];
    elseif order == 1
        b = 4*BT/(1 + 2*BT);
        b = b/(k0*kp);
        a = 1;
    end
end
