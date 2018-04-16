function [ farrow_output0, fin ] = finterp(xk,mu,init)

    % Grab initial Conditions
    ff_sa_prev0 = init(1);
    ff_sc_prev0 = init(2);
    ff_sf_prev0 = init(3);

    % Farrow Filter for a0 - cubic interpolation
    ff_sa0 = xk;
    ff_sb0 = ff_sa0*(1/6);
    ff_sc0 = ff_sa_prev0; % previous sa value
    ff_sd0 = ff_sc0*(-1/2);
    ff_se0 = ff_sb0 + ff_sd0;
    ff_sf0 = ff_sc_prev0; % previous sc value
    ff_sg0 = ff_sf0*(1/2);
    ff_sh0 = ff_se0 + ff_sg0;
    ff_si0 = ff_sf_prev0; % previous sf value
    ff_sj0 = ff_si0*(1/6);
    ff_sk0 = ff_sh0 + ff_sj0;
    ff_sl0 = mu;
    ff_sm0 = ff_sk0*ff_sl0;
    
    ff_sn0 = ff_sa_prev0;
    ff_so0 = ff_sn0*(-1/2);
    ff_sp0 = ff_sc_prev0;
    ff_sq0 = ff_so0 - ff_sp0;
    ff_sr0 = ff_sf_prev0;
    ff_ss0 = ff_sr0*(1/2);
    ff_st0 = ff_sq0 + ff_ss0;
    ff_su0 = ff_st0 + ff_sm0;
    
    ff_sv0 = ff_sl0*ff_su0;
    
    ff_sw0 = ff_sa0*(-1/6);
    ff_sx0 = ff_sa_prev0;
    ff_sy0 = ff_sx0 + ff_sw0;
    ff_sz0 = ff_sc_prev0;
    ff_saa0 = ff_sz0*(-1/2);
    ff_sac0 = ff_sy0 + ff_saa0;
    ff_sab0 = ff_sf_prev0;
    ff_sad0 = ff_sab0*(-1/3);
    ff_sae0 = ff_sac0 + ff_sad0;
    ff_sal0 = ff_sae0 + ff_sv0;
    ff_saf0 = ff_sl0*ff_sal0;
    
    ff_sag0 = ff_sa_prev0;
    ff_sah0 = ff_sc_prev0;
    farrow_output0 = ff_saf0 + ff_sah0;

    % Set final conditions
    fin = [ ff_sa0 ff_sc0 ff_sf0 ];
end