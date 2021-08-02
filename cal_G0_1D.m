function G0 = cal_G0_1D(z,w)
    %1D greens function, normal incidence, x-polarization
    dz = z(2) - z(1); 
    Nz = length(z);
    G_dz = @(dz) 1j * (w/2) * exp(1j * w * abs(dz));
    G0 = zeros(Nz,Nz);
    for i = 1:Nz
        dz_i = z(i) - z;
        G0(:,i) = G_dz(dz_i);
    end
    G0 = G0 * dz;

end