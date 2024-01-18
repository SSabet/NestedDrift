function imu = MU_inv(u,par)
    imu = u.^(-1/par.gamma);
end