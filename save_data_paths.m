% Tara Urner
% Jia Lab
% This script should be run from RayOptics_Main and sets where outputs of
% the program should be written

save_paths = struct();
m_p = image_params.master_path;

%%%% PSF CALIBRATION %%%%

if image_params.PSF_cal == 1
    save_paths.psf_save = [m_p filesep 'psf_cal_out']; 
    mkdir(save_paths.psf_save)
end    

%%%% MAGNIFICATION CALIBRATION %%%%

if image_params.MAG_cal == 1
    save_paths.mag_save = [m_p filesep 'mag_cal_out'];
    mkdir(save_paths.mag_save)
elseif image_params.Mag_Analysis ==1
    save_paths.mag_save = [m_p filesep 'mag_cal_out'];
end

%%%% RECONSTRUCTION %%%%

save_paths.recon_save = [m_p filesep 'recon_out'];
if ~exist(save_paths.recon_save, 'dir')
    mkdir(save_paths.recon_save)
end    

image_params.save_paths = save_paths;