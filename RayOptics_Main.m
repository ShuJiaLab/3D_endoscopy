%% Light Field Endoscope Image Reconstruction

% Tara Urner, Andrew Inman, Ben Lapid
% Fall 2019 - Spring 2020
% Jia Lab
close all
clear 


image_params = struct();

%% Runtime params
multiple_data = 0; % set to one and all data files in the 'data' folder will be reconstructed.
manual_data_path = 0; % sometimes having a 'data' folder is too restrictive


master_path = ''; % Path to your data here

if manual_data_path == 1
    data_path = [''];
end        

% Make sure to set this line to 1 when calibrating the PSF, 
% and 0 after it has been done once, 
% this step takes a while depending on the size of your images and PSF stack.
image_params.PSF_cal            = 1; % if PSF calibration step needs to be carried out
image_params.lens_cal           = 0; % if we want to use centroids and lens
image_params.verbose            = 0; % if enabled, additional information will be displayed during runtime
image_params.manual_lens_fit    = 0;
image_params.MAG_cal            = 0; % if MAG calibration step needs to be carried out
image_params.Mag_Analysis       = 0; % if magnification analysis should be performed
image_params.PSF_cal_in         = 0; % for tracking the blue PSF, a calibration dataset generated based on the better
image_params.recon_range        = [5:0.1:7];
image_params.gaussian_psf       = 0;
image_params.normshifts         = 0;

                                  % green or red psf is needed  
image_params.size_th            = 400000;
                                                             
                                 
% Manual thresholding
image_params.manual_thresh      = 1; %if auto lens detect fails, a manual threshold can be set
thresh = 10;

if image_params.PSF_cal == 1 && multiple_data == 1
    error('Running PSF calibration and multiple data at the same time is redundant and not recommended. To override this warning please comment out this line.')
end

if image_params.manual_thresh ==1
    image_params.thresh = thresh;
end


%% Calibration params
image_params.magnification      = 1.45;
image_params.K1                 = -.7*1e-6;
image_params.K2                 = -.7*1e-6;
image_params.K3                 = 0.1*1e-11;

%% Lens sizing params
image_params.Pixel_D            = 0.00185; %pixel diameter
image_params.lens_pitch         = 1.4;
image_params.GRIN_padding       = 100;
image_params.PSF_gaussbox       = 7;

%% Load Files

image_params.master_path = master_path
if manual_data_path == 1
    image_params.data_path = data_path
end

load_data_paths
save_data_paths

im_p = image_params; % for re-initialization for looped data processing
%%
for i=1:size(image_params.load_paths.data_file,1)

    %% Load datafile
    if multiple_data == 1
        image_params.datafile = image_params.load_paths.data_file(i).name;
        data = double(imread([image_params.load_paths.data_path filesep image_params.datafile]));
    else
        image_params.datafile = image_params.load_paths.data_file;
        data = double(imread([image_params.load_paths.data_path filesep image_params.datafile]));
    end
    
    %% Lens Detection and Image Rectification
    image_params = LDetRect(image_params);   

    
    %% Optional PSF Parametrization Step (Built in Barrel Distortion) %%
    if image_params.PSF_cal == 1
        [image_params,psf_overlay] = PSF_calibration(image_params);
        image_params = Shift_Caclulation(image_params,psf_overlay);
        
    %% Otherwise load PSF Parametrization file %%
    else
        load([image_params.load_paths.psf_pos_path filesep 'positionsxy.mat']);
        psf_overlay = imread([image_params.load_paths.psf_pos_path filesep 'psf_overlay.tif']);
        image_params.Positionsxy = Positionsxy;
        image_params = Shift_Caclulation(image_params,psf_overlay);
    end    

    
    %% Data Barrel Correction
    if image_params.manual_lens_fit == 1
        image_params = BarrelDistortionCorrectionDataManualShift(image_params,data);
    else
        image_params = BarrelDistortionCorrectionData(image_params,data);
    end
    
    %% Reconstruction %%
    [Recon, image_params] = RayOpticsRecon(image_params);
    
    if multiple_data == 1
        image_params = im_p;
        multiple_data = 1;
        clearvars -except image_params multiple_data im_p
        close all      
    end

end




