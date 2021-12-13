% Tara Urner
% Jia Lab
% This script should be called from RayOptics_Main and loads data from the 3D LF endoscope to be reconstructed
% based on a specified naming convention and file structure. The folder input
% format should be:
%
% "cali" : lens calibration image, i.e. bright image that shows outlines of
% the lenses
% I: This should be updated to check for an output file called
% "centers" if lens fitting has already been run for this dataset
%
% "data" : file/files to be reconstructed
% I: multiple data files will need to be supported, but in a smart way.
%
% "mag" : files for optional magnification calibration. Individual files
% rather than image stack.
%
% "psf" : files for optional shifts calibration or positions file for reconstruction. Individual files rather
% than image stack.
%
% "psf_cal_in": specifically added for the blue PSF, for tracking based on
% green or red high SNR data.
%
% "psf_cal_out" : where positions files written out by psf_calibration go.
% Loaded for reconstructions without reloading the PSF
%
% "mag_cal_out" : same as above but for magnification dataset

load_paths = struct();

m_p = image_params.master_path;

%%%% LENS CALIBRATION %%%%

%%% LENS CAL IN %%%

if image_params.lens_cal == 1
    cali_path = [m_p filesep 'lens_cal'];
    cd(cali_path) 
    f = dir('*.tiff');
    load_paths.cali_file = f.name;
    load_paths.cali_path = cali_path;
    [l1,~,~,t1] = islink(load_paths.cali_file);
    if l1
        load_paths.cali_file_link = t1;
    end
else
    cali_path = [m_p filesep 'cali'];
    cd(cali_path) 
    f = dir('*.tiff');
    load_paths.cali_path = cali_path;
    load_paths.cali_file = f.name; 
    [l1,~,~,t1] = islink(load_paths.cali_file);
    if l1
        load_paths.cali_file_link = t1;
    end
end

if length(f)>1
    error(['Too many image files in ' califolder ' Currently loading data corresponding to only one lens calibration file at a time is supported'])
end



%%%% DATA PATH %%%%
if isfield(image_params,'datapath')
    data_path = image_params.datapath;
else
    data_path = [m_p filesep 'data'];
end
cd(data_path)
f = [dir('*.tiff'); dir('*.tif')];
if multiple_data == 1
    load_paths.data_file = f;
    [l2,~,~,t2] = islink(load_paths.data_file(1).name);
    if l2
        load_paths.data_file_link = t2;
    end
elseif multiple_data ==0
    if length(f)>1
       [load_paths.data_file, ~] = uigetfile('*.tiff;*.tif', 'Please select which data file you would like to process');
    else
       load_paths.data_file = f.name; 
    end
    [l2,~,~,t2] = islink(load_paths.data_file);
    if l2
        load_paths.data_file_link = t2;
    end
end
load_paths.data_path = data_path;


%%%% PSF %%%%
if image_params.PSF_cal == 1
   load_paths.psf_path = [m_p filesep 'psf'];
else
   load_paths.psf_pos_path = [m_p filesep 'psf_cal_out'];
   load_paths.psf_pos_file = [load_paths.psf_pos_path filesep 'positionsxy.mat']; 
end
    
%%%% PSF IN %%%%

if image_params.PSF_cal_in
    load_paths.psf_cal_in_path = [m_p filesep 'psf_cal_in'];
    load_paths.psf_cal_in_file = [load_paths.psf_cal_in_path filesep 'positionsxy.mat'];
end

%%%% MAG %%%%
if image_params.MAG_cal == 1
    load_paths.mag_path = [m_p filesep 'mag'];
    load_paths.mag_pos_file = [load_paths.mag_path filesep 'rgb_shifts.mat'];
elseif image_params.Mag_Analysis == 1 && image_params.MAG_cal == 0
    load_paths.mag_pos_path = [m_p filesep 'mag_cal_out'];
    load_paths.mag_pos_file = [load_paths.mag_pos_path filesep 'positionsxyMAG.mat'];
end    

%%%% END FILE SELECTION %%%%

image_params.load_paths = load_paths;