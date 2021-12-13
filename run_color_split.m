your_path = '' % path to code folder

% Split the data
color_split([your_path filesep 'EndoscopeCode' filesep 'Example_Input' filesep 'data'])

% Split the calibration file
color_split([your_path filesep 'EndoscopeCode' filesep 'Example_Input' filesep 'cali'])

% Split the PSF
color_split([your_path filesep 'EndoscopeCode' filesep 'Example_Input' filesep 'psf'])
