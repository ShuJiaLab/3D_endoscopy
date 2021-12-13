clear
close all

fuse = struct()

%%
% Input paths to single color reconstructions
fuse.r_path = '' % Your path here: '/red/recon_out';
fuse.g_path = '' % Your path here: '/green/recon_out';
fuse.b_path = '' % Your path here:'/blue/recon_out';

% Save path
fuse.save_path = '' %Where you want to save'/RGB PSF/';
%%
mkdir(fuse.save_path);

% Find files in each path
cd(fuse.r_path)
fuse.reddir = dir('*.tif');
rf = {fuse.reddir.name};

cd(fuse.b_path)
fuse.bluedir = dir('*.tif');
bf = {fuse.bluedir.name};

cd(fuse.g_path)
fuse.greendir = dir('*.tif');
gf = {fuse.greendir.name};

c = {'R','G','B'};
f = {rf,gf,bf};
l = cellfun(@length, f);
longest = l==max(l);

if sum(longest) == 3
   color1 = f{1};
   color1_name = c{1};
   colors23 = [2,3];
else

    % Run based on the longest list of files
    color1 = f{longest};
    color1_name = c{longest};
    colors23 = find(~longest);
end
f_idx = 1;

for i=1:length(color1)
    
   disp(['Fusing ' num2str(i) ' of ' num2str(length(color1)) '...'])
   
   % Match colors based on filenames
   
   
   match1_f = color1{i};
   name = strsplit(color1{i},'_');
   id = strcat(name{3},'__',name{4})
   
   color2 = f{colors23(1)};
   color2_name = c{colors23(1)};
   match2 = strfind(color2,id);
   matches = find(not(cellfun('isempty',match2)));
   if isempty(find(cellfun('isempty',match2)==0))
       disp('File not found')
       f_idx = i-1;
       continue
   end
   match2_f = color2{matches}
   
   color3 = f{colors23(2)};
   color3_name = c{colors23(2)};
   match3 = strfind(color3,id);
   matches = find(not(cellfun('isempty',match3)));
   if isempty(find(cellfun('isempty',match2)==0))
       disp('File not found')
       f_idx = i-1;
       continue
   end
   match3_f = color3{matches}
   
   % Fuse colors back in RGB order
   colors_order = {color1_name, color2_name, color3_name}
   colors_files = {match1_f, match2_f, match3_f}
   
   t1 = strfind(colors_order,'R');
   R_idx = find(not(cellfun('isempty',t1)));
   R = colors_files{R_idx};
   
   t1 = strfind(colors_order,'G');
   G_idx = find(not(cellfun('isempty',t1)));
   G = colors_files{G_idx};  
   
   t1 = strfind(colors_order,'B');
   B_idx = find(not(cellfun('isempty',t1)));
   B = colors_files{B_idx};  
   
   %fuse(i).fileloadorder = 
   
   fuse.filenames{f_idx,1} = R; 
   fuse.filenames{f_idx,2} = G;
   fuse.filenames{f_idx,3} = B;
   
   f_idx = f_idx + 1;
end

% Determine images sizes (recons must all have been run with same psf)

% Red
cd(fuse.r_path)
m = dir('*.mat')
load([fuse.r_path filesep m(1).name])
size_recons
fuse.red_p = [re,cc,Recon_rows,Recon_cols]

% Green
cd(fuse.g_path)
m = dir('*.mat')
load([fuse.g_path filesep m(1).name])
size_recons
fuse.green_p = [re,cc,Recon_rows,Recon_cols]

% Blue
cd(fuse.b_path)
m = dir('*.mat')
load([fuse.b_path filesep m(1).name])
size_recons
fuse.blue_p = [re,cc,Recon_rows,Recon_cols]

% These parameters determine the overlap of the reconstructions.
% Everything is relative to the center of the image, but reconstructions
% can be different sizes depending on their shift values.
size_param = [fuse.red_p;fuse.green_p;fuse.blue_p;];
dim = size_param(size_param == max(max(size_param)));
% White balance and fuse

cd(fuse.save_path);

cf = fuse.filenames;

for i=1:size(cf,1)
    
    name = strsplit(cf{i,1},'.');
    tag = strsplit(name{1},'_');
    tag2 = strsplit(name{2},'_');
    filename = [strjoin(tag(1:4),'_') '_' strjoin(tag2(2:6),'_') '_fuse.tif'];
    filepath = [fuse.save_path filesep filename];
        
    t = imfinfo([fuse.r_path filesep cf{i,1}]);
    h = t.Height;
    w = t.Width;
    z = length(t);
    rgbImage = zeros(dim(1),dim(2));
    rgbImage = uint8(rgbImage);
    
    for k=1:length(t)
        if (k > 555)
            break;
        end
        disp(['Writing file ' num2str(k) ' of ' num2str(length(t)) ' in volume ' num2str(i) ' of ' num2str(length(cf))...
            ' to ' filename])

        % Red
        offsetr_start = (dim(1) - size_param(1,3))/2;
        offsetr_end = dim(1) - offsetr_start;
        rgbImage(19+offsetr_start:offsetr_end+18,19+offsetr_start:offsetr_end+18,1) = uint8(imread([fuse.r_path filesep cf{i,1}],k));
        % Green
        offsetg_start = (dim(1) - size_param(2,3))/2;
        offsetg_end = dim(1) - offsetg_start;
        rgbImage(1+offsetg_start:offsetg_end,1+offsetg_start:offsetg_end,2) = uint8(imread([fuse.g_path filesep cf{i,2}],k));
        % Blue
        offsetb_start = (dim(1) - size_param(3,3))/2;
        offsetb_end = dim(1) - offsetb_start;
        rgbImage(1+offsetb_start:offsetb_end,1+offsetb_start:offsetb_end,3) = uint8(imread([fuse.b_path filesep cf{i,3}],k));
        
        %BalancedRGB = WhiteBalanceRGBtoGray(rgbImage);
        
        imwrite(uint8(rgbImage),filepath,'compression','none','writemode','append');
    end
    

end