function color_split(path)
% takes all files in "path" and splits them into Red Green and Blue
% directories with sudirectories named the same as the input directory.
cd(path)
d = dir("*.tiff");
p = strsplit(path,filesep);
dirname = p{end};
parent_dir = strjoin(p(1:end-1),filesep);
r_d = [parent_dir filesep 'red' filesep]
g_d = [parent_dir filesep 'green' filesep]
b_d = [parent_dir filesep 'blue' filesep]
mkdir([r_d dirname])
mkdir([g_d dirname])
mkdir([b_d dirname])

for i = 1 : length(d)
   im = imread([path filesep d(i).name]);
   shortname = strsplit(d(i).name,'.');
   im2 = WhiteBalanceRGBtoGray(im);
   r = im2(:,:,1);
   g = im2(:,:,2);
   b = im2(:,:,3);
   imwrite(r,[r_d dirname filesep shortname{1} '_red.tiff'],'compression','none')
   imwrite(g, [g_d dirname filesep shortname{1} '_green.tiff'], 'compression','none')
   imwrite(b, [b_d dirname filesep shortname{1} '_blue.tiff'], 'compression','none')
end    

fprintf('color_split_done');
end