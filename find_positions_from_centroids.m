function p = find_positions_from_centroids(centroids)
    p.bc_centers = centroids;
    p.bc_rows = p.bc_centers(:,2);
    p.bc_cols = p.bc_centers(:,1);

    % find positions based on row and column values
    s = sortrows(p.bc_centers,2) % sorted on rows column
    top_lenses = sortrows(s(1:2,:),1);
    top_left = top_lenses(1,:);
    top_right = top_lenses(2,:);
    center_lenses = sortrows(s(3:5,:),1);
    center_left = center_lenses(1,:);
    center = center_lenses(2,:);
    center_right = center_lenses(3,:);
    bottom_lenses = sortrows(s(6:7,:),1);
    bottom_left = bottom_lenses(1,:);
    bottom_right = bottom_lenses(2,:);


    p.lens_positions = struct('top_left',top_left,'top_right',top_right,...
        'center_left',center_left,'center',center,'center_right',center_right,...
        'bottom_left',bottom_left,'bottom_right',bottom_right)
end