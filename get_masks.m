function [masks, region] = get_masks(glmodel, contrast, clusterFWEcorrect, extent, Num)

    EXPT = exploration_expt();

    % get ROI masks
    switch contrast
        case 'badre'
            % clusters = masks from paper
            masks = badre_2012_create_masks(false);
            %masks = masks(1); % TODO use all masks

        case 'dlpfc'
            % clusters = masks from paper
            masks = dlpfc_2012_create_masks(false);
            %masks = masks(1); % TODO use all masks

        case 'tommy'
            % clusters = masks from paper
            masks = tommy_2017_create_masks(false);

        otherwise


            if endsWith(contrast, '.img') || endsWith(contrast, '.nii')
                % it's a path to a .img or .nii mask
                masks = {contrast};

            elseif iscell(contrast)

                % the "contrast" is actually the cell array of mask paths
                masks = contrast;

            else

                % it's an actual contrast
                %

                % group-level settings
                p = 0.001;
                alpha = 0.05;
                Dis = 20;
                if ~exist('Num', 'var')
                    Num = 1; % # peak voxels per cluster; default in bspmview is 3
                end
                direct = '+';

                % TODO HACK FIXME for DV
                if ismember(glmodel, [47])
                    direct = '-';
                end
                % TODO HACK FIXME for RU
                if ismember(glmodel, [36]) && strcmp(contrast, 'RU') && clusterFWEcorrect
                    direct = '-';
                end

                [V, Y, C, CI, region, extent, stat, mni, cor, results_table] = ccnl_extract_clusters(EXPT, glmodel, contrast, p, direct, alpha, Dis, Num, clusterFWEcorrect, extent);

                r = 10 / 1.5; % 10 mm radius

                % create spherical masks around peak voxel of each cluster (intersected with cluster)
                %
                for c = 1:length(region)
                    masks{c} = sprintf('sphere_glm%d_%s_%d_%d_%d_r=%dmm.nii', glmodel, replace(contrast, ' ', '_'), mni(c,1), mni(c,2), mni(c,3), round(r * 1.5));
                    cmask = CI == CI(cor(c,1), cor(c,2), cor(c,3));
                    ccnl_create_spherical_mask(cor(c,1), cor(c,2), cor(c,3), r, masks{c}, cmask);
                end
            end

    end



    if ~exist('region', 'var')
        for c = 1:length(masks)
            mask = masks{c};
            [~, masknames{c}, ~] = fileparts(mask);
            region{c,:} = masknames{c};
        end
    end
