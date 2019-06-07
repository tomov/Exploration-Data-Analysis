function [masks, region] = get_masks_contra(glmodel, contrast, clusterFWEcorrect, extent, Num, sphere)

    % same as get_masks but contralateral ROIs
    %

    if ~exist('sphere', 'var')
        sphere = 10; % 10 mm radius by default
    end

    EXPT = exploration_expt();

    % get ROI masks
    switch contrast
        case 'badre'
            % clusters = masks from paper
            masks = badre_2012_create_masks(false, sphere);
            %masks = masks(1); % TODO use all masks
            assert(false, 'not supported');

        case 'dlpfc'
            % clusters = masks from paper
            masks = dlpfc_2012_create_masks(false, sphere);
            %masks = masks(1); % TODO use all masks
            assert(false, 'not supported');

        case 'tommy'
            % clusters = masks from paper
            masks = tommy_2017_create_masks(false, sphere);
            assert(false, 'not supported');

        case 'conj_3'
            % conjunction of RU (GLM 39), TU (GLM 40), and V (GLM 41)
            
            files = dir(fullfile('masks', 'ClusterMask_conj_3_*.nii'));
            for i = 1:length(files)
                masks{i} = fullfile('masks', files(i).name);
            end
            assert(false, 'not supported');

        case 'conj_36'
            % conjunction of RU (GLM 36), TU (GLM 36)
            
            files = dir(fullfile('masks', 'ClusterMask_conj_36_*.nii'));
            for i = 1:length(files)
                masks{i} = fullfile('masks', files(i).name);
            end
            assert(false, 'not supported');

        case 'conj_45'
            % conjunction of RU (GLM 36), TU (GLM 36)
            
            files = dir(fullfile('masks', 'ClusterMask_conj_45_*.nii'));
            for i = 1:length(files)
                masks{i} = fullfile('masks', files(i).name);
            end
            assert(false, 'not supported');

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
                direct = '+/-';

                % TODO HACK FIXME for DV
                if ismember(glmodel, [47, 29])
                    direct = '-';
                end
                % TODO HACK FIXME for RU
                if ismember(glmodel, [36]) && strcmp(contrast, 'RU') && clusterFWEcorrect
                    direct = '-';
                end

                [V, Y, C, CI, region, extent, stat, mni, cor, results_table] = ccnl_extract_clusters(EXPT, glmodel, contrast, p, direct, alpha, Dis, Num, clusterFWEcorrect, extent);

                [~, Vmask] = load_mask('masks/mask.nii');
                mni(:,1) = - mni(:,1); % flip horizontally
                cor = mni2cor(mni, Vmask.mat);

                r = sphere / 1.5; % 10 mm radius

                % create spherical masks around peak voxel of each cluster (intersected with cluster)
                %
                for c = 1:length(region)
                    masks{c} = sprintf('sphere_contra_glm%d_%s_%d_%d_%d_r=%dmm.nii', glmodel, replace(contrast, ' ', '_'), mni(c,1), mni(c,2), mni(c,3), round(r * 1.5));

                    % no intersect mask!
                    ccnl_create_spherical_mask(cor(c,1), cor(c,2), cor(c,3), r, masks{c});
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
