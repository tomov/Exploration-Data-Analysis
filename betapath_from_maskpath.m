function betapath = betapath_from_maskpath(mask, subj, suffix)

    if ~exist('suffix', 'var')
        suffix = '';
    else
        if ~startsWith(suffix, '_')
            suffix = ['_', suffix];
        end
    end

    [~, maskname, ~] = fileparts(mask);
    betapath = ['betas_', maskname, '_subj=', num2str(subj), suffix, '.mat'];
end
