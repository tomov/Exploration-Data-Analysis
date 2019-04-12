function [formula_both, formula_orig] = get_formula(regressor, do_orth, mixed_effects, intercept)

% notice we don't include decoded regressors in random effects 
% that's b/c we don't want to overparameterize the model; we want to be maximally conservative

switch regressor
    case 'RU'
        if mixed_effects
            if do_orth
                formula_both = 'C ~ -1 + V + RU + VTU + decRU_orth + (-1 + V + RU + VTU|S)';
            else
                formula_both = 'C ~ -1 + V + RU + VTU + decRU + (-1 + V + RU + VTU|S)';
            end
            formula_orig = 'C ~ -1 + V + RU + VTU + (-1 + V + RU + VTU|S)';
        else
            if do_orth
                formula_both = 'C ~ -1 + V + RU + VTU + decRU_orth';
            else
                formula_both = 'C ~ -1 + V + RU + VTU + decRU';
            end
            formula_orig = 'C ~ -1 + V + RU + VTU';
        end

    case 'V'
        if mixed_effects
            if do_orth
                formula_both = 'C ~ -1 + V + RU + VTU + decV_orth + (-1 + V + RU + VTU|S)';
            else
                formula_both = 'C ~ -1 + V + RU + VTU + decV + (-1 + V + RU + VTU|S)';
            end
            formula_orig = 'C ~ -1 + V + RU + VTU + (-1 + V + RU + VTU|S)';
        else
            if do_orth
                formula_both = 'C ~ -1 + V + RU + VTU + decV_orth';
            else
                formula_both = 'C ~ -1 + V + RU + VTU + decV';
            end
            formula_orig = 'C ~ -1 + V + RU + VTU';
        end

    case 'TU'
        if mixed_effects
            if do_orth
                formula_both = 'C ~ -1 + V + RU + VTU + VdecTU_orth + (-1 + V + RU + VTU|S)';
            else
                formula_both = 'C ~ -1 + V + RU + VTU + VdecTU + (-1 + V + RU + VTU|S)';
            end
            formula_orig = 'C ~ -1 + V + RU + VTU + (-1 + V + RU + VTU|S)';
        else
            if do_orth
                formula_both = 'C ~ -1 + V + RU + VTU + VdecTU_orth';
            else
                formula_both = 'C ~ -1 + V + RU + VTU + VdecTU';
            end
            formula_orig = 'C ~ -1 + V + RU + VTU';
        end

    case 'DV'
        if mixed_effects
            if do_orth
                formula_both = 'C ~ -1 + V + RU + VTU + decDV_orth + (-1 + V + RU + VTU|S)';
            else
                formula_both = 'C ~ -1 + V + RU + VTU + decDV + (-1 + V + RU + VTU|S)';
            end
            formula_orig = 'C ~ -1 + V + RU + VTU + (-1 + V + RU + VTU|S)';
        else
            if do_orth
                formula_both = 'C ~ -1 + V + RU + VTU + decDV_orth';
            else
                formula_both = 'C ~ -1 + V + RU + VTU + decDV';
            end
            formula_orig = 'C ~ -1 + V + RU + VTU';
        end

    case 'both'
        if mixed_effects
            if do_orth
                formula_both = 'C ~ -1 + V + RU + VTU + decRU_orth + VdecTU_orth + (-1 + V + RU + VTU|S)';
            else
                formula_both = 'C ~ -1 + V + RU + VTU + decRU + VdecTU + (-1 + V + RU + VTU|S)';
            end
            formula_orig = 'C ~ -1 + V + RU + VTU + (-1 + V + RU + VTU|S)';
        else
            if do_orth
                formula_both = 'C ~ -1 + V + RU + VTU + decRU_orth + VdecTU_orth';
            else
                formula_both = 'C ~ -1 + V + RU + VTU + decRU + VdecTU';
            end
            formula_orig = 'C ~ -1 + V + RU + VTU';
        end

    case 'three'
        if mixed_effects
            if do_orth
                formula_both = 'C ~ -1 + V + RU + VTU + decV_orth + decRU_orth + VdecTU_orth + (-1 + V + RU + VTU|S)';
            else
                formula_both = 'C ~ -1 + V + RU + VTU + decV + decRU + VdecTU + (-1 + V + RU + VTU|S)';
            end
            formula_orig = 'C ~ -1 + V + RU + VTU + (-1 + V + RU + VTU|S)';
        else
            if do_orth
                formula_both = 'C ~ -1 + V + RU + VTU + decV_orth + decRU_orth + VdecTU_orth';
            else
                formula_both = 'C ~ -1 + V + RU + VTU + decV + decRU + VdecTU';
            end
            formula_orig = 'C ~ -1 + V + RU + VTU';
        end

    otherwise
        assert(false);
end


if intercept
    formula_both = strrep(formula_both, '-1', '1');
    formula_orig = strrep(formula_orig, '-1', '1');
end
