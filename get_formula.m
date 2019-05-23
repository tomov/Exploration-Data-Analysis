function [formula_both, formula_orig] = get_formula(regressor, do_orth, mixed_effects, intercept)

switch regressor
    case 'RU'
        if mixed_effects
            if do_orth
                formula_both = 'C ~ -1 + V + RU + VTU + decRU_orth + (-1 + V + RU + VTU + decRU_orth|S)';
            else
                formula_both = 'C ~ -1 + V + RU + VTU + decRU + (-1 + V + RU + VTU + decRU|S)';
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
                formula_both = 'C ~ -1 + V + RU + VTU + decV_orth + (-1 + V + RU + VTU + decV_orth|S)';
            else
                formula_both = 'C ~ -1 + V + RU + VTU + decV + (-1 + V + RU + VTU + decV|S)';
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
                formula_both = 'C ~ -1 + V + RU + VTU + VdecTU_orth + (-1 + V + RU + VTU + VdecTU_orth|S)';
            else
                formula_both = 'C ~ -1 + V + RU + VTU + VdecTU + (-1 + V + RU + VTU + VdecTU|S)';
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

    case 'VTU'
        if mixed_effects
            if do_orth
                formula_both = 'C ~ -1 + V + RU + VTU + decVTU_orth + (-1 + V + RU + VTU + decVTU_orth|S)';
            else
                formula_both = 'C ~ -1 + V + RU + VTU + decVTU + (-1 + V + RU + VTU + decVTU|S)';
            end
            formula_orig = 'C ~ -1 + V + RU + VTU + (-1 + V + RU + VTU|S)';
        else
            if do_orth
                formula_both = 'C ~ -1 + V + RU + VTU + decVTU_orth';
            else
                formula_both = 'C ~ -1 + V + RU + VTU + decVTU';
            end
            formula_orig = 'C ~ -1 + V + RU + VTU';
        end

    case 'DV'
        if mixed_effects
            if do_orth
                formula_both = 'C ~ -1 + V + RU + VTU + decDV_orth + (-1 + V + RU + VTU + decDV_orth|S)';
            else
                formula_both = 'C ~ -1 + V + RU + VTU + decDV + (-1 + V + RU + VTU + decDV|S)';
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
                formula_both = 'C ~ -1 + V + RU + VTU + decRU_orth + VdecTU_orth + (-1 + V + RU + VTU + decRU_orth + VdecTU_orth|S)';
            else
                formula_both = 'C ~ -1 + V + RU + VTU + decRU + VdecTU + (-1 + V + RU + VTU + decRU + VdecTU|S)';
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
                formula_both = 'C ~ -1 + V + RU + VTU + decV_orth + decRU_orth + VdecTU_orth + (-1 + V + RU + VTU + decV_orth + decRU_orth + VdecTU_orth|S)';
            else
                formula_both = 'C ~ -1 + V + RU + VTU + decV + decRU + VdecTU + (-1 + V + RU + VTU + decV + decRU + VdecTU|S)';
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
