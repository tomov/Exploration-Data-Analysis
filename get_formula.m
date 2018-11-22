function [formula_both, formula_orig, formula_dec] = get_formula(regressor, do_orth, mixed_effects)

switch regressor
    case 'RU'
        if mixed_effects
            if do_orth
                formula_both = 'C ~ -1 + V + RU + VTU + decRU_orth + (-1 + V + RU + VTU + decRU_orth|S)';
            else
                formula_both = 'C ~ -1 + V + RU + VTU + decRU + (-1 + V + RU + VTU + decRU|S)';
            end
            formula_orig = 'C ~ -1 + V + RU + VTU + (-1 + V + RU + VTU|S)';
            formula_dec = 'C ~ -1 + V + decRU + VTU + (-1 + V + decRU + VTU|S)';
        else
            if do_orth
                formula_both = 'C ~ -1 + V + RU + VTU + decRU_orth';
            else
                formula_both = 'C ~ -1 + V + RU + VTU + decRU';
            end
            formula_orig = 'C ~ -1 + V + RU + VTU';
            formula_dec = 'C ~ -1 + V + decRU + VTU';
        end

    case 'TU'
        if mixed_effects
            if do_orth
                formula_both = 'C ~ -1 + V + RU + VTU + VdecTU_orth + (-1 + V + RU + VTU + VdecTU_orth|S)';
            else
                formula_both = 'C ~ -1 + V + RU + VTU + VdecTU + (-1 + V + RU + VTU + VdecTU|S)';
            end
            formula_orig = 'C ~ -1 + V + RU + VTU + (-1 + V + RU + VTU|S)';
            formula_dec = 'C ~ -1 + V + RU + VdecTU + (-1 + V + RU + VdecTU|S)';
        else
            if do_orth
                formula_both = 'C ~ -1 + V + RU + VTU + VdecTU_orth';
            else
                formula_both = 'C ~ -1 + V + RU + VTU + VdecTU';
            end
            formula_orig = 'C ~ -1 + V + RU + VTU';
            formula_dec = 'C ~ -1 + V + RU + VdecTU';
        end

    case 'both'
        if mixed_effects
            if do_orth
                formula_both = 'C ~ -1 + V + RU + VTU + decRU_orth + VdecTU_orth + (-1 + V + RU + VTU + decRU_orth + VdecTU_orth|S)';
            else
                formula_both = 'C ~ -1 + V + RU + VTU + decRU + VdecTU + (-1 + V + RU + VTU + decRU + VdecTU|S)';
            end
            formula_orig = 'C ~ -1 + V + RU + VTU + (-1 + V + RU + VTU|S)';
            formula_dec = 'C ~ -1 + V + decRU + VdecTU + (-1 + V + decRU + VdecTU|S)';
        else
            if do_orth
                formula_both = 'C ~ -1 + V + RU + VTU + decRU_orth + VdecTU_orth';
            else
                formula_both = 'C ~ -1 + V + RU + VTU + decRU + VdecTU';
            end
            formula_orig = 'C ~ -1 + V + RU + VTU';
            formula_dec = 'C ~ -1 + V + decRU + VdecTU';
        end

    otherwise
        assert(false);
end


