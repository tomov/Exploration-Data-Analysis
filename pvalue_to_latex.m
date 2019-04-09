function p_string = pvalue_to_latex(p)

    if p > 0.0001
        p_string = sprintf('p = %.4f', p);
    else
        log10p = log10(p);
        if log10p < -20
            log10p = -20;
        end
        p_string = sprintf('p < 10^{%.0f}', ceil(log10p));
    end
end
