function w = getEffects(results, fixed, order)

if ~exist('order', 'var')
    order = [3 1 2]; % in results_VTURU
end

[w_f, names_f] = fixedEffects(results);
[w_r, names_r] = randomEffects(results);

if fixed 
    w = w_f;
    for i = 1:3
        if order(i)
            new_w(i) = w(order(i));
        else
            new_w(i) = 0;
        end
    end
else
    for subj = 1:size(w_r,1)/3
        w(subj,:) = w_f + w_r((subj - 1) * 3 + 1 : subj * 3);
    end
    for i = 1:3
        if order(i)
            new_w(:,i) = w(:,order(i));
        else
            new_w(:,i) = zeros(size(w(:,1)));
        end
    end
end

w = new_w;

