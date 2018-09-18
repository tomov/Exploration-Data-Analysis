function w = getEffects(results, fixed)

[w_f, names_f] = fixedEffects(results);
[w_r, names_r] = randomEffects(results);

if fixed 
    w = w_f;
    w1 = w(3);
    w2 = w(1);
    w3 = w(2);
else
    for subj = 1:size(w_r,1)/3
        w(subj,:) = w_f + w_r((subj - 1) * 3 + 1 : subj * 3);
    end
    w1 = w(:,3);
    w2 = w(:,1);
    w3 = w(:,2);
end

w = [w1 w2 w3]; % reorder them
