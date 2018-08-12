for i = 1:31
    results_VTURU{i} = getModelAll(data,i);
end 

save('behavioral_weights.mat', 'results_VTURU');
