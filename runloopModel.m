for i = 1:20
    new_data{i} = getModelAll(data,i);
end 

save new_data 