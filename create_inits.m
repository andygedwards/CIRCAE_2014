function inits  = create_inits(S,cycle)
    load('statenames');
    if strcmp(cycle,'end')
        inits = cell2struct(num2cell(S{end}(end,:)),names.states,2);
    else
        inits = cell2struct(num2cell(S{cycle}(end,:)),names.states,2);
    end
end
