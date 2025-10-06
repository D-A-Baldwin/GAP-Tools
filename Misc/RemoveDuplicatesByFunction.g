RemoveDuplicatesByFunction := function(list, same_func)
	# This function takes as its inputs a list of objects and
	# a function that takes two inputs. Two elements x and y of 
	# the list are considered duplicates if same_func(x,y) = true
	# The output is the same list with all duplicates removed.
	
    local result, i, j, keep;
	
    result := [];
    for i in [1..Length(list)] do
        keep := true;
        for j in [1..Length(result)] do
            if same_func(list[i], result[j]) then
                keep := false;
                break;
            fi;
        od;
        if keep then
            Add(result, list[i]);
        fi;
    od;
    return result;
end;
