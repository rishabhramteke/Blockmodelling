function [logVal] = myLog(val)

if val == 0
    logVal = 0;
else
    logVal = log2(val);
end

end