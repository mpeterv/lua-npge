-- lua-npge, Nucleotide PanGenome explorer (Lua module)
-- Copyright (C) 2014-2015 Boris Nagaev
-- See the LICENSE file for terms of use.

return function(a, b)
    if #a ~= #b then
        return false
    end
    local n = #a
    for i = 1, n do
        if a[i] ~= b[i] then
            return false
        end
    end
    return true
end