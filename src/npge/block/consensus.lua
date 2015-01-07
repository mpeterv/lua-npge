local ATGC = {'A', 'T', 'G', 'C'}

local consensus_at_pos = function(rows, bp)
    local counts = {}
    for _, row in ipairs(rows) do
        local letter = row:sub(bp + 1, bp + 1)
        counts[letter] = (counts[letter] or 0) + 1
    end
    local best_letter = 'N'
    local best_count = 0
    for _, letter in ipairs(ATGC) do
        local count = counts[letter] or 0
        if count > best_count then
            best_letter = letter
            best_count = count
        end
    end
    return best_letter
end

return function(block)
    local rows = {}
    for fragment in block:iter_fragments() do
        table.insert(rows, block:text(fragment))
    end
    local result = {}
    for bp = 0, block:length() - 1 do
        local letter = consensus_at_pos(rows, bp)
        table.insert(result, letter)
    end
    return table.concat(result)
end
