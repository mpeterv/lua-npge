local has_pos

has_pos = function(fragment, index)
    if not fragment:parted() then
        if fragment:ori() == 1 then
            return index >= fragment:start() and
                index <= fragment:stop()
        else
            return index <= fragment:start() and
                index >= fragment:stop()
        end
    else
        local a, b = fragment:parts()
        return has_pos(a, index) or has_pos(b, index)
    end
end

return has_pos
