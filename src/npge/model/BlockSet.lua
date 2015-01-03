
local BlockSet = {}
local BlockSet_mt = {}
local bs_mt = {}

BlockSet_mt.__index = BlockSet_mt
bs_mt.__index = bs_mt

local is_prepangenome = function(seq2fragments)
    for seq, fragments in pairs(seq2fragments) do
        local lengths_sum = 0
        local prev
        for _, fragment in ipairs(fragments) do
            lengths_sum = lengths_sum + fragment:length()
            if prev and prev:common(fragment) > 0 then
                return false
            end
            prev = fragment
        end
        if lengths_sum ~= seq:length() then
            return false
        end
    end
    return true
end

local parent_or_fragment = function(self, f)
    local parent = self._parent_of_parts[f]
    return parent or f
end

BlockSet_mt.__call = function(self, sequences, blocks)
    local name2seq = {}
    local seq2fragments = {}
    for _, sequence in ipairs(sequences) do
        assert(sequence:type() == 'Sequence')
        assert(not name2seq[sequence:name()])
        name2seq[sequence:name()] = sequence
        seq2fragments[sequence] = {}
    end
    local parent_of_parts = {}
    local block_by_fragment = {}
    for _, block in ipairs(blocks) do
        for fragment in block:iter_fragments() do
            block_by_fragment[fragment] = block
            local seq = fragment:sequence()
            local name = seq:name()
            assert(name2seq[name],
                ("No sequence with name %q"):format(seq:name()))
            if not fragment:parted() then
                table.insert(seq2fragments[seq], fragment)
            else
                local a, b = fragment:parts()
                table.insert(seq2fragments[seq], a)
                table.insert(seq2fragments[seq], b)
                parent_of_parts[a] = fragment
                parent_of_parts[b] = fragment
            end
        end
    end
    for seq, fragments in pairs(seq2fragments) do
        table.sort(fragments)
    end
    local prepangenome = is_prepangenome(seq2fragments)
    local bs = {_name2seq=name2seq, _blocks=blocks,
        _seq2fragments=seq2fragments,
        _parent_of_parts=parent_of_parts,
        _block_by_fragment=block_by_fragment,
        _prepangenome=prepangenome}
    return setmetatable(bs, bs_mt)
end

bs_mt.type = function(self)
    return "BlockSet"
end

bs_mt.size = function(self)
    return #(self._blocks)
end

bs_mt.is_prepangenome = function(self)
    return self._prepangenome
end

bs_mt.blocks = function(self)
    local clone = require 'npge.util.clone'
    return clone.array(self._blocks)
end

bs_mt.iter_blocks = function(self)
    local it, t, index = ipairs(self._blocks)
    return function()
        index, block = it(t, index)
        return block
    end
end

bs_mt.fragments = function(self, sequence)
    local clone = require 'npge.util.clone'
    return clone.array_from_it(self:iter_fragments(sequence))
end

bs_mt.iter_fragments = function(self, sequence)
    -- iterate pairs (fragment, sub-fragment)
    local fragments = self._seq2fragments[sequence]
    assert(fragments, "Sequence not in blockset")
    local it, t, index = ipairs(fragments)
    return function()
        index, fragment = it(t, index)
        local parent = parent_or_fragment(self, fragment)
        return parent, fragment
    end
end

bs_mt.sequences = function(self)
    local seqs = {}
    for name, seq in pairs(self._name2seq) do
        table.insert(seqs, seq)
    end
    return seqs
end

bs_mt.iter_sequences = function(self)
    local name, seq
    return function()
        name, seq = next(self._name2seq, name)
        return seq
    end
end

bs_mt.has_sequence = function(self, seq)
    return self._seq2fragments[seq] ~= nil
end

bs_mt.sequence_by_name = function(self, name)
    return self._name2seq[name]
end

bs_mt.block_by_fragment = function(self, fragment)
    return self._block_by_fragment[fragment]
end

bs_mt.overlapping_fragments = function(self, fragment)
    local concat_arrays = require 'npge.util.concat_arrays'
    local unique = require 'npge.util.unique'
    if fragment:parted() then
        local a, b = fragment:parts()
        return unique(concat_arrays(
            self:overlapping_fragments(a),
            self:overlapping_fragments(b)))
    end
    local seq = fragment:sequence()
    local fragments = self._seq2fragments[seq]
    assert(fragments, "Sequence not in blockset")
    local result = {}
    local add_fragment_or_parent = function(f)
        table.insert(result, parent_or_fragment(self, f))
    end
    local upper = require('npge.util.binary_search').upper
    local index = upper(fragments, fragment)
    for i = index, #fragments do
        if fragment:common(fragments[i]) > 0 then
            add_fragment_or_parent(fragments[i])
        else
            break
        end
    end
    for i = index - 1, 1, -1 do
        if fragment:common(fragments[i]) > 0 then
            add_fragment_or_parent(fragments[i])
        else
            break
        end
    end
    return unique(result)
end

bs_mt.next = function(self, fragment)
    if fragment:parted() then
        local a, b = fragment:parts()
        local f = (a < b) and a or b
        return self:next(f)
    end
    local seq = fragment:sequence()
    local fragments = self._seq2fragments[seq]
    assert(fragments, "Sequence not in blockset")
    local lower = require('npge.util.binary_search').lower
    local index = lower(fragments, fragment)
    assert(fragments[index] == fragment)
    local f
    if index < #fragments then
        f = fragments[index + 1]
    elseif index == #fragments and seq:circular() then
        f = fragments[1]
    else
        return nil
    end
    return parent_or_fragment(self, f)
end

bs_mt.prev = function(self, fragment)
    if fragment:parted() then
        local a, b = fragment:parts()
        local f = (a < b) and b or a
        return self:prev(f)
    end
    local seq = fragment:sequence()
    local fragments = self._seq2fragments[seq]
    assert(fragments, "Sequence not in blockset")
    local lower = require('npge.util.binary_search').lower
    local index = lower(fragments, fragment)
    assert(fragments[index] == fragment)
    local f
    if index > 1 then
        f = fragments[index - 1]
    elseif index == 1 and seq:circular() then
        f = fragments[#fragments]
    else
        return nil
    end
    return parent_or_fragment(self, f)
end

return setmetatable(BlockSet, BlockSet_mt)

