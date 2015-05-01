-- lua-npge, Nucleotide PanGenome explorer (Lua module)
-- Copyright (C) 2014-2015 Boris Nagaev
-- See the LICENSE file for terms of use.

return function(blockset, prefix)
    local HasOverlap = require 'npge.algo.HasOverlap'
    assert(not HasOverlap(blockset))
    prefix = prefix or ''
    local seq2block = {}
    local sequences = {}
    for i, block in ipairs(blockset:blocks()) do
        local text
        if block:size() == 1 then
            local Fragment = require 'npge.model.Fragment'
            local fragment = block:fragments()[1]
            local seq0 = fragment:sequence()
            local f =  Fragment(seq0, 0, seq0:length() - 1, 1)
            if f == fragment then
                -- the only fragment of this block
                -- covers whole sequence
                text = seq0:text()
            end
        end
        if not text then
            local consensus = require 'npge.block.consensus'
            text = consensus(block)
        end
        local name = ('%sconsensus%06d'):format(prefix, i)
        local Sequence = require 'npge.model.Sequence'
        local seq = Sequence(name, text)
        seq2block[seq] = block
        table.insert(sequences, seq)
    end
    local BlockSet = require 'npge.model.BlockSet'
    return BlockSet(sequences, {}), seq2block
end
