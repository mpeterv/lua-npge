-- lua-npge, Nucleotide PanGenome explorer (Lua module)
-- Copyright (C) 2014-2016 Boris Nagaev
-- See the LICENSE file for terms of use.

return function(blockset)
    local alignLeft = require 'npge.block.alignLeft'
    local blocks = {}
    for block in blockset:iterBlocks() do
        table.insert(blocks, alignLeft(block))
    end
    local BlockSet = require 'npge.model.BlockSet'
    return BlockSet(blockset:sequences(), blocks)
end
