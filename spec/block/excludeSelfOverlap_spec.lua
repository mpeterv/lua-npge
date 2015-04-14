-- lua-npge, Nucleotide PanGenome explorer (Lua module)
-- Copyright (C) 2014-2015 Boris Nagaev
-- See the LICENSE file for terms of use.

describe("npge.block.excludeSelfOverlap", function()
    it("find non self-overlapping blocks", function()
        local model = require 'npge.model'
        local a = model.Sequence("a&chr1&c", "ATGC")
        local b = model.Sequence("b&chr1&c", "ATGC")
        local a_0_1_dir = model.Fragment(a, 0, 1, 1)
        local a_1_0_rev = model.Fragment(a, 1, 0, -1)
        local a_2_2_dir = model.Fragment(a, 2, 2, 1)
        local a_3_1_dir = model.Fragment(a, 3, 1, 1) -- parted
        local a_1_2_dir = model.Fragment(a, 1, 2, 1)
        local b_0_1_dir = model.Fragment(b, 0, 1, 1)
        local b_1_1_dir = model.Fragment(b, 1, 1, 1)
        local eso = require 'npge.block.excludeSelfOverlap'
        local Block = model.Block
        assert.same(eso(Block {a_0_1_dir, a_2_2_dir}), {
            Block {a_0_1_dir, a_2_2_dir},
        })
        assert.same(eso(Block {a_3_1_dir, a_2_2_dir}), {
            Block {a_3_1_dir, a_2_2_dir},
        })
        assert.same(eso(Block {a_0_1_dir, b_0_1_dir}), {
            Block {a_0_1_dir, b_0_1_dir},
        })
        assert.equal(type(eso(Block {a_3_1_dir, a_0_1_dir})),
            'table')
    end)
end)
