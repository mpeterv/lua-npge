describe("algo.Merge", function()
    it("merges two blocksets", function()
        local model = require 'npge.model'
        local s = model.Sequence("s", "ATAT")
        local f1 = model.Fragment(s, 0, 0, 1)
        local b1 = model.Block({f1})
        local bs1 = model.BlockSet({s}, {b1})
        local f2 = model.Fragment(s, 1, 1, 1)
        local b2 = model.Block({f2})
        local bs2 = model.BlockSet({s}, {b2})
        local Merge = require 'npge.algo.Merge'
        local sum = Merge(bs1, bs2)
        assert.truthy(sum:same_sequences(bs1))
        assert.same(sum:fragments(s), {f1, f2})
        assert.equal(sum:size(), 2)
    end)

    it("merges 3 blocksets", function()
        local model = require 'npge.model'
        local s = model.Sequence("s", "ATAT")
        local f1 = model.Fragment(s, 0, 0, 1)
        local b1 = model.Block({f1})
        local bs1 = model.BlockSet({s}, {b1})
        local f2 = model.Fragment(s, 1, 1, 1)
        local b2 = model.Block({f2})
        local bs2 = model.BlockSet({s}, {b2})
        local f3 = model.Fragment(s, 2, 2, 1)
        local b3 = model.Block({f3})
        local bs3 = model.BlockSet({s}, {b3})
        local Merge = require 'npge.algo.Merge'
        local sum = Merge(bs1, bs2, bs3)
        assert.equal(sum, model.BlockSet({s}, {b1, b2, b3}))
    end)

    it("merges 1 blocksets (return the argument)", function()
        local model = require 'npge.model'
        local s = model.Sequence("s", "ATAT")
        local f1 = model.Fragment(s, 0, 0, 1)
        local b1 = model.Block({f1})
        local bs1 = model.BlockSet({s}, {b1})
        local Merge = require 'npge.algo.Merge'
        local sum = Merge(bs1)
        assert.equal(sum, model.BlockSet({s}, {b1}))
    end)

    it("throws if called with no arguments", function()
        assert.has_error(function()
            local Merge = require 'npge.algo.Merge'
            local blockset = Merge()
        end)
    end)
end)