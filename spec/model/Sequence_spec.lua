-- lua-npge, Nucleotide PanGenome explorer (Lua module)
-- Copyright (C) 2014-2016 Boris Nagaev
-- See the LICENSE file for terms of use.

local Sequence = require 'npge.model.Sequence'

describe("npge.model.Sequence", function()
    it("sequence creation", function()
        local s = Sequence("test_name", "ATGC",
            "test description")
        assert.are.equal(s:name(), "test_name")
        assert.are.equal(s:text(), "ATGC")
        assert.are.equal(s:description(), "test description")
    end)

    it("throws on empty sequence", function()
        assert.has_error(function()
            Sequence('name', '')
        end)
    end)

    it("throws on unnamed sequence", function()
        assert.has_error(function()
            Sequence('', 'ATG')
        end)
    end)

    it("sequence creation with no description", function()
        local s = Sequence("test_name", "ATGC")
        assert.are.equal(s:name(), "test_name")
        assert.are.equal(s:text(), "ATGC")
        assert.are.equal(s:description(), '')
    end)

    it("sequence with lower text", function()
        local s = Sequence("test_name", "atgc")
        assert.are.equal(s:text(), "ATGC")
    end)

    it("sequence with gaps", function()
        local s = Sequence("test_name", "AT-GC")
        assert.are.equal(s:text(), "ATGC")
    end)

    it("sequence with new lines", function()
        local s = Sequence("test_name", "AT\nGC")
        assert.are.equal(s:text(), "ATGC")
    end)

    it("sequence with N", function()
        local s = Sequence("test_name", "ATNGC")
        assert.are.equal(s:text(), "ATNGC")
    end)

    it("sequence with Cornish-Bowden notation", function()
        local s = Sequence("test_name", "ATYGC")
        assert.are.equal(s:text(), "ATNGC")
    end)

    it("sets genome, chromosome and circular", function()
        local s = Sequence("genome&chromosome&c", "AAA")
        assert.are.equal(s:genome(), "genome")
        assert.are.equal(s:chromosome(), "chromosome")
        assert.is.truthy(s:circular())
    end)

    it("sets genome, chromosome and linear", function()
        local s = Sequence("genome&chromosome&l", "AAA")
        assert.are.equal(s:genome(), "genome")
        assert.are.equal(s:chromosome(), "chromosome")
        assert.is.falsy(s:circular())
    end)

    it("doesn't accept empty genome", function()
        local s = Sequence("&chr&c", "AAA")
        assert.are.equal(s:genome(), nil)
        assert.are.equal(s:chromosome(), nil)
        assert.is.falsy(s:circular())
    end)

    it("doesn't accept empty chromosome", function()
        local s = Sequence("genome&&c", "AAA")
        assert.are.equal(s:genome(), nil)
        assert.are.equal(s:chromosome(), nil)
        assert.is.falsy(s:circular())
    end)

    it("doesn't accept empty circularity", function()
        local s = Sequence("genome&chr&", "AAA")
        assert.are.equal(s:genome(), nil)
        assert.are.equal(s:chromosome(), nil)
        assert.is.falsy(s:circular())
    end)

    it("sets bad name", function()
        local s = Sequence("genome&chromosome&l&1", "AAA")
        assert.are.equal(s:genome(), nil)
        assert.are.equal(s:chromosome(), nil)
        assert.is.falsy(s:circular())
    end)

    it("sets bad circularity", function()
        local s = Sequence("genome&chromosome&zzz", "AAA")
        assert.are.equal(s:genome(), "genome")
        assert.are.equal(s:chromosome(), "chromosome")
        assert.is.falsy(s:circular())
    end)

    it("gets substrings from sequence's text", function()
        local s = Sequence("test_name", "ATGC")
        assert.are.equal(s:sub(0, 1), 'AT')
    end)

    it("throws if sub called with bad arguments", function()
        local s = Sequence("test_name", "ATGC")
        assert.has_error(function()
            assert.are.equal(s:sub(-1, 1), 'AT')
        end)
        assert.has_error(function()
            assert.are.equal(s:sub(0, 4), 'AT')
        end)
        assert.has_error(function()
            assert.are.equal(s:sub(3, 2), 'AT')
        end)
    end)

    it("gets length of text", function()
        local s = Sequence("test_name", "ATGC")
        assert.are.equal(s:length(), 4)
    end)

    it("has type Sequence", function()
        local s = Sequence("test_name", "ATGC")
        assert.are.equal(s:type(), "Sequence")
    end)

    it("compares sequences", function()
        local s1 = Sequence("test_name", 'ATGC')
        local s2 = Sequence("test_name", 'ATGC')
        local s3 = Sequence("test_name", 'ATGC', 'd')
        local s4 = Sequence("test_name", 'AT')
        local s5 = Sequence("test_name", 'TTTT')
        local s6 = Sequence("test", 'TTTT')
        assert.equal(s1, s1)
        assert.equal(s1, s2)
        assert.equal(s1, s3)
        assert.equal(s1, s4)
        assert.equal(s1, s5)
        assert.not_equal(s1, s6)
    end)

    it("makes string representation of sequence", function()
        local s1 = Sequence("test_name", 'ATGC')
        assert.truthy(tostring(s1))
    end)
end)
