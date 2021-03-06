-- lua-npge, Nucleotide PanGenome explorer (Lua module)
-- Copyright (C) 2014-2016 Boris Nagaev
-- See the LICENSE file for terms of use.

describe("npge.util.fromFasta", function()
    it("parses fasta representation", function()
        local fromFasta = require 'npge.util.fromFasta'
        local fasta = [[
>foo descr
ATGC

>bar several words
AAA
TTT

>nodescription
ACC
]]
        local textToIt = require 'npge.util.textToIt'
        local lines = textToIt(fasta)
        local parser = fromFasta(lines)
        local foo_name, foo_descr, foo_text = parser()
        assert.equal(foo_name, "foo")
        assert.equal(foo_descr, "descr")
        assert.equal(foo_text, "ATGC")
        local bar_name, bar_descr, bar_text = parser()
        assert.equal(bar_name, "bar")
        assert.equal(bar_descr, "several words")
        assert.equal(bar_text, "AAATTT")
        local nod_name, nod_descr, nod_text = parser()
        assert.equal(nod_name, "nodescription")
        assert.equal(nod_descr, "")
        assert.equal(nod_text, "ACC")
    end)

    it("parses fasta representation from generator (newlines)",
    function()
        local itFromArray = require 'npge.util.itFromArray'
        local fromFasta = require 'npge.util.fromFasta'
        local parser = fromFasta(itFromArray {
            ">foo bar\n",
            "AAAA\n",
            "TTTT\n",
        })
        local foo, bar, text = parser()
        assert.equal(foo, "foo")
        assert.equal(bar, "bar")
        assert.equal(text, "AAAATTTT")
    end)

    it("parses fasta representation from generator (space after >)",
    function()
        local itFromArray = require 'npge.util.itFromArray'
        local fromFasta = require 'npge.util.fromFasta'
        local parser = fromFasta(itFromArray {
            "> foo bar\n",
            "AAAA\n",
            "TTTT\n",
        })
        local foo, bar, text = parser()
        assert.equal(foo, "foo")
        assert.equal(bar, "bar")
        assert.equal(text, "AAAATTTT")
    end)

    it("parses fasta representation from multiline generator",
    function()
        local itFromArray = require 'npge.util.itFromArray'
        local fromFasta = require 'npge.util.fromFasta'
        local parser = fromFasta(itFromArray {
            [[>foo bar
            AAAA
            TTTT]]
        })
        local foo, bar, text = parser()
        assert.equal(foo, "foo")
        assert.equal(bar, "bar")
        assert.equal(text, "AAAATTTT")
    end)
end)
