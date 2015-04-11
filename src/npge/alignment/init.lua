-- lua-npge, Nucleotide PanGenome explorer (Lua module)
-- Copyright (C) 2014-2015 Boris Nagaev
-- See the LICENSE file for terms of use.

return {
    left = require 'npge.alignment.left',
    anchor = require 'npge.alignment.anchor',
    toAtgcn = require 'npge.alignment.toAtgcn',
    toAtgcnAndGap = require 'npge.alignment.toAtgcnAndGap',
    unwindRow = require 'npge.alignment.unwindRow',
    complement = require 'npge.alignment.complement',
    complementRows = require 'npge.alignment.complementRows',
    moveIdentical = require 'npge.alignment.moveIdentical',
    join = require 'npge.alignment.join',
    alignRows = require 'npge.alignment.alignRows',
    identity = require 'npge.alignment.identity',
}
