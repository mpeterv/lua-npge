-- lua-npge, Nucleotide PanGenome explorer (Lua module)
-- Copyright (C) 2014-2015 Boris Nagaev
-- See the LICENSE file for terms of use.

package = "lua-npge"
version = "dev-1"
source = {
    url = "git://github.com/npge/lua-npge.git"
}
description = {
    summary = "Nucleotide PanGenome explorer (Lua module)",
    homepage = "https://github.com/npge/lua-npge",
    license = "MIT",
}
dependencies = {
    "lua >= 5.1",
}
dependencies.platforms = {
    unix = {
        "alnbox",
        "luaposix",
    },
}
external_dependencies = {
    BOOST = {
        header = "boost/foreach.hpp"
    }
}
build = {
    type = "builtin",
    modules = {
        ['npge'] = 'src/npge/init.lua',
        ['npge.cpp'] = {
            sources = {
                "src/npge/cpp/lua_npge.cpp",
                "src/npge/cpp/model.cpp",
                "src/npge/cpp/throw_assert.cpp",
                "src/npge/cpp/strings.cpp",
                "src/npge/cpp/alignment.cpp",
                "src/npge/cpp/goodSlices.cpp",
                "src/npge/cpp/segmentTree.cpp",
            },
            libraries = {"stdc++"},
            incdirs = {"$(BOOST_INCDIR)"},
        },
        ['npge.config'] = 'src/npge/config.lua',
        ['npge.util'] = 'src/npge/util/init.lua',
        ['npge.util.arraysEqual'] = 'src/npge/util/arraysEqual.lua',
        ['npge.util.arraysLess'] = 'src/npge/util/arraysLess.lua',
        ['npge.util.asLines'] = 'src/npge/util/asLines.lua',
        ['npge.util.binary_search'] = 'src/npge/util/binary_search.lua',
        ['npge.util.clone'] = 'src/npge/util/clone.lua',
        ['npge.util.concatArrays'] = 'src/npge/util/concatArrays.lua',
        ['npge.util.endsWith'] = 'src/npge/util/endsWith.lua',
        ['npge.util.extractValue'] = 'src/npge/util/extractValue.lua',
        ['npge.util.fileExists'] = 'src/npge/util/fileExists.lua',
        ['npge.util.itFromArray'] = 'src/npge/util/itFromArray.lua',
        ['npge.util.loadstring'] = 'src/npge/util/loadstring.lua',
        ['npge.util.readFile'] = 'src/npge/util/readFile.lua',
        ['npge.util.readIt'] = 'src/npge/util/readIt.lua',
        ['npge.util.sandbox'] = 'src/npge/util/sandbox.lua',
        ['npge.util.split'] = 'src/npge/util/split.lua',
        ['npge.util.startsWith'] = 'src/npge/util/startsWith.lua',
        ['npge.util.timer'] = 'src/npge/util/timer.lua',
        ['npge.util.trim'] = 'src/npge/util/trim.lua',
        ['npge.util.unique'] = 'src/npge/util/unique.lua',
        ['npge.util.unpack'] = 'src/npge/util/unpack.lua',
        ['npge.util.writeIt'] = 'src/npge/util/writeIt.lua',
        ['npge.util.threads'] = 'src/npge/util/threads.lua',
        ['npge.util.mapItems'] = 'src/npge/util/mapItems.lua',
        ['npge.util.textToIt'] = 'src/npge/util/textToIt.lua',
        ['npge.model'] = 'src/npge/model/init.lua',
        ['npge.model.Block'] = 'src/npge/model/Block.lua',
        ['npge.model.BlockSet'] = 'src/npge/model/BlockSet.lua',
        ['npge.model.Fragment'] = 'src/npge/model/Fragment.lua',
        ['npge.model.Sequence'] = 'src/npge/model/Sequence.lua',
        ['npge.sequence'] = 'src/npge/sequence/init.lua',
        ['npge.sequence.fixPosition'] = 'src/npge/sequence/fixPosition.lua',
        ['npge.sequence.toFasta'] = 'src/npge/sequence/toFasta.lua',
        ['npge.fragment'] = 'src/npge/fragment/init.lua',
        ['npge.fragment.fragmentToSequence'] = 'src/npge/fragment/fragmentToSequence.lua',
        ['npge.fragment.hasPos'] = 'src/npge/fragment/hasPos.lua',
        ['npge.fragment.isSubfragmentOf'] = 'src/npge/fragment/isSubfragmentOf.lua',
        ['npge.fragment.reverse'] = 'src/npge/fragment/reverse.lua',
        ['npge.fragment.sequenceToFragment'] = 'src/npge/fragment/sequenceToFragment.lua',
        ['npge.fragment.subfragment'] = 'src/npge/fragment/subfragment.lua',
        ['npge.fragment.sub'] = 'src/npge/fragment/sub.lua',
        ['npge.block'] = 'src/npge/block/init.lua',
        ['npge.block.align'] = 'src/npge/block/align.lua',
        ['npge.block.consensus'] = 'src/npge/block/consensus.lua',
        ['npge.block.extend'] = 'src/npge/block/extend.lua',
        ['npge.block.goodSubblocks'] = 'src/npge/block/goodSubblocks.lua',
        ['npge.block.identity'] = 'src/npge/block/identity.lua',
        ['npge.block.isGood'] = 'src/npge/block/isGood.lua',
        ['npge.block.orient'] = 'src/npge/block/orient.lua',
        ['npge.block.reverse'] = 'src/npge/block/reverse.lua',
        ['npge.block.slice'] = 'src/npge/block/slice.lua',
        ['npge.block.unwind'] = 'src/npge/block/unwind.lua',
        ['npge.block.hasSelfOverlap'] = 'src/npge/block/hasSelfOverlap.lua',
        ['npge.block.excludeSelfOverlap'] = 'src/npge/block/excludeSelfOverlap.lua',
        ['npge.block.hasRepeats'] = 'src/npge/block/hasRepeats.lua',
        ['npge.block.genomes'] = 'src/npge/block/genomes.lua',
        ['npge.block.blockType'] = 'src/npge/block/blockType.lua',
        ['npge.block.giveName'] = 'src/npge/block/giveName.lua',
        ['npge.alignment'] = 'src/npge/alignment/init.lua',
        ['npge.alignment.alignRows'] = 'src/npge/alignment/alignRows.lua',
        ['npge.alignment.anchor'] = 'src/npge/alignment/anchor.lua',
        ['npge.alignment.complement'] = 'src/npge/alignment/complement.lua',
        ['npge.alignment.complementRows'] = 'src/npge/alignment/complementRows.lua',
        ['npge.alignment.identity'] = 'src/npge/alignment/identity.lua',
        ['npge.alignment.join'] = 'src/npge/alignment/join.lua',
        ['npge.alignment.left'] = 'src/npge/alignment/left.lua',
        ['npge.alignment.moveIdentical'] = 'src/npge/alignment/moveIdentical.lua',
        ['npge.alignment.toAtgcnAndGap'] = 'src/npge/alignment/toAtgcnAndGap.lua',
        ['npge.alignment.toAtgcn'] = 'src/npge/alignment/toAtgcn.lua',
        ['npge.alignment.unwindRow'] = 'src/npge/alignment/unwindRow.lua',
        ['npge.alignment.goodSlices'] = 'src/npge/alignment/goodSlices.lua',
        ['npge.algo'] = 'src/npge/algo/init.lua',
        ['npge.algo.AddGoodBlast'] = 'src/npge/algo/AddGoodBlast.lua',
        ['npge.algo.Align'] = 'src/npge/algo/Align.lua',
        ['npge.algo.BlastHits'] = 'src/npge/algo/BlastHits.lua',
        ['npge.algo.Blast'] = 'src/npge/algo/Blast.lua',
        ['npge.algo.BlockSetToLua'] = 'src/npge/algo/BlockSetToLua.lua',
        ['npge.algo.BlocksWithoutOverlaps'] = 'src/npge/algo/BlocksWithoutOverlaps.lua',
        ['npge.algo.ConsensusSequences'] = 'src/npge/algo/ConsensusSequences.lua',
        ['npge.algo.Cover'] = 'src/npge/algo/Cover.lua',
        ['npge.algo.FilterGoodBlocks'] = 'src/npge/algo/FilterGoodBlocks.lua',
        ['npge.algo.Genomes'] = 'src/npge/algo/Genomes.lua',
        ['npge.algo.GoodSubblocks'] = 'src/npge/algo/GoodSubblocks.lua',
        ['npge.algo.HasOverlap'] = 'src/npge/algo/HasOverlap.lua',
        ['npge.algo.HasSelfOverlap'] = 'src/npge/algo/HasSelfOverlap.lua',
        ['npge.algo.ExcludeSelfOverlap'] = 'src/npge/algo/ExcludeSelfOverlap.lua',
        ['npge.algo.Join'] = 'src/npge/algo/Join.lua',
        ['npge.algo.LoadFromLua'] = 'src/npge/algo/LoadFromLua.lua',
        ['npge.algo.Merge'] = 'src/npge/algo/Merge.lua',
        ['npge.algo.NonCovered'] = 'src/npge/algo/NonCovered.lua',
        ['npge.algo.Orient'] = 'src/npge/algo/Orient.lua',
        ['npge.algo.PangenomeMaker'] = 'src/npge/algo/PangenomeMaker.lua',
        ['npge.algo.PrimaryHits'] = 'src/npge/algo/PrimaryHits.lua',
        ['npge.algo.ReadFromBs'] = 'src/npge/algo/ReadFromBs.lua',
        ['npge.algo.ReadSequencesFromFasta'] = 'src/npge/algo/ReadSequencesFromFasta.lua',
        ['npge.algo.ReAlign'] = 'src/npge/algo/ReAlign.lua',
        ['npge.algo.UnwindBlocks'] = 'src/npge/algo/UnwindBlocks.lua',
        ['npge.algo.Workers'] = 'src/npge/algo/Workers.lua',
        ['npge.algo.WriteSequencesToFasta'] = 'src/npge/algo/WriteSequencesToFasta.lua',
        ['npge.algo.GiveNames'] = 'src/npge/algo/GiveNames.lua',
        ['npge.view'] = 'src/npge/view/init.lua',
        ['npge.view.BlockInConsole'] = 'src/npge/view/BlockInConsole.lua',
    },
    install = {
        bin = {
            "src/bin/MakePangenome",
        },
    },
}

