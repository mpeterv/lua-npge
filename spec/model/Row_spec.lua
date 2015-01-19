local model = require 'npge.model'
local Row = model.Row

describe("model.row", function()
    it("throws on empty string", function()
        assert.has_error(function()
            Row("")
        end)
    end)

    it("throws if called Row()", function()
        assert.has_error(function()
            Row()
        end)
    end)

    it("throws if Row(x) where x is not string", function()
        assert.has_error(function()
            Row({})
        end)
    end)

    it("uses row A-T-G-C #A_T_G_C", function()
        local r = Row("A-T-G-C")
        assert.are.equal(r:type(), "Row")
        assert.are.equal(r:length(), 7)
        assert.are.equal(r:fragment_length(), 4)
        assert.are.equal(r:text(), 'N-N-N-N')
        -- block2fragment
        assert.are.equal(r:block2fragment(0), 0)
        assert.are.equal(r:block2fragment(1), -1)
        assert.are.equal(r:block2fragment(2), 1)
        assert.are.equal(r:block2fragment(3), -1)
        assert.are.equal(r:block2fragment(4), 2)
        assert.are.equal(r:block2fragment(5), -1)
        assert.are.equal(r:block2fragment(6), 3)
        assert.has_error(function()
            r:block2fragment(-1)
        end)
        assert.has_error(function()
            r:block2fragment(-100)
        end)
        assert.has_error(function()
            r:block2fragment(7)
        end)
        assert.has_error(function()
            r:block2fragment(100)
        end)
        -- fragment2block
        assert.are.equal(r:fragment2block(0), 0)
        assert.are.equal(r:fragment2block(1), 2)
        assert.are.equal(r:fragment2block(2), 4)
        assert.are.equal(r:fragment2block(3), 6)
        assert.has_error(function()
            r:fragment2block(-1)
        end)
        assert.has_error(function()
            r:fragment2block(-100)
        end)
        assert.has_error(function()
            r:fragment2block(4)
        end)
        assert.has_error(function()
            r:fragment2block(100)
        end)
    end)

    local check_row = function(text)
        local to_atgcn = require 'npge.alignment.to_atgcn'
        local to_atgcn_and_gap =
            require 'npge.alignment.to_atgcn_and_gap'
        text = to_atgcn_and_gap(text)
        return function()
            local r = Row(text)
            assert.are.equal(r:length(), #text)
            local Sequence = require 'npge.model.Sequence'
            local ungapped = to_atgcn(text)
            assert.are.equal(r:fragment_length(), #ungapped)
            local text1 = text:gsub('[^-]', 'N')
            assert.are.equal(r:text(), text1)
            assert.are.equal(r:text(ungapped), text)
            local fp = 0
            for bp = 0, #text - 1 do
                local char = text:sub(bp + 1, bp + 1)
                if char ~= '-' then
                    assert.are.equal(r:block2fragment(bp), fp)
                    assert.are.equal(r:fragment2block(fp), bp)
                    assert.are.equal(r:block2left(bp), fp)
                    assert.are.equal(r:block2right(bp), fp)
                    assert.are.equal(r:block2nearest(bp), fp)
                    fp = fp + 1
                else
                    assert.are.equal(r:block2fragment(bp), -1)
                    local left = r:block2left(bp)
                    local right = r:block2right(bp)
                    local nearest = r:block2nearest(bp)
                    if fp >= 1 then
                        assert.are.equal(left, fp - 1)
                    else
                        assert.are.equal(left, -1)
                    end
                    if fp < #ungapped then
                        assert.are.equal(right, fp)
                    else
                        assert.are.equal(right, -1)
                    end
                    local nearest = r:block2nearest(bp)
                    if left ~= -1 and right ~= -1 then
                        local bp_l = r:fragment2block(left)
                        local dist_l = bp - bp_l
                        local bp_r = r:fragment2block(right)
                        local dist_r = bp_r - bp
                        if dist_l == dist_r then
                            assert(nearest == left or
                                   nearest == right)
                        elseif dist_l < dist_r then
                            assert.equal(nearest, left)
                        elseif dist_r < dist_l then
                            assert.equal(nearest, right)
                        end
                    elseif left ~= -1 then
                        assert.equal(nearest, left)
                    elseif right ~= -1 then
                        assert.equal(nearest, right)
                    else
                        assert.equal(nearest, -1)
                    end
                end
            end
            assert.has_error(function()
                r:block2fragment(-1)
            end)
            assert.has_error(function()
                r:block2fragment(-100)
            end)
            assert.has_error(function()
                r:block2fragment(#text)
            end)
            assert.has_error(function()
                r:block2fragment(2 * #text)
            end)
            assert.has_error(function()
                r:fragment2block(-1)
            end)
            assert.has_error(function()
                r:fragment2block(-100)
            end)
            assert.has_error(function()
                r:fragment2block(#ungapped)
            end)
            assert.has_error(function()
                r:fragment2block(#ungapped * 2)
            end)
        end
    end

    it("uses row A-T-G-C (2)", check_row("A-T-G-C"))

    it("uses row #ATGC", check_row("ATGC"))

    it("uses row A #row_of_one_letter", check_row("A"))

    it("uses row --A- #gaps_before_and_after",
        check_row("--A-"))

    it("uses row A----AT #gaps_between", check_row("A----AT"))

    it("uses #long row",
        check_row([[
TCACCATTATACAGTTATGGTATGAACTGGGTCTTCAT-AA------AA-AAAAATATTT
TTTTTTGTTTATGCCATCATAGTTGTTCAATTATGCTAGTTT-----------GAATACC
GAGCAAGAGCCACGTGCTTGAAAATCTTGCAAGCACTTTGAGGGGGAGCATTTTGAAAGC
TTAAGTTTGACTCAATAACTGCGATGGTTGAGGGTAAT----------------TT-ATG
-ATATATGACTTGCTTTCATCAAGTATGTCGCGTGATTACTGAAGCTTTCTCTGCCCTGC
ATAATGACCTATAATTATTC-----CAAAAAGCTTACTC
    ]]))

    it("compares rows", function()
        assert.equal(Row('A-T-G'), Row('A-T-G'))
        assert.equal(Row('A-T-G'), Row('N-N-N'))
        assert.not_equal(Row('A-T-G'), Row('ATG'))
        assert.not_equal(Row('A-T-G'), Row('-A-T-G'))
        assert.not_equal(Row('A-T-G'), Row('A-T-G-'))
        assert.not_equal(Row('A-T-G-'), Row('A-T-G--'))
        assert.not_equal(Row('-A-T-G-'), Row('--A-T-G-'))
    end)

    it("throws if Row:text() is applied not to string",
    function()
        assert.has_error(function()
            Row('T'):text({})
        end)
    end)

    it("throws in Row:text(x), where x of wrong length",
    function()
        assert.has_error(function()
            Row('TT'):text('T')
        end)
    end)

    it("throws in Row:block2fragment()", function()
        assert.has_error(function()
            Row('TT'):block2fragment()
        end)
    end)

    it("throws in Row:block2fragment(not number)", function()
        assert.has_error(function()
            Row('TT'):block2fragment({})
        end)
    end)

    it("throws in Row:block2fragment(-1)", function()
        assert.has_error(function()
            Row('TT'):block2fragment(-1)
        end)
    end)

    it("throws in Row:block2fragment(row length)", function()
        assert.has_error(function()
            Row('TT'):block2fragment(2)
        end)
    end)

    it("throws in Row:block2left()", function()
        assert.has_error(function()
            Row('TT'):block2left()
        end)
    end)

    it("throws in Row:block2left(not number)", function()
        assert.has_error(function()
            Row('TT'):block2left({})
        end)
    end)

    it("throws in Row:block2left(-1)", function()
        assert.has_error(function()
            Row('TT'):block2left(-1)
        end)
    end)

    it("throws in Row:block2left(row length)", function()
        assert.has_error(function()
            Row('TT'):block2left(2)
        end)
    end)

    it("throws in Row:block2right()", function()
        assert.has_error(function()
            Row('TT'):block2right()
        end)
    end)

    it("throws in Row:block2right(not number)", function()
        assert.has_error(function()
            Row('TT'):block2right({})
        end)
    end)

    it("throws in Row:block2right(-1)", function()
        assert.has_error(function()
            Row('TT'):block2right(-1)
        end)
    end)

    it("throws in Row:block2right(row length)", function()
        assert.has_error(function()
            Row('TT'):block2right(2)
        end)
    end)

    it("throws in Row:block2nearest()", function()
        assert.has_error(function()
            Row('TT'):block2nearest()
        end)
    end)

    it("throws in Row:block2nearest(not number)", function()
        assert.has_error(function()
            Row('TT'):block2nearest({})
        end)
    end)

    it("throws in Row:block2nearest(-1)", function()
        assert.has_error(function()
            Row('TT'):block2nearest(-1)
        end)
    end)

    it("throws in Row:block2nearest(row length)", function()
        assert.has_error(function()
            Row('TT'):block2nearest(2)
        end)
    end)

    it("throws in Row:fragment2block()", function()
        assert.has_error(function()
            Row('TT'):fragment2block()
        end)
    end)

    it("throws in Row:fragment2block(not number)", function()
        assert.has_error(function()
            Row('TT'):fragment2block({})
        end)
    end)

    it("throws in Row:fragment2block(-1)", function()
        assert.has_error(function()
            Row('TT'):fragment2block(-1)
        end)
    end)

    it("throws in Row:fragment2block(fragment length)",
    function()
        assert.has_error(function()
            Row('T-T'):fragment2block(2)
        end)
    end)
end)