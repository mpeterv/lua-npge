local model = require 'npge.model'

describe("model.fragment", function()
    it("creates fragment", function()
        local s = model.Sequence("test_name", "ATGC")
        local f = model.Fragment(s, 0, 3, 1)
        assert.are.equal(f:seq(), s)
        assert.are.equal(f:start(), 0)
        assert.are.equal(f:stop(), 3)
        assert.are.equal(f:ori(), 1)
    end)

    it("has type 'Fragment'", function()
        local s = model.Sequence("G&C&c", "ATGC")
        local f = model.Fragment(s, 1, 2, 1)
        assert.are.equal(f:type(), "Fragment")
    end)

    local Fragment = require 'npge.model.Fragment'
    local fragment_gen = function(seq, start, stop, ori)
        return function()
            return Fragment(seq, start, stop, ori)
        end
    end

    it("throws on bad fragment", function()
        -- linear
        local s = model.Sequence("genome&chromosome&l", "ATGC")
        assert.has.errors(fragment_gen(nil, 1, 2, 1))
        assert.has.errors(fragment_gen(s, 100, 110, 1))
        assert.has.errors(fragment_gen(s, 1, 2, 10))
        assert.has.errors(fragment_gen(s, 0, 4, 1))
        assert.has.errors(fragment_gen(s, 2, 1, 1))
        assert.has.errors(fragment_gen(s, 1, 2, -1))
    end)

    it("no throw on parted fragments on circular", function()
        -- circular
        local s = model.Sequence("genome&chromosome&c", "ATGC")
        assert.has_no.errors(fragment_gen(s, 2, 1, 1))
        assert.has_no.errors(fragment_gen(s, 1, 2, -1))
    end)

    it("is equal", function()
        local s = model.Sequence("genome&chromosome&c", "ATGC")
        local f = Fragment(s, 1, 2, 1)
        assert.are.equal(Fragment(s, 1, 2, 1), f)
        assert.are.equal(f, f)
        assert.are_not.equal(Fragment(s, 1, 2, -1), f)
    end)

    it("is less", function()
        local s1 = model.Sequence("ABC&chromosome&c", "ATGC")
        local s2 = model.Sequence("CDF&chromosome&c", "ATGC")
        assert(Fragment(s1, 0, 0, 1) < Fragment(s2, 0, 0, 1))
        assert(Fragment(s2, 0, 0, 1) > Fragment(s1, 0, 0, 1))
        assert(Fragment(s2, 0, 0, 1) >= Fragment(s1, 0, 0, 1))
        assert(Fragment(s1, 0, 0, 1) <= Fragment(s2, 0, 0, 1))
        assert(Fragment(s1, 0, 0, 1) <= Fragment(s1, 0, 0, 1))
        assert(Fragment(s1, 0, 0, 1) >= Fragment(s1, 0, 0, 1))
        assert(Fragment(s1, 0, 1, 1) < Fragment(s1, 1, 1, 1))
        assert(Fragment(s1, 0, 1, 1) < Fragment(s1, 0, 2, 1))
    end)

    it("throws in 'a < b' if a or a is parted", function()
        local s1 = model.Sequence("ABC&chromosome&c", "ATGC")
        local b
        assert.has_error(function()
            b = Fragment(s1, 0, 1, -1) < Fragment(s1, 0, 2, 1)
        end)
        assert.has_error(function()
            b = Fragment(s1, 0, 1, 1) < Fragment(s1, 2, 1, 1)
        end)
        assert.has_error(function()
            b = Fragment(s1, 0, 1, -1) < Fragment(s1, 2, 1, 1)
        end)
    end)

    it("has common positions with other fragment", function()
        local s1 = model.Sequence("ABC&chromosome&c", "ATGC")
        local F = Fragment
        assert.equal(F(s1, 0, 0, 1):common(F(s1, 0, 0, 1)), 1)
        assert.equal(F(s1, 0, 0, 1):common(F(s1, 0, 1, 1)), 1)
        assert.equal(F(s1, 0, 0, 1):common(F(s1, 1, 1, 1)), 0)
        assert.equal(F(s1, 1, 0, 1):common(F(s1, 2, 3, 1)), 2)
        assert.equal(F(s1, 0, 3, -1):common(F(s1, 0, 3, 1)), 2)
        assert.equal(F(s1, 0, 3, -1):common(F(s1, 2, 1, 1)), 2)
    end)

    it("gets id", function()
        local s = model.Sequence("G&C&c", "ATGC")
        local f = Fragment(s, 1, 2, 1)
        assert.are.equal(f:id(), "G&C&c_1_2_1")
    end)

    it("detects parted fragments", function()
        local s = model.Sequence("genome&chromosome&c", "ATGC")
        assert.is_true(Fragment(s, 1, 2, -1):parted())
        assert.is_false(Fragment(s, 1, 2, 1):parted())
        assert.is_true(Fragment(s, 2, 1, 1):parted())
        assert.is_false(Fragment(s, 1, 1, -1):parted())
    end)

    it("gets parts of parted fragment (positive)", function()
        local s = model.Sequence("genome&chromosome&c", "ATGC")
        local a, b = Fragment(s, 2, 0, 1):parts()
        assert.are.equal(a:start(), 2)
        assert.are.equal(a:stop(), 3)
        assert.are.equal(a:ori(), 1)
        assert.are.equal(b:start(), 0)
        assert.are.equal(b:stop(), 0)
        assert.are.equal(b:ori(), 1)
    end)

    it("gets parts of parted fragment (negative)", function()
        local s = model.Sequence("genome&chromosome&c", "ATGC")
        local a, b = Fragment(s, 1, 2, -1):parts()
        assert.are.equal(a:start(), 1)
        assert.are.equal(a:stop(), 0)
        assert.are.equal(a:ori(), -1)
        assert.are.equal(b:start(), 3)
        assert.are.equal(b:stop(), 2)
        assert.are.equal(b:ori(), -1)
    end)

    it("detects size of fragment", function()
        local s = model.Sequence("genome&chromosome&c", "ATGC")
        assert.are.equal(Fragment(s, 1, 2, -1):size(), 4)
        assert.are.equal(Fragment(s, 1, 3, -1):size(), 3)
        assert.are.equal(Fragment(s, 1, 2, 1):size(), 2)
        assert.are.equal(Fragment(s, 1, 3, 1):size(), 3)
        assert.are.equal(Fragment(s, 2, 1, 1):size(), 4)
        assert.are.equal(Fragment(s, 3, 1, 1):size(), 3)
        assert.are.equal(Fragment(s, 2, 1, -1):size(), 2)
        assert.are.equal(Fragment(s, 3, 1, -1):size(), 3)
    end)

    it("checks if fragment has sequence index", function()
        local s = model.Sequence("genome&chromosome&c", "ATGC")
        local f = Fragment(s, 3, 1, 1)
        assert.is.truthy(f:has(3))
        assert.is.truthy(f:has(0))
        assert.is.truthy(f:has(1))
        assert.is.falsy(f:has(2))
        local f = Fragment(s, 1, 3, -1)
        assert.is.truthy(f:has(3))
        assert.is.truthy(f:has(0))
        assert.is.truthy(f:has(1))
        assert.is.falsy(f:has(2))
    end)

    it("#is_subfragment", function()
        local s = model.Sequence("genome&chromosome&c", "ATGC")
        local f = Fragment(s, 3, 1, 1)
        assert.is.truthy(
            Fragment(s, 3, 1, 1):is_subfragment_of(f))
        assert.is.truthy(
            Fragment(s, 3, 0, 1):is_subfragment_of(f))
        assert.is.falsy(
            Fragment(s, 2, 2, 1):is_subfragment_of(f))
        assert.is.falsy(
            Fragment(s, 1, 3, 1):is_subfragment_of(f))
    end)

    it("#is_subfragment whole sequence", function()
        local s = model.Sequence("genome&chromosome&c", "ATGC")
        local f = Fragment(s, 3, 2, 1)
        assert.is.truthy(
            Fragment(s, 0, 3, 1):is_subfragment_of(f))
        assert.is.truthy(
            Fragment(s, 2, 1, 1):is_subfragment_of(f))
        assert.is.truthy(
            Fragment(s, 2, 0, 1):is_subfragment_of(f))
        assert.is.truthy(
            Fragment(s, 2, 0, -1):is_subfragment_of(f))
    end)

    it("#is_subfragment source parted", function()
        local s = model.Sequence("genome&chromosome&c", "ATGC")
        local f = Fragment(s, 3, 0, 1)
        assert.is.falsy(
            Fragment(s, 0, 3, 1):is_subfragment_of(f))
        assert.is.truthy(
            Fragment(s, 0, 3, -1):is_subfragment_of(f))
    end)

    it("#is_subfragment self parted", function()
        local s = model.Sequence("genome&chromosome&c", "ATGC")
        local f = Fragment(s, 1, 3, 1)
        assert.is.falsy(
            Fragment(s, 1, 3, -1):is_subfragment_of(f))
    end)

    it("gets text of fragment", function()
        local s = model.Sequence("genome&chromosome&c", "ATGC")
        assert.are.equal(Fragment(s, 0, 0, 1):text(), "A")
        assert.are.equal(Fragment(s, 0, 0, -1):text(), "T")
        assert.are.equal(Fragment(s, 0, 1, 1):text(), "AT")
        assert.are.equal(Fragment(s, 0, 1, -1):text(), "TGCA")
        assert.are.equal(Fragment(s, 1, 0, 1):text(), "TGCA")
    end)

    it("gets subfragment", function()
        local s = model.Sequence("genome&chromosome&c", "ATGC")
        local f = Fragment(s, 1, 2, 1)
        assert.are.equal(f:subfragment(0, 0, 1):text(), "T")
        assert.are.equal(f:subfragment(0, 0, -1):text(), "A")
        assert.are.equal(f:subfragment(0, 1, 1):text(), "TG")
        assert.are.equal(f:subfragment(1, 0, -1):text(), "CA")
        local f = Fragment(s, 2, 0, 1)
        assert.are.equal(f:subfragment(0, 0, 1):text(), "G")
        assert.are.equal(f:subfragment(0, 0, -1):text(), "C")
        assert.are.equal(f:subfragment(0, 2, 1):text(), "GCA")
        assert.are.equal(f:subfragment(2, 0, -1):text(), "TGC")
        local f = Fragment(s, 0, 2, -1)
        assert.are.equal(f:sub(1, 2, 1), "GC")
        assert.has_error(function()
            f:subfragment(1, 2, -1)
        end)
    end)

    it("gets char by index (positive)", function()
        local s = model.Sequence("genome&chromosome&c", "ATGC")
        local f = Fragment(s, 1, 2, 1)
        assert.are.equal(f:at(0), "T")
        assert.are.equal(f:at(1), "G")
    end)

    it("gets char by index (negative)", function()
        local s = model.Sequence("genome&chromosome&c", "ATGC")
        local f = Fragment(s, 1, 2, -1)
        assert.are.equal(f:at(0), "A")
        assert.are.equal(f:at(1), "T")
        assert.are.equal(f:at(2), "G")
        assert.are.equal(f:at(3), "C")
    end)
end)

