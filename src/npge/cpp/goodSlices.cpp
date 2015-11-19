/* lua-npge, Nucleotide PanGenome explorer (Lua module)
 * Copyright (C) 2014-2015 Boris Nagaev
 * See the LICENSE file for terms of use.
 */

#include <boost/foreach.hpp>

#include "npge.hpp"
#include "throw_assert.hpp"

namespace lnpge {

static int ssLength(const StartStop& ss) {
    return ss.second - ss.first + 1;
}

class GoodSlicer {
private:
    // prefix sums
    std::vector<int> score_sum_;

public:
    GoodSlicer(const Scores& score) {
        score_sum_.resize(score.size() + 1);
        score_sum_[0] = 0;
        for (int i = 0; i < score.size(); i++) {
            score_sum_[i + 1] = score_sum_[i] + score[i];
        }
    }

    int blockLength() const {
        return score_sum_.size() - 1;
    }

    int maxPos() const {
        return blockLength() - 1;
    }

    int countScore(int start, int stop) const {
        ASSERT_LTE(0, start);
        ASSERT_LTE(start, stop);
        ASSERT_LT(stop, blockLength());
        return score_sum_[stop + 1] - score_sum_[start - 1 + 1];
    }

    int countScore(const StartStop& region) const {
        return countScore(region.first, region.second);
    }

    bool testRegion(
        const StartStop& region,
        int min_identity
    ) const {
        int min_score = ssLength(region) * min_identity;
        return countScore(region) >= min_score;
    }

    // moves start, stop to inside to first MIN_END 100% col
    // returns if succeeded
    bool fixEnds(StartStop& region, int min_end) const {
        int& start = region.first;
        while (start + min_end <= region.second
               && !testRegion(
                    StartStop(start, start + min_end - 1),
               MAX_COLUMN_SCORE)) {
            start += 1;
        }
        int& stop = region.second;
        while (region.first + min_end <= stop
               && !testRegion(
                    StartStop(stop - min_end + 1, stop),
               MAX_COLUMN_SCORE)) {
            stop -= 1;
        }
        return
            testRegion(
                StartStop(start, start + min_end - 1),
                MAX_COLUMN_SCORE
            ) &&
            testRegion(
                StartStop(stop - min_end + 1, stop),
                MAX_COLUMN_SCORE
            );
    }

    // writes all continous good regions (score 100)
    void makeSimpleRegions(Coordinates& output) const {
        bool in_region = false;
        for (int i = 0; i < blockLength(); i++) {
            if (countScore(i, i) == MAX_COLUMN_SCORE) {
                if (!in_region) {
                    output.push_back(StartStop(i, i));
                    in_region = true;
                } else {
                    ASSERT_GT(output.size(), 0);
                    output.back().second = i;
                }
            } else {
                in_region = false;
            }
        }
    }

    void complement(
        Coordinates& output,
        const Coordinates& input
    ) const {
        int last_stop = -1;
        BOOST_FOREACH (const StartStop& region, input) {
            if (last_stop + 1 < region.first) {
                StartStop m(last_stop + 1, region.first - 1);
                output.push_back(m);
            }
            last_stop = region.second;
        }
        if (input.empty()) {
            output.push_back(StartStop(0, maxPos()));
        } else if (input.back().second < maxPos()) {
            StartStop last(input.back().second + 1, maxPos());
            output.push_back(last);
        }
    }

    void findLongRegions(
        Coordinates& output,
        const Coordinates& input,
        int frame
    ) const {
        BOOST_FOREACH (const StartStop& region, input) {
            if (ssLength(region) >= frame) {
                output.push_back(region);
            }
        }
    }

    void mergeContinuous(
        Coordinates& output,
        const Coordinates& input
    ) const {
        BOOST_FOREACH (const StartStop& region, input) {
            if (!output.empty() &&
                    output.back().second + 1 == region.first) {
                output.back().second = region.second;
            } else {
                output.push_back(region);
            }
        }
    }

    void findGood(
        Coordinates& output,
        const Coordinates& input,
        int frame,
        int min_identity, // 0-100
        int min_end
    ) const {
        BOOST_FOREACH (const StartStop& region, input) {
            StartStop t = region;
            t.first = std::max(0, t.first - frame);
            t.second = std::min(maxPos(), t.second + frame);
            if (!fixEnds(t, min_end)) {
                continue;
            }
            if (!testRegion(t, min_identity)) {
                continue;
            }
            output.push_back(region);
        }
    }

    void joinRegions(
        Coordinates& output,
        const Coordinates& input,
        int frame,
        int min_identity, // 0-100
        int min_end
    ) const {
        Coordinates good;
        findLongRegions(good, input, frame);
        Coordinates other;
        complement(other, good);
        findGood(good, other, frame, min_identity, min_end);
        std::sort(good.begin(), good.end());
        mergeContinuous(output, good);
    }

    void buildGoodRegions(
        Coordinates& output,
        int frame1,
        int frame2, // min_length
        int min_identity, // 0-100
        int min_end
    ) const {
        frame1 = std::min(frame1, min_end);
        frame2 = std::min(frame2, frame1);
        Coordinates s1, s2;
        makeSimpleRegions(s1);
        joinRegions(s2, s1, frame1, min_identity, min_end);
        joinRegions(output, s2, frame2, min_identity, min_end);
    }

    // see https://github.com/npge/npge/issues/26
    bool isGoodBlock(
        const Coordinates& good_regions,
        int min_length,
        int min_identity,
        int min_end
    ) const {
        // Asserts: are true if built with buildGoodRegions
        BOOST_FOREACH (const StartStop& region, good_regions) {
            ASSERT_TRUE(testRegion(region, min_identity));
            ASSERT_GTE(ssLength(region), min_length);
            StartStop copy = region;
            fixEnds(copy, min_end);
            ASSERT_TRUE(copy == region);
        }
        // global identity is good
        StartStop all(0, maxPos());
        if (!testRegion(all, min_identity)) {
            return false;
        }
        // find bad blocks
        Coordinates bad;
        complement(bad, good_regions);
        // no bad block >= MIN_LENGTH
        BOOST_FOREACH (const StartStop& bad_region, bad) {
            if (ssLength(bad_region) >= min_length) {
                return false;
            }
        }
        // first and last blocks are good
        if (!bad.empty()) {
            if (bad.front().first == 0) {
                return false;
            }
            if (bad.back().second == maxPos()) {
                return false;
            }
        }
        return false;
    }

    // O(N^2)
    void findGoodSlicesBrute(
        Coordinates& output,
        const Coordinates& good_regions,
        int min_length,
        int min_identity,
        int first_index,
        int last_index
    ) const {
        const StartStop invalid(-1, -1);
        StartStop best = invalid;
        for (int i = first_index; i <= last_index; i++) {
            for (int j = i; j <= last_index; j++) {
                StartStop region(
                    good_regions[i].first,
                    good_regions[j].second
                );
                if (ssLength(region) >= min_length
                    && testRegion(region, min_identity)) {
                    if (best == invalid
                        || ssLength(region) > ssLength(best)) {
                        best = region;
                    }
                }
            }
        }
        if (best != invalid) {
            output.push_back(best);
        }
    }

    void findGoodSlices(
        Coordinates& output,
        const Coordinates& good_regions,
        int min_length,
        int min_identity
    ) const {
        int n = good_regions.size();
        if (n == 0) {
            return;
        }
        int first_index = 0;
        for (int i = 0; i < n; i++) {
            if (i > 0) {
                StartStop bad_before(
                    good_regions[i - 1].second + 1,
                    good_regions[i].first - 1
                );
                if (ssLength(bad_before) >= min_length) {
                    // long bad region
                    int last_index = i - 1;
                    ASSERT_LTE(first_index, last_index);
                    findGoodSlicesBrute(
                        output,
                        good_regions,
                        min_length,
                        min_identity,
                        first_index,
                        last_index
                    );
                    first_index = i;
                }
            }
        }
        int last_index = n - 1;
        ASSERT_LTE(first_index, last_index);
        findGoodSlicesBrute(
            output,
            good_regions,
            min_length,
            min_identity,
            first_index,
            last_index
        );
    }
};

Coordinates goodSlices(const Scores& score,
                       int frame_length, int end_length,
                       int min_identity, int min_length) {
    GoodSlicer slicer(score);
    Coordinates good_regions;
    slicer.buildGoodRegions(
        good_regions,
        frame_length,
        min_length,
        min_identity,
        end_length
    );
    Coordinates good_slices;
    slicer.findGoodSlices(
        good_slices,
        good_regions,
        min_length,
        min_identity
    );
    return good_slices;
}

}
