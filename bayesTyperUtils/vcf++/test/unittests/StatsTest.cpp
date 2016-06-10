#include <map>

#include "catch.hpp"

#include "Variant.hpp"
#include "Allele.hpp"
#include "Sample.hpp"
#include "Attribute.hpp"

#include "Stats.hpp"

TEST_CASE("pop stats maf", "[Stats]") {

    string fmt_str("GT");
    map<string,Attribute::DetailedDescriptor> dscr_vec;

    Sample homozygot_ref_sample("0/0", fmt_str, dscr_vec, 2, false);
    Sample homozygot_alt_sample("1/1", fmt_str, dscr_vec, 2, false);
    Sample heterozygot_sample("0/1", fmt_str, dscr_vec, 2, false);

    REQUIRE(homozygot_ref_sample.ploidy() == Sample::Ploidy::Diploid);
    REQUIRE(homozygot_ref_sample.callStatus() == Sample::CallStatus::Complete);
    REQUIRE(homozygot_ref_sample.isInformative());

    REQUIRE(homozygot_alt_sample.ploidy() == Sample::Ploidy::Diploid);
    REQUIRE(homozygot_alt_sample.callStatus() == Sample::CallStatus::Complete);
    REQUIRE(homozygot_alt_sample.isInformative());

    REQUIRE(heterozygot_sample.ploidy() == Sample::Ploidy::Diploid);
    REQUIRE(heterozygot_sample.callStatus() == Sample::CallStatus::Complete);
    REQUIRE(heterozygot_sample.isInformative());

    Allele ref("A");
    vector<Allele> alts = {Allele("T")};
    AttributeSet info;

    Variant variant("chr1", 1, ref, alts, info);

    REQUIRE(variant.numAlls() == 2);

    for (uint i = 0; i < 25; i++) {

        REQUIRE(variant.addSample(to_string(i),homozygot_ref_sample));
    }

    for (uint i = 25; i < 50; i++) {

        REQUIRE(variant.addSample(to_string(i),homozygot_alt_sample));
    }

    for (uint i = 50; i < 100; i++) {

        REQUIRE(variant.addSample(to_string(i),heterozygot_sample));
    }

    // SECTION("allele_counter tests") {

    //     Stats::CompleteDiploidAlleleCounter allele_counter(2);
    //     auto num_parsed_samples = Stats::applyCounterToSamples(&variant, &allele_counter, regex(".+"));

    //     REQUIRE(allele_counter.alleleCountSum() == 200);
    //     REQUIRE(num_parsed_samples == 200);

    //     vector<float> correct_allele_freqs = {0.5, 0.5};
    //     REQUIRE(allele_counter.calcAlleleFreqs().first == correct_allele_freqs);
    //     REQUIRE(allele_counter.calcExpNumHeterozygotes().first == 50);
    // }

    // SECTION("inbreed_coef test") {

    //     auto inbreed_coef = Stats::calcInbreedingCoef(&variant);

    //     REQUIRE(inbreed_coef.second);
    //     REQUIRE(inbreed_coef.first == 0);
    // }
}
