
/*
Trio.cpp - This file is part of BayesTyper (v0.9)


The MIT License (MIT)

Copyright (c) 2016 Jonas Andreas Sibbesen and Lasse Maretty

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
*/


#include <algorithm>

#include "Trio.hpp"
#include "Auxiliaries.hpp"


Trio::Trio(Variant & cur_var, const TrioInfo & trio_info) {

    father = cur_var.getSample(trio_info.father);
    mother = cur_var.getSample(trio_info.mother);
    child = cur_var.getSample(trio_info.child);

    is_filtered = ((father.callStatus() != Sample::CallStatus::Complete) or (mother.callStatus() != Sample::CallStatus::Complete) or (child.callStatus() != Sample::CallStatus::Complete));
    is_diploid = ((father.ploidy() == Sample::Ploidy::Diploid) and (mother.ploidy() == Sample::Ploidy::Diploid) and (child.ploidy() == Sample::Ploidy::Diploid));

    assert((father.ploidy() != Sample::Ploidy::Zeroploid) or (mother.ploidy() != Sample::Ploidy::Zeroploid) or (child.ploidy() != Sample::Ploidy::Zeroploid));

    is_informative = true;
    is_concordant = false;
    is_exclusively_child_heterozygote = true;
    is_parents_bi_allelelic_heterozygote = false;
    is_reference_call = true;
    has_called_missing = false;

    if (!is_filtered) {

        auto father_genotype = father.genotypeEstimate();
        auto mother_genotype = mother.genotypeEstimate();
        auto child_genotype = child.genotypeEstimate();

        bool first_child_allele_in_father = false;
        bool second_child_allele_in_father = false;
        bool first_child_allele_in_mother = false;
        bool second_child_allele_in_mother = false;

        if (child_genotype.size() > 0) {

            assert(child.ploidy() != Sample::Ploidy::Zeroploid);

            first_child_allele_in_father = (find(father_genotype.begin(), father_genotype.end(), child_genotype.front()) != father_genotype.end());
            second_child_allele_in_father = (find(father_genotype.begin(), father_genotype.end(), child_genotype.back()) != father_genotype.end());
            first_child_allele_in_mother = (find(mother_genotype.begin(), mother_genotype.end(), child_genotype.front()) != mother_genotype.end());
            second_child_allele_in_mother = (find(mother_genotype.begin(), mother_genotype.end(), child_genotype.back()) != mother_genotype.end());
        }

        if (cur_var.chrom() == "chrX") {

            if (!is_diploid) {

                assert(father.ploidy() == Sample::Ploidy::Haploid);
                assert(father_genotype.size() == 1);

                assert(mother.ploidy() == Sample::Ploidy::Diploid);
                assert(mother_genotype.size() == 2);

                if (child.ploidy() == Sample::Ploidy::Diploid) {

                    assert(child_genotype.size() == 2);

                    is_concordant = ((first_child_allele_in_father and second_child_allele_in_mother) or (second_child_allele_in_father and first_child_allele_in_mother));

                } else {

                    assert(child.ploidy() == Sample::Ploidy::Haploid);
                    assert(child_genotype.size() == 1);

                    assert(first_child_allele_in_mother == second_child_allele_in_mother);

                    is_concordant = first_child_allele_in_mother;
                }

            } else {

                is_informative = false;
                is_concordant = false;
            }

        } else if (cur_var.chrom() == "chrY") {

            if (!is_diploid) {

                assert(father.ploidy() == Sample::Ploidy::Haploid);
                assert(father_genotype.size() == 1);

                assert(mother.ploidy() == Sample::Ploidy::Zeroploid);
                assert(mother_genotype.size() == 0);

                assert(!first_child_allele_in_mother and !second_child_allele_in_mother);
                assert(first_child_allele_in_father == second_child_allele_in_father);

                if (child.ploidy() == Sample::Ploidy::Haploid) {

                    assert(child_genotype.size() == 1);

                    is_concordant = first_child_allele_in_father;

                } else {

                    assert(child.ploidy() == Sample::Ploidy::Zeroploid);
                    assert(child_genotype.size() == 0);

                    assert(!first_child_allele_in_father);

                    is_informative = false;
                    is_concordant = false;
                }

            } else {

                is_informative = false;
                is_concordant = false;
            }

        } else {

            assert(is_diploid);

            assert(father_genotype.size() == 2);
            assert(mother_genotype.size() == 2);
            assert(child_genotype.size() == 2);

            is_concordant = ((first_child_allele_in_father and second_child_allele_in_mother) or (second_child_allele_in_father and first_child_allele_in_mother));

            if (cur_var.numAlls() == 2) {

                if ((father_genotype.front() != father_genotype.back()) and (mother_genotype.front() != mother_genotype.back())) {

                    is_parents_bi_allelelic_heterozygote = true;
                }
            }
        }

        if (child.ploidy() != Sample::Ploidy::Diploid) {

            is_exclusively_child_heterozygote = false;

        } else if ((father.ploidy() == Sample::Ploidy::Zeroploid) or (mother.ploidy() == Sample::Ploidy::Zeroploid)) {

            is_exclusively_child_heterozygote = false;

        } else if ((father_genotype.front() != father_genotype.back()) or (mother_genotype.front() != mother_genotype.back())) {

            is_exclusively_child_heterozygote = false;

        } else if (child_genotype.front() == child_genotype.back()) {

            is_exclusively_child_heterozygote = false;
        }

        auto all_called_alleles = father_genotype;
        all_called_alleles.insert(all_called_alleles.end(), mother_genotype.begin(), mother_genotype.end());
        all_called_alleles.insert(all_called_alleles.end(), child_genotype.begin(), child_genotype.end());

        assert(all_called_alleles.size() > 0);

        for (auto &allele_idx: all_called_alleles) {

            if (cur_var.allele(allele_idx).isMissing()) {

                has_called_missing = true;
            }

            if (allele_idx > 0) {

                is_reference_call = false;
            }
        }

        de_novo_event.first = false;

        if (is_diploid and !is_concordant and is_exclusively_child_heterozygote and !has_called_missing) {

            uint ancestral_allele_idx = 0;
            uint de_novo_allele_idx = 0;

            if ((!first_child_allele_in_father and !first_child_allele_in_mother) and (second_child_allele_in_father or second_child_allele_in_mother)) {

                de_novo_allele_idx = child_genotype.front();

            } else if ((!second_child_allele_in_father and !second_child_allele_in_mother) and (first_child_allele_in_father or first_child_allele_in_mother)) {

                de_novo_allele_idx = child_genotype.back();
            }

            if (de_novo_allele_idx == father_genotype.front()) {

                ancestral_allele_idx = mother_genotype.front();

            } else {

                ancestral_allele_idx = father_genotype.front();
            }

            if (de_novo_allele_idx > 0) {

                assert(de_novo_allele_idx != ancestral_allele_idx);

                de_novo_event.first = true;

                de_novo_event.second.child_id = trio_info.child;
                de_novo_event.second.ancestral_allele_idx = ancestral_allele_idx;
                de_novo_event.second.de_novo_allele_idx = de_novo_allele_idx;
            }
        }
    }
}

bool Trio::isFiltered() {

    return is_filtered;
}

bool Trio::isDiploid() {

    return is_diploid;
}

bool Trio::isInformative() {

    assert(!is_filtered);

    return is_informative;
}

bool Trio::isConcordant() {

    assert(!is_filtered);

    return is_concordant;
}

bool Trio::isExclusivelyChildHeterozygote() {

    assert(!is_filtered);

    return is_exclusively_child_heterozygote;
}

bool Trio::isParentsBiAllelicHeterozygote() {

    assert(!is_filtered);

    return is_parents_bi_allelelic_heterozygote;
}

bool Trio::isReferenceCall() {

    assert(!is_filtered);

    return is_reference_call;
}

bool Trio::hasCalledMissing() {

    assert(!is_filtered);

    return has_called_missing;
}

bool Trio::isDeNovo() {

    assert(!is_filtered);
    assert(is_diploid);

    return de_novo_event.first;
}

Trio::DeNovoEvent Trio::deNovoEvent() {

    assert(!is_filtered);
    assert(is_diploid);

    return de_novo_event.second;
}

vector<Trio::TrioInfo> Trio::parsePedigree(const VcfMetaData & meta_data, const string & trio_info_str) {

    vector<TrioInfo> trios;

    assert(!(trio_info_str.empty()));
    auto trio_info_str_split = Utils::splitString(trio_info_str, ':');

    uint suffix_length = 2;

    trios.reserve(trio_info_str_split.size());
    assert(trio_info_str_split.size() < 100);

    for (auto &trio_str: trio_info_str_split) {

        auto trio_str_split = Utils::splitString(trio_str, ',');
        assert(trio_str_split.size() == 3);

        TrioInfo trio_info;

        trio_info.id = to_string(trios.size() + 1);
        assert(trio_info.id.size() <= suffix_length);

        while (trio_info.id.size() < suffix_length) {

            trio_info.id.insert(0, "0");
        }

        trio_info.id.insert(0, "Trio");

        trio_info.father = trio_str_split.at(0);
        trio_info.mother = trio_str_split.at(1);
        trio_info.child = trio_str_split.at(2);

        meta_data.hasSampleId(trio_info.father);
        meta_data.hasSampleId(trio_info.mother);
        meta_data.hasSampleId(trio_info.child);

        trios.push_back(trio_info);
    }

    return trios;
}

vector<Trio::TrioInfo> Trio::parseGenomeDKPedigree(const VcfMetaData & meta_data) {

    vector<TrioInfo> trios;

    assert((meta_data.sampleIds().size() % 3) == 0);

    map<string, TrioInfo> temp_trios;

    for (auto &sample_id: meta_data.sampleIds()) {

        auto sample_id_split = Utils::splitString(sample_id, '-');
        assert(sample_id_split.size() == 2);

        auto family_member_id = stoi(sample_id_split.back());
        assert(family_member_id > 0);

        auto temp_trios_emplace = temp_trios.emplace(sample_id_split.front(), TrioInfo());
        temp_trios_emplace.first->second.id = sample_id_split.front();

        if (family_member_id == 1) {

            assert(temp_trios_emplace.first->second.father.empty());
            temp_trios_emplace.first->second.father = sample_id;

        } else if (family_member_id == 2) {

            assert(temp_trios_emplace.first->second.mother.empty());
            temp_trios_emplace.first->second.mother = sample_id;

        } else {

            assert(family_member_id <= 99);

            assert(temp_trios_emplace.first->second.child.empty());
            temp_trios_emplace.first->second.child = sample_id;
        }
    }

    assert((meta_data.sampleIds().size() / 3) == temp_trios.size());
    trios.reserve(temp_trios.size());

    for (auto &trio: temp_trios) {

        meta_data.hasSampleId(trio.second.father);
        meta_data.hasSampleId(trio.second.mother);
        meta_data.hasSampleId(trio.second.child);

        assert(trio.first == trio.second.id);
        assert(trio.first == Utils::splitString(trio.second.father, '-').front());
        assert(trio.first == Utils::splitString(trio.second.mother, '-').front());
        assert(trio.first == Utils::splitString(trio.second.child, '-').front());

        trios.push_back(trio.second);
    }

    return trios;
}
