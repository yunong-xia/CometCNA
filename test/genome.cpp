#include "genome.hpp"

#include <cassert>
#include <fstream>
#include <string>
#include <vector>

int main() {
    tumopp::Genome genome;
    assert(genome.is_viable());
    assert(!genome.is_viable(1.0, 8, 4.0, 0.2));
    assert(!genome.is_viable(8.0, 1, 4.0, 0.2));
    assert(!genome.is_viable(8.0, 8, 0.5, 0.2));
    assert(genome.is_viable(8.0, 8, 4.0, 0.0));

    for (int chr = 1; chr <= 22; ++chr) {
        for (const char arm_suffix : {'p', 'q'}) {
            const std::string arm = std::to_string(chr) + arm_suffix;
            assert(genome.breakpoints(arm).size() == 2);
            assert(genome.cn_states(arm).size() == 1);
            assert(genome.cn_states(arm).front() == 2);
        }
    }

    assert(genome.breakpoints("1p").front() == 0);
    assert(genome.breakpoints("1p").back() == 123400000);
    assert(genome.breakpoints("1q").front() == 123400000);
    assert(genome.breakpoints("1q").back() == 248956422);
    assert(genome.breakpoints("21p").front() == 0);
    assert(genome.breakpoints("21p").back() == 12000000);
    assert(genome.breakpoints("21q").front() == 12000000);
    assert(genome.breakpoints("21q").back() == 46709983);

    tumopp::Genome recorded_genome;
    tumopp::urbg_t record_engine(1);
    std::vector<std::string> cna_events;
    cna_events.reserve(20);
    cna_events.push_back(recorded_genome.mutate_wgd());
    for (int i = 1; i < 20; ++i) {
        cna_events.push_back(recorded_genome.mutate_cna_minussi_navins(record_engine));
    }

    std::ofstream cna_record("genome_cna_events.tsv");
    cna_record << "cna_event\tchr\tarm\tstart\tend\n";
    for (const auto& event : cna_events) {
        assert(!event.empty());
        assert(event.find("NA") == std::string::npos);
        cna_record << event;
        if (event.back() != '\n') {
            cna_record << '\n';
        }
    }
    cna_record.close();
    assert(cna_record);
    assert(cna_events.size() == 20);
    assert(cna_events.front() == "wgd\t\t\t\t");

    tumopp::urbg_t engine(1);
    for (int i = 0; i < 2000; ++i) {
        genome.mutate_cna_minussi_navins(engine);
    }

    for (int chr = 1; chr <= 22; ++chr) {
        for (const char arm_suffix : {'p', 'q'}) {
            const std::string arm = std::to_string(chr) + arm_suffix;
            const auto& breakpoints = genome.breakpoints(arm);
            const auto& cns = genome.cn_states(arm);

            assert(breakpoints.size() == cns.size() + 1);
            for (size_t j = 0; j + 1 < breakpoints.size(); ++j) {
                assert(breakpoints[j] < breakpoints[j + 1]);
                assert(cns[j] >= 0);
            }
        }
    }

    return 0;
}
