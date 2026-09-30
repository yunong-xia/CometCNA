#include "genome.hpp"
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>

void check(bool condition) {
    if (!condition) throw std::runtime_error("arm weight regression failed");
}

int main() {
    const std::string path = "test_arm_weights.tsv";
    const auto write = [&](const std::string& text) { std::ofstream out(path); out << text; };
    auto defaults = tumopp::Genome::load_arm_weights();
    for (double weight : *defaults) check(weight == 1.0);
    write("arm\tCNA_weight\n8q\t10\n1p\t0\n");
    auto weights = tumopp::Genome::load_arm_weights(path);
    auto named = tumopp::Genome::named_arm_weights(*weights);
    check(named.at("8q") == 10 && named.at("1p") == 0 && named.at("1q") == 1);
    for (const std::string bad : {
        "arm weight\n8q 10\n", "arm CNA_weight\n8q -1\n",
        "arm CNA_weight\n8q nan\n", "arm CNA_weight\n8q inf\n",
        "arm CNA_weight\n8q 2oops\n", "arm CNA_weight\n8q 2 extra\n",
        "arm CNA_weight\n8q 2\n8q 3\n", "arm CNA_weight\n23q 2\n", ""}) {
        write(bad);
        bool rejected = false;
        try { tumopp::Genome::load_arm_weights(path); }
        catch (const std::runtime_error&) { rejected = true; }
        check(rejected);
    }
    std::string table = "arm CNA_weight\n";
    for (int chr = 1; chr <= 22; ++chr)
        for (char arm : {'p', 'q'}) table += std::to_string(chr) + arm + " 0\n";
    write(table);
    bool rejected = false;
    try { tumopp::Genome::load_arm_weights(path); }
    catch (const std::runtime_error&) { rejected = true; }
    check(rejected);
    table.replace(table.find("8q 0"), 4, "8q 1");
    write(table);
    weights = tumopp::Genome::load_arm_weights(path);
    tumopp::urbg_t engine(42);
    int arm_events = 0, chromosome_events = 0, chromosome8 = 0;
    for (int i = 0; i < 6000; ++i) {
        tumopp::Genome genome(weights);
        std::istringstream record(genome.mutate_cna_minussi_navins(engine));
        std::string event, chr, arm;
        std::getline(record, event, '\t');
        std::getline(record, chr, '\t');
        std::getline(record, arm, '\t');
        if (event.find("chromosome_") == 0) {
            ++chromosome_events;
            if (chr == "8") ++chromosome8;
        } else {
            ++arm_events;
            check(arm == "8q");
        }
    }
    check(arm_events > 3000 && chromosome_events > 1500);
    check(chromosome8 > 30 && chromosome8 < 160);
    // Loading another configuration does not mutate existing shared weights.
    check((*weights)[15] == 1);
    for (double weight : *defaults) check(weight == 1);
}
