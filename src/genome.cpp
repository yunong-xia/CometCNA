/*! @file genome.cpp
    @brief Implementation of Genome class
*/

#include "genome.hpp"

#include <algorithm>
#include <limits>
#include <cmath>
#include <fstream>
#include <stdexcept>
#include <unordered_set>
#include <random>
#include <sstream>

namespace tumopp
{

    namespace
    {

        enum class CNAEvent
        {
            focal_gain,
            focal_loss,
            chromosome_gain,
            chromosome_loss,
            arm_gain,
            arm_loss,
        };

        // Shortest autosome p arm in hg38 cytoband table: chr21p is 12,000,000 bp.
        constexpr size_t MIN_FOCAL_CNA_SIZE = 12000000;

        std::string chromosome_from_arm(const char *arm)
        {
            std::string chromosome;
            for (const char *p = arm; *p >= '0' && *p <= '9'; ++p)
            {
                chromosome += *p;
            }
            return chromosome;
        }

    }

    // ArmInterval is defined in genome.hpp. Should only be used by genome.cpp
    const std::array<ArmInterval, 44> &Genome::arm_intervals()
    {
        static const std::array<ArmInterval, 44> intervals{{
            {"1p", 0, 123400000},
            {"1q", 123400000, 248956422},
            {"2p", 0, 93900000},
            {"2q", 93900000, 242193529},
            {"3p", 0, 90900000},
            {"3q", 90900000, 198295559},
            {"4p", 0, 50000000},
            {"4q", 50000000, 190214555},
            {"5p", 0, 48800000},
            {"5q", 48800000, 181538259},
            {"6p", 0, 59800000},
            {"6q", 59800000, 170805979},
            {"7p", 0, 60100000},
            {"7q", 60100000, 159345973},
            {"8p", 0, 45200000},
            {"8q", 45200000, 145138636},
            {"9p", 0, 43000000},
            {"9q", 43000000, 138394717},
            {"10p", 0, 39800000},
            {"10q", 39800000, 133797422},
            {"11p", 0, 53400000},
            {"11q", 53400000, 135086622},
            {"12p", 0, 35500000},
            {"12q", 35500000, 133275309},
            {"13p", 0, 17700000},
            {"13q", 17700000, 114364328},
            {"14p", 0, 17200000},
            {"14q", 17200000, 107043718},
            {"15p", 0, 19000000},
            {"15q", 19000000, 101991189},
            {"16p", 0, 36800000},
            {"16q", 36800000, 90338345},
            {"17p", 0, 25100000},
            {"17q", 25100000, 83257441},
            {"18p", 0, 18500000},
            {"18q", 18500000, 80373285},
            {"19p", 0, 26200000},
            {"19q", 26200000, 58617616},
            {"20p", 0, 28100000},
            {"20q", 28100000, 64444167},
            {"21p", 0, 12000000},
            {"21q", 12000000, 46709983},
            {"22p", 0, 15000000},
            {"22q", 15000000, 50818468},
        }};
        return intervals;
    }

    size_t Genome::arm_length(const ArmInterval &interval)
    {
        return interval.end - interval.start;
    }

    size_t Genome::chromosome_length(const size_t chr)
    {
        return arm_length(arm_intervals()[2 * (chr - 1)]) +
               arm_length(arm_intervals()[2 * (chr - 1) + 1]);
    }

    int Genome::bounded_cn(const int cn, const int delta)
    {
        if (cn == 0)
        {
            return 0;
        }
        return std::max(0, cn + delta);
    }

    // Omitted arms retain weight 1. Parse once, then share immutable weights across genomes.
    std::shared_ptr<const Genome::ArmWeights> Genome::load_arm_weights(const std::string &path)
    {

        // initialize arm CIN weights
        auto weights = std::make_shared<ArmWeights>();
        weights->fill(1.0);

        // if the path of the CIN arm weight configuration file is ""
        // each arm has equal weight = 1 (return the initialized weights)
        if (path.empty())
        {
            return weights;
        }

        // otherwise, try to read the configurtion file
        std::ifstream input(path);
        // if cannot open this file, throw error
        if (!input)
            throw std::runtime_error("Cannot open CNA arm weights: " + path);

        std::string line;
        size_t line_number = 0;
        bool header_seen = false;
        std::unordered_set<std::string> seen;
        while (std::getline(input, line))
        {
            ++line_number;
            std::istringstream row(line);
            std::string arm, value, extra;
            if (!(row >> arm))
                continue;
            const auto fail = [&]()
            {
                throw std::runtime_error("Invalid CNA arm weights at " + path + ":" +
                                         std::to_string(line_number) + ": " + line);
            };
            if (!(row >> value) || (row >> extra))
                fail();
            if (!header_seen)
            {
                if (arm != "arm" || value != "CNA_weight")
                    fail();
                header_seen = true;
                continue;
            }
            const auto &intervals = arm_intervals();
            const auto found = std::find_if(intervals.begin(), intervals.end(),
                                            [&](const ArmInterval &interval)
                                            { return arm == interval.arm; });
            if (found == intervals.end() || !seen.insert(arm).second)
                fail();
            double weight = 0.0;
            try
            {
                size_t consumed = 0;
                weight = std::stod(value, &consumed);
                if (consumed != value.size())
                    fail();
            }
            catch (const std::exception &)
            {
                fail();
            }
            if (!std::isfinite(weight) || weight < 0.0)
                fail();
            (*weights)[static_cast<size_t>(found - intervals.begin())] = weight;
        }
        if (input.bad())
            throw std::runtime_error("Cannot read CNA arm weights: " + path);
        if (!header_seen)
            throw std::runtime_error("Missing arm CNA_weight header: " + path);
        double total = 0.0;
        for (const double weight : *weights)
            total += weight;
        if (!std::isfinite(total) || total <= 0.0)
        {
            throw std::runtime_error("CNA arm weights must have a finite positive sum: " + path);
        }
        return weights;
    }

    std::unordered_map<std::string, double> Genome::named_arm_weights(const ArmWeights &weights)
    {
        std::unordered_map<std::string, double> result;
        for (size_t i = 0; i < weights.size(); ++i)
            result[arm_intervals()[i].arm] = weights[i];
        return result;
    }

    // load arm selection coefficients from a configuration file.
    // If the path is empty, return a shared pointer to an array of 0.0s.
    std::shared_ptr<const Genome::ArmSelectionCoefs> Genome::load_arm_selection_coefs(const std::string &path)
    {
        // initialize arm CIN weights
        auto selection_coefs = std::make_shared<ArmSelectionCoefs>();
        selection_coefs->fill(0.0);

        // if the string of path variable is ""
        // return neutral fitness to all arms
        if (path.empty())
        {
            return selection_coefs;
        }

        // otherwise, try to read the configurtion file
        std::ifstream input(path);
        // if cannot open this file, throw error
        if (!input)
            throw std::runtime_error("Cannot open CNA arm selection coefficients: " + path);

        // now read configuration file

        std::string line;
        size_t line_number = 0;
        bool header_seen = false;
        std::unordered_set<std::string> seen;
        while (std::getline(input, line))
        {                                  // read a line from the configuration file
            ++line_number;                 // increment line number
            std::istringstream row(line);  // create a string stream to parse the line
            std::string arm, value, extra; // parse the line into three parts: arm, value, and extra
            if (!(row >> arm))
                continue;             // if the line is empty, continue to the next line
            const auto fail = [&]() { // define a lambda function to throw an error if the line is invalid
                throw std::runtime_error("Invalid CNA arm selection coefficient at " + path + ":" +
                                         std::to_string(line_number) + ": " + line);
            };
            if (!(row >> value) || (row >> extra))
                fail(); // if the line does not have exactly two parts, throw an error
            if (!header_seen)
            { // if the header has not been seen yet, check if the first line is the header
                if (arm != "arm" || value != "selection_coef")
                    fail();
                header_seen = true;
                continue;
            }
            const auto &intervals = arm_intervals();                            // get the arm intervals from the genome
            const auto found = std::find_if(intervals.begin(), intervals.end(), // find the arm in the intervals
                                            [&](const ArmInterval &interval)
                                            { return arm == interval.arm; });
            if (found == intervals.end() || !seen.insert(arm).second)
                fail();                  // if the arm is not found or has already been seen, throw an error
            double selection_coef = 0.0; // initialize the weight to 0.0
            try
            {
                size_t consumed = 0;
                selection_coef = std::stod(value, &consumed);
                if (consumed != value.size())
                    fail();
            }
            catch (const std::exception &)
            {
                fail();
            }
            if (!std::isfinite(selection_coef) || selection_coef < 0.0)
                fail();
            (*selection_coefs)[static_cast<size_t>(found - intervals.begin())] = selection_coef; // set the selection coefficient for the arm
        }
        if (input.bad())
            throw std::runtime_error("Cannot read CNA arm selection coefficients: " + path);
        if (!header_seen)
            throw std::runtime_error("Missing arm CNA_selection_coef header: " + path);

        return selection_coefs;
    }
    
    // Convert the ArmSelectionCoefs array to a map of arm names to selection coefficients.
    std::unordered_map<std::string, double> Genome::named_arm_selection_coefs(const ArmSelectionCoefs &selection_coefs)
    {
        std::unordered_map<std::string, double> result;
        for (size_t i = 0; i < selection_coefs.size(); ++i)
        {
            result[arm_intervals()[i].arm] = selection_coefs[i];
        }
        return result;
    }


    double Genome::get_arm_selection_coef(const std::string &arm)
    {
        const auto &intervals = arm_intervals();

        const auto found = std::find_if(
            intervals.begin(), intervals.end(),
            [&](const ArmInterval &interval)
            {
                return arm == interval.arm;
            });

        if (found == intervals.end())
        {
            throw std::invalid_argument("Unknown chromosome arm: " + arm);
        }

        const size_t index =
            static_cast<size_t>(found - intervals.begin());

        return selection_coefs_->at(index);
    }

    // Bias arm selection for focal and arm-level events only.
    // Whole-chromosome events select uniformly among the 22 autosomes.
    const ArmInterval &Genome::select_arm_biased(urbg_t &engine4)
    {
        std::discrete_distribution<size_t> selector(arm_weights_->begin(), arm_weights_->end());
        return arm_intervals()[selector(engine4)];
    }

    // Each cell divions has a probability to generate a new CNA in one of the daughter cells. The new CNA can be a focal gain, focal loss, chromosome gain, chromosome loss, arm gain, or arm loss. The probabilities of each event are equal (1/6). This is based on the model from Minussi, Navins. et al. 2017. "Breast Tumor Heterogeneity: Source of Fitness, Hurdle for Therapy." Trends in Cancer 3 (5): 474–93. https://doi.org/10.1016/j.trecan.2017.04.001.
    // This function is just to choose which CNA event to occur,
    // and implement the effects of the CNA event on the genome.
    std::string Genome::mutate_cna_minussi_navins(urbg_t &engine4)
    {
        std::uniform_int_distribution<int> event_dist(0, 5);
        const int event = event_dist(engine4);
        std::ostringstream oss;

        // Focal gain: duplicate a short interval within one chromosome arm.
        if (event == static_cast<int>(CNAEvent::focal_gain))
        {
            const auto &arm_interval = select_arm_biased(engine4);
            const size_t length = arm_length(arm_interval);

            // first select the focal gain size, then place it within the arm
            std::uniform_int_distribution<size_t> size_selector(MIN_FOCAL_CNA_SIZE, length);
            const size_t event_size = size_selector(engine4);
            std::uniform_int_distribution<size_t> start_selector(
                arm_interval.start,
                arm_interval.end - event_size);
            const size_t start = start_selector(engine4);
            const size_t end = start + event_size;

            // apply the focal gain to the genome
            apply_focal_delta(arm_interval.arm, start, end, +1);

            oss << "focal_gain\t" << chromosome_from_arm(arm_interval.arm)
                << "\t" << arm_interval.arm << "\t" << start << "\t" << end << "\n";
            return oss.str();
            // Focal loss: delete a short interval within one chromosome arm.
        }
        else if (event == static_cast<int>(CNAEvent::focal_loss))
        {
            const auto &arm_interval = select_arm_biased(engine4);
            const size_t length = arm_length(arm_interval);

            // first select the focal loss size, then place it within the arm
            std::uniform_int_distribution<size_t> size_selector(MIN_FOCAL_CNA_SIZE, length);
            const size_t event_size = size_selector(engine4);
            std::uniform_int_distribution<size_t> start_selector(
                arm_interval.start,
                arm_interval.end - event_size);
            const size_t start = start_selector(engine4);
            const size_t end = start + event_size;

            // apply the focal loss to the genome
            apply_focal_delta(arm_interval.arm, start, end, -1);

            oss << "focal_loss\t" << chromosome_from_arm(arm_interval.arm)
                << "\t" << arm_interval.arm << "\t" << start << "\t" << end << "\n";
            return oss.str();
            // Chromosome gain: increase copy number across both p and q arms.
        }
        else if (event == static_cast<int>(CNAEvent::chromosome_gain))
        {
            // randomly select one chromosome
            std::uniform_int_distribution<size_t> chr_selector(1, 22);
            const size_t chr = chr_selector(engine4);

            // chr gain
            apply_chromosome_delta(chr, +1);
            oss << "chromosome_gain\t" << chr << "\t\t0\t" << chromosome_length(chr) << "\n";
            return oss.str();
            // Chromosome loss: decrease copy number across both p and q arms.
        }
        else if (event == static_cast<int>(CNAEvent::chromosome_loss))
        {
            // randomly select one chromosome
            std::uniform_int_distribution<size_t> chr_selector(1, 22);
            const size_t chr = chr_selector(engine4);

            // chr loss
            apply_chromosome_delta(chr, -1);
            oss << "chromosome_loss\t" << chr << "\t\t0\t" << chromosome_length(chr) << "\n";
            return oss.str();
            // Arm gain: increase copy number across one selected chromosome arm.
        }
        else if (event == static_cast<int>(CNAEvent::arm_gain))
        {
            const auto &arm_interval = select_arm_biased(engine4);

            // arm gain
            apply_arm_delta(arm_interval.arm, +1);
            oss << "arm_gain\t" << chromosome_from_arm(arm_interval.arm)
                << "\t" << arm_interval.arm << "\t" << arm_interval.start
                << "\t" << arm_interval.end << "\n";
            return oss.str();
            // Arm loss: decrease copy number across one selected chromosome arm.
        }
        else if (event == static_cast<int>(CNAEvent::arm_loss))
        {
            const auto &arm_interval = select_arm_biased(engine4);

            // arm loss
            apply_arm_delta(arm_interval.arm, -1);
            oss << "arm_loss\t" << chromosome_from_arm(arm_interval.arm)
                << "\t" << arm_interval.arm << "\t" << arm_interval.start
                << "\t" << arm_interval.end << "\n";
            return oss.str();
        }

        return "";
    }

    // this function is called when the cell wgd function is called.
    std::string Genome::mutate_wgd()
    {
        for (auto &pair : arm_cns)
        {
            auto &cns = pair.second;
            for (int &cn : cns)
            {
                cn *= 2;
            }
        }
        return "wgd\t\t\t\t";
    }

    bool Genome::is_viable(const double max_ploidy,
                           const int max_segment_cn,
                           const double max_normalized_segment_cn,
                           const double max_nullisomy_fraction) const
    {

        double ploidy = 0.0;
        double total_autosome_length = 0.0;
        int highest_segment_cn = 0;
        double highest_nullisomy_fraction = 0.0;

        for (const auto &arm_interval : arm_intervals())
        {
            const double current_arm_length = static_cast<double>(arm_length(arm_interval));
            total_autosome_length += current_arm_length;
            double current_arm_nullisomy_length = 0.0;

            const auto &bps = arm_breakpoints.at(arm_interval.arm);
            const auto &cns = arm_cns.at(arm_interval.arm);

            for (size_t i = 0; i < cns.size(); ++i)
            {
                highest_segment_cn = std::max(highest_segment_cn, cns[i]);

                const double segment_length = static_cast<double>(bps[i + 1] - bps[i]);

                // average ploidy:
                // sum total CN dosage across the autosomes, weighted by segment length.
                ploidy += segment_length * cns[i];

                // Nullisomy fraction of this arm.
                // count how much of this arm has total CN equal to 0.
                if (cns[i] == 0)
                {
                    current_arm_nullisomy_length += segment_length;
                }
            }

            highest_nullisomy_fraction = std::max(
                highest_nullisomy_fraction,
                current_arm_nullisomy_length / current_arm_length);
        }

        ploidy /= total_autosome_length;

        // CINner-style average ploidy checkpoint.
        if (ploidy <= 0.0 || ploidy > max_ploidy)
        {
            return false;
        }

        // CINner-style highest absolute CN checkpoint.
        if (highest_segment_cn > max_segment_cn)
        {
            return false;
        }

        // CINner-style nullisomy burden checkpoint, adapted from bins to arms.
        if (highest_nullisomy_fraction > max_nullisomy_fraction)
        {
            return false;
        }

        // CINner-style normalized-CN checkpoint:
        // If any segment has CN/ploidy > max_normalized_segment_cn, then inviable.
        for (const auto &arm_interval : arm_intervals())
        {
            const auto &cns = arm_cns.at(arm_interval.arm);
            for (const int cn : cns)
            {
                if (static_cast<double>(cn) / ploidy > max_normalized_segment_cn)
                {
                    return false;
                }
            }
        }

        return true;
    }

    void Genome::split_segment(const std::string &arm, size_t bp)
    {
        auto &bps = arm_breakpoints.at(arm);
        auto &cns = arm_cns.at(arm);

        auto pos = std::lower_bound(bps.begin(), bps.end(), bp);

        // Breakpoint already exists.
        if (pos != bps.end() && *pos == bp)
        {
            return;
        }

        size_t index = static_cast<size_t>(pos - bps.begin());

        // Do not insert outside the represented arm.
        if (index == 0 || index >= bps.size())
        {
            return;
        }

        bps.insert(pos, bp);

        // Duplicate the CN state of the segment that was split.
        cns.insert(cns.begin() + index, cns[index - 1]);
    }

    void Genome::apply_focal_delta(
        const std::string &arm,
        size_t start,
        size_t end,
        int delta)
    {
        split_segment(arm, start);
        split_segment(arm, end);

        auto &bps = arm_breakpoints.at(arm);
        auto &cns = arm_cns.at(arm);

        for (size_t i = 0; i < cns.size(); ++i)
        {
            if (bps[i] >= start && bps[i + 1] <= end)
            {
                cns[i] = bounded_cn(cns[i], delta);
            }
        }
    }

    void Genome::apply_arm_delta(const std::string &arm, const int delta)
    {
        auto &cns = arm_cns.at(arm);
        for (auto &cn : cns)
        {
            cn = bounded_cn(cn, delta);
        }
    }

    void Genome::apply_chromosome_delta(const size_t chr, const int delta)
    {
        apply_arm_delta(std::to_string(chr) + "p", delta);
        apply_arm_delta(std::to_string(chr) + "q", delta);
    }

}
