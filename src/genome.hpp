/*! @file genome.hpp
    @brief Interface of Genome class
*/
#pragma once
#ifndef TUMOPP_GENOME_HPP_
#define TUMOPP_GENOME_HPP_

#include "random.hpp"

#include <algorithm>
#include <array>
#include <cstddef>
#include <string>
#include <unordered_map>
#include <vector>


namespace tumopp {

/*! arm size data
*/
struct ArmInterval {
    const char* arm;
    size_t start;
    size_t end;
};



class Genome{
  public:
    Genome() {
        initialize_arms();
    }
                
    const std::vector<size_t>& breakpoints(const std::string& arm) const {
        return arm_breakpoints.at(arm);
    }

    const std::vector<int>& cn_states(const std::string& arm) const {
        return arm_cns.at(arm);
    }

    void add_breakpoint(const std::string& arm, size_t bp);


    std::string mutate_cna_minussi_navins(urbg_t& engine4);
    std::string mutate_wgd();

    bool is_viable(double max_ploidy = 8.0,
                   int max_segment_cn = 8,
                   double max_normalized_segment_cn = 4.0,
                   double max_nullisomy_fraction = 0.2) const;

  private:
        
    // key structures
    // arm_breakpoints: BPs on each arm in chromosome coordinates.
    // For example, 1p starts at 0, but 1q starts at the chr1 centromere coordinate.
    // arm_cns: CN levels on each arm. (initially 2 for diploid)
    std::unordered_map<std::string, std::vector<size_t>> arm_breakpoints;
    std::unordered_map<std::string, std::vector<int>> arm_cns;
        
    // hardcoded arm intervals, given the hg38 cytoband data from UCSC 
    static const std::array<ArmInterval, 44>& arm_intervals();
        
    // BELOW ARE PRIVATE HELPER FUNCTIONS

    // hg 38 arm length 
    static size_t arm_length(const ArmInterval& interval);
    static size_t chromosome_length(size_t chr);

    static int bounded_cn(int cn, int delta);
        
    // at the beginning of the simulation,
    // initialize the arm_breakpoints and arm_cns for each arm
    void initialize_arms() {

        // reserve 44 arms for human genome (22 autosomes)
        arm_breakpoints.reserve(arm_intervals().size());
        arm_cns.reserve(arm_intervals().size());
            

        // initialize each arm with a single segment (CN=2) and breakpoints
        // in chromosome coordinates
        for (const auto& interval : arm_intervals()) {
            const std::string arm(interval.arm);
            arm_breakpoints.emplace(arm, std::vector<size_t>{interval.start, interval.end});
            arm_cns.emplace(arm, std::vector<int>{2});
        }
    }
        
    void split_segment(const std::string& arm, size_t bp);
    void apply_focal_delta(const std::string& arm, size_t start, size_t end, int delta);
    void apply_arm_delta(const std::string& arm, int delta);
    void apply_chromosome_delta(size_t chr, int delta);
};



}

#endif // TUMOPP_GENOME_HPP_
