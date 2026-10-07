#include "genome.hpp"
#include "cell.hpp"
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <iostream>

namespace
{

    void check(bool condition,
               const std::string &error_message)
    {
        if (!condition)
            throw std::runtime_error(error_message);
    }

    void check_selection_coefs_initialization()
    {
        auto defaults_selection_coefs = tumopp::Genome::load_arm_selection_coefs();
        for (double selection_coef : *defaults_selection_coefs)
            check(selection_coef == 0.0,
                  "Default arm selection coefficients should be initialized to 0.0");
    }

    void check_selection_coefs_loading()
    {

        // write to a temporary file and load it back to check if the loading works correctly
        const std::string path = "test_arm_selection_coefs.tsv";

        const auto write = [&](const std::string &text)
        {
            std::ofstream out(path);
            check(out.is_open(), "Cannot create test configuration");
            out << text;
            out.close();
            check(!out.fail(), "Cannot write test configuration");
        };

        // 14q 0.5
        // 1p 0.0
        // 1q unspecified, should default to 0.0
        write("arm\tselection_coef\n14q\t0.5\n1p\t0.0\n");

        // load the selection coefficients from the temporary file
        auto selection_coefs = tumopp::Genome::load_arm_selection_coefs(path);

        // remove the temporary file after loading
        std::remove(path.c_str());

        // check if the loaded selection coefficients match the expected values
        auto named = tumopp::Genome::named_arm_selection_coefs(*selection_coefs);
        check(named.at("14q") == 0.5, "14q coefficient should be 0.5");
        check(named.at("1p") == 0.0, "Explicit 1p coefficient should be zero");
        check(named.at("1q") == 0.0, "Unspecified 1q coefficient should default to zero");
    }

    void check_selection_coefs_loading_invalid()
    {

        // create an invalid filepath
        const std::string path = "XXXXXX";

        // try to load from this invalid path
        bool rejected = false;
        try
        {
            tumopp::Genome::load_arm_selection_coefs(path);
        }
        catch (const std::runtime_error &)
        {
            rejected = true;
        }

        // check if the loading was rejected as expected
        check(rejected, "Invalid arm selection coefficients should be rejected");
    }

    void check_birth_rate_imposed_by_14q_gain()
    {
        
    }

}

int main()
{
    try
    {
        check_selection_coefs_initialization();
        check_selection_coefs_loading();
        check_selection_coefs_loading_invalid();
    }
    catch (const std::runtime_error &e)
    {
        std::cerr << e.what() << std::endl;
        return 1; // failure
    }

    std::cout << "All tests passed\n";
    return 0;
}
