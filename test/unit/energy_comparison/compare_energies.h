#include <iostream>
#include <fstream>
#include <string>

#include <limits> // epsilon for limit
#include <utility> // pair

#include <bitset> // TODO: Remove

namespace test_utils {
/**
 * @brief Helper function to write collective errors to file for further analysis
 *
 * @param errs The vector of all errors
 * @param field_per_line The number of values to write per file line
 */
void write_error_ouput( std::vector<double> errs, int field_per_line, std::string err_file_base_name)
{
    int counter = 0;
    std::ofstream outputFile(err_file_base_name);

    for (auto e : errs)
    {
        counter++;
        outputFile << counter << " " << e*100.0 << " "; // Convert to percent and dump
        if (counter % field_per_line == 0)
        {
            outputFile << std::endl;
        }
    }
    outputFile.close();
}

/**
 * @brief Helper function to compare numbers and calculate a absolute error
 *
 * @param A The first value to compare
 * @param B The second value to compare
 *
 * @return The calculated error
 */
double calculate_abs_error(double A, double B)
{
    return std::abs(A-B) / std::min(A,B);
}

/**
 * @brief Helper function to compare numbers and calculate a relative error
 *
 * @param A The first value to compare
 * @param B The second value to compare
 *
 * @return The calculated error
 */
double calculate_relative_error(double A, double B)
{
     return std::abs(A-B);
}

/**
 * @brief Function to compare errors to a given tolerance, and decide if it's within range
 *
 * @param A The first value to compare
 * @param B The second value to compare
 * @param relative_tolerance The relative tolerance to use when comparing
 *
 * @return A pair containing true/false if it's within tolerance, and the calculated error
 */
std::pair<bool, double> compare_error(double A, double B, double relative_tolerance)
{
    bool within_tol = false;
    double err = 0.0;

    // Right now this is pretty arbitrary..
    double abs_threshhold = 10 * std::numeric_limits<float>::epsilon();

    // Calculate if we're withing tolerances
    // If we're close to relative, do absolute
    if (std::abs(std::min(A,B)) < abs_threshhold)
    {
        err = calculate_relative_error(A, B);

        // Finding a relative error to 0 doesn't make much
        // sense, so lets do absolute error instead
        if ( err < std::numeric_limits<double>::epsilon() )
        {
            within_tol = true;
        }
        else {
            within_tol = false;
        }
    }
    else { // Do absolute error

        err = calculate_abs_error(A, B);

        if (err < relative_tolerance)
        {
            within_tol = true;
        }
        else {
            within_tol = false;
        }
    }
    return { within_tol, err };
}

enum FIELD_ENUM {
    Individual = 0, // Track each field individually
    Sum // Sum the masked fields
};

/**
 * @brief Function to compare the contents of two energy files
 *
 * @param file_a First file to compare
 * @param file_b Second file to compare
 * @param relative_tolerance Relative tolerance which is acceptable
 * @param field_mask A mask to specify which fields in the file to use
 * @param sum_mask A mask to specify which fields in the file to sum and compare
 * @param write_err_output If you should write the error output to a file
 * @param err_file_base_name Base filename for writing output
 * @param num_lines_to_skip The number of lines to skup into the file
 *
 * @NOTE A typical energy file is:
 * <step> <ex> <ey> <ez> <bx> <by> <bz> <particle energies...>
 * and the bit maps go accordingly with <step> being the LSB.
 * A mask for b fields only would be 0x000001110
 *
 * @NOTE We could * use bitsets for the masking but * they're generally slower
 *
 * @return True is they match (within tol), false if not
 */
bool compare_energies(
        const std::string file_a,
        const std::string file_b,
        const double relative_tolerance,
        const unsigned short field_mask = 0b1111111111111111, /// short has 16 bytes, assume all are true
        const FIELD_ENUM field_enum = FIELD_ENUM::Individual, /// short has 16 bytes, assume all are true
        const int write_err_ouput = 0, // If the run should dump the errors to disk
        const std::string err_file_base_name =  "err.out", // File name to write errors to
        const int num_lines_to_skip = 0 // Most energy files have 3 lines of padding
)
{
    // TODO: I could easily have a policy here based on the type of the field_mask
    std::vector<double> errs;

    const int DEFAULT_FILED_COUNT = 7;

    unsigned short agg_total = 0;
    unsigned short v = field_mask;
    // Count set bits
    for (agg_total = 0; v; agg_total++)
    {
        v &= v - 1; // clear the least significant bit set
    }

    try {

        bool match = true;

        std::string line1 = "";
        std::string line2 = "";

        std::ifstream f1 (file_a);
        std::ifstream f2 (file_b);

        double max_err = 0.0;
        int max_err_line = -1;

        // This is for counting the number of tokens on a line (changes
        // based on number of species). It can likely be done much better
        int line_token_count = 0;

        if (f1.is_open() && f2.is_open())
        {

            // Perform skipping
            for (int i = 0; i < num_lines_to_skip; i++)
            {
                getline(f1,line1);
                getline(f2,line2);
            }

            int counter = num_lines_to_skip;

            // Do processing
            while ( getline(f1,line1) )
            {
                getline(f2,line2);

                // Tokenize lines
                std::stringstream linestream1(line1);
                std::string item1;

                std::stringstream linestream2(line2);
                std::string item2;

                int used_line_token_count = 0;
                int total_line_token_count = 0;

                double sum_A = 0.0;
                double sum_B = 0.0;
                std::pair<bool, double> returned_err;
                returned_err.second = -1.0; // set a dummy value to show uninit

                int agg_count = 0;
                while (getline(linestream1, item1, ' '))
                {
                    bool write_this_err_ouput = write_err_ouput;
                    std::cout << "Setting write_this_err_ouput tp " << write_this_err_ouput << std::endl;

                    getline(linestream2, item2, ' ');
                    total_line_token_count++;

                    // Use this field
                    //std::cout << "this_line " << this_line_token_count << " mask " << field_mask << std::endl;

                    // Take the value one, and shift it to generate the mask to compare
                    unsigned short this_line_token_mask = 1 << (total_line_token_count - 1); // Set correct highest bit on
                    //this_line_token_mask |= this_line_token_mask-1; // Set lower bits on

                    // If this field is within our requested mask, use it
                    if (this_line_token_mask & field_mask)
                    {
                        used_line_token_count++;
                        std::cout << "Parsing field " << used_line_token_count << " val " << item1 << std::endl;

                        double A = std::stod(item1);
                        double B = std::stod(item2);

                        if (
                                (field_enum == FIELD_ENUM::Sum) && // Need to aggregate
                                (agg_count < agg_total) // Not done aggregating yet
                            )
                        {
                            // Need to aggregate..
                            sum_A += A;
                            sum_B += B;
                            agg_count++;

                            std::cout << "sum a " << sum_A << " += " << A << std::endl;
                            std::cout << "sum b " << sum_B << " += " << B << std::endl;

                            // Don't write this particular one
                            write_this_err_ouput = false;

                            if (agg_count == agg_total) { // final_aggregation
                                returned_err = compare_error(sum_A, sum_B, relative_tolerance);
                                write_this_err_ouput = true;
                            }
                        }
                        else // We can just compare this val
                        {

                            returned_err = compare_error(A, B, relative_tolerance);

                        }

                        if (returned_err.second != -1.0)  // Has some value set
                        {
                            bool returned_match = returned_err.first;

                            if (!returned_match) {
                                match = false;
                            }

                            double err = returned_err.second;

                            // Track max absolute error
                            if (err > max_err)
                            {
                                max_err = err;
                                max_err_line = counter;
                            }


                            // If we track the errors, track this one
                            if (write_this_err_ouput)
                            {
                                errs.push_back(err);
                            }
                        }
                    }
                    else {
                        std::cout << "Skipping field " << this_line_token_mask << " val " << item1 << std::endl;
                    }
                }
                line_token_count = used_line_token_count;
                counter++;
            }

            f1.close();
            f2.close();
        }
        else {
            std::cerr << "Unable to open file";
            return false;
        }

        std::cout << "Field mask : " << field_mask << std::endl;
        std::cout << "Fields used : " << line_token_count << std::endl;

        std::cout << "Max found err was " << max_err*100 << "% on line " << max_err_line << " (Threshold: " <<
            relative_tolerance*100 << "%)" << std::endl;

        if (write_err_ouput)
        {
            int err_per_line = line_token_count;
            if (field_enum == FIELD_ENUM::Sum) // Need to aggregate
            {
                err_per_line /= agg_total; // Reduce by aggregation factor
            }

            std::cout << "Writing error output " << errs.size() << std::endl;
            write_error_ouput( errs, err_per_line, err_file_base_name);
        }


        return match;
    }
    catch (const std::exception &exc) // Catching all is bad form, but OK for now..
    {
        // catch anything thrown within try block that derives from std::exception
        std::cerr << exc.what();
        return false;
    }

}

} // namespace
