#include <iostream>
#include <fstream>
#include <string>

#include <limits> // epsilon for limit
#include <utility> // pair

// TODO: add a namespace or similar

/**
 * @brief Helper function to write collective errors to file for further analysis
 *
 * @param errs The vector of all errors
 * @param field_per_line The number of values to write per file line
 */
void write_error_ouput( std::vector<double> errs, int field_per_line)
{
    int counter = 0;
    std::ofstream outputFile("err.out");

    for (auto e : errs)
    {
        outputFile << e*100.0 << " "; // Convert to percent and dump
        counter++;
        if (counter == field_per_line)
        {
            counter = 0;
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


/**
 * @brief Function to compare the contents of two energy files
 *
 * @param file_a First file to compare
 * @param file_b Second file to compare
 * @param relative_tolerance Relative tolerance which is acceptable
 * @param field_mask A mask to specify which fields in the file to use
 * @param sum_mask A mask to specify which fields in the file to sum and compare
 * @param num_lines_to_skip
 * @param write_err_ouput
 *
 * @return True is they match (within tol), false if not
 */
bool compare_energies(
        const std::string file_a,
        const std::string file_b,
        const double relative_tolerance,
        const unsigned short field_mask = 0b1111111111111111, /// short has 16 bytes, assume all are true
        const int num_lines_to_skip = 0, // Most energy files have 3 lines of padding
        const int write_err_ouput = 0 // If the run should dump the errors to disk
)
{
    std::vector<double> errs;

    const int DEFAULT_FILED_COUNT = 7;

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

                int this_line_token_count = 0;
                while (getline(linestream1, item1, ' '))
                {
                    this_line_token_count++;
                    getline(linestream2, item2, ' ');

                    double A = std::stod(item1);
                    double B = std::stod(item2);

                    std::pair<bool, double> returned_err = compare_error(A, B,
                            relative_tolerance);

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
                    if (write_err_ouput)
                    {
                        errs.push_back(err);
                    }
                }
                line_token_count = this_line_token_count;
                counter++;
            }

            f1.close();
            f2.close();
        }
        else {
            std::cerr << "Unable to open file";
            return false;
        }

        int num_species = line_token_count - DEFAULT_FILED_COUNT;
        std::cout << "Num species detected: " << num_species << std::endl;

        std::cout << "Max found err was " << max_err*100 << "% on line " << max_err_line << " (Threshold: " <<
            relative_tolerance*100 << "%)" << std::endl;

        if (write_err_ouput)
        {
            //int field_per_line = DEFAULT_FILED_COUNT + num_species;
            write_error_ouput( errs, line_token_count);
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

