#include <iostream>
#include <fstream>
#include <string>
#include <limits>

/**
 * @brief Function to compare the contents of two energy files
 *
 * @param file_a First file to compare
 * @param file_b Second file to compare
 * @param relative_tolerance Relative tolerance which is acceptable
 *
 * @return True is they match (within tol), false if not
 */

bool compare_energies(
        std::string file_a,
        std::string file_b,
        double relative_tolerance,
        int num_lines_to_skip = 3 // Most energy files have 3 lines of padding
)
{

    try {

        bool match = true;

        std::string line1 = "";
        std::string line2 = "";

        std::ifstream f1 (file_a);
        std::ifstream f2 (file_b);

        if (f1.is_open() && f2.is_open())
        {

            // Perform skipping
            for (int i = 0; i < num_lines_to_skip; i++)
            {
                getline(f1,line1);
                getline(f2,line2);
            }

            // Do processing
            while ( getline(f1,line1) )
            {
                getline(f2,line2);

                // Tokenize lines
                std::stringstream linestream1(line1);
                std::string item1;

                std::stringstream linestream2(line2);
                std::string item2;

                while (getline(linestream1, item1, ' '))
                {
                    getline(linestream2, item2, ' ');

                    double A = std::stod(item1);
                    double B = std::stod(item2);

                    bool within_tol = false;

                    std::cout <<
                        item1 << " (" << A << ") " <<
                        " vs " <<
                        item2 << " (" << B << ") " <<
                        std::endl;


                    // Calculate if we're withing tolerances

                    // Right now this is pretty arbitrary..
                    double abs_threshhold = 10 * std::numeric_limits<float>::epsilon();
                    // If we're close to relative, do absolute
                    if (std::abs(std::min(A,B)) < abs_threshhold)
                    {
                        // Finding a relative error to 0 doesn't make much
                        // sense, so lets do absolute error instead
                        if ( std::abs(A-B) < std::numeric_limits<double>::epsilon() )
                        {
                            within_tol = true;
                        }
                        else {
                            within_tol = false;
                        }
                        //break;
                    }

                    else if ((std::abs(A-B) / std::min(A,B)) < relative_tolerance)
                    {
                        within_tol = true;
                    }
                    else {
                        within_tol = false;
                    }

                    if (! within_tol ) {
                        match = false;
                        return match;
                    }
                    std::cout << item1 << std::endl;
                }

                std::cout << line1 << '\n';
                std::cout << line2 << '\n';
            }
            f1.close();
            f2.close();
        }
        else {
            std::cout << "Unable to open file";
            return false;
        }

        return match;
    }
    catch (...) { // Catching all is bad form, but OK for now..
        return false;
    }

}
