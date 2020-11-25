//------------------------------------------------------------------------------
//  The @header2, @begin_table, @item3 and @end_table commands are custom
//  defined commands in Doxygen.in. They are defined under ALIASES. For the page
//  created here, the 80 column limit is exceeded. Arguments of aliases are
//  separated by ','. If you intended ',' to be a string you must use an escaped
//  comma '\,'.
//
///  @page wout_diff_cl_parsing_sec Command Line Arguments
///
///  @tableofcontents
///
///  @section wout_diff_cl_parsing_intro Introduction
///  This contains a description of the command line arguments. All arguments
///  take the form of
///
///  @fixed_width{-arg=value}
///
///  @section wout_diff_cl_parsing_arg_sec Command Line Arguments
///  @header2{Argument, Takes Value, Discription}
///  @begin_table
///     @item3{@fixed_width{-h},          N, Displays the help text and exits the program.}
///     @item3{@fixed_width{-wout_file1}, Y, First wout file name.}
///     @item3{@fixed_width{-wout_file2}, Y, Second wout file name.}
///     @item3{@fixed_width{-quantity},   Y, Wout quantity to check.}
///     @item3{@fixed_width{-tol},        Y, Tolarance value.}
///  @end_table
///
///  @section wout_diff_cl_pasring_prog_ref_sec Programmers Reference
///  Reference material for the coding to implement command line parsing is found
///  in the @ref vmec_test::commandline_parser class.
//------------------------------------------------------------------------------
//******************************************************************************
///  @file wout_diff.cpp
///  @brief Utility to check the difference between two wout files.
//******************************************************************************

#include <iostream>

#include "vmec_test_commandline_parser.hpp"

//------------------------------------------------------------------------------
///  @brief Main test program.
///
///  @param[in] argc Number of commandline arguments.
///  @param[in] argv Array of arguments strings.
//------------------------------------------------------------------------------
int main(int argc, const char * argv[]) {
    const wout_test::commandline_parser args(argc, argv, []() -> void {
//                   "                                        ''                                      "
        std::cout << "                                                                                " << std::endl;
        std::cout << "                                     DIFF WOUT                                  " << std::endl;
        std::cout << "                                                                                " << std::endl;
        std::cout << "Usage: xdiff_wout [-arg][=option] ...                                           " << std::endl;
        std::cout << "                                                                                " << std::endl;
        std::cout << "Options:                                                                        " << std::endl;
        std::cout << "All options are displayes as [arg][takesoption][Discription]                    " << std::endl;
        std::cout << "  -h          N Display this information.                                       " << std::endl;
        std::cout << "                                                                                " << std::endl;
        std::cout << "  -wout_file1 Y First wout file name.                                           " << std::endl;
        std::cout << "                                                                                " << std::endl;
        std::cout << "  -wout_file2 Y Second wout file name.                                          " << std::endl;
        std::cout << "                                                                                " << std::endl;
        std::cout << "  -quantity   Y Wout quantity to check.                                         " << std::endl;
        std::cout << "                                                                                " << std::endl;
        std::cout << std::endl;

        exit(0);
    });
}
