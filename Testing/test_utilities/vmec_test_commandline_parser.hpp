//******************************************************************************
///  @file vmec_test_commandline_parser.hpp
///  @brief Contains classes to interpolate full and half grid quanities.
//******************************************************************************

#ifndef vmec_test_commandline_parser_hpp
#define vmec_test_commandline_parser_hpp

#include <map>
#include <string>
#include <sstream>
#include <functional>

namespace vmec_test {

///  Type for the argument map.
    typedef std::map<std::string, std::string> arg_map;
///  Type for a key value pair.
    typedef std::pair<std::string, std::string> arg_element;

//------------------------------------------------------------------------------
///  @brief A radial quantity.
//------------------------------------------------------------------------------
    class commandline_parser {
    public:
///  Parsed commands.
        const arg_map commands;
///  Help callback function.
        const std::function<void(void)> help;

//------------------------------------------------------------------------------
///  @brief Factory method to parse the commandline and produce the arguments.
///
///  @param[in] argc Number of commandline arguments.
///  @param[in] argv Commandline strings.
///  @param[in] help Call back function to display the help message.
///  @returns A constructed map of commandline argument key value pairs.
//------------------------------------------------------------------------------
        static arg_map parse_commands(const size_t argc,
                                      const char * argv[],
                                      const std::function<void(void)> help) {
            if (argc == 0) {
                help();
            }

            arg_map commands;

            for (size_t i = 1; i < argc; i++) {
                std::string arg(argv[i]);

                if (arg == "-h") {
                    help();
                } else {
                    size_t eqpos = arg.find('=');
                    if (eqpos != std::string::npos) {
                        std::string key = arg.substr(0, eqpos);
                        std::string value = arg.substr(eqpos + 1, std::string::npos);

                        commands.insert(arg_element(key, value));
                    } else {
                        commands.insert(arg_element(arg, ""));
                    }
                }
            }

            return commands;
        }

//------------------------------------------------------------------------------
///  @brief Construct a commandline_parser object by pasring command line
///  arguments.
///
///  @param[in] argc Number of commandline arguments.
///  @param[in] argv Commandline strings.
///  @param[in] help Call back function to display the help message.
//------------------------------------------------------------------------------
        commandline_parser(const size_t argc,
                           const char * argv[],
                           const std::function<void(void)> help) :
        commands(parse_commands(argc, argv)) {}

//------------------------------------------------------------------------------
///  @brief Check if command arg was set.
///
///  @param[in] key Commandline key to check.
///  @returns True if the key was set.
//------------------------------------------------------------------------------
        bool is_set(const std::string &key) const {
            return commands.find(key) != commands.end();
        }

//------------------------------------------------------------------------------
///  @brief Get the value of the agument.
///
///  @param[in] key Commandline key to check.
///  @returns Value of the argument.
//------------------------------------------------------------------------------
        template<typename TYPE>
        TYPE get(const std::string &key) const {
            if (!is_set(key)) {
                help();
            }

            std::stringstream value_stream(commands.at(key));
            TYPE temp;
            value_stream >> temp;

            return temp;
        }
    };
}

#endif
