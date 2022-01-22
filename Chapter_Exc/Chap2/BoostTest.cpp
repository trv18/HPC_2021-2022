#include <iostream>
#include <string>
#include <boost/program_options.hpp>
namespace po = boost::program_options;

int main(int argc, char* argv[]) {
    po::options_description opts(
    "Available options.");
    opts.add_options()
    ("start,s", po::value<int>()->default_value(0),
    "Starting value.")
    ("end, e", po::value<int>()->default_value(10),
    "Ending value.")
    ("Random", po::value<std::string>() -> default_value("Test"),
    "Testing input.")

    ("help", "Print help message.");

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, opts), vm);
    po::notify(vm);

    if (vm.count("help")) {
    std::cout << opts << std::endl;
    }
    const int start = vm["start"].as<int>();
    const int end = vm["end"].as<int>();

    std::cout << "Chosen Start Value = " << start << std::endl;
 }