#include <iostream>
#include <home/vdwti/HPC_21_22/Boost/program_options.hpp>
namespace po = boost::program_options;

int main(int argc, char* argv[]) {
    po::options_description opts(
    "Available options.");
    opts.add_options()
    ("start", po::value<int>()->default_value(0),
    "Starting value.")
    ("end", po::value<int>()->default_value(10),
    "Ending value.")
    ("help", "Print help message.");

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, opts), vm);
    po::notify(vm);

    if (vm.count("help")) {
    std::cout << opts << std::endl;
    }
    const int start = vm["start"].as<int>();
    const int end = vm["end"].as<int>();
 }