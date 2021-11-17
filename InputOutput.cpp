//
// Created by maximilian on 30.08.21.
//


#include "InputOutput.h"

/**
 * Parses all Options from command-args
 * @param argc - argc from main methpd
 * @param argv - argv from main method
 * @return Returns a full parameter-struct with all values parsed.
 */
Parameters getParametersFromProgramOptions(char **argc, int argv) {

    namespace po = boost::program_options;
    using namespace std;

    Parameters parameters;

    po::options_description desc("Available options", 100,50);
    desc.add_options()
    ("help,h", "prints this help message")
    ;

    po::options_description description_latticeOptions("Lattice options", 100,50);
    description_latticeOptions.add_options()
            ("sites,N", po::value<unsigned int>(&parameters.sites)->default_value(64), "set number of sites to use (one dimension)")
            ("dimension,D", po::value<unsigned int>(&parameters.dimension)->default_value(2), "set dimension (currently only D=2 is supported)")
            ;
    desc.add(description_latticeOptions);

    po::options_description description_simulationAndMeasurement("Simulation and Measurement", 100,50);
    description_simulationAndMeasurement.add_options()
            ("nos,A", po::value<unsigned int>(&parameters.number_of_simulations)->default_value(10), "Number of simulations")
            ("output,o", po::value<string>(&parameters.outFileName)->default_value("measurements"), "output file prefix to store results")
            ("printAllSteps,p", po::bool_switch(&parameters.flag_printEveryStep), "print observables for each iteration")
            ("measureTime", po::value<unsigned int>(&parameters.measure_time)->default_value(1), "Time-Steps for measurements after system is thermalized")
            ("seed", po::value<unsigned long long>(&parameters.seed)->default_value(std::chrono::high_resolution_clock::now().time_since_epoch().count(), "current time"), "Seed for the random number generator.")
            ;
    desc.add(description_simulationAndMeasurement);

    po::options_description description_thermalizationOptions("Thermalization options", 100,50);
    description_thermalizationOptions.add_options()
            ("timeDelta,t", po::value<TFloat>(&parameters.integrator_timeDelta)->default_value(0.05,"0.05"), "TimeStep for integration")
            ("thermTime1", po::value<unsigned int>(&parameters.thermalization_time_1)->default_value(1500), "Thermalization time-steps: sigma(ModelA), n(ModelA), (symmetryBreaking J active)")
            ("symmetryBreakingJ,S", po::value<TFloat>(&parameters.symmetryBreakingJ)->default_value(0), "External constant field for thermalization phase 1. (If set to zero, the value is +-0.1")
            ("thermTime2", po::value<unsigned int>(&parameters.thermalization_time_2)->default_value(25000), "Thermalization time-steps: sigma(ModelA), n(ModelA), (symmetryBreaking J = 0)")
            ("fastTherm", po::bool_switch(&parameters.fastThermalization), "Use larger time-steps for the first half of thermalization phase 2")
            ("dynamicChangeMode", po::value<unsigned int>(&parameters.dynamicChangeMode)->default_value(DYNAMIC_CHANGE_MODE_RE_INIT), "Specify how to treat nu-field on dynamic change\n0: Re-Init nu-field\n1: Subtract average\n2: Set field to zero")
            ("thermTime3", po::value<unsigned int>(&parameters.thermalization_time_3)->default_value(25000), "Thermalization time-steps: sigma(ModelA), n(ModelB)")
            ;
    desc.add(description_thermalizationOptions);

    po::options_description description_parameters("Physical parameters", 100,50);
    description_parameters.add_options()
            ("temperature,T", po::value<TFloat>(&parameters.temperature)->default_value(4), "physical parameter: Temperature")
            ("couplingConstant,C", po::value<TFloat>(&parameters.C)->default_value(0.1,"0.1"), "Coupling constant to couple both fields")
            ("parJ,J", po::value<TFloat>(&parameters.J)->default_value(0), "External constant field")
            ("msqA", po::value<TFloat>(&parameters.msq_A)->default_value(-1), "Square mass for Model-A field")
            ("lambdaA", po::value<TFloat>(&parameters.lambda_A)->default_value(1), "Lambda (quartic term) for Model-A field")
            ("gammaA", po::value<TFloat>(&parameters.gamma_A)->default_value(0.1, "0.1"), "Langevin-Damping for Model-A field")
            ("gammaB", po::value<TFloat>(&parameters.gamma_B)->default_value(0.1, "0.1"), "Langevin-Damping for Model-B field")
            ("msqB", po::value<TFloat>(&parameters.msq_B)->default_value(1), "Square mass for Model-B field")
            ("mu", po::value<TFloat>(&parameters.mu)->default_value(1), "Mobility coefficient for n-Field")
            ;
    desc.add(description_parameters);

    po::variables_map vm;
    po::store(po::parse_command_line(argv, argc, desc), vm);
    po::notify(vm);

    if( vm.count( "help" ) ) { std::cout << desc << std::endl << std::endl; exit(0); }

    if(parameters.dynamicChangeMode > 2) {std::cerr << "Error: Invalid dynamicChangeMode" << std::endl; exit(1);}

    parameters.volume = parameters.dimension == 3 ? POWER_THREE(parameters.sites) : POWER_TWO(parameters.sites);
    parameters.standard_deviation_A = std::sqrt(2* parameters.gamma_A * parameters.temperature);
    parameters.standard_deviation_B = std::sqrt(2* parameters.gamma_B * parameters.temperature * parameters.mu);
    parameters.integrator_timeDelta_SQRT = std::sqrt(parameters.integrator_timeDelta);

    return parameters;
}

/**
 * Pretty output for console
 * @param observables
 */
void printObservables(Observables observables) {
    std::cout

    << "sigma: "
            << std::setfill(' ') << std::setw(8)
            << std::fixed << std::setprecision(5)
    << observables.sigma_avg
    << "   sigma^2: "
            << std::setfill(' ') << std::setw(8)
            << std::fixed << std::setprecision(5)
    << observables.sigma_sq_avg
    << "      n: "
            << std::setfill(' ') << std::setw(8)
            << std::fixed << std::setprecision(5)
    << observables.n_avg
    << "   n^2: "
            << std::setfill(' ') << std::setw(8)
            << std::fixed << std::setprecision(5)
    << observables.n_sq_avg

    << std::endl;
}

/**
 * Pretty output for console
 * printed before start to show selected options
 * @param parameters
 */
void printParameters(const Parameters& parameters) {
    std::cout << "grid: ";
    for(int k = 0; k < parameters.dimension-1; k++) {
        std::cout << parameters.sites << "x";
    }
    std::cout << parameters.sites << "  Simulations: " << parameters.number_of_simulations
    << "  Seed: " << parameters.seed
    << std::endl;

    std::cout << "ThermTime (1): " << parameters.thermalization_time_1
    << "  ThermTime (2): " << parameters.thermalization_time_2
    << "  ThermTime (3): " << parameters.thermalization_time_3
    << "\nMeasureTime: " << parameters.measure_time << "  TimeStep: " << parameters.integrator_timeDelta
    << "  fastTherm: ";
    if(parameters.fastThermalization) {
        std::cout << "On";
    }else {
        std:: cout << "Off";
    };
    std::cout << std::endl;

    std::cout << "T=" << parameters.temperature << "  C=" << parameters.C
    << "  J=" << parameters.J << "  symmetryBreakingJ=";
    if(parameters.symmetryBreakingJ == 0) {
        std::cout << "random";
    }else {
        std::cout << parameters.symmetryBreakingJ;
    }
    std::cout << std::endl;
}

/**
 * Format time-millis
 * @param milliSeconds
 */
void printRuntime(unsigned long long milliSeconds) {
    if(milliSeconds < 10000) {
        std::cout << "Runtime: " << milliSeconds << " ms" << std::endl;
    }else if(milliSeconds < 1000 * 60) {
        std::cout << "Runtime: " << milliSeconds/1000 << " seconds" << std::endl;
    }else {
        unsigned long long hr = milliSeconds / 3600000;
        milliSeconds -= 3600000 * hr;
        unsigned long long min = milliSeconds / 60000;
        milliSeconds -= 60000 * min;
        unsigned long long sec = milliSeconds / 1000;

        std::cout << "Runtime: "
        << std::setfill('0') << std::setw(2) << hr
        << ":"
        << std::setfill('0') << std::setw(2) << min
        << ":"
        << std::setfill('0') << std::setw(2) << sec
        << std::endl;
    }
}

void writeListToFileLine(TFloat parameter_T, TFloat parameter_C, const std::string& name, const std::vector<TFloat>& floatList) {
    std::ofstream outfile;
    outfile.open(name, std::ios_base::app); // append instead of overwrite
    outfile << parameter_T << " " << parameter_C;
    for(auto obs : floatList) {
        outfile << " " << obs;
    }
    outfile << std::endl;
    outfile.close();
}

void writeParametersToFile(std::basic_ofstream<char> &outfile, const Parameters &parameters) {
    //Write parameters
    outfile << "\"N\": " << parameters.sites
            << ", \"T\": " << parameters.temperature
            << ", \"C\": " << parameters.C
            << ", \"D\": " << parameters.dimension
            << ", \"t\": " << parameters.integrator_timeDelta
            << ", \"msq_A\": " << parameters.msq_A
            << ", \"lambda_A\": " << parameters.lambda_A
            << ", \"gamma_A\": " << parameters.gamma_A
            << ", \"msq_B\": " << parameters.msq_B
            << ", \"gamma_B\": " << parameters.gamma_B
            << ", \"J\": " << parameters.J
            << ", \"mu\": " << parameters.mu
            << ", \"fastTherm\": " << (parameters.fastThermalization ? "true" : "false")
            << ", \"thermTime1\": " << parameters.thermalization_time_1
            << ", \"thermTime2\": " << parameters.thermalization_time_2
            << ", \"thermTime3\": " << parameters.thermalization_time_3
            << ", \"dynamicChangeMode\": " << (parameters.dynamicChangeMode == DYNAMIC_CHANGE_MODE_RE_INIT ?"\"ReInit\"" :
            parameters.dynamicChangeMode == DYNAMIC_CHANGE_MODE_ZERO ? "\"Zero\"" : "\"Subtract\"")
            << ", \"measure_time\": " << parameters.measure_time
            << ", \"p\": " << (parameters.flag_printEveryStep ? "true" : "false")
            << ", \"seed\": " << parameters.seed
            << ", \"A\": " << parameters.number_of_simulations
            << ", \"S\": " << parameters.symmetryBreakingJ;
}

//void writeAllSteps(const Parameters& parameters, const std::string& name, const std::vector<std::vector<TFloat>>& stepList) {
void writeAllSteps(const Parameters& parameters, const std::map<std::string, std::vector<TFloat>>& stepMap) {
//    writeListToFileLine(parameters.temperature, parameters.C, parameters.outFileName + "-allSteps_" + name + ".txt", floatList);
    std::ofstream outfile;
    outfile.open(parameters.outFileName + "-allSteps.yaml", std::ios_base::app); // append instead of overwrite

    //Begin Json-like object
    outfile << "{";

    writeParametersToFile(outfile, parameters);

    //Write observable arrays
    for (const auto& entry: stepMap) {
        outfile << ", \"" << entry.first << "\": [";
        for(auto obs : (entry.second)) {outfile << obs << ", ";}
        outfile << "]";
    }

    //Close Json-like object
    outfile << "}," << std::endl;
    outfile.close();
}

/**
 * Write observables to files.
 * The format is JSON-like.
 * Although you should use a yaml-parser to load the file in other programs
 * @param parameters
 * @param observablesList
 */
void writeToDatabase(Parameters parameters, const std::vector<Observables>& observablesList) {

    std::ofstream outfile;
    outfile.open(parameters.outFileName + ".yaml", std::ios_base::app); // append instead of overwrite

    //Begin Json-like object
    outfile << "{";

    parameters.flag_printEveryStep = false;
    parameters.measure_time = 1;
    writeParametersToFile(outfile, parameters);

    //Write observable arrays
    outfile << ",\"sigma\": [";
    for(auto obs : observablesList) {outfile << obs.sigma_avg << ", ";}
    outfile << "], \"n\": [";
    for(auto obs : observablesList) {outfile << obs.n_avg << ", ";}
    outfile << "], \"sigma_sq\": [";
    for(auto obs : observablesList) {outfile << obs.sigma_sq_avg << ", ";}
    outfile << "], \"n_sq\": [";
    for(auto obs : observablesList) {outfile << obs.n_sq_avg << ", ";}
    outfile << "],";

    //Close Json-like object
    outfile << "}," << std::endl;
    outfile.close();
}
