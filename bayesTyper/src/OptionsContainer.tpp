#include <string>
#include <vector>
#include <sstream>
#include <unordered_map>

#include "boost/algorithm/string/split.hpp"
#include "boost/algorithm/string.hpp"
#include "boost/lexical_cast.hpp"

#include "Utils.hpp"
#include "Regions.hpp"

OptionsContainer::OptionsContainer(std::string program_version, std::string start_time) {

    OptionValue<std::string> * option_value_version = new OptionValue<std::string>(program_version);
    assert(options.emplace("program-version", std::pair<OptionValueBase*, std::string>(option_value_version, program_version)).second);


    OptionValue<std::string> * option_value_time = new OptionValue<std::string>(start_time);
    assert(options.emplace("start-time", std::pair<OptionValueBase*, std::string>(option_value_time, start_time)).second);
}


OptionsContainer::~OptionsContainer() {

    for (auto &option: options) {

        delete option.second.first;
    }
}


template<typename ValueType>
void OptionsContainer::parseValue(std::string option, ValueType value) {

    OptionValue<ValueType> * option_value = new OptionValue<ValueType>(value);

    std::stringstream value_stream;
    value_stream << value;   

    assert(options.emplace(option, std::pair<OptionValueBase*, std::string>(option_value, value_stream.str())).second);
}


void OptionsContainer::parseFiles(std::string option, std::string files) {

    std::vector<std::string> file_names;
    split(file_names, files, boost::is_any_of(","));

    OptionValue<std::vector<std::string> > * option_value = new OptionValue<std::vector<std::string> >(file_names);
    
    assert(option_value->value.size() >= 1);
    assert(options.emplace(option, std::pair<OptionValueBase*, std::string>(option_value, files)).second);
}

void OptionsContainer::parseRegions(std::string option, std::string regions) {

    OptionValue<Regions> * option_value = new OptionValue<Regions>(Regions(regions));
    assert(options.emplace(option, std::pair<OptionValueBase*, std::string>(option_value, regions)).second);
}


void OptionsContainer::parsePriors(std::string option, std::string priors) {

    OptionValue<std::vector<std::pair<double,double> > > * option_value = new OptionValue<std::vector<std::pair<double,double> > >(parsePriorString(priors));

    assert(option_value->value.size() >= 1);
    assert(options.emplace(option, std::pair<OptionValueBase*, std::string>(option_value, priors)).second);
}


std::vector<std::pair<double,double> > OptionsContainer::parsePriorString(std::string & prior_string) {

    std::vector<std::pair<double,double> > prior_info;

    std::vector<std::string> prior_string_splittet;
    split(prior_string_splittet, prior_string, boost::is_any_of(":"));

    for (uint i = 0; i < prior_string_splittet.size(); i++) {

        std::vector<std::string> prior_parameters;
        boost::split(prior_parameters, prior_string_splittet.at(i), boost::is_any_of(","));
        assert(prior_parameters.size() == 2);
        prior_info.emplace_back(boost::lexical_cast<double>(prior_parameters.front()), boost::lexical_cast<double>(prior_parameters.back()));
    }  

    return prior_info;  
}


template<typename ValueType>
const ValueType & OptionsContainer::getValue(std::string option) const {

    return static_cast<OptionValue<ValueType> * >(options.at(option).first)->value;
}

std::string OptionsContainer::writeVCFFileHeader() {

    std::stringstream header_stream;

    header_stream << "##CommandLine=<ID=BayesTyper, CommandLineOptions=\"";

    auto oit = options.cbegin();

    header_stream << oit->first << "=" << "\"" << oit->second.second << "\"";
    oit++;

    while (oit != options.cend()) {

        header_stream << ", " << oit->first << "=" << "\"" << oit->second.second << "\"";
        oit++;
    }

    header_stream << "\">\n";

    return header_stream.str();
}



