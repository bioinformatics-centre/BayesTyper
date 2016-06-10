#include <string>
#include <unordered_map>

#include "VcfFile.hpp"
#include "JoiningString.hpp"

int main(int argc, char const *argv[]) {

	if (argc != 2) {

		std::cout << "USAGE: bayesTyperFilter <bayesTyperVariant.vcf>" << std::endl;
		return 1;
	}

	string vcf_filename(argv[1]);
	GenotypedVcfFileReader vcf_file(vcf_filename, true);

	VcfFileWriter output_vcf(vcf_filename + "_parsed.vcf", vcf_file.metaData(), true);
	Variant * current_variant;

	int vars = 0;

	while (vcf_file.getNextVariant(&current_variant)) {

		vars++;
		output_vcf.write(current_variant);

		if ((vars % 100000) == 0) {

			std::cout << "[" << Utils::getLocalTime() << "] Parsed " << vars << " variants" << endl;
		}

		delete current_variant;
	}

	return 0;
}
