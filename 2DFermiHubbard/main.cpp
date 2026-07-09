#include "includes/FermiHubbard_UI.hpp"


int main(const int argc, char* argv[]) {
	std::unique_ptr<FermiUI> intface = std::make_unique<FermiHubbard_UI::ui>(argc, argv);
	intface->make_sim();
	return 0;
}