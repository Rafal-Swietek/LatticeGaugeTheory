#include "includes/FermionicQLM_UI.hpp"


int main(const int argc, char* argv[]) {
	std::unique_ptr<FermiUI> intface = std::make_unique<FermionicQLM_UI::ui>(argc, argv);
	intface->make_sim();
	return 0;
}