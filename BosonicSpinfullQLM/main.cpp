#include "includes/BosonicSpinfullQLM_UI.hpp"


int main(const int argc, char* argv[]) {
	std::unique_ptr<BoseUI> intface = std::make_unique<BosonicSpinfullQLM_UI::ui>(argc, argv);
	intface->make_sim();
	return 0;
}