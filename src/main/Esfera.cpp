#include <iostream>


#include "../core/simul.h"


/** Main call */
int main( int args, char* argv[] )
{
	std::cout << "****************************************************" << std::endl;
	
	try
	{
		Simul simu;
		simu.initialize();
		simu.run();
	}
	catch( std::exception &e )
	{
		std::cerr << "\n !!!!!! Error: " << e.what() << std::endl;
		return EXIT_FAILURE;
	}

	std::cout << "****************************************************" << std::endl;
	return EXIT_SUCCESS;
}
