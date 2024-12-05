#ifndef INCLUDED_basic_assert_hh
#define INCLUDED_basic_assert_hh

#include <iostream>

namespace basic {

#define ALWAYS_ASSERT(_Expression) \
	if ( !(_Expression) ) { \
		std::cerr << std::endl; \
		std::cerr << "FAILED ASSERTION: " << std::endl << #_Expression << std::endl << std::endl; \
		std::cerr << "EXIT FROM: "        << std::endl << __FILE__ << " line: " << __LINE__ << std::endl << std::endl; \
		std::exit(-1); \
	}

#define ALWAYS_ASSERT_MSG(_Expression,_MSG) \
	if ( !(_Expression) ) { \
		std::cerr << std::endl; \
		std::cerr << "FAILED ASSERTION: " << std::endl << #_Expression << std::endl << std::endl; \
		std::cerr << "MESSAGE: "          << std::endl << _MSG << std::endl << std::endl; \
		std::cerr << "EXIT FROM: "        << std::endl << __FILE__ << " line: " << __LINE__ << std::endl << std::endl; \
		std::exit(-1); \
	}

#define EXIT_WITH_MSG(_MSG) \
    { \
		std::cerr << "EXIT WITH MESSAGE: "          << std::endl << _MSG << std::endl << std::endl; \
		std::cerr << "EXIT FROM: "        << std::endl << __FILE__ << " line: " << __LINE__ << std::endl << std::endl; \
		std::exit(-1); \
    }

}

#endif
