g++ -shared -std=c++11 -fPIC `python -m pybind11 --includes` pydssp.cpp -o dssp`python3-config --extension-suffix` dssp.cpp iocif.cpp mas.cpp primitives-3d.cpp structure.cpp utils.cpp -lboost_date_time -lboost_filesystem -lboost_iostreams -lboost_program_options  -lboost_system -lboost_thread -lpthread