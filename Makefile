CXXFLAGS = -O3

Mirage: main.cpp inverse_potts_model.cpp

	$(CXX) $(CXXFLAGS) -o Ipm main.cpp inverse_potts_model.cpp -std=c++11  
