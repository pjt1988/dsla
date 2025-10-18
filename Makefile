SRC_FILES := $(wildcard **/*.cpp)
CXX := g++-12
OBJ_FILES := $(patsubst %.cpp,%.o,$(SRC_FILES))
LDFLAGS := -fopenmp 
INC_DIR := inc
CPPFLAGS := -g -Wall -pedantic -std=c++20 -fopenmp -I$(INC_DIR)
CXXFLAGS := -O0


dsla: $(OBJ_FILES)
	$(CXX) $(LDFLAGS) -I$(INC_DIR) -o $@ $^

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c -o $@ $<

clean:
	rm -f $(OBJ_FILES) dsla
