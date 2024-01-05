SRC_FILES := $(wildcard **/*.cpp)
OBJ_FILES := $(patsubst %.cpp,%.o,$(SRC_FILES))
LDFLAGS := -fopenmp 
INC_DIR := inc
CPPFLAGS := -g -Wall -pedantic -std=c++20 -fopenmp -I$(INC_DIR)
CXXFLAGS := -O0


dsla: $(OBJ_FILES)
	g++-12 $(LDFLAGS) -I$(INC_DIR) -o $@ $^

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp
	g++-12 $(CPPFLAGS) $(CXXFLAGS) -c -o $@ $<

clean:
	rm -f $(OBJ_FILES) dsla

