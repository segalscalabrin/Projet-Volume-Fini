# Variables
G++ = g++
CXXFLAGS = -Wall -O3 -I src               # Options de compilation
SRC = $(wildcard src/**/*.cpp src/*.cpp)  # Recherche des fichiers sources
OBJ = $(patsubst src/%.cpp, build/%.o, $(SRC))  # Fichiers objets dans 'build'
OUT = levelset                            # Nom de l'exécutable principal
TEST_OUT = levelset_test                  # Nom de l'exécutable de tests

# Compilation principale
all: $(OUT)

$(OUT): $(OBJ)
	$(G++) $(CXXFLAGS) $(OBJ) -o $(OUT)

# Compilation des tests
test: $(TEST_OUT)

$(TEST_OUT): $(OBJ)
	$(G++) $(CXXFLAGS) $(wildcard tests/*.cpp) $(OBJ) -o $(TEST_OUT)

# Règle pour compiler les fichiers .o dans build/
build/%.o: src/%.cpp
	@mkdir -p $(dir $@)
	$(G++) $(CXXFLAGS) -c $< -o $@

# Nettoyage
clean:
	rm -rf build $(OUT) $(TEST_OUT)
