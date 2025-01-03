#!/bin/bash
# Script pour construire le projet

# Définir les chemins relatifs
ROOT_DIR=$(pwd)
BUILD_DIR="${ROOT_DIR}/builds"
SRC_DIR="${ROOT_DIR}/src"
LOG_DIR="${ROOT_DIR}/logs"

# Crée les répertoires nécessaires
mkdir -p "${BUILD_DIR}"
mkdir -p "${LOG_DIR}"

# Crée dynamiquement le fichier CMakeLists.txt dans le répertoire builds
cat <<EOL > "${BUILD_DIR}/CMakeLists.txt"
cmake_minimum_required(VERSION 3.1)

project(levelset)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

set(SRCS
    ${SRC_DIR}/levelset.cpp
    ${SRC_DIR}/main_loop.cpp
)

include_directories(${SRC_DIR})

find_package(NEOS REQUIRED)

add_executable(levelset \${SRCS})

target_link_libraries(levelset neos::neos)

set_target_properties(levelset PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${ROOT_DIR}
)
EOL

# Aller dans le répertoire builds
cd "${BUILD_DIR}" || exit

# Exécuter cmake et enregistrer les logs
cmake . > "${LOG_DIR}/cmake_log.txt" 2>&1

# Exécuter make et enregistrer les logs
make > "${LOG_DIR}/make_log.txt" 2>&1

# Revenir au répertoire principal
cd "${ROOT_DIR}" || exit

echo "Build terminé. Les logs sont dans ${LOG_DIR}."
