# Add hunter install directory to the find_package variables
if (EXISTS "${CMAKE_CURRENT_LIST_DIR}/_3rdParty/Hunter/install-root-dir")
  file(READ "${CMAKE_CURRENT_LIST_DIR}/_3rdParty/Hunter/install-root-dir" HunterInstall)
  list(APPEND CMAKE_FIND_ROOT_PATH "${HunterInstall}")
  list(APPEND CMAKE_PREFIX_PATH "${HunterInstall}")
endif()
