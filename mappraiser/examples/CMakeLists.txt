include(GNUInstallDirs)

# Install the modified toast_so_sim workflow in the prefix

install(PROGRAMS "toast_so_sim.py"
        DESTINATION bin)
