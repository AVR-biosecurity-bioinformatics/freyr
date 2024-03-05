##### this loads the listed R scripts 

# "projectDir" needs to be defined in the sys.source call

# script_list <-
#     c(
#         "jack_notes/.Rprofile",
#         "bin/functions.R",
#         "bin/themes.R",
#         "bin/_targets_packages.R"
#     )

# print(script_list)

# for (script in script_list) {
#     source(file.path(projectDir, script))
# }

source(file.path(projectDir, "jack_notes/.Rprofile"))
source(file.path(projectDir, "bin/functions.R"))
source(file.path(projectDir, "bin/themes.R"))
source(file.path(projectDir, "bin/_targets_packages.R"))