###############################################################################
################### MOOSE Application Standard Makefile #######################
###############################################################################
#
# Optional Environment variables
# MOOSE_DIR        - Root directory of the MOOSE project 
#
###############################################################################
# Use the MOOSE submodule if it exists and MOOSE_DIR is not set

MOOSE_SUBMODULE    := $(CURDIR)/moose
ifneq ($(wildcard $(MOOSE_SUBMODULE)/framework/Makefile),)
  MOOSE_DIR        ?= $(MOOSE_SUBMODULE)
else
  MOOSE_DIR        ?= $(shell dirname `pwd`)/moose
endif

# framework
FRAMEWORK_DIR      := $(MOOSE_DIR)/framework
include $(FRAMEWORK_DIR)/build.mk
include $(FRAMEWORK_DIR)/moose.mk

################################## MODULES ####################################
ALL_MODULES := yes
include $(MOOSE_DIR)/modules/modules.mk
###############################################################################

# dep apps
APPLICATION_DIR    := $(CURDIR)
APPLICATION_NAME   := dgosprey
BUILD_EXEC         := yes
DEP_APPS           := $(shell $(FRAMEWORK_DIR)/scripts/find_dep_apps.py $(APPLICATION_NAME) )
include            $(FRAMEWORK_DIR)/app.mk

###############################################################################
# Additional special case targets should be added here
#export EXTERNAL_FLAGS += -std=c++11 -fpermissive -Wunused-parameter -Wc++11-extensions -Wwrite-strings -Wconversion-null -Wunused-variable
#export libmesh_CXXFLAGS+= -std=c++11 -fpermissive -Wunused-parameter -Wc++11-extensions -Wwrite-strings -Wconversion-null -Wunused-variable

