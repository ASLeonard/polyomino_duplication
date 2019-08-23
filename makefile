MAKEFLAGS+="-j $(nproc)"

#Compiler and Linker
CXX         := g++-9

#The Target Binary Program
Du_TARGET   := NEW_DuplicationEvolution


#The Directories, Source, Includes, Objects, Binary and Resources
SRCDIR      := src
INCDIR      := includes
LIBDIR      := polyomino_core
BUILDDIR    := build
TARGETDIR   := bin
SRCEXT      := cpp
DEPEXT      := d
OBJEXT      := o

#VPATH=src:polyomino/src

#Flags, Libraries and Includes
CXXFLAGS    := -std=gnu++2a -Wall -Wextra -pedantic -pipe -march=haswell -no-pie $(cmdflag)
ifndef DEBUG
CXXFLAGS += -O3 -fopenmp -flto -flto-partition=none -ffunction-sections -fdata-sections
else
CXXFLAGS += -g3 -O0
endif

ifdef ILen
CXXFLAGS += -DGCC_INTERFACE_LENGTH=$(ILen)
Du_TARGET := $(Du_TARGET)_L$(ILen)
endif

ifdef FW
CXXFLAGS += -DFULL_WRITE=1
endif

INC         := -I$(INCDIR) -I$(LIBDIR)/$(INCDIR)
INCDEP      := -I$(INCDIR) -I$(LIBDIR)/$(INCDIR)

#---------------------------------------------------------------------------------
#DO NOT EDIT BELOW THIS LINE
#---------------------------------------------------------------------------------
Du_SOURCES := $(shell find $(SRCDIR) -type f -name duplication_*.$(SRCEXT))
Du_OBJECTS     := $(patsubst $(SRCDIR)/%,$(BUILDDIR)/%,$(Du_SOURCES:.$(SRCEXT)=.$(OBJEXT)))


#Default Make
all: Du

#Clean only Objects
clean:
	@$(RM) -rf $(BUILDDIR)

#Pull in dependency info for *existing* .o files

-include $(Du_OBJECTS:.$(OBJEXT)=.$(DEPEXT))
-include $(CORE_OBJECTS:.$(OBJEXT)=.$(DEPEXT))

Du: $(Du_OBJECTS) $(CORE_OBJECTS) 
	@mkdir -p $(TARGETDIR)
	$(CXX) $(CXXFLAGS) -Wl,--gc-sections -o $(TARGETDIR)/$(Du_TARGET) $^


#Compile
$(BUILDDIR)/%.$(OBJEXT): $(SRCDIR)/%.$(SRCEXT)
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) $(INC) -c -o $@ $<
	@$(CXX) $(CXXFLAGS) $(INCDEP) -MM $(SRCDIR)/$*.$(SRCEXT) > $(BUILDDIR)/$*.$(DEPEXT)
	@cp -f $(BUILDDIR)/$*.$(DEPEXT) $(BUILDDIR)/$*.$(DEPEXT).tmp
	@sed -e 's|.*:|$(BUILDDIR)/$*.$(OBJEXT):|' < $(BUILDDIR)/$*.$(DEPEXT).tmp > $(BUILDDIR)/$*.$(DEPEXT)
	@sed -e 's/.*://' -e 's/\\$$//' < $(BUILDDIR)/$*.$(DEPEXT).tmp | fmt -1 | sed -e 's/^ *//' -e 's/$$/:/' >> $(BUILDDIR)/$*.$(DEPEXT)
	@rm -f $(BUILDDIR)/$*.$(DEPEXT).tmp


#Non-File Targets
.PHONY: all clean Du check-and-reinit-submodules

check-and-reinit-submodules:
	@if git submodule status | egrep -q '^[-]|^[+]' ; then \
		echo "INFO: Need to reinitialize git submodules"; \
		git submodule update --init; \

	fi
