## no implicit rules
.SUFFIXES: 
FORTEXT:=f90

## definitions
FC=gfortran
COMPILE.f08 = $(FC) $(FCFLAGS) $(TARGET_ARCH) -c
MAKEMOD.f08 = $(FC) $(FCFLAGS) $(TARGET_ARCH) -fsyntax-only -c

## mod and smod directory
MODDIR := .mod
ifneq ($(MODDIR),)
  $(shell test -d $(MODDIR) || mkdir -p $(MODDIR))
  FCFLAGS+= -J $(MODDIR)
endif

# $(call source-to-extension,source-file-list,new-extension)
define source-to-extension
  $(strip \
    $(foreach ext,$(FORTEXT),\
      $(subst .$(ext),.$2,$(filter %.$(ext),$1))))
endef
# $(call modsource-pattern-rule,extension)
define modsource-pattern-rule
%.anc: %.$1
	$$(MAKEMOD.f08) $$<
	@touch $$@
endef

## sources for the program
SOURCES:=AnalyzerExcitonTransfer.f90 \
	 countNumberOfLine.f90 \
     indexage.f90 \
     parseTrajOrGeoxyz.f90\
     fileGenerator.f90 \
     paramAllocatableTab.f90 \
     sortGeoFile.f90 \
     parseIonsAndWC.f90 \
     InitGeoForLen.f90 \
     lenghtenMonomer.f90 \
     initGeoForStacking.f90 \
     translateGeoForStackV2.f90 \
     stackMonomer.f90 \
     mod_analyzeTraj.f90 \
     fillSetB.f90 \
     analyzeGeoTraj1.f90 \
     mod_analyzeGeoTraj.f90 \
     mod_excitonDiffusion.f90 \
     excitonDiff.f90
OBJECTS:=$(call source-to-extension,$(SOURCES),o)
ANCHORS:=$(call source-to-extension,$(SOURCES),anc)

## main and clean targets
#main: $(OBJECTS)
AnalyzerExcitonTransfer: $(OBJECTS)
	$(FC) -o $@ $+

.PHONY: clean
clean:
	-rm -rf *.mod *.smod $(OBJECTS) $(ANCHORS) main
	-test -d $(MODDIR) && rm -r $(MODDIR)

## compilation rules
$(foreach ext,$(FORTEXT),$(eval $(call modsource-pattern-rule,$(ext))))
%.o: %.anc
	$(COMPILE.f08) $(OUTPUT_OPTION) $(wildcard $(addprefix $*.,$(FORTEXT)))
	@touch $@

## target because -fsyntax-only also needs it)
AnalyzerExcitonTransfer.anc : countNumberOfLine.anc \
	indexage.o excitonDiff.anc \
	mod_analyzeGeoTraj.anc analyzeGeoTraj1.anc fileGenerator.anc stackMonomer.anc \
	lenghtenMonomer.anc parseIonsAndWC.anc paramAllocatableTab.anc
parseIonsAndWC.anc : countNumberOfLine.anc fileGenerator.anc
lenghtenMonomer.anc :InitGeoForLen.anc sortGeoFile.anc parseTrajOrGeoxyz.anc \
	fileGenerator.anc indexage.anc
stackMonomer.anc : initGeoForStacking.anc translateGeoForStackV2.anc \
	parseTrajOrGeoxyz.anc sortGeoFile.anc fileGenerator.anc indexage.anc
analyzeGeoTraj1.anc : indexage.anc mod_analyzeTraj.anc fillSetB.anc
excitonDiff.anc : mod_excitonDiffusion.anc fileGenerator.anc
#.PHONY: clean
#clean:
#  -$(RM) *.mod *.smod $(OBJECTS) $(ANCHORS) AnalyzerExcitonTransfer
#  -$(TEST) -d $(MODDIR) && $(RM) -r $(MODDIR)
#  -$(TEST) -d $(MODDIRTMP) && $(RM) -r $(MODDIRTMP)
