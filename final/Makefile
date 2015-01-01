SHELL := /bin/bash
.SUFFIXES:
.PHONY: all clean dist test test-% release report dist handin
T=ent
SRC := $(wildcard src/*.cpp)
OBJ := $(addprefix build/,$(subst src/,,$(SRC:.cpp=.o)))
CXXFLAGS += -g3 -Wall -Wextra -std=c++11

DEPEND := $(wildcard *.tex) img
TEX := xelatex -shell-escape -interaction=nonstopmode
COMMON_SRC := src/*.cc web
RELEASE_SRC := $(COMMON_SRC) report/report.pdf Makefile
DIST_SRC := $(COMMON_SRC) *.tex Makefile

NULL :=
SPACE := $(NULL) $(NULL)
remove-space = $(subst $(SPACE),,$1)
notify = notify-send xelatex $1

all: report/report.pdf
#all: $T

sinclude $(OBJ:.o=.d)

$T: $(OBJ)
	$(LINK.cc) $^ $(LDLIBS) -o $@

build/%.o: src/%.cpp | build
	g++ -std=c++11 -MM -MP -MT $@ -MF $(@:.o=.d) $<
	$(COMPILE.cc) $(OUTPUT_OPTION) $<

build:
	-mkdir $@

clean:
	$(RM) -r build

test-%: tests/%
	$(MAKE) -C $? test

report report/report.pdf: $(addprefix report/,report.tex 1.tex _title.tex)
	cd report && $(TEX) $(notdir $<)

report/report.html: report/1.tex
	pandoc -f latex --mathjax --base-header-level 2 <(sed '/\\title{/d' $<) -o $@

report/%.png:report/%.dot
	$(DOT) $< -o $@

handin:
	cd src && tar vcf /tmp/handin.tar *.{cpp,h} --exclude main.cpp

release: report/report.pdf
ifeq ($(words $(RELEASE_SRC)),1)
	cp -a $(RELEASE_SRC) "/tmp/$(call remove-space,$(notdir $(CURDIR)))-release.pdf"
else
	tar zcvf "/tmp/$(call remove-space,$(notdir $(CURDIR)))-release.tar.gz" $(RELEASE_SRC) --exclude '.*' --transform 's,report/report.pdf,report.pdf,' --transform s,^,final/,
endif

watch:
	(inotifywait -mqe CLOSE_WRITE report/*.tex | while read; do make report && $(call notify,success) || $(call notify,fail); done) &

help:
	@echo 'make COMMAND'
	@echo
	@echo 'COMMAND:'
	@echo '  all       build output/$(TARGET).pdf'
	@echo '  release   copy pdf to /tmp/ if there is only one target, otherwise make a tarball'
	@echo "  clean     remove \`output'"
	@echo '  dist      backup'
	@echo '  watch     open pdf with llpp and rebuild target upon changes'
