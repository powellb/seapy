# Build seapy modules, tag file, and documentation
.PHONY: etags
etags:
	rm -f TAGS
	find . '!' '('  -name __init__.py ')' -name "*\.py" | xargs etags -a --language=python

oalib.so: src/oalib.F
	cd src; $(MAKE) $(MFLAGS)

doc: force_build
	cd doc; $(MAKE) $(MFLAGS) html && ln -fs _build/html html

force_build:
	true

all: etags oalib.so doc

