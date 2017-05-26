# Build seapy modules, tag file, and documentation
.PHONY: gtags
etags:
	rm -f TAGS
	find . '!' '('  -name __init__.py ')' -name "*\.py" | xargs etags -a --language=python

gtags:
	gtags -i --gtagslabel=pygments --statistics

oalib.so: src/oalib.F
	cd src; $(MAKE) $(MFLAGS) oalib

hindices.so: src/hindices.F
	cd src; $(MAKE) $(MFLAGS) hindices

doc: force_build
	cd doc; $(MAKE) $(MFLAGS) html && ln -fs _build/html html

force_build:
	true

clean:
	rm -f TAGS GPATH GRTAGS GTAGS oalib.so model/hindices.so doc/html

all: oalib.so hindices.so doc

