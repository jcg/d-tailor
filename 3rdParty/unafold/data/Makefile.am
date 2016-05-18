pkgdata_DATA = $(rna) $(dna) $(dat) $(nul) $(nuld) $(nuldat)
CLEANFILES = $(rna) $(dna) $(dat) $(nul) $(nuld) $(nuldat)

rna = asint1x2.DG asint1x2.DH dangle.DG dangle.DH loop.DG loop.DH miscloop.DG miscloop.DH sint2.DG sint2.DH sint4.DG sint4.DH stack.DG stack.DH tstackh.DG tstackh.DH tstacki.DG tstacki.DH tstackm.DG tstackm.DH tstacke.DG tstacke.DH triloop.DG triloop.DH tloop.DG tloop.DH hexaloop.DG hexaloop.DH

dna = asint1x2.DGD asint1x2.DHD dangle.DGD dangle.DHD loop.DGD loop.DHD miscloop.DGD miscloop.DHD sint2.DGD sint2.DHD sint4.DGD sint4.DHD stack.DGD stack.DHD tstackh.DGD tstackh.DHD tstacki.DGD tstacki.DHD tstackm.DGD tstackm.DHD tstacke.DGD tstacke.DHD triloop.DGD triloop.DHD tloop.DGD tloop.DHD hexaloop.DGD hexaloop.DHD

dat = asint1x2.DAT dangle.DAT loop.DAT miscloop.DAT sint2.DAT sint4.DAT stack.DAT tstackh.DAT tstacki.DAT tstackm.DAT tstacke.DAT triloop.DAT tloop.DAT hexaloop.DAT extinction.DAT

nul = asint1x2.NUL dangle.NUL loop.NUL miscloop.NUL sint2.NUL sint4.NUL stack.NUL tstackh.NUL tstacki.NUL tstackm.NUL tstacke.NUL triloop.NUL tloop.NUL hexaloop.NUL

nuld = asint1x2.NULD dangle.NULD loop.NULD miscloop.NULD sint2.NULD sint4.NULD stack.NULD tstackh.NULD tstacki.NULD tstackm.NULD tstacke.NULD triloop.NULD tloop.NULD hexaloop.NULD

nuldat = asint1x2.NULDAT dangle.NULDAT loop.NULDAT miscloop.NULDAT sint2.NULDAT sint4.NULDAT stack.NULDAT tstackh.NULDAT tstacki.NULDAT tstackm.NULDAT tstacke.NULDAT triloop.NULDAT tloop.NULDAT hexaloop.NULDAT

asint1x2.DG : $(srcdir)/rules/asint1x2.dg
	$(PERL) $(top_srcdir)/machine.pl $(INFINITY) $? $@

asint1x2.DH : $(srcdir)/rules/asint1x2.dh
	$(PERL) $(top_srcdir)/machine.pl 0 $? $@

dangle.DG : $(srcdir)/rules/dangle.dg
	$(PERL) $(top_srcdir)/machine.pl $(INFINITY) $? $@

dangle.DH : $(srcdir)/rules/dangle.dh
	$(PERL) $(top_srcdir)/machine.pl 0 $? $@

hexaloop.DG : $(srcdir)/rules/hexaloop.dg
	$(PERL) $(top_srcdir)/machine.pl $(INFINITY) $? $@

hexaloop.DH : $(srcdir)/rules/hexaloop.dh
	$(PERL) $(top_srcdir)/machine.pl 0 $? $@

loop.DG : $(srcdir)/rules/loop.dg
	$(PERL) $(top_srcdir)/machine.pl $(INFINITY) $? $@

loop.DH : $(srcdir)/rules/loop.dh
	$(PERL) $(top_srcdir)/machine.pl 0 $? $@

miscloop.DG : $(srcdir)/rules/miscloop.dg
	$(PERL) $(top_srcdir)/machine.pl $(INFINITY) $? $@

miscloop.DH : $(srcdir)/rules/miscloop.dh
	$(PERL) $(top_srcdir)/machine.pl 0 $? $@

sint2.DG : $(srcdir)/rules/sint2.dg
	$(PERL) $(top_srcdir)/machine.pl --no-wc-mismatches --no-gu-mismatches $(INFINITY) $? $@

sint2.DH : $(srcdir)/rules/sint2.dh
	$(PERL) $(top_srcdir)/machine.pl --no-wc-mismatches --no-gu-mismatches 0 $? $@

sint4.DG : $(srcdir)/rules/sint4.dg
	$(PERL) $(top_srcdir)/machine.pl --no-wc-mismatches --no-gu-mismatches $(INFINITY) $? $@

sint4.DH : $(srcdir)/rules/sint4.dh
	$(PERL) $(top_srcdir)/machine.pl --no-wc-mismatches --no-gu-mismatches 0 $? $@

stack.DG : $(srcdir)/rules/stack.dg
	$(PERL) $(top_srcdir)/machine.pl $(INFINITY) $? $@

stack.DH : $(srcdir)/rules/stack.dh
	$(PERL) $(top_srcdir)/machine.pl 0 $? $@

tstackh.DG : $(srcdir)/rules/tstackh.dg
	$(PERL) $(top_srcdir)/machine.pl $(INFINITY) $? $@

tstackh.DH : $(srcdir)/rules/tstackh.dh
	$(PERL) $(top_srcdir)/machine.pl 0 $? $@

tstacki.DG : $(srcdir)/rules/tstacki.dg
	$(PERL) $(top_srcdir)/machine.pl $(INFINITY) $? $@

tstacki.DH : $(srcdir)/rules/tstacki.dh
	$(PERL) $(top_srcdir)/machine.pl 0 $? $@

tstackm.DG : $(srcdir)/rules/tstackm.dg
	$(PERL) $(top_srcdir)/machine.pl $(INFINITY) $? $@

tstackm.DH : $(srcdir)/rules/tstackm.dh
	$(PERL) $(top_srcdir)/machine.pl 0 $? $@

tstacke.DG : $(srcdir)/rules/tstacke.dg
	$(PERL) $(top_srcdir)/machine.pl $(INFINITY) $? $@

tstacke.DH : $(srcdir)/rules/tstacke.dh
	$(PERL) $(top_srcdir)/machine.pl 0 $? $@

tloop.DG : $(srcdir)/rules/tloop.dg
	$(PERL) $(top_srcdir)/machine.pl $(INFINITY) $? $@

tloop.DH : $(srcdir)/rules/tloop.dh
	$(PERL) $(top_srcdir)/machine.pl 0 $? $@

triloop.DG : $(srcdir)/rules/triloop.dg
	$(PERL) $(top_srcdir)/machine.pl $(INFINITY) $? $@

triloop.DH : $(srcdir)/rules/triloop.dh
	$(PERL) $(top_srcdir)/machine.pl 0 $? $@

asint1x2.DGD : $(srcdir)/rules/asint1x2.dgd
	$(PERL) $(top_srcdir)/machine.pl $(INFINITY) $? $@

asint1x2.DHD : $(srcdir)/rules/asint1x2.dhd
	$(PERL) $(top_srcdir)/machine.pl 0 $? $@

dangle.DGD : $(srcdir)/rules/dangle.dgd
	$(PERL) $(top_srcdir)/machine.pl $(INFINITY) $? $@

dangle.DHD : $(srcdir)/rules/dangle.dhd
	$(PERL) $(top_srcdir)/machine.pl 0 $? $@

hexaloop.DGD : $(srcdir)/rules/hexaloop.dgd
	$(PERL) $(top_srcdir)/machine.pl $(INFINITY) $? $@

hexaloop.DHD : $(srcdir)/rules/hexaloop.dhd
	$(PERL) $(top_srcdir)/machine.pl 0 $? $@

loop.DGD : $(srcdir)/rules/loop.dgd
	$(PERL) $(top_srcdir)/machine.pl $(INFINITY) $? $@

loop.DHD : $(srcdir)/rules/loop.dhd
	$(PERL) $(top_srcdir)/machine.pl 0 $? $@

miscloop.DGD : $(srcdir)/rules/miscloop.dgd
	$(PERL) $(top_srcdir)/machine.pl $(INFINITY) $? $@

miscloop.DHD : $(srcdir)/rules/miscloop.dhd
	$(PERL) $(top_srcdir)/machine.pl 0 $? $@

sint2.DGD : $(srcdir)/rules/sint2.dgd
	$(PERL) $(top_srcdir)/machine.pl --no-wc-mismatches --no-gu-mismatches --no-double-gu-close $(INFINITY) $? $@

sint2.DHD : $(srcdir)/rules/sint2.dhd
	$(PERL) $(top_srcdir)/machine.pl --no-wc-mismatches --no-gu-mismatches --no-double-gu-close 0 $? $@

sint4.DGD : $(srcdir)/rules/sint4.dgd
	$(PERL) $(top_srcdir)/machine.pl --no-wc-mismatches --no-gu-mismatches $(INFINITY) $? $@

sint4.DHD : $(srcdir)/rules/sint4.dhd
	$(PERL) $(top_srcdir)/machine.pl --no-wc-mismatches --no-gu-mismatches 0 $? $@

stack.DGD : $(srcdir)/rules/stack.dgd
	$(PERL) $(top_srcdir)/machine.pl $(INFINITY) $? $@

stack.DHD : $(srcdir)/rules/stack.dhd
	$(PERL) $(top_srcdir)/machine.pl 0 $? $@

tstackh.DGD : $(srcdir)/rules/tstackh.dgd
	$(PERL) $(top_srcdir)/machine.pl $(INFINITY) $? $@

tstackh.DHD : $(srcdir)/rules/tstackh.dhd
	$(PERL) $(top_srcdir)/machine.pl 0 $? $@

tstacki.DGD : $(srcdir)/rules/tstacki.dgd
	$(PERL) $(top_srcdir)/machine.pl $(INFINITY) $? $@

tstacki.DHD : $(srcdir)/rules/tstacki.dhd
	$(PERL) $(top_srcdir)/machine.pl 0 $? $@

tstackm.DGD : $(srcdir)/rules/tstackm.dgd
	$(PERL) $(top_srcdir)/machine.pl $(INFINITY) $? $@

tstackm.DHD : $(srcdir)/rules/tstackm.dhd
	$(PERL) $(top_srcdir)/machine.pl 0 $? $@

tstacke.DGD : $(srcdir)/rules/tstacke.dgd
	$(PERL) $(top_srcdir)/machine.pl $(INFINITY) $? $@

tstacke.DHD : $(srcdir)/rules/tstacke.dhd
	$(PERL) $(top_srcdir)/machine.pl 0 $? $@

tloop.DGD : $(srcdir)/rules/tloop.dgd
	$(PERL) $(top_srcdir)/machine.pl $(INFINITY) $? $@

tloop.DHD : $(srcdir)/rules/tloop.dhd
	$(PERL) $(top_srcdir)/machine.pl 0 $? $@

triloop.DGD : $(srcdir)/rules/triloop.dgd
	$(PERL) $(top_srcdir)/machine.pl $(INFINITY) $? $@

triloop.DHD : $(srcdir)/rules/triloop.dhd
	$(PERL) $(top_srcdir)/machine.pl 0 $? $@

asint1x2.DAT : $(srcdir)/rules/asint1x2.dat
	$(PERL) $(top_srcdir)/machine.pl $(INFINITY) $? $@

dangle.DAT : $(srcdir)/rules/dangle.dat
	$(PERL) $(top_srcdir)/machine.pl $(INFINITY) $? $@

hexaloop.DAT : $(srcdir)/rules/hexaloop.dat
	$(PERL) $(top_srcdir)/machine.pl $(INFINITY) $? $@

loop.DAT : $(srcdir)/rules/loop.dat
	$(PERL) $(top_srcdir)/machine.pl $(INFINITY) $? $@

miscloop.DAT : $(srcdir)/rules/miscloop.dat
	$(PERL) $(top_srcdir)/machine.pl $(INFINITY) $? $@

sint2.DAT : $(srcdir)/rules/sint2.dat
	$(PERL) $(top_srcdir)/machine.pl --no-wc-mismatches --no-gu-close $(INFINITY) $? $@

sint4.DAT : $(srcdir)/rules/sint4.dat
	$(PERL) $(top_srcdir)/machine.pl --no-wc-mismatches $(INFINITY) $? $@

stack.DAT : $(srcdir)/rules/stack.dat
	$(PERL) $(top_srcdir)/machine.pl --no-gu-gu-stacks $(INFINITY) $? $@

tstackh.DAT : $(srcdir)/rules/tstackh.dat
	$(PERL) $(top_srcdir)/machine.pl $(INFINITY) $? $@

tstacki.DAT : $(srcdir)/rules/tstacki.dat
	$(PERL) $(top_srcdir)/machine.pl $(INFINITY) $? $@

tstackm.DAT : $(srcdir)/rules/tstackm.dat
	$(PERL) $(top_srcdir)/machine.pl $(INFINITY) $? $@

tstacke.DAT : $(srcdir)/rules/tstacke.dat
	$(PERL) $(top_srcdir)/machine.pl $(INFINITY) $? $@

tloop.DAT : $(srcdir)/rules/tloop.dat
	$(PERL) $(top_srcdir)/machine.pl $(INFINITY) $? $@

triloop.DAT : $(srcdir)/rules/triloop.dat
	$(PERL) $(top_srcdir)/machine.pl $(INFINITY) $? $@

extinction.DAT : $(srcdir)/rules/extinction.dat
	$(PERL) $(top_srcdir)/machine.pl $(INFINITY) $? $@

asint1x2.NUL : $(srcdir)/rules/asint1x2.dg
	$(PERL) $(top_srcdir)/machine.pl --null-energies $(INFINITY) $? $@

dangle.NUL : $(srcdir)/rules/dangle.dg
	$(PERL) $(top_srcdir)/machine.pl --null-energies $(INFINITY) $? $@

hexaloop.NUL : $(srcdir)/rules/hexaloop.dg
	$(PERL) $(top_srcdir)/machine.pl --null-energies $(INFINITY) $? $@

loop.NUL : $(srcdir)/rules/loop.dg
	$(PERL) $(top_srcdir)/machine.pl --null-energies $(INFINITY) $? $@

miscloop.NUL : $(srcdir)/rules/miscloop.dg
	$(PERL) $(top_srcdir)/machine.pl --null-energies $(INFINITY) $? $@

sint2.NUL : $(srcdir)/rules/sint2.dg
	$(PERL) $(top_srcdir)/machine.pl --null-energies --no-wc-mismatches --no-gu-mismatches $(INFINITY) $? $@

sint4.NUL : $(srcdir)/rules/sint4.dg
	$(PERL) $(top_srcdir)/machine.pl --null-energies --no-wc-mismatches --no-gu-mismatches $(INFINITY) $? $@

stack.NUL : $(srcdir)/rules/stack.dg
	$(PERL) $(top_srcdir)/machine.pl --null-energies $(INFINITY) $? $@

tstackh.NUL : $(srcdir)/rules/tstackh.dg
	$(PERL) $(top_srcdir)/machine.pl --null-energies $(INFINITY) $? $@

tstacki.NUL : $(srcdir)/rules/tstacki.dg
	$(PERL) $(top_srcdir)/machine.pl --null-energies $(INFINITY) $? $@

tstackm.NUL : $(srcdir)/rules/tstackm.dg
	$(PERL) $(top_srcdir)/machine.pl --null-energies $(INFINITY) $? $@

tstacke.NUL : $(srcdir)/rules/tstacke.dg
	$(PERL) $(top_srcdir)/machine.pl --null-energies $(INFINITY) $? $@

tloop.NUL : $(srcdir)/rules/tloop.dg
	$(PERL) $(top_srcdir)/machine.pl --null-energies $(INFINITY) $? $@

triloop.NUL : $(srcdir)/rules/triloop.dg
	$(PERL) $(top_srcdir)/machine.pl --null-energies $(INFINITY) $? $@

asint1x2.NULD : $(srcdir)/rules/asint1x2.dgd
	$(PERL) $(top_srcdir)/machine.pl --null-energies $(INFINITY) $? $@

dangle.NULD : $(srcdir)/rules/dangle.dgd
	$(PERL) $(top_srcdir)/machine.pl --null-energies $(INFINITY) $? $@

hexaloop.NULD : $(srcdir)/rules/hexaloop.dgd
	$(PERL) $(top_srcdir)/machine.pl --null-energies $(INFINITY) $? $@

loop.NULD : $(srcdir)/rules/loop.dgd
	$(PERL) $(top_srcdir)/machine.pl --null-energies $(INFINITY) $? $@

miscloop.NULD : $(srcdir)/rules/miscloop.dgd
	$(PERL) $(top_srcdir)/machine.pl --null-energies $(INFINITY) $? $@

sint2.NULD : $(srcdir)/rules/sint2.dgd
	$(PERL) $(top_srcdir)/machine.pl --null-energies --no-wc-mismatches --no-gu-mismatches $(INFINITY) $? $@

sint4.NULD : $(srcdir)/rules/sint4.dgd
	$(PERL) $(top_srcdir)/machine.pl --null-energies --no-wc-mismatches --no-gu-mismatches $(INFINITY) $? $@

stack.NULD : $(srcdir)/rules/stack.dgd
	$(PERL) $(top_srcdir)/machine.pl --null-energies $(INFINITY) $? $@

tstackh.NULD : $(srcdir)/rules/tstackh.dgd
	$(PERL) $(top_srcdir)/machine.pl --null-energies $(INFINITY) $? $@

tstacki.NULD : $(srcdir)/rules/tstacki.dgd
	$(PERL) $(top_srcdir)/machine.pl --null-energies $(INFINITY) $? $@

tstackm.NULD : $(srcdir)/rules/tstackm.dgd
	$(PERL) $(top_srcdir)/machine.pl --null-energies $(INFINITY) $? $@

tstacke.NULD : $(srcdir)/rules/tstacke.dgd
	$(PERL) $(top_srcdir)/machine.pl --null-energies $(INFINITY) $? $@

tloop.NULD : $(srcdir)/rules/tloop.dgd
	$(PERL) $(top_srcdir)/machine.pl --null-energies $(INFINITY) $? $@

triloop.NULD : $(srcdir)/rules/triloop.dgd
	$(PERL) $(top_srcdir)/machine.pl --null-energies $(INFINITY) $? $@

asint1x2.NULDAT : $(srcdir)/rules/asint1x2.dat
	$(PERL) $(top_srcdir)/machine.pl --null-energies $(INFINITY) $? $@

dangle.NULDAT : $(srcdir)/rules/dangle.dat
	$(PERL) $(top_srcdir)/machine.pl --null-energies $(INFINITY) $? $@

hexaloop.NULDAT : $(srcdir)/rules/hexaloop.dat
	$(PERL) $(top_srcdir)/machine.pl --null-energies $(INFINITY) $? $@

loop.NULDAT : $(srcdir)/rules/loop.dat
	$(PERL) $(top_srcdir)/machine.pl --null-energies $(INFINITY) $? $@

miscloop.NULDAT : $(srcdir)/rules/miscloop.dat
	$(PERL) $(top_srcdir)/machine.pl --null-energies $(INFINITY) $? $@

sint2.NULDAT : $(srcdir)/rules/sint2.dat
	$(PERL) $(top_srcdir)/machine.pl --null-energies --no-wc-mismatches --no-gu-close $(INFINITY) $? $@

sint4.NULDAT : $(srcdir)/rules/sint4.dat
	$(PERL) $(top_srcdir)/machine.pl --null-energies --no-wc-mismatches $(INFINITY) $? $@

stack.NULDAT : $(srcdir)/rules/stack.dat
	$(PERL) $(top_srcdir)/machine.pl --null-energies --no-gu-gu-stacks $(INFINITY) $? $@

tstackh.NULDAT : $(srcdir)/rules/tstackh.dat
	$(PERL) $(top_srcdir)/machine.pl --null-energies $(INFINITY) $? $@

tstacki.NULDAT : $(srcdir)/rules/tstacki.dat
	$(PERL) $(top_srcdir)/machine.pl --null-energies $(INFINITY) $? $@

tstackm.NULDAT : $(srcdir)/rules/tstackm.dat
	$(PERL) $(top_srcdir)/machine.pl --null-energies $(INFINITY) $? $@

tstacke.NULDAT : $(srcdir)/rules/tstacke.dat
	$(PERL) $(top_srcdir)/machine.pl --null-energies $(INFINITY) $? $@

tloop.NULDAT : $(srcdir)/rules/tloop.dat
	$(PERL) $(top_srcdir)/machine.pl --null-energies $(INFINITY) $? $@

triloop.NULDAT : $(srcdir)/rules/triloop.dat
	$(PERL) $(top_srcdir)/machine.pl --null-energies $(INFINITY) $? $@
