bin_SCRIPTS = concentrations.pl concentrations-same.pl hybrid2.pl hybrid-2s.pl hybrid-ss-2s.pl hybrid-select.pl melt.pl plot2ann.pl h-num.pl ss-count.pl ct-energy-det.pl UNAFold.pl ITC-plot.pl vantHoff-plot.pl
DISTCLEANFILES = concentrations.pl concentrations-same.pl hybrid2.pl hybrid-2s.pl hybrid-ss-2s.pl hybrid-select.pl melt.pl plot2ann.pl h-num.pl ss-count.pl ct-energy-det.pl UNAFold.pl ITC-plot.pl vantHoff-plot.pl

install-exec-hook :
	cd $(DESTDIR)$(bindir) && $(LN_S) -f hybrid2.pl hybrid2-x.pl
	cd $(DESTDIR)$(bindir) && $(LN_S) -f hybrid2.pl hybrid2-min.pl
	cd $(DESTDIR)$(bindir) && $(LN_S) -f hybrid2.pl hybrid2-min-x.pl
	cd $(DESTDIR)$(bindir) && $(LN_S) -f hybrid2.pl hybrid2-2s.pl
	cd $(DESTDIR)$(bindir) && $(LN_S) -f hybrid2.pl hybrid2-2s-x.pl
