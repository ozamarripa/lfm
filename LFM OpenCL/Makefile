include make.def

lfm: ./Library/lfm.cpp
	$(CPPC) $^ $(OPENCLINC) $(OPENCLFLAGS) $(CPPFLAGS) $(LIBS) -o $@

clean:
	rm -f lfm
