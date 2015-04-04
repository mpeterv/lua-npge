# lua-npge, Nucleotide PanGenome explorer (Lua module)
# Copyright (C) 2014-2015 Boris Nagaev
# See the LICENSE file for terms of use.

BUSTED := $(wildcard /usr/local/lib/luarocks/rocks/busted/*/bin/busted)

clean:
	rm -f luacov.*.out
	rm -f *.gcov
	rm -f src/npge/*/*.gcov
	rm -f src/npge/*/*.gcda
	rm -f src/npge/*/*.gcno
	rm -f src/npge/*/*.o
	rm -f -r npge

test:
	busted -c
	gcov src/npge/*/*.c
