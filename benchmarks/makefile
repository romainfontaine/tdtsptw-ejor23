ALL = Ari18Vu20 Rif20 Constant

all: $(ALL)

Ari18Vu20: Ari18Vu20/.dirstamp
Rif20: Rif20/.dirstamp
Constant: dum95/.dirstamp

Ari18Vu20/.dirstamp:
	tar -xvf Ari18Vu20.tar.gz && touch $@

Rif20/.dirstamp:
	tar -xvf Rif20.tar.gz && touch $@

dum95/.dirstamp:
	tar -xvf constant.tar.gz && touch $@

.PHONY: clean

clean:
	rm -rf Ari18Vu20 Rif20 pot96 lan93 dum95 gen98 ohl07 asc96 pes98 das10

