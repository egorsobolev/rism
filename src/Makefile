include config.mk

all:	librism.a straj mkhist mkavgw awrism

clean:  clean_librism.a clean_straj clean_mkhist clean_mkavgw clean_awrism

librism.a:
	cd rism && make

straj:  librism.a
	cd straj && make

mkhist: librism.a
	cd dhist && make

mkavgw:	librism.a
	cd avgw && make

awrism: librism.a
	cd awrism && make

clean_librism.a:
	cd rism && make clean

clean_straj:
	cd straj && make clean

clean_mkhist:
	cd dhist && make clean

clean_mkavgw:
	cd avgw && make clean

clean_awrism:
	cd awrism && make clean
