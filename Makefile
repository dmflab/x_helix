all: erase compile

erase:
	rm -f x-helix 
compile:
	cc -o x-helix x_helix.c -lm

